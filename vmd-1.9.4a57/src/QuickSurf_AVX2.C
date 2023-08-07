/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2019 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: QuickSurf_AVX2.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.10 $      $Date: 2021/11/12 16:36:33 $
 *
 ***************************************************************************/
/**
 * \file QuickSurf_AVX2.C
 * \brief SIMD vectorized gaussian surface kernels for x86 CPUs with
 *        AVX2+FMA instructions.
 */

// The x86 AVX code path requires FMA and AVX2 integer instructions
// in order to achieve performance that actually beats SSE2.
#if VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__) || (defined(_MSC_VER) && (_MSC_VER >= 1916))
#include <immintrin.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "QuickSurf.h"
#include "Inform.h"
#include "utilities.h"
#include "WKFUtils.h"
#include "VolumetricData.h"

#include "VMDDisplayList.h"
#include "Displayable.h"
#include "DispCmds.h"

#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))


/*
 * David J. Hardy
 * 12 Dec 2008
 *
 * aexpfnx() - Approximate expf() for negative x.
 *
 * Assumes that x <= 0.
 *
 * Assumes IEEE format for single precision float, specifically:
 * 1 sign bit, 8 exponent bits biased by 127, and 23 mantissa bits.
 *
 * Interpolates exp() on interval (-1/log2(e), 0], then shifts it by
 * multiplication of a fast calculation for 2^(-N).  The interpolation
 * uses a linear blending of 3rd degree Taylor polynomials at the end
 * points, so the approximation is once differentiable.
 *
 * The error is small (max relative error per interval is calculated
 * to be 0.131%, with a max absolute error of -0.000716).
 *
 * The cutoff is chosen so as to speed up the computation by early
 * exit from function, with the value chosen to give less than the
 * the max absolute error.  Use of a cutoff is unnecessary, except
 * for needing to shift smallest floating point numbers to zero,
 * i.e. you could remove cutoff and replace by:
 *
 * #define MINXNZ  -88.0296919311130  // -127 * log(2)
 *
 *   if (x < MINXNZ) return 0.f;
 *
 * Use of a cutoff causes a discontinuity which can be eliminated
 * through the use of a switching function.
 *
 * We can obtain arbitrarily smooth approximation by taking k+1 nodes on
 * the interval and weighting their respective Taylor polynomials by the
 * kth order Lagrange interpolant through those nodes.  The wiggle in the
 * polynomial interpolation due to equidistant nodes (Runge's phenomenon)
 * can be reduced by using Chebyshev nodes.
 */

#if defined(__GNUC__) && ! defined(__INTEL_COMPILER)
#define __align(X)  __attribute__((aligned(X) ))
#if (__GNUC__ < 4)
#define MISSING_mm_cvtsd_f64
#endif
#else
#define __align(X) __declspec(align(X) )
#endif

#define MLOG2EF    -1.44269504088896f

/*
 * Interpolating coefficients for linear blending of the
 * 3rd degree Taylor expansion of 2^x about 0 and -1.
 */
#define SCEXP0     1.0000000000000000f
#define SCEXP1     0.6987082824680118f
#define SCEXP2     0.2633174272827404f
#define SCEXP3     0.0923611991471395f
#define SCEXP4     0.0277520543324108f

/* for single precision float */
#define EXPOBIAS   127
#define EXPOSHIFT   23


/* cutoff is optional, but can help avoid unnecessary work */
#define ACUTOFF    -10

typedef union flint_t {
  float f;
  int n;
} flint;

#if VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__)
// AVX variant of the 'flint' union above
typedef union AVXreg_t {
  __m256  f;  // 8x float (AVX)
  __m256i i;  // 8x 32-bit int (AVX)
  struct {
    float r0, r1, r2, r3, r4, r5, r6, r7;  // get the individual registers
  } floatreg;
} AVXreg;
#endif



void vmd_gaussdensity_avx2(int verbose,
                           int natoms, const float *xyzr,
                           const float *atomicnum,
                           const float *colors,
                           float *densitymap, float *voltexmap,
                           const int *numvoxels,
                           float radscale, float gridspacing,
                           float isovalue, float gausslim) {
  int i, x, y, z;
  int maxvoxel[3];
  maxvoxel[0] = numvoxels[0]-1;
  maxvoxel[1] = numvoxels[1]-1;
  maxvoxel[2] = numvoxels[2]-1;
  const float invgridspacing = 1.0f / gridspacing;

#if VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__)
  // Variables for AVX2 optimized inner loop
  __m256 gridspacing8_8 = _mm256_set1_ps(gridspacing * 8.0f);
  __attribute__((aligned(16))) float sxdelta8[8]; // 16-byte aligned for AVX2
  for (x=0; x<8; x++)
    sxdelta8[x] = ((float) x) * gridspacing;
#endif

  // compute colors only if necessary, since they are costly
  if (voltexmap != NULL) {
    float invisovalue = 1.0f / isovalue;
    // compute both density map and floating point color texture map
    for (i=0; i<natoms; i++) {
      if (verbose && ((i & 0x3fff) == 0)) {
        printf(".");
        fflush(stdout);
      }

      ptrdiff_t ind = i*4L;
      float scaledrad = xyzr[ind + 3] * radscale;

      // MDFF atomic number weighted density factor
      float atomicnumfactor = 1.0f;
      if (atomicnum != NULL) {
        atomicnumfactor = atomicnum[i];
      }

      // negate, precompute reciprocal, and change to base 2 from the outset
      float arinv = -(1.0f/(2.0f*scaledrad*scaledrad)) * MLOG2EF;
      float radlim = gausslim * scaledrad;
      float radlim2 = radlim * radlim; // cutoff test done in cartesian coords
      radlim *= invgridspacing;

#if VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__)
      __m256 atomicnumfactor_8 = _mm256_set1_ps(atomicnumfactor);
#if VMDUSESVMLEXP
      // Use of Intel's SVML requires changing the pre-scaling factor
      __m256 arinv_8 = _mm256_set1_ps(arinv * (2.718281828f/2.0f) / MLOG2EF);
#else
      __m256 arinv_8 = _mm256_set1_ps(arinv); // Use inlined exp approximation
#endif
#endif

      float tmp;
      tmp = xyzr[ind  ] * invgridspacing;
      int xmin = MAX((int) (tmp - radlim), 0);
      int xmax = MIN((int) (tmp + radlim), maxvoxel[0]);
      tmp = xyzr[ind+1] * invgridspacing;
      int ymin = MAX((int) (tmp - radlim), 0);
      int ymax = MIN((int) (tmp + radlim), maxvoxel[1]);
      tmp = xyzr[ind+2] * invgridspacing;
      int zmin = MAX((int) (tmp - radlim), 0);
      int zmax = MIN((int) (tmp + radlim), maxvoxel[2]);

      float dz = zmin*gridspacing - xyzr[ind+2];
      for (z=zmin; z<=zmax; z++,dz+=gridspacing) {
        float dy = ymin*gridspacing - xyzr[ind+1];
        for (y=ymin; y<=ymax; y++,dy+=gridspacing) {
          float dy2dz2 = dy*dy + dz*dz;

          // early-exit when outside the cutoff radius in the Y-Z plane
          if (dy2dz2 >= radlim2)
            continue;

          ptrdiff_t addr = ptrdiff_t(z * numvoxels[0]) * ptrdiff_t(numvoxels[1]) + ptrdiff_t(y * numvoxels[0]);
          float dx = xmin*gridspacing - xyzr[ind];
          x=xmin;

#if VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__)
          // Use AVX when we have a multiple-of-8 to compute
          // finish all remaining density map points with SSE or regular non-SSE loop
          {
            __align(16) AVXreg n;
            __align(16) AVXreg y;
            __m256 dy2dz2_8 = _mm256_set1_ps(dy2dz2);
            __m256 dx_8 = _mm256_add_ps(_mm256_set1_ps(dx), _mm256_load_ps(&sxdelta8[0]));

            for (; (x+7)<=xmax; x+=8,dx_8=_mm256_add_ps(dx_8, gridspacing8_8)) {
              __m256 r2 = _mm256_fmadd_ps(dx_8, dx_8, dy2dz2_8);
              __m256 d;
#if VMDUSESVMLEXP
              // use Intel's SVML exp2() routine
              y.f = _mm256_exp2_ps(_mm256_mul_ps(r2, arinv_8));
#else
              // use our (much faster) fully inlined exponential approximation
              y.f = _mm256_mul_ps(r2, arinv_8);         /* already negated and in base 2 */
              n.i = _mm256_cvttps_epi32(y.f);
              d = _mm256_cvtepi32_ps(n.i);
              d = _mm256_sub_ps(d, y.f);

              // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
              // Perform Horner's method to evaluate interpolating polynomial.
              y.f = _mm256_fmadd_ps(d, _mm256_set1_ps(SCEXP4), _mm256_set1_ps(SCEXP3));
              y.f = _mm256_fmadd_ps(y.f, d, _mm256_set1_ps(SCEXP2));
              y.f = _mm256_fmadd_ps(y.f, d, _mm256_set1_ps(SCEXP1));
              y.f = _mm256_fmadd_ps(y.f, d, _mm256_set1_ps(SCEXP0));

              // Calculate 2^N exactly by directly manipulating floating point exponent,
              // then use it to scale y for the final result.
              // We need AVX2 instructions to be able to operate on
              // 8-wide integer types efficiently.
              n.i = _mm256_sub_epi32(_mm256_set1_epi32(EXPOBIAS), n.i);
              n.i = _mm256_slli_epi32(n.i, EXPOSHIFT);
              y.f = _mm256_mul_ps(y.f, n.f);
              y.f = _mm256_mul_ps(y.f, atomicnumfactor_8); // MDFF density maps
#endif

              // At present, we do unaligned loads/stores since we can't guarantee
              // that the X-dimension is always a multiple of 8.
              float *ufptr = &densitymap[addr + x];
              d = _mm256_loadu_ps(ufptr);
              _mm256_storeu_ps(ufptr, _mm256_add_ps(d, y.f));

              // Accumulate density-weighted color to texture map.
              // Pre-multiply colors by the inverse isovalue we will extract
              // the surface on, to cause the final color to be normalized.
              d = _mm256_mul_ps(y.f, _mm256_set1_ps(invisovalue));
              ptrdiff_t caddr = (addr + x) * 3L;

#if 0 && VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__)
              // convert rgb3f AOS format to 8-element SOA vectors using shuffle instructions
              float *txptr = &voltexmap[caddr];
              __m128 *m = (__m128 *) txptr;
              // unaligned load of 8 consecutive rgb3f texture map texels
              __m256 m03 = _mm256_castps128_ps256(m[0]); // load lower halves
              __m256 m14 = _mm256_castps128_ps256(m[1]);
              __m256 m25 = _mm256_castps128_ps256(m[2]);
              m03  = _mm256_insertf128_ps(m03, m[3], 1); // load upper halves
              m14  = _mm256_insertf128_ps(m14, m[4], 1);
              m25  = _mm256_insertf128_ps(m25, m[5], 1);

              // upper Rs and Gs
              __m256 rg = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE(2,1,3,2));
              // lower Gs and Bs
              __m256 gb = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE(1,0,2,1));
              __m256 r  = _mm256_shuffle_ps(m03, rg , _MM_SHUFFLE(2,0,3,0));
              __m256 g  = _mm256_shuffle_ps(gb , rg , _MM_SHUFFLE(3,1,2,0));
              __m256 b  = _mm256_shuffle_ps(gb , m25, _MM_SHUFFLE(3,0,3,1));

              // accumulate density-scaled colors into texels
              r = _mm256_fmadd_ps(d, _mm256_set1_ps(colors[ind    ]), r);
              g = _mm256_fmadd_ps(d, _mm256_set1_ps(colors[ind + 1]), g);
              b = _mm256_fmadd_ps(d, _mm256_set1_ps(colors[ind + 2]), b);

              // convert 8-element SOA vectors to rgb3f AOS format using shuffle instructions
              __m256 rrg = _mm256_shuffle_ps(r, g, _MM_SHUFFLE(2,0,2,0));
              __m256 rgb = _mm256_shuffle_ps(g, b, _MM_SHUFFLE(3,1,3,1));
              __m256 rbr = _mm256_shuffle_ps(b, r, _MM_SHUFFLE(3,1,2,0));
              __m256 r03 = _mm256_shuffle_ps(rrg, rbr, _MM_SHUFFLE(2,0,2,0));
              __m256 r14 = _mm256_shuffle_ps(rgb, rrg, _MM_SHUFFLE(3,1,2,0));
              __m256 r25 = _mm256_shuffle_ps(rbr, rgb, _MM_SHUFFLE(3,1,3,1));

              // unaligned store of consecutive rgb3f texture map texels
              m[0] = _mm256_castps256_ps128( r03 );
              m[1] = _mm256_castps256_ps128( r14 );
              m[2] = _mm256_castps256_ps128( r25 );
              m[3] = _mm256_extractf128_ps( r03 ,1);
              m[4] = _mm256_extractf128_ps( r14 ,1);
              m[5] = _mm256_extractf128_ps( r25 ,1);
#else
              // color by atom colors
              float r, g, b;
              r = colors[ind    ];
              g = colors[ind + 1];
              b = colors[ind + 2];

              AVXreg tmp;
              tmp.f = d;
              float density;
              density = tmp.floatreg.r0;
              voltexmap[caddr     ] += density * r;
              voltexmap[caddr +  1] += density * g;
              voltexmap[caddr +  2] += density * b;

              density = tmp.floatreg.r1;
              voltexmap[caddr +  3] += density * r;
              voltexmap[caddr +  4] += density * g;
              voltexmap[caddr +  5] += density * b;

              density = tmp.floatreg.r2;
              voltexmap[caddr +  6] += density * r;
              voltexmap[caddr +  7] += density * g;
              voltexmap[caddr +  8] += density * b;

              density = tmp.floatreg.r3;
              voltexmap[caddr +  9] += density * r;
              voltexmap[caddr + 10] += density * g;
              voltexmap[caddr + 11] += density * b;

              density = tmp.floatreg.r4;
              voltexmap[caddr + 12] += density * r;
              voltexmap[caddr + 13] += density * g;
              voltexmap[caddr + 14] += density * b;

              density = tmp.floatreg.r5;
              voltexmap[caddr + 15] += density * r;
              voltexmap[caddr + 16] += density * g;
              voltexmap[caddr + 17] += density * b;

              density = tmp.floatreg.r6;
              voltexmap[caddr + 18] += density * r;
              voltexmap[caddr + 19] += density * g;
              voltexmap[caddr + 20] += density * b;

              density = tmp.floatreg.r7;
              voltexmap[caddr + 21] += density * r;
              voltexmap[caddr + 22] += density * g;
              voltexmap[caddr + 23] += density * b;
#endif
            }
          }
#endif

          // finish all remaining density map points with regular non-SSE loop
          for (; x<=xmax; x++,dx+=gridspacing) {
            float r2 = dx*dx + dy2dz2;

            // use our (much faster) fully inlined exponential approximation
            float mb = r2 * arinv;         /* already negated and in base 2 */
            int mbflr = (int) mb;          /* get int part, floor() */
            float d = mbflr - mb;          /* remaining exponent, -1 < d <= 0 */

            /* approx with linear blend of Taylor polys */
            float sy = SCEXP0 + d*(SCEXP1 + d*(SCEXP2 + d*(SCEXP3 + d*SCEXP4)));

            /* 2^(-mbflr) */
            flint scalfac;
            scalfac.n = (EXPOBIAS - mbflr) << EXPOSHIFT;

            // XXX assume we are never beyond the cutoff value in this loop
            float density = (sy * scalfac.f);

            density *= atomicnumfactor; // MDFF Cryo-EM atomic number density

            // accumulate density value to density map
            densitymap[addr + x] += density;

            // Accumulate density-weighted color to texture map.
            // Pre-multiply colors by the inverse isovalue we will extract
            // the surface on, to cause the final color to be normalized.
            density *= invisovalue;
            ptrdiff_t caddr = (addr + x) * 3L;

            // color by atom colors
            voltexmap[caddr    ] += density * colors[ind    ];
            voltexmap[caddr + 1] += density * colors[ind + 1];
            voltexmap[caddr + 2] += density * colors[ind + 2];
          }
        }
      }
    }
  } else {
    // compute density map only
    for (i=0; i<natoms; i++) {
      if (verbose && ((i & 0x3fff) == 0)) {
        printf(".");
        fflush(stdout);
      }

      ptrdiff_t ind = i*4L;
      float scaledrad = xyzr[ind+3] * radscale;

      // MDFF atomic number weighted density factor
      float atomicnumfactor = 1.0f;
      if (atomicnum != NULL) {
        atomicnumfactor = atomicnum[i];
      }

      // negate, precompute reciprocal, and change to base 2 from the outset
      float arinv = -(1.0f/(2.0f*scaledrad*scaledrad)) * MLOG2EF;
      float radlim = gausslim * scaledrad;
      float radlim2 = radlim * radlim; // cutoff test done in cartesian coords
      radlim *= invgridspacing;

#if VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__)
      __m256 atomicnumfactor_8 = _mm256_set1_ps(atomicnumfactor);
#if VMDUSESVMLEXP
      // Use of Intel's SVML requires changing the pre-scaling factor
      __m256 arinv_8 = _mm256_set1_ps(arinv * (2.718281828f/2.0f) / MLOG2EF);
#else
      __m256 arinv_8 = _mm256_set1_ps(arinv); // Use inlined exp approximation
#endif
#endif

      float tmp;
      tmp = xyzr[ind  ] * invgridspacing;
      int xmin = MAX((int) (tmp - radlim), 0);
      int xmax = MIN((int) (tmp + radlim), maxvoxel[0]);
      tmp = xyzr[ind+1] * invgridspacing;
      int ymin = MAX((int) (tmp - radlim), 0);
      int ymax = MIN((int) (tmp + radlim), maxvoxel[1]);
      tmp = xyzr[ind+2] * invgridspacing;
      int zmin = MAX((int) (tmp - radlim), 0);
      int zmax = MIN((int) (tmp + radlim), maxvoxel[2]);

      float dz = zmin*gridspacing - xyzr[ind+2];
      for (z=zmin; z<=zmax; z++,dz+=gridspacing) {
        float dy = ymin*gridspacing - xyzr[ind+1];
        for (y=ymin; y<=ymax; y++,dy+=gridspacing) {
          float dy2dz2 = dy*dy + dz*dz;

          // early-exit when outside the cutoff radius in the Y-Z plane
          if (dy2dz2 >= radlim2)
            continue;

          ptrdiff_t addr = ptrdiff_t(z * numvoxels[0]) * ptrdiff_t(numvoxels[1]) + ptrdiff_t(y * numvoxels[0]);
          float dx = xmin*gridspacing - xyzr[ind];
          x=xmin;

#if VMDUSEAVX2 && defined(__AVX__) && defined(__AVX2__)
          // Use AVX when we have a multiple-of-8 to compute
          // finish all remaining density map points with regular C++ loop
          {
            __align(16) AVXreg n;
            __align(16) AVXreg y;
            __m256 dy2dz2_8 = _mm256_set1_ps(dy2dz2);
            __m256 dx_8 = _mm256_add_ps(_mm256_set1_ps(dx), _mm256_load_ps(&sxdelta8[0]));

            for (; (x+7)<=xmax; x+=8,dx_8=_mm256_add_ps(dx_8, gridspacing8_8)) {
              __m256 r2 = _mm256_fmadd_ps(dx_8, dx_8, dy2dz2_8);
              __m256 d;
#if VMDUSESVMLEXP
              // use Intel's SVML exp2() routine
              y.f = _mm256_exp2_ps(_mm256_mul_ps(r2, arinv_8));
#else
              // use our (much faster) fully inlined exponential approximation
              y.f = _mm256_mul_ps(r2, arinv_8);         /* already negated and in base 2 */
              n.i = _mm256_cvttps_epi32(y.f);
              d = _mm256_cvtepi32_ps(n.i);
              d = _mm256_sub_ps(d, y.f);

              // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
              // Perform Horner's method to evaluate interpolating polynomial.
#if defined(__FMA__)
              // use FMA3 instructions when possible
              y.f = _mm256_fmadd_ps(d, _mm256_set1_ps(SCEXP4), _mm256_set1_ps(SCEXP3));
              y.f = _mm256_fmadd_ps(y.f, d, _mm256_set1_ps(SCEXP2));
              y.f = _mm256_fmadd_ps(y.f, d, _mm256_set1_ps(SCEXP1));
              y.f = _mm256_fmadd_ps(y.f, d, _mm256_set1_ps(SCEXP0));
#else
              y.f = _mm256_mul_ps(d, _mm256_set1_ps(SCEXP4));      /* for x^4 term */
              y.f = _mm256_add_ps(y.f, _mm256_set1_ps(SCEXP3));    /* for x^3 term */
              y.f = _mm256_mul_ps(y.f, d);
              y.f = _mm256_add_ps(y.f, _mm256_set1_ps(SCEXP2));    /* for x^2 term */
              y.f = _mm256_mul_ps(y.f, d);
              y.f = _mm256_add_ps(y.f, _mm256_set1_ps(SCEXP1));    /* for x^1 term */
              y.f = _mm256_mul_ps(y.f, d);
              y.f = _mm256_add_ps(y.f, _mm256_set1_ps(SCEXP0));    /* for x^0 term */
#endif

              // Calculate 2^N exactly by directly manipulating floating point exponent,
              // then use it to scale y for the final result.
              // We need AVX2 instructions to be able to operate on
              // 8-wide integer types efficiently.
              n.i = _mm256_sub_epi32(_mm256_set1_epi32(EXPOBIAS), n.i);
              n.i = _mm256_slli_epi32(n.i, EXPOSHIFT);
              y.f = _mm256_mul_ps(y.f, n.f);
              y.f = _mm256_mul_ps(y.f, atomicnumfactor_8); // MDFF density maps
#endif

              // At present, we do unaligned loads/stores since we can't guarantee
              // that the X-dimension is always a multiple of 8.
              float *ufptr = &densitymap[addr + x];
              d = _mm256_loadu_ps(ufptr);
              _mm256_storeu_ps(ufptr, _mm256_add_ps(d, y.f));
            }
          }
#endif

          // finish all remaining density map points with regular non-SSE loop
          for (; x<=xmax; x++,dx+=gridspacing) {
            float r2 = dx*dx + dy2dz2;

            // use our (much faster) fully inlined exponential approximation
            float mb = r2 * arinv;         /* already negated and in base 2 */
            int mbflr = (int) mb;          /* get int part, floor() */
            float d = mbflr - mb;          /* remaining exponent, -1 < d <= 0 */

            /* approx with linear blend of Taylor polys */
            float sy = SCEXP0 + d*(SCEXP1 + d*(SCEXP2 + d*(SCEXP3 + d*SCEXP4)));

            /* 2^(-mbflr) */
            flint scalfac;
            scalfac.n = (EXPOBIAS - mbflr) << EXPOSHIFT;

            // XXX assume we are never beyond the cutoff value in this loop
            float density = (sy * scalfac.f);

            density *= atomicnumfactor; // MDFF Cryo-EM atomic number density

            densitymap[addr + x] += density;
          }
        }
      }
    }
  }

  // prevent x86 AVX-SSE transition performance loss due to CPU state 
  // transition penalties or false dependence on upper register state
  _mm256_zeroupper();
}
 



