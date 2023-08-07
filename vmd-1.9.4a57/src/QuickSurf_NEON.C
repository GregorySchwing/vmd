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
 *      $RCSfile: QuickSurf_NEON.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $      $Date: 2021/11/12 16:36:33 $
 *
 ***************************************************************************/
/**
 * \file QuickSurf_NEON.C
 * \brief SIMD vectorized gaussian surface kernels for NEON instructions.
 */

// The NEON code path requires FMA instructions
// Due to differences in code generation between gcc/intelc/clang/msvc, we
// don't have to check for a defined(__NEON__)
#if defined(VMDCPUDISPATCH) && defined(VMDUSENEON)
#include <arm_neon.h>

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

// NEON variant of the 'flint' union above
typedef union NEONreg_t {
  float32x4_t f;  // 4x float (NEON)
  int32x4_t i;    // 4x 32-bit int (NEON)
  struct {
    float r0, r1, r2, r3;  // get the individual registers
  } floatreg;
} NEONreg;

void vmd_gaussdensity_neon(int verbose,
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

  // Variables for NEON optimized inner loop
  float32x4_t gridspacing4_4;
  __attribute__((aligned(16))) float sxdelta4[4]; // 16-byte aligned for 

  gridspacing4_4 = vdupq_n_f32(gridspacing * 4.0f);
  for (x=0; x<4; x++)
    sxdelta4[x] = ((float) x) * gridspacing;

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

      float32x4_t atomicnumfactor_4;
      float32x4_t arinv_4;
      atomicnumfactor_4 = vdupq_n_f32(atomicnumfactor);

      // Use our fully inlined exp approximation
      arinv_4 = vdupq_n_f32(arinv);

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

          // Use NEON when we have a multiple-of-4 to compute
          // finish all remaining density map points with regular C++ loop
          {
            __align(16) NEONreg n;
            __align(16) NEONreg y;
            float32x4_t dy2dz2_4 = vdupq_n_f32(dy2dz2);
            float32x4_t dx_4 = vaddq_f32(vdupq_n_f32(dx), vld1q_f32(&sxdelta4[0]));

            for (; (x+3)<=xmax; x+=4,dx_4=vaddq_f32(dx_4, gridspacing4_4)) {
              float32x4_t r2 = vfmaq_f32(dy2dz2_4, dx_4, dx_4);
              float32x4_t d;

              // use our (much faster) fully inlined exponential approximation
              y.f = vmulq_f32(r2, arinv_4);  // already negated and in base 2
              n.i = vcvtq_s32_f32(y.f);
              d = vcvtq_f32_s32(n.i);
              d = vsubq_f32(d, y.f);

              // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
              // Perform Horner's method to evaluate interpolating polynomial.
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP3), vdupq_n_f32(SCEXP4), d);
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP2), d, y.f);
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP1), d, y.f);
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP0), d, y.f);

              // Calculate 2^N exactly by directly manipulating floating point
              // exponent, then use it to scale y for the final result.
              n.i = vsubq_s32(vdupq_n_s32(EXPOBIAS), n.i);
              n.i = vshlq_s32(n.i, vdupq_n_s32(EXPOSHIFT));
              y.f = vmulq_f32(y.f, n.f);
              y.f = vmulq_f32(y.f, atomicnumfactor_4); // MDFF density maps

              float *ufptr = &densitymap[addr + x];
              d = vld1q_f32(ufptr); // ARM64 NEON isn't alignment-sensitive
              vst1q_f32(ufptr, vaddq_f32(d, y.f));

              // Accumulate density-weighted color to texture map.
              // Pre-multiply colors by the inverse isovalue we will extract
              // the surface on, to cause the final color to be normalized.
              d = vmulq_f32(y.f, vdupq_n_f32(invisovalue));
              ptrdiff_t caddr = (addr + x) * 3L;

#if 0
              // NEON AOS/SOA conversions
#else
              // color by atom colors
              float r, g, b;
              r = colors[ind    ];
              g = colors[ind + 1];
              b = colors[ind + 2];

              NEONreg tmp;
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
#endif
            }
          }

          // finish all remaining density map points with regular C++ loop
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

      float32x4_t atomicnumfactor_4 = vdupq_n_f32(atomicnumfactor);
      float32x4_t arinv_4= vdupq_n_f32(arinv); // Use inlined exp approximation

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

          // Use NEON when we have a multiple-of-4 to compute
          // finish all remaining density map points with C++ loop
          {
            __align(16) NEONreg n;
            __align(16) NEONreg y;
            float32x4_t dy2dz2_4 = vdupq_n_f32(dy2dz2);
            float32x4_t dx_4 = vaddq_f32(vdupq_n_f32(dx), vld1q_f32(&sxdelta4[0]));

            for (; (x+3)<=xmax; x+=4,dx_4=vaddq_f32(dx_4, gridspacing4_4)) {
              float32x4_t r2 = vfmaq_f32(dy2dz2_4, dx_4, dx_4);
              float32x4_t d;

              // use our (much faster) fully inlined exponential approximation
              y.f = vmulq_f32(r2, arinv_4);         /* already negated and in base 2 */
              n.i = vcvtq_s32_f32(y.f);
              d = vcvtq_f32_s32(n.i);
              d = vsubq_f32(d, y.f);

              // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
              // Perform Horner's method to evaluate interpolating polynomial.
#if __ARM_FEATURE_FMA
              // use FMA3 instructions when possible
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP3), vdupq_n_f32(SCEXP4), d);
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP2), d, y.f);
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP1), d, y.f);
              y.f = vfmaq_f32(vdupq_n_f32(SCEXP0), d, y.f);
#else
              y.f = vmulq_f32(d, vdupq_n_f32(SCEXP4));      /* for x^4 term */
              y.f = vaddq_f32(y.f, vdupq_n_f32(SCEXP3));    /* for x^3 term */
              y.f = vmulq_f32(y.f, d);
              y.f = vaddq_f32(y.f, vdupq_n_f32(SCEXP2));    /* for x^2 term */
              y.f = vmulq_f32(y.f, d);
              y.f = vaddq_f32(y.f, vdupq_n_f32(SCEXP1));    /* for x^1 term */
              y.f = vmulq_f32(y.f, d);
              y.f = vaddq_f32(y.f, vdupq_n_f32(SCEXP0));    /* for x^0 term */
#endif

              // Calculate 2^N exactly by directly manipulating floating point 
              // exponent, then use it to scale y for the final result.
              n.i = vsubq_s32(vdupq_n_s32(EXPOBIAS), n.i);
              n.i = vshlq_s32(n.i, vdupq_n_s32(EXPOSHIFT));
              y.f = vmulq_f32(y.f, n.f);
              y.f = vmulq_f32(y.f, atomicnumfactor_4); // MDFF density maps

              float *ufptr = &densitymap[addr + x];
              d = vld1q_f32(ufptr); // ARM64 NEON isn't alignment-sensitive
              vst1q_f32(ufptr, vaddq_f32(d, y.f));
            }
          }

          // finish all remaining density map points with regular C++ loop
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
}
 

#endif

