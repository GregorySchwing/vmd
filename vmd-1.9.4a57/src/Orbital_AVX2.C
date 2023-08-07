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
 *	$RCSfile: Orbital_AVX2.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 2021/02/16 17:24:26 $
 *
 ***************************************************************************/
/**
 * \file Orbital_AVX2.C
 * \brief SIMD vectorized MO kernels for AVX2 instructions.
 */

// Due to differences in code generation between gcc/intelc/clang/msvc, we
// don't have to check for a defined(__AVX2__)
#if defined(VMDCPUDISPATCH) && defined(VMDUSEAVX2) 

#include <immintrin.h>

#include <math.h>
#include <stdio.h>
#include "Orbital.h"
#include "DrawMolecule.h"
#include "utilities.h"
#include "Inform.h"
#include "WKFThreads.h"
#include "WKFUtils.h"
#include "ProfileHooks.h"

#define ANGS_TO_BOHR 1.88972612478289694072f

#if defined(__GNUC__) && ! defined(__INTEL_COMPILER)
#define __align(X)  __attribute__((aligned(X) ))
#else
#define __align(X) __declspec(align(X) )
#endif

#define MLOG2EF    -1.44269504088896f

#if 0
static void print_mm256_ps(__m256 v) {
  __attribute__((aligned(32))) float tmp[8]; // 32-byte aligned for AVX2
  _mm256_storeu_ps(&tmp[0], v);

  printf("mm256: ");
  int i;
  for (i=0; i<8; i++) 
    printf("%g ", tmp[i]);
  printf("\n");
}
#endif


//
// John Stone, January 2021
//
// aexpfnxavx2() - AVX2 version of aexpfnx().
//

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

typedef union AVX2reg_t {
  __m256  f;  // 8x float (AVX)
  __m256i i;  // 8x 32-bit int (AVX2)
} AVX2reg;

__m256 aexpfnxavx2(__m256 x) {
  __align(32) AVX2reg scal;
  scal.f = _mm256_cmp_ps(x, _mm256_set1_ps(ACUTOFF), _CMP_GE_OQ);  // Is x within cutoff?
  // If all x are outside of cutoff, return 0s.
  if (_mm256_movemask_ps(scal.f) == 0) {
    return _mm256_set1_ps(0.0f);
  }
  // Otherwise, scal.f contains mask to be ANDed with the scale factor

  /*
   * Convert base:  exp(x) = 2^(N-d) where N is integer and 0 <= d < 1.
   *
   * Below we calculate n=N and x=-d, with "y" for temp storage,
   * calculate floor of x*log2(e) and subtract to get -d.
   */
  __align(32) AVX2reg n;
  __m256 mb = _mm256_mul_ps(x, _mm256_set1_ps(MLOG2EF));
  n.i = _mm256_cvttps_epi32(mb);
  __m256 mbflr = _mm256_cvtepi32_ps(n.i);
  __m256 d = _mm256_sub_ps(mbflr, mb);

  // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
  // Perform Horner's method to evaluate interpolating polynomial.
  __m256 y;
  y = _mm256_fmadd_ps(d, _mm256_set1_ps(SCEXP4), _mm256_set1_ps(SCEXP3));
  y = _mm256_fmadd_ps(y, d, _mm256_set1_ps(SCEXP2));
  y = _mm256_fmadd_ps(y, d, _mm256_set1_ps(SCEXP1));
  y = _mm256_fmadd_ps(y, d, _mm256_set1_ps(SCEXP0));

  // Calculate 2^N exactly by directly manipulating floating point exponent,
  // then use it to scale y for the final result.
  n.i = _mm256_sub_epi32(_mm256_set1_epi32(EXPOBIAS), n.i);
  n.i = _mm256_slli_epi32(n.i, EXPOSHIFT);
  scal.f = _mm256_and_ps(scal.f, n.f);
  y = _mm256_mul_ps(y, scal.f);

  return y;
}


//
// AVX2 implementation for Xeons that don't have special fctn units
//
int evaluate_grid_avx2(int numatoms,
                          const float *wave_f, const float *basis_array,
                          const float *atompos,
                          const int *atom_basis,
                          const int *num_shells_per_atom,
                          const int *num_prim_per_shell,
                          const int *shell_types,
                          const int *numvoxels,
                          float voxelsize,
                          const float *origin,
                          int density,
                          float * orbitalgrid) {
  if (!orbitalgrid)
    return -1;

  int nx, ny, nz;
  __attribute__((aligned(32))) float sxdelta[8]; // 32-byte aligned for AVX2
  for (nx=0; nx<8; nx++) 
    sxdelta[nx] = ((float) nx) * voxelsize * ANGS_TO_BOHR;

  // Calculate the value of the orbital at each gridpoint and store in 
  // the current oribtalgrid array
  int numgridxy = numvoxels[0]*numvoxels[1];
  for (nz=0; nz<numvoxels[2]; nz++) {
    float grid_x, grid_y, grid_z;
    grid_z = origin[2] + nz * voxelsize;
    for (ny=0; ny<numvoxels[1]; ny++) {
      grid_y = origin[1] + ny * voxelsize;
      int gaddrzy = ny*numvoxels[0] + nz*numgridxy;
      for (nx=0; nx<numvoxels[0]; nx+=8) {
        grid_x = origin[0] + nx * voxelsize;

        // calculate the value of the wavefunction of the
        // selected orbital at the current grid point
        int at;
        int prim, shell;

        // initialize value of orbital at gridpoint
        __m256 value = _mm256_set1_ps(0.0f);

        // initialize the wavefunction and shell counters
        int ifunc = 0; 
        int shell_counter = 0;

        // loop over all the QM atoms
        for (at=0; at<numatoms; at++) {
          int maxshell = num_shells_per_atom[at];
          int prim_counter = atom_basis[at];

          // calculate distance between grid point and center of atom
          float sxdist = (grid_x - atompos[3*at  ])*ANGS_TO_BOHR;
          float sydist = (grid_y - atompos[3*at+1])*ANGS_TO_BOHR;
          float szdist = (grid_z - atompos[3*at+2])*ANGS_TO_BOHR;

          float sydist2 = sydist*sydist;
          float szdist2 = szdist*szdist;
          float yzdist2 = sydist2 + szdist2;

          __m256 xdelta = _mm256_load_ps(&sxdelta[0]); // aligned load
          __m256 xdist  = _mm256_set1_ps(sxdist);
          xdist = _mm256_add_ps(xdist, xdelta);
          __m256 ydist  = _mm256_set1_ps(sydist);
          __m256 zdist  = _mm256_set1_ps(szdist);
          __m256 xdist2 = _mm256_mul_ps(xdist, xdist);
          __m256 ydist2 = _mm256_mul_ps(ydist, ydist);
          __m256 zdist2 = _mm256_mul_ps(zdist, zdist);
          __m256 dist2  = _mm256_set1_ps(yzdist2); 
          dist2 = _mm256_add_ps(dist2, xdist2);
 
          // loop over the shells belonging to this atom
          // XXX this is maybe a misnomer because in split valence
          //     basis sets like 6-31G we have more than one basis
          //     function per (valence-)shell and we are actually
          //     looping over the individual contracted GTOs
          for (shell=0; shell < maxshell; shell++) {
            __m256 contracted_gto = _mm256_set1_ps(0.0f);

            // Loop over the Gaussian primitives of this contracted 
            // basis function to build the atomic orbital
            // 
            // XXX there's a significant opportunity here for further
            //     speedup if we replace the entire set of primitives
            //     with the single gaussian that they are attempting 
            //     to model.  This could give us another 6x speedup in 
            //     some of the common/simple cases.
            int maxprim = num_prim_per_shell[shell_counter];
            int shelltype = shell_types[shell_counter];
            for (prim=0; prim<maxprim; prim++) {
              // XXX pre-negate exponent value
              float exponent       = -basis_array[prim_counter    ];
              float contract_coeff =  basis_array[prim_counter + 1];

              // contracted_gto += contract_coeff * exp(-exponent*dist2);
              __m256 expval = _mm256_mul_ps(_mm256_set1_ps(exponent), dist2);
              // exp2f() equivalent required, use base-2 approximation
              __m256 retval = aexpfnxavx2(expval);
              contracted_gto = _mm256_fmadd_ps(_mm256_set1_ps(contract_coeff), retval, contracted_gto);

              prim_counter += 2;
            }

            /* multiply with the appropriate wavefunction coefficient */
            __m256 tmpshell = _mm256_set1_ps(0.0f);
            switch (shelltype) {
              // use FMADD instructions
              case S_SHELL:
                value = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), contracted_gto, value);
                break;

              case P_SHELL:
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), xdist, tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), ydist, tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), zdist, tmpshell);
                value = _mm256_fmadd_ps(tmpshell, contracted_gto, value);
                break;

              case D_SHELL:
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), xdist2, tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(xdist, ydist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), ydist2, tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(xdist, zdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(ydist, zdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), zdist2, tmpshell);
                value = _mm256_fmadd_ps(tmpshell, contracted_gto, value);
                break;

              case F_SHELL:
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(xdist2, xdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(xdist2, ydist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(ydist2, xdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(ydist2, ydist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(xdist2, zdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(_mm256_mul_ps(xdist, ydist), zdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(ydist2, zdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(zdist2, xdist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(zdist2, ydist), tmpshell);
                tmpshell = _mm256_fmadd_ps(_mm256_set1_ps(wave_f[ifunc++]), _mm256_mul_ps(zdist2, zdist), tmpshell);
                value = _mm256_fmadd_ps(tmpshell, contracted_gto, value);
                break;
 
#if 0
              default:
                // avoid unnecessary branching and minimize use of pow()
                int i, j; 
                float xdp, ydp, zdp;
                float xdiv = 1.0f / xdist;
                for (j=0, zdp=1.0f; j<=shelltype; j++, zdp*=zdist) {
                  int imax = shelltype - j; 
                  for (i=0, ydp=1.0f, xdp=pow(xdist, imax); i<=imax; i++, ydp*=ydist, xdp*=xdiv) {
                    tmpshell += wave_f[ifunc++] * xdp * ydp * zdp;
                  }
                }
                value += tmpshell * contracted_gto;
#endif
            } // end switch

            shell_counter++;
          } // end shell
        } // end atom

        // return either orbital density or orbital wavefunction amplitude
        if (density) {
          __m256 mask = _mm256_cmp_ps(value, _mm256_set1_ps(0.0f), _CMP_LT_OQ);
          __m256 sqdensity = _mm256_mul_ps(value, value);
          __m256 orbdensity = sqdensity;
          __m256 nsqdensity = _mm256_and_ps(sqdensity, mask);
          orbdensity = _mm256_sub_ps(orbdensity, nsqdensity);
          orbdensity = _mm256_sub_ps(orbdensity, nsqdensity);
          _mm256_storeu_ps(&orbitalgrid[gaddrzy + nx], orbdensity);
        } else {
          _mm256_storeu_ps(&orbitalgrid[gaddrzy + nx], value);
        }
      }
    }
  }

  // prevent x86 AVX2 clock rate limiting performance loss due to 
  // false dependence on upper vector register state for scalar or 
  // SSE instructions executing after an AVX2 instruction has written
  // an upper register.
  _mm256_zeroupper();

  return 0;
}

#endif


