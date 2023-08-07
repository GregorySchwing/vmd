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
 *	$RCSfile: Orbital_AVX512.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 2020/10/27 04:18:28 $
 *
 ***************************************************************************/
/**
 * \file Orbital_AVX512.C
 * \brief SIMD vectorized MO kernels for AVX-512F instructions.
 */

// Due to differences in code generation between gcc/intelc/clang/msvc, we
// don't have to check for a defined(__AVX512F__)
#if defined(VMDCPUDISPATCH) && defined(VMDUSEAVX512) 

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
static void print_mm512_ps(__m512 v) {
  __attribute__((aligned(64))) float tmp[16]; // 64-byte aligned for AVX512
  _mm512_storeu_ps(&tmp[0], v);

  printf("mm512: ");
  int i;
  for (i=0; i<16; i++) 
    printf("%g ", tmp[i]);
  printf("\n");
}
#endif


//
// John Stone, March 2017
//
// aexpfnxavx512f() - AVX-512F version of aexpfnx().
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

typedef union AVX512reg_t {
  __m512  f;  // 16x float (AVX-512F)
  __m512i i;  // 16x 32-bit int (AVX-512F)
} AVX512reg;

__m512 aexpfnxavx512f(__m512 x) {
  __mmask16 mask;
  mask = _mm512_cmpnle_ps_mask(_mm512_set1_ps(ACUTOFF), x); // Is x within cutoff?
#if 0
  // If all x are outside of cutoff, return 0s.
  if (_mm512_movemask_ps(scal.f) == 0) {
    return _mm512_set1_ps(0.0f);
  }
  // Otherwise, scal.f contains mask to be ANDed with the scale factor
#endif

  /*
   * Convert base:  exp(x) = 2^(N-d) where N is integer and 0 <= d < 1.
   *
   * Below we calculate n=N and x=-d, with "y" for temp storage,
   * calculate floor of x*log2(e) and subtract to get -d.
   */
  __align(64) AVX512reg n;
  __m512 mb = _mm512_mul_ps(x, _mm512_set1_ps(MLOG2EF));
  n.i = _mm512_cvttps_epi32(mb);
  __m512 mbflr = _mm512_cvtepi32_ps(n.i);
  __m512 d = _mm512_sub_ps(mbflr, mb);

  // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
  // Perform Horner's method to evaluate interpolating polynomial.
  __m512 y;
  y = _mm512_fmadd_ps(d, _mm512_set1_ps(SCEXP4), _mm512_set1_ps(SCEXP3));
  y = _mm512_fmadd_ps(y, d, _mm512_set1_ps(SCEXP2));
  y = _mm512_fmadd_ps(y, d, _mm512_set1_ps(SCEXP1));
  y = _mm512_fmadd_ps(y, d, _mm512_set1_ps(SCEXP0));

  // Calculate 2^N exactly by directly manipulating floating point exponent,
  // then use it to scale y for the final result.
  n.i = _mm512_sub_epi32(_mm512_set1_epi32(EXPOBIAS), n.i);
  n.i = _mm512_slli_epi32(n.i, EXPOSHIFT);
  n.f = _mm512_mask_mul_ps(n.f, mask, _mm512_set1_ps(0.0f), n.f);
  y = _mm512_mul_ps(y, n.f);
  return y;
}


//
// AVX-512F implementation for Xeons that don't have special fctn units
//
int evaluate_grid_avx512f(int numatoms,
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
  __attribute__((aligned(64))) float sxdelta[16]; // 64-byte aligned for AVX512
  for (nx=0; nx<16; nx++) 
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
      for (nx=0; nx<numvoxels[0]; nx+=16) {
        grid_x = origin[0] + nx * voxelsize;

        // calculate the value of the wavefunction of the
        // selected orbital at the current grid point
        int at;
        int prim, shell;

        // initialize value of orbital at gridpoint
        __m512 value = _mm512_set1_ps(0.0f);

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

          __m512 xdelta = _mm512_load_ps(&sxdelta[0]); // aligned load
          __m512 xdist  = _mm512_set1_ps(sxdist);
          xdist = _mm512_add_ps(xdist, xdelta);
          __m512 ydist  = _mm512_set1_ps(sydist);
          __m512 zdist  = _mm512_set1_ps(szdist);
          __m512 xdist2 = _mm512_mul_ps(xdist, xdist);
          __m512 ydist2 = _mm512_mul_ps(ydist, ydist);
          __m512 zdist2 = _mm512_mul_ps(zdist, zdist);
          __m512 dist2  = _mm512_set1_ps(yzdist2); 
          dist2 = _mm512_add_ps(dist2, xdist2);
 
          // loop over the shells belonging to this atom
          // XXX this is maybe a misnomer because in split valence
          //     basis sets like 6-31G we have more than one basis
          //     function per (valence-)shell and we are actually
          //     looping over the individual contracted GTOs
          for (shell=0; shell < maxshell; shell++) {
            __m512 contracted_gto = _mm512_set1_ps(0.0f);

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
              __m512 expval = _mm512_mul_ps(_mm512_set1_ps(exponent), dist2);
              // exp2f() equivalent required, use base-2 approximation
              __m512 retval = aexpfnxavx512f(expval);
              contracted_gto = _mm512_fmadd_ps(_mm512_set1_ps(contract_coeff), retval, contracted_gto);

              prim_counter += 2;
            }

            /* multiply with the appropriate wavefunction coefficient */
            __m512 tmpshell = _mm512_set1_ps(0.0f);
            switch (shelltype) {
              // use FMADD instructions
              case S_SHELL:
                value = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), contracted_gto, value);
                break;

              case P_SHELL:
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), xdist, tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), ydist, tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), zdist, tmpshell);
                value = _mm512_fmadd_ps(tmpshell, contracted_gto, value);
                break;

              case D_SHELL:
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), xdist2, tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(xdist, ydist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), ydist2, tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(xdist, zdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(ydist, zdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), zdist2, tmpshell);
                value = _mm512_fmadd_ps(tmpshell, contracted_gto, value);
                break;

              case F_SHELL:
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(xdist2, xdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(xdist2, ydist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(ydist2, xdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(ydist2, ydist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(xdist2, zdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(_mm512_mul_ps(xdist, ydist), zdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(ydist2, zdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(zdist2, xdist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(zdist2, ydist), tmpshell);
                tmpshell = _mm512_fmadd_ps(_mm512_set1_ps(wave_f[ifunc++]), _mm512_mul_ps(zdist2, zdist), tmpshell);
                value = _mm512_fmadd_ps(tmpshell, contracted_gto, value);
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
          __mmask16 mask = _mm512_cmplt_ps_mask(value, _mm512_set1_ps(0.0f));
          __m512 sqdensity = _mm512_mul_ps(value, value);
          __m512 orbdensity = _mm512_mask_mul_ps(sqdensity, mask, sqdensity,
                                                 _mm512_set1_ps(-1.0f));
          _mm512_storeu_ps(&orbitalgrid[gaddrzy + nx], orbdensity);
        } else {
          _mm512_storeu_ps(&orbitalgrid[gaddrzy + nx], value);
        }
      }
    }
  }

  // prevent x86 AVX-512 clock rate limiting performance loss due to 
  // false dependence on upper vector register state for scalar or 
  // SSE instructions executing after an AVX-512 instruction has written
  // an upper register.
  _mm256_zeroupper();

  return 0;
}

#endif


