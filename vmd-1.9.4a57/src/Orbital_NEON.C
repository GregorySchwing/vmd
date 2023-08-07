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
 *	$RCSfile: Orbital_NEON.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 2022/03/31 05:55:36 $
 *
 ***************************************************************************/
/**
 * \file Orbital_NEON.C
 * \brief SIMD vectorized MO kernels for NEON instructions.
 */

// Due to differences in code generation between gcc/intelc/clang/msvc, we
// don't have to check for a defined(__NEON__)
#if defined(VMDCPUDISPATCH) && defined(VMDUSENEON) 
#include <arm_neon.h>

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

// #if defined(__GNUC__) 
#define __align(X)  __attribute__((aligned(X) ))
// #endif

#define MLOG2EF    -1.44269504088896f

#if 0
static void print_float32x4_t(float32x4_t v) {
  __attribute__((aligned(16))) float tmp[4]; // 16-byte aligned for NEON
  vst1q_f32(&tmp[0], v);

  printf("print_float32x4_t: ");
  int i;
  for (i=0; i<4; i++)
    printf("%g ", tmp[i]);
  printf("\n");
}

static void print_int32x4_t(int32x4_t v) {
  __attribute__((aligned(16))) int tmp[4]; // 16-byte aligned for NEON
  vst1q_s32(&tmp[0], v);

  printf("print_int32x4_t: ");
  int i;
  for (i=0; i<4; i++)
    printf("%d ", tmp[i]);
  printf("\n");
}

static void print_hex32x4_t(int32x4_t v) {
  __attribute__((aligned(16))) int tmp[4]; // 16-byte aligned for NEON
  vst1q_s32(&tmp[0], v);

  printf("print_hex32x4_t: ");
  int i;
  for (i=0; i<4; i++)
    printf("%08x ", tmp[i]);
  printf("\n");
}
#endif


//
// John Stone, February 2021
//
// aexpfnxneon() - NEON version of aexpfnx().
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

typedef union NEONreg_t {
  float32x4_t f;  // 4x float (NEON)
  int32x4_t   i;  // 4x 32-bit int (NEON)
} NEONreg;

float32x4_t aexpfnxneon(float32x4_t x) {
  __align(16) NEONreg scal;

#if 1
  // NEON seems to lack a convenient way to test if any lane was true, so
  // we use a different approach, and perform a horizontal maximum, comparing
  // the result of that against the cutoff to see if all four values are 
  // below the the cutoff, and to early exit returning zeros in that case.
  // If all x are outside of cutoff, return 0s.  There may be a better 
  // early exit scheme here if we dig and find some more useful NEON 
  // instructions.  This block of code is reverse-ordered from the x86 
  // variants as a result of the different scheme used here.
  float32x2_t tmp;
  tmp = vpmax_f32(vget_low_f32(x), vget_high_f32(x));
  tmp = vpmax_f32(tmp, tmp);
  float vmax = vget_lane_f32(tmp, 0);
  if (vmax < ACUTOFF) {
    return vdupq_n_f32(0.0f);
  }
#endif
  // Otherwise, scal.f contains mask to be ANDed with the scale factor
  scal.f = vcvtq_f32_u32(vcgeq_f32(x, vdupq_n_f32(ACUTOFF)));  // Is x within cutoff?

  /*
   * Convert base:  exp(x) = 2^(N-d) where N is integer and 0 <= d < 1.
   *
   * Below we calculate n=N and x=-d, with "y" for temp storage,
   * calculate floor of x*log2(e) and subtract to get -d.
   */
  __align(16) NEONreg n;
  float32x4_t mb = vmulq_f32(x, vdupq_n_f32(MLOG2EF));
  n.i = vcvtq_s32_f32(mb);
  float32x4_t mbflr = vcvtq_f32_s32(n.i);
  float32x4_t d = vsubq_f32(mbflr, mb);

  // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
  // Perform Horner's method to evaluate interpolating polynomial.
  float32x4_t y;
#if __ARM_FEATURE_FMA
  y = vfmaq_f32(vdupq_n_f32(SCEXP3), vdupq_n_f32(SCEXP4), d);
  y = vfmaq_f32(vdupq_n_f32(SCEXP2), d, y);
  y = vfmaq_f32(vdupq_n_f32(SCEXP1), d, y);
  y = vfmaq_f32(vdupq_n_f32(SCEXP0), d, y);
#else
  y = vmulq_f32(d, vdupq_n_f32(SCEXP4));    /* for x^4 term */
  y = vaddq_f32(y, vdupq_n_f32(SCEXP3));    /* for x^3 term */
  y = vmulq_f32(y, d);
  y = vaddq_f32(y, vdupq_n_f32(SCEXP2));    /* for x^2 term */
  y = vmulq_f32(y, d);
  y = vaddq_f32(y, vdupq_n_f32(SCEXP1));    /* for x^1 term */
  y = vmulq_f32(y, d);
  y = vaddq_f32(y, vdupq_n_f32(SCEXP0));    /* for x^0 term */
#endif

  // Calculate 2^N exactly by directly manipulating floating point exponent,
  // then use it to scale y for the final result.
  n.i = vsubq_s32(vdupq_n_s32(EXPOBIAS), n.i);
  n.i = vshlq_s32(n.i, vdupq_n_s32(EXPOSHIFT));
  scal.i = vandq_s32(scal.i, n.i);
  y = vmulq_f32(y, scal.f);

  return y;
}


//
// NEON implementation for Xeons that don't have special fctn units
//
int evaluate_grid_neon(int numatoms,
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
  __attribute__((aligned(16))) float sxdelta[4]; // 16-byte aligned for NEON
  for (nx=0; nx<4; nx++) 
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
      for (nx=0; nx<numvoxels[0]; nx+=4) {
        grid_x = origin[0] + nx * voxelsize;

        // calculate the value of the wavefunction of the
        // selected orbital at the current grid point
        int at;
        int prim, shell;

        // initialize value of orbital at gridpoint
        float32x4_t value = vdupq_n_f32(0.0f);

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

          float32x4_t xdelta = vld1q_f32(&sxdelta[0]); // aligned load
          float32x4_t xdist  = vdupq_n_f32(sxdist);
          xdist = vaddq_f32(xdist, xdelta);
          float32x4_t ydist  = vdupq_n_f32(sydist);
          float32x4_t zdist  = vdupq_n_f32(szdist);
          float32x4_t xdist2 = vmulq_f32(xdist, xdist);
          float32x4_t ydist2 = vmulq_f32(ydist, ydist);
          float32x4_t zdist2 = vmulq_f32(zdist, zdist);
          float32x4_t dist2  = vdupq_n_f32(yzdist2); 
          dist2 = vaddq_f32(dist2, xdist2);
 
          // loop over the shells belonging to this atom
          // XXX this is maybe a misnomer because in split valence
          //     basis sets like 6-31G we have more than one basis
          //     function per (valence-)shell and we are actually
          //     looping over the individual contracted GTOs
          for (shell=0; shell < maxshell; shell++) {
            float32x4_t contracted_gto = vdupq_n_f32(0.0f);

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
              float32x4_t expval = vmulq_f32(vdupq_n_f32(exponent), dist2);
              // exp2f() equivalent required, use base-2 approximation
              float32x4_t retval = aexpfnxneon(expval);
              contracted_gto = vfmaq_f32(contracted_gto, retval, vdupq_n_f32(contract_coeff));

              prim_counter += 2;
            }

            /* multiply with the appropriate wavefunction coefficient */
            float32x4_t tmpshell = vdupq_n_f32(0.0f);
            switch (shelltype) {
              // use FMADD instructions
              case S_SHELL:
                value = vfmaq_f32(value, contracted_gto, vdupq_n_f32(wave_f[ifunc++]));
                break;

              case P_SHELL:
                tmpshell = vfmaq_f32(tmpshell, xdist, vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, ydist, vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, zdist, vdupq_n_f32(wave_f[ifunc++]));
                value = vfmaq_f32(value, contracted_gto, tmpshell);
                break;

              case D_SHELL:
                tmpshell = vfmaq_f32(tmpshell, xdist2, vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(xdist, ydist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, ydist2, vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(xdist, zdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(ydist, zdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, zdist2, vdupq_n_f32(wave_f[ifunc++]));
                value = vfmaq_f32(value, contracted_gto, tmpshell);
                break;

              case F_SHELL:
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(xdist2, xdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(xdist2, ydist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(ydist2, xdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(ydist2, ydist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(xdist2, zdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(vmulq_f32(xdist, ydist), zdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(ydist2, zdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(zdist2, xdist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(zdist2, ydist), vdupq_n_f32(wave_f[ifunc++]));
                tmpshell = vfmaq_f32(tmpshell, vmulq_f32(zdist2, zdist), vdupq_n_f32(wave_f[ifunc++]));
                value = vfmaq_f32(value, contracted_gto, tmpshell);
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
          float32x4_t mask = vcvtq_f32_u32(vcltq_f32(value, vdupq_n_f32(0.0f)));
          float32x4_t sqdensity = vmulq_f32(value, value);
          float32x4_t orbdensity = sqdensity;
          float32x4_t nsqdensity = vmulq_f32(sqdensity, mask);
          orbdensity = vsubq_f32(orbdensity, nsqdensity);
          orbdensity = vsubq_f32(orbdensity, nsqdensity);
          vst1q_f32(&orbitalgrid[gaddrzy + nx], orbdensity);
        } else {
          vst1q_f32(&orbitalgrid[gaddrzy + nx], value);
        }
      }
    }
  }

  return 0;
}

#endif


