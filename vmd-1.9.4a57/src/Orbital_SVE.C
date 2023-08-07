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
 *	$RCSfile: Orbital_SVE.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.9 $	$Date: 2022/04/16 03:04:51 $
 *
 ***************************************************************************/
/**
 * \file Orbital_SVE.C
 * \brief SIMD vectorized MO kernels for ARM SVE instructions.
 */

//
// Notes:  
//   Tests on the ORNL Wombat cluster reveal that several 
// compiler versions don't generate correct SVE code for 
// particular intrinsics, such as svcvt_s32_f32_x(),
// svcvt_s32_f32_z(), and likely others.
//
// Compilers found to generate bad code include:
//   armclang 21.1
//   armclang 22.0
//   llvm/clang 10.0.1
//
// Compilers generating correct code include:
//   armclang 20.3
//   llvm/clang 11.0.1
//

#if defined(VMDCPUDISPATCH) && defined(__ARM_FEATURE_SVE) 
#include <arm_sve.h>

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

#if defined(__GNUC__) 
#define __align(X)  __attribute__((aligned(X) ))
#else
#define __align(X) __declspec(align(X) )
#endif

#define MLOG2EF    -1.44269504088896f

#if 0
static void print_svfloat32_t(svfloat32_t v) {
  float tmp[svcntw()]; // vector-size aligned for SVE

  svbool_t pg = svptrue_b32();
  svst1_f32(pg, &tmp[0], v);

  printf("print_svfloat32_t: ");
  for (int i=0; i<svcntw(); i++)
    printf("%g ", tmp[i]);
  printf("\n");
}

static void print_svint32_t(svint32_t v) {
  int tmp[svcntw()]; // vector-size-aligned for SVE

  svbool_t pg = svptrue_b32();
  svst1_s32(pg, &tmp[0], v);

  printf("print_svint32_t: ");
  for (int i=0; i<svcntw(); i++)
    printf("%d ", tmp[i]);
  printf("\n");
}

static void print_svhex32_t(svint32_t v) {
  int tmp[4]; // vector-size-aligned for SVE

  svbool_t pg = svptrue_b32();
  svst1_s32(pg, &tmp[0], v);

  printf("print_svhex32_t: ");
  for (int i=0; i<svcntw(); i++)
    printf("%08x ", tmp[i]);
  printf("\n");
}
#endif


//
// John Stone, April 2022
//
// aexpfnxsve() - ARM SVE version of aexpfnx().
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

svfloat32_t aexpfnxarmsve(svfloat32_t x) {
  svbool_t pg = svptrue_b32();

  // If all x are outside of cutoff, return 0s.
  float fmv = svmaxv(pg, x);
  if (fmv < ACUTOFF) {
    return  svdup_f32(0.0f);
  }

  // if x < ACUTOFF, we return 0..
  pg = svcmpge_f32(pg, x, svdup_f32(ACUTOFF)); // Is x within cutoff?

  /*
   * Convert base:  exp(x) = 2^(N-d) where N is integer and 0 <= d < 1.
   *
   * Below we calculate n=N and x=-d, with "y" for temp storage,
   * calculate floor of x*log2(e) and subtract to get -d.
   */
  svfloat32_t mb = svmul_f32_x(pg, x, svdup_f32(MLOG2EF));
  svint32_t intflr = svcvt_s32_f32_x(pg, mb);
  svfloat32_t mbflr = svcvt_f32_s32_x(pg, intflr);
  svfloat32_t d = svsub_f32_x(pg, mbflr, mb);

  // Approximate 2^{-d}, 0 <= d < 1, by interpolation.
  // Perform Horner's method to evaluate interpolating polynomial.
  svfloat32_t y;
  y = svmad_f32_x(pg, d, svdup_f32(SCEXP4), svdup_f32(SCEXP3));
  y = svmad_f32_x(pg, y, d, svdup_f32(SCEXP2));
  y = svmad_f32_x(pg, y, d, svdup_f32(SCEXP1));
  y = svmad_f32_x(pg, y, d, svdup_f32(SCEXP0));

  // Calculate 2^N exactly by directly manipulating floating point exponent,
  // then use it to scale y for the final result.
  svint32_t flint = svsub_s32_x(pg, svdup_s32(EXPOBIAS), intflr);
  flint = svlsl_s32_x(pg, flint, svdup_u32(EXPOSHIFT));
  y = svmul_f32_z(pg, y, svreinterpret_f32_s32(flint));
  return y;
}


//
// ARM SVE implementation, uses exponential approximation
//
int evaluate_grid_sve(int numatoms,
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

  svbool_t pg = svptrue_b32();

  int nx, ny, nz;
  int vecsize = svcntw();
  float sxdelta[vecsize];
  for (nx=0; nx<vecsize; nx++) 
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
        svfloat32_t value = svdup_f32(0.0f);

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

          svfloat32_t xdelta = svld1_f32(pg, &sxdelta[0]); // aligned load
          svfloat32_t xdist  = svdup_f32(sxdist);
          xdist = svadd_f32_z(pg, xdist, xdelta);
          svfloat32_t ydist  = svdup_f32(sydist);
          svfloat32_t zdist  = svdup_f32(szdist);
          svfloat32_t xdist2 = svmul_f32_z(pg, xdist, xdist);
          svfloat32_t ydist2 = svmul_f32_z(pg, ydist, ydist);
          svfloat32_t zdist2 = svmul_f32_z(pg, zdist, zdist);
          svfloat32_t dist2  = svdup_f32(yzdist2); 
          dist2 = svadd_f32_z(pg, dist2, xdist2);
 
          // loop over the shells belonging to this atom
          // XXX this is maybe a misnomer because in split valence
          //     basis sets like 6-31G we have more than one basis
          //     function per (valence-)shell and we are actually
          //     looping over the individual contracted GTOs
          for (shell=0; shell < maxshell; shell++) {
            svfloat32_t contracted_gto = svdup_f32(0.0f);

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
              svfloat32_t expval = svmul_f32_z(pg, svdup_f32(exponent), dist2);
              // exp2f() equivalent required, use base-2 approximation
              svfloat32_t retval = aexpfnxarmsve(expval);
              contracted_gto = svmad_f32_z(pg, svdup_f32(contract_coeff), retval, contracted_gto);

              prim_counter += 2;
            }

            /* multiply with the appropriate wavefunction coefficient */
            svfloat32_t tmpshell = svdup_f32(0.0f);
            switch (shelltype) {
              // use FMADD instructions
              case S_SHELL:
                value = svmad_f32_z(pg, svdup_f32(wave_f[ifunc++]), contracted_gto, value);
                break;

              case P_SHELL:
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), xdist, tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), ydist, tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), zdist, tmpshell);
                value = svmad_f32_z(pg, tmpshell, contracted_gto, value);
                break;

              case D_SHELL:
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), xdist2, tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, xdist, ydist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), ydist2, tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, xdist, zdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, ydist, zdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), zdist2, tmpshell);
                value = svmad_f32_z(pg, tmpshell, contracted_gto, value);
                break;

              case F_SHELL:
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, xdist2, xdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, xdist2, ydist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, ydist2, xdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, ydist2, ydist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, xdist2, zdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, svmul_f32_x(pg, xdist, ydist), zdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, ydist2, zdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, zdist2, xdist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, zdist2, ydist), tmpshell);
                tmpshell = svmad_f32_x(pg, svdup_f32(wave_f[ifunc++]), svmul_f32_x(pg, zdist2, zdist), tmpshell);
                value = svmad_f32_z(pg, tmpshell, contracted_gto, value);
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
          pg = svcmplt_f32(pg, value, svdup_f32(0.0f));
          svfloat32_t sqdensity = svmul_f32_z(pg, value, value);
          svfloat32_t orbdensity = svsel_f32(pg, sqdensity, svmul_f32_z(pg, sqdensity, svdup_f32(-1.0f)));
          svst1_f32(pg, &orbitalgrid[gaddrzy + nx], orbdensity);
        } else {
          svst1_f32(pg, &orbitalgrid[gaddrzy + nx], value);
        }
      }
    }
  }

  return 0;
}

#else

int evaluate_grid_sve(int numatoms,
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
  return -1; // unimplemented by this compiler toolchain
}

#endif
