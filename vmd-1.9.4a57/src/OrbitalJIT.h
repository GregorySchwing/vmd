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
 *      $RCSfile: OrbitalJIT.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $      $Date: 2020/02/26 06:00:57 $
 *
 ***************************************************************************/
/**
 * \file OrbitalJIT.h
 * \brief Just-in-time (JIT) CUDA and OpenCL kernel generation for 
 *        computation of molecular orbitals on a uniformly spaced grid.
 */

#define ORBITAL_JIT_CUDA   0
#define ORBITAL_JIT_OPENCL 1

int orbital_jit_generate(int jitlanguage,
                         const char * srcfilename, int numatoms,
                         const float *wave_f, const float *basis_array,
                         const int *atom_basis,
                         const int *num_shells_per_atom,
                         const int *num_prim_per_shell,
                         const int *shell_types);
