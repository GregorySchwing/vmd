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
 *	$RCSfile: CUDAOrbital.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 2020/02/26 20:38:31 $
 *
 ***************************************************************************/
/**
 * \file CUDAOrbital.h
 * \brief CUDA Orbital kernels for computing molecular orbital amplitudes
 *        on uniform grids.
 *
 * This work is described in the following papers:
 *
 *  "NAMD goes quantum: an integrative suite for hybrid simulations"
 *  Marcelo C. R. Melo, Rafael C. Bernardi, Till Rudack, Maximilian Scheurer, 
 *  Christoph Riplinger, James C. Phillips, Julio D. C. Maia, Gerd B. Rocha, 
 *  João V. Ribeiro, John E. Stone, Frank Neese, Klaus Schulten, 
 *  and Zaida Luthey-Schulten.
 *  Nature Methods, 15: 351-354, 2018.
 *  https://www.nature.com/articles/nmeth.4638
 *
 *  "Early Experiences Porting the NAMD and VMD Molecular Simulation and 
 *   Analysis Software to GPU-Accelerated OpenPOWER Platforms"
 *  John E. Stone, Antti-Pekka Hynninen, James C. Phillips, and Klaus Schulten.
 *  International Workshop on OpenPOWER for HPC (IWOPH'16), 
 *  LNCS 9945, pp. 188-206, 2016.
 *  http://dx.doi.org/10.1007/978-3-319-46079-6_14                          
 *  
 *  "Evaluation of Emerging Energy-Efficient Heterogeneous Computing Platforms 
 *  for Biomolecular and Cellular Simulation Workloads"
 *  John E. Stone, Michael J. Hallock, James C. Phillips, Joseph R. Peterson, 
 *  Zaida Luthey-Schulten, and Klaus Schulten.
 *  25th International Heterogeneity in Computing Workshop, 
 *  2016 IEEE International Parallel and Distributed Processing Symposium 
 *  Workshops (IPDPSW), pp. 89-100, 2016.
 *  http://dx.doi.org/10.1109/IPDPSW.2016.130
 *
 *  "High Performance Computation and Interactive Display of Molecular Orbitals  *   on GPUs and Multi-core CPUs"
 *  John E. Stone, Jan Saam, David J. Hardy, Kirby L. Vandivort, 
 *  Wen-mei W. Hwu, and Klaus Schulten.
 *  In Proceedings of the 2nd Workshop on General-Purpose Processing 
 *  on Graphics Processing Units, ACM International Conference Proceeding 
 *  Series, volume 383, pp. 9-18, 2009.
 *  http://doi.acm.org/10.1145/1513895.1513897
 */

#ifndef CUDAORBITAL_H
#define CUDAORBITAL_H

int vmd_cuda_evaluate_orbital_grid(wkf_threadpool_t *devpool,
                       int numatoms,
                       const float *wave_f, int num_wave_f,
                       const float *basis_array, int num_basis,
                       const float *atompos,
                       const int *atom_basis,
                       const int *num_shells_per_atom,
                       const int *num_prim_per_shell,
                       const int *shell_types,
                       int num_shells,
                       const int *numvoxels,
                       float voxelsize,
                       const float *origin,
                       int density,
                       float *orbitalgrid);

#endif
