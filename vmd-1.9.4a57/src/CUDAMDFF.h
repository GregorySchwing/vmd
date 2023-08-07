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
 *      $RCSfile: CUDAMDFF.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.6 $      $Date: 2020/02/26 20:18:17 $
 *
 ***************************************************************************/
/**
 * \file CUDAMDFF.h
 * \brief CUDA kernels for computing simulated cryo-EM density maps,
 *        cross correlations, and other MDFF-related calculations.
 *
 * "GPU-Accelerated Analysis and Visualization of Large Structures
 *  Solved by Molecular Dynamics Flexible Fitting"
 *  John E. Stone, Ryan McGreevy, Barry Isralewitz, and Klaus Schulten.
 *  Faraday Discussions, 169:265-283, 2014.
 *  http://dx.doi.org/10.1039/C4FD00005F
 */

//
// Compute a simulated density map for a given atom selection and
// reference density map.
//
int vmd_cuda_calc_density(const AtomSel *sel, MoleculeList *mlist,  
                          int quality, float radscale, float gridspacing,
                          VolumetricData ** synthvol,
                          const VolumetricData * refmap,
                          VolumetricData ** diffvol,
                          int verbose);

//
// Fast single-pass algorithm for computing a synthetic density map
// for a given atom selection and reference map, and compute the 
// cross correlation, optionally saving a spatially localized
// cross correlation map, difference map, and simulated density map.
//
int vmd_cuda_compare_sel_refmap(const AtomSel *sel, MoleculeList *mlist,
                                int quality, float radscale, float gridspacing,
                                const VolumetricData * refmap,
                                VolumetricData **synthvol,
                                VolumetricData **diffvol,
                                VolumetricData **spatialccvol,
                                float *CC, float ccthreshdensity,
                                int verbose);

