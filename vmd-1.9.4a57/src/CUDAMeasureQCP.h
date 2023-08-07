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
 *      $RCSfile: CUDAMeasureQCP.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $       $Date: 2020/02/26 20:16:56 $
 *
 ***************************************************************************/
/**
 * \file CUDAMeasureQCP.h
 * \brief CUDA Quaternion Characteristic Polynomial calculation.
 *
 *  Compute RMSD values for unaligned structures without
 *  actually performing the alginment.  This is particularly useful for
 *  computing large dissimilarity matrices required for
 *  trajectory clustering analysis.
 */

int qcp_soa_gpu_ooc(wkf_threadpool_t *devpool, // VMD GPU worker thread pool
                    int nfiles, const char **trjfileset, const AtomSel *sel,
                    int first, int last, int step, float *rmsdmat);


