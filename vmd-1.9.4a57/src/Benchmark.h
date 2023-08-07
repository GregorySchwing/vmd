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
 *      $RCSfile: Benchmark.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $      $Date: 2020/07/24 19:26:42 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * Various CPU/memory subsystem benchmarking routines.
 * The peak performance numbers achieved within a VMD build can be 
 * used to determine how well the VMD build was optimized, the 
 * performance of the host CPU/memory systems, SMP scaling efficiency, etc.
 *
 * The streaming memory bandwidth tests are an alternative implementation 
 * of McCalpin's STREAM benchmark.
 *
 ***************************************************************************/


/// built-in minor variation of the STREAM benchmark
int stream_bench(int N, double *time, double *mbsec);

/// benchmarking routines for atom selections and SIMD-accelerated
/// array/vector/volume analysis routines
void vmdbench_minmax_1fv(int sz, int reps, double &runtime, double &bwmbsec);

void vmdbench_minmaxmean_1fv(int sz, int reps,
                             double &runtime, double &bwmbsec);

void vmdbench_minmax_3fv(int sz, int reps, double &runtime, double &bwmbsec);

void vmdbench_analyze_selection(int sz, int reps,
                                double &runtime, double &bwmbsec);


