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
 *      $RCSfile: Benchmark.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.13 $      $Date: 2020/07/28 08:21:19 $
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

#include <stdlib.h>
#include <string.h>
#include "WKFUtils.h"
#include "WKFThreads.h"
#include "utilities.h"

/*
 * On compilers that accept the C99 'restrict' keyword, we can give
 * the compiler additional help with optimization.  Since the caller is
 * contained within the same source file, this shouldn't be necessary
 * in the current case however. 
 */
#if 0
#define RESTRICT restrict
#else
#define RESTRICT 
#endif

/*
 * If we want, we can create compiler-specific vectorization 
 * helper macros to assist with achieving peak performance, though 
 * this really shouldn't be required.
 */
#if 0
#define VECTORIZEME _Pragma("vector always")
#else
#define VECTORIZEME 
#endif


/*
 * Double precision stream bandwidth tests
 */

void dstream_init(double * RESTRICT a, double * RESTRICT b,
                  double * RESTRICT c, int N) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++) {
    a[j] = 1.0;
    b[j] = 2.0;
    c[j] = 0.0;
  }
}

void dstream_copy(double * RESTRICT a, const double * RESTRICT b, 
                 int N, double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = b[j];

  *mbsize = (2L * sizeof(double) * N) / (1024.0 * 1024.0);
}

void dstream_scale(double * RESTRICT a, const double * RESTRICT b, 
                  double scalar, int N, double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = scalar * b[j];

  *mbsize = (2L * sizeof(double) * N) / (1024.0 * 1024.0);
}

void dstream_add(double * RESTRICT a, const double * RESTRICT b, 
                const double * RESTRICT c, int N, double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = b[j] + c[j];

  *mbsize = (3L * sizeof(double) * N) / (1024.0 * 1024.0);
}

void dstream_triad(double * RESTRICT a, const double * RESTRICT b, 
                  const double * RESTRICT c, double scalar, int N, 
                  double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = b[j] + scalar * c[j];

  *mbsize = (3L * sizeof(double) * N) / (1024.0 * 1024.0);
}



/*
 * Single precision stream bandwidth tests
 */

void fstream_init(float * RESTRICT a, float * RESTRICT b,
                  float * RESTRICT c, int N) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++) {
    a[j] = 1.0f;
    b[j] = 2.0f;
    c[j] = 0.0f;
  }
}

void fstream_copy(float * RESTRICT a, const float * RESTRICT b, 
                 int N, double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = b[j];

  *mbsize = (2L * sizeof(float) * N) / (1024.0 * 1024.0);
}

void fstream_scale(float * RESTRICT a, const float * RESTRICT b, 
                   float scalar, int N, double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = scalar * b[j];

  *mbsize = (2L * sizeof(float) * N) / (1024.0 * 1024.0);
}

void fstream_add(float * RESTRICT a, const float * RESTRICT b, 
                 const float * RESTRICT c, int N, double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = b[j] + c[j];

  *mbsize = (3L * sizeof(float) * N) / (1024.0 * 1024.0);
}

void fstream_triad(float * RESTRICT a, const float * RESTRICT b, 
                  const float * RESTRICT c, float scalar, int N, 
                  double *mbsize) {
  int j;
VECTORIZEME
  for (j=0; j<N; j++)
    a[j] = b[j] + scalar * c[j];

  *mbsize = (3L * sizeof(float) * N) / (1024.0 * 1024.0);
}


/*
 * run the benchmark
 */
int stream_bench(int N, double *time, double *mbsec) {
  double *da, *db, *dc;
  float *fa, *fb, *fc;
  wkf_timerhandle timer;
  int rc = 0;

  timer = wkf_timer_create();

  /*
   * run double precision benchmarks
   */
  da = (double *) malloc(N * sizeof(double));
  db = (double *) malloc(N * sizeof(double));
  dc = (double *) malloc(N * sizeof(double));

  if ((da != NULL) && (db != NULL) && (dc != NULL)) {
    double mbsz;

    dstream_init(da, db, dc, N);

    wkf_timer_start(timer);
    dstream_copy(da, db, N, &mbsz);
    wkf_timer_stop(timer);
    time[0] = wkf_timer_time(timer);
    mbsec[0] = mbsz / time[0];

    wkf_timer_start(timer);
    dstream_scale(da, db, 2.0, N, &mbsz);
    wkf_timer_stop(timer);
    time[1] = wkf_timer_time(timer);
    mbsec[1] = mbsz / time[1];

    wkf_timer_start(timer);
    dstream_add(da, db, dc, N, &mbsz);
    wkf_timer_stop(timer);
    time[2] = wkf_timer_time(timer);
    mbsec[2] = mbsz / time[2];

    wkf_timer_start(timer);
    dstream_triad(da, db, dc, 2.0, N, &mbsz);
    wkf_timer_stop(timer);
    time[3] = wkf_timer_time(timer);
    mbsec[3] = mbsz / time[3];
  } else {
    rc = -1;
  }

  if (da)
    free(da);
  if (db)
    free(db);
  if (dc)
    free(dc);

  if (rc) {
    wkf_timer_destroy(timer);
    return rc;
  }

  /*
   * run float precision benchmarks
   */
  fa = (float *) malloc(N * sizeof(float));
  fb = (float *) malloc(N * sizeof(float));
  fc = (float *) malloc(N * sizeof(float));

  if ((fa != NULL) && (fb != NULL) && (fc != NULL)) {
    double mbsz;

    fstream_init(fa, fb, fc, N);

    wkf_timer_start(timer);
    fstream_copy(fa, fb, N, &mbsz);
    wkf_timer_stop(timer);
    time[4] = wkf_timer_time(timer);
    mbsec[4] = mbsz / time[4];

    wkf_timer_start(timer);
    fstream_scale(fa, fb, 2.0, N, &mbsz);
    wkf_timer_stop(timer);
    time[5] = wkf_timer_time(timer);
    mbsec[5] = mbsz / time[5];

    wkf_timer_start(timer);
    fstream_add(fa, fb, fc, N, &mbsz);
    wkf_timer_stop(timer);
    time[6] = wkf_timer_time(timer);
    mbsec[6] = mbsz / time[6];

    wkf_timer_start(timer);
    fstream_triad(fa, fb, fc, 2.0, N, &mbsz);
    wkf_timer_stop(timer);
    time[7] = wkf_timer_time(timer);
    mbsec[7] = mbsz / time[7];
  } else {
    rc = -1;
  }

  if (fa)
    free(fa);
  if (fb)
    free(fb);
  if (fc)
    free(fc);

  wkf_timer_destroy(timer);

  return rc;
}



void vmdbench_minmax_1fv(int sz, int reps, double &runtime, double &bwmbsec) {
  int i;
  float minf=0, maxf=0;

  // generate trivial test array
  float *fv = (float *) malloc(sz * sizeof(float));
  for (i=0; i<sz; i++) {
    fv[i] = (float) i;
  }

  wkf_timerhandle timer;
  timer = wkf_timer_create();
  wkf_timer_start(timer);
  int r;
  for (r=0; r<reps; r++)
    minmax_1fv_aligned(fv, sz, &minf, &maxf);
  wkf_timer_stop(timer);
  runtime = wkf_timer_time(timer);

//  printf("minmax_1fv_aligned: %f\n", runtime);
//  printf("  min: %f max: %f\n", minf, maxf);

  bwmbsec = (reps * sz * sizeof(float) / (1024.0 * 1024.0)) / runtime;
//  printf("  BW: %.1f MB/sec\n", bwmbsec);

  free(fv);
  wkf_timer_destroy(timer);
}


void vmdbench_minmaxmean_1fv(int sz, int reps,
                             double &runtime, double &bwmbsec) {
  int i;
  float minf=0, maxf=0, meanf=0;

  // generate trivial test array
  float *fv = (float *) malloc(sz * sizeof(float));
  for (i=0; i<sz; i++) {
    fv[i] = (float) i;
  }

  wkf_timerhandle timer;
  timer = wkf_timer_create();
  wkf_timer_start(timer);
  int r;
  for (r=0; r<reps; r++)
    minmaxmean_1fv_aligned(fv, sz, &minf, &maxf, &meanf);
  wkf_timer_stop(timer);
  runtime = wkf_timer_time(timer);

//  printf("minmaxmean_1fv_aligned: %f\n", runtime);
//  printf("  min: %f max: %f mean: %f\n", minf, maxf, meanf);

  bwmbsec = (reps * sz * sizeof(float) / (1024.0 * 1024.0)) / runtime;
//  printf("  BW: %.1f MB/sec\n", bw);

  free(fv);
  wkf_timer_destroy(timer);
}


void vmdbench_minmax_3fv(int sz, int reps, double &runtime, double &bwmbsec) {
  int i;
  float minfv[3] = { 0 }, maxfv[3] = { 0 };

  // generate trivial test array
  float *fv = (float *) malloc(3L * sz * sizeof(float));
  for (i=0; i<sz * 3L; i++) {
    fv[i] = (float) i;
  }

  wkf_timerhandle timer;
  timer = wkf_timer_create();
  wkf_timer_start(timer);
  int r;
  for (r=0; r<reps; r++)
    minmax_3fv_aligned(fv, sz, minfv, maxfv);
  wkf_timer_stop(timer);
  runtime = wkf_timer_time(timer);

//  printf("minmax_3fv_aligned: %f\n", wkf_timer_time(timer));
//  int i;
//  for (i=0; i<3; i++) {
//    minf += minfv[i];
//    maxf += maxfv[i];
//  }
//  printf("  min: %f max: %f\n", minf, maxf);

  bwmbsec = (reps * 3L * sz * sizeof(float) / (1024.0 * 1024.0)) / runtime;
//  printf("  BW: %.1f MB/sec\n", bwmbsec);

  free(fv);
  wkf_timer_destroy(timer);
}


void vmdbench_analyze_selection(int sz, int reps,
                                double &runtime, double &bwmbsec) {
  int i;
  int first=0, last=-1, selected=0;
  int *on = (int *) calloc(1, sz * sizeof(int));

  // set one atom per group of 8 rotating through all lanes
  int lane=0;
  for (i=sz/2; i<(sz-7); i+=8) {
    on[i+lane] = 1;
    lane = (lane+1) & 0x7; // swizzle through lanes
  }

  wkf_timerhandle timer;
  timer = wkf_timer_create();
  wkf_timer_start(timer);
  int r;
  for (r=0; r<reps; r++)
    analyze_selection_aligned(sz, on, &first, &last, &selected);
  wkf_timer_stop(timer);
  runtime = wkf_timer_time(timer);

//  printf("selection stats: first %d  last %d  sel %d  time: %f\n",
//         first, last, selected, runtime);
  bwmbsec = (reps * sz * sizeof(int) / (1024.0 * 1024.0)) / runtime;
//  printf("  BW: %.1f MB/sec\n", bwmbsec);

  free(on);
  wkf_timer_destroy(timer);
}






