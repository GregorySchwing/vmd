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
 *      $RCSfile: util_simd_SVE.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.6 $        $Date: 2022/04/09 17:58:39 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * ARM SVE vector helper routines and vectorized loops for use via 
 * runtime CPU dispatch.
 *
 * Arm C Language Extensions for SVE documentation used to develop 
 * the kernels in this source file (document Arm_100987_0000_06_en,
 * version 00bet6, Copyright dates 2015-2020):
 *   https://developer.arm.com/docs/100987/latest
 *   https://developer.arm.com/documentation/100987/0000
 *
 * Earlier SVE docs:
 *   https://developer.arm.com/solutions/hpc/resources/hpc-white-papers/a-sneak-peek-into-sve-and-vla-programming 
 *
 ***************************************************************************/

#if defined(VMDCPUDISPATCH) && defined(__ARM_FEATURE_SVE)
#include <arm_sve.h>

#include "WKFThreads.h" // CPU capability flags
// #include <string.h>
// #include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

int arm_sve_vecsize_32bits(void) {
  return svcntw();
}

int arm_sve_vecsize_64bits(void) {
  return svcntd();
}


// Compute min/max/mean values for a an arbitrary array of floats
void minmaxmean_1fv_aligned_sve(const float *f, ptrdiff_t n,
                                float *fmin, float *fmax, float *fmean) {
  if (n < 1) {
    *fmin = 0.0f;
    *fmax = 0.0f;
    *fmean = 0.0f;
    return;
  }

  svbool_t pg = svptrue_b32();
  svfloat32_t minv = svdup_f32(f[0]);
  svfloat32_t maxv = minv;
  svfloat64_t meanv = svdup_f64(0.0);

  for (ptrdiff_t i=0; i<n; i+=svcntw()) {
    pg = svwhilelt_b32(i, n);
    svfloat32_t tmp = svld1(pg, (float32_t *) &f[i]);

    minv = svmin_m(pg, minv, tmp);
    maxv = svmax_m(pg, maxv, tmp);
    meanv = svadd_z(pg, meanv, svcvt_f64_z(pg, tmp));
  }

  pg = svptrue_b32();
  *fmin = svminv(pg, minv);
  *fmax = svmaxv(pg, maxv);
  *fmean = float(svaddv(pg, meanv) / n);
}


// Compute min/max values for an arbitrary array of floats
void minmax_1fv_aligned_sve(const float *f, ptrdiff_t n,
                            float *fmin, float *fmax) {
  if (n < 1)
    return;

  svbool_t pg = svptrue_b32();
  svfloat32_t minv = svdup_f32(f[0]);
  svfloat32_t maxv = minv;
  for (ptrdiff_t i=0; i<n; i+=svcntw()) {
    pg = svwhilelt_b32(i, n);
    svfloat32_t tmp = svld1(pg, (float32_t *) &f[i]);
    minv = svmin_m(pg, minv, tmp);
    maxv = svmax_m(pg, maxv, tmp);
  }

  pg = svptrue_b32();
  *fmin = svminv(pg, minv);
  *fmax = svmaxv(pg, maxv);
}


// Compute min/max values for an arbitrary array of float3s
// input value n3 is the number of 3-element vectors to process
void minmax_3fv_aligned_sve(const float *f, const ptrdiff_t n3,
                            float *fmin, float *fmax) {
  if (n3 < 1)
    return;

  svbool_t pg = svptrue_b32();
  svfloat32x3_t minv = svcreate3(svdup_f32(f[0]),
                                 svdup_f32(f[0]),
                                 svdup_f32(f[0]));

  svfloat32x3_t maxv = minv;
  int vlen = svcntw();
  int vlen3 = vlen*3;
  ptrdiff_t cnt, i;
  for (cnt=0,i=0; cnt<n3; cnt+=vlen,i+=vlen3) {
    pg = svwhilelt_b32(cnt, n3);
    svfloat32x3_t tmp = svld3(pg, (float32_t *) &f[i]);
    svset3(minv, 0, svmin_m(pg, svget3(minv, 0), svget3(tmp, 0)));
    svset3(maxv, 0, svmax_m(pg, svget3(maxv, 0), svget3(tmp, 0)));
    svset3(minv, 1, svmin_m(pg, svget3(minv, 1), svget3(tmp, 1)));
    svset3(maxv, 1, svmax_m(pg, svget3(maxv, 1), svget3(tmp, 1)));
    svset3(minv, 2, svmin_m(pg, svget3(minv, 2), svget3(tmp, 2)));
    svset3(maxv, 2, svmax_m(pg, svget3(maxv, 2), svget3(tmp, 2)));
  }

  pg = svptrue_b32();
  fmin[0] = svminv(pg, svget3(minv, 0));
  fmax[0] = svmaxv(pg, svget3(maxv, 0));
  fmin[1] = svminv(pg, svget3(minv, 1));
  fmax[1] = svmaxv(pg, svget3(maxv, 1));
  fmin[2] = svminv(pg, svget3(minv, 2));
  fmax[2] = svmaxv(pg, svget3(maxv, 2));
}


#else // CPUDISPATCH+SVE

int arm_sve_vecsize_32bits(void) {
  return -1;
}

int arm_sve_vecsize_64bits(void) {
  return -1;
}

#endif

