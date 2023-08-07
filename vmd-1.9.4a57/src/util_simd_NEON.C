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
 *      $RCSfile: util_simd_NEON.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $        $Date: 2022/03/31 16:20:47 $
 *
 ***************************************************************************/
/**
 * \file util_simd_NEON.C
 * \brief SIMD vectorized kernels for ARM64 NEON instructions.
 */

// Due to differences in code generation between gcc/intelc/clang/msvc, we
// don't have to check for a defined(__NEON__)
#if defined(VMDCPUDISPATCH) && defined(VMDUSENEON)
#include <arm_neon.h>

#include "WKFThreads.h" // CPU capability flags
// #include <string.h>
// #include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>


//
// Helper routine for use when coping with unaligned
// buffers returned by malloc() on many GNU systems:
//   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24261
//   http://www.sourceware.org/bugzilla/show_bug.cgi?id=206
//
// XXX until all compilers support uintptr_t, we have to do
//     dangerous and ugly things with pointer casting here...
//
#if defined(_WIN64) /* sizeof(size_t) == sizeof(void*) */
#define myintptrtype size_t
#elif 1 /* sizeof(unsigned long) == sizeof(void*) */
#define myintptrtype unsigned long
#else /* C99 */
#define myintptrtype uintptr_t
#endif


//
// Aligment test routines for vector instructions
//
static int test_alignment_Nbyte_powertwo(const void *ptr, unsigned int alignsz) {
  unsigned int alignmask = alignsz - 1;
  return (((myintptrtype) ptr) == (((myintptrtype) ptr) & (~alignmask)));
}


// helper routine to perform a min among all 4 elements of an __m128
static float fmin_f32x4(float32x4_t min4) {
  float *f1 = (float *) &min4;
  float min1 = f1[0];
  if (f1[1] < min1) min1 = f1[1];
  if (f1[2] < min1) min1 = f1[2];
  if (f1[3] < min1) min1 = f1[3];
  return min1;
}

static float fmax_f32x4(float32x4_t max4) {
  float *f1 = (float *) &max4;
  float max1 = f1[0];
  if (f1[1] > max1) max1 = f1[1];
  if (f1[2] > max1) max1 = f1[2];
  if (f1[3] > max1) max1 = f1[3];
  return max1;
}

static double hadd_f64x2(float64x2_t sum2) {
  double *d = (double *) &sum2;
  return d[0] + d[1]; 
}


// Compute min/max/mean values for a an arbitrary array of floats
void minmaxmean_1fv_aligned_neon(const float *f, ptrdiff_t n,
                                 float *fmin, float *fmax, float *fmean) {
  if (n < 1) {
    *fmin = 0.0f;
    *fmax = 0.0f;
    *fmean = 0.0f;
    return;
  }

  float32x4_t minv = vdupq_n_f32(f[0]);
  float32x4_t maxv = minv;
  float64x2_t meanv = vdupq_n_f64(0.0);

  for (ptrdiff_t i=0; i<n; i+=4) {
    float32x4_t tmp = vld1q_f32(&f[i]);
    minv = vminq_f32(minv, tmp);
    maxv = vmaxq_f32(maxv, tmp);
    meanv = vaddq_f64(meanv, vcvt_f64_f32(vget_high_f32(tmp)));
    meanv = vaddq_f64(meanv, vcvt_f64_f32(vget_high_f32(tmp)));
  }

  *fmin = fmin_f32x4(minv);
  *fmax = fmax_f32x4(maxv);
  *fmean = hadd_f64x2(meanv) / double(n);
}


// Compute min/max values for an arbitrary array of floats
// Compute min/max values for a 16-byte-aligned array of floats
void minmax_1fv_aligned_neon(const float *f, ptrdiff_t n, float *fmin, float *fmax) {
  if (n < 1)
    return;

  ptrdiff_t i=0;
  float min1 = f[0];
  float max1 = f[0];

  // roll up to the first 16-byte-aligned array index
  for (i=0; ((i<n) && !test_alignment_Nbyte_powertwo(&f[i], 16)); i++) {
    if (f[i] < min1) min1 = f[i];
    if (f[i] > max1) max1 = f[i];
  }

  // NEON vectorized min/max loop
  float32x4_t min4 = vdupq_n_f32(min1);
  float32x4_t max4 = vdupq_n_f32(max1);

  // do groups of 32 elements
  for (; i<(n-31); i+=32) {
    float32x4_t f4;
    f4 = vld1q_f32(&f[i   ]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
    f4 = vld1q_f32(&f[i+ 4]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
    f4 = vld1q_f32(&f[i+ 8]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
    f4 = vld1q_f32(&f[i+12]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);

    f4 = vld1q_f32(&f[i+16]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
    f4 = vld1q_f32(&f[i+20]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
    f4 = vld1q_f32(&f[i+24]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
    f4 = vld1q_f32(&f[i+28]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
  }

  // do groups of 4 elements
  for (; i<(n-3); i+=4) {
    float32x4_t f4 = vld1q_f32(&f[i]); // assume 16-byte aligned array!
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
  }

  // finish last elements off
  for (; i<n; i++) {
    float32x4_t f4 = vdupq_n_f32(f[i]);
    min4 = vminq_f32(min4, f4);
    max4 = vmaxq_f32(max4, f4);
  }

  // compute min/max among the final 4-element vectors by shuffling
  // and and reducing the elements within the vectors
  *fmin = fmin_f32x4(min4);
  *fmax = fmax_f32x4(max4);
}


// Compute min/max values for an arbitrary array of float3s
// input value n3 is the number of 3-element vectors to process
void minmax_3fv_aligned_neon(const float *f, const ptrdiff_t n3,
                             float *fmin, float *fmax) {
  float minx, maxx, miny, maxy, minz, maxz;
  const ptrdiff_t end = n3*3L;

  if (n3 < 1)
    return;

  ptrdiff_t i=0;
  minx=maxx=f[i  ];
  miny=maxy=f[i+1];
  minz=maxz=f[i+2];

  // Since we may not be on a 16-byte boundary when we start, we roll
  // through the first few items with plain C until we get to one.
  for (; i<end; i+=3L) {
    // exit if/when we reach a 16-byte boundary for both arrays
    if (test_alignment_Nbyte_powertwo(&f[i], 16)) {
      break;
    }

    float tmpx = f[i  ];
    if (tmpx < minx) minx = tmpx;
    if (tmpx > maxx) maxx = tmpx;

    float tmpy = f[i+1];
    if (tmpy < miny) miny = tmpy;
    if (tmpy > maxy) maxy = tmpy;

    float tmpz = f[i+2];
    if (tmpz < minz) minz = tmpz;
    if (tmpz > maxz) maxz = tmpz;
  }

  // initialize min/max values
  float32x4_t xmin4 = vdupq_n_f32(minx);
  float32x4_t xmax4 = vdupq_n_f32(maxx);
  float32x4_t ymin4 = vdupq_n_f32(miny);
  float32x4_t ymax4 = vdupq_n_f32(maxy);
  float32x4_t zmin4 = vdupq_n_f32(minz);
  float32x4_t zmax4 = vdupq_n_f32(maxz);

  for (; i<(end-11); i+=12) {
    // aligned load of four consecutive 3-element vectors into
    // three 4-element vectors with de-interleaving...
    float32x4x3_t soa = vld3q_f32(&f[i]);

    // compute mins and maxes
    xmin4 = vminq_f32(xmin4, soa.val[0]);
    xmax4 = vmaxq_f32(xmax4, soa.val[0]);
    ymin4 = vminq_f32(ymin4, soa.val[1]);
    ymax4 = vmaxq_f32(ymax4, soa.val[1]);
    zmin4 = vminq_f32(zmin4, soa.val[2]);
    zmax4 = vmaxq_f32(zmax4, soa.val[2]);
  }

  minx = fmin_f32x4(xmin4);
  miny = fmin_f32x4(ymin4);
  minz = fmin_f32x4(zmin4);

  maxx = fmax_f32x4(xmax4);
  maxy = fmax_f32x4(ymax4);
  maxz = fmax_f32x4(zmax4);

  // regular C code...
  for (; i<end; i+=3) {
    float tmpx = f[i  ];
    if (tmpx < minx) minx = tmpx;
    if (tmpx > maxx) maxx = tmpx;

    float tmpy = f[i+1];
    if (tmpy < miny) miny = tmpy;
    if (tmpy > maxy) maxy = tmpy;

    float tmpz = f[i+2];
    if (tmpz < minz) minz = tmpz;
    if (tmpz > maxz) maxz = tmpz;
  }

  fmin[0] = minx;
  fmax[0] = maxx;
  fmin[1] = miny;
  fmax[1] = maxy;
  fmin[2] = minz;
  fmax[2] = maxz;
}

#endif // CPUDISPATCH+NEON

