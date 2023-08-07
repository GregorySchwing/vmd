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
 *      $RCSfile: CUDAMeasureQCP.cu,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.37 $      $Date: 2021/12/20 23:49:52 $
 *
 ***************************************************************************/
/**
 * \file CUDAMeasureQCP.cu
 * \brief CUDA Quaternion Characteristic Polynomial calculation.
 *
 *  Compute RMSD values for unaligned structures without
 *  actually performing the alginment.  This is particularly useful for
 *  computing large dissimilarity matrices required for
 *  trajectory clustering analysis.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>

#include "molfile_plugin.h"

#include "Inform.h"
#include "utilities.h"
#include "WKFThreads.h"
#include "WKFUtils.h"
#include "CUDAKernels.h" 
#include "Measure.h"

//
// Restrict macro to make it easy to do perf tuning tests
//
#if 1
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

//#define VMDUSECUDAGDS 1
#if defined(VMDUSECUDAGDS)
#include </usr/local/gds-beta-0.7.1/lib/cufile.h>        // GPU-Direct Storage

// direct calls to JS plugin for devel/testing until the plugin manager
// and headers incorporate the new out-of-core GPU-direct I/O hooks.
#define VMDJSPLUGININCLUDESRC 1
#include "/home/johns/plugins/molfile_plugin/src/jsplugin.c"

#if 1
#define CUERR { cudaError_t err; \
  if ((err = cudaGetLastError()) != cudaSuccess) { \
  printf("CUDA error: %s, %s line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
  printf("Thread aborting...\n"); \
  return NULL; }}
#else
#define CUERR
#endif

typedef struct {
  int nfiles;
  const char **trjfileset;
  jshandle **jshandles;
  CUfileHandle_t *cfh;
  int devcount;
  int natoms;
  const AtomSel *sel;
  int first;
  int last;
  int step;
  float *rmsdmat;
} gpuqcprmsdmatoocthreadparms;
#else
typedef int fio_fd;
#endif

typedef struct {
  const AtomSel *sel;
  int first;
  int last;
  int step;
  float *rmsdmat;
  int padcnt;
  int framecrdsz;
  float *crds;
} gpuqcprmsdthreadparms;


//
// Device global variable for block count used in single-pass
// device-wide parallel reductions.
//
__device__ unsigned int glob_block_count = 0;


//
// Warp-wide sum reduction
//
template <typename T>
__inline__ __device__ T warp_sum_reduction(T v) {
  for (int offset = warpSize/2; offset > 0; offset >>=1 ) {
#if CUDART_VERSION >= 9000
    v += __shfl_down_sync(0xffffffff, v, offset);
#else
    v += __shfl_down(v, offset); 
#endif
  }
  return v;
}


//
// Block-wide sum reduction
// XXX would break with blockDim > 1024 (more than 32 warps) 
//
template <typename T>
__inline__ __device__ T block_sum_reduction(T v) {
  static __shared__ T shared[32];
  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;

  v = warp_sum_reduction(v);

  __syncthreads(); // needed when called multiple times in a row

  if (lane == 0)
    shared[wid] = v;

  __syncthreads();

  if (threadIdx.x < blockDim.x / warpSize) {
    v = shared[lane];
  } else {
    v = 0;
  }

  if (wid == 0)
    v = warp_sum_reduction(v);

  return v;
}


// 
// convert AOS (float3) coordinate layout to SOA (xxxx, yyyy, zzzz).
// 
__global__ static void 
vmd_float3_aos_to_soa(int natoms, float3 *xyz,
                      float *crdx, float *crdy, float *crdz) {
  unsigned int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int tcnt = gridDim.x * blockDim.x;
  int i;
  for (i=tid; i<natoms; i+=tcnt) {
    float3 crd = xyz[i];
    crdx[i] = crd.x;
    crdy[i] = crd.y;
    crdz[i] = crd.z;
  }
}


// 
// convert selected AOS (float3) coordinate layout to SOA (xxxx, yyyy, zzzz).
// 
__global__ static void 
vmd_float3_aos_to_soa_selected(int natoms, int *idxmap, float3 *xyz,
                               float *crdx, float *crdy, float *crdz) {
  unsigned int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int tcnt = gridDim.x * blockDim.x;
  int i;
  for (i=tid; i<natoms; i+=tcnt) {
    int dst = idxmap[i];
    if (dst > 0) {
      float3 crd = xyz[i];
      crdx[dst] = crd.x;
      crdy[dst] = crd.y;
      crdz[dst] = crd.z;
    }
  }
}


int * vmd_gpu_selection_indexlist(const AtomSel *sel) {
  int i, j;
  int *cpuidxlist = new int[sel->selected];

  for (j=0,i=sel->firstsel; i<=sel->lastsel; i++) {
    if (sel->on[i]) {
      cpuidxlist[j] = i;
      j++; 
    }
  } 

  int *gpuidxlist = NULL;
  cudaMalloc((void**) &gpuidxlist, sel->selected * sizeof(int));
  cudaMemcpy(gpuidxlist, cpuidxlist, sel->selected * sizeof(int),
             cudaMemcpyHostToDevice);
  delete [] cpuidxlist;

  return gpuidxlist;
}


//
// Device-wide QCP inner product kernel 
//   Notes:  This kernel is designed to compute a QCP inner product for
//   a single pairwise RMSD between two "large" structures, where "large"
//   means that we have sufficient parallelism (due to the total number of
//   atoms) to completely saturate GPU with work.  
//   Since the algorithm does linear work and only a relatively small number 
//   of FLOPS per atomic coordinate pair, this "large structure" case 
//   is inherently memory bandwidth bound, at least in isolation.  There may
//   an L2 cache performance benefit associated with scheduling back-to-back 
//   trajectory frame pair RMSD calculations where one of the frames is held
//   constant while looping over several others, but we could do much
//   better by writing a special-case kernel that would incorporate looping
//   internally within each thread block, thereby keeping one of the two 
//   coordinate sets entirely resident in machine registers, so long as the
//   problem was small enough to fit within the maximum CUDA 1-D grid size.
//
__global__ static void 
vmd_qcp_innerprod_soa_devicewide(double *pr,
                                 const float * RESTRICT crdx1, 
                                 const float * RESTRICT crdy1, 
                                 const float * RESTRICT crdz1,
                                 const float * RESTRICT crdx2, 
                                 const float * RESTRICT crdy2, 
                                 const float * RESTRICT crdz2,
                                 const int natoms, 
                                 const float * RESTRICT weight) {
  unsigned int tid  = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int tcnt = gridDim.x * blockDim.x;

  __shared__ int isLastBlock[1];
  if(threadIdx.x == 0) {
    isLastBlock[0] = 0;
  }
  __syncthreads();

  double G, a0, a1, a2, a3, a4, a5, a6, a7, a8;
  G=a0=a1=a2=a3=a4=a5=a6=a7=a8=0.0;

  int start = tid; 
  int stop  = natoms;
  if (weight != NULL) {
    for (int i=start; i<stop; i+=tcnt) {
      float fw = weight[i];
      float fx1 = crdx1[i];
      float fy1 = crdy1[i];
      float fz1 = crdz1[i];

      float fx2 = crdx2[i];
      float fy2 = crdy2[i];
      float fz2 = crdz2[i];

      G += fw * ((fx1*fx1 + fy1*fy1 + fz1*fz1) + (fx2*fx2 + fy2*fy2 + fz2*fz2));

      a0 += fx1 * fx2;
      a1 += fx1 * fy2;
      a2 += fx1 * fz2;

      a3 += fy1 * fx2;
      a4 += fy1 * fy2;
      a5 += fy1 * fz2;

      a6 += fz1 * fx2;
      a7 += fz1 * fy2;
      a8 += fz1 * fz2;
    }
  } else {
    for (int i=start; i<stop; i+=tcnt) {
      float fx1 = crdx1[i];
      float fy1 = crdy1[i];
      float fz1 = crdz1[i];

      float fx2 = crdx2[i];
      float fy2 = crdy2[i];
      float fz2 = crdz2[i];

      G += ((fx1*fx1 + fy1*fy1 + fz1*fz1) + (fx2*fx2 + fy2*fy2 + fz2*fz2));

      a0 += fx1 * fx2;
      a1 += fx1 * fy2;
      a2 += fx1 * fz2;

      a3 += fy1 * fx2;
      a4 += fy1 * fy2;
      a5 += fy1 * fz2;

      a6 += fz1 * fx2;
      a7 += fz1 * fy2;
      a8 += fz1 * fz2;
    }
  }

  __syncthreads();

  G *= 0.5; 

  // block-wide sum reduction of inner product
  double bG  = block_sum_reduction(G);
  double ba0 = block_sum_reduction(a0);
  double ba1 = block_sum_reduction(a1);
  double ba2 = block_sum_reduction(a2);
  double ba3 = block_sum_reduction(a3);
  double ba4 = block_sum_reduction(a4);
  double ba5 = block_sum_reduction(a5);
  double ba6 = block_sum_reduction(a6);
  double ba7 = block_sum_reduction(a7);
  double ba8 = block_sum_reduction(a8);
 
  __syncthreads();

  // thread 0 in each block writes out per-block partial sums
  if (threadIdx.x == 0) {
    pr[(0*gridDim.x)+blockIdx.x] = ba0;
    pr[(1*gridDim.x)+blockIdx.x] = ba1;
    pr[(2*gridDim.x)+blockIdx.x] = ba2;
    pr[(3*gridDim.x)+blockIdx.x] = ba3;
    pr[(4*gridDim.x)+blockIdx.x] = ba4;
    pr[(5*gridDim.x)+blockIdx.x] = ba5;
    pr[(6*gridDim.x)+blockIdx.x] = ba6;
    pr[(7*gridDim.x)+blockIdx.x] = ba7;
    pr[(8*gridDim.x)+blockIdx.x] = ba8;
    pr[(9*gridDim.x)+blockIdx.x] = bG;

    __threadfence(); // ensure all prior memory writes post before continuing

    // increment atomic counter of number of blocks that have finished
    // their work, so that the last block to finish knows to do the final
    // parallel reduction of all of the individual per-block partial sums
    unsigned int old_block_count = atomicInc(&glob_block_count, gridDim.x);
    isLastBlock[0] = ( old_block_count == (gridDim.x -1) );
  }

  __syncthreads();

  // the last block performs the final reduction of all of the individual
  // per-block partial sums
  if (isLastBlock[0]) {
    glob_block_count = 0;
    __threadfence(); // ensure block_count memory write posts before continuing

    // each thread loops over all 10 items doing the individual reductions
    for (int l=0; l < 10; ++l) {
      double *pr_a = pr+(l*gridDim.x);

      // each thread reduces all 10 items
      float sum = 0; 
      for (int j = threadIdx.x; j < gridDim.x; j += blockDim.x) {
        sum += pr_a[j];
      }

      sum = block_sum_reduction(sum);
      if (threadIdx.x==0)
        pr_a[0]=sum;
    }
  }
}


//
// Single-thread-block-per-structure-pair QCP inner product kernel 
//   Notes:  This kernel is designed to compute QCP inner products for
//   many pairwise RMSDs between "small" structures, where "small"
//   means that we have sufficient parallelism (due to the total number of
//   atoms) to supply a single CUDA thread block with sufficient work for
//   roughly 256 threads per thread block.
//   Since the algorithm does linear work and only a relatively small number 
//   of FLOPS per atomic coordinate pair, the "small structure" case 
//   is memory bandwidth bound in the worst case, but L1/L2 cache bandwidth
//   bound in the best cases.  
//   There is an opportunity to write a special case variant of this kernel
//   that stores the atomic coordinates for one of the structures in machine
//   registers such that they could be repeatedly reused, but this creates
//   some additional challenges in performing parallel sum reductions.
//
__global__ static void 
vmd_qcp_innerprod_soa_blockperpair(double *results,
                                   const float * RESTRICT crdx1, 
                                   const float * RESTRICT crdy1, 
                                   const float * RESTRICT crdz1,
                                   const float * RESTRICT crdx2, 
                                   const float * RESTRICT crdy2, 
                                   const float * RESTRICT crdz2,
                                   const int natoms, 
                                   const float * RESTRICT weight,
                                   const int num_structs) {
  unsigned int tid  = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int tcnt = gridDim.x * blockDim.x;

  double G, a0, a1, a2, a3, a4, a5, a6, a7, a8;

  // the grid of thread blocks loop over all of the structures
  for (int struct_id = blockIdx.x; struct_id < num_structs; struct_id+=gridDim.x) {
    G=a0=a1=a2=a3=a4=a5=a6=a7=a8=0.0;
    int struct_offset = struct_id * natoms * sizeof(float);
    int start = struct_offset + tid;
    int stop  = struct_offset + natoms;
    if (weight != NULL) {
      for (int i=start; i<stop; i+=tcnt) {
        float fw = weight[i];
        float fx1 = crdx1[i];
        float fy1 = crdy1[i];
        float fz1 = crdz1[i];

        float fx2 = crdx2[i];
        float fy2 = crdy2[i];
        float fz2 = crdz2[i];

        G += fw * ((fx1*fx1 + fy1*fy1 + fz1*fz1) + (fx2*fx2 + fy2*fy2 + fz2*fz2));

        a0 += fx1 * fx2;
        a1 += fx1 * fy2;
        a2 += fx1 * fz2;

        a3 += fy1 * fx2;
        a4 += fy1 * fy2;
        a5 += fy1 * fz2;

        a6 += fz1 * fx2;
        a7 += fz1 * fy2;
        a8 += fz1 * fz2;
      }
    } else {
      for (int i=start; i<stop; i+=tcnt) {
        float fx1 = crdx1[i];
        float fy1 = crdy1[i];
        float fz1 = crdz1[i];

        float fx2 = crdx2[i];
        float fy2 = crdy2[i];
        float fz2 = crdz2[i];

        G += ((fx1*fx1 + fy1*fy1 + fz1*fz1) + (fx2*fx2 + fy2*fy2 + fz2*fz2));

        a0 += fx1 * fx2;
        a1 += fx1 * fy2;
        a2 += fx1 * fz2;

        a3 += fy1 * fx2;
        a4 += fy1 * fy2;
        a5 += fy1 * fz2;

        a6 += fz1 * fx2;
        a7 += fz1 * fy2;
        a8 += fz1 * fz2;
      }
    }

    __syncthreads();

    G *= 0.5; 

    // block-wide sum reduction of inner product
    double bG  = block_sum_reduction(G);
    double ba0 = block_sum_reduction(a0);
    double ba1 = block_sum_reduction(a1);
    double ba2 = block_sum_reduction(a2);
    double ba3 = block_sum_reduction(a3);
    double ba4 = block_sum_reduction(a4);
    double ba5 = block_sum_reduction(a5);
    double ba6 = block_sum_reduction(a6);
    double ba7 = block_sum_reduction(a7);
    double ba8 = block_sum_reduction(a8);

    __syncthreads();

    // thread 0 in each block writes out per-block partial sums
    if (threadIdx.x == 0) {
      int addr = struct_id * 10;
      results[addr    ] = ba0;
      results[addr + 1] = ba1;
      results[addr + 2] = ba2;
      results[addr + 3] = ba3;
      results[addr + 4] = ba4;
      results[addr + 5] = ba5;
      results[addr + 6] = ba6;
      results[addr + 7] = ba7;
      results[addr + 8] = ba8;
      results[addr + 9] = bG;
    }
  }
}


// compute lower-triangular indices i,j from linear array index
static int idx2sub_tril(long N, long ind, long *J, long *I) {
  if (ind > (N*(N+1)/2)) {
    return -1; // out of bounds
  }

#if 1
  long i = long(N - 2 - floor(sqrt(-8*ind + 4*N*(N-1)-7)/2.0 - 0.5));
  long j = ind + i + 1 - N*(N-1)/2 + (N-i)*((N-i)-1)/2;
#elif 0
  long i, j;
  // j = ceil(0.5*(2*N+1 - sqrt(-8*ind + 1 + 4*N*(N+1))));
  // i = ind - (j-1)*N + j*(j-1)/2;

  // i = floor((2*N+1 - sqrt((2*N+1)*(2*N+1) - 8*ind)) / 2);
  // XXX deal with ambiguous types for sqrt() on Solaris
  double tmp2np1 = 2*N+1;
  i = floor((tmp2np1 - sqrt(tmp2np1*tmp2np1 - 8.0*ind)) / 2);
  j = ind - N*i + i*(i-1)/2 + i;
#endif

  *I = i;
  *J = j;

  return 0;
}


static void * measure_rmsdmat_qcp_thread(void *voidparms) {
  int threadid;
  gpuqcprmsdthreadparms *parms = NULL;
  wkf_threadpool_worker_getdata(voidparms, (void **) &parms);
  wkf_threadpool_worker_getid(voidparms, &threadid, NULL);

  //
  // copy in per-thread parameters
  //
  float *rmsdmat = parms->rmsdmat;

#if 0
  int framecrdsz = parms->framecrdsz;
  float *crds = parms->crds;
#endif
  const AtomSel *sel = parms->sel;
  int first  = parms->first;
  int last   = parms->last;
  int step   = parms->step;

#if 1
  printf("QCP GPU[%d] running...\n", threadid);
#endif

  int framecount = (last - first + 1) / step;

  wkf_tasktile_t tile;
  while (wkf_threadlaunch_next_tile(voidparms, 8, &tile) != WKF_SCHED_DONE) {
    long idx;

    for (idx=tile.start; idx<tile.end; idx++) {
      long i, j;

      // compute i,j from idx...
      // only compute off-diagonal elements, so we use (framecount-1)
      if (idx2sub_tril(framecount-1, idx, &i, &j)) {
        printf("qcpthread[%d]: work idx %ld out of triangle!\n", threadid, idx);
        break;
      }

      // calculate the (weighted) inner product of two structures
      double A[9];
      double E0 = 0;

#if 0
      float *xj = crds + (j * 3 * framecrdsz);
      float *yj = xj + framecrdsz;
      float *zj = xj + framecrdsz*2;

      float *xi = crds + (i * 3 * framecrdsz);
      float *yi = xi + framecrdsz;
      float *zi = xi + framecrdsz*2;

      E0 = InnerProductSOA(A, xj, yj, zj, xi, yi, zi,
                           sel->selected, NULL /* weight */);
#endif

      // calculate the RMSD & rotational matrix
      FastCalcRMSDAndRotation(NULL, A, &rmsdmat[j*framecount + i],
                              E0, sel->selected, -1);
    }
  }

  return NULL;
}


int qcp_soa_gpu(wkf_threadpool_t *devpool, // VMD GPU worker thread pool
                const AtomSel *sel, float *crds, int framecrdsz, int padcnt,
                int first, int last, int step, int framecount,
                float *rmsdmat) {
  //
  // copy in per-thread parameters
  //
  gpuqcprmsdthreadparms parms;
  memset(&parms, 0, sizeof(parms));

  parms.sel = sel;
  parms.rmsdmat = rmsdmat;
  parms.padcnt = padcnt;
  parms.framecrdsz = framecrdsz;
  parms.crds = crds;
  parms.first = first;
  parms.last = last;
  parms.step = step;

  // spawn child threads to do the work
  wkf_tasktile_t tile;
  tile.start=0;
  tile.end=(framecount-1)*(framecount-1)/2; // only compute off-diag elements

#if 1
  wkf_threadpool_sched_dynamic(devpool, &tile);
  wkf_threadpool_launch(devpool, measure_rmsdmat_qcp_thread, &parms, 1);
#else
  int i, j;
  for (j=0; j<framecount; j++) {
    float *xj = crds + (j * 3 * framecrdsz);
    float *yj = xj + framecrdsz;
    float *zj = xj + framecrdsz*2;
    for (i=0; i<j; i++) {
      // calculate the (weighted) inner product of two structures
      double A[9];

      float *xi = crds + (i * 3 * framecrdsz);
      float *yi = xi + framecrdsz;
      float *zi = xi + framecrdsz*2;

      double E0 = InnerProductSOA(A, xj, yj, zj, xi, yi, zi,
                                  sel->selected, NULL /* weight */);

      // calculate the RMSD & rotational matrix
      FastCalcRMSDAndRotation(NULL, A, &rmsdmat[j*framecount + i],
                              E0, sel->selected, -1);
    }
  }
#endif

  return 0;
}  



#if defined(VMDUSECUDAGDS)
#define VMDGDSMAXFRAMEBUF     8
static void * measure_rmsdmat_qcp_ooc_thread(void *voidparms) {
  int threadid;
  gpuqcprmsdmatoocthreadparms *parms = NULL;
  wkf_threadpool_worker_getdata(voidparms, (void **) &parms);
  wkf_threadpool_worker_getid(voidparms, &threadid, NULL);

  //
  // copy in per-thread parameters
  //
  float *rmsdmat = parms->rmsdmat;

  int nfiles = parms->nfiles;
  int natoms = parms->natoms;
  const AtomSel *sel = parms->sel;
  int first  = parms->first;
  int last   = parms->last;
  int step   = parms->step;

  int usecufile = 1;
  fio_fd *hostfds = NULL;
  int pinhostiobuffer = 1;
  if (getenv("VMDGDSHOSTNOPIN")) {
    pinhostiobuffer=0; 
  }

  if (getenv("VMDGDSUSEHOST")) {
    usecufile=0;
    hostfds = (fio_fd *) calloc(1, nfiles * sizeof(fio_fd));

    int hostusedirectio = 1;
    if (getenv("VMDGDSHOSTBUFFEREDIO") != NULL)
      hostusedirectio = 0;
 
    int openmode = FIO_READ;
    if (hostusedirectio) 
      openmode |= FIO_DIRECT;

    int i;
    for (i=0; i<nfiles; i++) {
      if (fio_open(parms->trjfileset[i], openmode, &hostfds[i]) < 0) {
        if (hostusedirectio) {
          printf("Thr[%d] direct I/O unavailable or can't open file '%s'\n", 
                 threadid, parms->trjfileset[i]);
        } else {
          printf("Thr[%d] can't open file '%s'\n", 
                 threadid, parms->trjfileset[i]);
        }
        return NULL;
      }
#if 0
      else {
        printf("Thr[%d]: fd[%d]: %d, '%s'\n", threadid, i, hostfds[i],
               parms->trjfileset[i]);
      }
#endif
    }
  }

  if (hostfds && usecufile) {
    printf("Inconsistent cufile/hostfds state, aborting!\n");
    return NULL;
  }

  /* ensure we have a large enough allocation so we can align */
  /* the starting pointer to a blocksz page boundary          */
  long blocksz = MOLFILE_DIRECTIO_MIN_BLOCK_SIZE;
  long sz = 3L*sizeof(float)*natoms + blocksz;

  /* pad the allocation to an even multiple of the block size */
  size_t blockpadsz = (sz + (blocksz - 1)) & (~(blocksz - 1));

  int framecount = (last - first + 1) / step;
  int framesperfile = framecount / nfiles;

  if (threadid == 0) {
    printf("Thr[%2d] running frame: %d natoms: %d selected: %d\n", 
           framecount, threadid, natoms, sel->selected);
  }

  cudaError_t crc;
  cudaStream_t qcpstream;
  double *pr=NULL;
  float *devptr=NULL;
  float *hostptr=NULL;
  float *hostptr_unaligned=NULL;

  float *crdx1=NULL, *crdy1=NULL, *crdz1=NULL;
  float *crdx2=NULL, *crdy2=NULL, *crdz2=NULL;
  int *gpuidxlist = NULL;
  int multiframeio = 0;
  if (getenv("VMDGDSMULTIFRAME"))
    multiframeio = atoi(getenv("VMDGDSMULTIFRAME"));
  if (multiframeio > VMDGDSMAXFRAMEBUF)
    multiframeio = VMDGDSMAXFRAMEBUF;

  // set block sizes and counts for QCP calcs
  dim3 QCPBsz = dim3(256, 1, 1);
  dim3 QCPGsz = dim3((natoms + QCPBsz.x - 1) / QCPBsz.x, 1, 1);

  if (parms->devcount > 0) {
    long gpuallocsz = (VMDGDSMAXFRAMEBUF+1) * blockpadsz;

    if (threadid == 0) {
      printf("Thr[%2d] Allocating QCP GPU timestep buf: %ld \n", 
             threadid, gpuallocsz);
    }
    crc = cudaMalloc((void**) &devptr, gpuallocsz);

    if (hostfds != NULL) {
      if (pinhostiobuffer) {
        crc = cudaMallocHost((void**) &hostptr, gpuallocsz);
        CUERR
      } else {
        hostptr = (float *) alloc_aligned_ptr(gpuallocsz, 4096, 
                                              (void**) &hostptr_unaligned);
      }
    }

    gpuidxlist = vmd_gpu_selection_indexlist(sel);
    long crdsz = sel->selected * sizeof(float);
    // 10x QCP sums
    crc = cudaMalloc((void**) &pr, 10 * sizeof(double) * (QCPGsz.x + 1));

    // atomic coord buffers
    crc = cudaMalloc((void**) &crdx1, crdsz);
    crc = cudaMalloc((void**) &crdy1, crdsz);
    crc = cudaMalloc((void**) &crdz1, crdsz);
    crc = cudaMalloc((void**) &crdx2, crdsz);
    crc = cudaMalloc((void**) &crdy2, crdsz);
    crc = cudaMalloc((void**) &crdz2, crdsz);
    if (crc != cudaSuccess) {
      printf("Thr[%2d], Failed to allocate GPU buffer!\n", threadid);
      return NULL; // XXX error handling needs to be done here
    }

    cudaStreamCreate(&qcpstream);

#if defined(VMDUSECUDAGDS)
    cuFileBufRegister(devptr, gpuallocsz, 0);
#endif
  }


#if 1
#define VMDGDSIOBENCHMARKONLY 1
  int verbose = (getenv("VMDGDSVERBOSE") != NULL) ? 1 : 0;

  wkf_tasktile_t tile;
  while (wkf_threadlaunch_next_tile(voidparms, VMDGDSMAXFRAMEBUF * 1, &tile) != WKF_SCHED_DONE) {
    //
    // simple I/O + compute benchmarking...
    //
    int idx;
    for (idx=tile.start; idx<tile.end; idx++) {
      int myfileidx = (threadid / 4) % nfiles;
      int fileframeidx = idx % framesperfile;

      //
      // compute multi-frame or single-frame I/O offsets and sizes
      //
      long startoffset, foffset, readlen;
      read_js_timestep_index_offsets(parms->jshandles[myfileidx], natoms, 
                                     fileframeidx, 
                                     0, natoms, NULL,
                                     &startoffset, &foffset, &readlen);
      if (multiframeio) {
        // multi-frame reads use the same starting offset, but the
        // read length is computed from the first and last frames
        // in the group.
        long multistartoffset, multifoffset, multireadlen;
        read_js_timestep_index_offsets(parms->jshandles[myfileidx], natoms, 
                                       fileframeidx+multiframeio-1, 
                                       0, natoms, NULL,
                                       &multistartoffset, &multifoffset, 
                                       &multireadlen);

        multireadlen = (multifoffset + multireadlen) - foffset;

        //printf("** readlen: %ld  multireadlen: %ld\n", readlen, multireadlen);
        readlen = multireadlen;
        idx+=multiframeio-1; // add in the required increment...
      }

      //
      // perform the required I/O via GDS or by host kernel I/O
      //
      long ret=0;
      if (usecufile) {
        ret = cuFileRead(parms->cfh[myfileidx], (char *) devptr, readlen, foffset, 0);
      } else if (hostfds) {
        foffset=0;
        ret=fio_fseek(hostfds[myfileidx], foffset, FIO_SEEK_SET);
        if (ret<0) { printf("fio_fseek() error!\n"); return NULL; }
        ret=fio_fread(hostptr, readlen, 1, hostfds[myfileidx]);
        if (ret<0) { printf("fio_fseek() error!\n"); return NULL; }
        cudaMemcpy(devptr, hostptr, readlen, cudaMemcpyHostToDevice);
        CUERR
      } else {
        printf("Inconsistent cufile/hostfds state, aborting!\n");
        return NULL;
      }

      // handle errors if they have occured
      if (ret < 0) {
        printf("Thr[%2d] Error: cuFileRead(): %ld\n", threadid, ret);
        return NULL; // XXX error handling needs to be done here
      }

#if 0
      float hostbuf[64];
      cudaMemcpy(hostbuf, devptr, sizeof(hostbuf), cudaMemcpyDeviceToHost);
      int l;
      for (l=0; l<10; l++) {
        printf("%.1f ", hostbuf[l]); 
      }
      printf("\n");
#endif

      float rmsd=0;
#if 0
      cudaDeviceSynchronize();
      CUERR // check for errors

      dim3 Bsz = dim3(256, 1, 1);
      dim3 Gsz = dim3((natoms + Bsz.x - 1) / Bsz.x, 1, 1);
      if (sel->selected == natoms) {
        vmd_float3_aos_to_soa<<<Gsz, Bsz, 0, qcpstream>>>(natoms, (float3 *) devptr, crdx1, crdy1, crdz1);
      } else {
        vmd_float3_aos_to_soa_selected<<<Gsz, Bsz, 0, qcpstream>>>(natoms, gpuidxlist, (float3 *) devptr, crdx1, crdy1, crdz1);
      }

      vmd_qcp_innerprod_soa_devicewide<<<QCPGsz, QCPBsz, 0, qcpstream>>>(pr, crdx1, crdy1, crdz1, crdx1, crdy1, crdz1, natoms, NULL);

      double pr_host[10];
      cudaMemcpy(pr_host, pr, 10 * sizeof(double), cudaMemcpyDeviceToHost);
#if 0
      printf("Thr[%2d] pr_host %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n", 
             threadid, pr_host[0], pr_host[1], pr_host[2],
             pr_host[3], pr_host[4], pr_host[5],
             pr_host[6], pr_host[7], pr_host[8], pr_host[9]);
#endif

      // calculate the RMSD & rotational matrix
      double E0 = pr_host[9];
      FastCalcRMSDAndRotation(NULL, pr_host, &rmsd, E0, sel->selected, -1);
      //rmsdmat[j*framecount + i] = rmsd;
#endif

      if (verbose) {
        printf("Thr[%2d]F[%d][tile: %d to %d] frame: %d RMSD: %4.1f cuFile len: %ld off: %ld\n", 
               threadid, myfileidx, tile.start, tile.end, idx, rmsd, 
               readlen, foffset); 
      }  
    }
  }

#else

  // 
  // Real QCP loops
  // 
  int xfileidx = -1;
  int yfileidx = -1;
  wkf_tasktile_t tile;
  long lastx = -1, lasty = -1;
  while (wkf_threadlaunch_next_tile(voidparms, 8, &tile) != WKF_SCHED_DONE) {
    int idx;
    for (idx=tile.start; idx<tile.end; idx++) {

      //
      // Perform I/O to get coordinates into the GPU
      //

      // compute i,j from idx...
      // only compute off-diagonal elements, so we use (framecount-1)
      long x=-1, y=-1, ret=0;
      if (idx2sub_tril(framecount-1, idx, &x, &y)) {
        printf("qcpthread[%d]: work idx %d out of triangle!\n", threadid, idx);
        break;
      }

      long startoffset, foffset, readlen;
//printf("Thr[%2d] pair(%3ld x %3ld) cuFileRead(): offset %ld  readlen: %ld\n", 
//       threadid, x, y, foffset, readlen);

      // if the "row" frame has changed, read in the new one
      if (lasty != y) {
        // XXX Hard-coded DGX-2 specific topology-file locality scheme.
        //     Topology determination should be entirely runtime-determined
        yfileidx = (threadid / 4) % nfiles;
        long fileframeidx = y % framesperfile;
        read_js_timestep_index_offsets(parms->jshandles[yfileidx], 
                                       natoms, fileframeidx, 
                                       0, natoms, NULL,
                                       &startoffset, &foffset, &readlen);
printf("Thr[%2d] pair(%3ld x %3ld) cuFileRead(Y): offset %ld  readlen: %ld\n", 
       threadid, x, y, foffset, readlen);
        lasty = y;
      }

      // if no errs and the "column" frame has changed, read in the new one
      if (ret >= 0 && lastx != x) {
        // XXX Hard-coded DGX-2 specific topology-file locality scheme.
        //     Topology determination should be entirely runtime-determined
        xfileidx = (threadid / 4) % nfiles;
        long fileframeidx = x % framesperfile;
        read_js_timestep_index_offsets(parms->jshandles[xfileidx], 
                                       natoms, fileframeidx, 
                                       0, natoms, NULL,
                                       &startoffset, &foffset, &readlen);
printf("Thr[%2d] pair(%3ld x %3ld) cuFileRead(X): offset %ld  readlen: %ld\n", 
       threadid, x, y, foffset, readlen);
        ret = cuFileRead(parms->cfh[xfileidx], (char *) devptr+blockpadsz, readlen, foffset);
          lastx = x;
      }

      // handle errors if they have occured
      if (ret < 0) {
        const char *descp = "unknown error code";
#if 0
        // XXX this requires linkage with libcuda.so, which is
        //     against convention, so we avoid it for now
        if (cuGetErrorName(ret, &descp) != CUDA_SUCCESS)
          descp = "unknown cuda error";
#endif

        printf("Thr[%2d] Error: cuFileRead(): %ld, '%s'\n", threadid, ret, descp);
        return NULL; // XXX error handling needs to be done here
      }

      //
      // Produce contiguous coordinates based on active atom selection
      //

      // 
      // Run QCP calculation(s)
      // 
      cudaDeviceSynchronize();
      CUERR // check for errors

      dim3 Bsz = dim3(256, 1, 1);
      dim3 Gsz = dim3((natoms + Bsz.x - 1) / Bsz.x, 1, 1);
      if (sel->selected == natoms) {
        vmd_float3_aos_to_soa<<<Gsz, Bsz, 0, qcpstream>>>(natoms, (float3 *) devptr, crdx1, crdy1, crdz1);
      } else {
        vmd_float3_aos_to_soa_selected<<<Gsz, Bsz, 0, qcpstream>>>(natoms, gpuidxlist, (float3 *) devptr, crdx1, crdy1, crdz1);
      }

      vmd_qcp_innerprod_soa_devicewide<<<QCPGsz, QCPBsz, 0, qcpstream>>>(pr, crdx1, crdy1, crdz1, crdx1, crdy1, crdz1, natoms, NULL);

#if 0
      double pr_host[10];
      cudaMemcpy(pr_host, pr, 10 * sizeof(double), cudaMemcpyDeviceToHost);

      // calculate the RMSD & rotational matrix
      float rmsd=0;
      double E0 = pr_host[9];
      FastCalcRMSDAndRotation(NULL, pr_host, &rmsd, E0, sel->selected, -1);
      //rmsdmat[j*framecount + i] = rmsd;
#endif
    }
  }
#endif

  cudaFree(gpuidxlist);
  cudaFree(pr);
  cudaFree(crdx1);
  cudaFree(crdy1);
  cudaFree(crdz1);
  cudaFree(crdx2);
  cudaFree(crdy2);
  cudaFree(crdz2);

#if defined(VMDUSECUDAGDS)
  if (usecufile) {
    cuFileBufDeregister(devptr);
  }
#endif

  if (hostfds != NULL) {
    int i;
    for (i=0; i<nfiles; i++) {
      fio_fclose(hostfds[i]);
    }
    free(hostfds);
  }

  if (hostptr != NULL) {
    if (pinhostiobuffer) {
      cudaFreeHost(hostptr);
    } else {
      free(hostptr_unaligned);
    }
  }

  return NULL;
}

#endif


int qcp_soa_gpu_ooc(wkf_threadpool_t *devpool, // VMD GPU worker thread pool
                    int nfiles, const char **trjfileset, 
                    const AtomSel *sel,
                    int first, int last, int step, float *rmsdmat) {
  printf("qcp_soa_gpu_ooc()\n");
  wkf_threadpool_t *bigpool = NULL;

#if defined(VMDUSECUDAGDS)
  int devcount;
  cudaError_t crc = cudaGetDeviceCount(&devcount);
  printf("  GPU device count: %d\n", devcount);
  if (devcount==0)
    printf("  No GPU devices, continuing with host only...\n");

  CUfileHandle_t * cfh = (CUfileHandle_t *) calloc(1, nfiles * sizeof(CUfileHandle_t));
  CUfileDescr_t * cfhdesc = (CUfileDescr_t *) calloc(1, nfiles * sizeof(CUfileDescr_t)); 
  memset(&cfh[0], 0, sizeof(cfh));
  memset(&cfhdesc[0], 0, sizeof(cfhdesc));

  int natomschk = 0;
  jshandle **jshandles = (jshandle **) calloc(1, nfiles * sizeof(jshandle *));
  fio_fd *directio_fds = (fio_fd *) calloc(1, nfiles * sizeof(fio_fd));

  int i;
  for (i=0; i<nfiles; i++) {
    const char *filename = trjfileset[i];
    printf("File[%d] GDS setup, opening '%s'\n", i, filename);
    jshandles[i] = (jshandle *) open_js_read(filename, "js", &natomschk);
    if (!jshandles[i]) {
      printf("File[%d] open_js_read failed for file %s\n", i, filename);
      return -1; // deal with error handling later
    }

#if vmdplugin_ABIVERSION > 17
    long blocksz = MOLFILE_DIRECTIO_MIN_BLOCK_SIZE;
    int filepgalignsz = 1;
    read_js_timestep_pagealign_size(jshandles[i], &filepgalignsz);
    if (filepgalignsz != blocksz) {
      printf("File[%d] Plugin-returned page alignment size mismatch!\n", i);
    } else {
      printf("File[%d] Page alignment size: %d\n", i, filepgalignsz);
    }
#endif

    read_js_timestep_index_offsets(jshandles[i], natomschk, 0, 0, 0, 
                                   &directio_fds[i], NULL, NULL, NULL);

    cfhdesc[i].handle.fd = directio_fds[i]; // typedef of Unix FD
    cfhdesc[i].type = CU_FILE_HANDLE_TYPE_OPAQUE_FD;
    CUfileError_t cferr = cuFileHandleRegister(&cfh[i], &cfhdesc[i]);

    if (cferr.err != CU_FILE_SUCCESS) {
      printf("File[%d] cuFileImportExternalFile on fd %d failed!\n", 
             i, cfhdesc[i].handle.fd);
      return -1; // XXX error handling needs to be done here
    }
  }


  //
  // copy in per-thread parameters
  //
  gpuqcprmsdmatoocthreadparms parms;
  memset(&parms, 0, sizeof(parms));
  parms.devcount = devcount;
  parms.nfiles = nfiles;
  parms.trjfileset = trjfileset;
  parms.jshandles = jshandles;
  parms.cfh = cfh;
  parms.natoms = sel->num_atoms;
  parms.sel = sel;
  parms.rmsdmat = rmsdmat;
  parms.first = first;
  parms.last = last;
  parms.step = step;

  int framecount = nfiles * (last / step);

  // create timers
  wkf_timerhandle timer;
  timer=wkf_timer_create();

  // spawn child threads to do the work
  wkf_tasktile_t tile;
  tile.start=0;
#if defined(VMDGDSIOBENCHMARKONLY)
  tile.end=framecount - 1; // first row only
#elif 1
  tile.end=framecount*16 - 1; // first 16 rows only
#else
  tile.end=(framecount-1)*(framecount-1)/2; // only compute off-diag elements
#endif

printf("** qcp_soa_gpu_ooc(): tile start: %d  end: %d\n", tile.start, tile.end);
 
  int gdsthreadspergpu = 1;
  if (getenv("VMDGDSTHREADSPERGPU") != NULL) 
    gdsthreadspergpu = atoi(getenv("VMDGDSTHREADSPERGPU"));

printf("** gdsthreadspergpu: %d\n", gdsthreadspergpu);

  if (gdsthreadspergpu > 1) {
    // XXX extra-large GPU device thread pool
    int workercount = devcount * gdsthreadspergpu;

    int *devlist = new int[workercount];
    int k;
    for (k=0; k<workercount; k++) {
      devlist[k] = k / gdsthreadspergpu; // XXX ignores VMD CUDA device masks
    }

    msgInfo << "Creating Multi-worker (" 
            << gdsthreadspergpu << " per-GPU) CUDA device pool..." << sendmsg;
    bigpool=wkf_threadpool_create(workercount, devlist);
    delete [] devlist;

    // associate each worker thread with a specific GPU
    if (getenv("VMDCUDAVERBOSE") != NULL)
      wkf_threadpool_launch(bigpool, vmd_cuda_devpool_setdeviceonly, (void*)"VMD CUDA Dev Init", 1);
    else
      wkf_threadpool_launch(bigpool, vmd_cuda_devpool_setdeviceonly, NULL, 1);

    // clear all available device memory on each of the GPUs
    wkf_threadpool_launch(bigpool, vmd_cuda_devpool_clear_device_mem, NULL, 1);

    // XXX override which GPU device pool we're going to use
    devpool = bigpool;
  }

  // XXX affinitize GPU worker threads for best perf 
  wkf_threadpool_launch(devpool, vmd_cuda_affinitize_threads, NULL, 1);

  wkf_threadpool_sched_dynamic(devpool, &tile);
  wkf_timer_start(timer);
  wkf_threadpool_launch(devpool, measure_rmsdmat_qcp_ooc_thread, &parms, 1);
  wkf_timer_stop(timer);

  double runtime = wkf_timer_time(timer); 
  double gbytes = sel->num_atoms * 12L * (tile.end+1) / (1024.0 * 1024.0 * 1024.0);

  printf("QCP natoms: %d, fsz: %ld, tsz: %ld\n", 
         sel->num_atoms, sel->num_atoms * 12L, 
         sel->num_atoms * 12L * (tile.end+1)); 

  int pinhostiobuffer = 1;
  if (getenv("VMDGDSHOSTNOPIN"))
    pinhostiobuffer=0; 

  int hostusedirectio = 1;
  if (getenv("VMDGDSHOSTBUFFEREDIO") != NULL)
    hostusedirectio = 0;

  int usecufile=1;
  if (getenv("VMDGDSUSEHOST"))
    usecufile=0;

  if (usecufile) {
    printf("OOC I/O via GDS + cuFile\n");
  } else {
    printf("OOC I/O via host, %s APIs, %s memory buffers\n",
           (hostusedirectio) ? "Direct I/O" : "Buffered I/O",
           (pinhostiobuffer) ? "pinned" : "unpinned");
  } 

  int multiframeio = 0;
  if (getenv("VMDGDSMULTIFRAME"))
    multiframeio = atoi(getenv("VMDGDSMULTIFRAME"));
  if (multiframeio > VMDGDSMAXFRAMEBUF)
    multiframeio = VMDGDSMAXFRAMEBUF;
  if (multiframeio) {
    printf("QCP GDS multi-frame read opt: %d frames per call, %ld bytes\n",
           multiframeio, 
           multiframeio * sel->num_atoms * 12L);
  }

  printf("QCP runtime: %.1f, %.2fGB/sec\n", runtime, gbytes/runtime);
         
  for (i=0; i<nfiles; i++) { 
#if defined(VMDUSECUDAGDS)
    cuFileHandleDeregister(cfh[i]);
#endif
    close_js_read(jshandles[i]);
  }
#endif

#if defined(VMDUSECUDAGDS)
  if (cfh != NULL)
    free(cfh);

  if (cfhdesc != NULL)
     free(cfhdesc);

  if (jshandles != NULL)
    free(jshandles);

  if (directio_fds != NULL)
    free(directio_fds);
#endif

  // if we created an extra large-thread-count-per-GPU thread pool, we 
  // need to destroy it here...
  if (bigpool != NULL) 
    wkf_threadpool_destroy(bigpool);

  return 0;
}  
