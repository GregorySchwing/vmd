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
 *      $RCSfile: CUDABench.cu,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.42 $      $Date: 2022/02/09 04:03:19 $
 *
 ***************************************************************************/
/**
 * \file CUDABench.cu
 * \brief Short benchmark kernels to measure GPU performance.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>

#include "Inform.h"
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
#endif


#define CUERR { cudaError_t err; \
  if ((err = cudaGetLastError()) != cudaSuccess) { \
  printf("CUDA error: %s, %s line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
  return -1; }}


//
// Benchmark peak Multiply-Add instruction performance, in GFLOPS
//

// FMADD16 macro contains a sequence of operations that the compiler
// won't optimize out, and will translate into a densely packed block
// of multiply-add instructions with no intervening register copies/moves
// or other instructions. 
#define FMADD16 \
    tmp0  = tmp0*tmp4+tmp7;     \
    tmp1  = tmp1*tmp5+tmp0;     \
    tmp2  = tmp2*tmp6+tmp1;     \
    tmp3  = tmp3*tmp7+tmp2;     \
    tmp4  = tmp4*tmp0+tmp3;     \
    tmp5  = tmp5*tmp1+tmp4;     \
    tmp6  = tmp6*tmp2+tmp5;     \
    tmp7  = tmp7*tmp3+tmp6;     \
    tmp8  = tmp8*tmp12+tmp15;   \
    tmp9  = tmp9*tmp13+tmp8;    \
    tmp10 = tmp10*tmp14+tmp9;   \
    tmp11 = tmp11*tmp15+tmp10;  \
    tmp12 = tmp12*tmp8+tmp11;   \
    tmp13 = tmp13*tmp9+tmp12;   \
    tmp14 = tmp14*tmp10+tmp13;  \
    tmp15 = tmp15*tmp11+tmp14;

// CUDA grid, thread block, loop, and MADD operation counts
#define GRIDSIZEX       6144  // number of 1-D thread blocks
#define BLOCKSIZEX      64    // number of threads per 1-D block
#define GLOOPS          2000  // iteration count (all threads)
#define FMADD16COUNT    32    // 32 reps
#define FLOPSPERFMADD16 32    // 16 MULs and 16 ADDs

// FLOP counting
#define FLOPSPERLOOP (FMADD16COUNT * FLOPSPERFMADD16)

__global__ static void madd_kernel(float *doutput) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  float tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
  float tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15;
  tmp0=tmp1=tmp2=tmp3=tmp4=tmp5=tmp6=tmp7=0.0f;
  tmp8=tmp9=tmp10=tmp11=tmp12=tmp13=tmp14=tmp15 = 0.0f;

  tmp15=tmp7 = blockIdx.x * 0.001f; // prevent compiler from optimizing out
  tmp1 = blockIdx.y * 0.001f;       // the body of the loop...

  int loop;
  for(loop=0; loop<GLOOPS; loop++){
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
    FMADD16
  }

  doutput[tid] = tmp0+tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7
                 +tmp8+tmp9+tmp10+tmp11+tmp12+tmp13+tmp14+tmp15;
}


static int cudamaddgflops(int cudadev, double *gflops, int testloops) {
  float *doutput = NULL;
  dim3 Bsz, Gsz;
  wkf_timerhandle timer;
  int i;

  cudaError_t rc;
  rc = cudaSetDevice(cudadev);
  if (rc != cudaSuccess) {
#if CUDART_VERSION >= 2010
    rc = cudaGetLastError(); // query last error and reset error state
    if (rc != cudaErrorSetOnActiveProcess)
      return -1; // abort and return an error
#else
    cudaGetLastError(); // just ignore and reset error state, since older CUDA
                        // revs don't have a cudaErrorSetOnActiveProcess enum
#endif
  }


  // setup CUDA grid and block sizes
  Bsz.x = BLOCKSIZEX;
  Bsz.y = 1;
  Bsz.z = 1;
  Gsz.x = GRIDSIZEX;
  Gsz.y = 1;
  Gsz.z = 1;

  // allocate output array
  cudaMalloc((void**)&doutput, BLOCKSIZEX * GRIDSIZEX * sizeof(float));
  CUERR // check and clear any existing errors

  // warmup run
  madd_kernel<<<Gsz, Bsz>>>(doutput);
  cudaDeviceSynchronize(); // wait for kernel to finish

  // benchmark run
  timer=wkf_timer_create();
  wkf_timer_start(timer);
  for (i=0; i<testloops; i++) { 
    madd_kernel<<<Gsz, Bsz>>>(doutput);
  }
  cudaDeviceSynchronize(); // wait for kernel to finish
  CUERR // check and clear any existing errors
  wkf_timer_stop(timer);

  double runtime = wkf_timer_time(timer);
  double gflop = ((double) GLOOPS) * ((double) FLOPSPERLOOP) *
                  ((double) BLOCKSIZEX) * ((double) GRIDSIZEX) * (1.0e-9) * testloops;
  
  *gflops = gflop / runtime;

  cudaFree(doutput);
  CUERR // check and clear any existing errors

  wkf_timer_destroy(timer);

  return 0;
}

typedef struct {
  int deviceid;
  int testloops;
  double gflops;
} maddthrparms;

static void * cudamaddthread(void *voidparms) {
  maddthrparms *parms = (maddthrparms *) voidparms;
  cudamaddgflops(parms->deviceid, &parms->gflops, parms->testloops);
  return NULL;
}

int vmd_cuda_madd_gflops(int numdevs, int *devlist, double *gflops,
                         int testloops) {
  maddthrparms *parms;
  wkf_thread_t * threads;
  int i;

  /* allocate array of threads */
  threads = (wkf_thread_t *) calloc(numdevs * sizeof(wkf_thread_t), 1);

  /* allocate and initialize array of thread parameters */
  parms = (maddthrparms *) malloc(numdevs * sizeof(maddthrparms));
  for (i=0; i<numdevs; i++) {
    if (devlist != NULL)
      parms[i].deviceid = devlist[i];
    else
      parms[i].deviceid = i;

    parms[i].testloops = testloops;
    parms[i].gflops = 0.0;
  }

#if defined(VMDTHREADS)
  /* spawn child threads to do the work */
  /* thread 0 must also be processed this way otherwise    */
  /* we'll permanently bind the main thread to some device */
  for (i=0; i<numdevs; i++) {
    wkf_thread_create(&threads[i], cudamaddthread, &parms[i]);
  }

  /* join the threads after work is done */
  for (i=0; i<numdevs; i++) {
    wkf_thread_join(threads[i], NULL);
  }
#else
  /* single thread does all of the work */
  cudamaddthread((void *) &parms[0]);
#endif

  for (i=0; i<numdevs; i++) {
    gflops[i] = parms[i].gflops; 
  }

  /* free thread parms */
  free(parms);
  free(threads);

  return 0;
}






//
// Host-GPU memcpy I/O bandwidth benchmark
//

#define BWITER 500
#define LATENCYITER 50000

static int cudabusbw(int cudadev, 
                     double *hdmbsec, double *hdlatusec, 
                     double *phdmbsec, double *phdlatusec, 
                     double *dhmbsec, double *dhlatusec,
                     double *pdhmbsec, double *pdhlatusec) {
  float *hdata = NULL;   // non-pinned DMA buffer
  float *phdata = NULL;  // pinned DMA buffer
  float *ddata = NULL;
  int i;
  double runtime;
  wkf_timerhandle timer;
  int memsz = 1024 * 1024 * sizeof(float);

  *hdmbsec = 0.0;
  *hdlatusec = 0.0;
  *dhmbsec = 0.0;
  *dhlatusec = 0.0;
  *phdmbsec = 0.0;
  *phdlatusec = 0.0;
  *pdhmbsec = 0.0;
  *pdhlatusec = 0.0;

  // attach to the selected device
  cudaError_t rc;
  rc = cudaSetDevice(cudadev);
  if (rc != cudaSuccess) {
#if CUDART_VERSION >= 2010
    rc = cudaGetLastError(); // query last error and reset error state
    if (rc != cudaErrorSetOnActiveProcess)
      return -1; // abort and return an error
#else
    cudaGetLastError(); // just ignore and reset error state, since older CUDA
                        // revs don't have a cudaErrorSetOnActiveProcess enum
#endif
  }

  // allocate non-pinned output array
  hdata = (float *) malloc(memsz); 

  // allocate pinned output array
  cudaMallocHost((void**) &phdata, memsz);
  CUERR // check and clear any existing errors

  // allocate device memory
  cudaMalloc((void**) &ddata, memsz);
  CUERR // check and clear any existing errors

  // create timer
  timer=wkf_timer_create();

  //
  // Host to device timings
  //

  // non-pinned bandwidth
  wkf_timer_start(timer);
  for (i=0; i<BWITER; i++) {
    cudaMemcpy(ddata, hdata, memsz,  cudaMemcpyHostToDevice);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *hdmbsec = ((double) BWITER) * ((double) memsz) / runtime / (1024.0 * 1024.0);

  // non-pinned latency
  wkf_timer_start(timer);
  for (i=0; i<LATENCYITER; i++) {
    cudaMemcpy(ddata, hdata, 1,  cudaMemcpyHostToDevice);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *hdlatusec = runtime * 1.0e6 / ((double) LATENCYITER);


  // pinned bandwidth
  wkf_timer_start(timer);
  for (i=0; i<BWITER; i++) {
    cudaMemcpy(ddata, phdata, memsz,  cudaMemcpyHostToDevice);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *phdmbsec = ((double) BWITER) * ((double) memsz) / runtime / (1024.0 * 1024.0);

  // pinned latency
  wkf_timer_start(timer);
  for (i=0; i<LATENCYITER; i++) {
    cudaMemcpy(ddata, phdata, 1,  cudaMemcpyHostToDevice);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *phdlatusec = runtime * 1.0e6 / ((double) LATENCYITER);

 
  //
  // Device to host timings
  //

  // non-pinned bandwidth
  wkf_timer_start(timer);
  for (i=0; i<BWITER; i++) {
    cudaMemcpy(hdata, ddata, memsz,  cudaMemcpyDeviceToHost);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *dhmbsec = ((double) BWITER) * ((double) memsz) / runtime / (1024.0 * 1024.0);

  // non-pinned latency
  wkf_timer_start(timer);
  for (i=0; i<LATENCYITER; i++) {
    cudaMemcpy(hdata, ddata, 1,  cudaMemcpyDeviceToHost);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *dhlatusec = runtime * 1.0e6 / ((double) LATENCYITER);


  // pinned bandwidth
  wkf_timer_start(timer);
  for (i=0; i<BWITER; i++) {
    cudaMemcpy(phdata, ddata, memsz,  cudaMemcpyDeviceToHost);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *pdhmbsec = ((double) BWITER) * ((double) memsz) / runtime / (1024.0 * 1024.0);

  // pinned latency
  wkf_timer_start(timer);
  for (i=0; i<LATENCYITER; i++) {
    cudaMemcpy(phdata, ddata, 1,  cudaMemcpyDeviceToHost);
  }
  wkf_timer_stop(timer);
  CUERR // check and clear any existing errors
  runtime = wkf_timer_time(timer);
  *pdhlatusec = runtime * 1.0e6 / ((double) LATENCYITER);
 
 
  cudaFree(ddata);
  CUERR // check and clear any existing errors
  cudaFreeHost(phdata);
  CUERR // check and clear any existing errors
  free(hdata);

  wkf_timer_destroy(timer);

  return 0;
}

typedef struct {
  int deviceid;
  double hdmbsec;
  double hdlatusec;
  double phdmbsec;
  double phdlatusec;
  double dhmbsec;
  double dhlatusec;
  double pdhmbsec;
  double pdhlatusec;
} busbwthrparms;

static void * cudabusbwthread(void *voidparms) {
  busbwthrparms *parms = (busbwthrparms *) voidparms;
  cudabusbw(parms->deviceid, 
            &parms->hdmbsec, &parms->hdlatusec,
            &parms->phdmbsec, &parms->phdlatusec,
            &parms->dhmbsec, &parms->dhlatusec,
            &parms->pdhmbsec, &parms->pdhlatusec);
  return NULL;
}

int vmd_cuda_bus_bw(int numdevs, int *devlist, 
                    double *hdmbsec, double *hdlatusec,
                    double *phdmbsec,double *phdlatusec,
                    double *dhmbsec, double *dhlatusec,
                    double *pdhmbsec, double *pdhlatusec) {
  busbwthrparms *parms;
  wkf_thread_t * threads;
  int i;

  /* allocate array of threads */
  threads = (wkf_thread_t *) calloc(numdevs * sizeof(wkf_thread_t), 1);

  /* allocate and initialize array of thread parameters */
  parms = (busbwthrparms *) malloc(numdevs * sizeof(busbwthrparms));
  for (i=0; i<numdevs; i++) {
    if (devlist != NULL)
      parms[i].deviceid = devlist[i];
    else
      parms[i].deviceid = i;
    parms[i].hdmbsec = 0.0;
    parms[i].hdlatusec = 0.0;
    parms[i].phdmbsec = 0.0;
    parms[i].phdlatusec = 0.0;
    parms[i].dhmbsec = 0.0;
    parms[i].dhlatusec = 0.0;
    parms[i].pdhmbsec = 0.0;
    parms[i].pdhlatusec = 0.0;
  }

#if defined(VMDTHREADS)
  /* spawn child threads to do the work */
  /* thread 0 must also be processed this way otherwise    */
  /* we'll permanently bind the main thread to some device */
  for (i=0; i<numdevs; i++) {
    wkf_thread_create(&threads[i], cudabusbwthread, &parms[i]);
  }

  /* join the threads after work is done */
  for (i=0; i<numdevs; i++) {
    wkf_thread_join(threads[i], NULL);
  }
#else
  /* single thread does all of the work */
  cudabusbwthread((void *) &parms[0]);
#endif

  for (i=0; i<numdevs; i++) {
    hdmbsec[i] = parms[i].hdmbsec; 
    hdlatusec[i] = parms[i].hdlatusec; 
    phdmbsec[i] = parms[i].phdmbsec; 
    phdlatusec[i] = parms[i].phdlatusec; 
    dhmbsec[i] = parms[i].dhmbsec; 
    dhlatusec[i] = parms[i].dhlatusec; 
    pdhmbsec[i] = parms[i].pdhmbsec; 
    pdhlatusec[i] = parms[i].pdhlatusec; 
  }

  /* free thread parms */
  free(parms);
  free(threads);

  return 0;
}



//
// GPU device global memory bandwidth benchmark
//
template <class T>
__global__ void gpuglobmemcpybw(T *dest, const T *src) {
  const unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
  dest[idx] = src[idx];
}

template <class T>
__global__ void gpuglobmemsetbw(T *dest, const T val) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  dest[idx] = val;
}

typedef float4 datatype;

static int cudaglobmembw(int cudadev, double *gpumemsetgbsec, double *gpumemcpygbsec) {
  int i;
  int len = 1 << 22; // one thread per data element
  int loops = 500;
  datatype *src, *dest;
  datatype val=make_float4(1.0f, 1.0f, 1.0f, 1.0f);

  // initialize to zero for starters
  float memsettime = 0.0f;
  float memcpytime = 0.0f;
  *gpumemsetgbsec = 0.0;
  *gpumemcpygbsec = 0.0;

  // attach to the selected device
  cudaError_t rc;
  rc = cudaSetDevice(cudadev);
  if (rc != cudaSuccess) {
#if CUDART_VERSION >= 2010
    rc = cudaGetLastError(); // query last error and reset error state
    if (rc != cudaErrorSetOnActiveProcess)
      return -1; // abort and return an error
#else
    cudaGetLastError(); // just ignore and reset error state, since older CUDA
                        // revs don't have a cudaErrorSetOnActiveProcess enum
#endif
  }

  cudaMalloc((void **) &src, sizeof(datatype)*len);
  CUERR
  cudaMalloc((void **) &dest, sizeof(datatype)*len);
  CUERR

  dim3 BSz(256, 1, 1);
  dim3 GSz(len / (BSz.x * BSz.y * BSz.z), 1, 1); 

  // do a warm-up pass
  gpuglobmemsetbw<datatype><<< GSz, BSz >>>(src, val);
  CUERR
  gpuglobmemsetbw<datatype><<< GSz, BSz >>>(dest, val);
  CUERR
  gpuglobmemcpybw<datatype><<< GSz, BSz >>>(dest, src);
  CUERR

  cudaEvent_t start, end;
  cudaEventCreate(&start);
  cudaEventCreate(&end);

  // execute the memset kernel
  cudaEventRecord(start, 0);
  for (i=0; i<loops; i++) {
    gpuglobmemsetbw<datatype><<< GSz, BSz >>>(dest, val);
  }
  CUERR
  cudaEventRecord(end, 0);
  CUERR
  cudaEventSynchronize(start);
  CUERR
  cudaEventSynchronize(end);
  CUERR
  cudaEventElapsedTime(&memsettime, start, end);
  CUERR

  // execute the memcpy kernel
  cudaEventRecord(start, 0);
  for (i=0; i<loops; i++) {
    gpuglobmemcpybw<datatype><<< GSz, BSz >>>(dest, src);
  }
  cudaEventRecord(end, 0);
  CUERR
  cudaEventSynchronize(start);
  CUERR
  cudaEventSynchronize(end);
  CUERR
  cudaEventElapsedTime(&memcpytime, start, end);
  CUERR

  cudaEventDestroy(start);
  CUERR
  cudaEventDestroy(end);
  CUERR

  *gpumemsetgbsec = (len * sizeof(datatype) / (1024.0 * 1024.0)) / (memsettime / loops);
  *gpumemcpygbsec = (2 * len * sizeof(datatype) / (1024.0 * 1024.0)) / (memcpytime / loops);
  cudaFree(dest);
  cudaFree(src);
  CUERR

  return 0;
}

typedef struct {
  int deviceid;
  double memsetgbsec;
  double memcpygbsec;
} globmembwthrparms;

static void * cudaglobmembwthread(void *voidparms) {
  globmembwthrparms *parms = (globmembwthrparms *) voidparms;
  cudaglobmembw(parms->deviceid, &parms->memsetgbsec, &parms->memcpygbsec);
  return NULL;
}

int vmd_cuda_globmem_bw(int numdevs, int *devlist, 
                        double *memsetgbsec, double *memcpygbsec) {
  globmembwthrparms *parms;
  wkf_thread_t * threads;
  int i;

  /* allocate array of threads */
  threads = (wkf_thread_t *) calloc(numdevs * sizeof(wkf_thread_t), 1);

  /* allocate and initialize array of thread parameters */
  parms = (globmembwthrparms *) malloc(numdevs * sizeof(globmembwthrparms));
  for (i=0; i<numdevs; i++) {
    if (devlist != NULL)
      parms[i].deviceid = devlist[i];
    else
      parms[i].deviceid = i;
    parms[i].memsetgbsec = 0.0;
    parms[i].memcpygbsec = 0.0;
  }

#if defined(VMDTHREADS)
  /* spawn child threads to do the work */
  /* thread 0 must also be processed this way otherwise    */
  /* we'll permanently bind the main thread to some device */
  for (i=0; i<numdevs; i++) {
    wkf_thread_create(&threads[i], cudaglobmembwthread, &parms[i]);
  }

  /* join the threads after work is done */
  for (i=0; i<numdevs; i++) {
    wkf_thread_join(threads[i], NULL);
  }
#else
  /* single thread does all of the work */
  cudaglobmembwthread((void *) &parms[0]);
#endif

  for (i=0; i<numdevs; i++) {
    memsetgbsec[i] = parms[i].memsetgbsec;
    memcpygbsec[i] = parms[i].memcpygbsec;
  }

  /* free thread parms */
  free(parms);
  free(threads);

  return 0;
}


//
// Benchmark latency for complete threadpool barrier wakeup/run/sleep cycle
//
static void * vmddevpoollatencythread(void *voidparms) {
  return NULL;
}

static void * vmddevpooltilelatencythread(void *voidparms) {
  int threadid=-1;
  int tilesize=1;
  void *parms=NULL;
  wkf_threadpool_worker_getid(voidparms, &threadid, NULL);
  wkf_threadpool_worker_getdata(voidparms, (void **) &parms);

  // grind through task tiles until none are left
  wkf_tasktile_t tile;
  while (wkf_threadpool_next_tile(voidparms, tilesize, &tile) != WKF_SCHED_DONE) {
    // do nothing but eat work units...
  }

  return NULL;
}


// no-op kernel for timing kernel launches
__global__ static void nopkernel(float * ddata) {
  unsigned int xindex  = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int yindex  = blockIdx.y * blockDim.y + threadIdx.y;
  unsigned int outaddr = gridDim.x * blockDim.x * yindex + xindex;

  if (ddata != NULL)
    ddata[outaddr] = outaddr;
}

// empty kernel for timing kernel launches
__global__ static void voidkernel(void) {
  return;
}

static void * vmddevpoolcudatilelatencythread(void *voidparms) {
  int threadid=-1;
  int tilesize=1;
  float *parms=NULL;
  wkf_threadpool_worker_getid(voidparms, &threadid, NULL);

  // XXX Note that we expect parms to be set to NULL or a valid CUDA
  //     global memory pointer for correct operation of the NOP kernel below
  wkf_threadpool_worker_getdata(voidparms, (void **) &parms);

#if 0
  // scale tile size by device performance
  tilesize=4; // GTX 280, Tesla C1060 starting point tile size
  wkf_threadpool_worker_devscaletile(voidparms, &tilesize);
#endif

  // grind through task tiles until none are left
  wkf_tasktile_t tile;
  dim3 Gsz(1,1,0);
  dim3 Bsz(8,8,1);
  while (wkf_threadpool_next_tile(voidparms, tilesize, &tile) != WKF_SCHED_DONE) {
    // launch a no-op CUDA kernel
    nopkernel<<<Gsz, Bsz, 0>>>(parms);
  }

  // wait for all GPU kernels to complete
  cudaDeviceSynchronize();

  return NULL;
}


int vmd_cuda_devpool_latency(wkf_threadpool_t *devpool, int tilesize,
                             double *kernlaunchlatency,
                             double *barlatency,
                             double *cyclelatency, 
                             double *tilelatency,
                             double *kernellatency) {
  int i;
  wkf_tasktile_t tile;
  wkf_timerhandle timer;
  int loopcount;

  timer=wkf_timer_create();

  // execute just a CUDA kernel launch and measure latency on whatever
  // GPU we get.
  loopcount = 15000;
  dim3 VGsz(1,1,0);
  dim3 VBsz(8,8,1);
  wkf_timer_start(timer);
  for (i=0; i<loopcount; i++) {
    voidkernel<<<VGsz, VBsz, 0>>>();
  }
  // wait for GPU kernels to complete
  cudaDeviceSynchronize();
  wkf_timer_stop(timer);
  *kernlaunchlatency = wkf_timer_time(timer) / ((double) loopcount);

  // execute just a raw barrier sync and measure latency
  loopcount = 15000;
  wkf_timer_start(timer);
  for (i=0; i<loopcount; i++) {
    wkf_threadpool_wait(devpool);
  }
  wkf_timer_stop(timer);
  *barlatency = wkf_timer_time(timer) / ((double) loopcount);

  // time wake-up, launch, and sleep/join of device pool doing a no-op
  loopcount = 5000;
  wkf_timer_start(timer);
  for (i=0; i<loopcount; i++) {
    tile.start=0;
    tile.end=0;
    wkf_threadpool_sched_dynamic(devpool, &tile);
    wkf_threadpool_launch(devpool, vmddevpoollatencythread, NULL, 1);
  }
  wkf_timer_stop(timer);
  *cyclelatency = wkf_timer_time(timer) / ((double) loopcount);

  // time wake-up, launch, and sleep/join of device pool eating tiles
  loopcount = 5000;
  wkf_timer_start(timer);
  for (i=0; i<loopcount; i++) {
    tile.start=0;
    tile.end=tilesize;
    wkf_threadpool_sched_dynamic(devpool, &tile);
    wkf_threadpool_launch(devpool, vmddevpooltilelatencythread, NULL, 1);
  }
  wkf_timer_stop(timer);
  *tilelatency = wkf_timer_time(timer) / ((double) loopcount);

  // time wake-up, launch, and sleep/join of device pool eating tiles
  loopcount = 2000;
  wkf_timer_start(timer);
  for (i=0; i<loopcount; i++) {
    tile.start=0;
    tile.end=tilesize;
    wkf_threadpool_sched_dynamic(devpool, &tile);
    wkf_threadpool_launch(devpool, vmddevpoolcudatilelatencythread, NULL, 1);
  }
  wkf_timer_stop(timer);
  *kernellatency = wkf_timer_time(timer) / ((double) loopcount);

  wkf_timer_destroy(timer);

#if 1
  vmd_cuda_measure_latencies(devpool);
#endif

  return 0;
}


//
// Benchmark CUDA kernel launch and memory copy latencies in isolation
//
typedef struct {
  int deviceid;
  int testloops;
  double kernlatency;
  double bcopylatency;
  double kbseqlatency;
} latthrparms;

static void * vmddevpoolcudalatencythread(void *voidparms) {
  int threadid=-1;
  latthrparms *parms=NULL;

  wkf_threadpool_worker_getid(voidparms, &threadid, NULL);
  wkf_threadpool_worker_getdata(voidparms, (void **) &parms);
  if (parms->deviceid == threadid) { 
    wkf_timerhandle timer;
    timer=wkf_timer_create();
    printf("Thread/device %d running...\n", threadid);
    cudaStream_t devstream;
    cudaStreamCreate(&devstream);

    char *hostbuf = (char *) calloc(1, 65536 * sizeof(char));
    char  *gpubuf = NULL;
    cudaMalloc((void**)&gpubuf, 65536 * sizeof(char));

    dim3 Gsz(1,1,0);
    dim3 Bsz(8,8,1);

    // measure back-to-back NULL kernel launches
    wkf_timer_start(timer);
    int i;
    for (i=0; i<parms->testloops; i++) {
      // launch a no-op CUDA kernel
      nopkernel<<<Gsz, Bsz, 0, devstream>>>(NULL);
    }
    // wait for all GPU kernels to complete
    cudaStreamSynchronize(devstream);
    wkf_timer_stop(timer);
    parms->kernlatency =  1000000 * wkf_timer_time(timer) / ((double) parms->testloops);

    // measure back-to-back round-trip 1-byte memcpy latencies
    wkf_timer_start(timer);
    for (i=0; i<parms->testloops; i++) {
      cudaMemcpyAsync(gpubuf, hostbuf, 1, cudaMemcpyHostToDevice, devstream);
      cudaMemcpyAsync(hostbuf, gpubuf, 1, cudaMemcpyDeviceToHost, devstream);
    }
    // wait for all GPU kernels to complete
    cudaStreamSynchronize(devstream);
    wkf_timer_stop(timer);
    parms->kernlatency =  1000000 * wkf_timer_time(timer) / ((double) parms->testloops);

    printf("NULL kernel launch latency (usec): %.2f\n", parms->kernlatency);

    cudaStreamDestroy(devstream);
    cudaFree(gpubuf);
    free(hostbuf);
    wkf_timer_destroy(timer);
  }

  return NULL;
}


int vmd_cuda_measure_latencies(wkf_threadpool_t *devpool) {
  latthrparms thrparms;
  int workers = wkf_threadpool_get_workercount(devpool);
  int i;
printf("vmd_cuda_measure_latencies()...\n");
  for (i=0; i<workers; i++) {
    memset(&thrparms, 0, sizeof(thrparms));
    thrparms.deviceid = i;
    thrparms.testloops = 2500;
    wkf_threadpool_launch(devpool, vmddevpoolcudalatencythread, &thrparms, 1);
  }

  return 0;
}


#if defined(VMDUSECUDAGDS)
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
} gpuoocbenchthreadparms;

#define VMDGDSMAXFRAMEBUF     8
static void * gpu_ooc_bench_thread(void *voidparms) {
  int threadid, numthreads;
  gpuoocbenchthreadparms *parms = NULL;
  wkf_threadpool_worker_getdata(voidparms, (void **) &parms);
  wkf_threadpool_worker_getid(voidparms, &threadid, &numthreads);

  //
  // copy in per-thread parameters
  //
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
    printf("Thr[%2d] %d frames total, natoms: %d selected: %d\n",
           threadid, framecount, natoms, sel->selected);
    printf("Thr[%2d] %d frames/file\n", threadid, framesperfile);
  }

  cudaError_t crc;
  cudaStream_t oocstream;
  float *devptr=NULL;
  float *hostptr=NULL;
  float *hostptr_unaligned=NULL;

  float *crdx1=NULL, *crdy1=NULL, *crdz1=NULL;
  float *crdx2=NULL, *crdy2=NULL, *crdz2=NULL;
  int multiframeio = 0;
  if (getenv("VMDGDSMULTIFRAME"))
    multiframeio = atoi(getenv("VMDGDSMULTIFRAME"));
  if (multiframeio > VMDGDSMAXFRAMEBUF)
    multiframeio = VMDGDSMAXFRAMEBUF;

  // set block sizes and counts for IO bench calcs
  dim3 IOBsz = dim3(256, 1, 1);
  dim3 IOGsz = dim3((natoms + IOBsz.x - 1) / IOBsz.x, 1, 1);

  if (parms->devcount > 0) {
    long gpuallocsz = (VMDGDSMAXFRAMEBUF+1) * blockpadsz;

    if (threadid == 0) {
      printf("Thr[%2d] Allocating GPU timestep I/O buf: %ld \n",
             threadid, gpuallocsz);
    }
    crc = cudaMalloc((void**) &devptr, gpuallocsz);

    if (hostfds != NULL) {
      if (pinhostiobuffer) {
        crc = cudaMallocHost((void**) &hostptr, gpuallocsz);
      } else {
        hostptr = (float *) alloc_aligned_ptr(gpuallocsz, 4096,
                                              (void**) &hostptr_unaligned);
        if (!hostptr) {
          printf("Thr[%d]: Failed allocation!\n", threadid);
          return NULL;
        }
      }
    }

    long crdsz = sel->selected * sizeof(float);

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

    cudaStreamCreate(&oocstream);

#if defined(VMDUSECUDAGDS)
    cuFileBufRegister(devptr, gpuallocsz, 0);
#endif
  }

  int verbose = (getenv("VMDGDSVERBOSE") != NULL) ? 1 : 0;
  
  int filestrategy = 0;
  if (getenv("VMDGDSFILESTRATEGY")) {
    filestrategy = atoi(getenv("VMDGDSFILESTRATEGY"));
  }
  if (threadid == 0) {
    printf("Thr[%2d] file strategy set to: %d\n", threadid, filestrategy);
  }

  wkf_tasktile_t tile;
  while (wkf_threadlaunch_next_tile(voidparms, VMDGDSMAXFRAMEBUF * 1, &tile) != WKF_SCHED_DONE) {
    //
    // simple I/O + compute benchmarking...
    //
    int idx;
    int threadspergpu;
    if (parms->devcount > 0)
      threadspergpu = numthreads / parms->devcount;
    else 
      threadspergpu = 1;

    for (idx=tile.start; idx<tile.end; idx++) {
      int myfileidx, fileframeidx;

      switch (filestrategy) {
        case 1:
          myfileidx = (idx / multiframeio) % nfiles;
          fileframeidx = idx % framesperfile;
          break;

        case 2:
          myfileidx = (idx / (multiframeio * threadspergpu)) % nfiles;
          fileframeidx = idx % framesperfile;
          break;

        case 3:
          myfileidx = (threadid / 4) % nfiles;
          fileframeidx = idx % framesperfile;
          break;

        case 0:
        default:
          myfileidx = (threadid / threadspergpu) % nfiles;
          fileframeidx = idx % framesperfile;
          break;
      }

      //
      // compute multi-frame or single-frame I/O offsets and sizes
      //
      long startoffset, foffset, readlen;
      read_js_timestep_index_offsets(parms->jshandles[myfileidx],
                                     natoms, fileframeidx, 0, natoms, NULL,
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
      } else {
        printf("Inconsistent cufile/hostfds state, aborting!\n");
        return NULL;
      }

      // handle errors if they have occured
      if (ret < 0) {
        printf("Thr[%2d] Error: cuFileRead(): %ld\n", threadid, ret);
        return NULL; // XXX error handling needs to be done here
      }

      if (verbose) {
        printf("Thr[%2d]F[%d][tile: %d to %d] frame: %d cuFile len: %ld off: %ld\n",
               threadid, myfileidx, tile.start, tile.end, idx, 
               readlen, foffset);
      }
    }
  }

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


int gpu_ooc_bench(wkf_threadpool_t *devpool, // VMD GPU worker thread pool
                  int nfiles, const char **trjfileset, const AtomSel *sel,
                  int first, int last, int step) {
  printf("gpu_ooc_bench()\n");
  wkf_threadpool_t *bigpool = NULL;

#if defined(VMDUSECUDAGDS)
  int devcount;
  cudaError_t crc = cudaGetDeviceCount(&devcount);
  printf("gpu_ooc_bench) GPU device count: %d\n", devcount);
  if (devcount==0)
    printf("gpu_ooc_bench)   No GPU devices, continuing with host only...\n");

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
    printf("gpu_ooc_bench) File[%d] GDS setup, opening '%s'\n", i, filename);
    jshandles[i] = (jshandle *) open_js_read(filename, "js", &natomschk);
    if (!jshandles[i]) {
      printf("gpu_ooc_bench) File[%d] open_js_read failed for file %s\n", i, filename);
      return -1; // deal with error handling later
    }

#if vmdplugin_ABIVERSION > 17
    long blocksz = MOLFILE_DIRECTIO_MIN_BLOCK_SIZE;
    int filepgalignsz = 1;
    read_js_timestep_pagealign_size(jshandles[i], &filepgalignsz);
    if (filepgalignsz != blocksz) {
      printf("gpu_ooc_bench) File[%d] Plugin-returned page alignment size mismatch!\n", i);
    } else {
      printf("gpu_ooc_bench) File[%d] Page alignment size: %d\n", i, filepgalignsz);
    }
#endif

    read_js_timestep_index_offsets(jshandles[i], natomschk, 0, 0, 0,
                                   &directio_fds[i], NULL, NULL, NULL);

    cfhdesc[i].handle.fd = directio_fds[i]; // typedef of Unix FD
    cfhdesc[i].type = CU_FILE_HANDLE_TYPE_OPAQUE_FD;
    CUfileError_t cferr = cuFileHandleRegister(&cfh[i], &cfhdesc[i]);

    if (cferr.err != CU_FILE_SUCCESS) {
      printf("gpu_ooc_bench) File[%d] cuFileImportExternalFile on fd %d failed!\n",
             i, cfhdesc[i].handle.fd);
      return -1; // XXX error handling needs to be done here
    }
  }


  //
  // copy in per-thread parameters
  //
  gpuoocbenchthreadparms parms;
  memset(&parms, 0, sizeof(parms));
  parms.devcount = devcount;
  parms.nfiles = nfiles;
  parms.trjfileset = trjfileset;
  parms.jshandles = jshandles;
  parms.cfh = cfh;
  parms.natoms = sel->num_atoms;
  parms.sel = sel;
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
  tile.end=framecount - 1; // first row only

  printf("gpu_ooc_bench) tile start: %d  end: %d\n", tile.start, tile.end);

  int gdsthreadspergpu = 1;
  if (getenv("VMDGDSTHREADSPERGPU") != NULL)
    gdsthreadspergpu = atoi(getenv("VMDGDSTHREADSPERGPU"));

  printf("gpu_ooc_bench) gdsthreadspergpu: %d\n", gdsthreadspergpu);

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
  wkf_threadpool_launch(devpool, gpu_ooc_bench_thread, &parms, 1);
  wkf_timer_stop(timer);

  double runtime = wkf_timer_time(timer);
  double gbytes = sel->num_atoms * 12L * (tile.end+1) / (1024.0 * 1024.0 * 1024.0);

  printf("gpu_ooc_bench) natoms: %d, fsz: %ld, tsz: %ld\n",
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
    printf("GDS multi-frame read opt: %d frames per call, %ld bytes\n",
           multiframeio,
           multiframeio * sel->num_atoms * 12L);
  }

  printf("OOC runtime: %.1f, %.2fGB/sec\n", runtime, gbytes/runtime);

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



