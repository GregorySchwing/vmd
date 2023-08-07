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
 *      $RCSfile: CUDAUtil.cu,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.56 $        $Date: 2020/02/26 19:26:47 $
 *
 ***************************************************************************/
/**
 * \file CUDAUtil.cu
 * \brief Lightweight CUDA helper routines for use by CUDAAccel
 */

#include <string.h>
#include <stdio.h>
#include "CUDAKernels.h"
#include "WKFThreads.h"
#include "intstack.h"
#include "ProfileHooks.h"

// to eliminate a bunch of old ifdef tests, we'll just check
// once here for a reasonably recent rev of CUDA.
#if CUDART_VERSION <= 6000
#error VMD requires at least CUDA 6.x or later
#endif

#if defined(__cplusplus)
extern "C" {
#endif

// report true if the driver is compatible with the runtime
static int vmd_cuda_drv_runtime_compatible() {
#if CUDART_VERSION >= 2020
  int cuda_driver_version=-1;
  int cuda_runtime_version=0;

  cudaDriverGetVersion(&cuda_driver_version);
  cudaRuntimeGetVersion(&cuda_runtime_version);

#if 0
  printf("CUDA driver version: %d\n", cuda_driver_version);
  printf("CUDA runtime version: %d\n", cuda_runtime_version);
#endif

  if (cuda_driver_version == 0) 
    return VMDCUDA_ERR_NODEVICES;

  if (cuda_driver_version < cuda_runtime_version) {
#if defined(ARCH_LINUXARM64)
    // XXX workaround for the first native CUDA compiler toolchain (5.5)
    //     having a newer rev than the driver (310.32, CUDA 5.0) reports
    if (cuda_driver_version == 10010 && cuda_runtime_version == 10020)
      return VMDCUDA_ERR_NONE;
#endif
#if defined(ARCH_LINUXCARMA)
    // XXX workaround for the first native CUDA compiler toolchain (5.5)
    //     having a newer rev than the driver (310.32, CUDA 5.0) reports
    if (cuda_driver_version == 5000 && cuda_runtime_version == 5050)
      return VMDCUDA_ERR_NONE;
#endif
    return VMDCUDA_ERR_DRVMISMATCH;
  }
#endif  

  return VMDCUDA_ERR_NONE;
}


int vmd_cuda_device_props(int dev, char *name, int namelen,
                          int *devmajor, int *devminor, 
                          unsigned long *memb, int *clockratekhz,
                          int *smcount, int *integratedgpu,
                          int *asyncenginecount, int *kerneltimeout,
                          int *canmaphostmem, int *computemode,
                          int *spdpfpperfratio, 
                          int *pageablememaccess, 
                          int *pageablememaccessuseshostpagetables) {
  cudaError_t rc;
  cudaDeviceProp deviceProp;

  // this extra paranoia check costs about 0.04 seconds on a DGX-2 w/ 16 GPUs
  int vercheck;
  if ((vercheck = vmd_cuda_drv_runtime_compatible()) != VMDCUDA_ERR_NONE) {
    return vercheck;
  }

  memset(&deviceProp, 0, sizeof(cudaDeviceProp));
  if ((rc=cudaGetDeviceProperties(&deviceProp, dev)) != cudaSuccess) {
    // printf("error: %s\n", cudaGetErrorString(rc));
    if (rc == cudaErrorNotYetImplemented)
      return VMDCUDA_ERR_EMUDEVICE;
    return VMDCUDA_ERR_GENERAL;
  }

  if (name)
    strncpy(name, deviceProp.name, namelen);
  if (devmajor)
    *devmajor = deviceProp.major;
  if (devminor)
    *devminor = deviceProp.minor;
  if (memb)
    *memb = deviceProp.totalGlobalMem;
  if (clockratekhz)
    *clockratekhz = deviceProp.clockRate;
  if (smcount)
    *smcount = deviceProp.multiProcessorCount;
  if (asyncenginecount)
    *asyncenginecount = deviceProp.asyncEngineCount;
  if (kerneltimeout)
    *kerneltimeout = (deviceProp.kernelExecTimeoutEnabled != 0);
  if (integratedgpu)
    *integratedgpu = (deviceProp.integrated != 0);
  if (canmaphostmem)
    *canmaphostmem = (deviceProp.canMapHostMemory != 0);
  if (computemode)
    *computemode = deviceProp.computeMode;
  if (spdpfpperfratio)
    *spdpfpperfratio = deviceProp.singleToDoublePrecisionPerfRatio;
  if (pageablememaccess)
    *pageablememaccess = deviceProp.pageableMemoryAccess;
  if (pageablememaccessuseshostpagetables)
#if CUDA_VERSION >= 10000
    *pageablememaccessuseshostpagetables = deviceProp.pageableMemoryAccessUsesHostPageTables;
#else
    *pageablememaccessuseshostpagetables = 0;
#endif

  return VMDCUDA_ERR_NONE;
}


int vmd_cuda_num_devices(int *numdev) {
  int devcount=0;
  *numdev = 0;

  // When this is the very first CUDA API call during VMD startup, 
  // there's a 2.0 second startup lag associated with it on the DGX-2, 
  // likely due to CUDA runtime library internal initialization overheads
  // across the 16 GPUs.
  int vercheck;
  if ((vercheck = vmd_cuda_drv_runtime_compatible()) != VMDCUDA_ERR_NONE) {
    return vercheck;
  }

  if (cudaGetDeviceCount(&devcount) != cudaSuccess) {
    return VMDCUDA_ERR_NODEVICES;
  }

  // Do a sanity check in case we get complete gibberish back,
  // but no error. This can occur if we have either a driver or 
  // CUDA runtime that's badly mismatched.
  if (devcount > 100 || devcount < 0)
    return VMDCUDA_ERR_DRVMISMATCH;

#if 0
  // XXX we can eliminate this code now since device emulation
  //     went away back in the very early days of CUDA.  This extra
  //     pass over all of the device info just slows down our initialization
  //     on machines like a DGX-2 that have a huge number of GPUs...

  // disregard emulation mode as unusable for our purposes
  int usabledevs=0;
  int i;
  for (i=0; i<devcount; i++) {
    int devmajor, devminor, rc;

    rc = vmd_cuda_device_props(i, NULL, 0, &devmajor, &devminor, NULL, 
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
                               NULL, NULL, NULL);

    if (rc == VMDCUDA_ERR_NONE) {
      // Check for emulation mode devices, and ignore if found
      if (((devmajor >= 1) && (devminor >= 0)) &&
          ((devmajor != 9999) && (devminor != 9999))) {
        usabledevs++;
      }
    } else if (rc != VMDCUDA_ERR_EMUDEVICE) {
      return VMDCUDA_ERR_SOMEDEVICES;
    }
  } 
  *numdev = usabledevs;
#else
  *numdev = devcount;
#endif

  return VMDCUDA_ERR_NONE;
}


int vmd_cuda_peer_matrix(int *numdev, 
                         int **p2pmat, 
                         int **p2psupp,
                         int **p2patomics,
                         int **p2parrays, 
                         int **perfmat,
                         int *p2plinkcount,
                         int *islandcount) {
  *p2plinkcount = 0;
  *numdev = 0;
  *p2pmat = NULL;
  *p2psupp = NULL;
  *p2patomics = NULL;
  *p2parrays = NULL;
  *perfmat = NULL;
  *islandcount = 0;

  int rc = vmd_cuda_num_devices(numdev);
  if (rc == VMDCUDA_ERR_NONE) {
    // check for bogus GPU count values...
    int N = (*numdev);
    if (N < 1 || N > 1024)
      return VMDCUDA_ERR_GENERAL; 

    int i, j;
    int sz = N*N;

    // allocate and compute GPU P2P connectivity/perf/etc matrices
    *p2pmat = (int*) calloc(1, sz * sizeof(int));
    *p2psupp = (int*) calloc(1, sz * sizeof(int));
    *p2patomics = (int*) calloc(1, sz * sizeof(int));
    *p2parrays = (int*) calloc(1, sz * sizeof(int));
    *perfmat = (int*) calloc(1, sz * sizeof(int));
    for (j=0; j<N; j++) {
      for (i=j+1; i<N; i++) {  // symmetric matrix...
        cudaError_t err = cudaSuccess;
        int canAccessPeer = 0;
        int supported = 0;
        int atomics = 0;
        int arrays = 0;
        int perfrank = -1;

#if 0
printf("Checking[%d][%d]", j, i);
#endif
        err = cudaDeviceCanAccessPeer(&canAccessPeer, i, j);
        err = cudaDeviceGetP2PAttribute(&supported, cudaDevP2PAttrAccessSupported, i, j);
        err = cudaDeviceGetP2PAttribute(&atomics, cudaDevP2PAttrNativeAtomicSupported, i, j);
#if CUDA_VERSION > 9000
        err = cudaDeviceGetP2PAttribute(&arrays, cudaDevP2PAttrCudaArrayAccessSupported, i, j);
#endif
        err = cudaDeviceGetP2PAttribute(&perfrank, cudaDevP2PAttrPerformanceRank, i, j);
        if (err != cudaSuccess) {
          free(*p2pmat);
          free(*p2psupp);
          free(*p2patomics);
          free(*p2parrays); 
          free(*perfmat);
          return VMDCUDA_ERR_GENERAL; 
        }

        // symmetric matrices...
        int idx = j*N + i;
        int symidx = i*N + j;
    
        (*p2pmat)[idx] = canAccessPeer;
        (*p2pmat)[symidx] = canAccessPeer; 

        (*p2psupp)[idx] = supported;
        (*p2psupp)[symidx] = supported;

        (*p2patomics)[idx] = atomics;
        (*p2patomics)[symidx] = atomics;

        (*p2parrays)[idx] = arrays;
        (*p2parrays)[symidx] = arrays;

        (*perfmat)[idx] = perfrank;
        (*perfmat)[symidx] = perfrank;

        // count up the total number of P2P links
        if (canAccessPeer && supported)
          (*p2plinkcount)++;
      }
    }

#if 0
    for (j=0; j<N; j++) {
      for (i=0; i<N; i++) {  // symmetric matrix...
        int idx = j*N + i;
        printf("P2P[%2d][%2d]: %d, sup %d, atomics %d, arrays %d, perf %d\n", 
               j, i, (*p2pmat)[idx], (*p2psupp)[idx], (*p2patomics)[idx], 
               (*p2parrays)[idx], (*perfmat)[idx]);
      }
    }
#endif

    // Once we have the GPU P2P connectivity matrix we
    // can compute the number of GPU P2P "islands" via DFS
    int *flags = (int*) calloc(1, N * sizeof(int));
    IntStackHandle s;
    s = intstack_create(N);
    int islands, l;
    for (islands=0,l=0; l<N; l++) {
      // visit GPUs that haven't yet been accounted for
      if (!flags[l]) { 
        intstack_push(s, l);

        // Loop over all GPUs we have P2P links to in the same island
        while (!intstack_pop(s, &j)) {
          flags[j] = 1; // mark as visited 

          // check connectivity with peers
          for (i=0; i<N; i++) {
            // find GPUs connected to index i
            int idx = j*N + i;
            if ((*p2pmat)[idx]==1 && (*p2psupp)[idx]==1) {
              // push unvisited GPUs onto the stack
              if (flags[i] == 0) 
                intstack_push(s, i);
            }
          }
        }

        islands++;
      } 
    } 
    intstack_destroy(s);
    free(flags);
    *islandcount = islands;
  }
  
  return rc;
}


void * vmd_cuda_devpool_setdeviceonly(void * voidparms) {
  int count, id, dev;
  char *mesg;

  wkf_threadpool_worker_getid(voidparms, &id, &count);
  wkf_threadpool_worker_getdata(voidparms, (void **) &mesg);
  wkf_threadpool_worker_getdevid(voidparms, &dev);

  // mark CPU-GPU management threads for display in profiling tools
  char threadname[1024];
  sprintf(threadname, "VMD GPU threadpool[%d]", id);
  PROFILE_NAME_THREAD(threadname);

  /* set active device */
  cudaSetDevice(dev);

  if (mesg != NULL)
    printf("devpool thread[%d / %d], device %d message: '%s'\n", id, count, dev, mesg);

  return NULL;
}


// Permanently bind worker threads to CPUs to achieve 
// peak performance and reducing run-to-run jitter for
// graphics related algorithms.  This code simply performs
// a one-to-one mapping of threads to CPUs at present, and
// doesn't handle other cases with any greater logic yet.
void * vmd_cuda_affinitize_threads(void *voidparms) {
  int tid, numthreads;
  wkf_threadpool_worker_getid(voidparms, &tid, &numthreads);
  int physcpus = wkf_thread_numphysprocessors();
  int setthreadaffinity=1; // enable affinity assignment by default

  int devcount;
  cudaError_t crc = cudaGetDeviceCount(&devcount);

  int devid = 0;
  wkf_threadpool_worker_getdevid(voidparms, &devid);

  int cpuspergpu = (physcpus/2) / devcount; // XXX HT/SMT hack
  int workerthreadspergpu = numthreads / devcount;
  
  int verbose = (getenv("VMDCUDACPUAFFINITYVERBOSE") != NULL);

  // XXX need to add some kind of check to skip trying to set affinity
  //     if the OS, batch system, or other external CPU set management
  //     system has already restricted or set of CPUs we can use

  if ((physcpus > 0) && setthreadaffinity) {
    int taffinity=0;

#if 1
    // On Intel/AMD hardware, physically distinct CPU cores are numbered 
    // consecutively, and additional SMT "cores" appear as additional multiples
    // of the physical CPU core count.  It is therefore sufficient to 
    // use consecutive CPU indices to spread worker threads fully across CPUs
    // and SMT hardware threads
    taffinity=(devid * cpuspergpu) + (tid % workerthreadspergpu); // set affinity of tid to same CPU id if numthreads == cpus
#endif

    if (verbose)
      printf("Affinitizing GPU[%d] worker thread[%d] to CPU[%d]...\n", 
             devid, tid, taffinity);

    wkf_thread_set_self_cpuaffinity(taffinity);
  }

  // mark CPU threads for display in profiling tools
  char threadname[1024];
  sprintf(threadname, "VMD GPU threadpool[%d]", tid);
  PROFILE_NAME_THREAD(threadname);

  return NULL;
}



/* XXX hard-coded device indexing scheme, needs to use abstractions */
void * vmd_cuda_devpool_enable_P2P(void * voidparms) {
  int threadid;
  wkf_threadpool_worker_getid(voidparms, &threadid, NULL);

  int i;
  if (threadid == 0)
    printf("Enabling all-to-all GPU Peer Access for DGX-2\n");

  int devcount;
  cudaError_t crc = cudaGetDeviceCount(&devcount);

  int devid = 0;
  wkf_threadpool_worker_getdevid(voidparms, &devid);

  for (i=0; i<devcount; i++) {
    if (i != devid) {
      cudaError_t cuerr = cudaDeviceEnablePeerAccess(i, 0);
      if (cuerr != cudaSuccess) {
        printf("Thr[%2d] Error enabling peer access to GPU %d\n", threadid, i);
      }
    }
  }

  return NULL;
}


void * vmd_cuda_devpool_setdevice(void * voidparms) {
  int count, id, dev;
  cudaDeviceProp deviceProp;
  char *mesg;
  char *d_membuf;
  cudaError_t err;

  wkf_threadpool_worker_getid(voidparms, &id, &count);
  wkf_threadpool_worker_getdata(voidparms, (void **) &mesg);
  wkf_threadpool_worker_getdevid(voidparms, &dev);

  // XXX This point within the CPU thread GPU device association 
  //     process is the right point at which to enforce CPU thread
  //     affinity or affinity masks that are based on the host platform's
  //     NUMA, PCIe, and NVLink interconnect topology.  We should assign
  //     a CPU thread affinity or affinity mask based on a combination of:
  //
  //     A) Incoming thread affinity masks set by the job scheduler and
  //        subsequently enforced by the OS (e.g. IBM 'jsrun' on ORNL Summit)
  //
  //     B) The CPU affinity mask associated with each GPU device as
  //        reported by the NVML API.
  //
  //     C) An indexing scheme for the target CPU's logical cores
  //        that distributes GPU host threads among CPU cores, while
  //        avoiding co-scheduling multiple GPU management host threads
  //        on the same physical CPU cores.  This is an extra portability
  //        issue even among different CPU architectures on Linux, as the CPU 
  //        SMT logical core indexing scheme as used by IBM and by Intel are
  //        completely incompatible with each other:  
  //        C1) Intel lists CPU physical cores consecutively, then followed 
  //            by their hyperthread logical cores  So, to use only physical
  //            CPU cores, where N is the physical core count, one would use
  //            logical core indices 0 to N-1.
  //        C2) IBM lists POWER CPU physical cores followed consecutively by
  //            their associated logical SMT cores, then followed by the
  //            next physical cores, etc.  So, to use only physical
  //            CPU cores, where N is the physical core count, and S is the
  //            max SMT depth of each CPU, one would use logical core indices
  //            0, (S), (2*S), (3S), ... ((N-1)*S).

  // mark CPU-GPU management threads for display in profiling tools
  char threadname[1024];
  sprintf(threadname, "VMD GPU threadpool[%d]", id);
  PROFILE_NAME_THREAD(threadname);

  /* set active device */
  cudaSetDevice(dev);

  /* Query SM count and clock rate, and compute a speed scaling value */
  /* the current code uses a GeForce GTX 280 / Tesla C1060 as the     */
  /* "1.0" reference value, with 30 SMs, and a 1.3 GHz clock rate     */ 
  memset(&deviceProp, 0, sizeof(cudaDeviceProp));
  if (cudaGetDeviceProperties(&deviceProp, dev) == cudaSuccess) {
    float smscale = ((float) deviceProp.multiProcessorCount) / 30.0f;
    double clockscale = ((double) deviceProp.clockRate) / 1295000.0;
    float speedscale = smscale * ((float) clockscale);

#if 0
    printf("clock rate: %lf\n", (double) deviceProp.clockRate);
    printf("scale: %.4f smscale: %.4f clockscale: %.4f\n", 
           speedscale, smscale, clockscale);  
#endif

    if (deviceProp.canMapHostMemory != 0) {
#if 0
      printf("Enabled mapped host memory on device[%d]\n", dev);
#endif

      /* 
       * set blocking/yielding API behavior and enable mapped host memory
       * If this fails, then either we've got a problematic device, or 
       * we have already set the device flags within this thread (shouldn't
       * ever happen), or the device we're accessing doesn't actually support
       * mapped host memory (shouldn't ever happen since we check for that).
       */

#if defined(VMDLIBOPTIX)
      // when compiled with OptiX enabled, we tell the CUDA runtime to 
      // maintain the peak local memory size that occured at runtime
      // to avoid thrashing with difficult scenes
      err = cudaSetDeviceFlags(cudaDeviceScheduleAuto | cudaDeviceMapHost | cudaDeviceLmemResizeToMax);
#else
      err = cudaSetDeviceFlags(cudaDeviceScheduleAuto | cudaDeviceMapHost);
#endif
      if (err != cudaSuccess) {
        printf("Warning) thread[%d] can't set GPU[%d] device flags\n", id, dev);
        printf("Warning) CUDA error: %s\n", cudaGetErrorString(err)); 
      }
    }

    wkf_threadpool_worker_setdevspeed(voidparms, speedscale);

    /* 
     * Do a small 1MB device memory allocation to ensure that our context
     * has actually been initialized by the time we return.
     * If this tiny allocation fails, then something is seriously wrong
     * and we should mark this device as unusable for the rest of 
     * this VMD session.
     */
    if ((err = cudaMalloc((void **) &d_membuf, 1*1024*1024)) == cudaSuccess) {
      cudaFree(d_membuf); 
    } else {
      printf("Warning) thread[%d] can't init GPU[%d] found by device query\n", id, dev); 
      printf("Warning) CUDA error: %s\n", cudaGetErrorString(err));
      /* 
       * XXX we should mark the device unusable here so that no other code
       *     touchies it, but have no mechanism for doing that yet...
       */
    }
  }

  if (mesg != NULL)
    printf("devpool thread[%d / %d], device %d message: '%s'\n", id, count, dev, mesg);

  return NULL;
}



#if defined(__cplusplus)
}
#endif


