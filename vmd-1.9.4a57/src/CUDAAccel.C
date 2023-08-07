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
 *	$RCSfile: CUDAAccel.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.70 $	$Date: 2022/02/13 05:34:21 $
 *
 ***************************************************************************/
/**
 * \file CUDAAccel.C
 * \brief VMD CUDA GPU and topology enumeration and management interface.
 *
 * Class to enumerate and initialize CUDA GPU accelerator devices, while 
 * providing an abstraction layer that allows a user-defined subset of
 * GPU devices to be enabled, while excluding others based on their 
 * hardware capabilities, software configuration properties, or to facilitate
 * conflict free sharing of GPUs among multiple application processes 
 * on the same workstation or compute node.  In order for the abstraction to
 * work properly, all of the CUDA code in VMD needs to make use of the 
 * abstracted query functions rather than the native CUDA functions, so
 * that device masks and selection criteria are always honored, and to 
 * avoid needless duplication of GPU device selection logic, and the like. 
 */

#include <stdio.h>
#include <stdlib.h>
#include "config.h"     // rebuild on config changes
#include "Inform.h"
#include "ResizeArray.h"
#include "CUDAAccel.h"
#include "CUDAKernels.h"
#include "WKFThreads.h"
#include "ProfileHooks.h"

CUDAAccel::CUDAAccel(void) {
  cudaavail = 0;
  numdevices = 0;
  numphysdevices = 0;

  nvmlh=NULL;
  cudapool=NULL;

  if (getenv("VMDNOCUDA") != NULL) {
    msgInfo << "VMDNOCUDA environment variable is set, CUDA support disabled."
            << sendmsg;
    return; 
  }

#if defined(VMDCUDA)
  PROFILE_PUSH_RANGE("CUDAAccel::CUDAAccel()", 0);

  unsigned int gpumask = 0xffffffff;
  const char *gpumaskstr = getenv("VMDCUDADEVICEMASK");
  if (gpumaskstr != NULL) {
    unsigned int tmp;
    if (sscanf(gpumaskstr, "%x", &tmp) == 1) {
      gpumask = tmp;
      msgInfo << "Using GPU device mask '"
              << gpumaskstr << "'" << sendmsg;
    } else {
      msgInfo << "Failed to parse CUDA GPU device mask string '" 
              << gpumaskstr << "'" << sendmsg;
    }
  }

  // This is the very first CUDA API call during VMD startup.
  // There's a >= 2.0 second startup lag associated with it on the DGX-2, 
  // likely due to CUDA runtime library internal initialization overheads
  // across the 16 GPUs.  The first internal call checks the CUDA runtime
  // and driver version compatibility.
  int usabledevices = 0;
  int rc = 0;
  if ((rc=vmd_cuda_num_devices(&numphysdevices)) != VMDCUDA_ERR_NONE) {
    numdevices = 0;
    numphysdevices = 0;

    // Only emit error messages when there are CUDA GPUs on the machine
    // but that they can't be used for some reason
    // XXX turning this off for the time being, as some people have 
    //     NVIDIA drivers installed on machines with no NVIDIA GPU, as can
    //     happen with some distros that package the drivers by default.
    switch (rc) {
      case VMDCUDA_ERR_NODEVICES:
      case VMDCUDA_ERR_SOMEDEVICES:
//        msgInfo << "No CUDA accelerator devices available." << sendmsg;
        break;

#if 0
      case VMDCUDA_ERR_SOMEDEVICES:
        msgWarn << "One or more CUDA accelerators may exist but are not usable." << sendmsg; 
        msgWarn << "Check to make sure that GPU drivers are up to date." << sendmsg;
        break;
#endif

      case VMDCUDA_ERR_DRVMISMATCH:
        msgWarn << "Detected a mismatch between CUDA runtime and GPU driver" << sendmsg; 
        msgWarn << "Check to make sure that GPU drivers are up to date." << sendmsg;
//        msgInfo << "No CUDA accelerator devices available." << sendmsg;
        break;
    }
   
    PROFILE_POP_RANGE();
    return;
  }

  // 
  // Runtime load of NVML shared library (packaged with CUDA driver) to
  // manually obtain function points for query of low-level host platform
  // and GPU hardware details such as the best CPU affinity mask associated 
  // with each GPU, taking into account the NUMA node, PCIe topology, and 
  // NVLink topology that exist on the system.
  nvmlh = wrap_nvml_create();


  // The following loop queries the individual GPU hardware and API 
  // compatibility properties and records their results for subsequent use.
  // This phase of startup costs about 0.05 seconds on a DGX-2 with 16 GPUs.
  if (numphysdevices > 0) {
    cudaavail = 1;

    int i;
    for (i=0; i<numphysdevices; i++) {
      cudadevprops dp;
      memset(&dp, 0, sizeof(dp));
      if (!vmd_cuda_device_props(i, dp.name, sizeof(dp.name),
                                &dp.major, &dp.minor,
                                &dp.membytes, &dp.clockratekhz, 
                                &dp.smcount, &dp.integratedgpu,
                                &dp.asyncenginecount, 
                                &dp.kernelexectimeoutenabled,
                                &dp.canmaphostmem, &dp.computemode,
                                &dp.spdpfpperfratio, 
                                &dp.pageablememaccess,
                                &dp.pageablememaccessuseshostpagetables)) {
        dp.deviceid=i; // save the device index

        // Check that each GPU device has not been excluded by virtue of 
        // being used for display, by a GPU device mask, or by the CUDA
        // device mode being set to a "prohibited" status.
        if (!(dp.kernelexectimeoutenabled && getenv("VMDCUDANODISPLAYGPUS")) &&
            (gpumask & (1 << i)) && 
            (dp.computemode != computeModeProhibited)) {
          devprops.append(dp);
          usabledevices++;
        }
      } else {
        msgWarn << "  Failed to retrieve properties for CUDA accelerator " << i << sendmsg; 
      }
    }
  }

  // assign the final usable device count as the number of available
  // CUDA devices (physical device count is maintained separately)
  numdevices=usabledevices;

  // This code creates a pool of CPU worker threads (one per GPU) that
  // are hereafter responsible for managing each device.  To ensure that
  // the GPUs are all actually usable, each worker thread allocates a 
  // few bytes of memory and executes a trivial kernel on it.
  // On a DGX-2, this phase of startup costs about 7.63 seconds on 16 GPUs.
  devpool_init();

  PROFILE_POP_RANGE();
#endif
}


// destructor
CUDAAccel::~CUDAAccel(void) {
  devpool_fini(); ///< destroy pool of GPU worker threads

#if defined(VMDCUDA)
  // destroy the live connection to NVML library
  if (nvmlh != NULL) {
    wrap_nvml_destroy(nvmlh);
  }
#endif
}


void CUDAAccel::devpool_init(void) {
  cudapool=NULL;

#if defined(VMDCUDA)
  PROFILE_PUSH_RANGE("CUDAAccel::devpool_init()", 0);

  // don't proceed any further if there are no devices or CUDA usage
  // has been disabled by the user
  if (!cudaavail || numdevices == 0 || getenv("VMDNOCUDA") != NULL)
    return;

  // only use as many GPUs as CPU cores we're allowed to use
  int workercount=numdevices;
  if (workercount > wkf_thread_numprocessors())
    workercount=wkf_thread_numprocessors();

  int *devlist = new int[workercount];
  int i;
  for (i=0; i<workercount; i++) {
    devlist[i]=device_index(i);
  }

  msgInfo << "Creating CUDA device pool and initializing hardware..." << sendmsg;
  cudapool=wkf_threadpool_create(workercount, devlist);
  delete [] devlist;

  // associate each worker thread with a specific GPU
  if (getenv("VMDCUDAVERBOSE") != NULL)
    wkf_threadpool_launch(cudapool, vmd_cuda_devpool_setdevice, (void*)"VMD CUDA Dev Init", 1);
  else
    wkf_threadpool_launch(cudapool, vmd_cuda_devpool_setdevice, NULL, 1);

  // clear all available device memory on each of the GPUs
  wkf_threadpool_launch(cudapool, vmd_cuda_devpool_clear_device_mem, NULL, 1);

  // XXX enable fully-connected NVLink peer-to-peer GPU memory access
  //     when requested (not fully generalized yet).  This is done only
  //     once per VMD process, per GPU, and never again.
  if (getenv("VMDCUDAP2PENABLE") != NULL) {
    msgInfo << "Enabling DGX-2 fully-connected NVLink GPU P2P..." << sendmsg;
    wkf_threadpool_launch(cudapool, vmd_cuda_devpool_enable_P2P, NULL, 1);
  }

  PROFILE_POP_RANGE();
#endif
}

void CUDAAccel::devpool_fini(void) {
  if (!cudapool)
    return;

#if defined(VMDCUDA)
  devpool_wait();
  wkf_threadpool_destroy(cudapool);
#endif
  cudapool=NULL;
}

int CUDAAccel::devpool_launch(void *fctn(void *), void *parms, int blocking) {
  if (!cudapool)
    return -1;

  return wkf_threadpool_launch(cudapool, fctn, parms, blocking); 
}

int CUDAAccel::devpool_wait(void) {
  if (!cudapool)
    return -1;

  return wkf_threadpool_wait(cudapool);
}

void CUDAAccel::print_cuda_devices(void) {
  if (getenv("VMDCUDANODISPLAYGPUS")) {
    msgInfo << "Ignoring CUDA-capable GPUs used for display" << sendmsg;
  }

  if (!cudaavail || numdevices == 0) {
    msgInfo << "No CUDA accelerator devices available." << sendmsg;
    return;
  }

  if (nvmlh == NULL) {
    msgInfo << "Unable to load NVML library, GPU-CPU affinity unavailable." << sendmsg;
  }

  // XXX GPU P2P hardware features need to be abstracted by CUDAAccel in the
  //     same way that usable CUDA devices are, so that VMDCUDADEVICEMASK 
  //     affects the record keeping and reporting of P2P connectivity etc.
  //     If the user selects a subset of GPUs, we should disinclude consideration
  //     of the P2P topology that connects GPUs that were masked out. 
  //     We should take into account the mask's impact on the number of P2P islands,
  //     and not count or report links to GPUs that were masked out. 
  //     Since the low-level peer matrix helper function doesn't know anything about
  //     GPU device masks or other control environment variables, CUDAAccel 
  //     should filter the output by copying only the P2P connectivity matrix elements 
  //     that correspond to links between GPUs that are enabled.  The final filtered
  //     and abstracted P2P matrix can then be used by the rest of VMD with appropriate
  //     accessor functions that take the potentially sparse GPU mapping into account.
  int p2plinkcount=0, p2pislands=0;
#if defined(VMDCUDA)
  int numdev=0;
  int *p2pmat=NULL;
  int *p2psupp=NULL;
  int *p2patomics=NULL;
  int *p2parrays=NULL;
  int *perfmat=NULL;

  if (vmd_cuda_peer_matrix(&numdev, &p2pmat, &p2psupp, &p2patomics, &p2parrays,
                           &perfmat, &p2plinkcount, &p2pislands) != VMDCUDA_ERR_NONE) {
    msgWarn << "Unable to ascertain GPU peer-to-peer connectivity" << sendmsg;
  }

  if (p2pmat)
    free(p2pmat);
  if (p2psupp)
    free(p2psupp);
  if (p2patomics)
    free(p2patomics);
  if (p2parrays)
    free(p2parrays);
  if (perfmat)
    free(perfmat);
#endif

  // Report detected GPU hardware and PCIe/NVLink P2P topology
  msgInfo << "Detected " << numdevices << " available CUDA " 
          << ((numdevices > 1) ? "accelerators" : "accelerator:");

  // XXX update to account for device masks...
  if (p2plinkcount > 0) {
    msgInfo << ", " 
            << p2plinkcount << ((p2plinkcount > 1) ? " P2P links, " : " P2P link, ")
            << p2pislands << ((p2pislands > 1) ? " islands" : " island");
  }

  msgInfo << ":" << sendmsg;


  char oldstr[1024], outstr[1024], gpustr[1024], idxprefix[1024];
  int idxrangecount=0,firstidx=-1, lastidx=-1;
  const char *idxfmtstring10gpus  = "[%d]";
  const char *idxfmtspaces10gpus  = "   ";
  const char *idxfmtstring100gpus = "[%2d]";
  const char *idxfmtspaces100gpus = "    ";
  const char *gpuidxfmtstring, *gpuidxfmtspaces;

#if 0
  int outputlineperdevice = 1;
#else
  int outputlineperdevice = (getenv("VMDCUDAOUTPUTLINEPERDEVICE") != NULL);
#endif

  // when enumerating large DGX-2 class hardware, we ensure columns line up
  // by choosing format strings to fit range of device IDs we got
  if (device_index(numdevices-1) > 10) {
    gpuidxfmtstring = idxfmtstring100gpus;
    gpuidxfmtspaces = idxfmtspaces100gpus;
  } else {
    gpuidxfmtstring = idxfmtstring10gpus;
    gpuidxfmtspaces = idxfmtspaces10gpus;
  }

  memset(oldstr, 0, sizeof(oldstr));
  memset(gpustr, 0, sizeof(gpustr));
  memset(idxprefix, 0, sizeof(idxprefix));

  int i;
  int shiftgpuidx=0;
  for (i=0; i<numdevices; i++) {
    memset(outstr, 0, sizeof(outstr));

    // list primary GPU device attributes
    const char *devname = device_name(i);
    sprintf(gpustr, " %-20s %2d SM_%d.%d %.1f GHz", 
            (devname) ? devname : "NULL Device Name!",
            (device_sm_count(i) > 0) ? device_sm_count(i) : 0,
            device_version_major(i), device_version_minor(i),
            device_clock_ghz(i));
    strcpy(outstr, gpustr);

    // list memory capacity 
    int gpumemmb = (device_membytes(i) / (1024 * 1024));
    if (gpumemmb < 1000) {
      sprintf(gpustr, ", %4dMB RAM", gpumemmb);
    } else if (gpumemmb < 10240) {
      sprintf(gpustr, ", %.1fGB RAM", gpumemmb / 1024.0);
    } else {
      // round up to nearest GB
      sprintf(gpustr, ", %dGB RAM", (gpumemmb + 512) / 1024);
    }
    strcat(outstr, gpustr);

    // list optional hardware features and configuration attributes here...
    if (device_computemode(i) == computeModeProhibited) {
      strcat(outstr, ", Compute Mode: Prohibited");
    } else {
      int sfpr = device_spdpfpperfratio(i);
      if (sfpr > 2) {
        sprintf(gpustr, " SP%d", sfpr);
        strcat(outstr, gpustr);
      }

      /// device is an integrated GPU
      if (device_integratedgpu(i)) {
        strcat(outstr, " IGPU");
      }

      /// GPU driver has the Kernel Timeout enabled, potentially 
      /// killing any long running kernels
      if (device_kerneltimeoutenabled(i)) {
        strcat(outstr, " KT");
      }

      /// Number of async DMA engines on the GPU
      if (device_asyncenginecount(i)) {
        sprintf(gpustr, " AE%d", device_asyncenginecount(i));
        strcat(outstr, gpustr);
      }

      /// GPU supports zero-copy direct access to host memory
      if (device_canmaphostmem(i))
        strcat(outstr, " ZC");

      /// GPU supports coherently accessing pageable memory 
      /// without calling cudaHostRegister on it
      if (device_pageablememaccess(i)) {
        /// GPU supports coherently accessing pageable memory 
        /// making direct use of the host's VM page tables
        if (device_pageablememaccessuseshostpagetables(i))
          strcat(outstr, " PMT");
        else 
          strcat(outstr, " PM");
      }
    }

    if (outputlineperdevice) {
      // emit a status line per-device despite any redundancy
      sprintf(idxprefix, gpuidxfmtstring, device_index(i));
      msgInfo << idxprefix << outstr << sendmsg; 
    } else {
      // if the current GPU description is the same as the last one,
      // we don't bother duplicating its listing, and instead we 
      // list the GPU index range(s) with matching descriptive strings
      int newidx = device_index(i);
      if (!strcmp(oldstr, outstr)) {
        // if we have a gap in GPU IDs, we emit the partial index range
        // and continue from the first index after the gap
        if ((newidx - lastidx) > 1) { 
          if (lastidx > firstidx) {
            sprintf(idxprefix, "%d-%d", firstidx, lastidx);
            shiftgpuidx=1;
          } else {
            sprintf(idxprefix, "%s%d", (shiftgpuidx) ? "  " : "", firstidx); 
          }
        
          msgInfo << ((idxrangecount == 0) ? "[" : ",") << idxprefix;
          idxrangecount++;
          firstidx = newidx;
          lastidx = newidx;
        }
        lastidx=newidx;
      } else {
        if (firstidx < 0) {
          firstidx = newidx;
          lastidx = newidx;
          strcpy(oldstr, outstr);
          continue; 
        }
       
        if (lastidx > firstidx) {
          sprintf(idxprefix, "%d-%d", firstidx, lastidx); 
          shiftgpuidx=1;
        } else {
          sprintf(idxprefix, "%s%d", (shiftgpuidx) ? "  " : "", firstidx); 
        }
        msgInfo << ((idxrangecount == 0) ? "[" : ",") << idxprefix;
        msgInfo << "]" << oldstr << sendmsg; 

        idxrangecount = 0;
        firstidx = newidx;
        lastidx = newidx;
        strcpy(oldstr, outstr);
        memset(outstr, 0, sizeof(outstr));
      }
    } 
  } // end of loop over devices 

  if (!outputlineperdevice) {
    if (lastidx > firstidx) {
      sprintf(idxprefix, "%d-%d", firstidx, lastidx); 
    } else {
      sprintf(idxprefix, "%s%d", (shiftgpuidx) ? "  " : "", firstidx); 
    }
    msgInfo << ((idxrangecount == 0) ? "[" : ",") << idxprefix;
    msgInfo << "]";
    if (idxrangecount > 2) {
      msgInfo << ":" << sendmsg;
      msgInfo << gpuidxfmtspaces; // shift to right to line up with column
    }
    msgInfo << oldstr << sendmsg; 
  }
}

int CUDAAccel::num_devices(void) {
  return numdevices;
}

int CUDAAccel::device_index(int dev) {
  return devprops[dev].deviceid;
}

const char *CUDAAccel::device_name(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return NULL;
  return devprops[dev].name; 
}

int CUDAAccel::device_version_major(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return 0; 
  return devprops[dev].major;
}

int CUDAAccel::device_version_minor(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return 0; 
  return devprops[dev].minor;
}

unsigned long CUDAAccel::device_membytes(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return 0; 
  return devprops[dev].membytes;
}

float CUDAAccel::device_clock_ghz(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return 0; 
  return (float) (devprops[dev].clockratekhz / 1000000.0);
}

int CUDAAccel::device_sm_count(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].smcount;
}

int CUDAAccel::device_integratedgpu(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].integratedgpu;
}

int CUDAAccel::device_asyncenginecount(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].asyncenginecount;
}

int CUDAAccel::device_kerneltimeoutenabled(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].kernelexectimeoutenabled;
}

int CUDAAccel::device_canmaphostmem(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].canmaphostmem;
}

int CUDAAccel::device_computemode(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].computemode;
}

int CUDAAccel::device_spdpfpperfratio(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].spdpfpperfratio;
}

int CUDAAccel::device_pageablememaccess(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].pageablememaccess;
}

int CUDAAccel::device_pageablememaccessuseshostpagetables(int dev) {
  if (!cudaavail || dev < 0 || dev >= numdevices)
    return -1; 
  return devprops[dev].pageablememaccessuseshostpagetables;
}

