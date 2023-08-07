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
 *	$RCSfile: CUDAAccel.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.28 $	$Date: 2020/02/26 04:22:39 $
 *
 ***************************************************************************/
/**
 * \file CUDAAccel.h
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

#ifndef CUDACCEL_H
#define CUDACCEL_H

#include "WKFThreads.h"
#include "CUDAWrapNVML.h"

typedef struct {
  int deviceid;
  char name[80];
  int major;
  int minor;
  unsigned long membytes;
  int clockratekhz;
  int smcount;
  int integratedgpu;
  int asyncenginecount;
  int kernelexectimeoutenabled;
  int canmaphostmem;
  int computemode;
  int spdpfpperfratio;
  int pageablememaccess;
  int pageablememaccessuseshostpagetables;
} cudadevprops;

/// manages enumeration and initialization of CUDA devices
class CUDAAccel {
private:
  int cudaavail;              ///< flag: CUDA runtime is operable for VMD
  int numdevices;             ///< number of CUDA GPU accelerators available
  int numphysdevices;         ///< number of CUDA GPUs accelerators that exist

  wrap_nvml_handle *nvmlh;    ///< Handle to NVML wrapper library

  ResizeArray<cudadevprops> devprops; ///< detailed per-GPU device properties

  wkf_threadpool_t *cudapool; ///< Pool of GPU management threads

  void devpool_init(void);    ///< init pool GPU devices flagging any that fail
  void devpool_fini(void);    ///< shutdown device pool and destroy conexts

  // convenience enum to match CUDA driver APIs
  enum { computeModeDefault=0, 
         computeModeExclusive=1,
         computeModeProhibited=2 }; // computeMode;
 
public:
  CUDAAccel(void);
  virtual ~CUDAAccel(void);

  // functions for enumerating CUDA GPU accelerator devices
  // and their attributes
  void print_cuda_devices(void);
  int num_devices(void);
  int device_index(int dev);
  const char *device_name(int dev);
  int device_version_major(int dev);
  int device_version_minor(int dev);
  unsigned long device_membytes(int dev);
  float device_clock_ghz(int dev);
  int device_sm_count(int dev);
  int device_integratedgpu(int dev);
  int device_asyncenginecount(int dev);
  int device_kerneltimeoutenabled(int dev);
  int device_canmaphostmem(int dev);
  int device_computemode(int dev);
  int device_spdpfpperfratio(int dev);
  int device_pageablememaccess(int dev);
  int device_pageablememaccessuseshostpagetables(int dev);

  // functions for operating on an open pool of CUDA devices
  int devpool_launch(void *fctn(void *), void *parms, int blocking);
  int devpool_wait(void);
  wkf_threadpool_t * get_cuda_devpool(void) { return cudapool; }

};

#endif



