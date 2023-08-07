/**
 * \file CUDAWrapNVML.cu
 * \brief Runtime loader and wrapper for the NVIDIA NVML GPU management library.
 *
 * A trivial little dlopen()-based wrapper library for the
 * NVIDIA NVML library, to allow runtime discovery of NVML on an
 * arbitrary system.  This is all very hackish and simple-minded, but
 * it serves my immediate needs in the short term until NVIDIA provides
 * a static NVML wrapper library themselves, hopefully in
 * CUDA 6.5 or maybe sometime shortly after.
 *
 * $Revision: 1.9 $       $Date: 2020/06/25 16:05:41 $
 *
 * \author John E. Stone - john.stone@gmail.com
 * \copyright 
 * This trivial code is made available under the "new" 3-clause BSD license,
 * and/or any of the GPL licenses you prefer.
 * Feel free to use the code and modify as you see fit.
 */

#include <stdio.h>
#include <stdlib.h>
#include "CUDAWrapNVML.h"
#include "cuda_runtime.h"

/*
 * Wrappers to emulate dlopen() on other systems like Windows
 */
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
#include <windows.h>
static void *wrap_dlopen(const char *filename) {
  return (void *)LoadLibrary(filename);
}
static void *wrap_dlsym(void *h, const char *sym) {
  return (void *)GetProcAddress((HINSTANCE)h, sym);
}
static int wrap_dlclose(void *h) {
  /* FreeLibrary returns nonzero on success */
  return (!FreeLibrary((HINSTANCE)h));
}
#else
/* assume we can use dlopen itself... */
#include <dlfcn.h>
static void *wrap_dlopen(const char *filename) {
  return dlopen(filename, RTLD_NOW);
}
static void *wrap_dlsym(void *h, const char *sym) {
  return dlsym(h, sym);
}
static int wrap_dlclose(void *h) {
  return dlclose(h);
}
#endif

#if defined(__cplusplus)
extern "C" {
#endif

wrap_nvml_handle * wrap_nvml_create() {
  int i=0;
  wrap_nvml_handle *nvmlh = NULL;

  /* 
   * We use hard-coded library installation locations for the time being...
   * No idea where or if libnvidia-ml.so is installed on MacOS X, a 
   * deep scouring of the filesystem on one of the Mac CUDA build boxes
   * I used turned up nothing, so for now it's not going to work on OSX.
   */
#if defined(_WIN64)
  /* 64-bit Windows */
  const char *plibnvidia_ml[] = {"c:/Program Files/NVIDIA Corporation/NVSMI/nvml.dll", NULL};
#elif defined(_WIN32) || defined(_MSC_VER)
  /* 32-bit Windows */
  const char *plibnvidia_ml[] = {"c:/Program Files (x86)/NVIDIA Corporation/NVSMI/nvml.dll", NULL};
#elif defined(__linux) && (defined(__i386__) || defined(__ARM_ARCH_7A__))
  /* 32-bit linux assumed */
  const char *plibnvidia_ml[] = {"/usr/lib/libnvidia-ml.so", 
                                 "/usr/lib/libnvidia-ml.so.1", 
                                 "/usr/lib/x86-linux-gnu/libnvidia-ml.so",
                                 "/usr/lib/x86-linux-gnu/libnvidia-ml.so.1",
                                 NULL};
#elif defined(__linux)
  /* 64-bit linux assumed */
  const char *plibnvidia_ml[] = {"/usr/lib64/libnvidia-ml.so", 
                                 "/usr/lib64/libnvidia-ml.so.1",
                                 "/usr/lib/aarch64-linux-gnu/libnvidia-ml.so",
                                 "/usr/lib/aarch64-linux-gnu/libnvidia-ml.so.1",
                                 "/usr/lib/x86_64-linux-gnu/libnvidia-ml.so",
                                 "/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1",
                                 NULL};
#else
  //#error "Unrecognized platform: need NVML DLL path for this platform..."
  const char *plibnvidia_ml[] = {NULL};
#endif


  void *nvml_dll = NULL;
  const char *libnvidia_ml = NULL;

  // allow explicit user override of path if needed
  if (getenv("LIBNVMLPATH") != NULL)
    nvml_dll = wrap_dlopen(getenv("LIBNVMLPATH"));

  int sopathidx = 0;
  libnvidia_ml = plibnvidia_ml[sopathidx];
  while ((nvml_dll == NULL) && (libnvidia_ml != NULL)) {
    nvml_dll = wrap_dlopen(libnvidia_ml);
    sopathidx++;
    libnvidia_ml = plibnvidia_ml[sopathidx];
  }
  if (nvml_dll == NULL)
    return NULL;

  nvmlh = (wrap_nvml_handle *) calloc(1, sizeof(wrap_nvml_handle));

  nvmlh->nvml_dll = nvml_dll;  

  nvmlh->nvmlInit = (wrap_nvmlReturn_t (*)(void)) 
    wrap_dlsym(nvmlh->nvml_dll, "nvmlInit");

  nvmlh->nvmlDeviceGetCount = (wrap_nvmlReturn_t (*)(int *)) 
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetCount_v2");

  nvmlh->nvmlDeviceGetHandleByIndex = (wrap_nvmlReturn_t (*)(int, wrap_nvmlDevice_t *)) 
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetHandleByIndex_v2");

  nvmlh->nvmlDeviceGetPciInfo = (wrap_nvmlReturn_t (*)(wrap_nvmlDevice_t, wrap_nvmlPciInfo_t *)) 
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetPciInfo");

  nvmlh->nvmlDeviceGetName = (wrap_nvmlReturn_t (*)(wrap_nvmlDevice_t, char *, int))
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetName");

  nvmlh->nvmlDeviceGetTemperature = (wrap_nvmlReturn_t (*)(wrap_nvmlDevice_t, int, unsigned int *))
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetTemperature");

  nvmlh->nvmlDeviceGetFanSpeed = (wrap_nvmlReturn_t (*)(wrap_nvmlDevice_t, unsigned int *))
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetFanSpeed");

  nvmlh->nvmlDeviceGetPowerUsage = (wrap_nvmlReturn_t (*)(wrap_nvmlDevice_t, unsigned int *))
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetPowerUsage");

  nvmlh->nvmlDeviceGetCpuAffinity = (wrap_nvmlReturn_t (*)(wrap_nvmlDevice_t, unsigned int, unsigned long *))
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceGetCpuAffinity");

  nvmlh->nvmlDeviceSetCpuAffinity = (wrap_nvmlReturn_t (*)(wrap_nvmlDevice_t))
    wrap_dlsym(nvmlh->nvml_dll, "nvmlDeviceSetCpuAffinity");

  nvmlh->nvmlShutdown = (wrap_nvmlReturn_t (*)()) 
    wrap_dlsym(nvmlh->nvml_dll, "nvmlShutdown");

  if (nvmlh->nvmlInit == NULL || 
      nvmlh->nvmlShutdown == NULL ||
      nvmlh->nvmlDeviceGetCount == NULL ||
      nvmlh->nvmlDeviceGetHandleByIndex == NULL || 
      nvmlh->nvmlDeviceGetPciInfo == NULL ||
      nvmlh->nvmlDeviceGetName == NULL ||
      nvmlh->nvmlDeviceGetTemperature == NULL ||
      nvmlh->nvmlDeviceGetFanSpeed == NULL ||
      nvmlh->nvmlDeviceGetPowerUsage == NULL ||
      nvmlh->nvmlDeviceGetCpuAffinity == NULL ||
      nvmlh->nvmlDeviceSetCpuAffinity == NULL
      ) {
#if 0
    printf("Failed to obtain all required NVML function pointers\n");
#endif
    wrap_dlclose(nvmlh->nvml_dll);
    free(nvmlh);
    return NULL;
  }

  nvmlh->nvmlInit();
  nvmlh->nvmlDeviceGetCount(&nvmlh->nvml_gpucount);

  /* Query CUDA device count, in case it doesn't agree with NVML, since  */
  /* CUDA will only report GPUs with compute capability greater than 1.0 */ 
  if (cudaGetDeviceCount(&nvmlh->cuda_gpucount) != cudaSuccess) {
#if 0
    printf("Failed to query CUDA device count!\n");
#endif
    wrap_dlclose(nvmlh->nvml_dll);
    free(nvmlh);
    return NULL;
  }

  nvmlh->devs = (wrap_nvmlDevice_t *) calloc(nvmlh->nvml_gpucount, sizeof(wrap_nvmlDevice_t));
  nvmlh->nvml_pci_domain_id = (unsigned int*) calloc(nvmlh->nvml_gpucount, sizeof(unsigned int));
  nvmlh->nvml_pci_bus_id = (unsigned int*) calloc(nvmlh->nvml_gpucount, sizeof(unsigned int));
  nvmlh->nvml_pci_device_id = (unsigned int*) calloc(nvmlh->nvml_gpucount, sizeof(unsigned int));
  nvmlh->nvml_cuda_device_id = (int*) calloc(nvmlh->nvml_gpucount, sizeof(int));
  nvmlh->cuda_nvml_device_id = (int*) calloc(nvmlh->cuda_gpucount, sizeof(int));

  /* Obtain GPU device handles we're going to need repeatedly... */
  for (i=0; i<nvmlh->nvml_gpucount; i++) {
    nvmlh->nvmlDeviceGetHandleByIndex(i, &nvmlh->devs[i]);
  } 

  /* Query PCI info for each NVML device, and build table for mapping of */
  /* CUDA device IDs to NVML device IDs and vice versa                   */
  for (i=0; i<nvmlh->nvml_gpucount; i++) {
    wrap_nvmlPciInfo_t pciinfo;
    nvmlh->nvmlDeviceGetPciInfo(nvmlh->devs[i], &pciinfo);
    nvmlh->nvml_pci_domain_id[i] = pciinfo.domain;
    nvmlh->nvml_pci_bus_id[i]    = pciinfo.bus;
    nvmlh->nvml_pci_device_id[i] = pciinfo.device;
  }

  /* build mapping of NVML device IDs to CUDA IDs */
  for (i=0; i<nvmlh->nvml_gpucount; i++) {
    nvmlh->nvml_cuda_device_id[i] = -1;
  } 
  for (i=0; i<nvmlh->cuda_gpucount; i++) {
    cudaDeviceProp props;
    nvmlh->cuda_nvml_device_id[i] = -1;

    if (cudaGetDeviceProperties(&props, i) == cudaSuccess) {
      int j;
      for (j=0; j<nvmlh->nvml_gpucount; j++) {
        if ((nvmlh->nvml_pci_domain_id[j] == props.pciDomainID) &&
            (nvmlh->nvml_pci_bus_id[j]    == props.pciBusID) &&
            (nvmlh->nvml_pci_device_id[j] == props.pciDeviceID)) {
#if 0
          printf("CUDA GPU[%d] matches NVML GPU[%d]\n", i, j);
#endif
          nvmlh->nvml_cuda_device_id[j] = i;
          nvmlh->cuda_nvml_device_id[i] = j;
        }
      }
    }
  }

  return nvmlh;
}


int wrap_nvml_destroy(wrap_nvml_handle *nvmlh) {
  nvmlh->nvmlShutdown();

  if (nvmlh->nvml_pci_domain_id != NULL)
    free(nvmlh->nvml_pci_domain_id);

  if (nvmlh->nvml_pci_bus_id != NULL)
    free(nvmlh->nvml_pci_bus_id);

  if (nvmlh->nvml_pci_device_id != NULL)
    free(nvmlh->nvml_pci_device_id);

  if (nvmlh->nvml_cuda_device_id != NULL)
    free(nvmlh->nvml_cuda_device_id);

  if (nvmlh->cuda_nvml_device_id != NULL)
    free(nvmlh->cuda_nvml_device_id);

  if (nvmlh->devs != NULL)
    free(nvmlh->devs);

  wrap_dlclose(nvmlh->nvml_dll);
  free(nvmlh);
  return 0;
}


int wrap_nvml_get_gpucount(wrap_nvml_handle *nvmlh, int *gpucount) {
  *gpucount = nvmlh->nvml_gpucount;
  return 0; 
}

int wrap_cuda_get_gpucount(wrap_nvml_handle *nvmlh, int *gpucount) {
  *gpucount = nvmlh->cuda_gpucount;
  return 0; 
}

int wrap_nvml_get_gpu_name(wrap_nvml_handle *nvmlh,
                           int cudaindex, 
                           char *namebuf,
                           int bufsize) {
  int gpuindex = nvmlh->cuda_nvml_device_id[cudaindex];
  if (gpuindex < 0 || gpuindex >= nvmlh->nvml_gpucount)
    return -1;

  if (nvmlh->nvmlDeviceGetName(nvmlh->devs[gpuindex], namebuf, bufsize) != WRAPNVML_SUCCESS)
    return -1; 

  return 0;
}


int wrap_nvml_get_tempC(wrap_nvml_handle *nvmlh,
                        int cudaindex, unsigned int *tempC) {
  wrap_nvmlReturn_t rc;
  int gpuindex = nvmlh->cuda_nvml_device_id[cudaindex];
  if (gpuindex < 0 || gpuindex >= nvmlh->nvml_gpucount)
    return -1;

  rc = nvmlh->nvmlDeviceGetTemperature(nvmlh->devs[gpuindex], 0u /* NVML_TEMPERATURE_GPU */, tempC);
  if (rc != WRAPNVML_SUCCESS) {
    return -1; 
  }

  return 0;
}


int wrap_nvml_get_fanpcnt(wrap_nvml_handle *nvmlh,
                          int cudaindex, unsigned int *fanpcnt) {
  wrap_nvmlReturn_t rc;
  int gpuindex = nvmlh->cuda_nvml_device_id[cudaindex];
  if (gpuindex < 0 || gpuindex >= nvmlh->nvml_gpucount)
    return -1;

  rc = nvmlh->nvmlDeviceGetFanSpeed(nvmlh->devs[gpuindex], fanpcnt);
  if (rc != WRAPNVML_SUCCESS) {
    return -1; 
  }

  return 0;
}


int wrap_nvml_get_power_usage(wrap_nvml_handle *nvmlh,
                              int cudaindex,
                              unsigned int *milliwatts) {
  int gpuindex = nvmlh->cuda_nvml_device_id[cudaindex];
  if (gpuindex < 0 || gpuindex >= nvmlh->nvml_gpucount)
    return -1;

  if (nvmlh->nvmlDeviceGetPowerUsage(nvmlh->devs[gpuindex], milliwatts) != WRAPNVML_SUCCESS)
    return -1; 

  return 0;
}


int wrap_nvml_get_cpu_affinity(wrap_nvml_handle *nvmlh,
                               int cudaindex,
                               unsigned int cpuSetSize,
                               unsigned long *cpuSet) {
  int gpuindex = nvmlh->cuda_nvml_device_id[cudaindex];
  if (gpuindex < 0 || gpuindex >= nvmlh->nvml_gpucount)
    return -1;

  if (nvmlh->nvmlDeviceGetCpuAffinity(nvmlh->devs[gpuindex], cpuSetSize, cpuSet) != WRAPNVML_SUCCESS)
    return -1; 

  return 0;
}


int wrap_nvml_set_cpu_affinity(wrap_nvml_handle *nvmlh,
                               int cudaindex) {
  int gpuindex = nvmlh->cuda_nvml_device_id[cudaindex];
  if (gpuindex < 0 || gpuindex >= nvmlh->nvml_gpucount)
    return -1;

  if (nvmlh->nvmlDeviceSetCpuAffinity(nvmlh->devs[gpuindex]) != WRAPNVML_SUCCESS)
    return -1; 
  
  return 0;
}


#if defined(__cplusplus)
}
#endif


