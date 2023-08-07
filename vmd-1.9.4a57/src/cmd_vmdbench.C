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
 *      $RCSfile: cmd_vmdbench.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.29 $       $Date: 2020/07/26 06:22:53 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   text commands for benchmarking hardware performance
 ***************************************************************************/

#include <tcl.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AtomSel.h"
#include "Benchmark.h"
#include "config.h"
#include "VMDApp.h"
#include "TclCommands.h"
#include "CUDAKernels.h"
#include "CUDAAccel.h"
#include "WKFThreads.h"

static void cmd_vmdbench_usage(Tcl_Interp *interp) {
  Tcl_AppendResult(interp,
      "usage: vmdbench <command> [args...]\n"
      "vmdbench stream        [N]       - built-in STREAM memory bandwidth test\n",
      "vmdbench cudamadd      [devices] - CUDA multiply-add arithmetic (*)\n",
      "vmdbench cudabusbw     [devices] - CUDA host/device bus bandwidth (*)\n",
      "vmdbench cudaglobmembw [devices] - CUDA global memory bandwidth (*)\n",
      "vmdbench cudadevpool   [N]       - CUDA threadpool run-cycle latency (*)\n",
      "(*) Only available in CUDA-enabled builds of VMD\n",
      NULL);
}

int text_cmd_vmdbench(ClientData cd, Tcl_Interp *interp, int argc, 
                      const char *argv[]) {
#if defined(VMDCUDA)
  VMDApp *app = (VMDApp *)cd; // need VMDApp ptr GPU threadpool access
#endif

  if (argc == 1) {
    cmd_vmdbench_usage(interp);
    return TCL_ERROR;
  }

  if (argc >= 2) {
    if (!strupncmp(argv[1], "stream", CMDLEN)) {
      double times[8], mbsec[8];
      int N = 1024*1024 * 16;

      if (argc == 3) {
        if (Tcl_GetInt(interp, argv[2], &N) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench stream", NULL);
          return TCL_ERROR;
        }
      }

      int rc = stream_bench(N, times, mbsec);
      if (rc) {
        Tcl_AppendResult(interp,
          "unable to complete stream benchmark, out of memory", NULL);
        return TCL_ERROR;
      }

      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      const char *benchnames[] = {
        "copy (double)",
        "scale (double)",
        "add (double)",
        "triad (double)",
        "copy (float)",
        "scale (float)",
        "add (float)",
        "triad (float)"
      };

      Tcl_Obj *colNameObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Test", -1)); 
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Time", -1)); 
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("MB/sec", -1)); 
      Tcl_ListObjAppendElement(interp, tcl_result, colNameObj);

      int i;     
      for (i=0; i<8; i++) {
        Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewStringObj(benchnames[i], -1)); 
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(times[i])); 
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(mbsec[i])); 
        Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);

      }
      Tcl_SetObjResult(interp, tcl_result);

      return TCL_OK;

    } else if (!strupncmp(argv[1], "minmaxmean_1fv", CMDLEN)) {
      double runtime, mbsec;
      int N = 1024*1024 * 16;
      int reps = 1;

      if (argc >= 3) {
        if (Tcl_GetInt(interp, argv[2], &N) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench minmaxmean_1fv", NULL);
          return TCL_ERROR;
        }
      }
      if (argc == 4) {
        if (Tcl_GetInt(interp, argv[3], &reps) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench minmaxmean_1fv", NULL);
          return TCL_ERROR;
        }
      }

      vmdbench_minmaxmean_1fv(N, reps, runtime, mbsec);

      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_Obj *colNameObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Time", -1)); 
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("MB/sec", -1)); 
      Tcl_ListObjAppendElement(interp, tcl_result, colNameObj);

      Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(runtime)); 
      Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(mbsec)); 

      Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);
      Tcl_SetObjResult(interp, tcl_result);

      return TCL_OK;

    } else if (!strupncmp(argv[1], "minmax_3fv", CMDLEN)) {
      double runtime, mbsec;
      int N = 1024*1024 * 16;
      int reps = 1;

      if (argc >= 3) {
        if (Tcl_GetInt(interp, argv[2], &N) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench minmax_3fv", NULL);
          return TCL_ERROR;
        }
      }
      if (argc == 4) {
        if (Tcl_GetInt(interp, argv[3], &reps) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench minmax_3fv", NULL);
          return TCL_ERROR;
        }
      }

      vmdbench_minmax_3fv(N, reps, runtime, mbsec);

      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_Obj *colNameObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Time", -1)); 
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("MB/sec", -1)); 
      Tcl_ListObjAppendElement(interp, tcl_result, colNameObj);

      Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(runtime)); 
      Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(mbsec)); 

      Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);
      Tcl_SetObjResult(interp, tcl_result);

      return TCL_OK;

    } else if (!strupncmp(argv[1], "analyze_selection", CMDLEN)) {
      double runtime, mbsec;
      int N = 1024*1024 * 16;
      int reps = 1;

      if (argc >= 3) {
        if (Tcl_GetInt(interp, argv[2], &N) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench analyze_selection", NULL);
          return TCL_ERROR;
        }
      }
      if (argc == 4) {
        if (Tcl_GetInt(interp, argv[3], &reps) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench analyze_selection", NULL);
          return TCL_ERROR;
        }
      }

      vmdbench_analyze_selection(N, reps, runtime, mbsec);

      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_Obj *colNameObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Time", -1)); 
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("MB/sec", -1)); 
      Tcl_ListObjAppendElement(interp, tcl_result, colNameObj);

      Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(runtime)); 
      Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(mbsec)); 

      Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);
      Tcl_SetObjResult(interp, tcl_result);

      return TCL_OK;

    } else if (!strupncmp(argv[1], "cudamadd", CMDLEN)) {
#if defined(VMDCUDA)
      int numdevs, physnumdevs;
      int *devlist = NULL;
      vmd_cuda_num_devices(&physnumdevs);
      numdevs = physnumdevs;
#if !defined(VMDTHREADS)
      numdevs = 1;
#endif

      // handle optional device list arguments
      if (argc > 2) {
        if ((argc-2) > numdevs) {
          Tcl_AppendResult(interp, "vmdbench: bad device argument", NULL);
          return TCL_ERROR;
        } else {
          numdevs = argc-2;
        }
        devlist = (int *) malloc(numdevs * sizeof(int));
        int arg, dev;
        for (arg=0; arg<numdevs; arg++) {
          if (Tcl_GetInt(interp, argv[arg+2], &dev) != TCL_OK) {
            Tcl_AppendResult(interp, "vmdbench: bad device argument", NULL);
            free(devlist);
            return TCL_ERROR;
          }
          if (dev < 0 || dev >= physnumdevs) {
            Tcl_AppendResult(interp, "vmdbench: device argument out of range", NULL);
            free(devlist);
            return TCL_ERROR;
          }
          devlist[arg] = dev;
        } 
      }

      double *gflops = (double *) malloc(numdevs * sizeof(double));
      int testloops=10;
      if (getenv("VMDMADDLOOPS") != NULL)
        testloops = atoi(getenv("VMDMADDLOOPS"));

      vmd_cuda_madd_gflops(numdevs, devlist, gflops, testloops);

      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_Obj *colNameObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Device", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("GFLOPS", -1));
      Tcl_ListObjAppendElement(interp, tcl_result, colNameObj);

      int i;
      for (i=0; i<numdevs; i++) {
        Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
        if (devlist != NULL) 
          Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewIntObj(devlist[i]));
        else
          Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewIntObj(i));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(gflops[i]));
        Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);
      }
      Tcl_SetObjResult(interp, tcl_result);

      if (devlist)
        free(devlist);

      return TCL_OK;
#else 
      Tcl_AppendResult(interp, "CUDA Acceleration not available in this build", NULL);
      return TCL_ERROR;
#endif
    } else if (!strupncmp(argv[1], "cudabusbw", CMDLEN)) {
#if defined(VMDCUDA)
      int numdevs, physnumdevs;
      int *devlist = NULL;
      vmd_cuda_num_devices(&physnumdevs);
      numdevs = physnumdevs;
#if !defined(VMDTHREADS)
      numdevs = 1;
#endif

      // handle optional device list arguments
      if (argc > 2) {
        if ((argc-2) > numdevs) {
          Tcl_AppendResult(interp, "vmdbench: bad device argument", NULL);
          return TCL_ERROR;
        } else {
          numdevs = argc-2;
        }
        devlist = (int *) malloc(numdevs * sizeof(int));
        int arg, dev;
        for (arg=0; arg<numdevs; arg++) {
          if (Tcl_GetInt(interp, argv[arg+2], &dev) != TCL_OK) {
            Tcl_AppendResult(interp, "vmdbench: bad device argument", NULL);
            free(devlist);
            return TCL_ERROR;
          }
          if (dev < 0 || dev >= physnumdevs) {
            Tcl_AppendResult(interp, "vmdbench: device argument out of range", NULL);
            free(devlist);
            return TCL_ERROR;
          }
          devlist[arg] = dev;
        } 
      }

      double *hdmbsec = (double *) malloc(numdevs * sizeof(double));
      double *hdlatusec = (double *) malloc(numdevs * sizeof(double));
      double *phdmbsec = (double *) malloc(numdevs * sizeof(double));
      double *phdlatusec = (double *) malloc(numdevs * sizeof(double));
      double *dhmbsec = (double *) malloc(numdevs * sizeof(double));
      double *dhlatusec = (double *) malloc(numdevs * sizeof(double));
      double *pdhmbsec = (double *) malloc(numdevs * sizeof(double));
      double *pdhlatusec = (double *) malloc(numdevs * sizeof(double));

      vmd_cuda_bus_bw(numdevs, devlist, 
                      hdmbsec, hdlatusec, phdmbsec, phdlatusec,
                      dhmbsec, dhlatusec, pdhmbsec, pdhlatusec);

      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_Obj *colNameObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Device", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Host-device bandwidth (MB/sec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Host-device latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Host-device pinned bandwidth (MB/sec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Host-device pinned latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Device-host bandwidth (MB/sec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Device-host latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Device-host pinned bandwidth (MB/sec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Device-host pinned latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, tcl_result, colNameObj);

      int i;
      for (i=0; i<numdevs; i++) {
        Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
        if (devlist != NULL) 
          Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewIntObj(devlist[i]));
        else
          Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewIntObj(i));

        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(hdmbsec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(hdlatusec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(phdmbsec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(phdlatusec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(dhmbsec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(dhlatusec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(pdhmbsec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(pdhlatusec[i]));
        Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);
      }
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
#else 
      Tcl_AppendResult(interp, "CUDA Acceleration not available in this build", NULL);
      return TCL_ERROR;
#endif
    } else if (!strupncmp(argv[1], "cudaglobmembw", CMDLEN)) {
#if defined(VMDCUDA)
      int numdevs, physnumdevs;
      int *devlist = NULL;
      vmd_cuda_num_devices(&physnumdevs);
      numdevs = physnumdevs;
#if !defined(VMDTHREADS)
      numdevs = 1;
#endif

      // handle optional device list arguments
      if (argc > 2) {
        if ((argc-2) > numdevs) {
          Tcl_AppendResult(interp, "vmdbench: bad device argument", NULL);
          return TCL_ERROR;
        } else {
          numdevs = argc-2;
        }
        devlist = (int *) malloc(numdevs * sizeof(int));
        int arg, dev;
        for (arg=0; arg<numdevs; arg++) {
          if (Tcl_GetInt(interp, argv[arg+2], &dev) != TCL_OK) {
            Tcl_AppendResult(interp, "vmdbench: bad device argument", NULL);
            free(devlist);
            return TCL_ERROR;
          }
          if (dev < 0 || dev >= physnumdevs) {
            Tcl_AppendResult(interp, "vmdbench: device argument out of range", NULL);
            free(devlist);
            return TCL_ERROR;
          }
          devlist[arg] = dev;
        } 
      }

      double *memsetgbsec = (double *) malloc(numdevs * sizeof(double));
      double *memcpygbsec = (double *) malloc(numdevs * sizeof(double));

      vmd_cuda_globmem_bw(numdevs, devlist, memsetgbsec, memcpygbsec);

      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_Obj *colNameObj = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Device", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Memory set bandwidth (GB/sec)", -1));
      Tcl_ListObjAppendElement(interp, colNameObj, Tcl_NewStringObj("Memory copy bandwidth (GB/sec)", -1));
      Tcl_ListObjAppendElement(interp, tcl_result, colNameObj);

      int i;
      for (i=0; i<numdevs; i++) {
        Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
        if (devlist != NULL) 
          Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewIntObj(devlist[i]));
        else
          Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewIntObj(i));

        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(memsetgbsec[i]));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(memcpygbsec[i]));
        Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);
      }
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
#else 
      Tcl_AppendResult(interp, "CUDA Acceleration not available in this build", NULL);
      return TCL_ERROR;
#endif
    } else if (!strupncmp(argv[1], "cudadevpool", CMDLEN)) {
#if defined(VMDCUDA)
      int N=1;
      if (argc == 3) {
        if (Tcl_GetInt(interp, argv[2], &N) != TCL_OK) {
          Tcl_AppendResult(interp, " in vmdbench cudadevpool", NULL);
          return TCL_ERROR;
        }
      }

      wkf_threadpool_t * devpool = app->cuda->get_cuda_devpool();
      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewStringObj("Empty kernel launch latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewStringObj("Device pool barrier latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewStringObj("Device pool empty run cycle latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewStringObj("Device pool tile run latency (usec)", -1));
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewStringObj("Device pool GPU kernel tile latency (usec)", -1));

      int i;
      double kernlaunchlatency, barlatency;
      double cyclelatency, tilelatency;
      double kernellatency;
      for (i=0; i<2; i++) {
        vmd_cuda_devpool_latency(devpool, N, &kernlaunchlatency,
                                 &barlatency, &cyclelatency, 
                                 &tilelatency, &kernellatency);

        // do one warmup pass before we report the benchmark numbers
        if (i < 1)
          continue;

        // report the results
        Tcl_Obj *rowListObj = Tcl_NewListObj(0, NULL);
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(kernlaunchlatency*1000000));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(barlatency*1000000));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(cyclelatency*1000000));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(tilelatency*1000000));
        Tcl_ListObjAppendElement(interp, rowListObj, Tcl_NewDoubleObj(kernellatency*1000000));
        Tcl_ListObjAppendElement(interp, tcl_result, rowListObj);
      }

      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
#else 
      Tcl_AppendResult(interp, "CUDA Acceleration not available in this build", NULL);
      return TCL_ERROR;
#endif
    } else if (!strupncmp(argv[1], "cudaoocgdsio", CMDLEN)) {
#if defined(VMDCUDA)
      int first = 0;  // start with first frame by default
      int last = -1;  // finish with last frame by default
      int step = 1;   // use all frames by default
      int i;
      int nfiles = 0;
      char **trjfileset = NULL;

      if (argc < 2) {
        cmd_vmdbench_usage(interp);
        return TCL_ERROR;
      }
      AtomSel *sel = tcl_commands_get_sel(interp, argv[2]);
      if (!sel) {
        Tcl_AppendResult(interp, "cudaoocgdsio: no atom selection", NULL);
        return TCL_ERROR;
      }

      for (i=3; i<argc; i+=2) {
        const char *argvcur = argv[i];
        if (!strupncmp(argvcur, "first", CMDLEN)) {
          first = atoi(argv[i+1]);
        } else if (!strupncmp(argvcur, "last", CMDLEN)) {
          last = atoi(argv[i+1]);
        } else if (!strupncmp(argvcur, "step", CMDLEN)) {
          step = atoi(argv[i+1]);
        } else if (!strupncmp(argvcur, "files", CMDLEN)) {
          int list_num;
          const char **list_strs;
          if (Tcl_SplitList(interp, argv[i+1], &list_num, &list_strs) != TCL_OK) {
            Tcl_AppendResult(interp, "cudaoocgdsio: bad trajectory file list", NULL);
            return TCL_ERROR;
          }
         
          int f;
          nfiles = list_num;
          trjfileset = (char **) calloc(1, nfiles * sizeof(const char *));
          for (f=0; f<nfiles; f++) {
            trjfileset[f] = strdup(list_strs[f]);
            printf("File[%d] '%s'\n", f, trjfileset[f]);
          }
          Tcl_Free((char *)list_strs);
        } else {
          Tcl_AppendResult(interp, "cudaoocgdsio: invalid syntax, no such keyword: ", argvcur, NULL);
          return TCL_ERROR;
        }
      }

      int ret_val = gpu_ooc_bench(app->cuda->get_cuda_devpool(), 
                                  nfiles, (const char **) trjfileset, sel,
                                  first, last, step);

      if (ret_val < 0) {
        Tcl_AppendResult(interp, "cudaoocgdsio: an error occured", NULL);
        return TCL_ERROR;
      }

      if (trjfileset != NULL) {
        int f;
        for (f=0; f<nfiles; f++) {
          free(trjfileset[f]);
        }
        free(trjfileset);
      }

      return TCL_OK;
#else 
      Tcl_AppendResult(interp, "CUDA Acceleration not available in this build", NULL);
      return TCL_ERROR;
#endif
    } else {
      cmd_vmdbench_usage(interp);
      return TCL_ERROR;
    }
  } else {
    cmd_vmdbench_usage(interp);
    return TCL_ERROR;
  }
  
  // if here, everything worked out ok
  return TCL_OK;
}


