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
*      $RCSfile: TclGraphLayout.C,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.5 $         $Date: 2020/10/15 16:07:31 $
*
***************************************************************************/
/**
 *  \file TclGraphLayout.C
 *  \brief Graph layout commands needed for interactive
 *  clustering analysis tools, and similar usage.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VMDApp.h"
#include "Inform.h"
#include "GraphLayout.h"
#include <tcl.h>
#include "TclCommands.h"
#include "config.h" // for CMDLEN

//#if defined(VMDCUDA)
//#include "CUDAGraphLayout.h"
//#endif

int layout_fr(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
  int i;
  if ((argc <3) || (argc > 12 )) {
    msgErr << "Usage: node_count iterations [-weights weightlist]" << sendmsg;
    return 0;
  }

  int n = 0;
  int iters = 0;
  float area = 1.0f;
  float kscale = 1.0e3f;
  float tempscale = 0.2f;
  float distance_epsilon = 1.0e-6f;
  float *weights = NULL;

  if (Tcl_GetIntFromObj(interp, objv[1], &n) != TCL_OK) {
    Tcl_AppendResult(interp, "\n node_count incorrectly specified",NULL);
    return TCL_ERROR;
  }
  if (n<1) {
    Tcl_AppendResult(interp, "\n node_count incorrectly specified",NULL);
    return TCL_ERROR;
  }

  if (Tcl_GetIntFromObj(interp, objv[2], &iters) != TCL_OK) {
    Tcl_AppendResult(interp, "\n iterations incorrectly specified",NULL);
    return TCL_ERROR;
  }
  if (iters<0) {
    Tcl_AppendResult(interp, "\n iterations incorrectly specified",NULL);
    return TCL_ERROR;
  }


  for (i=3; i < argc; i++) {
    char *opt = Tcl_GetStringFromObj(objv[i], NULL);
    if (!strcmp(opt, "-weights")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No weights specified",NULL);
        return TCL_ERROR;
      }

      int matlen = 0;
      Tcl_Obj **data;
      if (Tcl_ListObjGetElements(interp, objv[i+1], &matlen, &data) != TCL_OK) {
        return TCL_ERROR;
      }

      if (matlen != (n*n)) {
        Tcl_AppendResult(interp, "Incorrect weight matrix size specified",NULL);
        return TCL_ERROR;
      }

      weights = new float[n*n];
      for (i=0; i<matlen; i++) {
        double tmp;
        if (Tcl_GetDoubleFromObj(interp, data[i], &tmp) != TCL_OK) {
          delete [] weights;
          return TCL_ERROR;
        }
        weights[i] = float(tmp);
      }
    }

    if (!strcmp(opt, "-area")) {
      double tmp;
      if (Tcl_GetDoubleFromObj(interp, objv[i+1], &tmp) != TCL_OK) {
        return TCL_ERROR;
      }
      area = float(tmp);
    }

    if (!strcmp(opt, "-kscale")) {
      double tmp;
      if (Tcl_GetDoubleFromObj(interp, objv[i+1], &tmp) != TCL_OK) {
        return TCL_ERROR;
      }
      kscale = float(tmp);
    }

    if (!strcmp(opt, "-tempscale")) {
      double tmp;
      if (Tcl_GetDoubleFromObj(interp, objv[i+1], &tmp) != TCL_OK) {
        return TCL_ERROR;
      }
      tempscale = float(tmp);
    }

    if (!strcmp(opt, "-distance_epsilon")) {
      double tmp;
      if (Tcl_GetDoubleFromObj(interp, objv[i+1], &tmp) != TCL_OK) {
        return TCL_ERROR;
      }
      distance_epsilon = float(tmp);
    }
  }


//printf("layout_fr()\n");
  GraphLayout *g = new GraphLayout(n, 0);

//printf("init_positions()\n");
  g->init_positions_box();

  if (weights != NULL) {
//printf("add_weight_matrix()\n");
    g->add_weight_matrix(weights);
  }

//printf("compute()\n");
  g->compute(iters, area, kscale, tempscale, distance_epsilon);


//printf("get_vertex_ptrs()\n");
  int numverts=0;
  const float *posx, *posy;
  g->get_vertex_ptrs(numverts, posx, posy);

//printf("generating vertex positions resul list...\n");
  Tcl_Obj *vertexlist = Tcl_NewListObj(0, NULL);
  for (i=0; i<numverts; i++) {
    Tcl_Obj *vertex = Tcl_NewListObj(0, NULL);
    Tcl_ListObjAppendElement(interp, vertex, Tcl_NewDoubleObj(posx[i]));
    Tcl_ListObjAppendElement(interp, vertex, Tcl_NewDoubleObj(posy[i]));
    Tcl_ListObjAppendElement(interp, vertexlist, vertex);
  }
  Tcl_SetObjResult(interp, vertexlist);

//printf("delete g\n");
  delete g;

  if (weights) {
//printf("delete weights\n");
    delete [] weights;
  }

  return 0;
}


int obj_graphlayout(ClientData cd, Tcl_Interp *interp, int argc,
                     Tcl_Obj * const objv[]){
  if (argc < 2) {
    Tcl_SetResult(interp,
    (char *) "Usage: graphlayout <command> [args...]\n"
      "Commands:\n"
      "fr -- Perform Fruchterman-Reingold style spring layout\n"
      ,
      TCL_STATIC);
    return TCL_ERROR;
  }
  char *argv1 = Tcl_GetStringFromObj(objv[1],NULL);

  VMDApp *app = (VMDApp *)cd;
  if (!strupncmp(argv1, "fr", CMDLEN))
    return layout_fr(app, argc-1, objv+1, interp);

  Tcl_SetResult(interp, (char *) "Type 'graphlayout' for summary of usage\n", TCL_VOLATILE);
  return TCL_OK;
}


