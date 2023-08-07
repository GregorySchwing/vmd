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
 *	$RCSfile: CUDADispCmds.cu,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.11 $	$Date: 2020/02/26 04:22:39 $
 *
 ***************************************************************************/
/**
 * \file CUDADispCmds.cu
 * \brief GPU-native VMDDisplayList display command buffer packing routines
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Scene.h"
#include "DispCmds.h"
#include "utilities.h"
#include "Matrix4.h"
#include "VMDDisplayList.h"

//*************************************************************
// draw a mesh consisting of vertices, facets, colors, normals etc.
void DispCmdTriMesh::cuda_putdata(const float * vertices_d,
                                  const float * normals_d,
                                  const float * colors_d,
                                  int num_facets,
                                  VMDDisplayList * dobj) {
  // make a triangle mesh (no strips)
  DispCmdTriMesh *ptr;
  if (colors_d == NULL) {
    ptr = (DispCmdTriMesh *)
                (dobj->append(DTRIMESH_C3F_N3F_V3F, sizeof(DispCmdTriMesh) +
                              sizeof(float) * num_facets * 3 * 6));
  } else {
    ptr = (DispCmdTriMesh *)
                (dobj->append(DTRIMESH_C3F_N3F_V3F, sizeof(DispCmdTriMesh) +
                              sizeof(float) * num_facets * 3 * 9));
  }

  if (ptr == NULL)
    return;

  ptr->numverts=num_facets * 3;
  ptr->numfacets=num_facets;

  float *c=NULL, *n=NULL, *v=NULL;
  if (colors_d == NULL) {
    ptr->pervertexcolors=0;
    ptr->getpointers(n, v);
  } else {
    ptr->pervertexcolors=1;
    ptr->getpointers(c, n, v);
    cudaMemcpy(c, colors_d,   ptr->numverts * 3 * sizeof(float), cudaMemcpyDeviceToHost);
  }

  cudaMemcpy(n, normals_d,  ptr->numverts * 3 * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(v, vertices_d, ptr->numverts * 3 * sizeof(float), cudaMemcpyDeviceToHost);
}


//*************************************************************
// draw a mesh consisting of vertices, facets, colors, normals etc.
void DispCmdTriMesh::cuda_putdata(const float * vertices_d,
                                  const float * normals_d,
                                  const unsigned char * colors_d,
                                  int num_facets,
                                  VMDDisplayList * dobj) {
  // make a triangle mesh (no strips)
  DispCmdTriMesh *ptr;
  if (colors_d == NULL) {
    ptr = (DispCmdTriMesh *)
                (dobj->append(DTRIMESH_C3F_N3F_V3F, sizeof(DispCmdTriMesh) +
                              sizeof(float) * num_facets * 3 * 6));
  } else {
    ptr = (DispCmdTriMesh *)
                (dobj->append(DTRIMESH_C4U_N3F_V3F, sizeof(DispCmdTriMesh) +
                              4 * sizeof(unsigned char) * num_facets * 3 +
                              sizeof(float) * num_facets * 3 * 6));
  }

  if (ptr == NULL)
    return;

  ptr->numverts=num_facets * 3;
  ptr->numfacets=num_facets;

  unsigned char *c=NULL;
  float *n=NULL, *v=NULL;
  if (colors_d == NULL) {
    ptr->pervertexcolors=0;
    ptr->getpointers(n, v);
  } else {
    ptr->pervertexcolors=1;
    ptr->getpointers(c, n, v);
    cudaMemcpy(c, colors_d,   ptr->numverts * 4 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
  }

  cudaMemcpy(n, normals_d,  ptr->numverts * 3 * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(v, vertices_d, ptr->numverts * 3 * sizeof(float), cudaMemcpyDeviceToHost);
}


//*************************************************************
// draw a mesh consisting of vertices, facets, colors, normals etc.
void DispCmdTriMesh::cuda_putdata(const float * vertices_d,
                                  const char * normals_d,
                                  const unsigned char * colors_d,
                                  int num_facets,
                                  VMDDisplayList * dobj) {
  // make a triangle mesh (no strips)
  DispCmdTriMesh *ptr;
  if (colors_d == NULL) {
    ptr = (DispCmdTriMesh *)
                (dobj->append(DTRIMESH_C4U_N3B_V3F, sizeof(DispCmdTriMesh) +
                              sizeof(char) * num_facets * 3 * 3 +
                              sizeof(float) * num_facets * 3 * 3));

  } else {
    ptr = (DispCmdTriMesh *)
                (dobj->append(DTRIMESH_C4U_N3B_V3F, sizeof(DispCmdTriMesh) +
                              4 * sizeof(unsigned char) * num_facets * 3 +
                              sizeof(char) * num_facets * 3 * 3 +
                              sizeof(float) * num_facets * 3 * 3));
  }

  if (ptr == NULL)
    return;

  ptr->numverts=num_facets * 3;
  ptr->numfacets=num_facets;

  unsigned char *c=NULL;
  signed char *n=NULL;
  float *v=NULL;
  if (colors_d == NULL) {
    ptr->pervertexcolors=0;
    ptr->getpointers(n, v);
  } else {
    ptr->pervertexcolors=1;
    ptr->getpointers(c, n, v);
    cudaMemcpy(c, colors_d,   ptr->numverts * 4 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
  }

  cudaMemcpy(n, normals_d,  ptr->numverts * 3 * sizeof(char), cudaMemcpyDeviceToHost);
  cudaMemcpy(v, vertices_d, ptr->numverts * 3 * sizeof(float), cudaMemcpyDeviceToHost);
}



