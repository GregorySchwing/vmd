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
 *	$RCSfile: DrawMolItemVolume.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.183 $	$Date: 2021/10/28 21:12:15 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   A continuation of rendering types from DrawMolItem
 *
 *   This file only contains representations for visualizing volumetric data
 ***************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "DrawMolItem.h"
#include "DrawMolecule.h"
#include "BaseMolecule.h" // need volume data definitions
#include "Inform.h"
#include "Isosurface.h"
#include "MoleculeList.h"
#include "Scene.h"
#include "VolumetricData.h"
#include "VMDApp.h"
#include "WKFUtils.h"

#define MYSGN(a) (((a) > 0) ? 1 : -1)

int DrawMolItem::draw_volume_get_colorid(void) {
  int colorid = 0;
  int catindex;

  switch (atomColor->method()) {
    case AtomColor::MOLECULE:
    case AtomColor::COLORID:
      colorid = atomColor->color_index();
      break;

    default:
      catindex = scene->category_index("Display");
      colorid = scene->category_item_value(catindex, "Foreground"); 
      break;
  }

  return colorid;
}


void DrawMolItem::draw_volume_box_solid(VolumetricData * v) {
  float v0[3], v1[3], v2[3];
  float vorigin[3], vxaxis[3], vyaxis[3], vzaxis[3];
  int usecolor;

  int i;
  for (i=0; i<3; i++) {
    vorigin[i] = float(v->origin[i]);
    vxaxis[i] = float(v->xaxis[i]);
    vyaxis[i] = float(v->yaxis[i]);
    vzaxis[i] = float(v->zaxis[i]);
  }

  append(DMATERIALON); // enable lighting and shading
  usecolor=draw_volume_get_colorid();
  cmdColorIndex.putdata(usecolor, cmdList);

  // Draw XY plane
  v0[0] = vorigin[0];
  v0[1] = vorigin[1];
  v0[2] = vorigin[2];

  v1[0] = v0[0] + vxaxis[0]; 
  v1[1] = v0[1] + vxaxis[1]; 
  v1[2] = v0[2] + vxaxis[2]; 

  v2[0] = v1[0] + vyaxis[0];
  v2[1] = v1[1] + vyaxis[1];
  v2[2] = v1[2] + vyaxis[2];
  cmdSquare.putdata(v2, v1, v0, cmdList);

  // Draw XZ plane
  v2[0] = v1[0] + vzaxis[0];
  v2[1] = v1[1] + vzaxis[1];
  v2[2] = v1[2] + vzaxis[2];
  cmdSquare.putdata(v0, v1, v2, cmdList);

  // Draw YZ plane
  v1[0] = v0[0] + vyaxis[0]; 
  v1[1] = v0[1] + vyaxis[1]; 
  v1[2] = v0[2] + vyaxis[2]; 

  v2[0] = v1[0] + vzaxis[0];
  v2[1] = v1[1] + vzaxis[1];
  v2[2] = v1[2] + vzaxis[2];
  cmdSquare.putdata(v2, v1, v0, cmdList);

  // Draw XY +Z plane
  v0[0] = vorigin[0] + vzaxis[0];
  v0[1] = vorigin[1] + vzaxis[1];
  v0[2] = vorigin[2] + vzaxis[2];

  v1[0] = v0[0] + vxaxis[0]; 
  v1[1] = v0[1] + vxaxis[1]; 
  v1[2] = v0[2] + vxaxis[2]; 

  v2[0] = v1[0] + vyaxis[0];
  v2[1] = v1[1] + vyaxis[1];
  v2[2] = v1[2] + vyaxis[2];
  cmdSquare.putdata(v0, v1, v2, cmdList);

  // Draw XZ +Y plane
  v0[0] = vorigin[0] + vyaxis[0];
  v0[1] = vorigin[1] + vyaxis[1];
  v0[2] = vorigin[2] + vyaxis[2];

  v1[0] = v0[0] + vxaxis[0]; 
  v1[1] = v0[1] + vxaxis[1]; 
  v1[2] = v0[2] + vxaxis[2]; 

  v2[0] = v1[0] + vzaxis[0];
  v2[1] = v1[1] + vzaxis[1];
  v2[2] = v1[2] + vzaxis[2];
  cmdSquare.putdata(v2, v1, v0, cmdList);

  // Draw YZ +X plane
  v0[0] = vorigin[0] + vxaxis[0];
  v0[1] = vorigin[1] + vxaxis[1];
  v0[2] = vorigin[2] + vxaxis[2];

  v1[0] = v0[0] + vyaxis[0]; 
  v1[1] = v0[1] + vyaxis[1]; 
  v1[2] = v0[2] + vyaxis[2]; 

  v2[0] = v1[0] + vzaxis[0];
  v2[1] = v1[1] + vzaxis[1];
  v2[2] = v1[2] + vzaxis[2];
  cmdSquare.putdata(v0, v1, v2, cmdList);

  draw_volume_box_lines(v);
}

void DrawMolItem::draw_volume_box_lines(VolumetricData * v) {
  float start[3], end[3];
  int xaxiscolor, yaxiscolor, zaxiscolor, axiscolor;
  int catindex;
  float vorigin[3], vxaxis[3], vyaxis[3], vzaxis[3];

  int i;
  for (i=0; i<3; i++) {
    vorigin[i] = float(v->origin[i]);
    vxaxis[i] = float(v->xaxis[i]);
    vyaxis[i] = float(v->yaxis[i]);
    vzaxis[i] = float(v->zaxis[i]);
  }

  // get colors from the scene pointer
  catindex = scene->category_index("Axes");
  xaxiscolor = scene->category_item_value(catindex, "X"); 
  yaxiscolor = scene->category_item_value(catindex, "Y"); 
  zaxiscolor = scene->category_item_value(catindex, "Z"); 

  catindex = scene->category_index("Display");
  axiscolor = scene->category_item_value(catindex, "Foreground"); 
 
  append(DMATERIALOFF);
  cmdLineType.putdata(SOLIDLINE, cmdList);
  cmdLineWidth.putdata(3, cmdList);

  start[0] = vorigin[0];
  start[1] = vorigin[1];
  start[2] = vorigin[2];

  // Draw X axis of volume box
  cmdColorIndex.putdata(xaxiscolor, cmdList);
  end[0] = start[0] + vxaxis[0];
  end[1] = start[1] + vxaxis[1];
  end[2] = start[2] + vxaxis[2];
  cmdLine.putdata(start, end, cmdList);

  // Draw Y axis of volume box
  cmdColorIndex.putdata(yaxiscolor, cmdList);
  end[0] = start[0] + vyaxis[0];
  end[1] = start[1] + vyaxis[1];
  end[2] = start[2] + vyaxis[2];
  cmdLine.putdata(start, end, cmdList);

  // Draw Z axis of volume box
  cmdColorIndex.putdata(zaxiscolor, cmdList);
  end[0] = start[0] + vzaxis[0];
  end[1] = start[1] + vzaxis[1];
  end[2] = start[2] + vzaxis[2];
  cmdLine.putdata(start, end, cmdList);

  // Draw remaining outline of volume box
  cmdLineWidth.putdata(1, cmdList);
  cmdColorIndex.putdata(axiscolor, cmdList);
 
  start[0] = vorigin[0] + vxaxis[0];    
  start[1] = vorigin[1] + vxaxis[1];    
  start[2] = vorigin[2] + vxaxis[2];    

  end[0] = start[0] + vyaxis[0];
  end[1] = start[1] + vyaxis[1];
  end[2] = start[2] + vyaxis[2];
  cmdLine.putdata(start, end, cmdList);

  end[0] = start[0] + vzaxis[0];
  end[1] = start[1] + vzaxis[1];
  end[2] = start[2] + vzaxis[2];
  cmdLine.putdata(start, end, cmdList);

  start[0] = vorigin[0] + vyaxis[0];    
  start[1] = vorigin[1] + vyaxis[1];    
  start[2] = vorigin[2] + vyaxis[2];    

  end[0] = start[0] + vxaxis[0];
  end[1] = start[1] + vxaxis[1];
  end[2] = start[2] + vxaxis[2];
  cmdLine.putdata(start, end, cmdList);

  end[0] = start[0] + vzaxis[0];
  end[1] = start[1] + vzaxis[1];
  end[2] = start[2] + vzaxis[2];
  cmdLine.putdata(start, end, cmdList);


  start[0] = vorigin[0] + vzaxis[0];    
  start[1] = vorigin[1] + vzaxis[1];    
  start[2] = vorigin[2] + vzaxis[2];    

  end[0] = start[0] + vxaxis[0];
  end[1] = start[1] + vxaxis[1];
  end[2] = start[2] + vxaxis[2];
  cmdLine.putdata(start, end, cmdList);

  end[0] = start[0] + vyaxis[0];
  end[1] = start[1] + vyaxis[1];
  end[2] = start[2] + vyaxis[2];
  cmdLine.putdata(start, end, cmdList);
 

  start[0] = vorigin[0] + vxaxis[0] + vyaxis[0] + vzaxis[0];
  start[1] = vorigin[1] + vxaxis[1] + vyaxis[1] + vzaxis[1];
  start[2] = vorigin[2] + vxaxis[2] + vyaxis[2] + vzaxis[2];
 
  end[0] = start[0] - vxaxis[0];
  end[1] = start[1] - vxaxis[1];
  end[2] = start[2] - vxaxis[2];
  cmdLine.putdata(start, end, cmdList);

  end[0] = start[0] - vyaxis[0];
  end[1] = start[1] - vyaxis[1];
  end[2] = start[2] - vyaxis[2];
  cmdLine.putdata(start, end, cmdList);

  end[0] = start[0] - vzaxis[0];
  end[1] = start[1] - vzaxis[1];
  end[2] = start[2] - vzaxis[2];
  cmdLine.putdata(start, end, cmdList);
}


void DrawMolItem::draw_volume_isosurface_points(const VolumetricData * v,
                                float isovalue, int stepsize, int thickness) {
  int x,y,z;
  float xax[3], yax[3], zax[3];
  int usecolor;
  ResizeArray<float> centers;
  ResizeArray<float> colors;

  int i;
  float vorigin[3];
  for (i=0; i<3; i++) {
    vorigin[i] = float(v->origin[i]);
  }

  // calculate cell axes
  v->cell_axes(xax, yax, zax);
  ptrdiff_t zplanesz = v->xsize * v->ysize;
 
  append(DMATERIALOFF);
  usecolor = draw_volume_get_colorid();
  const float *cp = scene->color_value(usecolor);
  cmdColorIndex.putdata(usecolor, cmdList);

  int pointcount = 0;
  for (z=0; z<v->zsize; z+=stepsize) {
    // Pre-extend our output buffers so we don't have to call append().
    // As long as the output buffers are large enough to hold a full 
    // additional plane of points, we can avoid individually appending points.
    // Since this resizing scheme bypasses the normal bookkeeping done by
    // append() methods, we are responsible for setting the number of 
    // occupied array elements via set_size(), such that the results 
    // returned by the num() method are correct.
    ptrdiff_t newsz = (pointcount + zplanesz) * 3L;
    if (newsz >= centers.num()) {
      colors.extend(newsz);
      centers.extend(newsz);
    }

    // set starting position of output pointers
    int idx = pointcount * 3;
    float *pos=&centers[0];
    float *cols=&colors[0];
    for (y=0; y<v->ysize; y+=stepsize) {
      ptrdiff_t rowidx = z*zplanesz + y*v->xsize;
      const float *rowaddr = &(v->data[rowidx]);  

      // compute starting coordinate of row...
      float rowx = vorigin[0] + y * yax[0] + z * zax[0];
      float rowy = vorigin[1] + y * yax[1] + z * zax[1];
      float rowz = vorigin[2] + y * yax[2] + z * zax[2];

      // loop through xsize - 1 rather than the full range
      int lastgrid = v->xsize - 1;
      for (x=0; x<lastgrid; x+=stepsize) {
        // draw a point if the isovalue falls between neighboring X samples
        // minimal test criterion: 
        //   both terms have to be true, or both false, indicating that
        //   the isovalue lies in between the two nearest voxel values
        if (!((rowaddr[x] <= isovalue) ^ (isovalue < rowaddr[x+1]))) {
          // draw a point there
          pos[idx  ] = rowx + x * xax[0];
          pos[idx+1] = rowy + x * xax[1];
          pos[idx+2] = rowz + x * xax[2];

          // XXX per-vertex colors are redundant here and we could get a 2x 
          //     draw perf gain here by eliminating them...
          cols[idx  ] = cp[0];
          cols[idx+1] = cp[1];
          cols[idx+2] = cp[2];

          idx+=3;
          pointcount++;
        }
      } 
    } 

    // sync up ResizeArray occupied array size
    centers.set_size(idx);
    colors.set_size(idx);

    // Avoid integer overflow by submitting geometry batches after 
    // every 10M primitives, well before the point of int32 wraparound.
    if (pointcount > 10000000) {
      cmdPointArray.putdata((float *) &centers[0],
                            (float *) &colors[0],
                            (float) thickness,
                            pointcount,
                            cmdList);

      // reset tmp buffers
      centers.clear();
      colors.clear(); // XXX get rid of per-vertex colors
      pointcount=0;
    }
  }

  if (pointcount > 0) {
    cmdPointArray.putdata((float *) &centers[0],
                          (float *) &colors[0],
                          (float) thickness,
                          pointcount,
                          cmdList);
  }
}


void DrawMolItem::draw_volume_isosurface_lit_points(VolumetricData * v, 
                                 float isovalue, int stepsize, int thickness) {
  int x,y,z;
  float xax[3], yax[3], zax[3];
  ResizeArray<float> centers;
  ResizeArray<float> normals;
  ResizeArray<float> colors;

  int i;
  float vorigin[3];
  for (i=0; i<3; i++) {
    vorigin[i] = float(v->origin[i]);
  }

  // calculate cell axes
  v->cell_axes(xax, yax, zax);
  ptrdiff_t zplanesz = v->xsize * v->ysize;

  // get direct access to gradient data for speed, and force
  // generation of the gradients if they don't already exist
  const float *volgradient = v->access_volume_gradient();

  append(DMATERIALON);
  int usecolor = draw_volume_get_colorid();
  const float *cp = scene->color_value(usecolor);
  cmdColorIndex.putdata(usecolor, cmdList);

  int pointcount = 0;
  for (z=0; z<v->zsize; z+=stepsize) {
    // Pre-extend our output buffers so we don't have to call append().
    // As long as the output buffers are large enough to hold a full
    // additional plane of points, we can avoid individually appending points.
    // Since this resizing scheme bypasses the normal bookkeeping done by
    // append() methods, we are responsible for setting the number of 
    // occupied array elements via set_size(), such that the results 
    // returned by the num() method are correct.
    ptrdiff_t newsz = (pointcount + zplanesz) * 3L;
    if ((pointcount + zplanesz)*3L >= centers.num()) {
      colors.extend(newsz);
      normals.extend(newsz);
      centers.extend(newsz);
    }

    // set starting position of output pointers
    int idx = pointcount * 3;
    float *pos=&centers[0];
    float *cols=&colors[0];
    float *norms=&normals[0]; 
    for (y=0; y<v->ysize; y+=stepsize) {
      ptrdiff_t rowidx = z*zplanesz + y*v->xsize;
      float *rowaddr = &(v->data[rowidx]);  

      // compute starting coordinate of row...
      float rowx = vorigin[0] + y * yax[0] + z * zax[0];
      float rowy = vorigin[1] + y * yax[1] + z * zax[1];
      float rowz = vorigin[2] + y * yax[2] + z * zax[2];

      // loop through xsize - 1 rather than the full range
      int lastgrid = v->xsize - 1;
      for (x=0; x<lastgrid; x+=stepsize) {
        // draw a point if the isovalue falls between neighboring X samples
        // minimal test criterion:
        //   both terms have to be true, or both false, indicating that
        //   the isovalue lies in between the two nearest voxel values
        if (!((rowaddr[x] <= isovalue) ^ (isovalue < rowaddr[x+1]))) {
          // draw a point there
          pos[idx  ] = rowx + x * xax[0];
          pos[idx+1] = rowy + x * xax[1];
          pos[idx+2] = rowz + x * xax[2];

          ptrdiff_t volgradidx = (rowidx + x) * 3L;
          norms[idx  ] = volgradient[volgradidx    ];
          norms[idx+1] = volgradient[volgradidx + 1];
          norms[idx+2] = volgradient[volgradidx + 2];
 
          // XXX per-vertex colors are redundant here and we could get a 1.3x 
          //     draw perf gain here by eliminating them...
          cols[idx  ] = cp[0];
          cols[idx+1] = cp[1];
          cols[idx+2] = cp[2];

          idx+=3;
          pointcount++;
        }
      } 
    } 

    // sync up ResizeArray occupied array size
    centers.set_size(idx);
    normals.set_size(idx);
    colors.set_size(idx);

    // Avoid integer overflow by submitting geometry batches after 
    // every 10M primitives, well before the point of int32 wraparound.
    if (pointcount > 10000000) {
      DispCmdLitPointArray cmdLitPointArray;
      cmdLitPointArray.putdata((float *) &centers[0],
                               (float *) &normals[0],
                               (float *) &colors[0],
                               (float) thickness,
                               pointcount,
                               cmdList);

      // reset tmp buffers
      centers.clear();
      normals.clear();
      colors.clear(); // XXX get rid of per-vertex colors
      pointcount=0;
    }
  }

  if (pointcount > 0) {
    DispCmdLitPointArray cmdLitPointArray;
    cmdLitPointArray.putdata((float *) &centers[0],
                             (float *) &normals[0],
                             (float *) &colors[0],
                             (float) thickness,
                             pointcount,
                             cmdList);
  }
}


void DrawMolItem::draw_volume_isosurface_lines(VolumetricData * v, 
                                 float isovalue, int stepsize, int thickness) {
  IsoSurface *s = new IsoSurface;

  // XXX We must handle multi-billion triangle isosurface meshes -- needs work
  // Now that VMD can load multi-gigavoxel volumes it isn't a safe assumption
  // that we can do a single-pass MC isosurface extraction, as this could
  // result in billions of triangles, and that creates problems for triangle
  // mesh indexing on GPUs that can't do 64-bit indices.  So this code must 
  // now decompose the MC isosurface calculation into multiple passes that
  // produce less than 2 billion triangles per pass.  Since these are 
  // inherently data-parallel, there's also an opportunity to pick up speed
  // for large volumes with the use of a persistent CPU thread pool and 
  // dicing the problem up into multiple segments that won't break indexing
  // and allow us to exploit parallelism.
  s->clear();
  s->compute(v, isovalue, stepsize);

  if (s->numtriangles > 0) {
    append(DMATERIALOFF); // enable lighting and shading
    cmdLineType.putdata(SOLIDLINE, cmdList);
    cmdLineWidth.putdata(thickness, cmdList);

    int usecolor = draw_volume_get_colorid();
    cmdColorIndex.putdata(usecolor, cmdList);

    // draw triangles
    ptrdiff_t i;
    for (i=0; i<s->numtriangles; i++) {
      float *addr;
      addr = &(s->v[i * 9]); 
      cmdLine.putdata(&addr[0], &addr[3], cmdList);
      cmdLine.putdata(&addr[3], &addr[6], cmdList);
      cmdLine.putdata(&addr[6], &addr[0], cmdList);
    }
  }

  delete s; // we don't need this stuff after this point 
}



void DrawMolItem::draw_volume_isosurface_trimesh(VolumetricData * v, 
                                                 float isovalue, int stepsize,
                                                 const float *voltex) {
  IsoSurface s;

  // XXX We must handle multi-billion triangle isosurface meshes -- needs work
  // Now that VMD can load multi-gigavoxel volumes it isn't a safe assumption
  // that we can do a single-pass MC isosurface extraction, as this could
  // result in billions of triangles, and that creates problems for triangle
  // mesh indexing on GPUs that can't do 64-bit indices.  So this code must 
  // now decompose the MC isosurface calculation into multiple passes that
  // produce less than 2 billion triangles per pass.  Since these are 
  // inherently data-parallel, there's also an opportunity to pick up speed
  // for large volumes with the use of a persistent CPU thread pool and 
  // dicing the problem up into multiple segments that won't break indexing
  // and allow us to exploit parallelism.
  s.clear();                 // initialize isosurface data
  s.compute(v, isovalue, stepsize); // compute the isosurface
  s.vertexfusion(36, 36);    // identify and eliminate duplicated vertices
  s.normalize();             // normalize interpolated gradient/surface normals

#if 1
  if (s.numtriangles > 0) {
    append(DMATERIALON); // enable lighting and shading
    int usecolor = draw_volume_get_colorid();
    cmdColorIndex.putdata(usecolor, cmdList);

    if (voltex != NULL) {
      // assign per-vertex colors by a 3-D texture map
      s.set_color_voltex_rgb3fv(voltex); 
    } else {
      // use a single color for the entire mesh
      s.set_color_rgb3fv(scene->color_value(usecolor)); 
    }

    // Create a triangle mesh
    // XXX don't try to stripify it since this triggers a crash in ACTC for
    //     unknown reasons
    // XXX ACTC is undoubtably unable to deal with billion-triangle meshes
    // XXX The DispCmd code won't yet scale past 2 billion triangle meshes
    cmdTriMesh.putdata(&s.v[0], &s.n[0], &s.c[0], s.v.num() / 3,
                       &s.f[0], s.numtriangles, 0, cmdList);
  }
#else
  if (s.numtriangles > 0) {
    append(DMATERIALON); // enable lighting and shading
    int usecolor = draw_volume_get_colorid();
    cmdColorIndex.putdata(usecolor, cmdList);

    // draw surface with per-vertex normals using a vertex array
    float *c = new float[s.numtriangles * 9L];
    const float *fp = scene->color_value(usecolor);
    ptrdiff_t i;
    for (i=0; i<s.numtriangles; i++) { 
      int ind = i * 9;

      c[ind    ] = fp[0]; // Red
      c[ind + 1] = fp[1]; // Green
      c[ind + 2] = fp[2]; // Blue

      ind+=3;
      c[ind    ] = fp[0]; // Red
      c[ind + 1] = fp[1]; // Green
      c[ind + 2] = fp[2]; // Blue

      ind+=3;
      c[ind    ] = fp[0]; // Red
      c[ind + 1] = fp[1]; // Green
      c[ind + 2] = fp[2]; // Blue
    }

    // Create a triangle mesh
    // XXX don't try to stripify it since this triggers a crash in ACTC for
    //     unknown reasons
    // XXX ACTC is undoubtably unable to deal with billion-triangle meshes
    // XXX The DispCmd code won't yet scale past 2 billion triangle meshes
    cmdTriMesh.putdata(&s.v[0], &s.n[0], c, s.v.num() / 3,
                       &s.f[0], s.numtriangles, 0, cmdList);

    delete [] c;
  }
#endif
}


// recompute sliceTextureCoord, sliceVertexCoord, and sliceNormal
// using the given volumetric data set, axis, and offset
// sliceAxis and sliceOffset.  
static void prepare_texture_coordinates(const VolumetricData *v, 
    float loc, int axis, float *sliceNormal, float *sliceTextureCoords,
    float *sliceVertexes) {

  float t0[3], t1[3], t2[3], t3[3];
  float v0[3], v1[3], v2[3], v3[3];
  float normal[3];

  float vorigin[3], vxaxis[3], vyaxis[3], vzaxis[3];
  int i;
  for (i=0; i<3; i++) {
    vorigin[i] = float(v->origin[i]);
    vxaxis[i] = float(v->xaxis[i]);
    vyaxis[i] = float(v->yaxis[i]);
    vzaxis[i] = float(v->zaxis[i]);
  }


  if (loc < 0.0f)
      loc = 0.0f;

  if (loc > 1.0f)
      loc = 1.0f;

  switch (axis) {
    // X-Axis
    case 0:
    default:
      t0[0] = loc;
      t0[1] = 0.0f;
      t0[2] = 0.0f;

      t1[0] = loc;
      t1[1] = 0.0f;
      t1[2] = 1.0f;

      t2[0] = loc;
      t2[1] = 1.0f;
      t2[2] = 1.0f;

      t3[0] = loc;
      t3[1] = 1.0f;
      t3[2] = 0.0f;

      v0[0] = vorigin[0] + (vxaxis[0] * loc);
      v0[1] = vorigin[1] + (vxaxis[1] * loc);
      v0[2] = vorigin[2] + (vxaxis[2] * loc);

      v1[0] = v0[0] + vzaxis[0];
      v1[1] = v0[1] + vzaxis[1];
      v1[2] = v0[2] + vzaxis[2];

      v2[0] = v0[0] + vzaxis[0] + vyaxis[0];
      v2[1] = v0[1] + vzaxis[1] + vyaxis[1];
      v2[2] = v0[2] + vzaxis[2] + vyaxis[2];

      v3[0] = v0[0]             + vyaxis[0];
      v3[1] = v0[1]             + vyaxis[1];
      v3[2] = v0[2]             + vyaxis[2];

      normal[0] = vxaxis[0];
      normal[1] = vxaxis[1];
      normal[2] = vxaxis[2];
      vec_normalize(&normal[0]);
      break;

    // Y-Axis
    case 1:
      t0[0] = 0.0f;
      t0[1] = loc;
      t0[2] = 0.0f;

      t1[0] = 1.0f;
      t1[1] = loc;
      t1[2] = 0.0f;

      t2[0] = 1.0f;
      t2[1] = loc;
      t2[2] = 1.0f;

      t3[0] = 0.0f;
      t3[1] = loc;
      t3[2] = 1.0f;

      v0[0] = vorigin[0] + (vyaxis[0] * loc);
      v0[1] = vorigin[1] + (vyaxis[1] * loc);
      v0[2] = vorigin[2] + (vyaxis[2] * loc);

      v1[0] = v0[0] + vxaxis[0];
      v1[1] = v0[1] + vxaxis[1];
      v1[2] = v0[2] + vxaxis[2];

      v2[0] = v0[0] + vxaxis[0] + vzaxis[0];
      v2[1] = v0[1] + vxaxis[1] + vzaxis[1];
      v2[2] = v0[2] + vxaxis[2] + vzaxis[2];

      v3[0] = v0[0]             + vzaxis[0];
      v3[1] = v0[1]             + vzaxis[1];
      v3[2] = v0[2]             + vzaxis[2];

      normal[0] = vyaxis[0];
      normal[1] = vyaxis[1];
      normal[2] = vyaxis[2];
      vec_normalize(&normal[0]);
      break;

    // Z-Axis
    case 2:
      t0[0] = 0.0f;
      t0[1] = 0.0f;
      t0[2] = loc;

      t1[0] = 1.0f;
      t1[1] = 0.0f;
      t1[2] = loc;

      t2[0] = 1.0f;
      t2[1] = 1.0f;
      t2[2] = loc;

      t3[0] = 0.0f;
      t3[1] = 1.0f;
      t3[2] = loc;

      v0[0] = vorigin[0] + (vzaxis[0] * loc);
      v0[1] = vorigin[1] + (vzaxis[1] * loc);
      v0[2] = vorigin[2] + (vzaxis[2] * loc);

      v1[0] = v0[0] + vxaxis[0];
      v1[1] = v0[1] + vxaxis[1];
      v1[2] = v0[2] + vxaxis[2];

      v2[0] = v0[0] + vxaxis[0] + vyaxis[0];
      v2[1] = v0[1] + vxaxis[1] + vyaxis[1];
      v2[2] = v0[2] + vxaxis[2] + vyaxis[2];

      v3[0] = v0[0]             + vyaxis[0];
      v3[1] = v0[1]             + vyaxis[1];
      v3[2] = v0[2]             + vyaxis[2];

      normal[0] = vzaxis[0];
      normal[1] = vzaxis[1];
      normal[2] = vzaxis[2];
      vec_normalize(&normal[0]);
      break;
  }

  vec_copy(sliceTextureCoords  , t0);
  vec_copy(sliceTextureCoords+3, t1);
  vec_copy(sliceTextureCoords+6, t2);
  vec_copy(sliceTextureCoords+9, t3);
  vec_copy(sliceVertexes  , v0);
  vec_copy(sliceVertexes+3, v1);
  vec_copy(sliceVertexes+6, v2);
  vec_copy(sliceVertexes+9, v3);
  vec_copy(sliceNormal, normal);
}

void DrawMolItem::updateVolumeTexture() {
  float vmin, vmax;
  int volid = atomColor->volume_index();
  if (atomRep->method()==AtomRep::VOLSLICE) {
    // the volslice rep has its own volid specification
    volid = (int)atomRep->get_data(AtomRep::SPHERERES);
  }

  VolumetricData *v = mol->modify_volume_data(volid);
  if (v == NULL) {
    msgInfo << "No volume data loaded at index " << volid << sendmsg;
    return;
  } 

  // determine if a new 3D texture needs to be prepared and uploaded.  This
  // occurs if:
  // (1) the coloring method has changed
  // (2) the choice of volumetric data set has changed.
  atomColor->get_colorscale_minmax(&vmin, &vmax);
  if (!vmin && !vmax) {
    v->datarange(vmin, vmax);
  }

  if (volumeTexture.getTextureMap() && !(needRegenerate & COL_REGEN) &&
      volid == voltexVolid && atomColor->method() == voltexColorMethod &&
      voltexDataMin == vmin && voltexDataMax == vmax) {
    // nothing to do
    return;
  }

  voltexColorMethod = atomColor->method();
  voltexVolid = volid;
  voltexDataMin = vmin; 
  voltexDataMax = vmax;

  volumeTexture.setGridData(v);

  // update the volumeTexture instance
  switch (atomColor->method()) {
    case AtomColor::POS:
    case AtomColor::POSX:
    case AtomColor::POSY:
    case AtomColor::POSZ:
      volumeTexture.generatePosTexture();
      break;

    case AtomColor::INDEX:
      volumeTexture.generateIndexTexture();
      break;

    case AtomColor::CHARGE:
      volumeTexture.generateChargeTexture(vmin, vmax);
      break;

    case AtomColor::VOLUME:
      volumeTexture.generateColorScaleTexture(vmin, vmax, scene);
      break;

#if 0
    case AtomColor::STRUCTURE:
      // 3-D Contour lines
      volumeTexture.generateContourLineTexture(0.5, 0.5);
      break;
#endif
      
    case AtomColor::NAME:
    default:
      // HSV color ramp
      volumeTexture.generateHSVTexture(vmin, vmax);
      break;
  }
}


void DrawMolItem::draw_volslice(int volid, float slice, int axis, int texmode) {
  const VolumetricData *v = mol->modify_volume_data(volid);
  float sliceNormal[3];
  float sliceTextureCoords[12];
  float sliceVertexes[12];
  if (!volumeTexture.getTextureMap()) {
    msgErr << "draw_volslice: no texture map has been generated!" << sendmsg;
    return;
  }

  sprintf(commentBuffer, "Mol[%d] Rep[%d] VolSlice", mol->id(), repNumber);
  cmdCommentX.putdata(commentBuffer, cmdList);

  prepare_texture_coordinates(v, slice, axis, sliceNormal, sliceTextureCoords,
                              sliceVertexes);
  
  // Rescale the texture coordinates so that they include just the 
  // part of the texture map to which we mapped our volume data.
  // Clamp range to just below 1.0 to prevent
  // textures from being clamped to the
  // border color on some video cards.
  float tscale[3] = { float(v->xsize), float(v->ysize), float(v->zsize) };
  const int *size = volumeTexture.getTextureSize();
  for (int i=0; i<3; i++) {
    float rescale = (tscale[i] / (float)size[i]) * 0.99999f;
    sliceTextureCoords[i  ] *= rescale;
    sliceTextureCoords[i+3] *= rescale;
    sliceTextureCoords[i+6] *= rescale;
    sliceTextureCoords[i+9] *= rescale;
  }
  // add command to draw the slice itself.
  append(DMATERIALON); // enable lighting and shading

  // Pass instructions for the slice itself.
  cmdVolSlice.putdata(texmode, sliceNormal, sliceVertexes, sliceTextureCoords, cmdList);
}


void DrawMolItem::draw_isosurface(int volid, float isovalue, int drawbox, int style, int stepsize, int thickness) {
  VolumetricData * v = NULL;

  v = mol->modify_volume_data(volid);
  if (v == NULL) {
    msgInfo << "No volume data loaded at index " << volid << sendmsg;
    return;
  } 

  sprintf(commentBuffer, "Mol[%d] Rep[%d] Isosurface", mol->id(), repNumber);
  cmdCommentX.putdata(commentBuffer, cmdList);

  // Safety checks to prevent stepsize from cratering surface extraction code
  if (stepsize >= v->xsize)
    stepsize = (v->xsize - 1);
  if (stepsize >= v->ysize)
    stepsize = (v->ysize - 1);
  if (stepsize >= v->zsize)
    stepsize = (v->zsize - 1);
  if (stepsize < 2)
    stepsize = 1;

  if (drawbox > 0) {
    // don't texture the box if color by volume is active
    if (atomColor->method() == AtomColor::VOLUME) {
      append(DVOLTEXOFF);
    }
    // wireframe only?  or solid?
    if (style > 0 || drawbox == 2) {
      draw_volume_box_lines(v);
    } else {
      draw_volume_box_solid(v);
    }
    if (atomColor->method() == AtomColor::VOLUME) {
      append(DVOLTEXON);
    }
  } 

  if ((drawbox == 2) || (drawbox == 0)) {
    switch (style) {
      case 3:
        // shaded points isosurface looping over X-axis, 1 point max per voxel
        draw_volume_isosurface_lit_points(v, isovalue, stepsize, thickness);
        break;

      case 2:
        // points isosurface looping over X-axis, 1 point max per voxel
        draw_volume_isosurface_points(v, isovalue, stepsize, thickness);
        break;

      case 1:
        // lines implementation, max of 18 line per voxel (3-per triangle)
        draw_volume_isosurface_lines(v, isovalue, stepsize, thickness);
        break;

      case 0:
      default:
        // trimesh polygonalized surface, max of 6 triangles per voxel
        // XXX to enable per-vertex colors rather than 3-D texture mapping,
        //     we need to expose a setting for this purpose, for use when
        //     people need to export to packages incapable of 3-D texturing.
        draw_volume_isosurface_trimesh(v, isovalue, stepsize);
        break;
    }
  }
}


#if 0
// draw a 2-D contour line in a slice plane
void DrawMolItem::draw_volslice_contour_lines(int volid, float slice, int axis) {
  const VolumetricData * v = NULL;

  v = mol->get_volume_data(volid);

  if (v == NULL) {
    msgInfo << "No volume data loaded at index " << volid << sendmsg;
    return;
  }

  msgInfo << "Placeholder for contour slice not implemented yet..."
}
#endif


// calculate seed voxels for field lines based on gradient magnitude
int DrawMolItem::calcseeds_grid(VolumetricData * v, ResizeArray<float> *seeds, int maxseedcount) {
  int i;
  float vorigin[3];
  for (i=0; i<3; i++) {
    vorigin[i] = float(v->origin[i]);
  }

  // calculate cell axes
  float xax[3], yax[3], zax[3];
  v->cell_axes(xax, yax, zax);

  int seedcount = maxseedcount+1; // force loop to run once
  int stepsize = 1;

  // get direct access to gradient data for speed, and force
  // generation of the gradients if they don't already exist
  const float *volgradient = v->access_volume_gradient();

  // iterate if we generate more seeds than we can really use
  while (seedcount > maxseedcount) {
    seedcount=0;
    seeds->clear();

    int x,y,z;
    float pos[3];
    ptrdiff_t zplanesz = v->xsize*v->ysize;
    for (z=0; z<v->zsize && seedcount < maxseedcount; z+=stepsize) {
      for (y=0; y<v->ysize; y+=stepsize) {
        ptrdiff_t rowidx = z*zplanesz + y*v->xsize;

        // compute starting coordinate of row...
        float rowx = vorigin[0] + y * yax[0] + z * zax[0];
        float rowy = vorigin[1] + y * yax[1] + z * zax[1];
        float rowz = vorigin[2] + y * yax[2] + z * zax[2];

        for (x=0; x<v->xsize; x+=stepsize) {
          float grad[3];
          vec_copy(grad, &volgradient[(rowidx + x) * 3L]);

          pos[0] = rowx + x * xax[0];
          pos[1] = rowy + x * xax[1];
          pos[2] = rowz + x * xax[2];
  
          seedcount++;
          seeds->append3(&pos[0]);
        }
      }
    }
    stepsize++; // increase stepsize if we have to try again
  }
  
  return seedcount;
}


// calculate seed voxels for field lines based on gradient magnitude
int DrawMolItem::calcseeds_gradient_magnitude(VolumetricData * v, ResizeArray<float> *seeds, float seedmin, float seedmax, int maxseedcount) {
  float seedmin2 = seedmin*seedmin;
  float seedmax2 = seedmax*seedmax;

  int i;
  float vorigin[3];
  for (i=0; i<3; i++) {
    vorigin[i] = float(v->origin[i]);
  }

  // calculate cell axes
  float xax[3], yax[3], zax[3];
  v->cell_axes(xax, yax, zax);

  int seedcount = maxseedcount+1; // force loop to run once
  int stepsize = 1;

  // get direct access to gradient data for speed, and force
  // generation of the gradients if they don't already exist
  const float *volgradient = v->access_volume_gradient();

  // iterate if we generate more seeds than we can really use
  while (seedcount > maxseedcount) {
    seedcount=0;
    seeds->clear();

    int x,y,z;
    float pos[3];
    ptrdiff_t zplanesz = v->xsize*v->ysize;
    for (z=0; z<v->zsize && seedcount < maxseedcount; z+=stepsize) {
      for (y=0; y<v->ysize; y+=stepsize) {
        ptrdiff_t rowidx = z*zplanesz + y*v->xsize;

        // compute starting coordinate of row...
        float rowx = vorigin[0] + y * yax[0] + z * zax[0];
        float rowy = vorigin[1] + y * yax[1] + z * zax[1];
        float rowz = vorigin[2] + y * yax[2] + z * zax[2];

        for (x=0; x<v->xsize; x+=stepsize) {
          float grad[3];
          float gradmag2;
          vec_copy(grad, &volgradient[(rowidx + x) * 3L]);
          gradmag2 = dot_prod(grad, grad);

          if ((gradmag2 <= seedmax2) &&
              (gradmag2 >= seedmin2)) {
            pos[0] = rowx + x * xax[0];
            pos[1] = rowy + x * xax[1];
            pos[2] = rowz + x * xax[2];
  
            seedcount++;
            seeds->append3(&pos[0]);
          }
        }
      }
    }
    stepsize++; // increase stepsize if we have to try again
  }
  
  return seedcount;
}


// draw a 3-D field lines that follow the volume gradient
void DrawMolItem::draw_volume_field_lines(int volid, int seedusegrid, int maxseeds, 
                                          float seedval, float deltacell, 
                                          float minlen, float maxlen, 
                                          int drawstyle, int tuberes, float thickness) {
  VolumetricData * v = NULL;
  v = mol->modify_volume_data(volid);
  int printdonemesg=0;

  if (v == NULL) {
    msgInfo << "No volume data loaded at index " << volid << sendmsg;
    return;
  }

  int seedcount = 0;
  int pointcount = 0;
  int totalpointcount = 0;
  int usecolor;
  ResizeArray<float> seeds;

  sprintf(commentBuffer,"Mol[%d] Rep[%d] FieldLines", mol->id(), repNumber);
  cmdCommentX.putdata(commentBuffer, cmdList);


  DispCmdSphereRes cmdSphereRes;
  if (drawstyle != 0) {
    cmdSphereRes.putdata(tuberes, cmdList);
    append(DMATERIALON); // enable lighting and shading
    thickness *= 0.05f;  // XXX hack until we have a better GUI
  } else {
    append(DMATERIALOFF);
  }

  usecolor = draw_volume_get_colorid();
  cmdColorIndex.putdata(usecolor, cmdList);

  if (seedusegrid)
    seedcount = calcseeds_grid(v, &seeds, maxseeds);
  else
    seedcount = calcseeds_gradient_magnitude(v, &seeds, seedval*0.5f, seedval*1.5f, maxseeds);

  // Integrate field lines starting with each of the seeds to simulate
  // particle advection.
  // Uses Euler's approximation for solving the initial value problem.
  // We could get a more accurate solution using a fourth order Runge-Kutta
  // method, but with more math per iteration.  We may want to implement 
  // the integrator as a user selected option.

  // The choice of integration step size is currently arbitrary,
  // but will become a user-defined parameter, since it affects speed
  // and accuracy.  A good default might be 0.25 times the smallest
  // grid cell spacing axis.
  float lx, ly, lz;
  v->cell_lengths(&lx, &ly, &lz);
  float mincelllen=lx;
  mincelllen = (mincelllen < ly) ? mincelllen : ly;
  mincelllen = (mincelllen < lz) ? mincelllen : lz;
  float delta=mincelllen * deltacell; // delta per step (compensates gradient magnitude)

  // minimum gradient magnitude, before we consider that we've found
  // a critical point in the dataset.
  float mingmag =  0.0001f;

  // max gradient magnitude, before we consider it a source/sink
  float maxgmag = 5;

  ResizeArray<float> points;
  
  // ensure that the volume gradient has been computed prior to
  // sampling it for field line construction... (discard pointer)
  v->access_volume_gradient();

  // For each seed point, integrate in both positive and
  // negative directions for a field line length up to
  // the maxlen criterion.
  wkfmsgtimer *msgt = wkf_msg_timer_create(1);
  int seed;
  for (seed=0; seed < seedcount; seed++) {
    // emit UI messages as integrator runs, for long calculations...
    if (!(seed & 7) && wkf_msg_timer_timeout(msgt)) {
      char tmpbuf[128] = { 0 };
      sprintf(tmpbuf, "%6.2f %% complete", (100.0f * seed) / (float) seedcount);
      msgInfo << "integrating " << seedcount << " field lines: " << tmpbuf << sendmsg;
      printdonemesg=1;
    }
 
    int direction;
    for (direction=-1; direction != 1; direction=1) {
      float pos[3], comsum[3];
      vec_copy(pos, &seeds[seed*3]); // integration starting point is the seed

      // init the arrays
      points.clear();

      // main integration loop
      pointcount=0;
      totalpointcount++;
      float len=0;
      int iterations=0;
      float dir = (float) direction;

      vec_zero(comsum); // clear center of mass accumulator

      while ((len<maxlen) && (totalpointcount < 100000000)) {
        float grad[3];

        // sample gradient at the current position
        v->voxel_gradient_interpolate_from_coord(pos, grad);

        // Early-exit if we run out of bounds (gradient returned will
        // be a vector of NANs), run into a critical point (zero gradient)
        // or a huge gradient at a source/sink point in the dataset.
        // Since IEEE FP defines that all tests against NaN should return
        // false, we invert the the bound tests with so that a comparison
        // with a gmag of NaN will early-exit the advection loop.
        float gmag = norm(grad);
        if (!(gmag >= mingmag && gmag <= maxgmag))
           break;

        // Draw the current point only after the gradient value
        // has been checked, so we don't end up with out-of-bounds
        // vertices.
        // Only emit a fraction of integration points for display since
        // the integrator stepsize needs to be small for more numerical
        // accuracy, but the field lines themselves can be well 
        // represented with fewer sample points.
        if (!(iterations & 1)) {
          // Add a vertex for this field line
          points.append3(&pos[0]);

          vec_incr(comsum, pos);
 
          pointcount++;
          totalpointcount++;
        }

        // adjust integration stepsize so we never move more than 
        // the distance specified by delta at each step, to compensate
        // for varying gradient magnitude
        vec_scaled_add(pos, dir * delta / gmag, grad); // integrate position
        len += delta; // accumulate distance

        iterations++;
      }

      int drawfieldline = 1;

      // only draw the field line for this seed if we have enough points.
      // If we haven't reached the minimum field line length, we'll
      // drop the whole field line.
      if (pointcount < 2 || len < minlen)
        drawfieldline = 0;

      // only draw if bounding sphere diameter exceeds minlen
      if (drawfieldline) {       
        float com[3];
        vec_scale(com, 1.0f / (float) pointcount, comsum);
        float minlen2 = minlen*minlen;

        drawfieldline = 0;
        int p;
        for (p=0; p<pointcount; p++) {
          if ((2.0f * distance2(com, &points[p*3])) > minlen2) {
            drawfieldline = 1;
            break;
          }
        }
      }

      // if drawing style is tubes or spheres we enable the alternate
      // rendering path
      if (drawstyle != 0) {
        // only draw the field line if it met all selection criteria
        if (drawfieldline) {
          cmdColorIndex.putdata(usecolor, cmdList);
          DispCmdCylinder cmdCyl;
          int maxcylidx = (pointcount - 1) * 3;
          int p;
          if (drawstyle == 1) {
            for (p=0; p<maxcylidx; p+=3) {
              cmdCyl.putdata(&points[p], &points[p+3], thickness,
                             tuberes, 0, cmdList);
            }
            maxcylidx++;
          }

          for (p=0; p<maxcylidx; p+=3) {
            cmdSphere.putdata(&points[p], thickness, cmdList);
          }
        }
      } else if (drawfieldline) {
        // only draw the field line if it met all selection criteria
        cmdLineType.putdata(SOLIDLINE, cmdList);
        cmdLineWidth.putdata((int) thickness, cmdList);
        cmdColorIndex.putdata(usecolor, cmdList);
        DispCmdPolyLineArray cmdPolyLineArray;
        cmdPolyLineArray.putdata(&points[0], pointcount, cmdList);
      } 
    }
  }
  wkf_msg_timer_destroy(msgt);

  if (printdonemesg) 
    msgInfo << "field line integration complete." << sendmsg;
}


