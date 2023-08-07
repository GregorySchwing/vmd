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
*      $RCSfile: OSPRay2DisplayDevice.h,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.5 $         $Date: 2021/12/13 07:54:00 $
*
***************************************************************************/
/**
 *  \file OSPRayDisplayDevice.h
 *  \brief FileRenderer subclass for the OSPRay ray tracing engine.
 * 
 *  This work is briefly outlined in:
 *   "OSPRay - A CPU Ray Tracing Framework for Scientific Visualization"
 *    Ingo Wald, Gregory Johnson, Jefferson Amstutz, Carson Brownlee,
 *    Aaron Knoll, Jim Jeffers, Johannes Guenther, Paul Navratil
 *    IEEE Vis, 2016 (in-press)
 * 
 *  Portions of this code are derived from Tachyon:
 *    "An Efficient Library for Parallel Ray Tracing and Animation"
 *    John E. Stone.  Master's Thesis, University of Missouri-Rolla,
 *    Department of Computer Science, April 1998
 * 
 *    "Rendering of Numerical Flow Simulations Using MPI"
 *    John Stone and Mark Underwood.
 *    Second MPI Developers Conference, pages 138-141, 1996.
 *    http://dx.doi.org/10.1109/MPIDC.1996.534105
 */

#ifndef LIBOSPRAY2DISPLAYDEVICE
#define LIBOSPRAY2DISPLAYDEVICE

#include <stdio.h>
#include "Matrix4.h"
#include "FileRenderer.h"
#include "WKFUtils.h" // timers

/// forward declarations of OSPRay API types 
class OSPRay2Renderer;

/// FileRenderer subclass to exports VMD scenes to OSPRay
class OSPRay2DisplayDevice: public FileRenderer {
private:
  int isinteractive;
  OSPRay2Renderer * ort;
  wkf_timerhandle ort_timer;

  void reset_vars(void); ///< reset internal state variables
  void write_lights(void);
  void write_materials(void);
  void add_material(void);

  int inclipgroup;            ///< whether a clipping group is currently active
#if 0
  int involtex;               ///< volume texturing is enabled
  int voltexID;               ///< current volume texturing ID
  float xplaneeq[4];          ///< volumetric texture plane equations
  float yplaneeq[4];
  float zplaneeq[4];
#endif

#if 1
  // storage and state variables needed to aggregate lone cylinders
  // with common color and material state into a larger buffer for
  // transmission to OSPRay
  int cylinder_matindex;
  Matrix4 *cylinder_xform;
  float cylinder_radius_scalefactor;
  ResizeArray<float> cylinder_vert_buffer;
  ResizeArray<float> cylinder_radii_buffer;
  ResizeArray<float> cylinder_color_buffer;
  // cylinder end caps, made from rings
  ResizeArray<float> cylcap_vert_buffer;
  ResizeArray<float> cylcap_norm_buffer;
  ResizeArray<float> cylcap_radii_buffer;
  ResizeArray<float> cylcap_color_buffer;

  /// reset cylinder buffer to empty state
  void reset_cylinder_buffer() {
    cylinder_matindex = -1;
    cylinder_xform = NULL;
    cylinder_radius_scalefactor=1.0f;
    cylinder_vert_buffer.clear();
    cylinder_radii_buffer.clear();
    cylinder_color_buffer.clear();

    cylcap_vert_buffer.clear();
    cylcap_norm_buffer.clear();
    cylcap_radii_buffer.clear();
    cylcap_color_buffer.clear();
  };
#endif

  // storage and state variables needed to aggregate lone triangles
  // with common color and material state into a larger buffer for
  // transmission to OSPRay
  int triangle_cindex;
  int triangle_matindex;
  Matrix4 *triangle_xform;
  ResizeArray<float> triangle_vert_buffer;
  ResizeArray<float> triangle_norm_buffer;

  /// reset triangle buffer to empty state
  void reset_triangle_buffer() {
    triangle_cindex = -1;   
    triangle_matindex = -1; 
    triangle_xform = NULL;
    triangle_vert_buffer.clear();
    triangle_norm_buffer.clear();
  };

protected:
#if 1
  void send_cylinder_buffer(void);
  void cylinder(float *, float *, float rad, int filled);
  void text(float *pos, float size, float thickness, const char *str);
#endif

#if 0
  void sphere(float *spdata);
#endif
  void sphere_array(int num, int res, float *centers, 
                    float *radii, float *colors);
  void send_triangle_buffer(void);
  void triangle(const float *, const float *, const float *,
                const float *, const float *, const float *);
  void tricolor(const float * xyz1, const float * xyz2, const float * xyz3,
                const float * n1,   const float * n2,   const float * n3,
                const float * c1,   const float * c2,   const float * c3);
  void trimesh_c4u_n3b_v3f(unsigned char *c, signed char *n, float *v, int numfacets);
  void trimesh_c4u_n3f_v3f(unsigned char *c, float *n, float *v, int numfacets);
  void trimesh_c4n3v3(int numverts, float * cnv, int numfacets, int * facets);
  void trimesh_n3b_v3f(signed char *n, float *v, int numfacets);
  void trimesh_n3f_v3f(float *n, float *v, int numfacets);
#if 0
  void trimesh_n3fopt_v3f(float *n, float *v, int numfacets);
#endif
  void tristrip(int numverts, const float * cnv,
                int numstrips, const int *vertsperstrip,
                const int *facets);

#if 0
  // define a volumetric texture map
  virtual void define_volume_texture(int ID, int xs, int ys, int zs,
                                     const float *xplaneeq,
                                     const float *yplaneeq,
                                     const float *zplaneeq,
                                     unsigned char *texmap);

  // enable volumetric texturing, either in "replace" or "modulate" mode
  virtual void volume_texture_on(int texmode);

  // disable volumetric texturing
  virtual void volume_texture_off(void);

  // begin a group of objects to be clipped by the same set of
  // clipping planes
  void start_clipgroup(void);
  void end_clipgroup(void);
#endif

public: 
  static int OSPRay_Global_Init(void);      ///< global init, call ONCE
  static void OSPRay_Global_Shutdown(void); ///< global shutdown, call ONCE
  OSPRay2DisplayDevice(VMDApp *, int interactive);
  virtual ~OSPRay2DisplayDevice(void);
  void write_header(void); 
  void write_trailer(void);
}; 

#endif

