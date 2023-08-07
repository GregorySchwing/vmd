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
 *	$RCSfile: Vrml2DisplayDevice.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.31 $	$Date: 2020/07/01 06:09:05 $
 *
 ***************************************************************************/
/**
 *  \file Vrml2DisplayDevice.h
 *  \brief FileRenderer subclass to export molecular graphics to a VRML2 file.
 *
 * VRML2 / VRML97 specification:
 *   http://www.web3d.org/x3d/specifications/vrml/ISO-IEC-14772-VRML97/
 *
 * List of VRML viewers at NIST:
 *   http://cic.nist.gov/vrml/vbdetect.html
 */

#ifndef VRML2DISPLAYDEVICE_H
#define VRML2DISPLAYDEVICE_H

#include <stdio.h>
#include "FileRenderer.h"

/// FileRenderer subclass to export VMD scenes to VRML2/VRML97 scene format
class Vrml2DisplayDevice : public FileRenderer {
private:
  void write_cindexmaterial(int, int); // write colors, materials etc.
  void write_colormaterial(float *, int); // write colors, materials etc.

  void cylinder_noxfrm(float *a, float *b, float rad, int filled);

protected:
  // assorted graphics functions
  void comment(const char *);
  void cone    (float *a, float *b, float rad, int /* resolution */);
  void cylinder(float *a, float *b, float rad, int filled);
  void line(float *xyz1, float *xyz2);
  void point(float *xyz);
  void sphere(float *xyzr);
  void text(float *pos, float size, float thickness, const char *str);
  void triangle(const float *, const float *, const float *,
                const float *, const float *, const float *);
  void tricolor(const float * xyz1, const float * xyz2, const float * xyz3,
                const float * n1,   const float * n2,   const float * n3,
                const float *c1,    const float *c2,    const float *c3);
  virtual void trimesh_c4n3v3(int numverts, float * cnv, 
                              int numfacets, int * facets);
  virtual void trimesh_c4u_n3b_v3f(unsigned char *c, signed char *n,
                                   float *v, int numfacets);
  virtual void tristrip(int numverts, const float * cnv,
                        int numstrips, const int *vertsperstrip,
                        const int *facets);

  void load(const Matrix4& mat);       ///< load transofrmation matrix
  void multmatrix(const Matrix4& mat); ///< concatenate transformation
  void set_color(int color_index);     ///< set the colorID

public:
  Vrml2DisplayDevice(void);
  void write_header(void);
  void write_trailer(void);
}; 

#endif



