/***************************************************************************
 *cr
 *cr		(C) Copyright 1995-2019 The Board of Trustees of the
 *cr			    University of Illinois
 *cr			     All Rights Reserved
 *cr
 ***************************************************************************/
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WavefrontDisplayDevice.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.18 $	$Date: 2020/07/01 06:09:05 $
 *
 ***************************************************************************/
/**
 *  \file WavefrontDisplayDevice.h
 *  \brief FileRenderer subclass to export molecular graphics to
 *         Wavefront OBJ format.
 *
 *   Render molecular scene to a Wavefront Object or "OBJ" file and an
 *   associated .mtl material file.
 *   This file format is one of the most universally supported 3-D model
 *   file formats, particular for animation software such as Maya,
 *   3-D Studio, etc.  The file format is simple, but good enough for getting
 *   a lot of things done.  The Wavefront format only supports per-facet
 *   colors and materials.  This code currently generates a fixed size set of
 *   color/material entries and indexes into this set based on the active VMD
 *   color index.
 */

#ifndef WavefrontDISPLAYDEVICE
#define WavefrontDISPLAYDEVICE

#include <stdio.h>
#include "FileRenderer.h"

/// FileRenderer subclass to export VMD scenes to Wavefront "OBJ" format
class WavefrontDisplayDevice : public FileRenderer {
private:
  FILE *mtlfile;     ///< handle to material file
  char *mtlfilename; ///< name of the material file
  int oldColorIndex;
  int oldMaterialIndex;
  int oldMaterialState;
  void write_material_block(void);        ///< write full material table
  void write_cindexmaterial(int, int);    ///< write colors, materials etc.
  void write_colormaterial(float *, int); ///< write colors, materials etc.

protected:
  void beginrepgeomgroup(const char *);
  void comment(const char *);
  void line(float *xyz1, float *xyz2);
  void point(float *xyz);
  void triangle(const float *, const float *, const float *,
                const float *, const float *, const float *);
  virtual void trimesh_c4n3v3(int numverts, float * cnv, 
                              int numfacets, int * facets);
  virtual void trimesh_c4u_n3b_v3f(unsigned char *c, signed char *n,
                                   float *v, int numfacets);
  virtual void tristrip(int numverts, const float * cnv,
                        int numstrips, const int *vertsperstrip,
                        const int *facets);

public:
  WavefrontDisplayDevice(void);            // constructor
  virtual ~WavefrontDisplayDevice(void);   // destructor
  virtual int open_file(const char *filename);
  virtual void close_file(void);
  void write_header (void);
  void write_trailer(void);
};

#endif

