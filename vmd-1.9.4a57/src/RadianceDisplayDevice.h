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
 *	$RCSfile: RadianceDisplayDevice.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.34 $	$Date: 2020/02/26 07:21:45 $
 *
 ***************************************************************************/
/**
 *  \file RadianceDisplayDevice.h
 *  \brief FileRenderer subclass to export molecular graphics to Radiance.
 *
 * Writes the scene to Radiance format.  For more information about Radiance,
 * see http://radsite.lbl.gov/radiance/HOME.html
 * Radiance provides programs to convert from its native image format
 * to something other common formats (e.g., ra_ps, ra_tiff, and ra_gif).
 */

#ifndef RADIANCEDISPLAYDEVICE_H
#define RADIANCEDISPLAYDEVICE_H

#include <stdio.h>
#include "FileRenderer.h"

/// FileRenderer subclass to export VMD scenes to a Radiance scene file
class RadianceDisplayDevice : public FileRenderer {
private:
  /// manage the colors
  ResizeArray<float> red;
  ResizeArray<float> green;
  ResizeArray<float> blue;
  ResizeArray<float> trans;
  int cur_color;             ///< active colorID
  void reset_vars(void);     ///< reset internal state variables

protected:
  /// assorted graphics functions
  void comment(const char *);   
  void cone_trunc(float *, float *, float, float, int); 
  void cylinder(float *, float *, float,int);
  void line(float *, float *);
  void point(float *);
  void sphere(float *);
  void square(float *, float *, float *, float *, float *);
  void triangle(const float *, const float *, const float *,
                const float *, const float *, const float *);

  void set_color(int); ///< set current colorID
   
public: 
  RadianceDisplayDevice();
  virtual ~RadianceDisplayDevice(void);
  void write_header(void);
  void write_trailer(void);
}; 

#endif

