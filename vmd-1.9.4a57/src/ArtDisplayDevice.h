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
 *	$RCSfile: ArtDisplayDevice.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.31 $	$Date: 2020/02/26 03:51:30 $
 *
 ***************************************************************************/
/**
 *  \file ArtDisplayDevice.h
 *  \brief FileRenderer subclass to export molecular graphics 
 *  to the ART ray tracer.
 * 
 *  ART is available from gondwana.ecr.mu.oz.au
 *   as part of the vort package.  To see the output I suggest:
 *     art plot.scn 1000 1000
 *     vort2ppm plot.pix > plot.ppm
 *     fromppm plot.ppm plot.rgb
 *     ipaste plot.rgb
 * 
 *  NOTE: As of 2011, the original VORT package home page has vanished,
 *       but it can currently be found here:
 *         http://bund.com.au/~dgh/eric/
 *         http://bund.com.au/~dgh/eric/vort.tar.gz
 *         http://bund.com.au/~dgh/eric/artimages/index.html
 */

#ifndef ARTDISPLAYDEVICE_H
#define ARTDISPLAYDEVICE_H

#include <stdio.h>
#include "FileRenderer.h"

/// FileRenderer subclass to export VMD scenes to ART ray tracer scene format
class ArtDisplayDevice : public FileRenderer {
private:
  int Initialized;    ///< was the output file created?

protected:
  // assorted graphics functions
  void comment(const char *);
  void cone(float *, float *, float, int); 
  void cylinder(float *, float *, float,int filled);
  void line(float *, float *);
  void point(float *);
  void sphere(float *);
  void square(float *, float *, float *, float *, float *);
  void triangle(const float *, const float *, const float *,
                const float *, const float *, const float *);

public: 
  ArtDisplayDevice();
  virtual ~ArtDisplayDevice(void);
  void write_header(void);
  void write_trailer(void);
}; 

#endif

