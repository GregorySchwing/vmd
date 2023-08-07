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
 *	$RCSfile: SnapshotDisplayDevice.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.27 $	$Date: 2020/02/26 07:21:45 $
 *
 ***************************************************************************/
/**
 * \file SnapshotDisplayDevice.h
 * \brief FileRenderer subclass that saves a DisplayDevice subclass 
 *        screen shot to an image file.
 * 
 * SnapshotDisplayDevice has the screen parameters since it's a DisplayDevice,
 * obtaining the info from "display".  
 */

#ifndef SNAPSHOTDISPLAYDEVICE
#define SNAPSHOTDISPLAYDEVICE

#include "FileRenderer.h"

/// FileRenderer subclass to save VMD images in a supported image file format
class SnapshotDisplayDevice : public FileRenderer {
private:
  DisplayDevice *display;

public:
  /// set up the commands for grabbing images from the screen
  /// pass in display to grab image from
  SnapshotDisplayDevice(DisplayDevice *);
  virtual int open_file(const char *filename);   ///< open output
  virtual void render(const VMDDisplayList*) {}  ///< ignore renders
  virtual void close_file(void); ///< capture and save the image to a file
};  
#endif

