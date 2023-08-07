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
 *      $RCSfile: androidvmdstart.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.6 $        $Date: 2020/02/26 06:39:27 $
 *
 ***************************************************************************/
/**
 *  \file androidvmdstart.h
 *  \brief Android platform specific application startup code.
 */

extern "C" {
  void log_android(const char *prompt, const char * msg);
}

