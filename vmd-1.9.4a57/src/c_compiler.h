/***************************************************************************
 *cr
 *cr            (C) Copyright 2013-2022 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
* RCS INFORMATION:
*
*      $RCSfile: c_compiler.h,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.1 $         $Date: 2022/01/20 18:30:07 $
*
***************************************************************************/
/**
 *  \file c_compiler.c 
 *  \brief APIs to report compile-time C standard and/or other compiler info
 */

#ifndef VMD_C_COMPILER_H
#define VMD_C_COMPILER_H

#ifdef __cplusplus
extern "C" {
#endif

const char *c_compiler_std();

#ifdef __cplusplus
}
#endif

#endif


