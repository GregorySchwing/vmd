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
*      $RCSfile: c_compiler.c,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.2 $         $Date: 2022/02/08 19:12:19 $
*
***************************************************************************/
/**
 *  \file c_compiler.c 
 *  \brief APIs to report compile-time C standard and/or other compiler info
 */

// runtime query of compile-time C compiler language version,
// used by "vmdinfo compilers" command...
const char *c_compiler_std() {
#if   (__STDC_VERSION__ >= 201112L)
      const char *ccversion = "C11"; // C 2011
#elif   (__STDC_VERSION__ >= 199901L)
      const char *ccversion = "C99"; // C 1999 
#elif (__STDC_VERSION__ >= 199409L)
      const char *ccversion = "C90"; // C 1990
#elif defined(__STDC__)
      const char *ccversion = "C89"; // C 1989 - ANSI C
#else
      const char *ccversion = "Unknown C";
#endif

  return ccversion;
}



