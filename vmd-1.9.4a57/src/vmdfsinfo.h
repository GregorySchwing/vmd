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
 *      $RCSfile: vmdfsinfo.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $      $Date: 2019/05/01 20:17:21 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Thin multi-platform wrapper around file and
 *   filesystem local/remote determination queries for use by
 *   OS kernel-bypass I/O implementations in VMD and the molfile plugins.
 *
 * LICENSE:
 *   UIUC Open Source License 
 *   http://www.ks.uiuc.edu/Research/vmd/plugins/pluginlicense.html
 *
 ***************************************************************************/


#ifndef VMD_FSINFO_H
#define VMD_FSINFO_H

#ifdef __cplusplus
extern "C" {
#endif

#define VMD_FS_ERROR      -1  /**< an error occured during OS FS stat calls  */
#define VMD_FS_IS_LOCAL    1  /**< FS is not remote, therefore assumed local */
#define VMD_FS_IS_REMOTE   2  /**< FS is known to be a remote/network mount  */
#define VMD_FS_IS_UNKNOWN  4  /**< Platform/OS has no determination method   */

/**
 * vmd_fstype_locality: thin multi-platform wrapper around file and
 *            filesystem local/remote determination queries for use by
 *            kernel-bypass I/O implementations in VMD and the molfile plugins.
 */
int vmd_fstype_locality(const char *pathname);

#ifdef __cplusplus
}
#endif

#endif

