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
 *      $RCSfile: vmdfsinfo.c,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $      $Date: 2019/05/22 18:54:04 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Routines for file and filesystem locality determination.
 *
 * LICENSE:
 *   UIUC Open Source License
 *   http://www.ks.uiuc.edu/Research/vmd/plugins/pluginlicense.html
 *
 ***************************************************************************/

#include <stdio.h>
#include <string.h>

#if defined(__sun)
#include <sys/types.h>
#include <sys/statvfs.h>
#endif

#if defined(__linux__) || defined(__ANDROID__)
#include <sys/statfs.h>
#include <sys/vfs.h>
#endif

#if defined(__linux__) || defined(__ANDROID__)
/* Linux-specific constants from coreutils' src/fs.h */
#define  S_MAGIC_ACFS      0x61636673 
#define  S_MAGIC_AFS       0x5346414F
#define  S_MAGIC_AUFS      0x61756673 
#define  S_MAGIC_CEPH      0x00C36400
#define  S_MAGIC_CIFS      0xFF534D42
#define  S_MAGIC_CODA      0x73757245
#define  S_MAGIC_FHGFS     0x19830326
#define  S_MAGIC_FUSEBLK   0x65735546
#define  S_MAGIC_FUSECTL   0x65735543
#define  S_MAGIC_GFS       0x01161970
#define  S_MAGIC_GPFS      0x47504653
#define  S_MAGIC_IBRIX     0x013111A8
#define  S_MAGIC_KAFS      0x6B414653
#define  S_MAGIC_LUSTRE    0x0BD00BD0
#define  S_MAGIC_NCP       0x564C
#define  S_MAGIC_NFS       0x6969
#define  S_MAGIC_NFSD      0x6E667364
#define  S_MAGIC_OCFS2     0x7461636F
#define  S_MAGIC_OVERLAYFS 0x794C7630
#define  S_MAGIC_PANFS     0xAAD7AAEA
#define  S_MAGIC_PIPEFS    0x50495045
#define  S_MAGIC_PRL_FS    0x7C7C6673
#define  S_MAGIC_SMB       0x517B
#define  S_MAGIC_SMB2      0xFE534D42
#define  S_MAGIC_SNFS      0xBEEFDEAD
#define  S_MAGIC_VMHGFS    0xBACBACBC
#define  S_MAGIC_VXFS      0xA501FCF5
#endif


#define FS_ERROR       -1
#define FS_IS_LOCAL     1
#define FS_IS_REMOTE    2
#define FS_IS_UNKNOWN   4


/*
 * Determine locality of the specified file's filesystem.
 * Some systems have statfvs.f_basetype[FSTYPSZ] (AIX, HP-UX, and Solaris).
 * Others have statvfs.f_fstypename[_VFS_NAMELEN] (NetBSD 3.0).
 * Others have statfs.f_fstypename[MFSNAMELEN] (NetBSD 1.5.2).
 * Still others have neither and have to get by with f_type (GNU/Linux).
 * But f_type may only exist in statfs (Cygwin).  
 */
int vmd_fstype_locality(const char *pathname) {
  int fstype = FS_IS_UNKNOWN;

#if defined(__sun)
  struct statvfs svfs;

  /* check filesystem info for the listed file, determine */
  /* if the file is on a local or remote filesystem       */
  memset(&svfs, 0, sizeof(svfs));
  if (statvfs(pathname, &svfs) < 0)
    return FS_ERROR;

#if 0
  printf("VFS basetype: %s\n", svfs.f_basetype);
#endif
  if (!strcmp(svfs.f_basetype, "nfs"))
    fstype = FS_IS_REMOTE;
  else 
    fstype = FS_IS_LOCAL;
#endif

#if defined(__linux__) || defined(__ANDROID__)
  struct statfs sfs;

  // check filesystem info for the listed file, determine
  // if the file is on a local or remote filesystem
  memset(&sfs, 0, sizeof(sfs));
  if (statfs(pathname, &sfs) < 0)
    return FS_ERROR;

  switch (sfs.f_type) { 
#if defined(__linux__) || defined(__ANDROID__)
    /* Linux-specific constants from coreutils' src/fs.h */
    case S_MAGIC_ACFS:      /* 0x61636673 remote */
    case S_MAGIC_AFS:       /* 0x5346414F remote */
    case S_MAGIC_AUFS:      /* 0x61756673 remote */
    case S_MAGIC_CEPH:      /* 0x00C36400 remote */ 
    case S_MAGIC_CIFS:      /* 0xFF534D42 remote */
    case S_MAGIC_CODA:      /* 0x73757245 remote */
    case S_MAGIC_FHGFS:     /* 0x19830326 remote */
    case S_MAGIC_FUSEBLK:   /* 0x65735546 remote */
    case S_MAGIC_FUSECTL:   /* 0x65735543 remote */
    case S_MAGIC_GFS:       /* 0x01161970 remote */
    case S_MAGIC_GPFS:      /* 0x47504653 remote */
    case S_MAGIC_IBRIX:     /* 0x013111A8 remote */
    case S_MAGIC_KAFS:      /* 0x6B414653 remote */
    case S_MAGIC_LUSTRE:    /* 0x0BD00BD0 remote */
    case S_MAGIC_NCP:       /* 0x564C remote */
    case S_MAGIC_NFS:       /* 0x6969 remote */
    case S_MAGIC_NFSD:      /* 0x6E667364 remote */
    case S_MAGIC_OCFS2:     /* 0x7461636F remote */
    case S_MAGIC_OVERLAYFS: /* 0x794C7630 remote */
    case S_MAGIC_PANFS:     /* 0xAAD7AAEA remote */
    case S_MAGIC_PIPEFS:    /* 0x50495045 remote */ 
    case S_MAGIC_PRL_FS:    /* 0x7C7C6673 remote */
    case S_MAGIC_SMB:       /* 0x517B remote */
    case S_MAGIC_SMB2:      /* 0xFE534D42 remote */
    case S_MAGIC_SNFS:      /* 0xBEEFDEAD remote */
    case S_MAGIC_VMHGFS:    /* 0xBACBACBC remote */
    case S_MAGIC_VXFS:      /* 0xA501FCF5 remote */
#endif
#if __GNU__
    case FSTYPE_NFS:
    case FSTYPE_AFS:
    case FSTYPE_FTP:
    case FSTYPE_HTTP:
#endif
#if defined(__linux__) || defined(__ANDROID__) || defined (__GNU__)
      fstype = FS_IS_REMOTE;
#endif

    default:
      fstype = FS_IS_LOCAL;
  }
#endif

  return fstype;
}

