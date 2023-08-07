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
 *      $RCSfile: Win32ftp.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.14 $      $Date: 2020/02/26 14:58:21 $
 *
 ***************************************************************************/
/**
 *  \file Win32ftp.h
 *  \brief Windows FTP client used for old webpdb/mol pdbload implementations
 */

#define FTP_FAILURE  -1
#define FTP_SUCCESS   0

vmd_ftpclient(const char * site, const char * remotefile, const char * localfile);
