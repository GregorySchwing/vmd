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
 *      $RCSfile: VMDMPI.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $      $Date: 2020/07/23 04:10:10 $
 *
 ***************************************************************************/
#ifndef VMDMPI_INC
#define VMDMPI_INC 1

int vmd_mpi_init(int *argc, char ***argv);
int vmd_mpi_barrier();
int vmd_mpi_fini();
int vmd_mpi_nodeinfo(int *noderank, int *nodecount);
int vmd_mpi_nodescan(int *noderank, int *nodecount,
                     char *nodename, int maxnodenamelen,
                     int gpucount);

#endif
