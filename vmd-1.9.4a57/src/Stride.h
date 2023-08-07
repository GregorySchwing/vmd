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
 *      $RCSfile: Stride.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.14 $       $Date: 2019/04/12 04:41:53 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  Stride interface class.
 ***************************************************************************/

#ifndef STRIDE_H__
#define STRIDE_H__

class DrawMolecule;
#include "ResizeArray.h"

extern int write_ss_input_pdb(DrawMolecule *mol, const char *inputfilename,
                              ResizeArray<int>& residues);

extern int ss_from_stride(DrawMolecule *);
extern int ss_from_dssp(DrawMolecule *);

#endif
