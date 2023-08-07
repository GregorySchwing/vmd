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
 *	$RCSfile: Fragment.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.21 $	$Date: 2020/10/28 15:09:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * A Fragment contains a list of residues which are connected
 * each other, and to no one else.  This is at the residue
 * level, and not the atom level.  The residue numbers are the
 * unique_resid as assigned in BaseMolecule
 *
 ***************************************************************************/
#ifndef FRAGMENT_H
#define FRAGMENT_H

#include "ResizeArray.h"

/// A Fragment contains a list of residues which are connected
/// each other, and to no one else.  This is at the residue
/// level, and not the atom level.  The residue numbers are the
/// unique_resid as assigned in BaseMolecule
class Fragment {
public:
  ResizeArray<int> residues;
  
  Fragment(void) : residues(1) {
  }
  
  int num(void) { return int(residues.num()); }
  int operator [](int i) { return residues[i]; }
  void append(int i) { residues.append(i); }
};

#endif

