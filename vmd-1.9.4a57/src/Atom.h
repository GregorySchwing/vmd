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
 *        $RCSfile: Atom.h,v $
 *        $Author: johns $        $Locker:  $                $State: Exp $
 *        $Revision: 1.75 $        $Date: 2020/02/26 03:51:30 $
 *
 ***************************************************************************/
/**
 *  \file Atom.h
 *  \brief Class containing fundamental (required) atomic properties,
 *  topology information, and associated macros.
 */

#ifndef ATOM_H
#define ATOM_H

#include <string.h>
#include <stdlib.h>
#include "utilities.h"

// maximum number of bonds allowed to other atoms
// XXX It may be desirable to finally get around
//     to reimplementing the bond storage using the same scheme that is used
//     for bond types and bond orders, which are each based on dynamically
//     allocated auxilliary arrays associated with string keywords.
//     There should be low impact on structure traversal performance except 
//     for a minor increase in CPU TLB pressure due to traversal of multiple
//     arrays rather than just one, and the code would look almost the same.
#if 0
// Anyone can hack this macro for whatever max bond count they need to 
// support unusual models.  
#define MAXATOMBONDS      256L // example of a huge user-requested bond count
#elif defined(ARCH_BLUEWATERS) || defined(ARCH_CRAY_XC) || defined(ARCH_CRAY_XK)
#define MAXATOMBONDS      8L
#else
#define MAXATOMBONDS      12L
#endif

// Atom type flags
#define ATOMNORMAL      0
#define ATOMPROTEINBACK 1
#define ATOMNUCLEICBACK 2
#define ATOMHYDROGEN    3

// Residue type flags
#define RESNOTHING      0
#define RESPROTEIN      1
#define RESNUCLEIC      2
#define RESWATERS       3

/// class/struct which holds data for one atom
class MolAtom {
public:
  // XXX contents of the Atom structure are ordered specifically so 
  // that the compiler will pack it efficiently.
 
  // items that make this particular atom unique and are absolutely 
  // needed to link it up to the rest of the structure, or are speed-critical
  short nameindex;              ///< atom name string index
  short typeindex;              ///< atom type string index
  int uniq_resid;               ///< unique resid, since there can be dups
  int bondTo[MAXATOMBONDS];     ///< list of atoms to which this atom is bonded
  signed char bonds;            ///< how many bonds this atom has
  signed char atomicnumber;     ///< element atomic number
  signed char altlocindex;      ///< alternate location identifier index
  char insertionstr[2];         ///< for insertion codes (padded to 2 chars)

  // items which could potentially be moved into other data structures 
  // to save memory, but are presently kept here for extra simplicity or speed
  short chainindex;             ///< chain identifier
  short segnameindex;           ///< atom segment name string index
  int resid;                    ///< resid from original file
  short resnameindex;           ///< atom residue name string index

  // ATOMNORMAL, ATOMPROTEINBACK, ATOMNUCLEICBACK, ATOMHYDROGEN
  // XXX this should be converted to an unsigned bit-field to save memory
  signed char atomType;         ///< is this atom part of the backbone?

  /// used to tell me if this atom is part of some larger complete structure
  // RESNOTHING, RESPROTEIN, RESNUCLEIC, RESWATERS
  // XXX this should be converted to an unsigned bit-field to save memory
  signed char residueType;      ///< is this part of a larger component?
                                ///< for instance, is this CG atom in an 
                                ///< amino acid of some sort?

  void init(int n, int theresid, const char *insertion) { 
    uniq_resid = 0; // don't know yet, found in BaseMolecule
    bonds = 0;
    resid = theresid;
    strncpy(insertionstr, insertion, 2);        insertionstr[1] = '\0';
    nameindex = typeindex = resnameindex = segnameindex = altlocindex = (-1);
    atomicnumber = (-1);

    for (int i=0; i<MAXATOMBONDS; i++) {
      bondTo[i] = -1;
    }
    atomType = ATOMNORMAL;
    residueType = RESNOTHING;
  }
  
  /// add a bond into the atom.  Note that each bond will be stored twice, so
  /// make sure when using bonds to define whether the 'official' bond is from
  /// low # atom -> high #, or vice versa
  int add_bond(int a, int type) {
    if(bonds >= MAXATOMBONDS) // fail
      return -1;

    bondTo[(int) bonds] = a;
    if (type == ATOMPROTEINBACK || type == ATOMNUCLEICBACK)
      atomType = type;
    bonds++;
    return 0;
  }

  /// return TRUE if this atom is bonded to the specified atom.  Returns FALSE
  /// otherwise.
  int bonded(int a) {
    for (int i=0; i < bonds; i++) {
      if (bondTo[i] == a) {
        return TRUE;
      }
    }

    return FALSE;
  }

};

#endif

