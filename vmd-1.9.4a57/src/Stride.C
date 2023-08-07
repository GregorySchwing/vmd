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
 *      $RCSfile: Stride.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.51 $       $Date: 2021/11/18 07:33:43 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  Stride interface class.
 *  STRIDE source code:
 *    http://webclu.bio.wzw.tum.de/stride/
 *    http://webclu.bio.wzw.tum.de/stride/stride.tar.gz
 *    https://en.wikipedia.org/wiki/STRIDE
 * 
 *  DSSP/xDSSP/HSSP family; 
 *    https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
 *    https://github.com/cmbi/dssp
 *    https://github.com/cmbi/hssp/releases
 *
 *  Paper on evaluation of secondary structure assignment methods:
 *    https://pubmed.ncbi.nlm.nih.gov/16164759/
 *
 *  Secondary structure assignment for conformationally irregular peptides:
 *    https://pubmed.ncbi.nlm.nih.gov/25424660/
 * 
 ***************************************************************************/

#include <stdio.h>      // for tmpfile
#include "DrawMolecule.h"
#include "Timestep.h"
#include "Residue.h"
#include "Inform.h"
#include "Stride.h"
#include "utilities.h"  // for vmd_delete_file
#include "VMDDir.h"
#include "PeriodicTable.h"

// 
// Despite prolific compile-time and link-time warnings about
// the use of APIs such as tempnam(), we need to use this approach
// because VMD itself doesn't create the STRIDE output, so temp file
// API variants that return an open file handle aren't useful to us,
// and the potential race condition they're meant to protect against
// is inherent when using system() or exec() calls to run an external
// process on input/output files.
//
// Note: the caller must free the returned C string.
//
char * get_temporary_filename(void) {
  char * filename = NULL;
#if defined(ARCH_MACOSXARM64) || defined(ARCH_MACOSXX86) || defined(ARCH_MACOSXX86_64)
  char *tmpstr = tmpnam(NULL);
  if (tmpstr !=NULL)
    filename = strdup(tmpstr);
#else
  filename  = tempnam(NULL, NULL);
#endif
  return filename;
}


// Write a PDB of the given timestep, stripping out all of the non-protein 
// atoms since they are not used for secondary structure assignment.
// The PDB input can be used with STRIDE, DSSP, and likely other SS programs.
// Write the uniq_resid of the residues written to the stride input file
// into the residues argument.
int write_ss_input_pdb(DrawMolecule *mol, const char *inputfilename,
                       ResizeArray<int>& residues) {
  const float *pos = mol->current()->pos;
  char name[6], resname[5], atomicelement[4], chain[2];
 
  residues.clear();

  FILE *inputfile = fopen(inputfilename,"w");
  if (!inputfile) {
    msgErr << "Unable to open input file '" << inputfilename 
           << "' for writing." << sendmsg;
    return 1;
  }

  int prev = -1; // previous uniq_resid
  int atomcount = 0;

  for (int i=0; i < mol->nAtoms; i++) {
    MolAtom *atom = mol->atom(i);
    // skip if this atom isn't part of protein
    if (atom->residueType != RESPROTEIN)
      continue;

    strncpy(name, (mol->atomNames).name(atom->nameindex), 4);
    name[4]='\0';

    strncpy(resname,(mol->resNames).name(atom->resnameindex),4);  
    resname[4] = '\0';

    strncpy(atomicelement, get_pte_label(atom->atomicnumber), 4);  
    atomicelement[1] = ' ';
    atomicelement[2] = ' ';
    atomicelement[3] = '\0';

    chain[0] = (mol->chainNames).name(atom->chainindex)[0];
    chain[1] = '\0';

    if (atom->uniq_resid != prev) {
      prev = atom->uniq_resid;
      residues.append(prev);
    }

    const float *loc=pos + 3L*i;
    const char *emptycols = "                       ";
    if (fprintf(inputfile, "ATOM  %5d %-4s %-4s%c%4ld    %8.3f%8.3f%8.3f%s%3s\n",
                ++atomcount, name, resname, chain[0], long(residues.num()-1),
                loc[0], loc[1], loc[2], emptycols, atomicelement) < 0) {
      msgErr << "Error writing line in Stride/DSSP input file for atom " << i << sendmsg;
      return 1;
    }
  }
  fprintf(inputfile, "TER                                                                             \n");

  fclose(inputfile);
  return 0;
}


//
// Here's the sprintf command stride uses on output in report.c:
//    sprintf(Tmp,"ASG  %3s %c %4s %4d    %c   %11s   %7.2f   %7.2f   %7.1f",
//            p->ResType,SpaceToDash(Chain[Cn]->Id),p->PDB_ResNumb,i+1,
//            p->Prop->Asn,Translate(p->Prop->Asn),p->Prop->Phi,
//            p->Prop->Psi,p->Prop->Solv);
//
//   Here's a typical line:
//ASG  THR -    9    9    B        Bridge   -113.85    157.34      21.9      ~~~~
// 

static int read_stride_record( DrawMolecule *mol, 
                               const char *stridefile,
                               const ResizeArray<int> &residues ) { 
  FILE *in = fopen(stridefile, "rt");
  if (!in) {
    msgErr << "Unable to find Stride output file: "
           << stridefile << sendmsg;
    return 1;
  }

  const int BUF_LEN = 90;
  char buf[BUF_LEN], resname[4], chain[1], uniq_residstr[5], ss;
  int resid=0;

  while (fgets(buf, BUF_LEN, in)) {
    if (strncmp(buf, "ASG", 3))
      continue;

    sscanf(buf,"ASG  %3s %c %4s %4d    %c",
           resname, chain, uniq_residstr, &resid, &ss);

    int index = atoi(uniq_residstr);
    if (index < 0 || index >= residues.num()) {
      msgErr << "invalid resid found in stride output file!" << sendmsg;
      msgErr << "Error found in the following line: " << sendmsg;
      msgErr << buf << sendmsg;
      fclose(in);
      return -1;
    }

    int uniq_resid = residues[index];
    switch (ss) {
      case 'H': // H: helix
	mol->residueList[uniq_resid]->sstruct = SS_HELIX_ALPHA;
	break;

      case 'G': // G: 3-10 helix
	mol->residueList[uniq_resid]->sstruct = SS_HELIX_3_10;
        break; 

      case 'I': // I: PI helix
	mol->residueList[uniq_resid]->sstruct = SS_HELIX_PI;
	break;

      case 'E': // E: Beta sheet
	mol->residueList[uniq_resid]->sstruct = SS_BETA;
	break;

      case 'B': // B: bridge
      case 'b': // b: bridge
	mol->residueList[uniq_resid]->sstruct = SS_BRIDGE;
	break;

      case 'T': // T: turn
        mol->residueList[uniq_resid]->sstruct = SS_TURN;
	break;

      default:
	// msgErr << "Internal error in read_stride_record\n" << sendmsg;
	mol->residueList[uniq_resid]->sstruct = SS_COIL;
    }
  }

  fclose(in);
  return 0;
}  


int ss_from_stride(DrawMolecule *mol) {
  int rc = 0;
  
  char *stridebin = getenv("STRIDE_BIN");
  if (!stridebin) {
    msgErr << "No STRIDE binary found; please set STRIDE_BIN environment variable" << sendmsg;
    msgErr << "to the location of the STRIDE binary." << sendmsg;
    return 1;
  }

  // check to see if the executable exists
  if (!vmd_file_is_executable(stridebin)) {
    msgErr << "STRIDE binary " << stridebin << " cannot be run; check permissions." << sendmsg;
    return 1;
  }

  char *infilename  = get_temporary_filename();
  char *outfilename = get_temporary_filename();
  if (infilename == NULL || outfilename == NULL) {
    msgErr << "Unable to create temporary files for STRIDE." << sendmsg;
    return 1;
  }

  char *stridecall = new char[strlen(stridebin)
                              +strlen(infilename)
                              +strlen(outfilename)
                              +16];

  sprintf(stridecall,"\"%s\" %s -f%s", stridebin, infilename, outfilename);

  ResizeArray<int> residues;

  if (write_ss_input_pdb(mol,infilename,residues)) {
    msgErr << "write_ss_input_pdb(): unable "
           << "to write input file for Stride\n" << sendmsg;
    rc = 1;
  }

  if (!rc) {
    vmd_system(stridecall);

    if (read_stride_record(mol,outfilename,residues)) {
      msgErr << "Stride::read_stride_record: unable "
             << "to read output file from Stride\n" << sendmsg;
      rc = 1;
    }
  }

  delete [] stridecall;
  vmd_delete_file(outfilename);
  vmd_delete_file(infilename);
  free(outfilename);
  free(infilename);

  return rc;
} 



//
// DSSP implementation
//

/* 
  Dssp output format:
return (kDSSPResidueLine % residue.GetNumber() % ca.mResSeq % ca.mICode % chainChar %     code % ss % helix[0] % helix[1] % helix[2] % bend % chirality % bridgelabel[0] %   bridgelabel[1] % bp[0] % bp[1] % sheet % floor(residue.Accessibility() + 0.5) %     NHO[0] % ONH[0] % NHO[1] % ONH[1] % residue.TCO() % residue.Kappa() % alpha %     residue.Phi() % residue.Psi() % ca.mLoc.mX % ca.mLoc.mY % ca.mLoc.mZ %       long_ChainID1 % long_ChainID2).str();

  This is the header line for the residue lines in a DSSP file:

  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA           CHAIN AUTHCHAIN
*/



static int read_dssp_record(DrawMolecule *mol, 
                            const char *dsspfile,
                            const ResizeArray<int> &residues ) { 
  FILE *in = fopen(dsspfile, "rt");
  if (!in) {
    msgErr << "Unable to find DSSP output file: "
           << dsspfile << sendmsg;
    return 1;
  }

  const int BUF_LEN = 1024;
  char buf[BUF_LEN], atom_num[6], chain[2], uniq_residstr[5], ss;
  bool start_flag = false;

  while (fgets(buf, BUF_LEN, in)) {
    // check if the header line is reached
    if (!start_flag){
      if (strncmp(buf, "  #  RESIDUE", 12)) {
        continue;
      } else {
        if (fgets(buf, BUF_LEN, in) == NULL) {
          fclose(in);
          return -1; // failed 
        }
        start_flag = true;
      }
    }

    sscanf(buf,"%5s %4s %1s %c",
           atom_num, uniq_residstr, chain, &ss);

    ss = buf[16];
 
    int index = atoi(uniq_residstr);
    if (!index) 
      continue; // for the gap between chains in DSSP output

    if (index < 0 || index >= residues.num()) {
      msgErr << "invalid resid found in DSSP output file!" << sendmsg;
      msgErr << "Error found in the following line: " << sendmsg;
      msgErr << buf << sendmsg;
      fclose(in);
      return -1;
    }

    int uniq_resid = residues[index];
#if 0
    msgErr << atom_num << " " << uniq_residstr << " " << chain <<  " " << ss << sendmsg;
    msgErr << buf << sendmsg;
#endif

    switch (ss) {
      case 'H': // H
        mol->residueList[uniq_resid]->sstruct = SS_HELIX_ALPHA;
        break;

      case 'B': // B
        mol->residueList[uniq_resid]->sstruct = SS_BRIDGE;
        break; 

      case 'E': // E
        mol->residueList[uniq_resid]->sstruct = SS_BETA;
        break;

      case 'G': // G
        mol->residueList[uniq_resid]->sstruct = SS_HELIX_3_10;
        break;

      case 'I': // I
        mol->residueList[uniq_resid]->sstruct = SS_HELIX_PI;
        break;

      case 'T': // T
        mol->residueList[uniq_resid]->sstruct = SS_TURN;
        break;

      case 'S': // S
                // this is the case where vmd doesn't cover
                // structure code "bent" assigned by DSSP  

      default:
        msgErr << "Internal error in read_dssp_record\n" << sendmsg;
//        printf("DSSP uniq_resid: %d  SS: %c\n", uniq_resid, ss);
//        printf("DSSP buf: '%s'\n\n", buf);
        mol->residueList[uniq_resid]->sstruct = SS_COIL;
    }
  }

  fclose(in);
  return 0;
}  




int ss_from_dssp(DrawMolecule *mol) {
  int rc = 0;
  
  char *dsspbin = getenv("DSSP_BIN");
  if (!dsspbin) {
    msgErr << "No DSSP/XDSSP binary found; please set DSSP_BIN environment variable" << sendmsg;
    msgErr << "to the location of the DSSP binary." << sendmsg;
    return 1;
  }

  // check to see if the executable exists
  if (!vmd_file_is_executable(dsspbin)) {
    msgErr << "DSSP binary " << dsspbin << " cannot be run; check permissions." << sendmsg;
    return 1;
  }

  char *infilename  = get_temporary_filename();
  char *outfilename = get_temporary_filename();
  if (infilename == NULL || outfilename == NULL) {
    msgErr << "Unable to create temporary files for DSSP." << sendmsg;
    return 1;
  }

  char *dsspcall = new char[strlen(dsspbin)
                            +strlen(infilename)
                            +strlen(outfilename)
                            +16];

  sprintf(dsspcall,"\"%s\" -i %s -o %s -v", dsspbin, infilename, outfilename);

  ResizeArray<int> residues;

  if (write_ss_input_pdb(mol,infilename,residues)) {
    msgErr << "write_ss_input_pdb(): unable "
           << "to write input file for DSSP\n" << sendmsg;
    rc = 1;
  }

  if (!rc) {
    vmd_system(dsspcall);

    if (read_dssp_record(mol,outfilename,residues)) {
      msgErr << "Dssp::read_dssp_record: unable "
             << "to read output file from DSSP\n" << sendmsg;
      rc = 1;
    }
  }

  delete [] dsspcall;
  vmd_delete_file(outfilename);
  vmd_delete_file(infilename);
  free(outfilename);
  free(infilename);
  return rc;
} 



