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
 *      $RCSfile: MeasureSr.C,v $
 *      $Author: gregs $        $Locker:  $             $State: Exp $
 *      $Revision: 1.0 $       $Date: 2023/08/07 21:21:00 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Code to compute the distance dependence of the Kirkwood factor, GK(r), 
 *   and the radial distribution of the dipole ordering, s(r), functions for MD trajectories.
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Measure.h"
#include "AtomSel.h"
#include "utilities.h"
#include "ResizeArray.h"
#include "MoleculeList.h"
#include "Inform.h"
#include "Timestep.h"
#include "VMDApp.h"
#include "WKFThreads.h"
#include "WKFUtils.h"
#include "CUDAAccel.h"
#include "CUDAKernels.h"

#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))

/*! the volume of a spherical cap is defined as:
 * <pre> pi / 9 * h^2 * (3 * r  - h) </pre>
 * with h = height of cap = radius - boxby2.
 * \brief the volume of a sperical cap. */
static inline double spherical_cap(const double &radius, const double &boxby2) {
  return (VMD_PI / 3.0 * (radius - boxby2) * (radius - boxby2)
          * ( 2.0 * radius + boxby2));
}

void dipole_cpu(int natoms1,  // array of the number of atoms in
                              // selection 1 in each frame. 
             int natoms2,     // array of the number of atoms in
                              // selection 2 in each frame. 
             int natoms3,     // array of the number of atoms in
                              // selection 3 in each frame. 
             float* xyz3,     // coordinates of selection 3.
                              // [natoms2][3]
             int natoms4,     // array of the number of atoms in
                              // selection 4 in each frame. 
             float* xyz4,     // coordinates of selection 4.
                              // [natoms2][3])
             float* q3,       // charges of atoms in selection3.
                              // [natoms])
             float* m3,       // masses of atoms in selection3.
                              // [natoms])
             float* rvec3,    // r vector of dipoles in selection3.
                              // [natoms3/natoms1][3])
             float* qrvec3,   // qr vector of dipoles in selection4.
                              // [natoms3/natoms1][3])
             float* mrvec3,   // mr vector of dipoles in selection4.
                              // [natoms3/natoms1][3])
             float* totalQ3,  // total charge of each molecule in selection3.
                              // [natoms3/natoms1])
             float* totalM3,   // total mass of each molecule in selection3.
                              // [natoms3/natoms1])
             float* q4,       // charges of atoms in selection4.
                              // [natoms])
             float* m4,       // masses of atoms in selection4.
                              // [natoms])
             float* rvec4,    // r vector of dipoles in selection4.
                              // [natoms4/natoms2][3])
             float* qrvec4,   // qr vector of dipoles in selection4.
                              // [natoms4/natoms2][3])
             float* mrvec4,   // mr vector of dipoles in selection4.
                              // [natoms4/natoms2][3])
             float* totalQ4,  // total charge of each molecule in selection4.
                              // [natoms4/natoms2])
             float* totalM4,   // total mass of each molecule in selection4.
                              // [natoms3/natoms1])
             float* dipoles3, // dipoles of selection 3.
                              // [natoms3/natoms1][3])
             float* dipoles4, // dipoles of selection 4.
                              // [natoms4/natoms2][3])
             int usecenter)   // Use geometric mean
{
  int iatom;
  int moleculeSize3 = natoms3/natoms1;
  for (iatom=0; iatom<natoms3; iatom++) {
    long addr = 3L * iatom;
    long dipoleaddr = 3L * iatom/moleculeSize3;

    rvec3[dipoleaddr    ] += xyz3[addr    ];
    rvec3[dipoleaddr + 1] += xyz3[addr + 1];
    rvec3[dipoleaddr + 2] += xyz3[addr + 2];

    qrvec3[dipoleaddr    ] += xyz3[addr    ] * q3[iatom];
    qrvec3[dipoleaddr + 1] += xyz3[addr + 1] * q3[iatom];
    qrvec3[dipoleaddr + 2] += xyz3[addr + 2] * q3[iatom];

    mrvec3[dipoleaddr    ] += xyz3[addr    ] * m3[iatom];
    mrvec3[dipoleaddr + 1] += xyz3[addr + 1] * m3[iatom];
    mrvec3[dipoleaddr + 2] += xyz3[addr + 2] * m3[iatom];

    totalQ3[iatom/moleculeSize3]+=q3[iatom];
    totalM3[iatom/moleculeSize3]+=m3[iatom];

  }
  msgInfo << "Finished sel3 dipoles..." << sendmsg;

  int moleculeSize4 = natoms4/natoms2;
  for (iatom=0; iatom<natoms4; iatom++) {
    long addr = 3L * iatom;
    long dipoleaddr = 3L * iatom/moleculeSize4;

    rvec4[dipoleaddr    ] += xyz4[addr    ];
    rvec4[dipoleaddr + 1] += xyz4[addr + 1];
    rvec4[dipoleaddr + 2] += xyz4[addr + 2];

    qrvec4[dipoleaddr    ] += xyz4[addr    ] * q4[iatom];
    qrvec4[dipoleaddr + 1] += xyz4[addr + 1] * q4[iatom];
    qrvec4[dipoleaddr + 2] += xyz4[addr + 2] * q4[iatom];

    mrvec4[dipoleaddr    ] += xyz4[addr    ] * m4[iatom];
    mrvec4[dipoleaddr + 1] += xyz4[addr + 1] * m4[iatom];
    mrvec4[dipoleaddr + 2] += xyz4[addr + 2] * m4[iatom];

    totalQ4[iatom/moleculeSize4]+=q4[iatom];
    totalM4[iatom/moleculeSize4]+=m4[iatom];

  }
  msgInfo << "Finished sel4 dipoles..." << sendmsg;

  switch (usecenter) {
      case 1:
      {
          for (iatom=0; iatom<natoms3/moleculeSize3; iatom++) {
            double rscale = totalQ3[iatom] / moleculeSize3; 
            long addr = 3L * iatom;
            dipoles3[addr    ] = (float) (qrvec3[addr    ] - (rvec3[addr    ] * rscale)); 
            dipoles3[addr + 1] = (float) (qrvec3[addr + 1] - (rvec3[addr + 1] * rscale)); 
            dipoles3[addr + 2] = (float) (qrvec3[addr + 2] - (rvec3[addr + 2] * rscale)); 
          }
          break;
      }

      case -1: 
      {
          for (iatom=0; iatom<natoms3/moleculeSize3; iatom++) {
            double mscale = totalQ3[iatom] / totalM3[iatom] / moleculeSize3; 
            long addr = 3L * iatom;
            dipoles3[addr    ] = (float) (qrvec3[addr    ] - (mrvec3[addr    ] * mscale)); 
            dipoles3[addr + 1] = (float) (qrvec3[addr + 1] - (mrvec3[addr + 1] * mscale)); 
            dipoles3[addr + 2] = (float) (qrvec3[addr + 2] - (mrvec3[addr + 2] * mscale)); 
          }
          break;
      }
      
      case 0: // fallthrough
      default: 
      {
          for (iatom=0; iatom<natoms3/moleculeSize3; iatom++) {
            long addr = 3L * iatom;
            dipoles3[addr    ] = qrvec3[addr    ]; 
            dipoles3[addr + 1] = qrvec3[addr + 1]; 
            dipoles3[addr + 2] = qrvec3[addr + 2]; 
          }
          break;
      }
  }

  switch (usecenter) {
      case 1:
      {
          for (iatom=0; iatom<natoms4/moleculeSize4; iatom++) {
            double rscale = totalQ4[iatom] / moleculeSize4; 
            long addr = 3L * iatom;
            dipoles4[addr    ] = (float) (qrvec4[addr    ] - (rvec4[addr    ] * rscale)); 
            dipoles4[addr + 1] = (float) (qrvec4[addr + 1] - (rvec4[addr + 1] * rscale)); 
            dipoles4[addr + 2] = (float) (qrvec4[addr + 2] - (rvec4[addr + 2] * rscale)); 
          }
          break;
      }

      case -1: 
      {
          for (iatom=0; iatom<natoms4/moleculeSize4; iatom++) {
            double mscale = totalQ4[iatom] / totalM4[iatom] / moleculeSize4; 
            long addr = 3L * iatom;
            dipoles4[addr    ] = (float) (qrvec4[addr    ] - (mrvec4[addr    ] * mscale)); 
            dipoles4[addr + 1] = (float) (qrvec4[addr + 1] - (mrvec4[addr + 1] * mscale)); 
            dipoles4[addr + 2] = (float) (qrvec4[addr + 2] - (mrvec4[addr + 2] * mscale)); 
          }
          break;
      }
      
      case 0: // fallthrough
      default: 
      {
          for (iatom=0; iatom<natoms4/moleculeSize4; iatom++) {
            long addr = 3L * iatom;
            dipoles4[addr    ] = qrvec4[addr    ]; 
            dipoles4[addr + 1] = qrvec4[addr + 1]; 
            dipoles4[addr + 2] = qrvec4[addr + 2]; 
          }
          break;
      }
  }
}


/*! \brief Calculate 1.0/sqrt(x) in double precision, but single range
 *
 * \param  x  Positive value to calculate inverse square root for, must be
 *            in the input domain valid for single precision.
 *
 * For now this is implemented with std::sqrt(x). However, we might
 * decide to use instrinsics or compiler-specific functions in the future, and
 * then we want to have the freedom to do the first step in single precision.
 *
 * \return 1.0/sqrt(x)
 */
static inline double invsqrt(double x)
{
    return 1.0 / std::sqrt(x);
}

/* WARNING:
 * Do _not_ use these routines to calculate the angle between two vectors
 * as acos(cos_angle(u,v)). While it might seem obvious, the acos function
 * is very flat close to -1 and 1, which will lead to accuracy-loss.
 * Instead, use the new gmx_angle() function directly.
 */
static inline float cos_angle(const float * a, const float * b)
{
    /*
     *                  ax*bx + ay*by + az*bz
     * cos-vec (a,b) =  ---------------------
     *                      ||a|| * ||b||
     */
    float   cosval;
    int    m;
    double aa, bb, ip, ipa, ipb, ipab; /* For accuracy these must be double! */

    ip = ipa = ipb = 0.0;
    for (m = 0; (m < 3); m++) /* 18 */
    {
        aa = a[m];
        bb = b[m];
        ip += aa * bb;
        ipa += aa * aa;
        ipb += bb * bb;
    }
    ipab = ipa * ipb;
    if (ipab > 0)
    {
        cosval = static_cast<float>(ip * invsqrt(ipab)); /*  7 */
    }
    else
    {
        cosval = 1;
    }
    /* 25 TOTAL */
    if (cosval > 1.0)
    {
        return 1.0;
    }
    if (cosval < -1.0)
    {
        return -1.0;
    }

    return cosval;
}

void sr_cpu(int natoms1,     // array of the number of atoms in
                              // selection 1 in each frame. 
             float* xyz,      // coordinates of first selection.
                              // [natoms1][3]
             int natoms2,     // array of the number of atoms in
                              // selection 2 in each frame. 
             float* xyz2,     // coordinates of selection 2.
                              // [natoms2][3]
             float* dipoles3,     // dipoles of selection 3.
                              // [natoms2][3]
             float* dipoles4,     // dipoles of selection 4.
                              // [natoms2][3]
             float* cell,     // the cell x y and z dimensions [3]
             float* hist,     // the histograms, 1 per block
                              // [ncudablocks][maxbin]
             float* hist_d,   // total projection per block
                              // [ncudablocks][maxbin]
             int maxbin,      // the number of bins in the histogram
             float rmin,      // the minimum value of the first bin
             float delr)      // the width of each bin
{
  int iatom, jatom, ibin;
  float rij, rxij, rxij2, x1, y1, z1, x2, y2, z2;
  float cellx, celly, cellz;
  float cos_theta;
  for (ibin=0; ibin<maxbin; ibin++) {
    hist[ibin]=0;
    hist_d[ibin]=0;
  }

  cellx = cell[0];
  celly = cell[1];
  cellz = cell[2];

  for (iatom=0; iatom<natoms1; iatom++) {
    long addr = 3L * iatom;
    xyz[addr    ] = fmodf(xyz[addr    ], cellx);
    xyz[addr + 1] = fmodf(xyz[addr + 1], celly);
    xyz[addr + 2] = fmodf(xyz[addr + 2], cellz);
  }

  for (iatom=0; iatom<natoms2; iatom++) {
    long addr = 3L * iatom;
    xyz2[addr    ] = fmodf(xyz2[addr    ], cellx);
    xyz2[addr + 1] = fmodf(xyz2[addr + 1], celly);
    xyz2[addr + 2] = fmodf(xyz2[addr + 2], cellz);
  }

  for (iatom=0; iatom<natoms1; iatom++) {
    x1 = xyz[3L * iatom    ];
    y1 = xyz[3L * iatom + 1];
    z1 = xyz[3L * iatom + 2];
    for (jatom=0;jatom<natoms2;jatom++) {
      x2 = xyz2[3L * jatom    ];
      y2 = xyz2[3L * jatom + 1];
      z2 = xyz2[3L * jatom + 2];

      rxij = fabsf(x1 - x2);
      rxij2 = cellx - rxij;
      rxij = MIN(rxij, rxij2);
      rij = rxij * rxij;

      rxij = fabsf(y1 - y2);
      rxij2 = celly - rxij;
      rxij = MIN(rxij, rxij2);
      rij += rxij * rxij;

      rxij = fabsf(z1 - z2);
      rxij2 = cellz - rxij;
      rxij = MIN(rxij, rxij2);
      rij += rxij * rxij;

      rij = sqrtf(rij);

      ibin = (int)floorf((rij-rmin)/delr);
      cos_theta = cos_angle(&dipoles3[3L*iatom],&dipoles4[3L*jatom]);
      if (ibin<maxbin && ibin>=0) {
        ++hist[ibin];
        hist_d[ibin]+=cos_theta;
      }
    }
  }

}




int measure_sr(VMDApp *app,
                AtomSel *sel1, AtomSel *sel2, 
                AtomSel *sel3, AtomSel *sel4,
                MoleculeList *mlist,
                const int count_h, double *gofr, 
                double *numint, double *histog,
                double *Gkr, double *avgcos, 
                const float delta, int first, int last, int step, 
                int *framecntr, int usepbc, int selupdate) {
  int i, j, frame;

  float a, b, c, alpha, beta, gamma;
  float pbccell[3];
  int isortho=0;     // orthogonal unit cell not assumed by default.
  int duplicates=0;  // counter for duplicate atoms in both selections.
  float rmin = 0.0f; // min distance to histogram

  // initialize a/b/c/alpha/beta/gamma to arbitrary defaults to please the compiler.
  a=b=c=9999999.0;
  alpha=beta=gamma=90.0;

  // reset counter for total, skipped, and _orth processed frames.
  framecntr[0]=framecntr[1]=framecntr[2]=0;

  // First round of sanity checks.
  // neither list can be undefined
  if (!sel1 || !sel2 ) {
    return MEASURE_ERR_NOSEL;
  }

  // make sure that both selections are from the same molecule
  // so that we know that PBC unit cell info is the same for both
  if (sel2->molid() != sel1->molid()) {
    return MEASURE_ERR_MISMATCHEDMOLS;
  }

  Molecule *mymol = mlist->mol_from_id(sel1->molid());
  int maxframe = mymol->numframes() - 1;
  int nframes = 0;

  if (last == -1)
    last = maxframe;

  if ((last < first) || (last < 0) || (step <=0) || (first < -1)
      || (maxframe < 0) || (last > maxframe)) {
      msgErr << "measure sr: bad frame range given."
             << " max. allowed frame#: " << maxframe << sendmsg;
    return MEASURE_ERR_BADFRAMERANGE;
  }

  // test for non-orthogonal PBC cells, zero volume cells, etc.
  if (usepbc) {
    for (isortho=1, nframes=0, frame=first; frame <=last; ++nframes, frame += step) {
      const Timestep *ts;

      if (first == -1) {
        // use current frame only. don't loop.
        ts = sel1->timestep(mlist);
        frame=last;
      } else {
        ts = mymol->get_frame(frame);
      }
      // get periodic cell information for current frame
      a = ts->a_length;
      b = ts->b_length;
      c = ts->c_length;
      alpha = ts->alpha;
      beta = ts->beta;
      gamma = ts->gamma;

      // check validity of PBC cell side lengths
      if (fabsf(a*b*c) < 0.0001) {
        msgErr << "measure sr: unit cell volume is zero." << sendmsg;
        return MEASURE_ERR_GENERAL;
      }

      // check PBC unit cell shape to select proper low level algorithm.
      if ((alpha != 90.0) || (beta != 90.0) || (gamma != 90.0)) {
        isortho=0;
      }
    }
  } else {
    // initialize a/b/c/alpha/beta/gamma to arbitrary defaults
    isortho=1;
    a=b=c=9999999.0;
    alpha=beta=gamma=90.0;
  }

  // until we can handle non-orthogonal periodic cells, this is fatal
  if (!isortho) {
    msgErr << "measure sr: only orthorhombic cells are supported (for now)." << sendmsg;
    return MEASURE_ERR_GENERAL;
  }

  // clear the result arrays
  for (i=0; i<count_h; ++i) {
    gofr[i] = Gkr[i] = avgcos[i] = numint[i] = histog[i] = 0.0;
  }
  const float *q = mymol->charge();
  const float *m = mymol->mass();
  // pre-allocate coordinate buffers of the max size we'll
  // ever need, so we don't have to reallocate if/when atom
  // selections are updated on-the-fly
  float *sel1coords = new float[3L*sel1->num_atoms];
  float *sel2coords = new float[3L*sel2->num_atoms];
  float *sel3coords = new float[3L*sel3->num_atoms];
  float *sel4coords = new float[3L*sel4->num_atoms];
  float *sel3q = new float[3L*sel3->num_atoms];
  float *sel4q = new float[3L*sel4->num_atoms];
  float *sel3m = new float[3L*sel3->num_atoms];
  float *sel4m = new float[3L*sel4->num_atoms];
  float *sel3rvec = new float[3L*sel3->num_atoms];
  float *sel3qrvec = new float[3L*sel4->num_atoms];
  float *sel3mrvec = new float[3L*sel3->num_atoms];
  float *sel3totalq = new float[sel3->num_atoms];
  float *sel3totalm = new float[sel4->num_atoms];
  float *sel4rvec = new float[3L*sel4->num_atoms];
  float *sel4qrvec = new float[3L*sel3->num_atoms];
  float *sel4mrvec = new float[3L*sel4->num_atoms];
  float *sel4totalq = new float[sel3->num_atoms];
  float *sel4totalm = new float[sel4->num_atoms];
  float *sel3dipoles = new float[sel3->num_atoms];
  float *sel4dipoles = new float[sel4->num_atoms];
  float *lhist = new float[count_h];
  float *lhist_dipoles = new float[count_h];

  for (nframes=0,frame=first; frame <=last; ++nframes, frame += step) {
    const Timestep *ts1, *ts2, *ts3, *ts4;

    if (frame  == -1) {
      // use current frame only. don't loop.
      ts1 = sel1->timestep(mlist);
      ts2 = sel2->timestep(mlist);
      ts3 = sel3->timestep(mlist);
      ts4 = sel4->timestep(mlist);
      frame=last;
    } else {
      sel1->which_frame = frame;
      sel2->which_frame = frame;
      sel3->which_frame = frame;
      sel4->which_frame = frame;
      ts1 = ts2 = ts3 = ts4 = mymol->get_frame(frame); // requires sels from same mol
    }

    if (usepbc) {
      // get periodic cell information for current frame
      a     = ts1->a_length;
      b     = ts1->b_length;
      c     = ts1->c_length;
      alpha = ts1->alpha;
      beta  = ts1->beta;
      gamma = ts1->gamma;
    }

    // compute half periodic cell size
    float boxby2[3];
    boxby2[0] = 0.5f * a;
    boxby2[1] = 0.5f * b;
    boxby2[2] = 0.5f * c;

    // update the selections if the user desires it
    if (selupdate) {
      if (sel1->change(NULL, mymol) != AtomSel::PARSE_SUCCESS)
        msgErr << "measure sr: failed to evaluate atom selection update";
      if (sel2->change(NULL, mymol) != AtomSel::PARSE_SUCCESS)
        msgErr << "measure sr: failed to evaluate atom selection update";
      if (sel3->change(NULL, mymol) != AtomSel::PARSE_SUCCESS)
        msgErr << "measure sr: failed to evaluate atom selection update";
      if (sel4->change(NULL, mymol) != AtomSel::PARSE_SUCCESS)
        msgErr << "measure sr: failed to evaluate atom selection update";
    }

    // check for duplicate atoms in the two lists, as these will have
    // to be subtracted back out of the first histogram slot
    if (sel2->molid() == sel1->molid()) {
      int i;
      for (i=0, duplicates=0; i<sel1->num_atoms; ++i) {
        if (sel1->on[i] && sel2->on[i])
          ++duplicates;
      }
    }

    // copy selected atoms to the two coordinate lists
    // requires that selections come from the same molecule
    const float *framepos = ts1->pos;
    for (i=0, j=0; i<sel1->num_atoms; ++i) {
      if (sel1->on[i]) {
        long a = i*3L;
        sel1coords[j    ] = framepos[a    ];
        sel1coords[j + 1] = framepos[a + 1];
        sel1coords[j + 2] = framepos[a + 2];
        j+=3;
      }
    }
    framepos = ts2->pos;
    for (i=0, j=0; i<sel2->num_atoms; ++i) {
      if (sel2->on[i]) {
        long a = i*3L;
        sel2coords[j    ] = framepos[a    ];
        sel2coords[j + 1] = framepos[a + 1];
        sel2coords[j + 2] = framepos[a + 2];
        j+=3;
      }
    }
    framepos = ts3->pos;
    for (i=0, j=0; i<sel3->num_atoms; ++i) {
      if (sel3->on[i]) {
        long a = i*3L;
        sel3coords[j    ] = framepos[a    ];
        sel3coords[j + 1] = framepos[a + 1];
        sel3coords[j + 2] = framepos[a + 2];
        j+=3;
        sel3q[j]=q[i];
        sel3m[j]=m[i];
      }
    }
    framepos = ts4->pos;
    for (i=0, j=0; i<sel4->num_atoms; ++i) {
      if (sel4->on[i]) {
        long a = i*3L;
        sel4coords[j    ] = framepos[a    ];
        sel4coords[j + 1] = framepos[a + 1];
        sel4coords[j + 2] = framepos[a + 2];
        j+=3;
        sel4q[j]=q[i];
        sel4m[j]=m[i];
      }
    }
    // copy unit cell information
    pbccell[0]=a;
    pbccell[1]=b;
    pbccell[2]=c;

    // clear the histogram for this frame
    memset(lhist, 0, count_h * sizeof(float));
    memset(lhist_dipoles, 0, count_h * sizeof(float));

    // Clear intermediate arrays
    memset(sel3rvec, 0, 3L*sel3->num_atoms * sizeof(float));
    memset(sel3qrvec, 0, 3L*sel3->num_atoms * sizeof(float));
    memset(sel3mrvec, 0, 3L*sel3->num_atoms * sizeof(float));
    memset(sel3totalq, 0, sel3->num_atoms * sizeof(float));
    memset(sel3totalm, 0, sel3->num_atoms * sizeof(float));

    memset(sel4rvec, 0, 3L*sel4->num_atoms * sizeof(float));
    memset(sel4qrvec, 0, 3L*sel4->num_atoms * sizeof(float));
    memset(sel4mrvec, 0, 3L*sel4->num_atoms * sizeof(float));
    memset(sel4totalq, 0, sel4->num_atoms * sizeof(float));
    memset(sel4totalm, 0, sel4->num_atoms * sizeof(float));

    memset(sel3dipoles, 0, sel3->num_atoms * sizeof(float));
    memset(sel4dipoles, 0, sel4->num_atoms * sizeof(float));

    if (isortho && sel1->selected && sel2->selected) {
      // do the sr calculation for orthogonal boxes.
      // XXX. non-orthogonal box not supported yet. detected and handled above.
      int rc=-1;
#if defined(VMDCUDA)
      if (!getenv("VMDNOCUDA") && (app->cuda != NULL)) {
        msgInfo << "Running multi-GPU sr calc..." << sendmsg;
        rc=rdf_gpu(app->cuda->get_cuda_devpool(),
                   usepbc,
                   sel1->selected, sel1coords,
                   sel2->selected, sel2coords, 
                   pbccell,
                   lhist,
                   count_h,
                   rmin,
                   delta);
      } 
#endif
      if (rc != 0) {
        msgInfo << "Running single-core CPU sr calc..." << sendmsg;
        dipole_cpu(sel1->selected, sel2->selected,
                    sel3->selected, sel3coords,
                    sel4->selected, sel4coords,
                    sel3q,sel3m,
                    sel3rvec,sel3qrvec,sel3mrvec,
                    sel3totalq,sel3totalm,
                    sel4q,sel4m,
                    sel4rvec,sel4qrvec,sel4mrvec,
                    sel4totalq,sel4totalm,
                    sel3dipoles,sel4dipoles,1);
        msgInfo << "Return from dipoles_cpu..." << sendmsg;
        sr_cpu(sel1->selected, sel1coords,
                sel2->selected, sel2coords,
                sel3dipoles,sel4dipoles,
                pbccell,
                lhist,
                lhist_dipoles,
                count_h,
                rmin,
                delta);
        lhist[0] -= duplicates;
        msgInfo << "Return from sr_cpu..." << sendmsg;
      }

      ++framecntr[2]; // frame processed with sr algorithm
    } else {
      ++framecntr[1]; // frame skipped
    }
    ++framecntr[0];   // total frames.

#if 0
    // XXX elimination of duplicates is now handled within the 
    //     GPU kernels themselves, so we do not need to subtract them
    //     off during the histogram normalization calculations.
    // Correct the first histogram slot for the number of atoms that are
    // present in both lists. they'll end up in the first histogram bin.
    // we subtract only from the first thread histogram which is always defined.
    lhist[0] -= duplicates;
#endif

    // in case of going 'into the edges', we should cut
    // off the part that is not properly normalized to
    // not confuse people that don't know about this.
    int h_max=count_h;
    float smallside=a;
    if (isortho && usepbc) {
      if(b < smallside) {
        smallside=b;
      }
      if(c < smallside) {
        smallside=c;
      }
      h_max=(int) (sqrtf(0.5f)*smallside/delta) +1;
      if (h_max > count_h) {
        h_max=count_h;
      }
    }

    // compute normalization function.
    double all=0.0;
    double pair_dens = 0.0;
    
    if (sel1->selected && sel2->selected) {
      if (usepbc) {
        pair_dens = a * b * c / ((double)sel1->selected * (double)sel2->selected - (double)duplicates);
      } else { // assume a particle volume of 30 \AA^3 (~ 1 water).
        pair_dens = 30.0 * (double)sel1->selected /
          ((double)sel1->selected * (double)sel2->selected - (double)duplicates);
      }
    }
    msgInfo << "Accumulating results..." << sendmsg;
    // XXX for orthogonal boxes, we can reduce this to rmax < sqrt(0.5)*smallest side
    double GkrSum = 0.0;
    for (i=0; i<h_max; ++i) {
      // radius of inner and outer sphere that form the spherical slice
      double r_in  = delta * (double)i;
      double r_out = delta * (double)(i+1);
      double slice_vol = 4.0 / 3.0 * VMD_PI
        * ((r_out * r_out * r_out) - (r_in * r_in * r_in));

      if (isortho && usepbc) {
        // add correction for 0.5*box < r <= sqrt(0.5)*box
        if (r_out > boxby2[0]) {
          slice_vol -= 2.0 * spherical_cap(r_out, boxby2[0]);
        }
        if (r_out > boxby2[1]) {
          slice_vol -= 2.0 * spherical_cap(r_out, boxby2[1]);
        }
        if (r_out > boxby2[2]) {
          slice_vol -= 2.0 * spherical_cap(r_out, boxby2[2]);
        }
        if (r_in > boxby2[0]) {
          slice_vol += 2.0 * spherical_cap(r_in, boxby2[0]);
        }
        if (r_in > boxby2[1]) {
          slice_vol += 2.0 * spherical_cap(r_in, boxby2[1]);
        }
        if (r_in > boxby2[2]) {
          slice_vol += 2.0 * spherical_cap(r_in, boxby2[2]);
        }
      }

      double normf = pair_dens / slice_vol;
      double histv = (double) lhist[i];

      gofr[i]   += normf * histv;
      GkrSum    += lhist_dipoles[i];
      Gkr[i]    += GkrSum;
      avgcos[i] += normf * lhist_dipoles[i];
      all       += histv;
      if (sel1->selected) {
        numint[i] += all / (double)(sel1->selected);
      }
      histog[i] += histv;
    }
    msgInfo << "Finished accumulating results..." << sendmsg;

  }
  msgInfo << "Deleting arrays..." << sendmsg;
  msgInfo << "Deleting sel1coords..." << sendmsg;

  delete [] sel1coords;
    msgInfo << "Deleting sel2coords..." << sendmsg;

  delete [] sel2coords;
      msgInfo << "Deleting sel3coords..." << sendmsg;

  delete [] sel3coords;
        msgInfo << "Deleting sel4coords..." << sendmsg;

  delete [] sel4coords;
        msgInfo << "Deleting sel3q..." << sendmsg;

  delete [] sel3q;
        msgInfo << "Deleting sel4q..." << sendmsg;

  delete [] sel4q;
        msgInfo << "Deleting sel3m..." << sendmsg;

  delete [] sel3m;
        msgInfo << "Deleting sel4m..." << sendmsg;

  delete [] sel4m;
        msgInfo << "Deleting sel3rvec..." << sendmsg;

  delete [] sel3rvec;
        msgInfo << "Deleting sel3qrvec..." << sendmsg;

  delete [] sel3qrvec;
        msgInfo << "Deleting sel3mrvec..." << sendmsg;

  delete [] sel3mrvec;
        msgInfo << "Deleting sel3totalq..." << sendmsg;

  delete [] sel3totalq;
        msgInfo << "Deleting sel3totalm..." << sendmsg;

  delete [] sel3totalm;
        msgInfo << "Deleting sel4rvec..." << sendmsg;

  delete [] sel4rvec;
        msgInfo << "Deleting sel4qrvec..." << sendmsg;

  delete [] sel4qrvec;
        msgInfo << "Deleting sel4mrvec..." << sendmsg;

  delete [] sel4mrvec;
        msgInfo << "Deleting sel4totalq..." << sendmsg;

  delete [] sel4totalq;
        msgInfo << "Deleting sel4totalm..." << sendmsg;

  delete [] sel4totalm;
        msgInfo << "Deleting sel3dipoles..." << sendmsg;

  delete [] sel3dipoles;
        msgInfo << "Deleting sel4dipoles..." << sendmsg;

  delete [] sel4dipoles;
        msgInfo << "Deleting lhist..." << sendmsg;

  delete [] lhist;
        msgInfo << "Deleting lhist_dipoles..." << sendmsg;

  delete [] lhist_dipoles;
  msgInfo << "Finished deleting arrays..." << sendmsg;

  int ngrp = sel1->num_atoms;
  double norm = 1.0 / (double) nframes;
  double normMol = 1.0 / ((double) nframes * (double) ngrp);
  msgInfo << "Normalizing results..." << sendmsg;

  for (i=0; i<count_h; ++i) {
    gofr[i]   *= norm;
    numint[i] *= norm;
    histog[i] *= norm;
    Gkr[i] *= normMol;
    avgcos[i] *= norm;
  }
  msgInfo << "Returning results..." << sendmsg;

  return MEASURE_NOERR;
}

