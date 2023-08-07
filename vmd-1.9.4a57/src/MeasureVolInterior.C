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
 *      $RCSfile: MeasureVolInterior.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.17 $      $Date: 2020/12/13 07:41:55 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Method for measuring the interior volume in a vesicle or capsid based    
 *   on an externally-provided simulated density map (e.g., from QuickSurf)
 *   and a threshold isovalue for inside/outside tests based on density.
 *   The approach computes and counts interior/exterior voxels based on 
 *   a parallel ray casting approach on an assumed orthorhombic grid.
 *   Juan R. Perilla - 2018 
 *   
 ***************************************************************************/

#include <tcl.h>
#include "TclCommands.h"
#include "AtomSel.h"
#include "VMDApp.h"
#include "MoleculeList.h"
#include "Molecule.h"
#include "VolumetricData.h"
#include "VolMapCreate.h"
#include "QuickSurf.h"
#include <math.h>
#include "MeasureVolInterior.h"
#include "utilities.h"

#include <sstream>
#include <string>

#define MIN(X,Y) (((X)<(Y))? (X) : (Y))

typedef struct{
  const VolumetricData *volmapA;
  VolumetricData *targetVol;
  VolumetricData *pmap;
  float *targetPmap;
  float *targetMap;
  const float *testMap;
  const int *numvoxels;
  float isovalue;
  wkf_mutex_t mtx;
  const float *rayDir;
  long Nout;
} volinparms;


static void * volinthread_prob(void *voidparms) {
  wkf_tasktile_t tile;
  volinparms *parms = NULL;
  wkf_threadlaunch_getdata(voidparms, (void **) &parms);

  int gx,gy,gz;
  long Nout=0l;
  while (wkf_threadlaunch_next_tile(voidparms, 16384, &tile) != WKF_SCHED_DONE) { 
    int myx;
    int xsize = parms->numvoxels[0];
    int ysize = parms->numvoxels[1];
    int xysize = ysize*xsize; // precompute to save math later
    for (myx=tile.start; myx<tile.end; myx++) {
      gz = myx / (xysize);           // int divide/mod are SLOW (20-100+ clocks)
      gy = (myx % (xysize)) / xsize; // maybe make outer loop and inner loop
      gx = myx % xsize;              // so we iterate over consecutive X...

      // skip protein voxels (-5.0f) 
      float voxelA = parms->targetMap[myx];
      if (voxelA != PROTEINVOXEL) {
        int coord[3] =  {gx, gy, gz};
        int step[3];
        float rayOrig[3] = {gx + 0.5f, gy + 0.5f, gz + 0.5f};
        float deltaDist[3];
        float next[3];
    
        for (int ii = 0; ii < 3; ii++) {
          // XXX floating point division is slow, we should pre-compute the
          // rayDir[] / rayDir[..] fractions ahead of the main loop over voxels
          // in a separate tiny table, so we can just use them over and over 
          // without actually performing the divide, which is terribly slow on
          // most CPUs...  Then the ternary expression just becomes a select or
          // conditional move instruction that is way way faster if in-cache.
          // Further, it may be possible to precompute deltaDist via table too,
          // rather than having to call sqrtf(), which is also very slow.
          const float x = (parms->rayDir[ii] != 0.0f) ? (parms->rayDir[0] / parms->rayDir[ii]) : parms->rayDir[0];
          const float y = (parms->rayDir[ii] != 0.0f) ? (parms->rayDir[1] / parms->rayDir[ii]) : parms->rayDir[1];
          const float z = (parms->rayDir[ii] != 0.0f) ? (parms->rayDir[2] / parms->rayDir[ii]) : parms->rayDir[2];
          deltaDist[ii] = sqrtf(x*x + y*y + z*z);
          if (parms->rayDir[ii] < 0.f) {
            step[ii] = -1;
            next[ii] = (rayOrig[ii] - coord[ii]) * deltaDist[ii];
          } else {
            step[ii] = 1;
            next[ii] = (coord[ii] + 1.f - rayOrig[ii]) * deltaDist[ii];
          }
        }
    
        int hit=0;
        while(hit==0) {
          // Perform DDA
          int side=0;
          for (int ii=1; ii < 3; ii++) {
            if (next[side] > next[ii]) {
              side=ii;
            }
          }
          next[side]  += deltaDist[side];
          coord[side] += step[side];

          // Check if out of bounds
          if (coord[side] < 0 || coord[side] >= parms->numvoxels[side] ) {
            //parms->targetMap[myx]=EXTERIORVOXEL;
            Nout+=1l;
            break; //XXX
          }

          // Check if ray has hit a wall
          if (parms->testMap[coord[2]*xysize + coord[1]*xsize + coord[0]] >= parms->isovalue) {
            hit=1;
            parms->targetPmap[myx] += 1.0f;
          }
        } // wall detection loop
      }
    }
  }

  // atomic update of total Nout value
  // XXX we should consider using atomic increments instead
  wkf_mutex_lock(&parms->mtx); 
  parms->Nout += Nout;
  wkf_mutex_unlock(&parms->mtx);

  return NULL;
} // volinthread


static void * volinthread(void *voidparms) {
  wkf_tasktile_t tile;
  volinparms *parms = NULL;
  wkf_threadlaunch_getdata(voidparms, (void **) &parms);

  int gx,gy,gz;
  long Nout=0l;
  while (wkf_threadlaunch_next_tile(voidparms, 16384, &tile) != WKF_SCHED_DONE) { 
    int myx;
    int xsize = parms->numvoxels[0];
    int ysize = parms->numvoxels[1];
    int xysize = ysize*xsize; // precompute to save math later
    for (myx=tile.start; myx<tile.end; myx++) {
      gz = myx / (xysize);           // int divide/mod are SLOW (20-100+ clocks)
      gy = (myx % (xysize)) / xsize; // maybe make outer loop and inner loop
      gx = myx % xsize;              // so we iterate over consecutive X...

      // Find voxels inside but skip voxels that we already know to be 
      // inside from previous rounds
      float voxelA = parms->targetMap[myx];
      if (voxelA != EXTERIORVOXEL) {
        int coord[3] = {gx, gy, gz};
        int step[3];
        float rayOrig[3] = {gx + 0.5f, gy + 0.5f, gz + 0.5f} ;
        float deltaDist[3];
        float next[3];
    
        for (int ii = 0; ii < 3; ii++) {
          // XXX floating point division is slow, we should pre-compute the
          // rayDir[] / rayDir[..] fractions ahead of the main loop over voxels
          // in a separate tiny table, so we can just use them over and over 
          // without actually performing the divide, which is terribly slow on
          // most CPUs...  Then the ternary expression just becomes a select or
          // conditional move instruction that is way way faster if in-cache.
          // Further, it may be possible to precompute deltaDist via table too,
          // rather than having to call sqrtf(), which is also very slow.
          const float x = (parms->rayDir[ii] != 0.0f) ? (parms->rayDir[0] / parms->rayDir[ii]) : parms->rayDir[0];
          const float y = (parms->rayDir[ii] != 0.0f) ? (parms->rayDir[1] / parms->rayDir[ii]) : parms->rayDir[1];
          const float z = (parms->rayDir[ii] != 0.0f) ? (parms->rayDir[2] / parms->rayDir[ii]) : parms->rayDir[2];
          deltaDist[ii] = sqrtf(x*x + y*y + z*z);
          if (parms->rayDir[ii] < 0.f) {
            step[ii] = -1;
            next[ii] = (rayOrig[ii] - coord[ii]) * deltaDist[ii];
          } else {
            step[ii] = 1;
            next[ii] = (coord[ii] + 1.f - rayOrig[ii]) * deltaDist[ii];
          }
        }
    
        int hit=0;
        while(hit==0) {
          // Perform DDA
          int side=0;
          for (int ii=1; ii < 3; ii++) {
            if (next[side] > next[ii]) {
              side=ii;
            }
          }
          next[side]  += deltaDist[side];
          coord[side] += step[side];

          // Check if out of bounds
          if (coord[side] < 0 || coord[side] >= parms->numvoxels[side] ) {
            parms->targetMap[myx]=EXTERIORVOXEL;
            Nout+=1l;
            break;
          }

          // Check if ray has hit a wall
          if (parms->testMap[coord[2]*xysize + coord[1]*xsize + coord[0]] >= parms->isovalue) {
            hit=1;
          }
        } // wall detection loop
      }
    }
  }

  // atomic update of total Nout value
  // XXX we should consider using atomic increments instead
  wkf_mutex_lock(&parms->mtx); 
  parms->Nout += Nout;
  wkf_mutex_unlock(&parms->mtx);

  return NULL;
} // volinthread


long volin_threaded_prob(const VolumetricData *volmapA, VolumetricData *targetVol, 
                    VolumetricData *targetPvol, float _isovalue, float *rayDir) {
  volinparms parms;
  memset(&parms, 0, sizeof(parms));
  parms.testMap = volmapA->data;
  parms.targetMap = targetVol->data;
  parms.targetPmap = targetPvol->data;
  int numvoxels [] = {volmapA->xsize, volmapA->ysize, volmapA->zsize};
  parms.numvoxels = numvoxels;
  parms.targetVol = targetVol;
  parms.volmapA = volmapA;
  parms.isovalue = _isovalue;
  parms.rayDir = rayDir;
  parms.Nout=0;

  long Nout;
  int physprocs = wkf_thread_numprocessors();
  int maxprocs = physprocs;
  float *voltexmap = NULL;

  // We can productively use only a few cores per socket due to the
  // limited memory bandwidth per socket. Also, hyperthreading
  // actually hurts performance.  These two considerations combined
  // with the linear increase in memory use prevent us from using large
  // numbers of cores with this simple approach, so if we've got more 
  // than 8 CPU cores, we'll iteratively cutting the core count in 
  // half until we're under 20 cores.
  while (maxprocs > 40) 
    maxprocs /= 2;

  // Limit the number of CPU cores used so we don't run the 
  // machine out of memory during surface computation.
  // Use either a dynamic or hard-coded heuristic to limit the
  // number of CPU threads we will spawn so that we don't run
  // the machine out of memory.  
  long volsz = numvoxels[0] * numvoxels[1] * numvoxels[2];
  long volmemsz = sizeof(float) * volsz;
  long volmemszkb = volmemsz / 1024;
  long volmemtexszkb = volmemszkb + ((voltexmap != NULL) ? 3*volmemszkb : 0);

  // Platforms that don't have a means of determining available
  // physical memory will return -1, in which case we fall back to the
  // simple hard-coded 2GB-max-per-core heuristic.
  long vmdcorefree = -1;

#if defined(ARCH_BLUEWATERS) || defined(ARCH_CRAY_XC) || defined(ARCH_CRAY_XK) || defined(ARCH_LINUXAMD64) || defined(ARCH_SOLARIS2_64) || defined(ARCH_SOLARISX86_64) || defined(ARCH_AIX6_64) || defined(ARCH_MACOSXARM64) || defined(ARCH_MACOSXX86_64) 
  // XXX The core-free query scheme has one weakness in that we might have a 
  // 32-bit version of VMD running on a 64-bit machine, where the available
  // physical memory may be much larger than is possible for a 
  // 32-bit VMD process to address.  To do this properly we must therefore
  // use conditional compilation safety checks here until we  have a better
  // way of determining this with a standardized helper routine.
  vmdcorefree = vmd_get_avail_physmem_mb();
#endif

  if (vmdcorefree >= 0) {
    // Make sure QuickSurf uses no more than a fraction of the free memory
    // as an upper bound alternative to the hard-coded heuristic.
    // This should be highly preferable to the fixed-size heuristic
    // we had used in all cases previously.
    while ((volmemtexszkb * maxprocs) > (1024*vmdcorefree/4)) {
      maxprocs /= 2;
    }
  } else {
    // Set a practical per-core maximum memory use limit to 2GB, for all cores
    while ((volmemtexszkb * maxprocs) > (2 * 1024 * 1024))
      maxprocs /= 2;
  }

  if (maxprocs < 1) 
    maxprocs = 1;
  int numprocs = maxprocs; // ever the optimist
  wkf_mutex_init(&parms.mtx);
  wkf_tasktile_t tile;
  tile.start = 0;
  tile.end = volsz;
  wkf_threadlaunch(numprocs, &parms, volinthread_prob, &tile);
  wkf_mutex_destroy(&parms.mtx);
  Nout=parms.Nout;

  return Nout;
}


long volin_threaded(const VolumetricData *volmapA, VolumetricData *targetVol, 
                    float _isovalue, float *rayDir) {
  volinparms parms;
  memset(&parms, 0, sizeof(parms));
  parms.testMap = volmapA->data;
  parms.targetMap = targetVol->data;
  int numvoxels [] = {volmapA->xsize, volmapA->ysize, volmapA->zsize};
  parms.numvoxels = numvoxels;
  parms.targetVol = targetVol;
  parms.volmapA = volmapA;
  parms.isovalue = _isovalue;
  parms.rayDir = rayDir;
  parms.Nout=0;

  long Nout;
  int physprocs = wkf_thread_numprocessors();
  int maxprocs = physprocs;
  float *voltexmap = NULL;

  // We can productively use only a few cores per socket due to the
  // limited memory bandwidth per socket. Also, hyperthreading
  // actually hurts performance.  These two considerations combined
  // with the linear increase in memory use prevent us from using large
  // numbers of cores with this simple approach, so if we've got more 
  // than 8 CPU cores, we'll iteratively cutting the core count in 
  // half until we're under 20 cores.
  while (maxprocs > 40) 
    maxprocs /= 2;

  // Limit the number of CPU cores used so we don't run the 
  // machine out of memory during surface computation.
  // Use either a dynamic or hard-coded heuristic to limit the
  // number of CPU threads we will spawn so that we don't run
  // the machine out of memory.  
  long volsz = numvoxels[0] * numvoxels[1] * numvoxels[2];
  long volmemsz = sizeof(float) * volsz;
  long volmemszkb = volmemsz / 1024;
  long volmemtexszkb = volmemszkb + ((voltexmap != NULL) ? 3*volmemszkb : 0);

  // Platforms that don't have a means of determining available
  // physical memory will return -1, in which case we fall back to the
  // simple hard-coded 2GB-max-per-core heuristic.
  long vmdcorefree = -1;

#if defined(ARCH_BLUEWATERS) || defined(ARCH_CRAY_XC) || defined(ARCH_CRAY_XK) || defined(ARCH_LINUXAMD64) || defined(ARCH_SOLARIS2_64) || defined(ARCH_SOLARISX86_64) || defined(ARCH_AIX6_64) || defined(MARCOSXARM64) || defined(ARCH_MACOSXX86_64) 
  // XXX The core-free query scheme has one weakness in that we might have a 
  // 32-bit version of VMD running on a 64-bit machine, where the available
  // physical memory may be much larger than is possible for a 
  // 32-bit VMD process to address.  To do this properly we must therefore
  // use conditional compilation safety checks here until we  have a better
  // way of determining this with a standardized helper routine.
  vmdcorefree = vmd_get_avail_physmem_mb();
#endif

  if (vmdcorefree >= 0) {
    // Make sure QuickSurf uses no more than a fraction of the free memory
    // as an upper bound alternative to the hard-coded heuristic.
    // This should be highly preferable to the fixed-size heuristic
    // we had used in all cases previously.
    while ((volmemtexszkb * maxprocs) > (1024*vmdcorefree/4)) {
      maxprocs /= 2;
    }
  } else {
    // Set a practical per-core maximum memory use limit to 2GB, for all cores
    while ((volmemtexszkb * maxprocs) > (2 * 1024 * 1024))
      maxprocs /= 2;
  }

  if (maxprocs < 1) 
    maxprocs = 1;
  int numprocs = maxprocs; // ever the optimist
  wkf_mutex_init(&parms.mtx);
  wkf_tasktile_t tile;
  tile.start = 0;
  tile.end = volsz;
  wkf_threadlaunch(numprocs, &parms, volinthread, &tile);
  wkf_mutex_destroy(&parms.mtx);
  Nout=parms.Nout;

  return Nout;
}


VolumetricData* CreateEmptyGrid(const VolumetricData *volmapA) {
  VolumetricData *targetVol;
  float *volmap;
  // Calculate dimension of the grid
  long volsz=volmapA->xsize * volmapA->ysize * volmapA->zsize;
  // Allocate new grid
  volmap = new float[volsz];
  // Initialize new grid to 0
  memset(volmap, 0, sizeof(float)*volsz);  

  // Setup Volumetric data
  targetVol = new VolumetricData("inout map", volmapA->origin, 
                                 volmapA->xaxis, volmapA->yaxis, volmapA->zaxis,
                                 volmapA->xsize, volmapA->ysize, volmapA->zsize,
                                 volmap);
  return targetVol;
}


VolumetricData* CreateProbGrid(const VolumetricData *volmapA) {
  VolumetricData *targetPvol;
  float *volmap;
  // Calculate dimension of the grid
  long volsz=volmapA->xsize * volmapA->ysize * volmapA->zsize;
  // Allocate new grid
  volmap = new float[volsz];
  // Initialize new grid to 0
  memset(volmap, 0, sizeof(float)*volsz);  

  // Setup Volumetric data
  targetPvol = new VolumetricData("P(in) map", volmapA->origin, 
                                 volmapA->xaxis, volmapA->yaxis, volmapA->zaxis,
                                 volmapA->xsize, volmapA->ysize, volmapA->zsize,
                                 volmap);
  return targetPvol;
}


VolumetricData* normalize_pmap(const VolumetricData *volmapA, int nrays) {
  long gridsz = volmapA->gridsize();
  const float* map = volmapA->data;
  float *volmap;

  VolumetricData *targetPmap;

  long volsz = volmapA->xsize * volmapA->ysize * volmapA->zsize;
  volmap = new float[volsz];
  memset(volmap, 0, sizeof(float)*volsz);

  float nrays_inv = 1.0f / float(nrays);
  for (long l=0; l<gridsz; l++) {
    volmap[l] = map[l] * nrays_inv;
  }

  targetPmap = new VolumetricData("P(in) map", volmapA->origin,
                                 volmapA->xaxis, volmapA->yaxis, volmapA->zaxis,
                                 volmapA->xsize, volmapA->ysize, volmapA->zsize,
                                 volmap);
  return targetPmap;
}


void VolIn_CleanGrid(VolumetricData *volmapA) {
  long volsz=volmapA->xsize * volmapA->ysize * volmapA->zsize;
  memset(volmapA->data, 0, sizeof(float)*volsz);
}


long countIsoGrids(const VolumetricData *volmapA, const float _isovalue) {
  long volInIso=0;
  long gridsize = volmapA->gridsize();
  const float* map=volmapA->data;

  for (long l=0; l<gridsize; l++) {
    if (map[l] == _isovalue) volInIso += 1;
  }

  return volInIso;
}


long markIsoGrid(const VolumetricData *volmapA, VolumetricData *targetVol, 
                 const float _isovalue) {
  long nVoxels=0;

  long gridsize = volmapA->gridsize();
  const float *map=volmapA->data;
  float* targetmap=targetVol->data;

  for (long l=0; l<gridsize; l++) {
    if (map[l] >= _isovalue) {
      targetmap[l]=PROTEINVOXEL;
      nVoxels+=1;
    }
  }

  return nVoxels;
}


long RaycastGrid(const VolumetricData *volmapA, VolumetricData *targetVol,
                 float _isovalue, float *rayDir) {
  int volDimens[3]; // Volume dimensions

  // Calculate dimension of the grid
  float *map=targetVol->data;
  volDimens[0]=volmapA->xsize;
  volDimens[1]=volmapA->ysize; 
  volDimens[2]=volmapA->zsize;

  float isovalue=_isovalue; // Set isovalue for wall hit detection

  // Ray moves in gridspace 
  for (int k=0; k < volDimens[2]; k++) {
    for (int j=0; j < volDimens[1]; j++) {
      long voxrowaddr = k*volDimens[0]*volDimens[1] + j*volDimens[0];
      for (int i=0; i < volDimens[0]; i++) {
        long voxeladdr = voxrowaddr + i;
        if (map[voxeladdr] != EXTERIORVOXEL) {
          float deltaDist[3];
          float next[3];
          int coord[3] =  { i,j,k };
          int step[3];
          float rayOrig[3] = { i + 0.5f , j + 0.5f , k + 0.5f } ;
    
          for (int ii = 0; ii < 3; ii++) {
            const float x = (rayDir[0] / rayDir[ii]);
            const float y = (rayDir[1] / rayDir[ii]);
            const float z = (rayDir[2] / rayDir[ii]);
            deltaDist[ii] = sqrtf(x*x + y*y + z*z);
            if (rayDir[ii] < 0.f) {
              step[ii] = -1;
              next[ii] = (rayOrig[ii] - coord[ii]) * deltaDist[ii];
            } else {
              step[ii] = 1;
              next[ii] = (coord[ii] + 1.f - rayOrig[ii]) * deltaDist[ii];
            }
          }
    
          int hit=0;
          while (hit==0) {
            // Perform DDA
            int side=0;
            for (int ii=1; ii<3; ii++) {
              if (next[side] > next[ii]) {
                side=ii;
              }
            }
            next[side]  += deltaDist[side];
            coord[side] += step[side];
  
            // Check if out of bounds
            if (coord[side] < 0 || coord[side] >= volDimens[side] ) {
              map[voxeladdr]=EXTERIORVOXEL;
              break;
            }

            // Check if ray has hit a wall
            if (volmapA->data[coord[2]*volDimens[0]*volDimens[1] + coord[1]*volDimens[0] + coord[0]] > isovalue) {
              hit=1;
            }
          } // wall detection loop
        } // conditional
      } // loop over i
    } // loop over j
  } //loop ver k

  return 0;
}


VolumetricData* process_pmap (const VolumetricData *pmap, float cutoff) {
  long gridsz=pmap->gridsize();
  const float* map=pmap->data;
  float *volmap;

  VolumetricData *targetPmap;

  long volsz = pmap->xsize*pmap->ysize*pmap->zsize;
  volmap = new float[volsz];
  memset(volmap, 0, sizeof(float)*volsz);

  for (long l=0; l<gridsz; l++) {
    if (map[l]==0.0f) {
      volmap[l]=PROTEINVOXEL;
    } else if ((map[l]-cutoff) >= VOLMAPTOLERANCE) {
      volmap[l]=INTERIORVOXEL;
    } else if ((map[l]-cutoff) <= VOLMAPTOLERANCE) {
      volmap[l]=EXTERIORVOXEL;
    }
  }

  targetPmap = new VolumetricData("ProcMap", pmap->origin,
                                 pmap->xaxis, pmap->yaxis, pmap->zaxis,
                                 pmap->xsize, pmap->ysize, pmap->zsize,
                                 volmap);

  return targetPmap;
}


long vol_probability(const VolumetricData* probmap, float cutoff, float gspace) {
  long Nvox=0;
  long gridsz=probmap->gridsize();
  const float* map=probmap->data;

  for (long l=0; l<gridsz; l++) {
    if (map[l] >= cutoff) {
      Nvox++;
    }
  }
  return Nvox;
}


bool isfloat(char* opt) {
  std::istringstream iss(opt);
  float f;
  iss >> std::noskipws >> f;
  return iss.eof() && !iss.fail();
}


