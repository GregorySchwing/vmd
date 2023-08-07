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
 *      $RCSfile: MeasureVolInterior.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $      $Date: 2019/09/27 00:40:24 $
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
VolumetricData* CreateEmptyGrid(const VolumetricData *);  
void VolInterior_CleanGrid(VolumetricData *);
long RaycastGrid(const VolumetricData *, VolumetricData *, float, float *);  
long volin_threaded(const VolumetricData *, VolumetricData *, float, float *);
long countIsoGrids(const VolumetricData *, const float);
long markIsoGrid(const VolumetricData *, VolumetricData *, const float);
VolumetricData* CreateProbGrid(const VolumetricData *);
VolumetricData* normalize_pmap(const VolumetricData *, int);
long volin_threaded_prob(const VolumetricData *, VolumetricData *, VolumetricData *, float, float *);
long vol_probability(const VolumetricData*,float,float);
bool isfloat(char*);
VolumetricData* process_pmap (const VolumetricData*, float);

#define EXTERIORVOXEL 5.0f
#define INTERIORVOXEL 0.0f
#define PROTEINVOXEL -5.0f
#define VOLMAPTOLERANCE 0.000000000001f

