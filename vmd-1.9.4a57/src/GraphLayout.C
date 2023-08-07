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
*      $RCSfile: GraphLayout.C,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.8 $         $Date: 2020/10/15 16:07:31 $
*
***************************************************************************/
/**
 *  \file GraphLayout.C
 *  \brief Algorithms for 2-D and 3-D graph layouts needed for interactive
 *  clustering analysis tools, and similar usage.  Currently based on a 
 *  modified variation of the Fructherman-Reingold layout algorithm.
 */

// Enable the fully optimized code path, although it and the less
// optimized path have different strengths and weaknesses from a 
// numerical perspective.
#define VMDFROPTIMIZED 1

// Fruchterman-Reingold spring force graph layout  
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include "GraphLayout.h"
#include "utilities.h"   // needed for sincosf() on MacOS X

GraphLayout::GraphLayout(int numverts, int numedges) {
  pos_x.set_size(numverts);
  pos_y.set_size(numverts);
  dsp_x.set_size(numverts);
  dsp_y.set_size(numverts);

  if (numedges > 0) {
    edges_v0.set_size(numedges);
    edges_v1.set_size(numedges);
    weights.set_size(numedges);
  }

  weight_matrix = NULL;
}


GraphLayout::~GraphLayout() {
  pos_x.clear();
  pos_y.clear();
  dsp_x.clear();
  dsp_y.clear();

  edges_v0.clear();
  edges_v1.clear();
  weights.clear();
}


void GraphLayout::add_edge(int v0, int v1, float weight) {
  // If the two indices are the same, or if the weight is zero,
  // do nothing (no link)
  if (v0 == v1 || weight == 0.f) {
    return;
  }

  // auto-extend the vertex array if we get an edge with a larger
  // index than the current size.  This is easy-to-use but inefficient
  // compared to initializing the array sizes at construction time.
  long maxvert = std::max(v0, v1);
  if (maxvert > pos_x.num()) {
    long extra = maxvert - pos_x.num(); 
    pos_x.extend(extra);
    pos_y.extend(extra);
    dsp_x.extend(extra);
    dsp_y.extend(extra);
  }

  edges_v0.append(v0);
  edges_v1.append(v1);
  weights.append(weight);
}


void GraphLayout::add_weight_matrix(const float *weight_matrix) {
  this->weight_matrix = weight_matrix;
}


// place nodes randomly in a square centered at origin
void GraphLayout::init_positions_box() {
  int numverts = pos_x.num();
  float rm_inv = 1.0f / float(RAND_MAX);

  srand(314159); // make reproducible
  int v;
  for (v=0; v<numverts; v++) {
    pos_x[v] = rm_inv * rand() - 0.5f;
    pos_y[v] = rm_inv * rand() - 0.5f;
  }

  // clear displacement vectors
  memset(&dsp_x[0], 0, numverts * sizeof(float));
  memset(&dsp_y[0], 0, numverts * sizeof(float));
}


// place nodes randomly in a circle centered at origin
void GraphLayout::init_positions_circle() {
  int numverts = pos_x.num();
  float radscale = 0.5f;
  float angle = 0.0f;
  float delta_angle = 2.0f * VMD_PIF / numverts;

  int v;
  for (v=0; v<numverts; v++) {
    float sa, ca;
    sincosf(angle, &sa, &ca);
    pos_x[v] = radscale * ca;
    pos_y[v] = radscale * sa;
    angle += delta_angle;
  }
}


// A variant of the Fruchterman-Reingold layout algorithm, 
// with modifications to add support for weights, implicit N^2 all-to-all
// node/vertex connectivity, and to allow efficient CPU vectorization 
// on modern processors.  The original Fruchterman-Reingold algorithm
// is described in:
//  Graph Drawing by Force-Directed Placement. 
//  Fruchterman, T. M. J., and Reingold, E. M.,
//  Software: Practice and Experience, 21(11), 1991. 
//
void GraphLayout::compute(int iters, float area, float kscale, 
                          float tempscale, float distance_epsilon) {
  int numverts = pos_x.num();
  int numedges = edges_v0.num();

printf("GraphLayout::compute(iters: %d  verts: %d  edges: %d  area: %g  kscale: %g  tempscale: %g\n", 
       iters, numverts, numedges, area, kscale, tempscale);

  if (iters < 1) 
    return;

//  float area = 1.0f;
//  float kscale = 1.0e3f;
//  float tempscale = 0.2f;
//  float distance_epsilon = 1.0e-06f; // prevent division by zero

  float k2 = area / numverts;
  float k = sqrtf(k2) * kscale;
  float k_inv = 1.0f / k;
  float iters_inv = 1.0f / iters; 

  for (int i=0; i<iters; i++) {
printf("iter: %d\r", i); fflush(stdout);

    // Temperature cool down; starts at 1, ends at 0
    // (other formulas can be used for the cooling)
    float temp = (1.f -  i * iters_inv);

    // Rather than a purely linear temperature cooldown, we use a
    // quadratic curve, but it might be useful to allow a user-defined
    // exponent here, since it's an insignificant cost.
    // Higher exponents yield a cooldown with more iterations having
    // a lower temperature (smaller max displacement per step).
    float tempsquared = (temp * temp) * tempscale;

    //
    // calc repulsive force: Fr(d) = -k^2/d
    //
    int v0, v1;
    for (v0=0; v0<numverts; v0++) {
      float posx = pos_x[v0];
      float posy = pos_y[v0];
      float dispx = 0.0f;
      float dispy = 0.0f;

      // avoid per-iteration if test for v0!=v1, or (dx!=0, dy!=0)
      // by splitting the loop into two runs...
      for (v1=0; v1<v0; v1++) {
        float dx = posx - pos_x[v1];
        float dy = posy - pos_y[v1];
        float dxy2 = dx*dx + dy*dy + distance_epsilon; 
#if defined(VMDFROPTIMIZED)
        float repulsion = k2 / dxy2;
        dispx += dx * repulsion;
        dispy += dy * repulsion;
#else
        float d_1 = 1.0f / sqrtf(dxy2);
        float repulsion = k2 * d_1;
        dispx += dx * d_1 * repulsion;
        dispy += dy * d_1 * repulsion;
#endif
      }

      // avoid per-iteration if test for v0!=v1, or (dx!=0, dy!=0)
      // by splitting the loop into two runs...
      for (v1=v0+1; v1<numverts; v1++) {
        float dx = posx - pos_x[v1];
        float dy = posy - pos_y[v1];
        float dxy2 = dx*dx + dy*dy + distance_epsilon;
#if defined(VMDFROPTIMIZED)
        float repulsion = k2 / dxy2;
        dispx += dx * repulsion;
        dispy += dy * repulsion;
#else
        float d_1 = 1.0f / sqrtf(dxy2);
        float repulsion = k2 * d_1;
        dispx += dx * d_1 * repulsion;
        dispy += dy * d_1 * repulsion;
#endif
      }

      // write new displacements back...
      dsp_x[v0] = dispx;
      dsp_y[v0] = dispy;
    }

    //
    // calc attractive forces, only for connected vertices
    //   Fa(d) = d^2/k * weight -- weight factor is optional
    //
    // if the edges list is empty, assume a fully-connected graph
    if (numedges == 0) {
      if (weight_matrix == NULL) {
        // N^2/2 vertex combinations in fully-connected graph, no weights
        for (v0=0; v0<numverts; v0++) {
          float posx = pos_x[v0];
          float posy = pos_y[v0];
          float dispx = dsp_x[v0];
          float dispy = dsp_y[v0];

          for (v1=0; v1<v0; v1++) {
            float dx = posx - pos_x[v1];
            float dy = posy - pos_y[v1];
            float dxy2 = dx*dx + dy*dy;
#if defined(VMDFROPTIMIZED)
            float attraction = sqrtf(dxy2) * k_inv;
            float dxa = dx * attraction;
            float dya = dy * attraction;
#else
            float attraction = dxy2 * k_inv;
            float d_1 = 1.0f / sqrtf(dxy2 + distance_epsilon);
            float dxa = dx * d_1 * attraction;
            float dya = dy * d_1 * attraction;
#endif

            // memory scatter here is bad for performance...
            dispx     -= dxa;
            dispy     -= dya;
            dsp_x[v1] += dxa;
            dsp_y[v1] += dya;
          }

          // write new displacements back...
          dsp_x[v0] = dispx;
          dsp_y[v0] = dispy;
        }
      } else {
        // N^2/2 vertex combinations in fully-connected graph
        for (v0=0; v0<numverts; v0++) {
          float posx = pos_x[v0];
          float posy = pos_y[v0];
          float dispx = dsp_x[v0];
          float dispy = dsp_y[v0];

          for (v1=0; v1<v0; v1++) {
            // NxN weight matrix
            float weight = weight_matrix[v1*numverts + v0];

            float dx = posx - pos_x[v1];
            float dy = posy - pos_y[v1];
            float dxy2 = dx*dx + dy*dy;
#if defined(VMDFROPTIMIZED)
            float attraction = sqrtf(dxy2) * k_inv * weight;
            float dxa = dx * attraction;
            float dya = dy * attraction;
#else
            float attraction = dxy2 * k_inv * weight;
            float d_1 = 1.0f / sqrtf(dxy2 + distance_epsilon);
            float dxa = dx * d_1 * attraction;
            float dya = dy * d_1 * attraction;
#endif

            // memory scatter here is bad for performance...
            dispx     -= dxa;
            dispy     -= dya;
            dsp_x[v1] += dxa;
            dsp_y[v1] += dya;
          }

          // write new displacements back...
          dsp_x[v0] = dispx;
          dsp_y[v0] = dispy;
        }
      }
    } else {
      int e;
      for (e=0; e<numedges; e++) {
        int v0 = edges_v0[e];
        int v1 = edges_v1[e];
        float weight = weights[e];

        float dx = pos_x[v0] - pos_x[v1];
        float dy = pos_y[v0] - pos_y[v1];
        float dxy2 = dx*dx + dy*dy;
#if defined(VMDFROPTIMIZED)
        float attraction = sqrtf(dxy2) * k_inv * weight;
        float dxa = dx * attraction;
        float dya = dy * attraction;
#else
        float attraction = dxy2 * k_inv * weight;
        float d_1 = 1.0f / sqrtf(dxy2 + distance_epsilon);
        float dxa = dx * d_1 * attraction;
        float dya = dy * d_1 * attraction;
#endif

        // memory scatter here is bad for performance...
        dsp_x[v0] -= dxa;
        dsp_y[v0] -= dya;
        dsp_x[v1] += dxa;
        dsp_y[v1] += dya;
      }
    }

    //
    // Integration: update graph vertex positions with
    // the accumulated displacement vectors, with maximum
    // displacement limited by the current squared "temperature", 
    // per a very simple simulated annealing scheme.
    //
    for (v0=0; v0<numverts; v0++) {
      float dx = dsp_x[v0];
      float dy = dsp_y[v0];

      float dxy2 = dx*dx + dy*dy;

      // limit maximum displacement according to current 
      // temperature factor, by rescaling the displacement magnitude 
      // to match the squared temperature magnitude.
      if (dxy2 > tempsquared * tempsquared) {
        float d = sqrtf(dxy2);
        float rescale = tempsquared / d;
        dx *= rescale;
        dy *= rescale;
      }

      pos_x[v0] += dx;
      pos_y[v0] += dy;
    } 
  }
printf("\n");

}



