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
*      $RCSfile: GraphLayout.h,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.3 $         $Date: 2020/05/29 16:48:55 $
*
***************************************************************************/
/**
 *  \file GraphLayout.h
 *  \brief Algorithms for 2-D and 3-D graph layouts needed for interactive
 *  clustering analysis tools, and similar usage.  Currently based on a
 *  modified variation of the Fructherman-Reingold layout algorithm.
 */

#include "ResizeArray.h"

class GraphLayout {
private:
  ResizeArray<float> pos_x;    // vertex x coord
  ResizeArray<float> pos_y;    // vertex y coord

  ResizeArray<float> dsp_x;    // vertex x displacement
  ResizeArray<float> dsp_y;    // vertex y displacement

  ResizeArray<int> edges_v0;   // explicit edge list
  ResizeArray<int> edges_v1;   // explicit edge list
  ResizeArray<float> weights;  // explicit edge weights

  const float *weight_matrix;  // external N^2 weight matrix

public:
  GraphLayout(int nverts, int nedges); // init size hints
  ~GraphLayout();   

  void add_edge(int v0, int v1, float weight);
  void add_weight_matrix(const float *weight_matrix);

  void init_positions_box();
  void init_positions_circle();

  void compute(int iters, float area, float kscale, 
               float tempscale, float distance_epsilon);

  void get_vertex_ptrs(int &numverts, const float * &px, const float * &py) {
    numverts = pos_x.num(); 
    px = &pos_x[0];
    py = &pos_y[0];
  }

};



