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
 *	$RCSfile: CaveScene.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.55 $	$Date: 2020/02/26 15:51:10 $
 *
 ***************************************************************************/
/**
 *  \file CaveScene.C
 *  \brief Scene subclass that maintains a list of Displayable objects 
 *         in shared memory and draws them in parallel on a CaveDisplayDevice.
 *
 * The CAVE specific Scene that accesses all of its information from a
 * shared memory arena, since the display lists are shared amoung
 * multiple concurrently running renderer processes or threads.
 * The Scene has a list of Displayable objects and display commands.
 * The command lists are used to draw the objects, the Displayable
 * objects to prepare and update objects for drawing.
 */

#include <stdlib.h>
#include <cave_ogl.h>

#include "CaveScene.h"
#include "CaveRoutines.h"
#include "DisplayDevice.h"
#include "Inform.h"
#include "utilities.h"
#include "VMDApp.h" // for VMDexit();

#define CAVEUSERWLOCK 1

////////////////////////////  constructor  
CaveScene::CaveScene(VMDApp *vmdapp) : app(vmdapp) {
#if defined(CAVEUSERWLOCK)
  msgInfo << "Creating R/W lock for shared mem access sync.\n" << sendmsg;
  draw_rwlock = CAVENewLock();   // allocate a new CAVE lock
  CAVESetWriteLock(draw_rwlock); // start out locked for writing
#else 
  msgInfo << "Creating draw barrier for " <<  CAVEConfig->ActiveWalls
          << " active CAVE walls and master process.\n" << sendmsg; 

  // setup the drawing barrier used to control slave rendering processes
  draw_barrier =(wkf_barrier_t *) malloc_from_CAVE_memory(sizeof(wkf_barrier_t));
  if (draw_barrier == NULL) {
    msgErr << "CANNOT ALLOCATE SHARED MEMORY FOR CaveScene CLASS !!!"
           << sendmsg;
  }
  memset(draw_barrier, 0, sizeof(wkf_barrier_t));

  if (getenv("VMDCAVEDRAWPROCS") != NULL) {
    wkf_thread_barrier_init_proc_shared(draw_barrier, atoi(getenv("VMDCAVEDRAWPROCS")) + 1);
  } else {
    wkf_thread_barrier_init_proc_shared(draw_barrier, CAVEConfig->ActiveWalls + 1);
  }
#endif 
}

////////////////////////////  destructor  
CaveScene::~CaveScene(void) {
#if defined(CAVEUSERWLOCK)
  CAVEFreeLock(draw_rwlock); // delete the drawing rwlock
#else
  // free things allocated from shared memory
  free_to_CAVE_memory(draw_barrier);
#endif

}

// Called by CAVElib, on each of the child rendering processes
void CaveScene::draw(DisplayDevice *display) {
#if defined(CAVEUSERWLOCK)
  if (display->is_renderer_process()) {
    CAVESetReadLock(draw_rwlock);
    Scene::draw(display);                // draw the scene
    CAVEUnsetReadLock(draw_rwlock);
  } else {
    CAVEUnsetWriteLock(draw_rwlock);
    vmd_msleep(1); // give 'em 1 millisecond to draw before re-locking 
    CAVESetWriteLock(draw_rwlock);
  }
#else
  // Barrier synchronization for all drawing processes and the master
  wkf_thread_barrier(draw_barrier, 0); // wait for start of drawing 
  Scene::draw(display);                // draw the scene
  wkf_thread_barrier(draw_barrier, 0); // wait for end of drawing.
#endif
}

// called in VMDupdate, this updates the values of numDisplayable[23]D
// this is called by the parent!
int CaveScene::prepare() {
  // check if the ESC key is being pressed; if so, quit
  // note: THIS IS A HACK, just for Bill Sherman ;-)
  if(CAVEgetbutton(CAVE_ESCKEY)) {
    CAVEExit();
    app->VMDexit("Exiting due to Cave escape key being pressed.", 10, 4);
  }

  return Scene::prepare(); // call regular scene prepare method
}

void *CaveScene::operator new(size_t s) {
  return malloc_from_CAVE_memory(s);
}

void CaveScene::operator delete(void *p, size_t) {
  free_to_CAVE_memory(p);
}

