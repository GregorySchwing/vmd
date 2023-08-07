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
*      $RCSfile: ANARIRenderer.C,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.26 $         $Date: 2021/12/18 07:27:39 $
*
***************************************************************************/
/**
 *  \file ANARIRenderer.C
 *  \brief FileRenderer subclass for the Khronos ANARI rendering interface.
 *
 *  This code is based on early, incomplete, developmental versions of 
 *  ANARI header files and specification, and should not be used as a
 *  reference for developing ANARI applications.  This warning text will
 *  be removed when ANARI is finalized.
 * 
 *  Portions of this code are derived from Tachyon:
 *    "An Efficient Library for Parallel Ray Tracing and Animation"
 *    John E. Stone.  Master's Thesis, University of Missouri-Rolla,
 *    Department of Computer Science, April 1998
 * 
 *    "Rendering of Numerical Flow Simulations Using MPI"
 *    John Stone and Mark Underwood.
 *    Second MPI Developers Conference, pages 138-141, 1996.
 *    http://dx.doi.org/10.1109/MPIDC.1996.534105
 */

#define VMDANARIREPGROUPING 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(__linux)
#include <unistd.h>   // needed for symlink() in movie recorder
#endif

#include "Inform.h"
#include "ImageIO.h"
#include "ANARIRenderer.h"
#include "Matrix4.h"
#include "utilities.h"
#include "WKFUtils.h"

// #if !(ANARI_VERSION_MAJOR >= 1 && ANARI_VERSION_MINOR >= 2)
// #error VMD requires ANARI >= 1.2.0 for correct transparent AO shading
// // VMD requires ANARI >= 1.1.2 for correct rendering of cylinders
// #endif

// enable the interactive ray tracing capability
#if defined(VMDANARI_INTERACTIVE_OPENGL)
#if (defined(WIN32) || defined(_WIN64)) && defined(_MSC_VER)
#include <windows.h> // must include windows.h prior to GL
#endif

#include <GL/gl.h>
#endif

#if 0
#define DBG() 
#else
#define DBG() printf("ANARIRenderer) ++ %s\n", __func__);
#endif

#include <algorithm>

static void vmd_anari_status_callback(void *userData, 
                                      ANARIDevice dev,
                                      ANARIObject src,
                                      ANARIDataType srctype,
                                      ANARIStatusSeverity sev,
                                      ANARIStatusCode code,
                                      const char *message) {
  switch (sev) {
    case ANARI_SEVERITY_FATAL_ERROR:
      printf("ANARIRenderer) FATAL: %s\n", message);
      break;

    case ANARI_SEVERITY_ERROR:
      printf("ANARIRenderer) ERROR: %s\n", message);
      break;

    case ANARI_SEVERITY_WARNING:
      printf("ANARIRenderer) WARN: %s\n", message);
      break;

    case ANARI_SEVERITY_PERFORMANCE_WARNING:
      printf("ANARIRenderer) PERF: %s\n", message);
      break;

    case ANARI_SEVERITY_INFO:
      printf("ANARIRenderer) INFO: %s\n", message);
      break;

    default:
      printf("ANARIRenderer) STATUS: %s\n", message);
      break;
  }
}


// Global ANARI initialization routine -- call it only ONCE...
void ANARIRender::ANARI_Global_Init(void) {
//  DBG();
}


// Global ANARI initialization routine -- call it only ONCE...
void ANARIRender::ANARI_Global_Shutdown(void) {
}

/// prototypes for GL fctns
#if defined(VMDOPENGL)
extern "C" { 
typedef void ( *__GLXextFuncPtr)(void);
__GLXextFuncPtr glXGetProcAddress (const GLubyte *procName);
}
#endif

/// constructor ... initialize some variables
ANARIRender::ANARIRender(void) {
  DBG();

  anr_timer = wkf_timer_create(); // create and initialize timer
  wkf_timer_start(anr_timer);

  dev = NULL;
  rendererworkarounds = ANARI_NONE;
  memset(lastcommentstring, 0, sizeof(lastcommentstring));
  seqgrp = 0;

  // set ANARI state handles/variables to NULL
  anariRenderer = NULL;
  anariFrameBuffer = NULL;
  anariCamera = NULL;
  anariWorld = NULL;
  anariLightData = NULL;

  // zero out last rep counters
  lastrepmesh=0;
  lastrepspheres=0;
  lastrepcyls=0;

  lasterror = 0;               // begin with no error state set
  context_created = 0;         // no context yet
  buffers_allocated = 0;       // flag no buffer allocated yet
  scene_created = 0;           // scene has been created

  // clear timers
  time_ctx_create = 0.0;
  time_ctx_setup = 0.0;
  time_ctx_validate = 0.0;
  time_ctx_AS_build = 0.0;
  time_ray_tracing = 0.0;
  time_image_io = 0.0;

  // set default scene background state
  scene_background_mode = RT_BACKGROUND_TEXTURE_SOLID;
  memset(scene_bg_color, 0, sizeof(scene_bg_color));
  memset(scene_bg_grad_top, 0, sizeof(scene_bg_grad_top));
  memset(scene_bg_grad_bot, 0, sizeof(scene_bg_grad_bot));
  memset(scene_gradient, 0, sizeof(scene_gradient));
  scene_gradient_topval = 1.0f;
  scene_gradient_botval = 0.0f;
  // XXX this has to be recomputed prior to rendering..
  scene_gradient_invrange = 1.0f / (scene_gradient_topval - scene_gradient_botval);

  cam_zoom = 1.0f;
  cam_stereo_eyesep = 0.06f;
  cam_stereo_convergence_dist = 2.0f;

  headlight_enabled = 0;     // VR HMD headlight disabled by default

  shadows_enabled = RT_SHADOWS_OFF; // disable shadows by default 
  aa_samples = 0;            // no AA samples by default

  ao_samples = 0;            // no AO samples by default
  ao_direct = 0.3f;          // AO direct contribution is 30%
  ao_ambient = 0.7f;         // AO ambient contribution is 70%

  dof_enabled = 0;           // disable DoF by default
  cam_dof_focal_dist = 2.0f;
  cam_dof_fnumber = 64.0f;

  fog_mode = RT_FOG_NONE;    // fog/cueing disabled by default
  fog_start = 0.0f;
  fog_end = 10.0f;
  fog_density = 0.32f;

  rendererworkarounds = ANARI_NONE;
  anari_rendermode = ANARI_PATHTRACER;
  anari_matclass = ANARI_MATTE;

// XXX debugging attach delay
printf("debugging PID: %d\n", getpid());
if (getenv("VMDANARISLEEP")) {
  int sleepsecs = atoi(getenv("VMDANARISLEEP"));
  sleep(sleepsecs);
}


  if (getenv("VMDANARIUSD"))  {
    printf("ANARIRenderer) Attempting to load library: 'usd'\n");
    lib = anariLoadLibrary("usd", vmd_anari_status_callback, NULL);
    dev = anariNewDevice(lib, "usd");

    if (!dev) {
      printf("ANARIRenderer: failed to load USD ANARI library! Exiting!\n");
      exit(-1);
    }

    int usdConnLogVerbosity = 0;
    anariSetParameter(dev, dev, "usd::connection.logverbosity", ANARI_INT32, &usdConnLogVerbosity);

    int usdOutputOmniverse = (getenv("VMDUSEOMNIVERSE") != NULL);
    if (usdOutputOmniverse) {
      if (!getenv("VMDOMNIVERSEHOST")) {
        anariSetParameter(dev, dev, "usd::serialize.hostname", ANARI_STRING, "localhost");
      } else { 
        anariSetParameter(dev, dev, "usd::serialize.hostname", ANARI_STRING, getenv("VMDOMNIVERSEHOST"));
     }
    }

    if (!getenv("VMDUSDOUTPUTPATH")) {
      anariSetParameter(dev, dev, "usd::serialize.outputpath", ANARI_STRING, "vmdanariusd");
    } else {
      anariSetParameter(dev, dev, "usd::serialize.outputpath", ANARI_STRING, getenv("VMDUSDOUTPUTPATH"));
    }

    int usdOutputBinary = (getenv("VMDANARIUSDASCII") != NULL);
    anariSetParameter(dev, dev, "usd::serialize.outputbinary", ANARI_BOOL, &usdOutputBinary);

    anariCommit(dev, dev);

    rendererworkarounds = ANARI_USD;
  }


  // catch-all that reverts to the ref device if we do no better
  if (!dev) {
    const char *userdev = getenv("VMDANARIDEVICE");
    if (userdev) {
      printf("ANARIRenderer) attempting to load ANARI device '%s'\n", userdev); 
      lib = anariLoadLibrary(userdev, vmd_anari_status_callback, NULL);

    } else {
      printf("ANARIRenderer) No library loaded, trying 'example'\n"); 
      lib = anariLoadLibrary("example", vmd_anari_status_callback, NULL);
    }

    // query available subtypes for this device
    const char **devices = anariGetDeviceSubtypes(lib);
    if (!devices) {
      printf("ANARIRenderer) No device subtypes returned.\n");
    } else {
      printf("ANARIRenderer) Available devices:\n");
      for (const char **d = devices; *d != NULL; d++)
        printf("ANARIRenderer)   %s\n", *d);
    }

    // query available renderers
    int havepathtracer = 0;
    const char **renderers = anariGetObjectSubtypes(lib, "default", ANARI_RENDERER);
    if (renderers) {
      printf("ANARIRenderer) Available renderers:\n");
      for (const char **r = renderers; *r != NULL; r++) {
        printf("ANARIRenderer)   %s\n", *r);
        if (strcmp(*r, "pathtracer") == 0)
          havepathtracer = 1;
      }
    } else {
      printf("ANARIRenderer) No renderers available!\n");
    }

    const char *rendererstr = "default";
    if (havepathtracer)
      rendererstr = "pathtracer"; 

    // inspect pathtracer renderer parameters
    const ANARIParameter *rendparams = anariGetObjectParameters(lib, "default", "pathtracer", ANARI_RENDERER);

    if (!rendparams) {
      printf("ANARIRenderer) Renderer '%s' has no parameters.\n", rendererstr);
    } else {
      printf("ANARIRenderer) Parameters of '%s':\n", rendererstr);
      for (const ANARIParameter *p = rendparams; p->name != NULL; p++) {
        const char *desc = (const char *) anariGetParameterInfo(lib,
            "default",
            "pathtracer",
            ANARI_RENDERER,
            p->name,
            p->type,
            "description",
            ANARI_STRING);

        const int *required = (const int *) anariGetParameterInfo(lib,
            "default",
            "pathtracer",
            ANARI_RENDERER,
            p->name,
            p->type,
            "required",
            ANARI_BOOL);

        printf("ANARIRenderer)   [%d] %s, %s: %s\n",
               int(p->type), p->name, required && *required ? "REQ" : "OPT", desc);
      }
    }

    dev = anariNewDevice(lib, "default");

    if (!dev) {
      printf("ANARIRenderer: failed to load ANARI library! Exiting!\n");
      exit(-1);
    }

    anariCommit(dev, dev);

    rendererworkarounds = ANARI_NONE;
    anari_rendermode = ANARI_PATHTRACER;
  }


  if (getenv("VMDANARITRACE")) {
    // Trace output will be written to filename specified by env variable:
    //   ANARI_DEBUG_TRACE  
    anariLoadLibrary("debug", vmd_anari_status_callback, NULL); // enable ANARI API tracing
  }



#if 0
#if 1
  if (getenv("VMDANARITRACE")) {
    // Trace output will be written to filename specified by env variable:
    //   ANARI_DEBUG_TRACE  
    lib = anariLoadLibrary("debug", vmd_anari_status_callback, NULL); // enable ANARI API tracing
  }

  if (getenv("VMDANARINVGL"))  {
    printf("ANARIRenderer) Attempting to load library: 'nvgl'\n");
    lib = anariLoadLibrary("nvgl", vmd_anari_status_callback, NULL);
    dev = anariNewDevice(lib, "nvgl");
    rendererworkarounds = ANARI_NVGL;
#if defined(VMDOPENGL)
    const void *ptr = (const void*) glXGetProcAddress;
    anariSetParameter(dev, dev, "oglGetProcAddress", ANARI_VOID_PTR, &ptr);
#endif
  } else if (getenv("VMDANARIUSD"))  {
    printf("ANARIRenderer) Attempting to load library: 'usd'\n");
    lib = anariLoadLibrary("usd", vmd_anari_status_callback, NULL);
    dev = anariNewDevice(lib, "usd");
    anariSetParameter(dev, dev, "serialize.location", ANARI_STRING, "/tmp/foobar");
    rendererworkarounds = ANARI_USD;
  } else if (getenv("VMDANARIOSPRAY"))  {
    printf("ANARIRenderer) Attempting to load library: 'ospray'\n");
    lib = anariLoadLibrary("ospray", vmd_anari_status_callback, NULL);
    dev = anariNewDevice(lib, "ospray");
    rendererworkarounds = ANARI_OSPRAY;
  }
#endif

  if (!dev) {
    printf("ANARIRenderer) No devices loaded, trying 'example'\n"); 
    lib = anariLoadLibrary("example", vmd_anari_status_callback, NULL);
    dev = anariNewDevice(lib, "example");
  }

  if (getenv("VMDANARIOWL") == NULL) {
    int loglevel = ANARI_LOG_INFO;
    loglevel = ANARI_LOG_DEBUG;
    printf("ANARIRenderer) setting device log level property...\n");
    anariSetParameter(dev, dev, "logLevel", ANARI_INT32, &loglevel);
  }
  anariCommit(dev, dev);

  printf("ANARIRenderer) setting renderer mode\n");
  anari_rendermode = ANARI_PATHTRACER;
  if (getenv("VMDANARISCIVIS")) {
    printf("ANARIRenderer) Renderer mode set to 'scivis'\n");
    anari_rendermode = ANARI_SCIVIS;
  }
  if (getenv("VMDANARIPATHTRACER")) {
    printf("ANARIRenderer) Renderer mode set to 'pathtracer'\n");
    anari_rendermode = ANARI_PATHTRACER;
  }
#endif



  // verbose = RT_VERB_MIN;    // quiet console except for perf/debugging cases
  verbose = RT_VERB_DEBUG;  // extensive console output
  check_verbose_env();      // see if the user has overridden verbose flag

  printf("ANARIRenderer) clear/init scene data structures\n");
  destroy_scene();          // clear/init geometry vectors
  anariInstances.clear();   // clear instance list

  init_materials();

  directional_lights.clear();
  positional_lights.clear();

  // clear all primitive lists
  trimesh_v3f_n3f_c3f.clear();
  spheres_color.clear();
  cylinders_color.clear();
  surfbufs.clear();

  double starttime = wkf_timer_timenow(anr_timer);

  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating context...\n");

  //
  // create the main renderer object needed early on for 
  // instantiation of materials, lights, etc.
  // 
  const char *rstr;
  switch (anari_rendermode) {
    case ANARI_PATHTRACER: rstr = "pathtracer";        break;
    case ANARI_DEFAULT:    rstr = "default";           break;
    case ANARI_SCIVIS:     rstr = "scivis";            break;
    case ANARI_AO:         rstr = "ambientocclusion";  break;
              default:     rstr = "default";           break;
  }

  if (rendererworkarounds == ANARI_NVGL) {
    rstr = "default";
    anari_rendermode = ANARI_DEFAULT;
    printf("ANARIRenderer) INFO: NVGL back-end requires default renderer\n");
  }

  if ((anariRenderer = anariNewRenderer(dev, rstr)) == NULL) {
    printf("ANARIRenderer) Failed to load renderer '%s'!\n", rstr);
  }
  if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) created renderer '%s'\n", rstr);
  }

  time_ctx_create = wkf_timer_timenow(anr_timer) - starttime;
  
  if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) context creation time: %.2f\n", time_ctx_create);
  }

  context_created = 1;
}

        
/// destructor
ANARIRender::~ANARIRender(void) {
  DBG();

  destroy_scene();

  if (context_created && (anariRenderer != NULL))
    anariRelease(dev, anariRenderer);

  anariRenderer=NULL;
  context_created=0;

  if (dev != NULL)
    anariRelease(dev, dev);
  dev=NULL;

  // how do we free a library? or do we?
  lib=NULL;

  materialcache.clear();

  directional_lights.clear();
  positional_lights.clear();

  wkf_timer_destroy(anr_timer);
}


void ANARIRender::check_verbose_env() {
  DBG();

  char *verbstr = getenv("VMDANARIVERBOSE");
  if (verbstr != NULL) {
//    printf("ANARIRenderer) verbosity config request: '%s'\n", verbstr);
    if (!strupcmp(verbstr, "MIN")) {
      verbose = RT_VERB_MIN;
      printf("ANARIRenderer) verbose setting: minimum\n");
    } else if (!strupcmp(verbstr, "TIMING")) {
      verbose = RT_VERB_TIMING;
      printf("ANARIRenderer) verbose setting: timing data\n");
    } else if (!strupcmp(verbstr, "DEBUG")) {
      verbose = RT_VERB_DEBUG;
      printf("ANARIRenderer) verbose setting: full debugging data\n");
    }
  }
}


void ANARIRender::setup_context(int w, int h) {
  DBG();
  double starttime = wkf_timer_timenow(anr_timer);
  time_ctx_setup = 0;

  lasterror = 0; /* XXX SUCCESS; */ // clear any error state
  width = w;
  height = h;

  if (!context_created)
    return;

  check_verbose_env(); // update verbose flag if changed since last run

  // maxPathLength -- supported by all renderers
  if (getenv("VMDANARIMAXDEPTH")) {
    int maxdepth = atoi(getenv("VMDANARIMAXDEPTH"));
    if (maxdepth > 0 && maxdepth <= 20) {
      printf("ANARIRenderer) Setting maxdepth to %d...\n", maxdepth);
      anariSetParameter(dev, anariRenderer, "maxPathLength", ANARI_INT32, &maxdepth);  
    } else {
      printf("ANARIRenderer) ignoring out-of-range maxdepth: %d...\n", maxdepth);
    }
  } else {
    int maxdepth = 20;
    anariSetParameter(dev, anariRenderer, "maxPathLength", ANARI_INT32, &maxdepth);  
  }
 
#if 0
  // XXX -- not implemented in ANARI 2.x presently
  // Implore ANARI to correctly handle lighting through transparent 
  // surfaces when AO is enabled
  const int one = 1;
  anariSetParameter(dev, anariRenderer, "aoTransparencyEnabled", ANARI_INT32, &one);
#endif

  time_ctx_setup = wkf_timer_timenow(anr_timer) - starttime;
}


void ANARIRender::destroy_scene() {
  DBG();

  double starttime = wkf_timer_timenow(anr_timer);
  time_ctx_destroy_scene = 0;

  // zero out all object counters
  cylinder_array_cnt = 0;
  cylinder_array_color_cnt = 0;
  ring_array_color_cnt = 0;
  sphere_array_cnt = 0;
  sphere_array_color_cnt = 0;
  tricolor_cnt = 0;
  trimesh_c4u_n3b_v3f_cnt = 0;
  trimesh_n3b_v3f_cnt = 0;
  trimesh_n3f_v3f_cnt = 0;
  trimesh_v3f_cnt = 0;

  // zero out last rep counters
  lastrepmesh=0;
  lastrepspheres=0;
  lastrepcyls=0;

  // frame depends on lots of other objects,
  // so we destroy it early
  framebuffer_destroy();

  // clear lists of primitives
  int i;
  for (i=0; i<trimesh_v3f_n3f_c3f.num(); i++) {
    free(trimesh_v3f_n3f_c3f[i].v);
    trimesh_v3f_n3f_c3f[i].v = NULL;
    free(trimesh_v3f_n3f_c3f[i].n);
    trimesh_v3f_n3f_c3f[i].n = NULL;
    free(trimesh_v3f_n3f_c3f[i].c);
    trimesh_v3f_n3f_c3f[i].c = NULL;
    free(trimesh_v3f_n3f_c3f[i].f);
    trimesh_v3f_n3f_c3f[i].f = NULL;
  }
  trimesh_v3f_n3f_c3f.clear();

  for (i=0; i<spheres_color.num(); i++) {
    free(spheres_color[i].xyz);
    spheres_color[i].xyz = NULL;
    free(spheres_color[i].radii);
    spheres_color[i].radii = NULL;
    free(spheres_color[i].colors);
    spheres_color[i].colors = NULL;
  }
  spheres_color.clear();

  for (i=0; i<cylinders_color.num(); i++) {
    free(cylinders_color[i].cyls);
    cylinders_color[i].cyls = NULL;
    free(cylinders_color[i].rads);
    cylinders_color[i].rads = NULL;
    free(cylinders_color[i].ind);
    cylinders_color[i].ind = NULL;
    free(cylinders_color[i].cols);
    cylinders_color[i].cols = NULL;
  }
  cylinders_color.clear();

  for (i=0; i<surfbufs.num(); i++) {
    free(surfbufs[i]);
    surfbufs[i] = NULL;
  }
  surfbufs.clear();

  anariInstances.clear();

  for (i=0; i<materialcache.num(); i++) {
    if (materialcache[i].isvalid)
      anariRelease(dev, materialcache[i].mat);
  }
  materialcache.clear();

  int lcnt = anariLights.num();
  for (i = 0; i < lcnt; ++i) {
    anariRelease(dev, anariLights[i]);
  }
  anariLights.clear();


  if (anariCamera != NULL) {
    anariRelease(dev, anariCamera);
    anariCamera = NULL;
  }

  if (anariWorld != NULL) {
    anariRelease(dev, anariWorld);
    anariWorld = NULL;
  }


  double endtime = wkf_timer_timenow(anr_timer);
  time_ctx_destroy_scene = endtime - starttime;

  scene_created = 0; // scene has been destroyed
}


void ANARIRender::update_rendering_state(int interactive) {
  DBG();
  if (!context_created)
    return;

  wkf_timer_start(anr_timer);

  // Set interactive/progressive rendering flag so that we wire up
  // the most appropriate renderer for the task.  For batch rendering
  // with AO, we would choose the largest possible sample batch size,
  // but for interactive we will always choose a batch size of 1 or maybe 2
  // to yield the best interactivity.
  interactive_renderer = interactive;

  // XXX set ANARI rendering state

  long totaltris = tricolor_cnt + trimesh_c4u_n3b_v3f_cnt + 
                   trimesh_n3b_v3f_cnt + trimesh_n3f_v3f_cnt + trimesh_v3f_cnt;

  if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) cyl %ld, ring %ld, sph %ld, tri %ld, tot: %ld  lt %ld\n",
           cylinder_array_cnt + cylinder_array_color_cnt,
           ring_array_color_cnt,
           sphere_array_cnt + sphere_array_color_cnt,
           totaltris,
           cylinder_array_cnt +  cylinder_array_color_cnt + ring_array_color_cnt + sphere_array_cnt + sphere_array_color_cnt + totaltris,
           directional_lights.num() + positional_lights.num());
  }

  if (verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) using fully general shader and materials.\n");
  }

  // XXX set ANARI background color

  if (verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) scene bg mode: %d\n", scene_background_mode);

    printf("ANARIRenderer) scene bgsolid: %.2f %.2f %.2f\n", 
           scene_bg_color[0], scene_bg_color[1], scene_bg_color[2]);

    printf("ANARIRenderer) scene bggradT: %.2f %.2f %.2f\n", 
           scene_bg_grad_top[0], scene_bg_grad_top[1], scene_bg_grad_top[2]);

    printf("ANARIRenderer) scene bggradB: %.2f %.2f %.2f\n", 
           scene_bg_grad_bot[0], scene_bg_grad_bot[1], scene_bg_grad_bot[2]);
  
    printf("ANARIRenderer) bg gradient: %f %f %f  top: %f  bot: %f\n",
           scene_gradient[0], scene_gradient[1], scene_gradient[2],
           scene_gradient_topval, scene_gradient_botval);
  }

  // update in case the caller changed top/bottom values since last recalc
  scene_gradient_invrange = 1.0f / (scene_gradient_topval - scene_gradient_botval);
  // XXX set ANARI background gradient

  // XXX set ANARI fog mode

  if (verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) adding lights: dir: %ld  pos: %ld\n", 
           directional_lights.num(), positional_lights.num());
  }

  // XXX set ANARI lights

  if (verbose == RT_VERB_DEBUG) 
    printf("ANARIRenderer) Finalizing ANARI scene graph...\n");

  // create group to hold instances

  // XXX we should create an acceleration object the instance shared
  //     by multiple PBC images


  // XXX ANARI AS builder initialization if there's any customization...

  // do final state variable updates before rendering begins
  if (verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) cam zoom factor %f\n", cam_zoom);
    printf("ANARIRenderer) cam stereo eye separation  %f\n", cam_stereo_eyesep);
    printf("ANARIRenderer) cam stereo convergence distance %f\n", 
           cam_stereo_convergence_dist);
    printf("ANARIRenderer) cam DoF focal distance %f\n", cam_dof_focal_dist);
    printf("ANARIRenderer) cam DoF f/stop %f\n", cam_dof_fnumber);
  }

  // define all of the standard camera params
  // XXX set ANARI camera state

  // define stereoscopic camera parameters
  // XXX set ANARI camera state

  // define camera DoF parameters
  // XXX set ANARI camera state

  // XXX set ANARI AO sample counts and light scaling factors

  if (verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) setting sample counts:  AA %d  AO %d\n", aa_samples, ao_samples);
    printf("ANARIRenderer) setting AO factors:  AOA %f  AOD %f\n", ao_ambient, ao_direct);
  }

  //
  // Handle AA samples either internally with loops internal to 
  // each ray launch point thread, or externally by iterating over
  // multiple launches, adding each sample to an accumulation buffer,
  // or a hybrid combination of the two.  
  //
#if 1
  ext_aa_loops = std::max(aa_samples, 1) * std::max(ao_samples, 1);
#else
  ext_aa_loops = 1;
  if (ao_samples > 0 || (aa_samples > 4)) {
    // if we have too much work for a single-pass rendering, we need to 
    // break it up into multiple passes of the right counts in each pass
    ext_aa_loops = 1 + aa_samples;
    // XXX set ANARI sample counts per launch...
  } else { 
    // if the scene is simple, e.g. no AO rays and AA sample count is small,
    // we can run it in a single pass and get better performance
    // XXX set ANARI sample counts per launch...
  }
  // XXX set ANARI accum buf normalization scaling factors
#endif

  if (verbose == RT_VERB_DEBUG) {
    if (ext_aa_loops > 1)
      printf("ANARIRenderer) Running ANARI multi-pass: %d loops\n", ext_aa_loops);
    else
      printf("ANARIRenderer) Running ANARI single-pass: %d total samples\n", 1+aa_samples);
  }

  // set the ray generation program to the active camera code...
  // XXX set ANARI camera mode and clear accum buf
  // set the active color accumulation ray gen program based on the 
  // camera/projection mode, stereoscopic display mode, 
  // and depth-of-field state
  // XXX set ANARI camera mode and accum buf mode
  // XXX set ANARI "miss" shading mode (solid or gradient)
}


void ANARIRender::framebuffer_config(int fbwidth, int fbheight) {
  DBG();
  if (!context_created)
    return;

  width = fbwidth;
  height = fbheight;

  // allocate and resize buffers to match request
  if (buffers_allocated) {
    // if the buffers already exist and match the current 
    // progressive/non-progressive rendering mode, just resize them
    if (verbose == RT_VERB_DEBUG) {
      printf("ANARIRenderer) resizing framebuffer\n");
    }
    framebuffer_resize(width, height);
  } else {
    // (re)allocate framebuffer and associated accumulation buffers if they
    // don't already exist or if they weren't bound properly for
    // current progressive/non-progressive rendering needs.
    if (verbose == RT_VERB_DEBUG) {
      printf("ANARIRenderer) creating framebuffer and accum. buffer\n");
    }

    // XXX for ANARI framebuffer setup is completely deferred to render-time

    buffers_allocated = 1;
  }
}


void ANARIRender::framebuffer_resize(int fbwidth, int fbheight) {
  DBG();
  if (!context_created)
    return;

  width = fbwidth;
  height = fbheight;

  if (buffers_allocated) {
    if (verbose == RT_VERB_DEBUG) 
      printf("ANARIRenderer) framebuffer_resize(%d x %d)\n", width, height);
    framebuffer_destroy();
  }

  // XXX for ANARI framebuffer setup is completely deferred to render-time

  buffers_allocated = 1;
}


void ANARIRender::framebuffer_destroy() {
  DBG();
  if (!context_created)
    return;

//  if (buffers_allocated) {
    // XXX for ANARI framebuffer teardown is done just after unmapping
//  }
  buffers_allocated = 0;
}



void ANARIRender::commit_rep() {
#if defined(VMDANARIREPGROUPING)
  int i;
  int nmeshbufs = trimesh_v3f_n3f_c3f.num() - lastrepmesh;
  int nspherebufs = spheres_color.num() - lastrepspheres;
  int ncylbufs = cylinders_color.num() - lastrepcyls;
  int numsurfs = nmeshbufs + nspherebufs + ncylbufs;

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Committing rep, %d surfs: %d meshbufs, %d spbufs, %d cylbufs\n", numsurfs, nmeshbufs, nspherebufs, ncylbufs);

  if (numsurfs < 1)
    return;

  ANARISurface *surfs = (ANARISurface *) calloc(1, numsurfs * sizeof(ANARISurface));
  surfbufs.append(surfs);
  int cursurf=0;

  // commit triangle mesh geometry after assigning materials
  for (i=lastrepmesh; i<trimesh_v3f_n3f_c3f.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("ANARIRenderer) Adding triangle mesh[%d]: %d tris ...\n", 
             i, trimesh_v3f_n3f_c3f[i].num);
 
    surfs[cursurf++] = trimesh_v3f_n3f_c3f[i].surf;
  } 
  lastrepmesh=trimesh_v3f_n3f_c3f.num();


  // commit sphere geometry after assigning materials
  for (i=lastrepspheres; i<spheres_color.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("ANARIRenderer) Adding sphere_color array [%d]: %d spheres ...\n",
             i, spheres_color[i].num);

    surfs[cursurf++] = spheres_color[i].surf;
  } 
  lastrepspheres=spheres_color.num();


  // commit cylinder geometry after assigning materials
  for (i=lastrepcyls; i<cylinders_color.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("ANARIRenderer) Adding cylinders_color array [%d]: %d cylinders...\n",
             i, cylinders_color[i].num);

    surfs[cursurf++] = cylinders_color[i].surf;
  }
  lastrepcyls=cylinders_color.num();


  // add all of the prims to the same group/instance
  ANARIGroup group = anariNewGroup(dev);
#if 1
printf("ANARIRenderer) Assigning group name.\n");
  // if using the USD back-end, assign the "name" tag for the geom...
  if (rendererworkarounds == ANARI_USD && (strlen(lastcommentstring) > 0)) {
    char strbuf[2048];
    sprintf(strbuf, "%s_seqgrp%d", lastcommentstring, seqgrp);
    anariSetParameter(dev, group, "name", ANARI_STRING, strbuf);
    seqgrp++; 
  }
#endif

printf("ANARIRenderer) Generating surface array...\n");
  ANARIArray1D surfobj = anariNewArray1D(dev, &surfs[0], 0, 0, ANARI_SURFACE, numsurfs, 0);
  anariCommit(dev, surfobj);

printf("ANARIRenderer) Generating group.\n");
  anariSetParameter(dev, group, "surface", ANARI_ARRAY1D, &surfobj);

printf("ANARIRenderer) Committing group.\n");
  anariCommit(dev, group); 
  anariRelease(dev, surfobj);

  for (i=0; i<numsurfs; i++) {
    anariRelease(dev, surfs[i]);
  }

  ANARIInstance instance = anariNewInstance(dev);
  anariSetParameter(dev, instance, "group", ANARI_GROUP, &group);
  anariCommit(dev, instance);
  anariRelease(dev, group);

  anariInstances.append(instance);
#endif
}


void ANARIRender::render_compile_and_validate(void) {
  int i;

  DBG();
  if (!context_created)
    return;

  //
  // finalize context validation, compilation, and AS generation 
  //
  double startctxtime = wkf_timer_timenow(anr_timer);

  // XXX any last ANARI state updates/checks

  // commit the final bit of accumulated rep geometry if not already...
  commit_rep();

  if ((anariWorld = anariNewWorld(dev)) == NULL) {
    printf("ANARIRenderer) Failed to create new world!\n");
  }

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) num spheres = %ld\n", spheres_color.num());


  // 
  // Set camera parms
  // 
  float cam_pos_orig[3] = {0.0f, 0.0f, 2.0f};
  float cam_U_orig[3] = {1.0f, 0.0f, 0.0f};
  float cam_V_orig[3] = {0.0f, 1.0f, 0.0f};
  float cam_W_orig[3] = {0.0f, 0.0f, -1.0f};

  float cam_pos[3], cam_U[3], cam_V[3], cam_W[3];
  vec_copy(cam_pos, cam_pos_orig);
  vec_copy(cam_U, cam_U_orig);
  vec_copy(cam_V, cam_V_orig);
  vec_copy(cam_W, cam_W_orig);

  if (camera_projection == ANARIRender::RT_ORTHOGRAPHIC) {
    if(!anariCamera) anariCamera = anariNewCamera(dev, "orthographic");

    float orthoheight = 2.0f * cam_zoom;
    anariSetParameter(dev, anariCamera, "height", ANARI_FLOAT32, &orthoheight);

    if (dof_enabled) {
      msgWarn << "ANARIRenderer) DoF not implemented for orthographic camera!" << sendmsg;
    }
  } else {
    if(!anariCamera) anariCamera = anariNewCamera(dev, "perspective");

    float camfovy = 2.0f*180.0f*(atanf(cam_zoom)/M_PI);
    anariSetParameter(dev, anariCamera, "fovy", ANARI_FLOAT32, &camfovy);

    if (dof_enabled) {
      anariSetParameter(dev, anariCamera, "focusDistance", ANARI_FLOAT32, &cam_dof_focal_dist);
      float dofaprad = cam_dof_focal_dist / (2.0f * cam_zoom * cam_dof_fnumber);
      anariSetParameter(dev, anariCamera, "apertureRadius", ANARI_FLOAT32, &dofaprad);
    } else {
      float dofaprad = 0.0f;
      anariSetParameter(dev, anariCamera, "apertureRadius", ANARI_FLOAT32, &dofaprad);
    }
  }

  if (anariCamera) {
    float camaspect = width / ((float) height); 
    anariSetParameter(dev, anariCamera, "aspect", ANARI_FLOAT32, &camaspect);
    anariSetParameter(dev, anariCamera, "position",  ANARI_FLOAT32_VEC3, cam_pos);
    anariSetParameter(dev, anariCamera, "direction", ANARI_FLOAT32_VEC3, cam_W);
    anariSetParameter(dev, anariCamera, "up",        ANARI_FLOAT32_VEC3, cam_V);
    anariCommit(dev, anariCamera);
  }

  // 
  // Set framebuffer 
  // 
  framebuffer_config(width, height);

  //
  // Set all lights
  //
  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) setting lights...\n");

  // The direct lighting scaling factor all of the other lights.
  float lightscale = 1.0f;
  if (ao_samples != 0)
    lightscale = ao_direct;

  for (i = 0; i < directional_lights.num(); ++i) {
    ANARILight light = anariNewLight(dev, "directional");
    anariSetParameter(dev, light, "color", ANARI_FLOAT32_VEC3, directional_lights[i].color);

//      // The direct lighting scaling factor is applied to the lights here.
//      anariSetParameter(dev, light, "intensity", ANARI_FLOAT32, &lightscale);

    // ANARI uses a light direction vector opposite to VMD and Tachyon 
    float lightDir[3];
    vec_negate(lightDir, directional_lights[i].dir);
    vec_normalize(lightDir); // just for good measure
    anariSetParameter(dev, light, "direction", ANARI_FLOAT32_VEC3, lightDir);
    anariCommit(dev, light);
    anariLights.append(light);
  }

  // 
  // update renderer state
  //

  // AO scaling factor is applied to a special ambient light.
  if (ao_samples != 0) {
    // ANARI spec puts AO on renderers, not as a light object
    if (getenv("VMDANARIDEVICE") && !strcmp(getenv("VMDANARIDEVICE"), "visrtx")) {
      float tmp = ao_ambient * 1.5f;
      anariSetParameter(dev, anariRenderer, "ambientLight", ANARI_FLOAT32, &tmp);
    } else {
      ANARILight light = anariNewLight(dev, "ambient");
  
      // AO scaling factor is applied to the special ambient light
      anariSetParameter(dev, light, "intensity", ANARI_FLOAT32, &ao_ambient);
 
      float whitecol[] = { 1.0f, 1.0f, 1.0f };
      anariSetParameter(dev, light, "color", ANARI_FLOAT32_VEC3, whitecol);
  
      anariCommit(dev, light);
      anariLights.append(light); // add AO ambient light
    }
  } 

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) setting sample counts...\n");

  // ANARI uses VEC4F only for background colors
  float bgcoltmp[4];
  vec_copy(bgcoltmp, scene_bg_color);
  bgcoltmp[3] = 1.0f;
  anariSetParameter(dev, anariRenderer, "backgroundColor", ANARI_FLOAT32_VEC4, bgcoltmp);

  if (rendererworkarounds != ANARI_USD) {
    if (ao_samples && interactive_renderer) {
      const int one = 1;
      anariSetParameter(dev, anariRenderer, "pixelSamples", ANARI_INT32, &one); // all renderers
      if (anari_rendermode == ANARI_SCIVIS)
        anariSetParameter(dev, anariRenderer, "aoSamples", ANARI_INT32, &one);    // scivis-only
    } else {
      anariSetParameter(dev, anariRenderer, "pixelSamples", ANARI_INT32, &aa_samples); // all renderers
      if (anari_rendermode == ANARI_SCIVIS)
        anariSetParameter(dev, anariRenderer, "aoSamples", ANARI_INT32, &ao_samples);    // scivis-only
    }

    if (getenv("VMDANARIAOMAXDIST")) {
      float tmp = atof(getenv("VMDANARIAOMAXDIST"));
      if (verbose == RT_VERB_DEBUG) {
        printf("ANARIRenderer) setting AO maxdist: %f\n", tmp);
      }
      anariSetParameter(dev, anariRenderer, "aoRadius", ANARI_FLOAT32, &tmp); // scivis-only
    }
  }

  if (rendererworkarounds != ANARI_USD) {
#if 1
    // XXX ANARI doesn't support rendering w/o shadows presently
    // render with/without shadows
    msgInfo << "Shadow rendering enabled." << sendmsg;
#else
    if (shadows_enabled || ao_samples) {
      if (shadows_enabled && !ao_samples)
        msgInfo << "Shadow rendering enabled." << sendmsg;

      const int one = 1;
      anariSetParameter(dev, anariRenderer, "shadowsEnabled", ANARI_INT32, &one);
    } else {
      const int zero = 0;
      anariSetParameter(dev, anariRenderer, "shadowsEnabled", ANARI_INT32, &zero);
    }
#endif
  }

  // render with ambient occlusion, but only if shadows are also enabled
  if (ao_samples) {
    msgInfo << "Ambient occlusion enabled." << sendmsg;
//    msgInfo << "Shadow rendering enabled." << sendmsg;
  }

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Committing geometry buffers...\n");

#if !defined(VMDANARIREPGROUPING)
  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Committing %ld trimesh buffers...\n", trimesh_v3f_n3f_c3f.num());

  // commit triangle mesh geometry after assigning materials
  for (i=0; i<trimesh_v3f_n3f_c3f.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("ANARIRenderer) Adding triangle mesh[%d]: %d tris ...\n", 
             i, trimesh_v3f_n3f_c3f[i].num);
 
    ANARIGroup group = anariNewGroup(dev);
#if 1
    // if using the USD back-end, assign the "name" tag for the geom...
    if (rendererworkarounds == ANARI_USD && (strlen(lastcommentstring) > 0)) {
 //     printf("ANARIRenderer)    Mesh: '%s'\n", lastcommentstring);

      char strbuf[2048];
      sprintf(strbuf, "%s_seqgrp%d", lastcommentstring, seqgrp);
      anariSetParameter(dev, group, "name", ANARI_STRING, strbuf);
      seqgrp++; 
    }
#endif
    ANARIArray1D surfs = anariNewArray1D(dev, &trimesh_v3f_n3f_c3f[i].surf, 0, 0, ANARI_SURFACE, 1, 0);
    anariCommit(dev, surfs);
    anariRelease(dev, trimesh_v3f_n3f_c3f[i].surf);
    anariSetParameter(dev, group, "surface", ANARI_ARRAY1D, &surfs);
    anariCommit(dev, group);
    anariRelease(dev, surfs);

    ANARIInstance instance = anariNewInstance(dev);
    anariSetParameter(dev, instance, "group", ANARI_GROUP, &group);
    anariCommit(dev, instance);
    anariRelease(dev, group);

    anariInstances.append(instance);
  } 

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Committing %ld sphere buffers...\n", spheres_color.num());

  // commit sphere geometry after assigning materials
  for (i=0; i<spheres_color.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("ANARIRenderer) Adding sphere_color array [%d]: %d spheres ...\n",
             i, spheres_color[i].num);

    ANARIGroup group = anariNewGroup(dev);
#if 1
    // if using the USD back-end, assign the "name" tag for the geom...
    if (rendererworkarounds == ANARI_USD && (strlen(lastcommentstring) > 0)) {
//      printf("ANARIRenderer)    SphereArray: '%s'\n", lastcommentstring);

      char strbuf[2048];
      sprintf(strbuf, "%s_seqgrp%d", lastcommentstring, seqgrp);
      anariSetParameter(dev, group, "name", ANARI_STRING, strbuf);
      seqgrp++; 
    }
#endif
    ANARIArray1D surfs = anariNewArray1D(dev, &spheres_color[i].surf, 0, 0, ANARI_SURFACE, 1, 0);
    anariCommit(dev, surfs);
    anariRelease(dev, spheres_color[i].surf);
    anariSetParameter(dev, group, "surface", ANARI_ARRAY1D, &surfs);
    anariCommit(dev, group); 
    anariRelease(dev, surfs);

    ANARIInstance instance = anariNewInstance(dev);
    anariSetParameter(dev, instance, "group", ANARI_GROUP, &group);
    anariCommit(dev, instance);
    anariRelease(dev, group);

    anariInstances.append(instance);
  } 

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Committing %ld cylinder buffers...\n", cylinders_color.num());

  // commit cylinder geometry after assigning materials
  for (i=0; i<cylinders_color.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("ANARIRenderer) Adding cylinders_color array [%d]: %d cylinders...\n",
             i, cylinders_color[i].num);

    ANARIGroup group = anariNewGroup(dev);
#if 1
    // if using the USD back-end, assign the "name" tag for the geom...
    if (rendererworkarounds == ANARI_USD && (strlen(lastcommentstring) > 0)) {
//      printf("ANARIRenderer)    CylinderArray: '%s'\n", lastcommentstring);

      char strbuf[2048];
      sprintf(strbuf, "%s_seqgrp%d", lastcommentstring, seqgrp);
      anariSetParameter(dev, group, "name", ANARI_STRING, strbuf);
      seqgrp++; 
    }
#endif
    ANARIArray1D surfs = anariNewArray1D(dev, &cylinders_color[i].surf, 0, 0, ANARI_SURFACE, 1, 0);
    anariCommit(dev, surfs);
    anariRelease(dev, cylinders_color[i].surf);
    anariSetParameter(dev, group, "surface", ANARI_ARRAY1D, &surfs);
    anariCommit(dev, group); 
    anariRelease(dev, surfs);

    ANARIInstance instance = anariNewInstance(dev);
    anariSetParameter(dev, instance, "group", ANARI_GROUP, &group);
    anariCommit(dev, instance);
    anariRelease(dev, group);

    anariInstances.append(instance);
  }
#endif


  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Attaching instances to scene...\n");

  // attach all instances to the scene...
  ANARIArray1D instances = anariNewArray1D(dev, &anariInstances[0], 0, 0, ANARI_INSTANCE, anariInstances.num(), 0);
  anariCommit(dev, instances);
  anariSetParameter(dev, anariWorld, "instance", ANARI_ARRAY1D, &instances);
  anariRelease(dev, instances);
  for (i=0; i<anariInstances.num(); i++) {
    anariRelease(dev, anariInstances[i]);
  }
  anariInstances.clear();

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Attaching %ld lights to scene...\n", anariLights.num());

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Committing anariWorld...\n");

  // if using the USD back-end, assign the "name" tag for the geom...
  if (rendererworkarounds == ANARI_USD && (strlen(lastcommentstring) > 0)) {
    anariSetParameter(dev, anariWorld, "name", ANARI_STRING, "VMD World");
  }

  if (anariLights.num() > 0) {
    // attach all lights to the scene...
    ANARIArray1D lights = anariNewArray1D(dev, &anariLights[0], 0, 0, ANARI_LIGHT, anariLights.num(), 0);
    anariCommit(dev, lights);
    anariSetParameter(dev, anariWorld, "light", ANARI_ARRAY1D, &lights);
    anariCommit(dev, anariWorld); // commit the completed scene
    anariRelease(dev, lights);
  } else {
    anariCommit(dev, anariWorld); // commit the completed scene
  }

  // print out world bounds
  float worldBounds[6] = {};
  if (anariGetProperty(dev, anariWorld, "bounds", ANARI_FLOAT32_BOX3, 
                       worldBounds, sizeof(worldBounds), ANARI_WAIT)) {
    printf("ANARIRenderer) world bounds: ({%f, %f, %f}, {%f, %f, %f}\n\n",
           worldBounds[0], worldBounds[1], worldBounds[2],
           worldBounds[3], worldBounds[4], worldBounds[5]);
  }

  if (verbose == RT_VERB_DEBUG)
    printf("ANARIRenderer) Committing anariRenderer...\n");

  anariCommit(dev, anariRenderer);

  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) Finalizing ANARI rendering kernels...\n");
  // XXX any last ANARI state updates/checks

  double contextinittime = wkf_timer_timenow(anr_timer);
  time_ctx_validate = contextinittime - startctxtime;

  //
  // Force ANARI to build the acceleration structure _now_, so we can time it
  //
  // XXX No way to force-build ANARI AS for timing?

  time_ctx_AS_build = wkf_timer_timenow(anr_timer) - contextinittime;
  if (verbose == RT_VERB_DEBUG) {
    printf("ANARIRenderer) launching render: %d x %d\n", width, height);
  }
}


#if defined(VMDANARI_INTERACTIVE_OPENGL)

static void *createanariraywindow(const char *wintitle, int width, int height) {
  printf("ANARIRenderer) Creating ANARI window: %d x %d...\n", width, height);

  void *win = glwin_create(wintitle, width, height);
  while (glwin_handle_events(win, GLWIN_EV_POLL_NONBLOCK) != 0);

  glDrawBuffer(GL_BACK);
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glClearColor(0.0, 0.0, 0.0, 1.0); /* black */
  glViewport(0, 0, width, height);
  glClear(GL_COLOR_BUFFER_BIT);

  glShadeModel(GL_FLAT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, width, height, 0.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glDrawBuffer(GL_BACK);
  glClear(GL_COLOR_BUFFER_BIT);

  glwin_swap_buffers(win);

  return win;
}


static void interactive_viewer_usage(void *win) {
  printf("ANARIRenderer) VMD TachyonL-ANARI Interactive Ray Tracer help:\n");
  printf("ANARIRenderer) ================================================\n");

  // check for Spaceball/SpaceNavigator/Magellan input devices
  int havespaceball = ((glwin_spaceball_available(win)) && (getenv("VMDDISABLESPACEBALLXDRV") == NULL));
  printf("ANARIRenderer) Spaceball/SpaceNavigator/Magellan: %s\n",
         (havespaceball) ? "Available" : "Not available");

  // check for stereo-capable display
  int havestereo, havestencil;
  glwin_get_wininfo(win, &havestereo, &havestencil);
  printf("ANARIRenderer) Stereoscopic display: %s\n",
         (havestereo) ? "Available" : "Not available");

  // check for vertical retrace sync
  int vsync=0, rc=0;
  if ((rc = glwin_query_vsync(win, &vsync)) == GLWIN_SUCCESS) {
    printf("ANARIRenderer) Vert retrace sync: %s\n", (vsync) ? "On" : "Off");
  } else {
    printf("ANARIRenderer) Vert retrace sync: indeterminate\n");
  }

  printf("ANARIRenderer)\n");
  printf("ANARIRenderer) General controls:\n");
  printf("ANARIRenderer)   space: save numbered snapshot image\n");
  printf("ANARIRenderer)       =: reset to initial view\n");
  printf("ANARIRenderer)       h: print this help info\n");
  printf("ANARIRenderer)       p: print current rendering parameters\n");
  printf("ANARIRenderer)   ESC,q: quit viewer\n");
  printf("ANARIRenderer)\n");
  printf("ANARIRenderer) Display controls\n");
  printf("ANARIRenderer)      F1: override shadows on/off (off=AO off too)\n");
  printf("ANARIRenderer)      F2: override AO on/off\n");
  printf("ANARIRenderer)      F3: override DoF on/off\n");
  printf("ANARIRenderer)      F4: override Depth cueing on/off\n");
// Not currently applicable to ANARI
// #ifdef USE_REVERSE_SHADOW_RAYS
//   printf("ANARIRenderer)      F5: enable/disable shadow ray optimizations\n");
// #endif
  printf("ANARIRenderer)     F12: toggle full-screen display on/off\n");
  printf("ANARIRenderer)   1-9,0: override samples per update auto-FPS off\n");
  printf("ANARIRenderer)      Up: increase DoF focal distance\n");
  printf("ANARIRenderer)    Down: decrease DoF focal distance\n");
  printf("ANARIRenderer)    Left: decrease DoF f/stop\n");
  printf("ANARIRenderer)   Right: increase DoF f/stop\n");
  printf("ANARIRenderer)       S: toggle stereoscopic display on/off (if avail)\n");
  printf("ANARIRenderer)       a: toggle AA/AO auto-FPS tuning on/off (on)\n");
  printf("ANARIRenderer)       g: toggle gradient sky xforms on/off (on)\n");
  printf("ANARIRenderer)       l: toggle light xforms on/off (on)\n");
  printf("ANARIRenderer)\n");
  printf("ANARIRenderer) Mouse controls:\n");
  printf("ANARIRenderer)       f: mouse depth-of-field mode\n");
  printf("ANARIRenderer)       r: mouse rotation mode\n");
  printf("ANARIRenderer)       s: mouse scaling mode\n");
  printf("ANARIRenderer)       t: mouse translation mode\n");

  int movie_recording_enabled = (getenv("VMDANARILIVEMOVIECAPTURE") != NULL);
  if (movie_recording_enabled) {
    printf("ANARIRenderer)\n");
    printf("ANARIRenderer) Movie recording controls:\n");
    printf("ANARIRenderer)       R: start/stop movie recording\n");
    printf("ANARIRenderer)       F: toggle movie FPS (24, 30, 60)\n");
  }
}


void ANARIRender::render_to_glwin(const char *filename) {
  DBG();
  int i;

  if (!context_created)
    return;

  enum RtMouseMode { RTMM_ROT=0, RTMM_TRANS=1, RTMM_SCALE=2, RTMM_DOF=3 };
  enum RtMouseDown { RTMD_NONE=0, RTMD_LEFT=1, RTMD_MIDDLE=2, RTMD_RIGHT=3 };
  RtMouseMode mm = RTMM_ROT;
  RtMouseDown mousedown = RTMD_NONE;

  // flags to interactively enable/disable shadows, AO, DoF
  int gl_shadows_on=(shadows_enabled) ? RT_SHADOWS_ON : RT_SHADOWS_OFF;

  int gl_fs_on=0; // fullscreen window state
  int owsx=0, owsy=0; // store last win size before fullscreen
  int gl_ao_on=(ao_samples > 0);
  int gl_dof_on, gl_dof_on_old;
  gl_dof_on=gl_dof_on_old=dof_enabled; 
  int gl_fog_on=(fog_mode != RT_FOG_NONE);

  // Enable live recording of a session to a stream of image files indexed
  // by their display presentation time, mapped to the nearest frame index
  // in a fixed-frame-rate image sequence (e.g. 24, 30, or 60 FPS), 
  // to allow subsequent encoding into a standard movie format.
  // XXX this feature is disabled by default at present, to prevent people
  //     from accidentally turning it on during a live demo or the like
  int movie_recording_enabled = (getenv("VMDANARILIVEMOVIECAPTURE") != NULL);
  int movie_recording_on = 0;
  double movie_recording_start_time = 0.0;
  int movie_recording_fps = 30;
  int movie_framecount = 0;
  int movie_lastframeindex = 0;
  const char *movie_recording_filebase = "vmdlivemovie.%05d.tga";
  if (getenv("VMDANARILIVEMOVIECAPTUREFILEBASE"))
    movie_recording_filebase = getenv("VMDANARILIVEMOVIECAPTUREFILEBASE");

  // Enable/disable Spaceball/SpaceNavigator/Magellan input 
  int spaceballenabled=(getenv("VMDDISABLESPACEBALLXDRV") == NULL) ? 1 : 0;
  int spaceballmode=0;       // default mode is rotation/translation
  int spaceballflightmode=0; // 0=moves object, 1=camera fly
  if (getenv("VMDANARISPACEBALLFLIGHT"))
    spaceballflightmode=1;


  // total AA/AO sample count
  int totalsamplecount=0;

  // counter for snapshots of live image...
  int snapshotcount=0;

  // flag to enable automatic AO sample count adjustment for FPS rate control
  int autosamplecount=1;

  // flag to enable transformation of lights and gradient sky sphere, 
  // so that they track camera orientation as they do in the VMD OpenGL display
  int xformlights=1, xformgradientsphere=1;

  //
  // allocate or reconfigure the framebuffer, accumulation buffer, 
  // and output streams required for progressive rendering, either
  // using the new progressive APIs, or using our own code.
  //
  // Unless overridden by environment variables, we use the incoming
  // window size parameters from VMD to initialize the RT image dimensions.
  // If image size is overridden, often when using HMDs, the incoming 
  // dims are window dims are used to size the GL window, but the image size
  // is set independently.
  int wsx=width, wsy=height;
  const char *imageszstr = getenv("VMDANARIIMAGESIZE");
  if (imageszstr) {
    if (sscanf(imageszstr, "%d %d", &width, &height) != 2) {
      width=wsx;
      height=wsy;
    } 
  } 
  framebuffer_config(width, height);

  // prepare the majority of ANARI rendering state before we go into 
  // the interactive rendering loop
  update_rendering_state(1);
  render_compile_and_validate();

  // make a copy of state we're going to interactively manipulate,
  // so that we can recover to the original state on-demand
  int samples_per_pass = 1;
  int cur_aa_samples = aa_samples;
  int cur_ao_samples = ao_samples;
  float cam_zoom_orig = cam_zoom;
  float scene_gradient_orig[3] = {0.0f, 1.0f, 0.0f};
  vec_copy(scene_gradient_orig, scene_gradient);

  float cam_pos_orig[3] = {0.0f, 0.0f, 2.0f};
  float cam_U_orig[3] = {1.0f, 0.0f, 0.0f};
  float cam_V_orig[3] = {0.0f, 1.0f, 0.0f};
  float cam_W_orig[3] = {0.0f, 0.0f, -1.0f};
  float cam_pos[3], cam_U[3], cam_V[3], cam_W[3];
  float hmd_U[3], hmd_V[3], hmd_W[3];

  vec_copy(cam_pos, cam_pos_orig);
  vec_copy(cam_U, cam_U_orig);
  vec_copy(cam_V, cam_V_orig);
  vec_copy(cam_W, cam_W_orig);

  // copy light directions
  anr_directional_light *cur_dlights = (anr_directional_light *) calloc(1, directional_lights.num() * sizeof(anr_directional_light));
  for (i=0; i<directional_lights.num(); i++) {
    vec_copy((float*)&cur_dlights[i].color, directional_lights[i].color);
    vec_copy((float*)&cur_dlights[i].dir, directional_lights[i].dir);
    vec_normalize((float*)&cur_dlights[i].dir);
  }

  // create the display window
  void *win = createanariraywindow("VMD TachyonL-ANARI Interactive Ray Tracer", width, height);
  interactive_viewer_usage(win);
  
  // check for stereo-capable display
  int havestereo=0, havestencil=0;
  int stereoon=0, stereoon_old=0;
  glwin_get_wininfo(win, &havestereo, &havestencil);

  // Override AA/AO sample counts since we're doing progressive rendering.
  // Choosing an initial AO sample count of 1 will give us the peak progressive 
  // display update rate, but we end up wasting time on re-tracing many
  // primary rays.  The automatic FPS optimization scheme below will update
  // the number of samples per rendering pass and assign the best values for
  // AA/AO samples accordingly.
  cur_aa_samples = samples_per_pass;
  if (cur_ao_samples > 0) {
    cur_aa_samples = 1;
    cur_ao_samples = samples_per_pass;
  }

  const char *statestr = "|/-\\.";
  int done=0, winredraw=1, accum_count=0;
  int state=0, mousedownx=0, mousedowny=0;
  float cur_cam_zoom = cam_zoom_orig;

  double fpsexpave=0.0; 
  
  double oldtime = wkf_timer_timenow(anr_timer);
  while (!done) { 
    int winevent=0;

    while ((winevent = glwin_handle_events(win, GLWIN_EV_POLL_NONBLOCK)) != 0) {
      int evdev, evval;
      char evkey;

      glwin_get_lastevent(win, &evdev, &evval, &evkey);
      glwin_get_winsize(win, &wsx, &wsy);

      if (evdev == GLWIN_EV_WINDOW_CLOSE) {
        printf("ANARIRenderer) display window closed, exiting...\n");
        done = 1;
        winredraw = 0;
      } else if (evdev == GLWIN_EV_KBD) {
        switch (evkey) {
          case  '1': autosamplecount=0; samples_per_pass=1; winredraw=1; break;
          case  '2': autosamplecount=0; samples_per_pass=2; winredraw=1; break;
          case  '3': autosamplecount=0; samples_per_pass=3; winredraw=1; break;
          case  '4': autosamplecount=0; samples_per_pass=4; winredraw=1; break;
          case  '5': autosamplecount=0; samples_per_pass=5; winredraw=1; break;
          case  '6': autosamplecount=0; samples_per_pass=6; winredraw=1; break;
          case  '7': autosamplecount=0; samples_per_pass=7; winredraw=1; break;
          case  '8': autosamplecount=0; samples_per_pass=8; winredraw=1; break;
          case  '9': autosamplecount=0; samples_per_pass=9; winredraw=1; break;
          case  '0': autosamplecount=0; samples_per_pass=10; winredraw=1; break;

          case  '=': /* recover back to initial state */
            vec_copy(scene_gradient, scene_gradient_orig);
            cam_zoom = cam_zoom_orig;
            vec_copy(cam_pos, cam_pos_orig);
            vec_copy(cam_U, cam_U_orig);
            vec_copy(cam_V, cam_V_orig);
            vec_copy(cam_W, cam_W_orig);

            // restore original light directions
            for (i=0; i<directional_lights.num(); i++) {
              vec_copy((float*)&cur_dlights[i].dir, directional_lights[i].dir);
              vec_normalize((float*)&cur_dlights[i].dir);
            }
            winredraw = 1;
            break;
 
          case  ' ': /* spacebar saves current image with counter */
            {
              char snapfilename[256];
              sprintf(snapfilename, "vmdsnapshot.%04d.tga", snapshotcount);
              const unsigned char *FB = (const unsigned char*)anariMapFrame(dev, anariFrameBuffer, "color");
              if (write_image_file_rgb4u(snapfilename, FB, width, height)) {
                printf("ANARIRenderer) Failed to write output image!\n");
              } else {
                printf("ANARIRenderer) Saved snapshot to '%s'             \n",
                       snapfilename);
              }
              anariUnmapFrame(dev, anariFrameBuffer, "color");
              snapshotcount++; 
            }
            break;

          case  'a': /* toggle automatic sample count FPS tuning */
            autosamplecount = !(autosamplecount);
            printf("\nANARIRenderer) Automatic AO sample count FPS tuning %s\n",
                   (autosamplecount) ? "enabled" : "disabled");
            break;

          case  'f': /* DoF mode */
            mm = RTMM_DOF;
            printf("\nANARIRenderer) Mouse DoF aperture and focal dist. mode\n");
            break;

          case  'g': /* toggle gradient sky sphere xforms */
            xformgradientsphere = !(xformgradientsphere);
            printf("\nANARIRenderer) Gradient sky sphere transformations %s\n",
                   (xformgradientsphere) ? "enabled" : "disabled");
            break;

          case  'h': /* print help message */
            printf("\n");
            interactive_viewer_usage(win);
            break;

          case  'l': /* toggle lighting xforms */
            xformlights = !(xformlights);
            printf("\nANARIRenderer) Light transformations %s\n",
                   (xformlights) ? "enabled" : "disabled");
            break;

          case  'p': /* print current RT settings */
            printf("\nANARIRenderer) Current Ray Tracing Parameters:\n"); 
            printf("ANARIRenderer) -------------------------------\n"); 
            printf("ANARIRenderer) Camera zoom: %f\n", cur_cam_zoom);
            printf("ANARIRenderer) Shadows: %s  Ambient occlusion: %s\n",
                   (gl_shadows_on) ? "on" : "off",
                   (gl_ao_on) ? "on" : "off");
            printf("ANARIRenderer) Antialiasing samples per-pass: %d\n",
                   cur_aa_samples);
            printf("ANARIRenderer) Ambient occlusion samples per-pass: %d\n",
                   cur_ao_samples);
            printf("ANARIRenderer) Depth-of-Field: %s f/num: %.1f  Foc. Dist: %.2f\n",
                   (gl_dof_on) ? "on" : "off", 
                   cam_dof_fnumber, cam_dof_focal_dist);
            printf("ANARIRenderer) Image size: %d x %d\n", width, height);
            break;

          case  'r': /* rotate mode */
            mm = RTMM_ROT;
            printf("\nANARIRenderer) Mouse rotation mode\n");
            break;

          case  's': /* scaling mode */
            mm = RTMM_SCALE;
            printf("\nANARIRenderer) Mouse scaling mode\n");
            break;

          case  'F': /* toggle live movie recording FPS (24, 30, 60) */
            if (movie_recording_enabled) {
              switch (movie_recording_fps) {
                case 24: movie_recording_fps = 30; break;
                case 30: movie_recording_fps = 60; break;
                case 60:
                default: movie_recording_fps = 24; break;
              }
              printf("\nANARIRenderer) Movie recording FPS rate: %d\n", 
                     movie_recording_fps);
            } else {
              printf("\nANARIRenderer) Movie recording not available.\n");
            }
            break;

          case  'R': /* toggle live movie recording mode on/off */
            if (movie_recording_enabled) {
              movie_recording_on = !(movie_recording_on);
              printf("\nANARIRenderer) Movie recording %s\n",
                     (movie_recording_on) ? "STARTED" : "STOPPED");
              if (movie_recording_on) {
                movie_recording_start_time = wkf_timer_timenow(anr_timer);
                movie_framecount = 0;
                movie_lastframeindex = 0;
              } else {
                printf("ANARIRenderer) Encode movie with:\n");
                printf("ANARIRenderer)   ffmpeg -f image2 -i vmdlivemovie.%%05d.tga -c:v libx264 -profile:v baseline -level 3.0 -pix_fmt yuv420p -b:v 15000000 output.mp4\n");
              }
            } else {
              printf("\nANARIRenderer) Movie recording not available.\n");
            }
            break;

          case  'S': /* toggle stereoscopic display mode */
            if (havestereo) {
              stereoon = (!stereoon);
              printf("\nANARIRenderer) Stereoscopic display %s\n",
                     (stereoon) ? "enabled" : "disabled");
              winredraw = 1;
            } else {
              printf("\nANARIRenderer) Stereoscopic display unavailable\n");
            }
            break;
 
          case  't': /* translation mode */
            mm = RTMM_TRANS;
            printf("\nANARIRenderer) Mouse translation mode\n");
            break;
            
          case  'q': /* 'q' key */
          case  'Q': /* 'Q' key */
          case 0x1b: /* ESC key */
            printf("\nANARIRenderer) Exiting on user input.               \n");
            done=1; /* exit from interactive RT window */
            break;
        }
      } else if (evdev != GLWIN_EV_NONE) {
        switch (evdev) {
          case GLWIN_EV_KBD_F1: /* turn shadows on/off */
            gl_shadows_on=(!gl_shadows_on) ? RT_SHADOWS_ON : RT_SHADOWS_OFF;
            // gl_shadows_on = (!gl_shadows_on);
            printf("\n");
            printf("ANARIRenderer) Shadows %s\n",
                   (gl_shadows_on) ? "enabled" : "disabled");
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_F2: /* turn AO on/off */
            gl_ao_on = (!gl_ao_on); 
            printf("\n");
            printf("ANARIRenderer) Ambient occlusion %s\n",
                   (gl_ao_on) ? "enabled" : "disabled");
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_F3: /* turn DoF on/off */
            gl_dof_on = (!gl_dof_on);
            printf("\n");
            if ((camera_projection == RT_ORTHOGRAPHIC) && gl_dof_on) {
              gl_dof_on=0; 
              printf("ANARIRenderer) Depth-of-field not available in orthographic mode\n");
            }
            printf("ANARIRenderer) Depth-of-field %s\n",
                   (gl_dof_on) ? "enabled" : "disabled");
            winredraw = 1;
            break;

          case GLWIN_EV_KBD_F4: /* turn fog/depth cueing on/off */
            gl_fog_on = (!gl_fog_on); 
            printf("\n");
            printf("ANARIRenderer) Depth cueing %s\n",
                   (gl_fog_on) ? "enabled" : "disabled");
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_F12: /* toggle full-screen window on/off */
            gl_fs_on = (!gl_fs_on);
            printf("\nANARIRenderer) Toggling fullscreen window %s\n",
                   (gl_fs_on) ? "on" : "off");
            if (gl_fs_on) { 
              if (glwin_fullscreen(win, gl_fs_on, 0) == 0) {
                owsx = wsx;
                owsy = wsy;
                glwin_get_winsize(win, &wsx, &wsy);
              } else {
                printf("ANARIRenderer) Fullscreen mode note available\n");
              }
            } else {
              glwin_fullscreen(win, gl_fs_on, 0);
              glwin_resize(win, owsx, owsy);
            }
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_UP: /* change depth-of-field focal dist */
            cam_dof_focal_dist *= 1.02f; 
            printf("\nANARIRenderer) DoF focal dist: %f\n", cam_dof_focal_dist);
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_DOWN: /* change depth-of-field focal dist */
            cam_dof_focal_dist *= 0.96f; 
            if (cam_dof_focal_dist < 0.02f) cam_dof_focal_dist = 0.02f;
            printf("\nANARIRenderer) DoF focal dist: %f\n", cam_dof_focal_dist);
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_RIGHT: /* change depth-of-field f/stop number */
            cam_dof_fnumber += 1.0f; 
            printf("\nANARIRenderer) DoF f/stop: %f\n", cam_dof_fnumber);
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_LEFT: /* change depth-of-field f/stop number */
            cam_dof_fnumber -= 1.0f; 
            if (cam_dof_fnumber < 1.0f) cam_dof_fnumber = 1.0f;
            printf("\nANARIRenderer) DoF f/stop: %f\n", cam_dof_fnumber);
            winredraw = 1; 
            break;

          case GLWIN_EV_MOUSE_MOVE:
            if (mousedown != RTMD_NONE) {
              int x, y;
              glwin_get_mousepointer(win, &x, &y);

              float zoommod = 2.0f*cur_cam_zoom/cam_zoom_orig;
              float txdx = (x - mousedownx) * zoommod / wsx;
              float txdy = (y - mousedowny) * zoommod / wsy;
              if (mm != RTMM_SCALE) {
                mousedownx = x;
                mousedowny = y;
              }

              if (mm == RTMM_ROT) {
                Matrix4 rm;
                if (mousedown == RTMD_LEFT) {
                  // when zooming in further from the initial view, we
                  // rotate more slowly so control remains smooth
                  rm.rotate_axis(cam_V, -txdx);
                  rm.rotate_axis(cam_U, -txdy);
                } else if (mousedown == RTMD_MIDDLE || 
                           mousedown == RTMD_RIGHT) {
                  rm.rotate_axis(cam_W, txdx);
                }
                rm.multpoint3d(cam_pos, cam_pos);
                rm.multnorm3d(cam_U, cam_U);
                rm.multnorm3d(cam_V, cam_V);
                rm.multnorm3d(cam_W, cam_W);

                if (xformgradientsphere) {
                  rm.multnorm3d(scene_gradient, scene_gradient);
                }
 
                if (xformlights) {
                  // update light directions (comparatively costly)
                  for (i=0; i<directional_lights.num(); i++) {
                    rm.multnorm3d((float*)&cur_dlights[i].dir, (float*)&cur_dlights[i].dir);
                  }
                }

                winredraw = 1;
              } else if (mm == RTMM_TRANS) {
                if (mousedown == RTMD_LEFT) {
                  float dU[3], dV[3];
                  vec_scale(dU, -txdx, cam_U);
                  vec_scale(dV,  txdy, cam_V);
                  vec_add(cam_pos, cam_pos, dU); 
                  vec_add(cam_pos, cam_pos, dV); 
                } else if (mousedown == RTMD_MIDDLE || 
                           mousedown == RTMD_RIGHT) {
                  float dW[3];
                  vec_scale(dW, txdx, cam_W);
                  vec_add(cam_pos, cam_pos, dW); 
                } 
                winredraw = 1;
              } else if (mm == RTMM_SCALE) {
                float txdx = (x - mousedownx) * 2.0 / wsx;
                float zoominc = 1.0 - txdx;
                if (zoominc < 0.01) zoominc = 0.01;
                cam_zoom = cur_cam_zoom * zoominc;
                winredraw = 1;
              } else if (mm == RTMM_DOF) {
                cam_dof_fnumber += txdx * 20.0f;
                if (cam_dof_fnumber < 1.0f) cam_dof_fnumber = 1.0f;
                cam_dof_focal_dist += -txdy; 
                if (cam_dof_focal_dist < 0.01f) cam_dof_focal_dist = 0.01f;
                winredraw = 1;
              }
            }
            break;

          case GLWIN_EV_MOUSE_LEFT:
          case GLWIN_EV_MOUSE_MIDDLE:
          case GLWIN_EV_MOUSE_RIGHT:
            if (evval) {
              glwin_get_mousepointer(win, &mousedownx, &mousedowny);
              cur_cam_zoom = cam_zoom;

              if (evdev == GLWIN_EV_MOUSE_LEFT) mousedown = RTMD_LEFT;
              else if (evdev == GLWIN_EV_MOUSE_MIDDLE) mousedown = RTMD_MIDDLE;
              else if (evdev == GLWIN_EV_MOUSE_RIGHT) mousedown = RTMD_RIGHT;
            } else {
              mousedown = RTMD_NONE;
            }
            break;

          case GLWIN_EV_MOUSE_WHEELUP:
            cam_zoom /= 1.1f; winredraw = 1; break;

          case GLWIN_EV_MOUSE_WHEELDOWN:
            cam_zoom *= 1.1f; winredraw = 1; break;
        }
      }
    }


    //
    // Support for Spaceball/Spacenavigator/Magellan devices that use
    // X11 ClientMessage protocol....
    //
    //
    // Support for Spaceball/Spacenavigator/Magellan devices that use
    // X11 ClientMessage protocol....
    //
    if (spaceballenabled) {
      // Spaceball/Spacenavigator/Magellan event state variables
      int tx=0, ty=0, tz=0, rx=0, ry=0, rz=0, buttons=0;
      if (glwin_get_spaceball(win, &rx, &ry, &rz, &tx, &ty, &tz, &buttons)) {
        // negate directions if we're in flight mode...
        if (spaceballflightmode) {
          rx= -rx;
          ry= -ry;
          rz= -rz;

          tx= -tx;
          ty= -ty;
          tz= -tz;
        }

        // check for button presses to reset the view
        if (buttons & 1) {
          printf("ANARIRenderer) spaceball button 1 pressed: reset view\n");
          vec_copy(scene_gradient, scene_gradient_orig);
          cam_zoom = cam_zoom_orig;
          vec_copy(cam_pos, cam_pos_orig);
          vec_copy(cam_U, cam_U_orig);
          vec_copy(cam_V, cam_V_orig);
          vec_copy(cam_W, cam_W_orig);

          // restore original light directions
          for (i=0; i<directional_lights.num(); i++) {
            vec_copy((float*)&cur_dlights[i].dir, directional_lights[i].dir);
            vec_normalize((float*)&cur_dlights[i].dir);
          }
          winredraw = 1;
        }

        // check for button presses to toggle spaceball mode
        if (buttons & 2) {
          spaceballmode = !(spaceballmode);
          printf("ANARIRenderer) spaceball mode: %s                       \n",
                 (spaceballmode) ? "scaling" : "rotation/translation");
        }

        // rotation/translation mode
        if (spaceballmode == 0) {
          float zoommod = 2.0f*cam_zoom/cam_zoom_orig;
          float divlen = sqrtf(wsx*wsx + wsy*wsy) * 50;

          // check for rotation and handle it...
          if (rx != 0 || ry !=0 || rz !=0) {
            Matrix4 rm;
            rm.rotate_axis(cam_U, -rx * zoommod / divlen);
            rm.rotate_axis(cam_V, -ry * zoommod / divlen);
            rm.rotate_axis(cam_W, -rz * zoommod / divlen);

            rm.multpoint3d(cam_pos, cam_pos);
            rm.multnorm3d(cam_U, cam_U);
            rm.multnorm3d(cam_V, cam_V);
            rm.multnorm3d(cam_W, cam_W);

            if (xformgradientsphere) {
              rm.multnorm3d(scene_gradient, scene_gradient);
            }

            if (xformlights) {
              // update light directions (comparatively costly)
              for (i=0; i<directional_lights.num(); i++) {
                rm.multnorm3d((float*)&cur_dlights[i].dir, (float*)&cur_dlights[i].dir);
              }
            }
            winredraw = 1;
          }

          // check for translation and handle it...
          if (tx != 0 || ty !=0 || tz !=0) {
            float dU[3], dV[3], dW[3];
            vec_scale(dU, -tx * zoommod / divlen, cam_U);
            vec_scale(dV, -ty * zoommod / divlen, cam_V);
            vec_scale(dW, -tz * zoommod / divlen, cam_W);
            vec_add(cam_pos, cam_pos, dU);
            vec_add(cam_pos, cam_pos, dV);
            vec_add(cam_pos, cam_pos, dW);
            winredraw = 1;
          }
        }
    
        // scaling mode
        if (spaceballmode == 1) {
          const float sbscale = 1.0f / (1024.0f * 8.0f);
          float zoominc = 1.0f - (rz * sbscale);
          if (zoominc < 0.01) zoominc = 0.01;
            cam_zoom *= zoominc;
            winredraw = 1;
        }

      }
    }


    // if there is no HMD, we use the camera orientation directly  
    vec_copy(hmd_U, cam_U);
    vec_copy(hmd_V, cam_V);
    vec_copy(hmd_W, cam_W);

    // XXX HMD handling goes here

    //
    // handle window resizing, stereoscopic mode changes,
    // destroy and recreate affected ANARI buffers
    //
    int resize_buffers=0;

    if (wsx != width) {
      width = wsx;
      resize_buffers=1;
    }
 
    if (wsy != height || (stereoon != stereoon_old)) {
      if (stereoon) {
        if (height != wsy * 2) {
          height = wsy * 2; 
          resize_buffers=1;
        }
      } else {
        height = wsy;
        resize_buffers=1;
      }
    }


    // check if stereo mode or DoF mode changed, both cases
    // require changing the active color accumulation ray gen program
    if ((stereoon != stereoon_old) || (gl_dof_on != gl_dof_on_old)) {
      // when stereo mode changes, we have to regenerate the
      // the RNG, accumulation buffer, and framebuffer
      if (stereoon != stereoon_old) {
        resize_buffers=1;
      }

      // update stereo and DoF state
      stereoon_old = stereoon;
      gl_dof_on_old = gl_dof_on;

      // set the active color accumulation ray gen mode based on the 
      // camera/projection mode, stereoscopic display mode, 
      // and depth-of-field state
      winredraw=1;
    }

    if (resize_buffers) {
      framebuffer_resize(width, height);

      // when movie recording is enabled, print the window size as a guide
      // since the user might want to precisely control the size or 
      // aspect ratio for a particular movie format, e.g. 1080p, 4:3, 16:9
      if (movie_recording_enabled) {
        printf("\rANARIRenderer) Window resize: %d x %d                               \n", width, height);
      }

      winredraw=1;
    }

    int frame_ready = 1; // Default to true
    unsigned int subframe_count = 1;
    if (!done) {
      //
      // If the user interacted with the window in a meaningful way, we
      // need to update the ANARI rendering state, recompile and re-validate
      // the context, and then re-render...
      //
      if (winredraw) {
        // update camera parameters
        anariSetParameter(dev, anariCamera, "position",  ANARI_FLOAT32_VEC3, cam_pos);
        anariSetParameter(dev, anariCamera, "direction", ANARI_FLOAT32_VEC3,   hmd_W);
        anariSetParameter(dev, anariCamera, "up",        ANARI_FLOAT32_VEC3,   hmd_V);
        float camaspect = width / ((float) height);
        anariSetParameter(dev, anariCamera, "aspect", ANARI_FLOAT32, &camaspect);

        float camfovy = 2.0f*180.0f*(atanf(cam_zoom)/M_PI);
        anariSetParameter(dev, anariCamera, "fovy", ANARI_FLOAT32, &camfovy);
 
        // update shadow state 
        // anariSetParameter(dev, anariRenderer, "shadowsEnabled", ANARI_INT32, &gl_shadows_on);

        // update AO state 
        if (gl_shadows_on && gl_ao_on) {
          const int one = 1;
          if (anari_rendermode == ANARI_SCIVIS)
            anariSetParameter(dev, anariRenderer, "aoSamples", ANARI_INT32, &one);
        } else {
          const int zero = 0;
          if (anari_rendermode == ANARI_SCIVIS)
            anariSetParameter(dev, anariRenderer, "aoSamples", ANARI_INT32, &zero);
        }

        // update depth cueing state
        // XXX update ANARI depth cueing state
 
        // update/recompute DoF values 
        // XXX ANARI only implements DoF for the perspective
        //     camera at the present time
        if (camera_projection == ANARIRender::RT_PERSPECTIVE) {
          if (gl_dof_on) {
            anariSetParameter(dev, anariCamera, "focusDistance", ANARI_FLOAT32, &cam_dof_focal_dist);
            float camaprad = cam_dof_focal_dist / (2.0f * cam_zoom * cam_dof_fnumber);
            anariSetParameter(dev, anariCamera, "apertureRadius", ANARI_FLOAT32, &camaprad);
          } else {
            float camaprad = 0.0f;
            anariSetParameter(dev, anariCamera, "apertureRadius", ANARI_FLOAT32, &camaprad);
          }
        }

        // commit camera updates once they're all done...
        anariCommit(dev, anariCamera);

        //
        // Update light directions in the ANARI light buffers
        //
        if (xformlights) {
          // AO scaling factor is applied at the renderer level, but
          // we apply the direct lighting scaling factor to the lights.
          float lightscale = 1.0f;
          if (ao_samples != 0)
            lightscale = ao_direct;

          // XXX assumes the only contents in the first part of the 
          //     light list are directional lights.  The new AO "ambient"
          //     light is the last light in the list now, so we can get
          //     away with this, but refactoring is still needed here.
          for (i=0; i<directional_lights.num(); i++) {
            anariSetParameter(dev, anariLights[i], "intensity", ANARI_FLOAT32, &lightscale);
            anariSetParameter(dev, anariLights[i], "color", ANARI_FLOAT32_VEC3, cur_dlights[i].color);

            float ltmp[3];
            vec_negate(ltmp, cur_dlights[i].dir);
            anariSetParameter(dev, anariLights[i], "direction", ANARI_FLOAT32_VEC3, ltmp);
            anariCommit(dev, anariLights[i]);
          }
        }

        // commit pending changes...
        anariCommit(dev, anariRenderer);

        // reset accumulation buffer 
        accum_count=0;
        totalsamplecount=0;
        if (anariFrameBuffer != NULL) {
//          anariResetAccumulation(dev, anariFrameBuffer);
        }

        // 
        // Sample count updates and ANARI state must always remain in 
        // sync, so if we only update sample count state during redraw events,
        // that's the only time we should recompute the sample counts, since
        // they also affect normalization factors for the accumulation buffer
        //

        // Update sample counts to achieve target interactivity
        if (autosamplecount) {
          if (fpsexpave > 37)
            samples_per_pass++;
          else if (fpsexpave < 30) 
            samples_per_pass--;
    
          // clamp sample counts to a "safe" range
          if (samples_per_pass > 14)
            samples_per_pass=14;
          if (samples_per_pass < 1)
            samples_per_pass=1;
        } 

        // split samples per pass either among AA and AO, depending on
        // whether DoF and AO are enabled or not. 
        if (gl_shadows_on && gl_ao_on) {
          if (gl_dof_on) {
            if (samples_per_pass < 4) {
              cur_aa_samples=samples_per_pass;
              cur_ao_samples=1;
            } else {
              int s = (int) sqrtf(samples_per_pass);
              cur_aa_samples=s;
              cur_ao_samples=s;
            }
          } else {
            cur_aa_samples=1;
            cur_ao_samples=samples_per_pass;
          }
        } else {
          cur_aa_samples=samples_per_pass;
          cur_ao_samples=0;
        }

        // update the current AA/AO sample counts since they may be changing if
        // FPS autotuning is enabled...
        // XXX update ANARI AA sample counts

        // observe latest AO enable/disable flag, and sample count
        if (gl_shadows_on && gl_ao_on) {
          // XXX update ANARI AA/AO sample counts
        } else {
          cur_ao_samples = 0;
          // XXX update ANARI AA/AO sample counts
        }
      } 


      // The accumulation buffer normalization factor must be updated
      // to reflect the total accumulation count before the accumulation
      // buffer is drawn to the output framebuffer
      // XXX update ANARI accum buf normalization factor

      // The accumulation buffer subframe index must be updated to ensure that
      // the RNGs for AA and AO get correctly re-seeded
      // XXX update ANARI accum subframe count

      // Force context compilation/validation
      // render_compile_and_validate();

      anariSetParameter(dev, anariFrameBuffer, "renderer", ANARI_RENDERER, &anariRenderer);
      anariSetParameter(dev, anariFrameBuffer, "camera", ANARI_CAMERA, &anariCamera);
      anariSetParameter(dev, anariFrameBuffer, "world", ANARI_WORLD, &anariWorld);
      anariCommit(dev, anariFrameBuffer);


      //
      // run the renderer 
      //
      frame_ready = 1; // Default to true
      subframe_count = 1;
      if (lasterror == 0 /* XXX SUCCESS */) {
        if (winredraw) {
//          anariResetAccumulation(dev, anariFrameBuffer);
          winredraw=0;
        }

        // iterate, adding to the accumulation buffer...
        anariRenderFrame(dev, anariFrameBuffer);
        anariFrameReady(dev, anariFrameBuffer, ANARI_WAIT);

        subframe_count++; // increment subframe index
        totalsamplecount += samples_per_pass;
        accum_count += cur_aa_samples;

        // copy the accumulation buffer image data to the framebuffer and
        // perform type conversion and normaliztion on the image data...
        // XXX launch ANARI accum copy/norm/finish

        if (lasterror == 0 /* XXX SUCCESS */) {
          if (frame_ready) {
            // display output image
            const unsigned char * img;
            img = (const unsigned char*)anariMapFrame(dev, anariFrameBuffer, "color");

#if 0
            glwin_draw_image_tex_rgb3u(win, (stereoon!=0)*GLWIN_STEREO_OVERUNDER, width, height, img);
#else
            glwin_draw_image_rgb3u(win, (stereoon!=0)*GLWIN_STEREO_OVERUNDER, width, height, img);
#endif
            anariUnmapFrame(dev, anariFrameBuffer, "color");

            // if live movie recording is on, we save every displayed frame
            // to a sequence sequence of image files, with each file numbered
            // by its frame index, which is computed by the multiplying image
            // presentation time by the image sequence fixed-rate-FPS value.
            if (movie_recording_enabled && movie_recording_on) {
              char moviefilename[2048];

              // compute frame number from wall clock time and the
              // current fixed-rate movie playback frame rate
              double now = wkf_timer_timenow(anr_timer);
              double frametime = now - movie_recording_start_time;
              int fidx = frametime * movie_recording_fps;

              // always force the first recorded frame to be 0
              if (movie_framecount==0)
                fidx=0;
              movie_framecount++;

#if defined(__linux)
              // generate symlinks for frame indices between the last written
              // frame and the current one so that video encoders such as
              // ffmpeg and mencoder can be fed the contiguous frame sequence
              // at a fixed frame rate, as they require
              sprintf(moviefilename, movie_recording_filebase,
                      movie_lastframeindex);
              int symidx;
              for (symidx=movie_lastframeindex; symidx<fidx; symidx++) {
                char symlinkfilename[2048];
                sprintf(symlinkfilename, movie_recording_filebase, symidx);
                symlink(moviefilename, symlinkfilename);
              }
#endif

              // write the new movie frame
              sprintf(moviefilename, movie_recording_filebase, fidx);
              const unsigned char *FB = (const unsigned char*)anariMapFrame(dev, anariFrameBuffer, "color");
              if (write_image_file_rgb4u(moviefilename, FB, width, height)) {
                movie_recording_on = 0;
                printf("\n");
                printf("ANARIRenderer) ERROR during writing image during movie recording!\n");
                printf("ANARIRenderer) Movie recording STOPPED\n");
              }
              anariUnmapFrame(dev, anariFrameBuffer, "color");

              movie_lastframeindex = fidx; // update last frame index written
            }
          }
        } else {
          printf("ANARIRenderer) An error occured during rendering. Rendering is aborted.\n");
          done=1;
          break;
        }
      } else {
        printf("ANARIRenderer) An error occured in AS generation. Rendering is aborted.\n");
        done=1;
        break;
      }
    }

    if (!done && frame_ready) {
      double newtime = wkf_timer_timenow(anr_timer);
      double frametime = (newtime-oldtime) + 0.00001f;
      oldtime=newtime;

      // compute exponential moving average for exp(-1/10)
      double framefps = 1.0f/frametime;
      fpsexpave = (fpsexpave * 0.90) + (framefps * 0.10);

      printf("ANARIRenderer) %c AA:%2d AO:%2d, %4d tot RT FPS: %.1f  %.4f s/frame sf: %d  \r",
             statestr[state], cur_aa_samples, cur_ao_samples, 
             totalsamplecount, fpsexpave, frametime, subframe_count);

      fflush(stdout);
      state = (state+1) & 3;
    }

  } // end of per-cycle event processing

  printf("\n");

  // write the output image upon exit...
  if (lasterror == 0 /* XXX SUCCESS */) {
    wkf_timer_start(anr_timer);
    // write output image
    const unsigned char *FB = (const unsigned char*)anariMapFrame(dev, anariFrameBuffer, "color");
    if (write_image_file_rgb4u(filename, FB, width, height)) {
      printf("ANARIRenderer) Failed to write output image!\n");
    }
    anariUnmapFrame(dev, anariFrameBuffer, "color");
    wkf_timer_stop(anr_timer);

    if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
      printf("ANARIRenderer) image file I/O time: %f secs\n", wkf_timer_time(anr_timer));
    }
  }

  glwin_destroy(win);
}

#endif


void ANARIRender::render_to_file(const char *filename) {
  DBG();
  if (!context_created)
    return;

  // Unless overridden by environment variables, we use the incoming
  // window size parameters from VMD to initialize the RT image dimensions.
  int wsx=width, wsy=height;
  const char *imageszstr = getenv("VMDANARIIMAGESIZE");
  if (imageszstr) {
    if (sscanf(imageszstr, "%d %d", &width, &height) != 2) {
      width=wsx;
      height=wsy;
    } 
  } 

  // config/allocate framebuffer and accumulation buffer
  framebuffer_config(width, height);

  update_rendering_state(0);
  render_compile_and_validate();
  double starttime = wkf_timer_timenow(anr_timer);

  // XXX extra world commit for the benefit of USD...
  if (rendererworkarounds == ANARI_USD) {
    printf("ANARIRenderer) *** extra anariWorld commit for USD...\n");
    anariCommit(dev, anariWorld);
  }


  //
  // XXX for ANARI we currently defer FB setup to just before rendering, 
  //     because the frame object serves as our main sync point and
  //     we have to have already committed our renderer, camera, world, etc,
  //     prior to committing the frame for the call to anariRenderFrame().
  // 


  anariFrameBuffer = anariNewFrame(dev);
  int imgsz[2];
  imgsz[0] = width;
  imgsz[1] = height;
  anariSetParameter(dev, anariFrameBuffer, "size", ANARI_UINT32_VEC2, &imgsz);

  // create intermediate output and accumulation buffers
  // One of: UFIXED8_VEC4, UFIXED8_RGBA_SRGB, FLOAT32_VEC4
//  ANARIDataType format = ANARI_UFIXED8_RGBA_SRGB;
  ANARIDataType format = ANARI_UFIXED8_VEC4;
  anariSetParameter(dev, anariFrameBuffer, "color", ANARI_DATA_TYPE, &format);

  anariSetParameter(dev, anariFrameBuffer, "renderer", ANARI_RENDERER, &anariRenderer);
  anariSetParameter(dev, anariFrameBuffer, "camera", ANARI_CAMERA, &anariCamera);
  anariSetParameter(dev, anariFrameBuffer, "world", ANARI_WORLD, &anariWorld);
  anariCommit(dev, anariFrameBuffer);


  //
  // run the renderer 
  //
  if (lasterror == 0 /* XXX SUCCESS */) {
    // clear the accumulation buffer
//    anariResetAccumulation(dev, anariFrameBuffer);

    // Render to the accumulation buffer for the required number of passes
    if (getenv("VMDANARINORENDER") == NULL) {
      if (rendererworkarounds == ANARI_USD) {
        anariRenderFrame(dev, anariFrameBuffer);
        anariFrameReady(dev, anariFrameBuffer, ANARI_WAIT);
      } else {
        int accum_sample;
        for (accum_sample=0; accum_sample<ext_aa_loops; accum_sample++) {
          // The accumulation subframe count must be updated to ensure that
          // any custom RNGs for AA and AO get correctly re-seeded
          anariRenderFrame(dev, anariFrameBuffer);
          anariFrameReady(dev, anariFrameBuffer, ANARI_WAIT);
        }
      }
    }

    // copy the accumulation buffer image data to the framebuffer and perform
    // type conversion and normaliztion on the image data...
    double rtendtime = wkf_timer_timenow(anr_timer);
    time_ray_tracing = rtendtime - starttime;

    if (rendererworkarounds != ANARI_USD) {
      if (lasterror == 0 /* XXX SUCCESS */) {
        // write output image to a file unless we are benchmarking
        if (getenv("VMDANARINOSAVE") == NULL) {
          const unsigned char *FB = (const unsigned char*)anariMapFrame(dev, anariFrameBuffer, "color");
          if (write_image_file_rgb4u(filename, FB, width, height)) {
            printf("ANARIRenderer) Failed to write output image!\n");
          }
          anariUnmapFrame(dev, anariFrameBuffer, "color");
        }
        time_image_io = wkf_timer_timenow(anr_timer) - rtendtime;
      } else {
        printf("ANARIRenderer) Error during rendering.  Rendering aborted.\n");
      }
    }

    if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
      printf("ANARIRenderer) ctx setup %.2f  valid %.2f  AS %.2f  RT %.2f io %.2f\n", time_ctx_setup, time_ctx_validate, time_ctx_AS_build, time_ray_tracing, time_image_io);
    }
  } else {
    printf("ANARIRenderer) Error during AS generation.  Rendering aborted.\n");
  }
}


void ANARIRender::add_material(int matindex,
                               float ambient, float diffuse, 
                               float specular,
                               float shininess, float reflectivity,
                               float opacity, 
                               float outline, float outlinewidth,
                               int transmode) {
if (dev == NULL) {
  printf("add_material() ANARI device is NULL!!!\n");
  return;
}
  

printf("ANARI Add mat[%d]\n", matindex);
  int oldmatcount = materialcache.num();
  if (oldmatcount <= matindex) {
    anr_material m;
    memset(&m, 0, sizeof(m));

    // XXX do something noticable so we see that we got a bad entry...
    m.ambient = 0.5f;
    m.diffuse = 0.7f;
    m.specular = 0.0f;
    m.shininess = 10.0f;
    m.reflectivity = 0.0f;
    m.opacity = 1.0f;
    m.transmode = 0;

    materialcache.appendN(m, matindex - oldmatcount + 1);
  }
 
  if (materialcache[matindex].isvalid) {
    return;
  } else {
    if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) Adding material[%d]\n", matindex);

    materialcache[matindex].ambient      = ambient;
    materialcache[matindex].diffuse      = diffuse; 
    materialcache[matindex].specular     = specular;
    materialcache[matindex].shininess    = shininess;
    materialcache[matindex].reflectivity = reflectivity;
    materialcache[matindex].opacity      = opacity;
    materialcache[matindex].outline      = outline;
    materialcache[matindex].outlinewidth = outlinewidth;
    materialcache[matindex].transmode    = transmode;

    // create an ANARI material object too...
    float mtmp[3];
    ANARIMaterial anariMat;

    if (rendererworkarounds == ANARI_USD) {
      anari_matclass = ANARI_MATTE;
    }


    // choose the right material depending on the active material class
    // and selected renderer type
    switch (anari_matclass) {
      case ANARI_MATTE: 
        {
          if (opacity < 1.0f) {
            // partial cut-out transparency
            anariMat = anariNewMaterial(dev, "transparentMatte");
          } else {
            anariMat = anariNewMaterial(dev, "matte");
          }

          mtmp[0] = mtmp[1] = mtmp[2] = materialcache[matindex].diffuse;
          anariSetParameter(dev, anariMat, "kd", ANARI_FLOAT32_VEC3, mtmp);
          anariSetParameter(dev, anariMat, "d", ANARI_FLOAT32, &materialcache[matindex].opacity);
        }
        break;

      default:
      case ANARI_OBJ: 
        {
          anariMat = anariNewMaterial(dev, "obj");

          mtmp[0] = mtmp[1] = mtmp[2] = materialcache[matindex].diffuse;
          anariSetParameter(dev, anariMat, "kd", ANARI_FLOAT32_VEC3, mtmp);
          anariSetParameter(dev, anariMat, "d", ANARI_FLOAT32, &materialcache[matindex].opacity);

          mtmp[0] = mtmp[1] = mtmp[2] = materialcache[matindex].specular;
          anariSetParameter(dev, anariMat, "ks", ANARI_FLOAT32_VEC3, mtmp);

          anariSetParameter(dev, anariMat, "ns", ANARI_FLOAT32, &materialcache[matindex].shininess);
        }
        break;
    }


    if (rendererworkarounds == ANARI_USD) {
      int usetrue = 1;
      anariSetParameter(dev, anariMat, "usevertexcolors", ANARI_BOOL, &usetrue);

      // if using the USD back-end, assign the "name" tag for the material...
      if (rendererworkarounds == ANARI_USD && (strlen(lastcommentstring) > 0)) {
        char strbuf[2048];
        sprintf(strbuf, "VMD material %d", matindex);
        anariSetParameter(dev, anariMat, "name", ANARI_STRING, strbuf);
      }
    }

    anariCommit(dev, anariMat);
    materialcache[matindex].mat = anariMat;
    materialcache[matindex].isvalid      = 1;
  }
}


// record the most recent comment token for use by ANARI object "name" tags
void ANARIRender::comment(const char *s) {
  commit_rep();  

  printf("ANARIRenderer) comment: '%s'\n", s);
  strncpy(lastcommentstring, s, sizeof(lastcommentstring) - 1);
  lastcommentstring[sizeof(lastcommentstring)-1] = '\0';
}


void ANARIRender::init_materials() {
  DBG();
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) init_materials()\n");

  materialcache.clear();
}


void ANARIRender::set_material(ANARISurface &surf, int matindex, float *uniform_color) {
  if (!context_created)
    return;

  if (verbose == RT_VERB_DEBUG) 
    printf("ANARIRenderer)   setting material %d\n", matindex);
  anariSetParameter(dev, surf, "material", ANARI_MATERIAL, &materialcache[matindex].mat);
}


void ANARIRender::attach_mesh(int numverts, int numfacets, int matindex,
                              anr_trimesh_v3f_n3f_c3f &mesh) {
  mesh.matindex = matindex;
  mesh.verts = anariNewArray1D(dev, mesh.v, 0, 0, ANARI_FLOAT32_VEC3, numverts, 0);
  anariCommit(dev, mesh.verts);
  mesh.norms = anariNewArray1D(dev, mesh.n, 0, 0, ANARI_FLOAT32_VEC3, numverts, 0);
  anariCommit(dev, mesh.norms);
  mesh.cols  = anariNewArray1D(dev, mesh.c, 0, 0, ANARI_FLOAT32_VEC4, numverts, 0);
  anariCommit(dev, mesh.cols);
  mesh.ind   = anariNewArray1D(dev, mesh.f, 0, 0, ANARI_UINT32_VEC3, numfacets, 0);
  anariCommit(dev, mesh.ind);

  mesh.geom  = anariNewGeometry(dev, "triangle");

  anariSetParameter(dev, mesh.geom, "vertex.position", ANARI_ARRAY1D, &mesh.verts);
  anariRelease(dev, mesh.verts);

  anariSetParameter(dev, mesh.geom, "vertex.normal",   ANARI_ARRAY1D, &mesh.norms);
  anariRelease(dev, mesh.norms);

  anariSetParameter(dev, mesh.geom, "primitive.index", ANARI_ARRAY1D, &mesh.ind);
  anariRelease(dev, mesh.ind);

  anariSetParameter(dev, mesh.geom, "vertex.color",    ANARI_ARRAY1D, &mesh.cols);
  anariRelease(dev, mesh.cols);

  anariCommit(dev, mesh.geom);

  mesh.surf = anariNewSurface(dev); 
  anariSetParameter(dev, mesh.surf, "geometry", ANARI_GEOMETRY, &mesh.geom);
  set_material(mesh.surf, matindex, NULL);

  anariCommit(dev, mesh.surf);
  anariRelease(dev, mesh.geom);
  trimesh_v3f_n3f_c3f.append(mesh); 
}


void ANARIRender::attach_sphere_array(int numsp, int matindex,
                                      anr_sphere_array_color &sparray) {
  sparray.matindex = matindex;

  sparray.cents = anariNewArray1D(dev, sparray.xyz, 0, 0, ANARI_FLOAT32_VEC3, numsp, 0);
  anariCommit(dev, sparray.cents);
  sparray.rads = anariNewArray1D(dev, sparray.radii, 0, 0, ANARI_FLOAT32, numsp, 0);
  anariCommit(dev, sparray.rads);
  sparray.cols = anariNewArray1D(dev, sparray.colors, 0, 0, ANARI_FLOAT32_VEC4, numsp, 0);
  anariCommit(dev, sparray.cols);

  sparray.geom  = anariNewGeometry(dev, "sphere");
  anariSetParameter(dev, sparray.geom, "vertex.position", ANARI_ARRAY1D, &sparray.cents);
  anariSetParameter(dev, sparray.geom, "vertex.radius",   ANARI_ARRAY1D, &sparray.rads);
  anariSetParameter(dev, sparray.geom, "vertex.color", ANARI_ARRAY1D, &sparray.cols);
  anariCommit(dev, sparray.geom);
  anariRelease(dev, sparray.cents);
  anariRelease(dev, sparray.rads);
  anariRelease(dev, sparray.cols);

  sparray.surf = anariNewSurface(dev);
  anariSetParameter(dev, sparray.surf, "geometry", ANARI_GEOMETRY, &sparray.geom);
  set_material(sparray.surf, matindex, NULL);
  anariCommit(dev, sparray.surf);
  anariRelease(dev, sparray.geom);

  spheres_color.append(sparray);
}


void ANARIRender::attach_cylinder_array(int numcyl, int matindex,
                                        anr_cylinder_array_color &cylarray) {
  cylarray.matindex = matindex;
  cylarray.cyls = anariNewArray1D(dev, cylarray.verts, 0, 0, ANARI_FLOAT32_VEC3, numcyl * 2, 0);
  anariCommit(dev, cylarray.cyls);
  cylarray.rads = anariNewArray1D(dev, cylarray.radii, 0, 0, ANARI_FLOAT32, numcyl, 0);
  anariCommit(dev, cylarray.rads);
  cylarray.cols = anariNewArray1D(dev, cylarray.colors, 0, 0, ANARI_FLOAT32_VEC4, numcyl, 0);
  anariCommit(dev, cylarray.cols);
  cylarray.ind  = anariNewArray1D(dev, cylarray.indices, 0, 0, ANARI_UINT32, numcyl * 1, 0);
  anariCommit(dev, cylarray.ind);

  cylarray.geom  = anariNewGeometry(dev, "curve");
  anariSetParameter(dev, cylarray.geom, "vertex.position", ANARI_ARRAY1D, &cylarray.cyls);
  anariSetParameter(dev, cylarray.geom, "vertex.radius",   ANARI_ARRAY1D, &cylarray.rads);
  anariSetParameter(dev, cylarray.geom, "vertex.index",    ANARI_ARRAY1D, &cylarray.ind);

  cylarray.surf = anariNewSurface(dev);
  anariSetParameter(dev, cylarray.surf, "geometry", ANARI_GEOMETRY, &cylarray.geom);
  anariSetParameter(dev, cylarray.surf, "primitive.color", ANARI_ARRAY1D, &cylarray.cols);
  set_material(cylarray.surf, matindex, NULL);
  anariCommit(dev, cylarray.surf);
  anariRelease(dev, cylarray.geom);

  free(cylarray.cyls);
  cylarray.cyls = NULL;
  free(cylarray.rads);
  cylarray.rads = NULL;
  free(cylarray.ind);
  cylarray.ind = NULL;
  free(cylarray.cols);
  cylarray.cols = NULL;

  cylinders_color.append(cylarray);
}



void ANARIRender::add_directional_light(const float *dir, const float *color) {
  DBG();
  anr_directional_light l;
  vec_copy(l.dir, dir);
  vec_copy(l.color, color);

  directional_lights.append(l);
}


void ANARIRender::add_positional_light(const float *pos, const float *color) {
  DBG();
  anr_positional_light l;
  vec_copy(l.pos, pos);
  vec_copy(l.color, color);

  positional_lights.append(l);
}


void ANARIRender::cylinder_array(Matrix4 *wtrans, float radius,
                                 float *uniform_color,
                                 int cylnum, float *points, int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating cylinder array: %d...\n", cylnum);

  cylinder_array_cnt += cylnum;
  
  anr_cylinder_array_color ca;
  memset(&ca, 0, sizeof(ca));
  ca.num = cylnum;
  ca.verts = (float *) calloc(1, cylnum * 6 * sizeof(float));
  ca.radii = (float *) calloc(1, cylnum * 1 * sizeof(float));
  ca.colors = (float *) calloc(1, cylnum * 4 * sizeof(float));
  ca.indices = (unsigned int *) calloc(1, cylnum * 2 * sizeof(unsigned int));

  int i,ind4,ind6;
  if (wtrans == NULL) {
    for (i=0,ind4=0,ind6=0; i<cylnum; i++,ind4+=4,ind6+=6) {
      vec_copy(&ca.verts[ind6  ], &points[ind6  ]);
      vec_copy(&ca.verts[ind6+3], &points[ind6+3]);
      ca.radii[i] = radius;
      vec_copy(&ca.colors[ind4], &uniform_color[0]);
      ca.colors[ind4 + 3] = 1.0f;
      ca.indices[i] = i*2;
    }
  } else {
    for (i=0,ind4=0,ind6=0; i<cylnum; i++,ind4+=4,ind6+=6) {
      // apply transforms on points, radii
      wtrans->multpoint3d(&points[ind6  ], &ca.verts[ind6  ]);
      wtrans->multpoint3d(&points[ind6+3], &ca.verts[ind6+3]);
      ca.radii[i] = radius;
      vec_copy(&ca.colors[ind4], &uniform_color[0]);
      ca.colors[ind4 + 3] = 1.0f;
      ca.indices[i] = i*2;
    }
  }

  attach_cylinder_array(cylnum, matindex, ca);
}


void ANARIRender::cylinder_array_color(Matrix4 & wtrans, float rscale,
                                       int cylnum, float *points,
                                       float *radii, float *colors,
                                       int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating color cylinder array: %d...\n", cylnum);
  cylinder_array_color_cnt += cylnum;

  anr_cylinder_array_color cac;
  memset(&cac, 0, sizeof(cac));
  cac.num = cylnum;
  cac.verts = (float *) calloc(1, cylnum * 6 * sizeof(float));
  cac.radii = (float *) calloc(1, cylnum * 1 * sizeof(float));
  cac.colors = (float *) calloc(1, cylnum * 4 * sizeof(float));
  cac.indices = (unsigned int *) calloc(1, cylnum * 2 * sizeof(unsigned int));

  int i, ind3, ind4, ind6;
  for (i=0,ind3=0,ind4=0,ind6=0; i<cylnum; i++,ind3+=3,ind4+=4,ind6+=6) {
    // apply transforms on points, radii
    wtrans.multpoint3d(&points[ind6  ], &cac.verts[ind6  ]);
    wtrans.multpoint3d(&points[ind6+3], &cac.verts[ind6+3]);
    cac.radii[i] = radii[i] * rscale; // radius
    vec_copy(&cac.colors[ind4], &colors[ind3]);
    cac.colors[ind4 + 3] = 1.0f;
    cac.indices[i] = i*2;
  }

  attach_cylinder_array(cylnum, matindex, cac);
}

#if 0
void ANARIRender::ring_array_color(Matrix4 & wtrans, float rscale,
                                   int rnum, float *centers,
                                   float *norms, float *radii, 
                                   float *colors, int matindex) {
}
#endif


void ANARIRender::sphere_array(Matrix4 *wtrans, float rscale,
                               float *uniform_color,
                               int numsp, float *centers,
                               float *radii,
                               int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating sphere array: %d...\n", numsp);
  sphere_array_cnt += numsp;

  const rgba c = { uniform_color[0], uniform_color[1], uniform_color[2], 1.0f};

  anr_sphere_array_color sp;
  memset(&sp, 0, sizeof(sp));
  sp.num = numsp;
  sp.xyz = (float *) calloc(1, numsp * 3*sizeof(float));
  sp.radii = (float *) calloc(1, numsp * sizeof(float));
  sp.colors = (float *) calloc(1, numsp * 4*sizeof(float));

  int i, ind3, ind4;
  if (wtrans == NULL) {
    if (radii == NULL) {
      for (i=0,ind3=0,ind4=0; i<numsp; i++,ind3+=3,ind4+=4) {
        // transform to eye coordinates
        vec_copy((float*) &sp.xyz[ind3], &centers[ind3]);
        sp.radii[i] = rscale;
        memcpy((float*) &sp.colors[ind4], &c, 4*sizeof(float));
      }
    } else {
      for (i=0,ind3=0,ind4=0; i<numsp; i++,ind3+=3,ind4+=4) {
        // transform to eye coordinates
        vec_copy((float*) &sp.xyz[ind3], &centers[ind3]);
        sp.radii[i] = radii[i] * rscale;
        memcpy((float*) &sp.colors[ind4], &c, 4*sizeof(float));
      }
    }
  } else {
    if (radii == NULL) {
      for (i=0,ind3=0,ind4=0; i<numsp; i++,ind3+=3,ind4+=4) {
        wtrans->multpoint3d(&centers[ind3], &sp.xyz[ind3]);
        sp.radii[i] = rscale;
        memcpy((float*) &sp.colors[ind4], &c, 4*sizeof(float));
      }
    } else {
      for (i=0,ind3=0,ind4=0; i<numsp; i++,ind3+=3,ind4+=4) {
        // transform to eye coordinates
        wtrans->multpoint3d(&centers[ind3], &sp.xyz[ind3]);
        sp.radii[i] = radii[i] * rscale;
        memcpy((float*) &sp.colors[ind4], &c, 4*sizeof(float));
      }
    }
  }

  attach_sphere_array(numsp, matindex, sp);
}


void ANARIRender::sphere_array_color(Matrix4 & wtrans, float rscale,
                                     int numsp, float *centers,
                                     float *radii, float *colors,
                                     int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating sphere array color: %d...\n", numsp);
  sphere_array_color_cnt += numsp;

  anr_sphere_array_color sp;
  memset(&sp, 0, sizeof(sp));
  sp.num = numsp;
  sp.xyz = (float *) calloc(1, numsp * 3*sizeof(float));
  sp.radii = (float *) calloc(1, numsp * sizeof(float));
  sp.colors = (float *) calloc(1, numsp * 4*sizeof(float));

  int i, ind3, ind4;
  for (i=0,ind3=0,ind4=0; i<numsp; i++,ind3+=3,ind4+=4) {
    wtrans.multpoint3d(&centers[ind3], &sp.xyz[ind3]);
    sp.radii[i] = radii[i] * rscale;
    vec_copy((float*) &sp.colors[ind4], &colors[ind3]);
    sp.colors[ind4 + 3] = 1.0f;
  }

  attach_sphere_array(numsp, matindex, sp);
}


void ANARIRender::tricolor_list(Matrix4 & wtrans, int numtris, float *vnc,
                                int matindex) {
  if (!context_created) return;
//if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating tricolor list: %d...\n", numtris);
  tricolor_cnt += numtris;

  // create and fill the ANARI trimesh memory buffer
  anr_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numtris;
  mesh.v = (float *) calloc(1, numtris * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numtris * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numtris * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numtris * 3*sizeof(int));
  
  float alpha = 1.0f;
//  alpha = materialcache[matindex].opacity;

  int i, ind, ind9, ind12;
  for (i=0,ind=0,ind9=0,ind12=0; i<numtris; i++,ind+=27,ind9+=9,ind12+=12) {
    // transform to eye coordinates
    wtrans.multpoint3d(&vnc[ind    ], (float*) &mesh.v[ind9    ]);
    wtrans.multpoint3d(&vnc[ind + 3], (float*) &mesh.v[ind9 + 3]);
    wtrans.multpoint3d(&vnc[ind + 6], (float*) &mesh.v[ind9 + 6]);

    wtrans.multnorm3d(&vnc[ind +  9], (float*) &mesh.n[ind9    ]);
    wtrans.multnorm3d(&vnc[ind + 12], (float*) &mesh.n[ind9 + 3]);
    wtrans.multnorm3d(&vnc[ind + 15], (float*) &mesh.n[ind9 + 6]);

    vec_copy(&mesh.c[ind12    ], &vnc[ind + 18]);
    mesh.c[ind12 +  3] = alpha;
    vec_copy(&mesh.c[ind12 + 4], &vnc[ind + 21]);
    mesh.c[ind12 +  7] = alpha;
    vec_copy(&mesh.c[ind12 + 8], &vnc[ind + 24]);
    mesh.c[ind12 + 11] = alpha;

    mesh.f[i*3  ] = i*3;
    mesh.f[i*3+1] = i*3 + 1;
    mesh.f[i*3+2] = i*3 + 2;
  }

  attach_mesh(numtris * 3, numtris, matindex, mesh);
}


void ANARIRender::trimesh_c4n3v3(Matrix4 & wtrans, int numverts,
                                 float *cnv, int numfacets, int * facets,
                                 int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating trimesh_c4n3v3: %d...\n", numfacets);
  trimesh_c4u_n3b_v3f_cnt += numfacets;

  // create and fill the ANARI trimesh memory buffer
  anr_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numfacets;
  mesh.v = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numfacets * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numfacets * 3*sizeof(int));
 
  float alpha = 1.0f;
//  alpha = materialcache[matindex].opacity;

  // XXX we are currently converting to triangle soup for ease of
  // initial implementation, but this is clearly undesirable long-term
  int i, ind, ind9, ind12;
  for (i=0,ind=0,ind9=0,ind12=0; i<numfacets; i++,ind+=3,ind9+=9,ind12+=12) {
    int v0 = facets[ind    ] * 10;
    int v1 = facets[ind + 1] * 10;
    int v2 = facets[ind + 2] * 10;

    // transform to eye coordinates
    wtrans.multpoint3d(cnv + v0 + 7, (float*) &mesh.v[ind9    ]);
    wtrans.multpoint3d(cnv + v1 + 7, (float*) &mesh.v[ind9 + 3]);
    wtrans.multpoint3d(cnv + v2 + 7, (float*) &mesh.v[ind9 + 6]);

    wtrans.multnorm3d(cnv + v0 + 4, (float*) &mesh.n[ind9    ]);
    wtrans.multnorm3d(cnv + v1 + 4, (float*) &mesh.n[ind9 + 3]);
    wtrans.multnorm3d(cnv + v2 + 4, (float*) &mesh.n[ind9 + 6]);

    vec_copy(&mesh.c[ind12    ], cnv + v0);
    mesh.c[ind12 +  3] = alpha;
    vec_copy(&mesh.c[ind12 + 4], cnv + v1);
    mesh.c[ind12 +  7] = alpha;
    vec_copy(&mesh.c[ind12 + 8], cnv + v2);
    mesh.c[ind12 + 11] = alpha;

    mesh.f[i*3  ] = i*3;
    mesh.f[i*3+1] = i*3 + 1;
    mesh.f[i*3+2] = i*3 + 2;
  }

  attach_mesh(numfacets * 3, numfacets, matindex, mesh);
}



// 
// This implementation translates from the most-compact host representation
// to the best that ANARI allows
//
void ANARIRender::trimesh_c4u_n3b_v3f(Matrix4 & wtrans, unsigned char *c, 
                                      signed char *n, float *v, 
                                      int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating trimesh_c4u_n3b_v3f: %d...\n", numfacets);
  trimesh_n3b_v3f_cnt += numfacets;

  // create and fill the ANARI trimesh memory buffer
  anr_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numfacets;
  mesh.v = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numfacets * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numfacets * 3*sizeof(int));
 
  float alpha = 1.0f;
//  alpha = materialcache[matindex].opacity;

  // XXX we are currently converting to triangle soup for ease of
  // initial implementation, but this is clearly undesirable long-term
  int i, ind, ind9, ind12;

  const float ci2f = 1.0f / 255.0f;
  const float cn2f = 1.0f / 127.5f;
  for (i=0,ind=0,ind9=0,ind12=0; i<numfacets; i++,ind+=3,ind9+=9,ind12+=12) {
    float norm[9];

    // conversion from GLbyte format, Table 2.6, p. 44 of OpenGL spec 1.2.1
    // float = (2c+1)/(2^8-1)
    norm[0] = n[ind9    ] * cn2f + ci2f;
    norm[1] = n[ind9 + 1] * cn2f + ci2f;
    norm[2] = n[ind9 + 2] * cn2f + ci2f;
    norm[3] = n[ind9 + 3] * cn2f + ci2f;
    norm[4] = n[ind9 + 4] * cn2f + ci2f;
    norm[5] = n[ind9 + 5] * cn2f + ci2f;
    norm[6] = n[ind9 + 6] * cn2f + ci2f;
    norm[7] = n[ind9 + 7] * cn2f + ci2f;
    norm[8] = n[ind9 + 8] * cn2f + ci2f;

    // transform to eye coordinates
    wtrans.multpoint3d(v + ind9    , (float*) &mesh.v[ind9    ]);
    wtrans.multpoint3d(v + ind9 + 3, (float*) &mesh.v[ind9 + 3]);
    wtrans.multpoint3d(v + ind9 + 6, (float*) &mesh.v[ind9 + 6]);

    wtrans.multnorm3d(norm    , (float*) &mesh.n[ind9    ]);
    wtrans.multnorm3d(norm + 3, (float*) &mesh.n[ind9 + 3]);
    wtrans.multnorm3d(norm + 6, (float*) &mesh.n[ind9 + 6]);

    float col[9];

    // conversion from GLubyte format, Table 2.6, p. 44 of OpenGL spec 1.2.1
    // float = c/(2^8-1)
    col[0] = c[ind12     ] * ci2f;
    col[1] = c[ind12 +  1] * ci2f;
    col[2] = c[ind12 +  2] * ci2f;
    col[3] = c[ind12 +  4] * ci2f;
    col[4] = c[ind12 +  5] * ci2f;
    col[5] = c[ind12 +  6] * ci2f;
    col[6] = c[ind12 +  8] * ci2f;
    col[7] = c[ind12 +  9] * ci2f;
    col[8] = c[ind12 + 10] * ci2f;

    vec_copy(&mesh.c[ind12    ], col    );
    mesh.c[ind12 +  3] = alpha;
    vec_copy(&mesh.c[ind12 + 4], col + 3);
    mesh.c[ind12 +  7] = alpha;
    vec_copy(&mesh.c[ind12 + 8], col + 6);
    mesh.c[ind12 + 11] = alpha;

    mesh.f[i*3  ] = i*3;
    mesh.f[i*3+1] = i*3 + 1;
    mesh.f[i*3+2] = i*3 + 2;
  }

  attach_mesh(numfacets * 3, numfacets, matindex, mesh);
}



void ANARIRender::trimesh_c4u_n3f_v3f(Matrix4 & wtrans, unsigned char *c, 
                                      float *n, float *v, 
                                      int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating trimesh_c4u_n3f_v3f: %d...\n", numfacets);
  tricolor_cnt += numfacets;

  // create and fill the ANARI trimesh memory buffer
  anr_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numfacets;
  mesh.v = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numfacets * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numfacets * 3*sizeof(int));
 
  float alpha = 1.0f;
//  alpha = materialcache[matindex].opacity;

  // XXX we are currently converting to triangle soup for ease of
  // initial implementation, but this is clearly undesirable long-term
  int i, ind, ind9, ind12;

  const float ci2f = 1.0f / 255.0f;
  for (i=0,ind=0,ind9=0,ind12=0; i<numfacets; i++,ind+=3,ind9+=9,ind12+=12) {
    // transform to eye coordinates
    wtrans.multpoint3d(v + ind9    , (float*) &mesh.v[ind9    ]);
    wtrans.multpoint3d(v + ind9 + 3, (float*) &mesh.v[ind9 + 3]);
    wtrans.multpoint3d(v + ind9 + 6, (float*) &mesh.v[ind9 + 6]);

    wtrans.multnorm3d(n + ind9    , &mesh.n[ind9    ]);
    wtrans.multnorm3d(n + ind9 + 3, &mesh.n[ind9 + 3]);
    wtrans.multnorm3d(n + ind9 + 6, &mesh.n[ind9 + 3]);

    // conversion from GLubyte format, Table 2.6, p. 44 of OpenGL spec 1.2.1
    // float = c/(2^8-1)
    float col[9];
    col[0] = c[ind12     ] * ci2f;
    col[1] = c[ind12 +  1] * ci2f;
    col[2] = c[ind12 +  2] * ci2f;
    col[3] = c[ind12 +  4] * ci2f;
    col[4] = c[ind12 +  5] * ci2f;
    col[5] = c[ind12 +  6] * ci2f;
    col[6] = c[ind12 +  8] * ci2f;
    col[7] = c[ind12 +  9] * ci2f;
    col[8] = c[ind12 + 10] * ci2f;

    vec_copy(&mesh.c[ind12    ], col    );
    mesh.c[ind12 +  3] = alpha;
    vec_copy(&mesh.c[ind12 + 4], col + 3);
    mesh.c[ind12 +  7] = alpha;
    vec_copy(&mesh.c[ind12 + 8], col + 6);
    mesh.c[ind12 + 11] = alpha;

    mesh.f[i*3  ] = i*3;
    mesh.f[i*3+1] = i*3 + 1;
    mesh.f[i*3+2] = i*3 + 2;
  }

  attach_mesh(numfacets * 3, numfacets, matindex, mesh);
}


void ANARIRender::trimesh_n3b_v3f(Matrix4 & wtrans, float *uniform_color, 
                                  signed char *n, float *v, 
                                  int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating trimesh_n3b_v3f: %d...\n", numfacets);
  trimesh_n3b_v3f_cnt += numfacets;

  // create and fill the ANARI trimesh memory buffer
  anr_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numfacets;
  mesh.v = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numfacets * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numfacets * 3*sizeof(int));
 
  float alpha = 1.0f;

  // XXX we are currently converting to triangle soup for ease of
  // initial implementation, but this is clearly undesirable long-term
  int i, ind, ind9, ind12;

  const float ci2f = 1.0f / 255.0f;
  const float cn2f = 1.0f / 127.5f;
  for (i=0,ind=0,ind9=0,ind12=0; i<numfacets; i++,ind+=3,ind9+=9,ind12+=12) {
    float norm[9];

    // conversion from GLbyte format, Table 2.6, p. 44 of OpenGL spec 1.2.1
    // float = (2c+1)/(2^8-1)
    norm[0] = n[ind9    ] * cn2f + ci2f;
    norm[1] = n[ind9 + 1] * cn2f + ci2f;
    norm[2] = n[ind9 + 2] * cn2f + ci2f;
    norm[3] = n[ind9 + 3] * cn2f + ci2f;
    norm[4] = n[ind9 + 4] * cn2f + ci2f;
    norm[5] = n[ind9 + 5] * cn2f + ci2f;
    norm[6] = n[ind9 + 6] * cn2f + ci2f;
    norm[7] = n[ind9 + 7] * cn2f + ci2f;
    norm[8] = n[ind9 + 8] * cn2f + ci2f;

    // transform to eye coordinates
    wtrans.multpoint3d(v + ind9    , (float*) &mesh.v[ind9    ]);
    wtrans.multpoint3d(v + ind9 + 3, (float*) &mesh.v[ind9 + 3]);
    wtrans.multpoint3d(v + ind9 + 6, (float*) &mesh.v[ind9 + 6]);

    wtrans.multnorm3d(norm    , (float*) &mesh.n[ind9    ]);
    wtrans.multnorm3d(norm + 3, (float*) &mesh.n[ind9 + 3]);
    wtrans.multnorm3d(norm + 6, (float*) &mesh.n[ind9 + 6]);

    vec_copy(&mesh.c[ind12    ], uniform_color);
    mesh.c[ind12 +  3] = alpha;
    vec_copy(&mesh.c[ind12 + 4], uniform_color);
    mesh.c[ind12 +  7] = alpha;
    vec_copy(&mesh.c[ind12 + 8], uniform_color);
    mesh.c[ind12 + 11] = alpha;

    mesh.f[i*3  ] = i*3;
    mesh.f[i*3+1] = i*3 + 1;
    mesh.f[i*3+2] = i*3 + 2;
  }

  attach_mesh(numfacets * 3, numfacets, matindex, mesh);
}


// XXX At present we have to build/populate a per-vertex color arrays,
//     but that should go away as soon as ANARI allows it.
void ANARIRender::trimesh_n3f_v3f(Matrix4 & wtrans, float *uniform_color, 
                                  float *n, float *v, int numfacets, 
                                  int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating trimesh_n3f_v3f: %d...\n", numfacets);
  trimesh_n3f_v3f_cnt += numfacets;
  // create and fill the ANARI trimesh memory buffer
  anr_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numfacets;
  mesh.v = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numfacets * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numfacets * 3*sizeof(int));

  float alpha = 1.0f;

  // create and fill the ANARI trimesh memory buffer
  int i, ind, ind9, ind12;

  for (i=0,ind=0,ind9=0,ind12=0; i<numfacets; i++,ind+=3,ind9+=9,ind12+=12) {
    // transform to eye coordinates
    wtrans.multpoint3d(v + ind9    , (float*) &mesh.v[ind9    ]);
    wtrans.multpoint3d(v + ind9 + 3, (float*) &mesh.v[ind9 + 3]);
    wtrans.multpoint3d(v + ind9 + 6, (float*) &mesh.v[ind9 + 6]);

    wtrans.multnorm3d(n + ind9    , (float*) &mesh.n[ind9    ]);
    wtrans.multnorm3d(n + ind9 + 3, (float*) &mesh.n[ind9 + 3]);
    wtrans.multnorm3d(n + ind9 + 6, (float*) &mesh.n[ind9 + 6]);

    vec_copy(&mesh.c[ind12    ], uniform_color);
    mesh.c[ind12 +  3] = alpha;
    vec_copy(&mesh.c[ind12 + 4], uniform_color);
    mesh.c[ind12 +  7] = alpha;
    vec_copy(&mesh.c[ind12 + 8], uniform_color);
    mesh.c[ind12 + 11] = alpha;

    mesh.f[i*3  ] = i*3;
    mesh.f[i*3+1] = i*3 + 1;
    mesh.f[i*3+2] = i*3 + 2;
  }

  attach_mesh(numfacets * 3, numfacets, matindex, mesh);
}


#if 0
void ANARIRender::trimesh_v3f(Matrix4 & wtrans, float *uniform_color, 
                              float *v, int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating trimesh_v3f: %d...\n", numfacets);
  trimesh_v3f_cnt += numfacets;

  set_material(geom, matindex, NULL);
  append_objects(buf, geom, instance);
}

#endif



void ANARIRender::tristrip(Matrix4 & wtrans, int numverts, const float * cnv,
                           int numstrips, const int *vertsperstrip,
                           const int *facets, int matindex) {
  if (!context_created) return;
  int i;
  int numfacets = 0;
  for (i=0; i<numstrips; i++) 
    numfacets += (vertsperstrip[i] - 2);  

  if (verbose == RT_VERB_DEBUG) printf("ANARIRenderer) creating tristrip: %d...\n", numfacets);
  tricolor_cnt += numfacets;

  // create and fill the ANARI trimesh memory buffer
  anr_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numfacets;
  mesh.v = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numfacets * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numfacets * 3*sizeof(int));

  float alpha = 1.0f;
//  alpha = materialcache[matindex].opacity;

  // XXX we are currently converting to triangle soup for ease of
  // initial implementation, but this is clearly undesirable long-term

  // render triangle strips one triangle at a time
  // triangle winding order is:
  //   v0, v1, v2, then v2, v1, v3, then v2, v3, v4, etc.
  int strip, t, v = 0;
  int stripaddr[2][3] = { {0, 1, 2}, {1, 0, 2} };

  // loop over all of the triangle strips
  i=0; // set triangle index to 0
  int ind9, ind12;
  for (strip=0,ind9=0,ind12=0; strip < numstrips; strip++) {
    // loop over all triangles in this triangle strip
    for (t = 0; t < (vertsperstrip[strip] - 2); t++) {
      // render one triangle, using lookup table to fix winding order
      int v0 = facets[v + (stripaddr[t & 0x01][0])] * 10;
      int v1 = facets[v + (stripaddr[t & 0x01][1])] * 10;
      int v2 = facets[v + (stripaddr[t & 0x01][2])] * 10;

      // transform to eye coordinates
      wtrans.multpoint3d(cnv + v0 + 7, (float*) &mesh.v[ind9    ]);
      wtrans.multpoint3d(cnv + v1 + 7, (float*) &mesh.v[ind9 + 3]);
      wtrans.multpoint3d(cnv + v2 + 7, (float*) &mesh.v[ind9 + 6]);

      wtrans.multnorm3d(cnv + v0 + 4, (float*) &mesh.n[ind9    ]);
      wtrans.multnorm3d(cnv + v1 + 4, (float*) &mesh.n[ind9 + 3]);
      wtrans.multnorm3d(cnv + v2 + 4, (float*) &mesh.n[ind9 + 6]);

      vec_copy(&mesh.c[ind12    ], cnv + v0);
      mesh.c[ind12 +  3] = alpha;
      vec_copy(&mesh.c[ind12 + 4], cnv + v1);
      mesh.c[ind12 +  7] = alpha;
      vec_copy(&mesh.c[ind12 + 8], cnv + v2);
      mesh.c[ind12 + 11] = alpha;

      mesh.f[i*3  ] = i*3;
      mesh.f[i*3+1] = i*3 + 1;
      mesh.f[i*3+2] = i*3 + 2;

      v++; // move on to next vertex
      i++; // next triangle
      ind9+=9;
      ind12+=12;
    }
    v+=2; // last two vertices are already used by last triangle
  }

  attach_mesh(numfacets * 3, numfacets, matindex, mesh);
}



