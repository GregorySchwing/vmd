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
*      $RCSfile: OSPRay2Renderer.C,v $
*      $Author: johns $      $Locker:  $               $State: Exp $
*      $Revision: 1.35 $         $Date: 2021/12/21 22:53:29 $
*
***************************************************************************/
/**
 *  \file OSPRay2Renderer.C
 *  \brief VMD built-in Tachyon/OSPRay ray tracing engine.
 * 
 *  This work is briefly outlined in:
 *   "OSPRay - A CPU Ray Tracing Framework for Scientific Visualization"
 *    Ingo Wald, Gregory Johnson, Jefferson Amstutz, Carson Brownlee,
 *    Aaron Knoll, Jim Jeffers, Johannes Guenther, Paul Navratil
 *    IEEE Vis, 2016 (in-press)
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(__linux)
#include <unistd.h>   // needed for symlink() in movie recorder
#endif

#include "Inform.h"
#include "ImageIO.h"
#include "OSPRay2Renderer.h"
// #include "OSPRay2Shaders.ih" /// ISPC code at some point?
#include "Matrix4.h"
#include "utilities.h"
#include "WKFUtils.h"

// #if !(OSPRAY_VERSION_MAJOR >= 1 && OSPRAY_VERSION_MINOR >= 2)
// #error VMD requires OSPRay >= 1.2.0 for correct transparent AO shading
// // VMD requires OSPRay >= 1.1.2 for correct rendering of cylinders
// #endif

// enable the interactive ray tracing capability
#if defined(VMDOSPRAY_INTERACTIVE_OPENGL)
#if (defined(WIN32) || defined(_WIN64)) && defined(_MSC_VER)
#include <windows.h> // must include windows.h prior to GL
#endif

#include <GL/gl.h>
#endif

#if 0
#define DBG() 
#else
#define DBG() printf("OSPRay2Renderer) %s\n", __func__);
#endif

static void vmd_ospray2_error_callback(void *userData, OSPError err, const char *detailstr) {
  printf("OSPRay2Renderer) ERROR: %s\n", detailstr);
}

static void vmd_ospray2_status_callback(void *userData, const char *detailstr) {
  printf("OSPRay2Renderer) STATUS: %s", detailstr);
}


// Global OSPRay initialization routine -- call it only ONCE...
int OSPRay2Renderer::OSPRay_Global_Init(void) {
  DBG();

  // initialize OSPRay itself
  const char *ospraynormalargs[] = {"vmd", "--osp:mpi"};
  const char *ospraydebugargs[] = {"vmd", "--osp:debug", "--osp:mpi"};
  const char **osprayargs = ospraynormalargs; 
  int argcount = 1;

  if (getenv("VMDOSPRAYDEBUG") != NULL) {
    osprayargs = ospraydebugargs;
    argcount=2;    
  }

  // only pass in the second "--osp:mpi" flag if the user has
  // requested that the MPI renderer back-end be enabled through
  // environment variable flags
  if (getenv("VMDOSPRAYMPI") || getenv("VMD_OSPRAY_MPI")) {
    msgInfo << "OSPRay2Renderer) Initializing OSPRay in MPI mode" << sendmsg;
    argcount++;
  }
 
  OSPError osprc = ospInit(&argcount, osprayargs);

  if (osprc != OSP_NO_ERROR)
    return -1; // return indication of failure

  return 0;
}

// Global OSPRay initialization routine -- call it only ONCE...
void OSPRay2Renderer::OSPRay_Global_Shutdown(void) {
  DBG();
  ospShutdown();
}

/// constructor ... initialize some variables
OSPRay2Renderer::OSPRay2Renderer(void) {
  DBG();

#if 1
printf("debugging PID: %d\n", getpid());
if (getenv("VMDOSPRAYSLEEP")) {
  int sleepsecs = atoi(getenv("VMDOSPRAYSLEEP"));
  sleep(sleepsecs);
}
#endif

  osp_timer = wkf_timer_create(); // create and initialize timer
  wkf_timer_start(osp_timer);

  osp_rendermode = RT_PATHTRACER;
  if (getenv("VMDOSPRAYSCIVIS")) {
    printf("OSPRay2Renderer) Renderer mode set to 'scivis'\n");
    osp_rendermode = RT_SCIVIS;
  }
  if (getenv("VMDOSPRAYPATHTRACER")) {
    printf("OSPRay2Renderer) Renderer mode set to 'pathtracer'\n");
    osp_rendermode = RT_PATHTRACER;
  }

  // set OSPRay state handles/variables to NULL
  ospRenderer = NULL;
  ospFrameBuffer = NULL;
  ospCamera = NULL;
  ospWorld = NULL;
  ospLightData = NULL;

  lasterror = 0;               // begin with no error state set
  context_created = 0;         // no context yet
  buffers_allocated = 0;       // flag no buffer allocated yet
  scene_created = 0;           // scene has been created

  destroy_scene();             // clear/init geometry vectors

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

  verbose = RT_VERB_MIN;  // keep console quiet except for perf/debugging cases
  check_verbose_env();    // see if the user has overridden verbose flag

  ospInstances.clear();   // clear instance list

  // clear all primitive lists
  trimesh_v3f_n3f_c3f.clear();
  spheres_color.clear();
  cylinders_color.clear();

  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating context...\n");

  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) setting error / status callbacks...\n");
  OSPDevice dev = ospGetCurrentDevice();
  ospDeviceSetErrorCallback(dev, vmd_ospray2_error_callback, NULL);
  ospDeviceSetStatusCallback(dev, vmd_ospray2_status_callback, NULL);
  int loglevel = OSP_LOG_INFO;
//  loglevel = OSP_LOG_DEBUG;
  ospDeviceSetParam(dev, "logLevel", OSP_INT, &loglevel);
  ospDeviceCommit(dev);
  ospDeviceRelease(dev);

  double starttime = wkf_timer_timenow(osp_timer);


  //
  // create the main renderer object needed early on for 
  // instantiation of materials, lights, etc.
  // 
  const char *rstr = (osp_rendermode == RT_SCIVIS) ? "scivis" : "pathtracer";
  if ((ospRenderer = ospNewRenderer(rstr)) == NULL) {
    printf("OSPRay2Renderer) Failed to load OSPRay renderer '%s'!\n", rstr);
  } 
  if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) created renderer '%s'\n", rstr);
  }

  // load and initialize all of the materials
  init_materials();

  time_ctx_create = wkf_timer_timenow(osp_timer) - starttime;
  
  if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) context creation time: %.2f\n", time_ctx_create);
  }

  context_created = 1;
}

        
/// destructor
OSPRay2Renderer::~OSPRay2Renderer(void) {
  DBG();

  destroy_scene();

  if (context_created && (ospRenderer != NULL))
    ospRelease(ospRenderer);

//  if (dev != NULL)
//    ospDeviceRelease(dev);

  wkf_timer_destroy(osp_timer);
}


void OSPRay2Renderer::check_verbose_env() {
  DBG();

  char *verbstr = getenv("VMDOSPRAYVERBOSE");
  if (verbstr != NULL) {
//    printf("OSPRay2Renderer) verbosity config request: '%s'\n", verbstr);
    if (!strupcmp(verbstr, "MIN")) {
      verbose = RT_VERB_MIN;
      printf("OSPRay2Renderer) verbose setting: minimum\n");
    } else if (!strupcmp(verbstr, "TIMING")) {
      verbose = RT_VERB_TIMING;
      printf("OSPRay2Renderer) verbose setting: timing data\n");
    } else if (!strupcmp(verbstr, "DEBUG")) {
      verbose = RT_VERB_DEBUG;
      printf("OSPRay2Renderer) verbose setting: full debugging data\n");
    }
  }
}


void OSPRay2Renderer::setup_context(int w, int h) {
  DBG();
  double starttime = wkf_timer_timenow(osp_timer);
  time_ctx_setup = 0;

  lasterror = 0; /* XXX SUCCESS; */ // clear any error state
  width = w;
  height = h;

  if (!context_created)
    return;

  check_verbose_env(); // update verbose flag if changed since last run

  // maxPathLength -- supported by all renderers
  if (getenv("VMDOSPRAYMAXDEPTH")) {
    int maxdepth = atoi(getenv("VMDOSPRAYMAXDEPTH"));
    if (maxdepth > 0 && maxdepth <= 20) {
      printf("OSPRay2Renderer) Setting maxdepth to %d...\n", maxdepth);
      ospSetParam(ospRenderer, "maxPathLength", OSP_INT, &maxdepth);  
    } else {
      printf("OSPRay2Renderer) ignoring out-of-range maxdepth: %d...\n", maxdepth);
    }
  } else {
    int maxdepth = 20;
    ospSetParam(ospRenderer, "maxPathLength", OSP_INT, &maxdepth);  
  }

#if 0
  // XXX -- not implemented in OSPRay 2.x presently
  // Implore OSPRay to correctly handle lighting through transparent 
  // surfaces when AO is enabled
  const int one = 1;
  ospSetParam(ospRenderer, "aoTransparencyEnabled", OSP_INT, &one);
#endif

  time_ctx_setup = wkf_timer_timenow(osp_timer) - starttime;
}


void OSPRay2Renderer::destroy_scene() {
  DBG();

  double starttime = wkf_timer_timenow(osp_timer);
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
    free(cylinders_color[i].ind);
    cylinders_color[i].ind = NULL;
    free(cylinders_color[i].cols);
    cylinders_color[i].cols = NULL;
  }
  cylinders_color.clear();

  ospInstances.clear();

  for (i=0; i<materialcache.num(); i++) {
    ospRelease(materialcache[i].mat);
  }
  materialcache.clear();

  int lcnt = ospLights.num();
  for (i = 0; i < lcnt; ++i) {
    ospRelease(ospLights[i]);
  }
  ospLights.clear();

  if (ospCamera != NULL) {
    ospRelease(ospCamera);
    ospCamera = NULL;
  }

  if (context_created && (ospWorld != NULL)) {
    ospRelease(ospWorld);
    ospWorld = NULL;
  }

  framebuffer_destroy();

  double endtime = wkf_timer_timenow(osp_timer);
  time_ctx_destroy_scene = endtime - starttime;

  scene_created = 0; // scene has been destroyed
}


void OSPRay2Renderer::update_rendering_state(int interactive) {
  DBG();
  if (!context_created)
    return;

  wkf_timer_start(osp_timer);

  // Set interactive/progressive rendering flag so that we wire up
  // the most appropriate renderer for the task.  For batch rendering
  // with AO, we would choose the largest possible sample batch size,
  // but for interactive we will always choose a batch size of 1 or maybe 2
  // to yield the best interactivity.
  interactive_renderer = interactive;

  // XXX set OSPRay rendering state

  long totaltris = tricolor_cnt + trimesh_c4u_n3b_v3f_cnt + 
                   trimesh_n3b_v3f_cnt + trimesh_n3f_v3f_cnt + trimesh_v3f_cnt;

  if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) cyl %ld, ring %ld, sph %ld, tri %ld, tot: %ld  lt %ld\n",
           cylinder_array_cnt + cylinder_array_color_cnt,
           ring_array_color_cnt,
           sphere_array_cnt + sphere_array_color_cnt,
           totaltris,
           cylinder_array_cnt +  cylinder_array_color_cnt + ring_array_color_cnt + sphere_array_cnt + sphere_array_color_cnt + totaltris,
           directional_lights.num() + positional_lights.num());
  }

  if (verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) using fully general shader and materials.\n");
  }

  // XXX set OSPRay background color

  if (verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) scene bg mode: %d\n", scene_background_mode);

    printf("OSPRay2Renderer) scene bgsolid: %.2f %.2f %.2f\n", 
           scene_bg_color[0], scene_bg_color[1], scene_bg_color[2]);

    printf("OSPRay2Renderer) scene bggradT: %.2f %.2f %.2f\n", 
           scene_bg_grad_top[0], scene_bg_grad_top[1], scene_bg_grad_top[2]);

    printf("OSPRay2Renderer) scene bggradB: %.2f %.2f %.2f\n", 
           scene_bg_grad_bot[0], scene_bg_grad_bot[1], scene_bg_grad_bot[2]);
  
    printf("OSPRay2Renderer) bg gradient: %f %f %f  top: %f  bot: %f\n",
           scene_gradient[0], scene_gradient[1], scene_gradient[2],
           scene_gradient_topval, scene_gradient_botval);
  }

  // update in case the caller changed top/bottom values since last recalc
  scene_gradient_invrange = 1.0f / (scene_gradient_topval - scene_gradient_botval);
  // XXX set OSPRay background gradient

  // XXX set OSPRay fog mode

  if (verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) adding lights: dir: %ld  pos: %ld\n", 
           directional_lights.num(), positional_lights.num());
  }

  // XXX set OSPRay lights

  if (verbose == RT_VERB_DEBUG) 
    printf("OSPRay2Renderer) Finalizing OSPRay scene graph...\n");

  // create group to hold instances

  // XXX we should create an acceleration object the instance shared
  //     by multiple PBC images


  // XXX OSPRay AS builder initialization if there's any customization...

  // do final state variable updates before rendering begins
  if (verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) cam zoom factor %f\n", cam_zoom);
    printf("OSPRay2Renderer) cam stereo eye separation  %f\n", cam_stereo_eyesep);
    printf("OSPRay2Renderer) cam stereo convergence distance %f\n", 
           cam_stereo_convergence_dist);
    printf("OSPRay2Renderer) cam DoF focal distance %f\n", cam_dof_focal_dist);
    printf("OSPRay2Renderer) cam DoF f/stop %f\n", cam_dof_fnumber);
  }

  // define all of the standard camera params
  // XXX set OSPRay camera state

  // define stereoscopic camera parameters
  // XXX set OSPRay camera state

  // define camera DoF parameters
  // XXX set OSPRay camera state

  // XXX set OSPRay AO sample counts and light scaling factors

  if (verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) setting sample counts:  AA %d  AO %d\n", aa_samples, ao_samples);
    printf("OSPRay2Renderer) setting AO factors:  AOA %f  AOD %f\n", ao_ambient, ao_direct);
  }

  //
  // Handle AA samples either internally with loops internal to 
  // each ray launch point thread, or externally by iterating over
  // multiple launches, adding each sample to an accumulation buffer,
  // or a hybrid combination of the two.  
  //
#if 1
  ext_aa_loops = 1;
#else
  ext_aa_loops = 1;
  if (ao_samples > 0 || (aa_samples > 4)) {
    // if we have too much work for a single-pass rendering, we need to 
    // break it up into multiple passes of the right counts in each pass
    ext_aa_loops = 1 + aa_samples;
    // XXX set OSPRay sample counts per launch...
  } else { 
    // if the scene is simple, e.g. no AO rays and AA sample count is small,
    // we can run it in a single pass and get better performance
    // XXX set OSPRay sample counts per launch...
  }
  // XXX set OSPRay accum buf normalization scaling factors
#endif

  if (verbose == RT_VERB_DEBUG) {
    if (ext_aa_loops > 1)
      printf("OSPRay2Renderer) Running OSPRay multi-pass: %d loops\n", ext_aa_loops);
    else
      printf("OSPRay2Renderer) Running OSPRay single-pass: %d total samples\n", 1+aa_samples);
  }

  // set the ray generation program to the active camera code...
  // XXX set OSPRay camera mode and clear accum buf
  // set the active color accumulation ray gen program based on the 
  // camera/projection mode, stereoscopic display mode, 
  // and depth-of-field state
  // XXX set OSPRay camera mode and accum buf mode
  // XXX set OSPRay "miss" shading mode (solid or gradient)
}


void OSPRay2Renderer::framebuffer_config(int fbwidth, int fbheight) {
  if (!context_created)
    return;

  width = fbwidth;
  height = fbheight;

  // allocate and resize buffers to match request
  if (buffers_allocated) {
    // if the buffers already exist and match the current 
    // progressive/non-progressive rendering mode, just resize them
    if (verbose == RT_VERB_DEBUG) {
      printf("OSPRay2Renderer) resizing framebuffer\n");
    }
    framebuffer_resize(width, height);
  } else {
    // (re)allocate framebuffer and associated accumulation buffers if they
    // don't already exist or if they weren't bound properly for
    // current progressive/non-progressive rendering needs.
    if (verbose == RT_VERB_DEBUG) {
      printf("OSPRay2Renderer) creating framebuffer and accum. buffer\n");
    }

    // create intermediate output and accumulation buffers
    ospFrameBuffer = ospNewFrameBuffer(width, height, OSP_FB_RGBA8, OSP_FB_COLOR | OSP_FB_ACCUM);
    ospCommit(ospFrameBuffer);
    ospResetAccumulation(ospFrameBuffer);

    buffers_allocated = 1;
  }
}


void OSPRay2Renderer::framebuffer_resize(int fbwidth, int fbheight) {
  if (!context_created)
    return;

  width = fbwidth;
  height = fbheight;

  if (buffers_allocated) {
    if (verbose == RT_VERB_DEBUG) 
      printf("OSPRay2Renderer) framebuffer_resize(%d x %d)\n", width, height);
    framebuffer_destroy();
  }

  ospFrameBuffer = ospNewFrameBuffer(width, height, OSP_FB_RGBA8, OSP_FB_COLOR | OSP_FB_ACCUM);
  ospCommit(ospFrameBuffer);
  ospResetAccumulation(ospFrameBuffer);
  buffers_allocated = 1;
}


void OSPRay2Renderer::framebuffer_destroy() {
  if (!context_created)
    return;

  if (buffers_allocated) {
    if (ospFrameBuffer)
      ospRelease(ospFrameBuffer);
  }
  buffers_allocated = 0;
}


void OSPRay2Renderer::render_compile_and_validate(void) {
  int i;

  DBG();
  if (!context_created)
    return;

  //
  // finalize context validation, compilation, and AS generation 
  //
  double startctxtime = wkf_timer_timenow(osp_timer);

  // XXX any last OSPRay state updates/checks

  if ((ospWorld = ospNewWorld()) == NULL) {
    printf("OSPRay2Renderer) Failed to create new world!\n");
  }

  if (verbose == RT_VERB_DEBUG)
    printf("OSPRayReenderer) num spheres = %ld\n", spheres_color.num());


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

  if (camera_projection == OSPRay2Renderer::RT_ORTHOGRAPHIC) {
    if(!ospCamera) ospCamera = ospNewCamera("orthographic");

    float camaspect = width / ((float) height); 
    ospSetParam(ospCamera, "aspect", OSP_FLOAT, &camaspect);
    float orthoheight = 2.0f * cam_zoom;
    ospSetParam(ospCamera, "height", OSP_FLOAT, &orthoheight);

    if (dof_enabled) {
      msgWarn << "OSPRay2Renderer) DoF not implemented for orthographic camera!" << sendmsg;
    }
  } else {
    if(!ospCamera) ospCamera = ospNewCamera("perspective");

    float camaspect = width / ((float) height); 
    ospSetParam(ospCamera, "aspect", OSP_FLOAT, &camaspect);
    float camfovy = 2.0f*180.0f*(atanf(cam_zoom)/float(M_PI));
    ospSetParam(ospCamera, "fovy", OSP_FLOAT, &camfovy);

    if (dof_enabled) {
      ospSetParam(ospCamera, "focusDistance", OSP_FLOAT, &cam_dof_focal_dist);
      float camaprad = cam_dof_focal_dist / (2.0f * cam_zoom * cam_dof_fnumber);
      ospSetParam(ospCamera, "apertureRadius", OSP_FLOAT, &camaprad);
    } else {
      const float zero = 0.0f;
      ospSetParam(ospCamera, "apertureRadius", OSP_FLOAT, &zero);
    }
  }

  if (ospCamera) {
    ospSetParam(ospCamera, "position",  OSP_VEC3F, cam_pos);
    ospSetParam(ospCamera, "direction", OSP_VEC3F, cam_W);
    ospSetParam(ospCamera, "up",        OSP_VEC3F, cam_V);
    ospCommit(ospCamera);
  }

  // 
  // Set framebuffer 
  // 
  framebuffer_config(width, height);

  //
  // Set all lights
  //

  // The direct lighting scaling factor all of the other lights.
  float lightscale = 1.0f;
#if 0
  if (ao_samples != 0)
    lightscale = ao_direct;
#endif

  for (i = 0; i < directional_lights.num(); ++i) {
    OSPLight light = ospNewLight("distant");

    // The direct lighting scaling factor is applied to the lights here.
    ospSetParam(light, "intensity", OSP_FLOAT, &lightscale);
    ospSetParam(light, "color", OSP_VEC3F, directional_lights[i].color);

    // OSPRay uses a light direction vector opposite to VMD and Tachyon 
    float lightDir[3];
    vec_negate(lightDir, directional_lights[i].dir);
    vec_normalize(lightDir); // just for good measure
    ospSetParam(light, "direction", OSP_VEC3F, lightDir);
    ospCommit(light);
    ospLights.append(light);
  }

  // AO scaling factor is applied to a special ambient light.
  if (ao_samples != 0) {
    OSPLight light = ospNewLight("ambient");

    // AO scaling factor is applied to the special ambient light
    ospSetParam(light, "intensity", OSP_FLOAT, &ao_ambient);
    float whitecol[] = { 1.0f, 1.0f, 1.0f };
    ospSetParam(light, "color", OSP_VEC3F, whitecol);
    ospCommit(light);
    ospLights.append(light); // add AO ambient light
  } 

  // 
  // update renderer state
  //
  ospSetParam(ospRenderer, "backgroundColor", OSP_VEC3F, scene_bg_color);

  if (ao_samples && interactive_renderer) {
    const int one = 1;
    ospSetParam(ospRenderer, "pixelSamples", OSP_INT, &one); // all renderers
    if (osp_rendermode == RT_SCIVIS)
      ospSetParam(ospRenderer, "aoSamples", OSP_INT, &one);  // scivis-only
  } else {
    ospSetParam(ospRenderer, "pixelSamples", OSP_INT, &aa_samples); // all renderers
    if (osp_rendermode == RT_SCIVIS)
      ospSetParam(ospRenderer, "aoSamples", OSP_INT, &ao_samples); // scivis-only
  }

  if (getenv("VMDOSPRAYAOMAXDIST")) {
    float tmp = float(atof(getenv("VMDOSPRAYAOMAXDIST")));
    if (verbose == RT_VERB_DEBUG) {
      printf("OSPRay2Renderer) setting AO maxdist: %f\n", tmp);
    }
    ospSetParam(ospRenderer, "aoRadius", OSP_FLOAT, &tmp); // scivis-only
  }

#if 1
  // XXX OSPRay 2.x doesn't support rendering w/o shadows presently
  // render with/without shadows
    msgInfo << "Shadow rendering enabled." << sendmsg;
#else
  if (shadows_enabled || ao_samples) {
    if (shadows_enabled && !ao_samples)
      msgInfo << "Shadow rendering enabled." << sendmsg;

    const int one = 1;
//    ospSetParam(ospRenderer, "shadowsEnabled", OSP_INT, &one);
  } else {
    const int zero = 0;
//    ospSetParam(ospRenderer, "shadowsEnabled", OSP_INT, &zero);
  }
#endif

  // render with ambient occlusion, but only if shadows are also enabled
  if (ao_samples) {
    msgInfo << "Ambient occlusion enabled." << sendmsg;
//    msgInfo << "Shadow rendering enabled." << sendmsg;
  }

  // commit triangle mesh geometry after assigning materials
  for (i=0; i<trimesh_v3f_n3f_c3f.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("OSPRay2Renderer) Adding triangle mesh[%d]: %d tris ...\n", 
             i, trimesh_v3f_n3f_c3f[i].num);

    OSPGroup group = ospNewGroup();
    OSPData geometricModels = ospNewSharedData(&trimesh_v3f_n3f_c3f[i].model, OSP_GEOMETRIC_MODEL, 1, 0);
    ospCommit(geometricModels);
    ospRelease(trimesh_v3f_n3f_c3f[i].model);
    ospSetParam(group, "geometry", OSP_DATA, &geometricModels);
    ospCommit(group);
    ospRelease(geometricModels);

    OSPInstance instance = ospNewInstance(group);
    ospCommit(instance); 
    ospRelease(group);

    ospInstances.append(instance);
  } 

  // commit sphere geometry after assigning materials
  for (i=0; i<spheres_color.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("OSPRay2Renderer) Adding sphere_color array [%d]: %d spheres ...\n",
             i, spheres_color[i].num);

    OSPGroup group = ospNewGroup();
    OSPData geometricModels = ospNewSharedData(&spheres_color[i].model, OSP_GEOMETRIC_MODEL, 1, 0);
    ospCommit(geometricModels);
    ospRelease(spheres_color[i].model);
    ospSetParam(group, "geometry", OSP_DATA, &geometricModels);
    ospRelease(geometricModels);
    ospCommit(group);

    OSPInstance instance = ospNewInstance(group);
    ospCommit(instance); 
    ospRelease(group);

    ospInstances.append(instance);
  } 

  // commit cylinder geometry after assigning materials
  for (i=0; i<cylinders_color.num(); i++) {
    if (verbose == RT_VERB_DEBUG)
      printf("OSPRay2Renderer) Adding cylinders_color array [%d]: %d cyls...\n",
             i, cylinders_color[i].num);

    OSPGroup group = ospNewGroup();
    OSPData geometricModels = ospNewSharedData(&cylinders_color[i].model, OSP_GEOMETRIC_MODEL, 1, 0);
    ospCommit(geometricModels);
    ospRelease(cylinders_color[i].model);
    ospSetParam(group, "geometry", OSP_DATA, &geometricModels);
    ospRelease(geometricModels);
    ospCommit(group);

    OSPInstance instance = ospNewInstance(group);
    ospCommit(instance); 
    ospRelease(group);

    ospInstances.append(instance);
  }

  // attach all instances to the scene...
  OSPData instances = ospNewSharedData(&ospInstances[0], OSP_INSTANCE, ospInstances.num(), 0);
  ospCommit(instances);
  ospSetParam(ospWorld, "instance", OSP_DATA, &instances);
  ospRelease(instances);
  for (i=0; i<ospInstances.num(); i++) {
    ospRelease(ospInstances[i]);
  }
  ospInstances.clear();

  if (ospLights.num() > 0) {
    // attach all lights to the scene...
    OSPData lights = ospNewSharedData(&ospLights[0], OSP_LIGHT, ospLights.num(), 0);
    ospCommit(lights);
    ospSetParam(ospWorld, "light", OSP_DATA, &lights);
    ospCommit(ospWorld); // commit the completed scene
    ospRelease(lights);
  } else {
    ospCommit(ospWorld); // commit the completed scene
  }

  // print out world bounds
  OSPBounds worldBounds = ospGetBounds(ospWorld);
  printf("OSPRay2Renderer) world bounds: ({%f, %f, %f}, {%f, %f, %f}\n\n",
         worldBounds.lower[0], worldBounds.lower[1], worldBounds.lower[2],
         worldBounds.upper[0], worldBounds.upper[1], worldBounds.upper[2]);

  ospCommit(ospRenderer);


  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) Finalizing OSPRay rendering kernels...\n");
  // XXX any last OSPRay state updates/checks

  double contextinittime = wkf_timer_timenow(osp_timer);
  time_ctx_validate = contextinittime - startctxtime;

  //
  // Force OSPRay to build the acceleration structure _now_, so we can time it
  //
  // XXX No way to force-build OSPRay AS for timing?

  time_ctx_AS_build = wkf_timer_timenow(osp_timer) - contextinittime;
  if (verbose == RT_VERB_DEBUG) {
    printf("OSPRay2Renderer) launching render: %d x %d\n", width, height);
  }
}


#if defined(VMDOSPRAY_INTERACTIVE_OPENGL)

static void *createospraywindow(const char *wintitle, int width, int height) {
  printf("OSPRay2Renderer) Creating OSPRay window: %d x %d...\n", width, height);

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
  printf("OSPRay2Renderer) VMD TachyonL-OSPRay Interactive Ray Tracer help:\n");
  printf("OSPRay2Renderer) ================================================\n");

  // check for Spaceball/SpaceNavigator/Magellan input devices
  int havespaceball = ((glwin_spaceball_available(win)) && (getenv("VMDDISABLESPACEBALLXDRV") == NULL));
  printf("OSPRay2Renderer) Spaceball/SpaceNavigator/Magellan: %s\n",
         (havespaceball) ? "Available" : "Not available");

  // check for stereo-capable display
  int havestereo, havestencil;
  glwin_get_wininfo(win, &havestereo, &havestencil);
  printf("OSPRay2Renderer) Stereoscopic display: %s\n",
         (havestereo) ? "Available" : "Not available");

  // check for vertical retrace sync
  int vsync=0, rc=0;
  if ((rc = glwin_query_vsync(win, &vsync)) == GLWIN_SUCCESS) {
    printf("OSPRay2Renderer) Vert retrace sync: %s\n", (vsync) ? "On" : "Off");
  } else {
    printf("OSPRay2Renderer) Vert retrace sync: indeterminate\n");
  }

  printf("OSPRay2Renderer)\n");
  printf("OSPRay2Renderer) General controls:\n");
  printf("OSPRay2Renderer)   space: save numbered snapshot image\n");
  printf("OSPRay2Renderer)       =: reset to initial view\n");
  printf("OSPRay2Renderer)       h: print this help info\n");
  printf("OSPRay2Renderer)       p: print current rendering parameters\n");
  printf("OSPRay2Renderer)   ESC,q: quit viewer\n");
  printf("OSPRay2Renderer)\n");
  printf("OSPRay2Renderer) Display controls\n");
  printf("OSPRay2Renderer)      F1: override shadows on/off (off=AO off too)\n");
  printf("OSPRay2Renderer)      F2: override AO on/off\n");
  printf("OSPRay2Renderer)      F3: override DoF on/off\n");
  printf("OSPRay2Renderer)      F4: override Depth cueing on/off\n");
// Not currently applicable to OSPRay
// #ifdef USE_REVERSE_SHADOW_RAYS
//   printf("OSPRay2Renderer)      F5: enable/disable shadow ray optimizations\n");
// #endif
  printf("OSPRay2Renderer)     F12: toggle full-screen display on/off\n");
  printf("OSPRay2Renderer)   1-9,0: override samples per update auto-FPS off\n");
  printf("OSPRay2Renderer)      Up: increase DoF focal distance\n");
  printf("OSPRay2Renderer)    Down: decrease DoF focal distance\n");
  printf("OSPRay2Renderer)    Left: decrease DoF f/stop\n");
  printf("OSPRay2Renderer)   Right: increase DoF f/stop\n");
  printf("OSPRay2Renderer)       S: toggle stereoscopic display on/off (if avail)\n");
  printf("OSPRay2Renderer)       a: toggle AA/AO auto-FPS tuning on/off (on)\n");
  printf("OSPRay2Renderer)       g: toggle gradient sky xforms on/off (on)\n");
  printf("OSPRay2Renderer)       l: toggle light xforms on/off (on)\n");
  printf("OSPRay2Renderer)\n");
  printf("OSPRay2Renderer) Mouse controls:\n");
  printf("OSPRay2Renderer)       f: mouse depth-of-field mode\n");
  printf("OSPRay2Renderer)       r: mouse rotation mode\n");
  printf("OSPRay2Renderer)       s: mouse scaling mode\n");
  printf("OSPRay2Renderer)       t: mouse translation mode\n");

  int movie_recording_enabled = (getenv("VMDOSPRAYLIVEMOVIECAPTURE") != NULL);
  if (movie_recording_enabled) {
    printf("OSPRay2Renderer)\n");
    printf("OSPRay2Renderer) Movie recording controls:\n");
    printf("OSPRay2Renderer)       R: start/stop movie recording\n");
    printf("OSPRay2Renderer)       F: toggle movie FPS (24, 30, 60)\n");
  }
}


void OSPRay2Renderer::render_to_glwin(const char *filename) {
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
  int movie_recording_enabled = (getenv("VMDOSPRAYLIVEMOVIECAPTURE") != NULL);
  int movie_recording_on = 0;
  double movie_recording_start_time = 0.0;
  int movie_recording_fps = 30;
  int movie_framecount = 0;
  int movie_lastframeindex = 0;
  const char *movie_recording_filebase = "vmdlivemovie.%05d.tga";
  if (getenv("VMDOSPRAYLIVEMOVIECAPTUREFILEBASE"))
    movie_recording_filebase = getenv("VMDOSPRAYLIVEMOVIECAPTUREFILEBASE");

  // Enable/disable Spaceball/SpaceNavigator/Magellan input 
  int spaceballenabled=(getenv("VMDDISABLESPACEBALLXDRV") == NULL) ? 1 : 0;
  int spaceballmode=0;       // default mode is rotation/translation
  int spaceballflightmode=0; // 0=moves object, 1=camera fly
  if (getenv("VMDOSPRAYSPACEBALLFLIGHT"))
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
  const char *imageszstr = getenv("VMDOSPRAYIMAGESIZE");
  if (imageszstr) {
    if (sscanf(imageszstr, "%d %d", &width, &height) != 2) {
      width=wsx;
      height=wsy;
    } 
  } 
  framebuffer_config(width, height);

  // prepare the majority of OSPRay rendering state before we go into 
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
  osp_directional_light *cur_dlights = (osp_directional_light *) calloc(1, directional_lights.num() * sizeof(osp_directional_light));
  for (i=0; i<directional_lights.num(); i++) {
    vec_copy((float*)&cur_dlights[i].color, directional_lights[i].color);
    vec_copy((float*)&cur_dlights[i].dir, directional_lights[i].dir);
    vec_normalize((float*)&cur_dlights[i].dir);
  }

  // create the display window
  void *win = createospraywindow("VMD TachyonL-OSPRay Interactive Ray Tracer", width, height);
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
  
  double oldtime = wkf_timer_timenow(osp_timer);
  while (!done) { 
    int winevent=0;

    while ((winevent = glwin_handle_events(win, GLWIN_EV_POLL_NONBLOCK)) != 0) {
      int evdev, evval;
      char evkey;

      glwin_get_lastevent(win, &evdev, &evval, &evkey);
      glwin_get_winsize(win, &wsx, &wsy);

      if (evdev == GLWIN_EV_WINDOW_CLOSE) {
        printf("OSPRay2Renderer) display window closed, exiting...\n");
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
              const unsigned char *FB = (const unsigned char*)ospMapFrameBuffer(ospFrameBuffer, OSP_FB_COLOR);
              if (write_image_file_rgb4u(snapfilename, FB, width, height)) {
                printf("OSPRay2Renderer) Failed to write output image!\n");
              } else {
                printf("OSPRay2Renderer) Saved snapshot to '%s'             \n",
                       snapfilename);
              }
              ospUnmapFrameBuffer(FB, ospFrameBuffer);
              snapshotcount++; 
            }
            break;

          case  'a': /* toggle automatic sample count FPS tuning */
            autosamplecount = !(autosamplecount);
            printf("\nOSPRay2Renderer) Automatic AO sample count FPS tuning %s\n",
                   (autosamplecount) ? "enabled" : "disabled");
            break;

          case  'f': /* DoF mode */
            mm = RTMM_DOF;
            printf("\nOSPRay2Renderer) Mouse DoF aperture and focal dist. mode\n");
            break;

          case  'g': /* toggle gradient sky sphere xforms */
            xformgradientsphere = !(xformgradientsphere);
            printf("\nOSPRay2Renderer) Gradient sky sphere transformations %s\n",
                   (xformgradientsphere) ? "enabled" : "disabled");
            break;

          case  'h': /* print help message */
            printf("\n");
            interactive_viewer_usage(win);
            break;

          case  'l': /* toggle lighting xforms */
            xformlights = !(xformlights);
            printf("\nOSPRay2Renderer) Light transformations %s\n",
                   (xformlights) ? "enabled" : "disabled");
            break;

          case  'p': /* print current RT settings */
            printf("\nOSPRay2Renderer) Current Ray Tracing Parameters:\n"); 
            printf("OSPRay2Renderer) -------------------------------\n"); 
            printf("OSPRay2Renderer) Camera zoom: %f\n", cur_cam_zoom);
            printf("OSPRay2Renderer) Shadows: %s  Ambient occlusion: %s\n",
                   (gl_shadows_on) ? "on" : "off",
                   (gl_ao_on) ? "on" : "off");
            printf("OSPRay2Renderer) Antialiasing samples per-pass: %d\n",
                   cur_aa_samples);
            printf("OSPRay2Renderer) Ambient occlusion samples per-pass: %d\n",
                   cur_ao_samples);
            printf("OSPRay2Renderer) Depth-of-Field: %s f/num: %.1f  Foc. Dist: %.2f\n",
                   (gl_dof_on) ? "on" : "off", 
                   cam_dof_fnumber, cam_dof_focal_dist);
            printf("OSPRay2Renderer) Image size: %d x %d\n", width, height);
            break;

          case  'r': /* rotate mode */
            mm = RTMM_ROT;
            printf("\nOSPRay2Renderer) Mouse rotation mode\n");
            break;

          case  's': /* scaling mode */
            mm = RTMM_SCALE;
            printf("\nOSPRay2Renderer) Mouse scaling mode\n");
            break;

          case  'F': /* toggle live movie recording FPS (24, 30, 60) */
            if (movie_recording_enabled) {
              switch (movie_recording_fps) {
                case 24: movie_recording_fps = 30; break;
                case 30: movie_recording_fps = 60; break;
                case 60:
                default: movie_recording_fps = 24; break;
              }
              printf("\nOSPRay2Renderer) Movie recording FPS rate: %d\n", 
                     movie_recording_fps);
            } else {
              printf("\nOSPRay2Renderer) Movie recording not available.\n");
            }
            break;

          case  'R': /* toggle live movie recording mode on/off */
            if (movie_recording_enabled) {
              movie_recording_on = !(movie_recording_on);
              printf("\nOSPRay2Renderer) Movie recording %s\n",
                     (movie_recording_on) ? "STARTED" : "STOPPED");
              if (movie_recording_on) {
                movie_recording_start_time = wkf_timer_timenow(osp_timer);
                movie_framecount = 0;
                movie_lastframeindex = 0;
              } else {
                printf("OSPRay2Renderer) Encode movie with:\n");
                printf("OSPRay2Renderer)   ffmpeg -f image2 -i vmdlivemovie.%%05d.tga -c:v libx264 -profile:v baseline -level 3.0 -pix_fmt yuv420p -b:v 15000000 output.mp4\n");
              }
            } else {
              printf("\nOSPRay2Renderer) Movie recording not available.\n");
            }
            break;

          case  'S': /* toggle stereoscopic display mode */
            if (havestereo) {
              stereoon = (!stereoon);
              printf("\nOSPRay2Renderer) Stereoscopic display %s\n",
                     (stereoon) ? "enabled" : "disabled");
              winredraw = 1;
            } else {
              printf("\nOSPRay2Renderer) Stereoscopic display unavailable\n");
            }
            break;
 
          case  't': /* translation mode */
            mm = RTMM_TRANS;
            printf("\nOSPRay2Renderer) Mouse translation mode\n");
            break;
            
          case  'q': /* 'q' key */
          case  'Q': /* 'Q' key */
          case 0x1b: /* ESC key */
            printf("\nOSPRay2Renderer) Exiting on user input.               \n");
            done=1; /* exit from interactive RT window */
            break;
        }
      } else if (evdev != GLWIN_EV_NONE) {
        switch (evdev) {
          case GLWIN_EV_KBD_F1: /* turn shadows on/off */
            gl_shadows_on=(!gl_shadows_on) ? RT_SHADOWS_ON : RT_SHADOWS_OFF;
            // gl_shadows_on = (!gl_shadows_on);
            printf("\n");
            printf("OSPRay2Renderer) Shadows %s\n",
                   (gl_shadows_on) ? "enabled" : "disabled");
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_F2: /* turn AO on/off */
            gl_ao_on = (!gl_ao_on); 
            printf("\n");
            printf("OSPRay2Renderer) Ambient occlusion %s\n",
                   (gl_ao_on) ? "enabled" : "disabled");
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_F3: /* turn DoF on/off */
            gl_dof_on = (!gl_dof_on);
            printf("\n");
            if ((camera_projection == RT_ORTHOGRAPHIC) && gl_dof_on) {
              gl_dof_on=0; 
              printf("OSPRay2Renderer) Depth-of-field not available in orthographic mode\n");
            }
            printf("OSPRay2Renderer) Depth-of-field %s\n",
                   (gl_dof_on) ? "enabled" : "disabled");
            winredraw = 1;
            break;

          case GLWIN_EV_KBD_F4: /* turn fog/depth cueing on/off */
            gl_fog_on = (!gl_fog_on); 
            printf("\n");
            printf("OSPRay2Renderer) Depth cueing %s\n",
                   (gl_fog_on) ? "enabled" : "disabled");
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_F12: /* toggle full-screen window on/off */
            gl_fs_on = (!gl_fs_on);
            printf("\nOSPRay2Renderer) Toggling fullscreen window %s\n",
                   (gl_fs_on) ? "on" : "off");
            if (gl_fs_on) { 
              if (glwin_fullscreen(win, gl_fs_on, 0) == 0) {
                owsx = wsx;
                owsy = wsy;
                glwin_get_winsize(win, &wsx, &wsy);
              } else {
                printf("OSPRay2Renderer) Fullscreen mode note available\n");
              }
            } else {
              glwin_fullscreen(win, gl_fs_on, 0);
              glwin_resize(win, owsx, owsy);
            }
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_UP: /* change depth-of-field focal dist */
            cam_dof_focal_dist *= 1.02f; 
            printf("\nOSPRay2Renderer) DoF focal dist: %f\n", cam_dof_focal_dist);
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_DOWN: /* change depth-of-field focal dist */
            cam_dof_focal_dist *= 0.96f; 
            if (cam_dof_focal_dist < 0.02f) cam_dof_focal_dist = 0.02f;
            printf("\nOSPRay2Renderer) DoF focal dist: %f\n", cam_dof_focal_dist);
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_RIGHT: /* change depth-of-field f/stop number */
            cam_dof_fnumber += 1.0f; 
            printf("\nOSPRay2Renderer) DoF f/stop: %f\n", cam_dof_fnumber);
            winredraw = 1; 
            break;

          case GLWIN_EV_KBD_LEFT: /* change depth-of-field f/stop number */
            cam_dof_fnumber -= 1.0f; 
            if (cam_dof_fnumber < 1.0f) cam_dof_fnumber = 1.0f;
            printf("\nOSPRay2Renderer) DoF f/stop: %f\n", cam_dof_fnumber);
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
                float txdx = (x - mousedownx) * 2.0f / wsx;
                float zoominc = 1.0f - txdx;
                if (zoominc < 0.01f) zoominc = 0.01f;
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
          printf("OSPRay2Renderer) spaceball button 1 pressed: reset view\n");
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
          printf("OSPRay2Renderer) spaceball mode: %s                       \n",
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
          if (zoominc < 0.01f) zoominc = 0.01f;
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
    // destroy and recreate affected OSPRay buffers
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
        printf("\rOSPRay2Renderer) Window resize: %d x %d                               \n", width, height);
      }

      winredraw=1;
    }

    int frame_ready = 1; // Default to true
    unsigned int subframe_count = 1;
    if (!done) {
      //
      // If the user interacted with the window in a meaningful way, we
      // need to update the OSPRay rendering state, recompile and re-validate
      // the context, and then re-render...
      //
      if (winredraw) {
        // update camera parameters
        ospSetParam(ospCamera, "position",  OSP_VEC3F, cam_pos);
        ospSetParam(ospCamera, "direction", OSP_VEC3F,   hmd_W);
        ospSetParam(ospCamera, "up",        OSP_VEC3F,   hmd_V);
        float camaspect = width / ((float) height);
        ospSetParam(ospCamera, "aspect", OSP_FLOAT, &camaspect);
        float camfovy = 2.0f*180.0f*(atanf(cam_zoom)/float(M_PI));
        ospSetParam(ospCamera, "fovy", OSP_FLOAT, &camfovy);
 
        // update shadow state 
//        ospSetParam(ospRenderer, "shadowsEnabled", OSP_INT, &gl_shadows_on);

        // update AO state 
        if (gl_shadows_on && gl_ao_on) {
          const int one = 1;
          if (osp_rendermode == RT_SCIVIS)
            ospSetParam(ospRenderer, "aoSamples", OSP_INT, &one);
        } else {
          const int zero = 0;
          if (osp_rendermode == RT_SCIVIS)
            ospSetParam(ospRenderer, "aoSamples", OSP_INT, &zero);
        }

        // update depth cueing state
        // XXX update OSPRay depth cueing state
 
        // update/recompute DoF values 
        // XXX OSPRay only implements DoF for the perspective
        //     camera at the present time
        if (camera_projection == OSPRay2Renderer::RT_PERSPECTIVE) {
          if (gl_dof_on) {
            ospSetParam(ospCamera, "focusDistance", OSP_FLOAT, &cam_dof_focal_dist);
            float camaprad = cam_dof_focal_dist / (2.0f * cam_zoom * cam_dof_fnumber);
            ospSetParam(ospCamera, "apertureRadius", OSP_FLOAT, &camaprad);
          } else {
            const float zero = 0.0f;
            ospSetParam(ospCamera, "apertureRadius", OSP_FLOAT, &zero);
          }
        }

        // commit camera updates once they're all done...
        ospCommit(ospCamera);

        //
        // Update light directions in the OSPRay light buffers
        //
        if (xformlights) {
          // AO scaling factor is applied at the renderer level, but
          // we apply the direct lighting scaling factor to the lights.
          float lightscale = 1.0f;
#if 0
          if (ao_samples != 0)
            lightscale = ao_direct;
#endif

          // XXX assumes the only contents in the first part of the 
          //     light list are directional lights.  The new AO "ambient"
          //     light is the last light in the list now, so we can get
          //     away with this, but refactoring is still needed here.
          for (i=0; i<directional_lights.num(); i++) {
            ospSetParam(ospLights[i], "intensity", OSP_FLOAT, &lightscale);
            ospSetParam(ospLights[i], "color", OSP_VEC3F, cur_dlights[i].color);
            float ltmp[3];
            vec_negate(ltmp, cur_dlights[i].dir);
            ospSetParam(ospLights[i], "direction", OSP_VEC3F, ltmp);
            ospCommit(ospLights[i]);
          }
        }

        // commit pending changes...
        ospCommit(ospRenderer);

        // reset accumulation buffer 
        accum_count=0;
        totalsamplecount=0;
        if (ospFrameBuffer != NULL) {
          ospResetAccumulation(ospFrameBuffer);
        }

        // 
        // Sample count updates and OSPRay state must always remain in 
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
        // XXX update OSPRay AA sample counts

        // observe latest AO enable/disable flag, and sample count
        if (gl_shadows_on && gl_ao_on) {
          // XXX update OSPRay AA/AO sample counts
        } else {
          cur_ao_samples = 0;
          // XXX update OSPRay AA/AO sample counts
        }
      } 


      // The accumulation buffer normalization factor must be updated
      // to reflect the total accumulation count before the accumulation
      // buffer is drawn to the output framebuffer
      // XXX update OSPRay accum buf normalization factor

      // The accumulation buffer subframe index must be updated to ensure that
      // the RNGs for AA and AO get correctly re-seeded
      // XXX update OSPRay accum subframe count

      // Force context compilation/validation
      // render_compile_and_validate();

      //
      // run the renderer 
      //
      frame_ready = 1; // Default to true
      subframe_count = 1;
      if (lasterror == 0 /* XXX SUCCESS */) {
        if (winredraw) {
          ospResetAccumulation(ospFrameBuffer);
          winredraw=0;
        }

        // iterate, adding to the accumulation buffer...
//printf("OSPRay2Renderer) ospRenderFrameBlocking(): [%d] ...\n", accum_sample);
#if 1
        OSPFuture fut;
        fut = ospRenderFrame(ospFrameBuffer, ospRenderer, ospCamera, ospWorld);
        ospWait(fut, OSP_FRAME_FINISHED);
        ospRelease(fut);
#else
        ospRenderFrameBlocking(ospFrameBuffer, ospRenderer, ospCamera, ospWorld);
#endif
        subframe_count++; // increment subframe index
        totalsamplecount += samples_per_pass;
        accum_count += cur_aa_samples;

        // copy the accumulation buffer image data to the framebuffer and
        // perform type conversion and normaliztion on the image data...
        // XXX launch OSPRay accum copy/norm/finish

        if (lasterror == 0 /* XXX SUCCESS */) {
          if (frame_ready) {
            // display output image
            const unsigned char * img;
            img = (const unsigned char*)ospMapFrameBuffer(ospFrameBuffer, OSP_FB_COLOR);

#if 0
            glwin_draw_image_tex_rgb3u(win, (stereoon!=0)*GLWIN_STEREO_OVERUNDER, width, height, img);
#else
            glwin_draw_image_rgb3u(win, (stereoon!=0)*GLWIN_STEREO_OVERUNDER, width, height, img);
#endif
            ospUnmapFrameBuffer(img, ospFrameBuffer);

            // if live movie recording is on, we save every displayed frame
            // to a sequence sequence of image files, with each file numbered
            // by its frame index, which is computed by the multiplying image
            // presentation time by the image sequence fixed-rate-FPS value.
            if (movie_recording_enabled && movie_recording_on) {
              char moviefilename[2048];

              // compute frame number from wall clock time and the
              // current fixed-rate movie playback frame rate
              double now = wkf_timer_timenow(osp_timer);
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
              const unsigned char *FB = (const unsigned char*)ospMapFrameBuffer(ospFrameBuffer, OSP_FB_COLOR);
              if (write_image_file_rgb4u(moviefilename, FB, width, height)) {
                movie_recording_on = 0;
                printf("\n");
                printf("OSPRay2Renderer) ERROR during writing image during movie recording!\n");
                printf("OSPRay2Renderer) Movie recording STOPPED\n");
              }
              ospUnmapFrameBuffer(FB, ospFrameBuffer);

              movie_lastframeindex = fidx; // update last frame index written
            }
          }
        } else {
          printf("OSPRay2Renderer) An error occured during rendering. Rendering is aborted.\n");
          done=1;
          break;
        }
      } else {
        printf("OSPRay2Renderer) An error occured in AS generation. Rendering is aborted.\n");
        done=1;
        break;
      }
    }

    if (!done && frame_ready) {
      double newtime = wkf_timer_timenow(osp_timer);
      double frametime = (newtime-oldtime) + 0.00001f;
      oldtime=newtime;

      // compute exponential moving average for exp(-1/10)
      double framefps = 1.0f/frametime;
      fpsexpave = (fpsexpave * 0.90) + (framefps * 0.10);

      printf("OSPRay2Renderer) %c AA:%2d AO:%2d, %4d tot RT FPS: %.1f  %.4f s/frame sf: %d  \r",
             statestr[state], cur_aa_samples, cur_ao_samples, 
             totalsamplecount, fpsexpave, frametime, subframe_count);

      fflush(stdout);
      state = (state+1) & 3;
    }

  } // end of per-cycle event processing

  printf("\n");

  // write the output image upon exit...
  if (lasterror == 0 /* XXX SUCCESS */) {
    wkf_timer_start(osp_timer);
    // write output image
    const unsigned char *FB = (const unsigned char*)ospMapFrameBuffer(ospFrameBuffer, OSP_FB_COLOR);
    if (write_image_file_rgb4u(filename, FB, width, height)) {
      printf("OSPRay2Renderer) Failed to write output image!\n");
    }
    ospUnmapFrameBuffer(FB, ospFrameBuffer);
    wkf_timer_stop(osp_timer);

    if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
      printf("OSPRay2Renderer) image file I/O time: %f secs\n", wkf_timer_time(osp_timer));
    }
  }

  glwin_destroy(win);
}

#endif


void OSPRay2Renderer::render_to_file(const char *filename) {
  DBG();
  if (!context_created)
    return;

  // Unless overridden by environment variables, we use the incoming
  // window size parameters from VMD to initialize the RT image dimensions.
  int wsx=width, wsy=height;
  const char *imageszstr = getenv("VMDOSPRAYIMAGESIZE");
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
  double starttime = wkf_timer_timenow(osp_timer);

  //
  // run the renderer 
  //
  if (lasterror == 0 /* XXX SUCCESS */) {
    // clear the accumulation buffer
    ospResetAccumulation(ospFrameBuffer);

    // Render to the accumulation buffer for the required number of passes
    if (getenv("VMDOSPRAYNORENDER") == NULL) {
      int accum_sample;
      for (accum_sample=0; accum_sample<ext_aa_loops; accum_sample++) {
        // The accumulation subframe count must be updated to ensure that
        // any custom RNGs for AA and AO get correctly re-seeded
#if 1
        OSPFuture fut;
        fut = ospRenderFrame(ospFrameBuffer, ospRenderer, ospCamera, ospWorld);
        ospWait(fut, OSP_FRAME_FINISHED);
        ospRelease(fut);
#else
        ospRenderFrameBlocking(ospFrameBuffer, ospRenderer, ospCamera, ospWorld);
#endif
      }
    }

    // copy the accumulation buffer image data to the framebuffer and perform
    // type conversion and normaliztion on the image data...
    double rtendtime = wkf_timer_timenow(osp_timer);
    time_ray_tracing = rtendtime - starttime;

    if (lasterror == 0 /* XXX SUCCESS */) {
      // write output image to a file unless we are benchmarking
      if (getenv("VMDOSPRAYNOSAVE") == NULL) {
        const unsigned char *FB = (const unsigned char*)ospMapFrameBuffer(ospFrameBuffer, OSP_FB_COLOR);
        if (write_image_file_rgb4u(filename, FB, width, height)) {
          printf("OSPRay2Renderer) Failed to write output image!\n");
        }
        ospUnmapFrameBuffer(FB, ospFrameBuffer);
      }
      time_image_io = wkf_timer_timenow(osp_timer) - rtendtime;
    } else {
      printf("OSPRay2Renderer) Error during rendering.  Rendering aborted.\n");
    }

    if (verbose == RT_VERB_TIMING || verbose == RT_VERB_DEBUG) {
      printf("OSPRay2Renderer) ctx setup %.2f  valid %.2f  AS %.2f  RT %.2f io %.2f\n", time_ctx_setup, time_ctx_validate, time_ctx_AS_build, time_ray_tracing, time_image_io);
    }
  } else {
    printf("OSPRay2Renderer) Error during AS generation.  Rendering aborted.\n");
  }
}


void OSPRay2Renderer::add_material(int matindex,
                                   float ambient, float diffuse, 
                                   float specular,
                                   float shininess, float reflectivity,
                                   float opacity, 
                                   float outline, float outlinewidth,
                                   int transmode) {
  int oldmatcount = materialcache.num();
  if (oldmatcount <= matindex) {
    osp_material m;
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
    if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) Adding material[%d]\n", matindex);

    materialcache[matindex].ambient      = ambient;
    materialcache[matindex].diffuse      = diffuse; 
    materialcache[matindex].specular     = specular;
    materialcache[matindex].shininess    = shininess;
    materialcache[matindex].reflectivity = reflectivity;
    materialcache[matindex].opacity      = opacity;
    materialcache[matindex].outline      = outline;
    materialcache[matindex].outlinewidth = outlinewidth;
    materialcache[matindex].transmode    = transmode;

    // create an OSPRay material object too...
    const char *rstr = (osp_rendermode == RT_SCIVIS) ? "scivis" : "pathtracer";

    OSPMaterial ospMat;
    if (ambient > 0.15f && osp_rendermode == RT_PATHTRACER) {
      float lumvalue = materialcache[matindex].ambient;
      lumvalue *= 8.0f; 

printf("***\n***\n*** LUMINOUS MATERIAL: %f\n***\n***\n", lumvalue);
      // 
      // Luminous material type
      //
      ospMat = ospNewMaterial(rstr, "luminous");
      float mtmp[3] = {1.0f, 1.0f, 1.0f};
      ospSetParam(ospMat, "color", OSP_VEC3F, mtmp);
      ospSetParam(ospMat, "intensity", OSP_FLOAT, &lumvalue);
    } else {

#if 1
      //
      // Path traced material type
      //
      ospMat = ospNewMaterial(rstr, "principled");
      float mtmp[3];

#if 0
      if (osp_rendermode == RT_PATHTRACER) {
        mtmp[0] = mtmp[1] = mtmp[2] = materialcache[matindex].ambient;
        ospSetParam(ospMat, "ka", OSP_VEC3F, mtmp);
      }
#endif

      mtmp[0] = mtmp[1] = mtmp[2] = 1.0f;
      ospSetParam(ospMat, "baseColor", OSP_VEC3F, mtmp);
  
      ospSetParam(ospMat, "diffuse", OSP_FLOAT, &materialcache[matindex].diffuse);
      ospSetParam(ospMat, "specular", OSP_FLOAT, &materialcache[matindex].specular);
//      ospSetParam(ospMat, "ns", OSP_FLOAT, &materialcache[matindex].shininess);
//      ospSetParam(ospMat, "transmission", OSP_FLOAT, &materialcache[matindex].opacity);
      ospSetParam(ospMat, "opacity", OSP_FLOAT, &materialcache[matindex].opacity);

#else

      //
      // Simple sci-vis OBJ material type
      //
      ospMat = ospNewMaterial(rstr, "obj");
      float mtmp[3];

#if 0
      if (osp_rendermode == RT_PATHTRACER) {
        mtmp[0] = mtmp[1] = mtmp[2] = materialcache[matindex].ambient;
        ospSetParam(ospMat, "ka", OSP_VEC3F, mtmp);
      }
#endif

      mtmp[0] = mtmp[1] = mtmp[2] = materialcache[matindex].diffuse;
      ospSetParam(ospMat, "kd", OSP_VEC3F, mtmp);
  
      if (osp_rendermode == RT_PATHTRACER) {
        mtmp[0] = mtmp[1] = mtmp[2] = materialcache[matindex].specular;
        ospSetParam(ospMat, "ks", OSP_VEC3F, mtmp);
      }

      ospSetParam(ospMat, "d", OSP_FLOAT, &materialcache[matindex].opacity);

      if (osp_rendermode == RT_PATHTRACER) {
        ospSetParam(ospMat, "ns", OSP_FLOAT, &materialcache[matindex].shininess);
      }
#endif
    }


    /// XXX The OSPRay path tracer supports filtered transparency 
    ///     with a "Tf" material value, but there are noteworthy
    ///     restrictions about energy conservation etc to worry about:
    /// https://www.ospray.org/documentation.html#materials

    ospCommit(ospMat);
    materialcache[matindex].mat = ospMat;

    materialcache[matindex].isvalid      = 1;
  }
}


void OSPRay2Renderer::init_materials() {
  DBG();
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer: init_materials()\n");

}


void OSPRay2Renderer::set_material(OSPGeometricModel &model, int matindex, float *uniform_color) {
  if (!context_created)
    return;

  if (verbose == RT_VERB_DEBUG) 
    printf("OSPRay2Renderer)   setting material %d\n", matindex);
  ospSetParam(model, "material", OSP_MATERIAL, &materialcache[matindex].mat);
}


void OSPRay2Renderer::attach_mesh(int numverts, int numfacets, int matindex,
                                  osp_trimesh_v3f_n3f_c3f &mesh) {
  mesh.matindex = matindex;
  mesh.verts = ospNewSharedData(mesh.v, OSP_VEC3F, numverts, 0);
  ospCommit(mesh.verts);
  mesh.norms = ospNewSharedData(mesh.n, OSP_VEC3F, numverts, 0);
  ospCommit(mesh.norms);
  mesh.cols  = ospNewSharedData(mesh.c, OSP_VEC4F, numverts, 0);
  ospCommit(mesh.cols);
  mesh.ind   = ospNewSharedData(mesh.f, OSP_VEC3UI, numfacets, 0);
  ospCommit(mesh.ind);

  mesh.geom  = ospNewGeometry("mesh");
  ospSetParam(mesh.geom, "vertex.position", OSP_DATA, &mesh.verts);
  ospSetParam(mesh.geom, "vertex.normal",   OSP_DATA, &mesh.norms);
  ospSetParam(mesh.geom, "index",           OSP_DATA, &mesh.ind);
  ospSetParam(mesh.geom, "vertex.color",    OSP_DATA, &mesh.cols);
  ospCommit(mesh.geom);
  ospRelease(mesh.verts);
  ospRelease(mesh.norms);
  ospRelease(mesh.cols);
  ospRelease(mesh.ind);

  mesh.model = ospNewGeometricModel(mesh.geom);
  set_material(mesh.model, matindex, NULL);
  ospCommit(mesh.model);
  ospRelease(mesh.geom);
  trimesh_v3f_n3f_c3f.append(mesh); 
}


void OSPRay2Renderer::attach_sphere_array(int numsp, int matindex,
                                          osp_sphere_array_color &sparray) {
  sparray.matindex = matindex;
  sparray.cents = ospNewSharedData(sparray.xyz, OSP_VEC3F, numsp, 0);
  ospCommit(sparray.cents);
  sparray.rads = ospNewSharedData(sparray.radii, OSP_FLOAT, numsp, 0);
  ospCommit(sparray.rads);
  sparray.cols = ospNewSharedData(sparray.colors, OSP_VEC4F, numsp, 0);
  ospCommit(sparray.cols);

  sparray.geom  = ospNewGeometry("sphere");
  ospSetParam(sparray.geom, "sphere.position", OSP_DATA, &sparray.cents);
  ospSetParam(sparray.geom, "sphere.radius",   OSP_DATA, &sparray.rads);
  ospCommit(sparray.geom);
  ospRelease(sparray.cents);
  ospRelease(sparray.rads);

  sparray.model = ospNewGeometricModel(sparray.geom);
  ospSetParam(sparray.model, "color",    OSP_DATA, &sparray.cols);
  set_material(sparray.model, matindex, NULL);
  ospCommit(sparray.model);
  ospRelease(sparray.geom);
  ospRelease(sparray.cols);

  spheres_color.append(sparray);
}


void OSPRay2Renderer::attach_cylinder_array(int numcyl, int matindex,
                                            osp_cylinder_array_color &cylarray) {
  cylarray.matindex = matindex;
  cylarray.cyls = ospNewSharedData(cylarray.vertsrads, OSP_VEC4F, numcyl * 2, 0);
  ospCommit(cylarray.cyls);
  cylarray.cols = ospNewSharedData(cylarray.colors, OSP_VEC4F, numcyl, 0);
  ospCommit(cylarray.cols);
  cylarray.ind  = ospNewSharedData(cylarray.indices, OSP_UINT, numcyl, 0);
  ospCommit(cylarray.ind);

  cylarray.geom  = ospNewGeometry("curve");
  ospSetParam(cylarray.geom, "vertex.position_radius", OSP_DATA, &cylarray.cyls);
  ospSetParam(cylarray.geom, "index",    OSP_DATA, &cylarray.ind);
  unsigned char type = OSP_DISJOINT;
  ospSetParam(cylarray.geom, "type",     OSP_UCHAR, &type);
  unsigned char basis = OSP_LINEAR;
  ospSetParam(cylarray.geom, "basis",    OSP_UCHAR, &basis);
  ospCommit(cylarray.geom);
  ospRelease(cylarray.cyls);
  ospRelease(cylarray.ind);

  cylarray.model = ospNewGeometricModel(cylarray.geom);
  ospSetParam(cylarray.model, "color",   OSP_DATA, &cylarray.cols);
  set_material(cylarray.model, matindex, NULL);
  ospCommit(cylarray.model);
  ospRelease(cylarray.geom);
  ospRelease(cylarray.cols);

  cylinders_color.append(cylarray);
}



void OSPRay2Renderer::add_directional_light(const float *dir, const float *color) {
  DBG();
  osp_directional_light l;
  vec_copy(l.dir, dir);
  vec_copy(l.color, color);

  directional_lights.append(l);
}


void OSPRay2Renderer::add_positional_light(const float *pos, const float *color) {
  DBG();
  osp_positional_light l;
  vec_copy(l.pos, pos);
  vec_copy(l.color, color);

  positional_lights.append(l);
}


void OSPRay2Renderer::cylinder_array(Matrix4 *wtrans, float radius,
                                   float *uniform_color,
                                   int cylnum, float *points, int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating cylinder array: %d...\n", cylnum);

  cylinder_array_cnt += cylnum;
  
  osp_cylinder_array_color ca;
  memset(&ca, 0, sizeof(ca));
  ca.num = cylnum;
  ca.vertsrads = (float *) calloc(1, cylnum * 8 * sizeof(float));
  ca.colors    = (float *) calloc(1, cylnum * 4 * sizeof(float));
  ca.indices = (unsigned int *) calloc(1, cylnum * 1 * sizeof(unsigned int));

  int i,ind4,ind6,ind8;
  if (wtrans == NULL) {
    for (i=0,ind4,ind6=0,ind8=0; i<cylnum; i++,ind4+=4,ind6+=6,ind8+=8) {
      vec_copy(&ca.vertsrads[ind8  ], &points[ind6  ]);
      ca.vertsrads[ind8+3] = radius;
      vec_copy(&ca.vertsrads[ind8+4], &points[ind6+3]);
      ca.vertsrads[ind8+7] = radius;

      vec_copy(&ca.colors[ind4], &uniform_color[0]);
      ca.colors[ind4+3] = 1.0f;

      ca.indices[i] = i*2;
    }
  } else {
    for (i=0,ind4=0,ind6=0,ind8=0; i<cylnum; i++,ind4+=4,ind6+=6,ind8+=8) {
      // apply transforms on points, radii
      wtrans->multpoint3d(&points[ind6  ], &ca.vertsrads[ind8  ]);
      ca.vertsrads[ind8+3] = radius;
      wtrans->multpoint3d(&points[ind6+3], &ca.vertsrads[ind8+4]);
      ca.vertsrads[ind8+7] = radius;

      vec_copy(&ca.colors[ind4], &uniform_color[0]);
      ca.colors[ind4+3] = 1.0f;

      ca.indices[i] = i*2;
    }
  }

  attach_cylinder_array(cylnum, matindex, ca);
}


void OSPRay2Renderer::cylinder_array_color(Matrix4 & wtrans, float rscale,
                                         int cylnum, float *points,
                                         float *radii, float *colors,
                                         int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating color cylinder array: %d...\n", cylnum);
  cylinder_array_color_cnt += cylnum;

  osp_cylinder_array_color cac;
  memset(&cac, 0, sizeof(cac));
  cac.num = cylnum;
  cac.vertsrads = (float *) calloc(1, cylnum * 8 * sizeof(float));
  cac.colors    = (float *) calloc(1, cylnum * 4 * sizeof(float));
  cac.indices = (unsigned int *) calloc(1, cylnum * 1 * sizeof(unsigned int));

  int i, ind3, ind4, ind6, ind8;
  for (i=0,ind3=0,ind4=0,ind6=0,ind8=0; i<cylnum; i++,ind3+=3,ind4+=4,ind6+=6,ind8+=8) {
    // apply transforms on points, radii
    wtrans.multpoint3d(&points[ind6  ], &cac.vertsrads[ind8  ]);
    cac.vertsrads[ind8+3] = radii[i] * rscale; // radius
    wtrans.multpoint3d(&points[ind6+3], &cac.vertsrads[ind8+4]);
    cac.vertsrads[ind8+7] = radii[i] * rscale; // radius

    vec_copy(&cac.colors[ind4], &colors[ind3]);
    cac.colors[ind4+3] = 1.0f;

    cac.indices[i] = i*2;
  }

  attach_cylinder_array(cylnum, matindex, cac);
}

#if 0
void OSPRay2Renderer::ring_array_color(Matrix4 & wtrans, float rscale,
                                     int rnum, float *centers,
                                     float *norms, float *radii, 
                                     float *colors, int matindex) {
}
#endif


void OSPRay2Renderer::sphere_array(Matrix4 *wtrans, float rscale,
                                 float *uniform_color,
                                 int numsp, float *centers,
                                 float *radii,
                                 int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating sphere array: %d...\n", numsp);
  sphere_array_cnt += numsp;

  const rgba c = { uniform_color[0], uniform_color[1], uniform_color[2], 1.0f};

  osp_sphere_array_color sp;
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


void OSPRay2Renderer::sphere_array_color(Matrix4 & wtrans, float rscale,
                                       int numsp, float *centers,
                                       float *radii, float *colors,
                                       int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating sphere array color: %d...\n", numsp);
  sphere_array_color_cnt += numsp;

  osp_sphere_array_color sp;
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


void OSPRay2Renderer::tricolor_list(Matrix4 & wtrans, int numtris, float *vnc,
                                  int matindex) {
  if (!context_created) return;
//if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating tricolor list: %d...\n", numtris);
  tricolor_cnt += numtris;

  // create and fill the OSPRay trimesh memory buffer
  osp_trimesh_v3f_n3f_c3f mesh;
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


void OSPRay2Renderer::trimesh_c4n3v3(Matrix4 & wtrans, int numverts,
                                   float *cnv, int numfacets, int * facets,
                                   int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating trimesh_c4n3v3: %d...\n", numfacets);
  trimesh_c4u_n3b_v3f_cnt += numfacets;

  // create and fill the OSPRay trimesh memory buffer
  osp_trimesh_v3f_n3f_c3f mesh;
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
// to the best that OSPRay allows
//
void OSPRay2Renderer::trimesh_c4u_n3b_v3f(Matrix4 & wtrans, unsigned char *c, 
                                         signed char *n, float *v, 
                                         int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating trimesh_c4u_n3b_v3f: %d...\n", numfacets);
  trimesh_n3b_v3f_cnt += numfacets;

  // create and fill the OSPRay trimesh memory buffer
  osp_trimesh_v3f_n3f_c3f mesh;
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



void OSPRay2Renderer::trimesh_c4u_n3f_v3f(Matrix4 & wtrans, unsigned char *c, 
                                         float *n, float *v, 
                                         int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating trimesh_c4u_n3f_v3f: %d...\n", numfacets);
  tricolor_cnt += numfacets;

  // create and fill the OSPRay trimesh memory buffer
  osp_trimesh_v3f_n3f_c3f mesh;
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


void OSPRay2Renderer::trimesh_n3b_v3f(Matrix4 & wtrans, float *uniform_color, 
                                     signed char *n, float *v, 
                                     int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating trimesh_n3b_v3f: %d...\n", numfacets);
  trimesh_n3b_v3f_cnt += numfacets;

  // create and fill the OSPRay trimesh memory buffer
  osp_trimesh_v3f_n3f_c3f mesh;
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
//     but that should go away as soon as OSPRay allows it.
void OSPRay2Renderer::trimesh_n3f_v3f(Matrix4 & wtrans, float *uniform_color, 
                                    float *n, float *v, int numfacets, 
                                    int matindex) {
  DBG();
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating trimesh_n3f_v3f: %d...\n", numfacets);
  trimesh_n3f_v3f_cnt += numfacets;
  // create and fill the OSPRay trimesh memory buffer
  osp_trimesh_v3f_n3f_c3f mesh;
  memset(&mesh, 0, sizeof(mesh));
  mesh.num = numfacets;
  mesh.v = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.n = (float *) calloc(1, numfacets * 9*sizeof(float));
  mesh.c = (float *) calloc(1, numfacets * 12*sizeof(float));
  mesh.f = (int *) calloc(1, numfacets * 3*sizeof(int));

  float alpha = 1.0f;

  // create and fill the OSPRay trimesh memory buffer
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
void OSPRay2Renderer::trimesh_v3f(Matrix4 & wtrans, float *uniform_color, 
                                float *v, int numfacets, int matindex) {
  if (!context_created) return;
  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating trimesh_v3f: %d...\n", numfacets);
  trimesh_v3f_cnt += numfacets;

  set_material(geom, matindex, NULL);
  append_objects(buf, geom, instance);
}

#endif



void OSPRay2Renderer::tristrip(Matrix4 & wtrans, int numverts, const float * cnv,
                             int numstrips, const int *vertsperstrip,
                             const int *facets, int matindex) {
  if (!context_created) return;
  int i;
  int numfacets = 0;
  for (i=0; i<numstrips; i++) 
    numfacets += (vertsperstrip[i] - 2);  

  if (verbose == RT_VERB_DEBUG) printf("OSPRay2Renderer) creating tristrip: %d...\n", numfacets);
  tricolor_cnt += numfacets;

  // create and fill the OSPRay trimesh memory buffer
  osp_trimesh_v3f_n3f_c3f mesh;
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



