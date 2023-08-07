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
 *	$RCSfile: Scene.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.100 $	$Date: 2021/12/21 22:56:37 $
 *
 ***************************************************************************/
/**
 *  \file Scene.C
 *  \brief Scene maintains a list of Displayable objects and draws them 
 *         to the active DisplayDevice.
 *
 * The Scene has a list of Displayable objects and display commands.
 * The command lists are used to draw the objects, the Displayable
 * objects to prepare and update objects for drawing.
 */

#include "Scene.h"
#include "DisplayDevice.h"
#include "Inform.h"
#include "DispCmds.h"
#include "utilities.h"
#include "FileRenderList.h"
#include "FileRenderer.h"
#include "ColorScaleTables.h"

static const int num_scalemethods = 36;
static const ColorScale defScales[] = {
 { {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 1.0}, "RWB"},
 { {0.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 0.0, 0.0}, "BWR"},
 { {1.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}, "RGryB"},
 { {0.0, 0.0, 1.0}, {0.5, 0.5, 0.5}, {1.0, 0.0, 0.0}, "BGryR"},
 { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, "RGB"},
 { {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, "BGR"},
 { {1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 0.0}, "RWG"},
 { {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, {1.0, 0.0, 0.0}, "GWR"},
 { {0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 1.0}, "GWB"},
 { {0.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 0.0}, "BWG"},
 { {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {1.0, 1.0, 1.0}, "BlkW"},
 { {1.0, 1.0, 1.0}, {0.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, "WBlk"},

 // sequential
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "cividis"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "viridis"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "magma"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "plasma"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "inferno"},

 // sequential-linear
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L3"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L8"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L9"}, 
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L16"}, 
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L17"}, 
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L18"}, 
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L19"}, 
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_L20"}, 

 // cyclic
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_C2"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_C4"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_C6"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_C7"},

 // isoluminant
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_I1"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_I2"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_I3"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_D11"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_D12"},

 // rainbow
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "turbo"},
 { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, "CET_R2"} 
};

static const char *colScaleMenuNames[] = {
 "RWB",
 "BWR",
 "RGryB",
 "BGryR",
 "RGB",
 "BGR",
 "RWG",
 "GWR",
 "GWB",
 "BWG",
 "BlkW",
 "WBlk",

 "Sequential/cividis",
 "Sequential/viridis",
 "Sequential/magma",
 "Sequential/plasma",
 "Sequential/inferno",

 "Linear/CET_L3",
 "Linear/CET_L8",
 "Linear/CET_L9",
 "Linear/CET_L16",
 "Linear/CET_L17",
 "Linear/CET_L18",
 "Linear/CET_L19",
 "Linear/CET_L20",

 "Cyclic/CET_C2",
 "Cyclic/CET_C4",
 "Cyclic/CET_C6",
 "Cyclic/CET_C7",

 "Isoluminant/CET_I1",
 "Isoluminant/CET_I2",
 "Isoluminant/CET_I3",
 "Isoluminant/CET_D11",
 "Isoluminant/CET_D12",

 "Rainbow/Turbo",
 "Rainbow/CET_R2"
};


enum colorScaleType {
  SCALE_TYPE_TABULATED=0,
  SCALE_TYPE_TABULATED_SRGB=1,
  SCALE_TYPE_DIVERGENT=2
};


/// color gradient/ramp used by value-based coloring methods
class ColorScaleTypeInfo {
public:
  enum colorScaleType type;
  const float *colortable;
  int colortabsz;

#if 0
  int operator==(const ColorScale c) {
    return !memcmp(&c, this, sizeof(ColorScale));
  }
#endif
};


// hard-coded color scale types
static const ColorScaleTypeInfo scaleTypes[] = {
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "RWB"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "BWR"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "RGryB"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "BGryR"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "RGB"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "BGR"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "RWG"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "GWR"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "GWB"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "BWG"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "BlkW"
  {SCALE_TYPE_DIVERGENT,            NULL,   0}, // "WBlk"

  {SCALE_TYPE_TABULATED,    cividis_data, 255}, // "cividis"
  {SCALE_TYPE_TABULATED,    viridis_data, 255}, // "viridis"
  {SCALE_TYPE_TABULATED,      magma_data, 255}, // "plasma"
  {SCALE_TYPE_TABULATED,     plasma_data, 255}, // "inferno"
  {SCALE_TYPE_TABULATED,    inferno_data, 255}, // "magma"

  {SCALE_TYPE_TABULATED,     CET_L3_data, 255}, // linear
  {SCALE_TYPE_TABULATED,     CET_L8_data, 255}, // linear
  {SCALE_TYPE_TABULATED,     CET_L9_data, 255}, // linear
  {SCALE_TYPE_TABULATED,    CET_L16_data, 255}, // linear
  {SCALE_TYPE_TABULATED,    CET_L17_data, 255}, // linear
  {SCALE_TYPE_TABULATED,    CET_L18_data, 255}, // linear
  {SCALE_TYPE_TABULATED,    CET_L19_data, 255}, // linear
  {SCALE_TYPE_TABULATED,    CET_L20_data, 255}, // linear

  {SCALE_TYPE_TABULATED,     CET_C2_data, 255}, // cyclic
  {SCALE_TYPE_TABULATED,     CET_C4_data, 255}, // cyclic 
  {SCALE_TYPE_TABULATED,     CET_C6_data, 255}, // cyclic 
  {SCALE_TYPE_TABULATED,     CET_C7_data, 255}, // cyclic 

  {SCALE_TYPE_TABULATED,     CET_I1_data, 255}, // isoluminant
  {SCALE_TYPE_TABULATED,     CET_I2_data, 255}, // isoluminant
  {SCALE_TYPE_TABULATED,     CET_I3_data, 255}, // isoluminant
  {SCALE_TYPE_TABULATED,    CET_D11_data, 255}, // isoluminant
  {SCALE_TYPE_TABULATED,    CET_D12_data, 255}, // isoluminant

  {SCALE_TYPE_TABULATED, turbo_srgb_data, 255}, // rainbow
  {SCALE_TYPE_TABULATED,     CET_R2_data, 255}  // rainbow
};


// maybe VMD should be using a color wheel scheme
// to make it easier to pick contrasting and complementary 
// colors.
static const char *defColorNames[REGCLRS] = {
  "blue",    "red",    "gray",    "orange",  
  "yellow",  "tan",    "silver",  "green",  
  "white",   "pink",   "cyan",   "purple",
  "lime",    "mauve",  "ochre",  "iceblue", 

  // XXX black is expected to be at the end of the list
  "black"

#if (REGCLRS > 17)
 ,"yellow2",  "yellow3",   "green2",    "green3",
  "cyan2",    "cyan3",     "blue2",     "blue3",
  "violet",   "violet2",   "magenta",   "magenta2",  
  "red2",     "red3",      "orange2",   "orange3"
#endif

};

const float Scene::defaultColor[] = {
   0.0f,   0.0f, 1.00f,  1.0f,  0.0f,  0.0f, // BLUE, RED
   0.35f, 0.35f, 0.35f,  1.0f,  0.5f,  0.0f, // GREY, ORANGE
   1.0f,   1.0f,  0.0f,  0.5f,  0.5f,  0.2f, // YELLOW, TAN
   0.6f,   0.6f,  0.6f,  0.0f,  1.0f,  0.0f, // SILVER, GREEN
   1.0f,   1.0f,  1.0f,  1.0f,  0.6f,  0.6f, // WHITE, PINK
   0.25f, 0.75f, 0.75f, 0.65f,  0.0f, 0.65f, // CYAN, PURPLE
   0.5f,   0.9f,  0.4f,  0.9f,  0.4f,  0.7f, // LIME, MAUVE
   0.5f,   0.3f,  0.0f,  0.5f,  0.5f, 0.75f, // OCHRE, ICEBLUE

   // XXX black is expected to be at the end of the list
   0.0f,   0.0f,  0.0f                       // BLACK
  
#if (REGCLRS > 17)
  ,0.88f, 0.97f, 0.02f,  0.55f, 0.90f, 0.02f, // yellow
   0.00f, 0.90f, 0.04f,  0.00f, 0.90f, 0.50f, // green
   0.00f, 0.88f, 1.00f,  0.00f, 0.76f, 1.00f, // cyan
   0.02f, 0.38f, 0.67f,  0.01f, 0.04f, 0.93f, // blue
   0.27f, 0.00f, 0.98f,  0.45f, 0.00f, 0.90f, // violet
   0.90f, 0.00f, 0.90f,  1.00f, 0.00f, 0.66f, // magenta
   0.98f, 0.00f, 0.23f,  0.81f, 0.00f, 0.00f, // red
   0.89f, 0.35f, 0.00f,  0.96f, 0.72f, 0.00f  // orange
#endif

};


/// Displayable subclass containing the background color information
class DisplayColor : public Displayable {
private:
  int colorCat;
  int dcindex;
  int dccolor;
  int changed;
 
protected:
  void do_color_changed(int cat) {
    if (cat == colorCat) {
      dccolor = scene->category_item_value(colorCat, dcindex);
      changed = 1;
    }
  } 
  void do_color_rgb_changed(int color) {
    if (color == dccolor) 
      changed = 1;
  }
   
public:
  DisplayColor(Displayable *d, const char *coloritemname, int colorindex)
  : Displayable(d), changed(0) {
    // Query the index of the Display color category, or add it if necessary
    colorCat = scene->category_index("Display");
    if (colorCat == -1) {
      colorCat = scene->add_color_category("Display");
    } 
    dcindex = scene->add_color_item(colorCat, coloritemname, colorindex);
    do_color_changed(colorCat);
  }
  int color_changed() const { return changed; }
  void clear_changed() { changed = 0; }
  int color_id() const { return dccolor; }
};


///  constructor 
Scene::Scene() : root(this) {
  set_background_mode(0);
  reset_lights();

  // initialize color names and values
  int i;
  for (i=BEGREGCLRS; i<REGCLRS; i++) {
    colorNames.add_name(defColorNames[i], i);
    set_color_value(i, defaultColor + 3L*i);
  }

  // initialize color scale parameter info
  for (i=0; i<num_scalemethods; i++) {
    colorScales.append(defScales[i]);
    colorScaleMenuNames.append(colScaleMenuNames[i]);
  }

  scaleActive = 0;
  scaleMin = 0.06f; // default color scale "offset" value
  scaleMid = 0.50f; // default color scale "midpoint" value
  scaleMax = 1.00f;
  scalePosterize = 0;
  scaleReverse = 0;
  create_colorscale();

  // background color 
  background = new DisplayColor(&root, "Background", REGBLACK);
  background_color_changed = 0;
  background_color_id = 0;

  // background gradient colors
  backgradtop = new DisplayColor(&root, "BackgroundTop", REGBLACK);
  backgradtop_color_changed = 0;
  backgradtop_color_id = 0;
  backgradbot = new DisplayColor(&root, "BackgroundBot", REGBLUE2);
  backgradbot_color_changed = 0;
  backgradbot_color_id = 0;

  // foreground color
  foreground = new DisplayColor(&root, "Foreground", REGWHITE);
  foreground_color_changed = 0;
  foreground_color_id = 0;
}

///////////////////////////  destructor
Scene::~Scene(void) {
  for (int i=0; i<categories.num(); i++) delete categories.data(i);
}

///////////////////////////  public routines 

// Set background drawing mode 
void Scene::set_background_mode(int mode) {
  backgroundmode = mode;
  backgroundmode_changed = 1;
}

// Query background drawing mode 
int Scene::background_mode(void) {
  return backgroundmode;
}


//
// routines to deal with light sources
//
  
// default light data
static const float def_light_color[3] = { 1.0, 1.0, 1.0 };
static const float def_light_pos[DISP_LIGHTS][3] = {
        { -0.1f, 0.1f,  1.0f }, {  1.0f,  2.0f, 0.5f },
        { -1.0f, 2.0f, -1.0f }, { -1.0f, -1.0f, 0.0f }
};
static const float def_adv_light_pos[DISP_LIGHTS][3] = {
        { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, 
        { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }
};

static const int def_light_on[DISP_LIGHTS] = { TRUE, TRUE, FALSE, FALSE };

void Scene::define_light(int n, const float *color, const float *position) {
  if (n < 0 || n >= DISP_LIGHTS)
    return;
  
  for (int i=0; i < 3; i++) {
    lightState[n].color[i] = color[i];
    lightState[n].pos[i] = position[i];
  }
  light_changed = 1;
}

void Scene::activate_light(int n, int turnon) {
  if (n < 0 || n >= DISP_LIGHTS )
    return;
  lightState[n].on = turnon;
  light_changed = 1;
}

void Scene::rotate_light(int n, float theta, char axis) {
  if (n < 0 || n >= DISP_LIGHTS)
    return;
  Matrix4 mat;
  mat.rot(theta,axis);
  mat.multpoint3d(lightState[n].pos,lightState[n].pos);
  light_changed = 1;
}

void Scene::move_light(int n, const float *p) {
  if (n < 0 || n >= DISP_LIGHTS) return;
  for (int i=0; i<3; i++) lightState[n].pos[i] = p[i];
  light_changed = 1;
}

const float *Scene::light_pos(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return lightState[n].pos;
}

const float *Scene::light_pos_default(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return def_light_pos[n];
}

const float *Scene::light_color(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return lightState[n].color;
}

const float *Scene::light_color_default(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return def_light_color;
}


void Scene::define_adv_light(int n, const float *color, 
                             const float *position,
                             float constant, float linear, float quad,
                             float *spotdir, 
                             float fallstart, float fallend, int spoton) {
  if (n < 0 || n >= DISP_LIGHTS)
    return;
  
  for (int i=0; i < 3; i++) {
    advLightState[n].color[i] = color[i];
    advLightState[n].pos[i] = position[i];
    advLightState[n].spotdir[i] = spotdir[i];
  }
  advLightState[n].constfactor = constant;
  advLightState[n].linearfactor = linear;
  advLightState[n].quadfactor = quad;
  advLightState[n].fallstart = fallstart;
  advLightState[n].fallend = fallend;
  advLightState[n].spoton = spoton;
  adv_light_changed = 1;
}

void Scene::activate_adv_light(int n, int turnon) {
  if (n < 0 || n >= DISP_LIGHTS )
    return;
  advLightState[n].on = turnon;
  adv_light_changed = 1;
}

void Scene::move_adv_light(int n, const float *p) {
  if (n < 0 || n >= DISP_LIGHTS) return;
  for (int i=0; i<3; i++) 
    advLightState[n].pos[i] = p[i];
  adv_light_changed = 1;
}

const float *Scene::adv_light_pos(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return advLightState[n].pos;
}

const float *Scene::adv_light_pos_default(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return def_adv_light_pos[n];
}

const float *Scene::adv_light_color(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return advLightState[n].color;
}

const float *Scene::adv_light_color_default(int n) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  return def_light_color;
}

void Scene::adv_light_attenuation(int n, float constant, float linear,
                                  float quad) {
  if (n < 0 || n >= DISP_LIGHTS) return;
  advLightState[n].constfactor = constant;
  advLightState[n].linearfactor = linear;
  advLightState[n].quadfactor = quad;
}

void Scene::adv_light_get_attenuation(int n, float &constant, float &linear,
                                      float &quad) const {
  if (n < 0 || n >= DISP_LIGHTS) return;
  constant = advLightState[n].constfactor;
  linear   = advLightState[n].linearfactor;
  quad     = advLightState[n].quadfactor;
}

void Scene::adv_light_spotlight(int n, float *spotdir, float fallstart, 
                                float fallend, int spoton) {
  if (n < 0 || n >= DISP_LIGHTS) return;
  advLightState[n].fallstart = fallstart;
  advLightState[n].fallend = fallend;
  advLightState[n].spoton = spoton;
  for (int i=0; i<3; i++) 
    advLightState[n].spotdir[i] = spotdir[i];
}

const float *Scene::adv_light_get_spotlight(int n, float &fallstart, 
                                            float &fallend, int &spoton) const {
  if (n < 0 || n >= DISP_LIGHTS) return NULL;
  fallstart = advLightState[n].fallstart;
  fallend   = advLightState[n].fallend;
  spoton    = advLightState[n].spoton;
  return advLightState[n].spotdir;
}

void Scene::reset_lights(void) {
  int i;

  // standard directional lights
  for (i=0; i<DISP_LIGHTS; i++) { 
    define_light(i, def_light_color, def_light_pos[i]);
    activate_light(i, def_light_on[i]);
  }
  light_changed = 1;

  // advanced lights
  for (i=0; i<DISP_LIGHTS; i++) { 
    float spotdir[] = { 0.0f, 0.0f, 1.0f };
    define_adv_light(i, def_light_color, def_light_pos[i],
                     1.0f, 0.0f, 0.0f,
                     spotdir, 0.3f, 0.7f, 0);
    activate_adv_light(i, 0);
  }
  adv_light_changed = 1;
}


// prepare all registered Displayables
int Scene::prepare() {
  background_color_changed = background->color_changed();
  background_color_id = background->color_id();
  background->clear_changed();

  backgradtop_color_changed = backgradtop->color_changed();
  backgradtop_color_id = backgradtop->color_id();
  backgradtop->clear_changed();

  backgradbot_color_changed = backgradbot->color_changed();
  backgradbot_color_id = backgradbot->color_id();
  backgradbot->clear_changed();

  foreground_color_changed = foreground->color_changed();
  foreground_color_id = foreground->color_id();
  foreground->clear_changed();

  return root.draw_prepare() || backgroundmode_changed || light_changed || 
         background_color_changed || 
         backgradtop_color_changed || backgradbot_color_changed || 
         foreground_color_changed;
}

// draws the scene to the given DisplayDevice
// this is the only Scene which tells the display to do graphics commands
//
// XXX implementation of multi-pass rendering, accumulation buffers,
//     and better snapshot behavior will require modifications to this
//     mechanism.  The current system is very simple minded, and doesn't
//     allow for property-based sorting of geometry, occlusion culling 
//     algorithms, etc.  This needs to be changed.
//
// XXX note, this method should really be a 'const' method since
//     it is run concurrently by several processes that share memory, but
//     we can't actually write it that way since the locking routines do
//     indeed write to lock variables.  The code in draw() should be written
//     as though it were a const method however, at least in terms of what
//     state is changed in the Scene class or subclass.
//
void Scene::draw(DisplayDevice *display) {
  if (!display)
    return;

  if (!display->is_renderer_process()) // master process doesn't draw
    return;

  // check background rendering mode
  if (backgroundmode_changed) {
    display->set_background_mode(backgroundmode);
  }

  // XXX Make a Displayable for lights, as we do for background color, and
  //     have a DispCmd for turning lights on and off.
  if (light_changed) {
    // set up the lights
    for (int i=0; i<DISP_LIGHTS; i++) {
      display->do_define_light(i, lightState[i].color, lightState[i].pos);
      display->do_activate_light(i, lightState[i].on);
    }
  }

#if 0
  // XXX advanced lights are not yet implemented by displays
  // XXX Make a Displayable for lights, as we do for background color, and
  //     have a DispCmd for turning lights on and off.
  if (adv_light_changed) {
    // set up the lights
    for (int i=0; i<DISP_LIGHTS; i++) {
      display->do_define_light(i, lightState[i].color, lightState[i].pos);
      display->do_activate_light(i, lightState[i].on);
    }
  }
#endif

  // update colors
  display->use_colors(colorData);

  // set the background color if any color definitions changed, or if the
  // the color assigned to background changed.  
  // XXX Make a DispCmd for background color so that background color can
  //     be handled automatically by the renderers, rather than special-cased
  //     as we do here.  This should work because the background displayable is
  //     the first item traversed in the scene graph. 
  if (background_color_changed) {
    display->set_background(color_value(background_color_id));
  }

  // set the background gradient colors if any definitions changed, or if
  // the color assigned to the gradient changed.
  if (backgradtop_color_changed || backgradbot_color_changed) {
    display->set_backgradient(color_value(backgradtop->color_id()),
                              color_value(backgradbot->color_id()));
  }

  // clear the display and set the viewport correctly
  display->prepare3D(TRUE);

  // on some machines with broken stereo implementations, we have
  // to draw in stereo even when VMD is set for mono mode
  // XXX this is all very ugly, and while this trick works for drawing
  // in normal mode, it doesn't fix the other single-buffer "stereo"
  // modes, since they don't fall into this test case.  Scene doesn't have
  // intimate knowledge of the details of the display subclass's 
  // stereo mode implementations, so a future rewrite will probably have
  // to move all of this code into the DisplayDevice class, with Scene
  // simply passing in a bunch of parameters so the display subclass
  // can deal with these problems for itself.
  if (display->forced_stereo_draws() && !(display->stereo_mode())) {
    display->left();  // this will be done as normal() since stereo is off
                      // but we'll actually draw to GL_BACK_LEFT
    display->prepareOpaque();
    root.draw(display);
    if (display->prepareTrans()) {
      root.draw(display);
    }

    display->right(); // this will be done as normal() since stereo is off
                      // but we'll actually draw to GL_BACK_RIGHT
    display->prepareOpaque();
    root.draw(display);
    if (display->prepareTrans()) {
      root.draw(display);
    }
  } else {
    // draw left eye first, if stereo
    if (display->stereo_mode())
      display->left();
    else
      display->normal();

    // draw all the objects for normal, or left eye position
    display->prepareOpaque();
    root.draw(display);

    if (display->prepareTrans()) {
      root.draw(display);
    }
   
    // now draw right eye, if in stereo
    if (display->stereo_mode()) {
      display->right();
      display->prepareOpaque();
      root.draw(display);

      if (display->prepareTrans()) {
        root.draw(display);
      }
    } 
  }  

  display->render_done();
  
  // update the display
  display->update(TRUE);
}

// draw the scene to a file in a given format, trying to match the
// view of the given DisplayDevice as closely as possible
// returns TRUE if successful, FALSE if not
int Scene::filedraw(FileRenderer *render, const char *filename, 
                    DisplayDevice *display) {
  int i;

  // Copy all relevant info from the current main DisplayDevice
  (*((DisplayDevice *)render)) = (*display);

  // set up the lights
  for (i=0; i<DISP_LIGHTS; i++) {
    if (lightState[i].on)
      render->do_define_light(i, lightState[i].color, lightState[i].pos);
    render->do_activate_light(i, lightState[i].on);
  }

  // set up the advanced lights
  for (i=0; i<DISP_LIGHTS; i++) {
    if (advLightState[i].on)
      render->do_define_adv_light(i, 
                                  advLightState[i].color,
                                  advLightState[i].pos,
                                  advLightState[i].constfactor,
                                  advLightState[i].linearfactor,
                                  advLightState[i].quadfactor,
                                  advLightState[i].spotdir,
                                  advLightState[i].fallstart,
                                  advLightState[i].fallend,
                                  advLightState[i].spoton);
    render->do_activate_adv_light(i, advLightState[i].on);
  }

  // set up the colors
  render->use_colors(colorData);
  render->set_background(color_value(background->color_id()));
  render->set_backgradient(color_value(backgradtop->color_id()),
                           color_value(backgradbot->color_id()));

  // returns FALSE if not able to open file
  if (!render->open_file(filename)) {
    return FALSE;
  }
  // writes the header
  render->prepare3D(TRUE);

  // draw the scene to the file
  root.draw(render);

  // write the trailer and close
  render->render_done();
  render->update(TRUE);

  // rendering was successful
  return TRUE;
}

// state updates have to be done outside of the draw routine to prevent
// race conditions when multiple slave renderers execute the draw() code.
void Scene::draw_finished() {
  // only master process changes state flags
  backgroundmode_changed=0; // reset background mode change flag
  light_changed = 0;        // reset light changing state cache
  adv_light_changed = 0;    // reset light changing state cache
}

// create the color scale based on current settings
    
////////////////
// definitions for scales with 3 colors
// 1 at val = 0; 0 at mid and greater
static float slope3_lower(float mid, float val) {
  if (val > mid) return 0.0f;
  if (mid == 0) return 1.0f;
  return (1.0f - val / mid);
}
  
// 1 at val = 1; 0 at mid and lesser
static float slope3_upper(float mid, float val) {
  if (val < mid) return 0.0f;
  if (mid == 1) return 1.0f;
  return (1.0f - (1.0f - val) / (1.0f - mid));
} 
    
// 1 at val = mid; 0 at val = 0 an 1
static float slope3_middle(float mid, float val) {
  if (val < mid) {
    return val / mid;
  }
  if (val > mid) { 
    return (1-val) / (1-mid);
  } 
  return 1.0;
}
  
    
// given RGB values at (0,mid,1), return the correct RGB triplet
static void scale_color(const float *lowRGB, const float *midRGB, 
                        const float *highRGB,
                        float mid, float val, float scaleMin, float *rgb) {
  float w1 = slope3_lower(mid, val);
  float w2 = slope3_middle(mid, val);
  float w3 = slope3_upper(mid, val);
  float wsum = (w1 + w2 + w3);

  rgb[0] = (lowRGB[0]*w1 + midRGB[0]*w2 + highRGB[0]*w3) / wsum + scaleMin;
  rgb[1] = (lowRGB[1]*w1 + midRGB[1]*w2 + highRGB[1]*w3) / wsum + scaleMin;
  rgb[2] = (lowRGB[2]*w1 + midRGB[2]*w2 + highRGB[2]*w3) / wsum + scaleMin;

  clamp_color(rgb);
}

void Scene::create_colorscale() {
  const ColorScale &scale = colorScales[scaleActive];

  if (scalePosterize <= 1) 
    scalePosterize = 0;

  int maxidx = MAPCLRS - 1;
  float fmaxidx_inv = 1.0f / float(maxidx);
  float fpost_inv = 1.0f / float(scalePosterize - 1);

  int i;
  switch (scaleTypes[scaleActive].type) {
    case SCALE_TYPE_DIVERGENT:
      for (i=0; i<=maxidx; i++) {
        int idx = (scaleReverse) ? (maxidx-i) : i;
        float *rcol = colorData + 3L*(BEGMAP + idx);
        float relpos = float(i) * fmaxidx_inv;

        relpos = (scalePosterize) ? (int(relpos * scalePosterize) * fpost_inv) : relpos;

        scale_color(scale.min, scale.mid, scale.max, 
                    scaleMid, relpos, scaleMin, rcol);
      }
      break;

    case SCALE_TYPE_TABULATED:
    case SCALE_TYPE_TABULATED_SRGB:
      const float *ctab = scaleTypes[scaleActive].colortable;
      if (ctab == NULL) {
        printf("color tab NULL!\n");
      }

      for (i=0; i<=maxidx; i++) {
        int idx = (scaleReverse) ? (maxidx-i) : i;
        float *rcol = colorData + 3L*(BEGMAP + idx);
        float relpos = float(i) / float(MAPCLRS-1);

        relpos = (scalePosterize) ? (int(relpos * scalePosterize) * fpost_inv) : relpos;

        int tabidx = int(scaleTypes[scaleActive].colortabsz * relpos);
        if (tabidx > scaleTypes[scaleActive].colortabsz)
          tabidx = scaleTypes[scaleActive].colortabsz;

        rcol[0] = ctab[3*tabidx    ];
        rcol[1] = ctab[3*tabidx + 1];
        rcol[2] = ctab[3*tabidx + 2];
      }

      // on-the-fly conversion of tabulated color scales from 
      // sRGB gamma color space to linear gamma RGB color space,
      if (scaleTypes[scaleActive].type == SCALE_TYPE_TABULATED_SRGB) {
        for (i=0; i<=maxidx; i++) {
          float *rcol = colorData + 3L*(BEGMAP + i);
          rcol[0] = powf(rcol[0], 2.2f);
          rcol[1] = powf(rcol[1], 2.2f);
          rcol[2] = powf(rcol[2], 2.2f);
        }
      }
      break;
  }
    
  root.color_scale_changed();
}


int Scene::colorscale_type(int index) {
  if (index < 0 || index >= num_scalemethods)
    return 0;
  return (scaleTypes[scaleActive].type != SCALE_TYPE_TABULATED);
}


int Scene::nearest_index(float r, float g, float b) const {
   const float *rcol = color_value(BEGREGCLRS);  // get the solid colors
   float lsq = r - rcol[0]; lsq *= lsq;
   float tmp = g - rcol[1]; lsq += tmp * tmp;
         tmp = b - rcol[2]; lsq += tmp * tmp;
   float best = lsq;
   int bestidx = BEGREGCLRS;
   for (int n= BEGREGCLRS+1; n < (BEGREGCLRS + REGCLRS); n++) {
      rcol = color_value(n); 
      lsq = r - rcol[0]; lsq *= lsq;
      tmp = g - rcol[1]; lsq += tmp * tmp;
      tmp = b - rcol[2]; lsq += tmp * tmp;
      if (lsq < best) {
       best = lsq;
       bestidx = n;
      }
   }
   return bestidx;
}

int Scene::get_colorscale_colors(int whichScale, 
      float min[3], float mid[3], float max[3]) {
  if (whichScale < 0 || whichScale >= colorScales.num())
    return FALSE;
  const ColorScale &scale = colorScales[whichScale];
  for (int i=0; i<3; i++) {
    min[i] = scale.min[i];
    mid[i] = scale.mid[i];
    max[i] = scale.max[i];
  }
  return TRUE;
}

int Scene::set_colorscale_colors(int whichScale, 
      const float min[3], const float mid[3], const float max[3]) {
  if (whichScale < 0 || whichScale >= colorScales.num())
    return FALSE;
  ColorScale &scale = colorScales[whichScale];
  for (int i=0; i<3; i++) {
    scale.min[i] = min[i];
    scale.mid[i] = mid[i];
    scale.max[i] = max[i];
  }
  create_colorscale();
  return TRUE;
}

