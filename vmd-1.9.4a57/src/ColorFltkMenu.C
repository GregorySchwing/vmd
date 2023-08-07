/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2019 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

#include "ColorFltkMenu.h"

#include <stdio.h>
#include <math.h>

#include <FL/fl_draw.H>
#include <FL/forms.H>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Toggle_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Input_Choice.H>
#if defined(VMD_USE_FLTK_COLORCHOOSER)
#include <FL/Fl_Color_Chooser.H>
#endif

#include "Command.h"
#include "Scene.h"
#include "VMDApp.h"
#include "Inform.h"
#include "utilities.h"

/// class to maintain a GUI-usable image of a ColorScale
class ColorscaleImage : public Fl_Box {
  VMDApp *app;
  int grating;
  unsigned char *data;

public:
  ColorscaleImage(int myx, int myy, int myw, int myh, 
                  int grateon, VMDApp *vmdapp) : Fl_Box(myx, myy, myw, myh), app(vmdapp), grating(grateon) {
    data = new unsigned char[3*w()*h()];
  }
  ~ColorscaleImage() { delete [] data; }

protected:
  virtual void draw() {
    if (grating) {
      // To validate color map usability we draw draw a 
      // Kovesi style 2-D color scale test w/ sinusoidal "grating" effect, 
      // by adding high-frequency 10% sinusoid to the value
      // See:  https://arxiv.org/abs/1509.03700
      const float *rgb = app->scene->color_value(REGCLRS);
      int xs = w(), ys = h();
      for (int ix=0; ix<xs; ix++) {
        for (int iy=0; iy<ys; iy++) {
          float vertscale = (0.05f * (ys-iy-1) / float(ys));
          float singrate =  sinf((xs / 5.0f) * VMD_TWOPIF * ix / float(xs));
          float csidx = (vertscale * singrate) + (ix/float(xs));

          // map to integer color scale index range   
          int idx = int(csidx * MAPCLRS);
          idx = (idx < 0)        ?           0 : idx;
          idx = (idx >= MAPCLRS) ? (MAPCLRS-1) : idx;

          int idx3 = idx*3;
          int img3 = (xs*iy + ix) * 3;

          unsigned char r = (unsigned char)(255*rgb[idx3 + 0]);
          unsigned char g = (unsigned char)(255*rgb[idx3 + 1]);
          unsigned char b = (unsigned char)(255*rgb[idx3 + 2]);

          data[img3    ] = r;
          data[img3 + 1] = g;
          data[img3 + 2] = b;
        }
      }
    } else {
      // just draw the plain colorscale data itself with no perturbation
      for (int i=0; i<w(); i++) {
        const float *rgb = app->scene->color_value(i*MAPCLRS/w()+REGCLRS);
        unsigned char r = (unsigned char)(255*rgb[0]);
        unsigned char g = (unsigned char)(255*rgb[1]);
        unsigned char b = (unsigned char)(255*rgb[2]);
        for (int j=0; j<h(); j++) {
          data[3*(w()*j + i)    ] = r;
          data[3*(w()*j + i) + 1] = g;
          data[3*(w()*j + i) + 2] = b;
        }
      }

    }

    fl_draw_image(data, x(), y(), w(), h());
  }

};


//#define VMD_FLCHART_WORKAROUND 1
#if defined(VMD_FLCHART_WORKAROUND)

// Use Fltk's definition of rint, since Windows sees fit to omit rint
// from its math library.  This definition isn't correct for negative
// arguments, but it's what's used in the original Fl_Chart class,
// from which this code was derived, so I'll keep doing it there way.
static double fltk_rint(double v) {return floor(v+.5);}

/// A replacement for Fl_Chart, which has had a bug that causes it 
/// to sometimes draw outside its bounding box.
/// Furthermore, by writing our own, we can customize its behavior 
/// significantly.
class myFl_Chart : public Fl_Widget {
private:
  int num;
  float *values;
  float min, max;
  int imin, imax;

public:
  myFl_Chart(int x, int y, int w, int h, const char *l=0) : Fl_Widget(x,y,w,h,l) {
    box(FL_BORDER_BOX);
    align(FL_ALIGN_BOTTOM);
    num = 0;
    values = NULL;
    min = max = 0;
    imin = imax = 0;
  }

  ~myFl_Chart() {
    delete [] values;
  }

  void clear() {
    delete [] values;
    values = NULL;
    num = 0;
    redraw();
  }

  void set_data(const float *data, int n) {
    if (n < 1) {
      clear();
      return;
    }
    delete [] values;
    values = new float[n];
    memcpy(values, data, n*sizeof(float));
    num = n;
    min = max = data[0];
    imin = imax = 0;
    for (int i=1; i<n; i++) {
      if (min > data[i]) { min = data[i]; imin = i; }
      if (max < data[i]) { max = data[i]; imax = i; }
    }
    redraw();
  }

  void bounds(float newmin, float newmax) {
    min = newmin;
    max = newmax;
  }

protected:
  void draw() {
    int xx, yy, ww, hh;
    if (!num) return;
    xx = x()+9;
    yy = y()+9;
    ww = w()-2*9;
    hh = h()-2*9;

    draw_box();

    double lh = fl_height(); // compensate for text height?
    double incr;
    int zeroh;
    if (min > 0) {
      incr = (hh-2*lh)/max;
      zeroh = yy+hh-9;
    } else if (max < 0) {
      incr = (hh-2*lh)/min;
      zeroh = yy-9;
    } else {
      incr = (hh-2*lh)/(max-min);
      zeroh = yy+hh+(int)fltk_rint(min*incr) - 9;
    }
    double bwidth = ww/double(num);

    for (int i=1; i<num; i++) {
      int x0 = xx    + (int)fltk_rint((i-.5)*bwidth);
      int x1 = xx    + (int)fltk_rint((i+.5)*bwidth);
      int y0 = zeroh - (int)fltk_rint(values[i-1]*incr);
      int y1 = zeroh - (int)fltk_rint(values[i]*incr);
      int color = FL_GREEN;

      if (i == imin) color = FL_RED;
      else if (i == imax) color = FL_BLUE;
      fl_color(color);
      if ((values[i-1]>0.0)!=(values[i]>0.0)) {
        double ttt = values[i-1]/(values[i-1]-values[i]);
        int xt = xx + (int)fltk_rint((i-.5+ttt)*bwidth);
        fl_polygon(x0,zeroh, x0,y0, xt,zeroh);
        fl_polygon(xt,zeroh, x1,y1, x1,zeroh);
      } else {
        fl_polygon(x0, zeroh, x0, y0, x1, y1, x1, zeroh);
      }
      fl_color(FL_BLACK);
      fl_line(x0,y0,x1,y1);
    }
    fl_line(xx, zeroh, xx+ww, zeroh);

#if 0
    char buf[30] = { 0 };
    sprintf(buf, "%d: %3.2f", imin, min);
    fl_draw(buf, xx+(int)fltk_rint((imin+.5)*bwidth), zeroh-(int)fltk_rint(min*incr),0,0,
         min >= 0 ? FL_ALIGN_BOTTOM : FL_ALIGN_TOP);
    sprintf(buf, "%d: %3.2f", imax, max);
    fl_draw(buf, xx+(int)fltk_rint((imax+.5)*bwidth), zeroh-(int)fltk_rint(max*incr),0,0,
         max >= 0 ? FL_ALIGN_BOTTOM : FL_ALIGN_TOP);
#endif
  }

};

#endif // VMD_FLCHART_WORKAROUND


#if defined(VMD_FLCHART_WORKAROUND)
#define CHARTCLASS myFl_Chart 
#else
#define CHARTCLASS Fl_Chart
#endif

/// class to maintain a GUI-usable image of a ColorScale
class ColorscaleLumaChart : public CHARTCLASS {
  VMDApp *app;
  ResizeArray<float> values;

public:
  ColorscaleLumaChart(int myx, int myy, int myw, int myh, VMDApp *vmdapp) : CHARTCLASS(myx, myy, myw, myh), app(vmdapp) {
  }
  ~ColorscaleLumaChart() { }

protected:
  virtual void draw() {
    values.clear();
    // just draw the plain colorscale data itself with no perturbation
    int i;
    for (i=0; i<w(); i++) {
      const float *rgb = app->scene->color_value(i*MAPCLRS/w()+REGCLRS);
      float r = rgb[0];
      float g = rgb[1];
      float b = rgb[2];

      // convert linear gamma RGB to CIEXYZ
//      float Xd65 = r*0.41239f + g*0.35758f + b*0.18048f;
      float Yd65 = r*0.21264f + g*0.71517f + b*0.07219f;
//      float Zd65 = r*0.01933f + g*0.11919f + b*0.95053f;

      // convert CIEXYZ to CIELAB
//      float fXd65 = (Xd65 > 0.008856) ? powf(Xd65, 1.0f/3.0f) : Xd65 * 7.787f + 16.0f/116.0f;
      float fYd65 = (Yd65 > 0.008856) ? powf(Yd65, 1.0f/3.0f) : Yd65 * 7.787f + 16.0f/116.0f;
//      float fZd65 = (Zd65 > 0.008856) ? powf(Zd65, 1.0f/3.0f) : Zd65 * 7.787f + 16.0f/116.0f;
      float labL = (116.0f * fYd65) - 16.0f;
//      float laba = 500.0f * (fXd65 - fYd65);
//      float labb = 500.0f * (fYd65 - fZd65);
      values.append(labL);
    }

#if defined(VMD_FLCHART_WORKAROUND)
    // XXX workaround Fl_Chart bug that could cause out-of-bounds drawing
    myFl_Chart::set_data(&(values[0]), values.num());
    myFl_Chart::draw();
#else 
    clear();
    type(FL_LINE_CHART);
    bounds(0.0, 100.0);
    for (i=0; i<values.num(); i++) {
      add(values[i]);
    }
    Fl_Chart::draw();
#endif // VMD_FLCHART_WORKAROUND
  }

};



void ColorFltkMenu::make_window() {
  size(475, 495);
  { 
    { Fl_Tabs* o = new Fl_Tabs(0, 0, 475, 495);
#if defined(VMDMENU_WINDOW)
      o->color(VMDMENU_WINDOW, FL_GRAY);
      o->selection_color(VMDMENU_WINDOW);
#endif

    //
    // Color editing tabs
    //

#if defined(VMD_USE_FLTK_COLORCHOOSER)
      // use FLTK-provided color chooser
      { Fl_Group* o = new Fl_Group(0, 20, 475, 425, "Color Definitions");
#if defined(VMDMENU_WINDOW)
        o->color(VMDMENU_WINDOW, FL_GRAY);
        o->selection_color(VMDMENU_WINDOW);
#endif
        //o->hide();
        { Fl_Hold_Browser* o = colordefbrowser = new Fl_Hold_Browser(10, 195, 110, 205);
          o->labeltype(FL_NO_LABEL);
          o->color(VMDMENU_BROWSER_BG, VMDMENU_BROWSER_SEL);
          o->callback(colordef_cb, this);
          VMDFLTKTOOLTIP(o, "Select color name to adjust RGB color definition")
        }
        { Fl_Color_Chooser *o = colchooser = new Fl_Color_Chooser(125, 195, 280, 205);
          colchooser->rgb(1.0, 1.0, 1.0);
          o->callback(colchooser_cb, this);
        }
        { Fl_Box *o = colchooserbox = new Fl_Box(FL_DOWN_BOX, 410, 195, 55, 205, "");
          o->color(fl_rgb_color(255, 255, 255));
        }
        colchoosedefaultbutton = new Fl_Button(230, 405, 85, 25, "Default");
#if defined(VMDMENU_WINDOW)
        colchoosedefaultbutton->color(VMDMENU_WINDOW, FL_GRAY);
#endif
        colchoosedefaultbutton->callback(default_cb, this);
        VMDFLTKTOOLTIP(colchoosedefaultbutton, "Reset to original RGB color")
#endif


#if defined(VMD_USE_VMD_COLORCHOOSER)
      { Fl_Group* o = new Fl_Group(0, 20, 475, 425, "Color Definitions");
#if defined(VMDMENU_WINDOW)
        o->color(VMDMENU_WINDOW, FL_GRAY);
        o->selection_color(VMDMENU_WINDOW);
#endif
        { Fl_Hold_Browser* o = colordefbrowser = new Fl_Hold_Browser(15, 195, 135, 120);
          o->labeltype(FL_NO_LABEL);
          o->color(VMDMENU_BROWSER_BG, VMDMENU_BROWSER_SEL);
          o->callback(colordef_cb, this);
          VMDFLTKTOOLTIP(o, "Select color name to adjust RGB color definition")
        }
        { Fl_Value_Slider* o = redscale = new Fl_Value_Slider(160, 195, 300, 20);
          o->type(FL_HORIZONTAL);
          o->color(VMDMENU_COLOR_RSLIDER);
          o->callback(rgb_cb, this);
          VMDFLTKTOOLTIP(o, "Adjust slider to change RGB color definition")
        }
        { Fl_Value_Slider* o = greenscale = new Fl_Value_Slider(160, 215, 300, 20);
          o->type(FL_HORIZONTAL);
          o->color(VMDMENU_COLOR_GSLIDER);
          o->callback(rgb_cb, this);
          VMDFLTKTOOLTIP(o, "Adjust slider to change RGB color definition")
        }
        { Fl_Value_Slider* o = bluescale = new Fl_Value_Slider(160, 235, 300, 20);
          o->type(FL_HORIZONTAL);
          o->color(VMDMENU_COLOR_BSLIDER);
          o->callback(rgb_cb, this);
          VMDFLTKTOOLTIP(o, "Adjust slider to change RGB color definition")
        }
        { Fl_Button* o = grayscalebutton = new Fl_Button(165, 265, 85, 25, "Grayscale");
          o->type(FL_TOGGLE_BUTTON);
#if defined(VMDMENU_WINDOW)
          o->color(VMDMENU_WINDOW, FL_GRAY);
#endif
          VMDFLTKTOOLTIP(o, "Lock sliders for grayscale color")
        }
        defaultbutton = new Fl_Button(290, 265, 85, 25, "Default");
#if defined(VMDMENU_WINDOW)
        defaultbutton->color(VMDMENU_WINDOW, FL_GRAY);
#endif
        defaultbutton->callback(default_cb, this);
        VMDFLTKTOOLTIP(defaultbutton, "Reset to original RGB color")
#endif

        // color category/name/color browsers common to both code versions
        { Fl_Hold_Browser* o = categorybrowser = new Fl_Hold_Browser(10, 75, 125, 100, "Categories");
          o->align(FL_ALIGN_TOP);
          o->color(VMDMENU_BROWSER_BG, VMDMENU_BROWSER_SEL);
          o->callback(category_cb, this);
          VMDFLTKTOOLTIP(o, "Select color category then name to set active color")
        }
        { Fl_Hold_Browser* o = itembrowser = new Fl_Hold_Browser(140, 75, 210, 100, "Names");
          o->align(FL_ALIGN_TOP);
          o->color(VMDMENU_BROWSER_BG, VMDMENU_BROWSER_SEL);
          o->callback(item_cb, this);
          VMDFLTKTOOLTIP(o, "Select color category then name to set active color")
        }
        { Fl_Hold_Browser* o = colorbrowser = new Fl_Hold_Browser(355, 75, 110, 100, "Colors");
          o->align(FL_ALIGN_TOP);
          o->color(VMDMENU_BROWSER_BG, VMDMENU_BROWSER_SEL);
          o->callback(color_cb, this);
          VMDFLTKTOOLTIP(o, "Select color category then name to set active color")
        }
        new Fl_Box(10, 30, 190, 25, "Assign colors to categories:");

        o->end();
      }

      { Fl_Group* o = new Fl_Group(0, 20, 475, 475, "Color Scale");
#if defined(VMDMENU_WINDOW)
        o->color(VMDMENU_WINDOW, FL_GRAY);
        o->selection_color(VMDMENU_WINDOW);
#endif
        o->hide();
        { Fl_Choice* o = scalemethod = new Fl_Choice(15, 40, 110, 25, "Method");
          o->color(VMDMENU_CHOOSER_BG, VMDMENU_CHOOSER_SEL);
          o->down_box(FL_BORDER_BOX);
          o->align(FL_ALIGN_TOP);
          o->callback(scalemethod_cb, this);
        }
        offsetvalue = new Fl_Value_Slider(200, 25, 265, 20, "Offset");
        offsetvalue->type(FL_HORIZONTAL);
        offsetvalue->color(VMDMENU_SLIDER_BG, VMDMENU_SLIDER_FG);
        offsetvalue->align(FL_ALIGN_LEFT);
        offsetvalue->range(-1.0, 1.0);
        offsetvalue->callback(scalesettings_cb, this);
        { Fl_Value_Slider* o = midpointvalue = new Fl_Value_Slider(200, 55, 265, 20, "Midpoint");
          o->type(FL_HORIZONTAL);
          midpointvalue->align(FL_ALIGN_LEFT);
          midpointvalue->color(VMDMENU_SLIDER_BG, VMDMENU_SLIDER_FG);
          o->range(0.0, 1.0);
          o->callback(scalesettings_cb, this);
        }
        { Fl_Toggle_Button *o = scalereversebutton = new Fl_Toggle_Button(15, 85, 110, 25, "Reverse");
          o->callback(scalesettings_cb, this);
        }
        { Fl_Input_Choice *o = scaleposterize = new Fl_Input_Choice(200, 85, 100, 25, "Posterize");
          o->add("Off");
          o->add("2");
          o->add("3");
          o->add("4");
          o->add("5");
          o->add("6");
          o->add("7");
          o->add("8");
          o->add("9");
          o->add("10");
          o->add("11");
          o->add("12");
          o->add("13");
          o->add("14");
          o->add("15");
          o->add("16");
          o->value("Off");
          o->callback(scalesettings_cb, this);
        }


        new Fl_Box(10, 135, 425, 20, "High spatial frequency color scale test grating");
        new Fl_Box(435, 155, 30, 20, "10%");
        new Fl_Box(435, 220, 30, 20, "5%");
        new Fl_Box(435, 285, 30, 20, "0%");
        image = new ColorscaleImage(10, 155, 425, 150, 1, app);
        new Fl_Box(10, 315, 425, 20, "CIELAB L* perceptual lightness");
        scalelumachart = new ColorscaleLumaChart(10, 335, 425, 150, app);
        new Fl_Box(435, 345, 30, 20, "100"); // down 10px to 100% line
        new Fl_Box(435, 455, 30, 20, "0"); // up 10px to Fl_Chart zero line

        o->end();
      }
      o->end();
    }
    end();
  }
}


ColorFltkMenu::ColorFltkMenu(VMDApp *vmdapp) 
: VMDFltkMenu("color", "Color Controls", vmdapp) {
  make_window();

  command_wanted(Command::COLOR_SCALE_METHOD);
  command_wanted(Command::COLOR_SCALE_SETTINGS);
  command_wanted(Command::COLOR_SCALE_COLORS);
  command_wanted(Command::COLOR_CHANGE);
  command_wanted(Command::COLOR_NAME);
  command_wanted(Command::MOL_RENAME);
  command_wanted(Command::COLOR_ADD_ITEM);

  // set up color scale values
  for (int j=0; j<vmdapp->num_colorscale_methods(); j++)
    scalemethod->add(app->colorscale_method_menuname(j));

  reset_color_categories();
  reset_color_names();
  reset_color_scale();
}

void ColorFltkMenu::reset_color_categories() {
  categorybrowser->clear();
  int n = app->num_color_categories();
  for (int j=0; j<n; j++)
    categorybrowser->add(app->color_category(j));
  categorybrowser->value(0);  //nothing selected
}

void ColorFltkMenu::reset_color_names() {
  colorbrowser->clear();
  colordefbrowser->clear();
  int n = app->num_regular_colors();
  for (int j=0; j<n; j++) {
    char buf[128] = { 0 };
    sprintf(buf, "%d %s", j, app->color_name(j));
    colorbrowser->add(buf);
    colordefbrowser->add(buf);
  }
  colorbrowser->value(0);    // nothing selected
  colordefbrowser->value(0); // nothing selected
}

void ColorFltkMenu::reset_color_scale() {
  int idx = app->colorscale_method_current();
  const char *menuname = app->colorscale_method_menuname(idx);
  set_chooser_from_string(menuname, scalemethod);

  float mid, min, max;
  int rev, posterize;
  app->colorscale_params(&mid, &min, &max, &rev, &posterize);
  offsetvalue->value(min);
  midpointvalue->value(mid);
  scalereversebutton->value(rev);

  // make controls active/inactive according to
  // the type of colorscale that is active
  if (app->scene->colorscale_type(idx)) {
    offsetvalue->activate();
    midpointvalue->activate(); 
  } else {
    offsetvalue->deactivate();
    midpointvalue->deactivate(); 
  }

  image->redraw();
  scalelumachart->redraw();
}

void ColorFltkMenu::update_chosen_color() {
  int catval = categorybrowser->value();
  int itemval = itembrowser->value();
  if (!catval || !itemval) {
    colordefbrowser->value(0);
    return;
  }
  const char *category = categorybrowser->text(catval);
  const char *item = itembrowser->text(itemval);
  const char *color = app->color_mapping(category, item);
  int index = app->color_index(color);
  if (index < 0) {
    msgErr << "ColorFltkMenu::update_chosen_color: invalid color" << sendmsg;
    return;
  }
  colorbrowser->value(index+1);
  colorbrowser->visible(index+1);
  colordefbrowser->value(index+1);
  colordefbrowser->visible(index+1);
  update_color_definition();
}

void ColorFltkMenu::update_color_definition() {
  float r, g, b;
  const char *colorname = app->color_name(colordefbrowser->value()-1);
  if (!app->color_value(colorname, &r, &g, &b)) return;
#if defined(VMD_USE_FLTK_COLORCHOOSER)
  colchooser->rgb(r, g, b);
#endif
#if defined(VMD_USE_VMD_COLORCHOOSER)
  redscale->value(r);
  greenscale->value(g);
  bluescale->value(b);
#endif
}

int ColorFltkMenu::act_on_command(int type, Command *) {
  switch (type) {
    case Command::COLOR_SCALE_METHOD:
    case Command::COLOR_SCALE_SETTINGS:
    case Command::COLOR_SCALE_COLORS:
      reset_color_scale();
      break;
    case Command::COLOR_CHANGE:
      update_color_definition();
      break;
    case Command::COLOR_NAME:
    case Command::MOL_RENAME:
      update_chosen_color();
      break;
    case Command::COLOR_ADD_ITEM:
      category_cb(NULL, this);
      break;
    default:
      ;
  }
  return FALSE;
}

void ColorFltkMenu::category_cb(Fl_Widget *, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  int val = self->categorybrowser->value();
  if (!val) return;
  const char *category = self->categorybrowser->text(val);
  int n = self->app->num_color_category_items(category);
  self->itembrowser->clear();
  for (int i=0; i<n; i++)
    self->itembrowser->add(self->app->color_category_item(category, i));
  self->itembrowser->value(0);
  self->colorbrowser->value(0);
}

void ColorFltkMenu::item_cb(Fl_Widget *, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  self->update_chosen_color();
}

void ColorFltkMenu::color_cb(Fl_Widget *, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  int catval = self->categorybrowser->value();
  int itemval = self->itembrowser->value();
  int colorval = self->colorbrowser->value();
  if (!catval || !itemval || !colorval) return;
  const char *category = self->categorybrowser->text(catval);
  const char *item = self->itembrowser->text(itemval);
  const char *color = self->app->color_name(colorval-1);
  self->app->color_change_name(category, item, color);
  self->update_chosen_color();
}

void ColorFltkMenu::colordef_cb(Fl_Widget *, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  self->update_color_definition();
}


#if defined(VMD_USE_FLTK_COLORCHOOSER)
// callback for FLTK-provided color chooser
void ColorFltkMenu::colchooser_cb(Fl_Widget *w, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;

  float r, g, b;
  r = (float)self->colchooser->r();
  g = (float)self->colchooser->g();
  b = (float)self->colchooser->b();
  self->colchooserbox->color(fl_rgb_color(uchar(255*r), 
                                          uchar(255*g), 
                                          uchar(255*b)));
  self->colchooserbox->redraw();

  int val = self->colordefbrowser->value();
  if (!val) return; // don't proceed further unless a color item is selected

  const char *color = self->colordefbrowser->text(val);
  self->app->color_change_rgb(color, r, g, b);
}
#endif


#if defined(VMD_USE_VMD_COLORCHOOSER)
// callback for classic VMD color sliders
void ColorFltkMenu::rgb_cb(Fl_Widget *w, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  int val = self->colordefbrowser->value();
  if (!val) return;
  const char *color = self->colordefbrowser->text(val);
  float r, g, b;
  if (self->grayscalebutton->value()) {
    r = g = b = (float)((Fl_Value_Slider *)w)->value();
  } else {
    r = (float)self->redscale->value();
    g = (float)self->greenscale->value();
    b = (float)self->bluescale->value();
  }
  self->app->color_change_rgb(color, r, g, b);
}
#endif


void ColorFltkMenu::default_cb(Fl_Widget *, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  int val = self->colordefbrowser->value();
  if (!val) return;
  const char *color = self->colordefbrowser->text(val);
  float r, g, b;
  if (!self->app->color_default_value(color, &r, &g, &b)) {
    msgErr << "ColorFltkMenu::default_cb(): invalid color" << sendmsg;
    return;
  }
  self->app->color_change_rgb(color, r, g, b);
}

void ColorFltkMenu::scalemethod_cb(Fl_Widget *w, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  Fl_Choice *choice = (Fl_Choice *)w;

  int idx = choice->value();
  const char *mn = choice->text();

  //
  // Find the index of a string that matches the name of the given menu
  // item.  XXX Only checks leaf node menu names, not full pathnames currently
  //
  // FLTK 1.1.4
  int i, strindex;
  int num = self->app->num_colorscale_methods();
  for (strindex=-1, i=0; i < num; i++) {
    // find leaf menu name from full menu path
    const char *colstr;
    if ((colstr = strrchr(self->app->colorscale_method_menuname(i), '/')) == NULL)
      colstr = self->app->colorscale_method_menuname(i);
    else
      colstr++;

    // compare leaf submenu item name against left color mode name
    if (!strcmp(colstr, mn)) {
      strindex=i;
      break;
    }
  }
  idx = strindex;

  self->app->colorscale_setmethod(idx);
}

void ColorFltkMenu::scalesettings_cb(Fl_Widget *, void *v) {
  ColorFltkMenu *self = (ColorFltkMenu *)v;
  float mid, min, max;
  int rev, posterize;
  self->app->colorscale_params(&mid, &min, &max, &rev, &posterize);
  mid = (float)self->midpointvalue->value();
  min = (float)self->offsetvalue->value();
  rev = self->scalereversebutton->value();
  const char *posterizestr = self->scaleposterize->value();
  posterize = 0;
  if (!strupcmp(posterizestr, "off")) {
    posterize = 0;
  } else {
    if (sscanf(posterizestr, "%d", &posterize) != 1) {
      posterize = 0;
    }
  }

  self->app->colorscale_setparams(mid, min, max, rev, posterize);
}



