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
 *	$RCSfile: MoleculeGraphics.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.57 $	$Date: 2020/11/06 07:19:59 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * This simply stores and retrieves graphics objects (nominally for use
 *  by the text interface).
 * 
 ***************************************************************************/
#include "DrawMolecule.h"
#include "DispCmds.h"
#include "PickMode.h"

// maximum number of floats required to describe a shape
#define MGMAXDATASIZE 30

class JString;
struct Material;

/// Displayable subclass to generate graphics directly from scripts and plugins
class MoleculeGraphics : public Displayable {
public:
  /// enumeration list of all supported geometric entitities
  enum Shapes {NONE, POINT, PICKPOINT, TRIANGLE, TRINORM, TRICOLOR, 
               LINE, CYLINDER, CONE, SPHERE, TEXT, 
               SPHERETUBE, COLOR, MATERIALS, MATERIAL};

private:
  /// private class used to hold geometric data for a single shape
  struct ShapeClass {
    int id;
    Shapes shape;
    float data[MGMAXDATASIZE];
    int numdata;
    void *extradata; // optional buffer for large primitive data arrays...

    ShapeClass(Shapes s=NONE, int ndata=0, int newid=-1)
              : id(newid), shape(s), numdata(ndata), extradata(NULL) { } 

    void clear() { 
      shape = NONE; 
      if (extradata != NULL)
        free(extradata); 
    }

    int operator==(const ShapeClass&) {return 1;}  // needed for ResizeArray
  };

  /// A ResizeArray containing all shapes in the molecule
  ResizeArray<ShapeClass> shapes;

  /// ShapeClass instances that use text keep a copy of a pointer that's stored
  /// in this array.  MoleculeGraphics takes care of freeing the pointer
  /// when the molecule is deleted.
  ResizeArray<char *> shapetext;

  char graphics_info[250];            ///< text command to regen graphics
  JString *graphics_info_xl;          ///< too-large text gfx regen commands

  int colorID;                        ///< most recently-set draw color index
  int molid;                          ///< id of my parent molecule
  int max_id; 
  int next_id;
  int next_index;                     ///< the end, unless replace was called
  int delete_count;                   ///< total number of deleted elements
  int needRegenerate;                 ///< whether cmdlist needs to be rebuilt
  int added(void);                    ///< actions after adding a new shape
  float cov_scale, cov_pos[3];        ///< center of volume scale factors

  void find_bounds(void);             ///< find bounding box for geometry
  void delete_shapetext();            ///< clear out the shapetext array
  virtual void create_cmdlist(void);  ///< regenerate the draw list

public:
  MoleculeGraphics(DrawMolecule *d) : 
    Displayable(d) {
      molid = d->id();
      max_id = 0;
      next_id = 0;
      next_index = 0;
      needRegenerate = 1;
      delete_count = 0;
      colorID = 0;
      graphics_info_xl = NULL;
  }
  virtual ~MoleculeGraphics(void){
    delete_shapetext();          ///< this is where I delete the data elements
  }
  virtual void prepare() { 
    if (needRegenerate) create_cmdlist(); 
  }

  /// center of volume, and scaling factor
  virtual void cov(float &x, float &y, float &z) {
    find_bounds(); x = cov_pos[0]; y = cov_pos[1]; z = cov_pos[2];
  }
  virtual float scale_factor(void) {
    find_bounds(); return cov_scale;
  }

  // manipulate the data values
  int add_point(const float *x);
  int add_pickpoint(const float *x);
  int add_triangle(const float *x1, const float *x2, const float *x3);
  int add_trinorm(const float *x1, const float *x2, const float *x3,
                  const float *nx1, const float *nx2, const float *nx3);
  int add_tricolor(const float *x1, const float *x2, const float *x3,
                   const float *nx1, const float *nx2, const float *nx3,
      int c1, int c2, int c3);
  int add_line(const float *x, const float *y, int line_style, int width);
  int add_cylinder(const float *x, const float *y, 
                   float radius, int res, int filled);
  int add_cone(const float *x, const float *y, 
               float radius, float radius2, int res);
  int add_sphere(const float *x, float r, int res);
  int add_text(const float *x, const char *text, float size, float thickness);
  int add_spheretube(const int numcoords, const float *xyz3fv, 
                     const int numradii,  const float *radii1fv, 
                     const int numcolors,
                     const float *rgb3fv, const int *colorids,
                     int drawtubes, int res);
  int use_materials(int yes_no);
  int use_color(int index);
  int use_material(const Material *);
 
  void delete_id(int id);   ///< delete an entry
  void delete_all(void);    ///< delete everything
  int replace_id(int id);   ///< same as delete, but next addition goes here
  int index_id(int id);     ///< returns -1 if doesn't exist, else the index
  int num_elements(void){ return int(shapes.num()); }
  int element_id(int index) { return shapes[index].shape != NONE ? shapes[index].id : -1; }

  /// return the command string required to re-make the shape; 
  /// class-internal string storage is set when called and returned
  /// to the caller, and will remain valid until the next such call.
  /// returns NULL if doesn't exist or is NONE
  const char *info_id(int id);

  virtual void pick_start(PickMode *pm, DisplayDevice *d, 
                          int btn, int tag, 
                          const int *cell /* [3] */,
                          int /* dim */, const float *) {

    pm->pick_graphics(molid, tag, btn, d);
  }

};

