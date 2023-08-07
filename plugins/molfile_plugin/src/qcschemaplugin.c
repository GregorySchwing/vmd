/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/* This is a plugin that will read input from a JSON file
** generated from QCSchema MolSSI efforts
** some more details will go here soon 
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include <sys/stat.h>
#include "qcschema_json.c"

#if defined(_AIX)
#include <strings.h>
#endif

#if defined(WIN32) || defined(WIN64)
#define strcasecmp stricmp
#endif

#include "molfile_plugin.h"
#include "unit_conversion.h"
#include "periodic_table.h"
#include "qmplugin.h"


#define ALLOCATE(array, type, size) \
  array = (type *)calloc(size, sizeof(type)); \
  if (array == NULL) { \
    fprintf(stderr, "qcschemaplugin) Memory allocation for %s failed!\n", #array); \
    return FALSE; \
  }

#define GET_LINE(x,y) if (!fgets(x, sizeof(x), y)) return FALSE

/* Read the basis set data */
static int get_basis (qmdata_t *);
static int shelltype_int(char *type);
static int fill_basis_arrays(qmdata_t *data);

/*static int read_geom_block(qmdata_t *data);*/
static int read_molecular_orbitals(qmdata_t *data);
static int read_wave_coeffs(FILE *file, qm_wavefunction_t *wave);
static int count_orbitals(qmdata_t *data);

typedef struct {
  json_value* geom_array;   /* store geometry array */
  json_value* symbol_array; /* store symbol array */
  json_value* model;        /* store method, basis, etc */
  json_value* keywords;     /* store keywords */
  json_value* extras;       /* store extras */
  json_value* provenance;   /* store creator program, version, routine, etc */
  json_value* properties;   /* store properties */
  json_value* return_result;   /* store result */
  int coordsonly;
} jsondata_t;

/*********************************************************
 *
 * Open file and fill corresponding json_values.
 * After, file is not used. Returns natoms.
 *
 *********************************************************/
static void *open_qcschema_read(const char *filename,
                              const char *filetype,
                              int *natoms) {
  FILE *fd;
  qmdata_t *data = NULL;
  char buffer[1024];
  char keystring[20];
  json_char* json;
  json_value* value;
  json_value* aux_value;
  char* file_contents;
  struct stat filestatus;
  int file_size;

  if ( stat(filename, &filestatus) != 0) {
    fprintf(stderr, "File %s not found\n", filename);
    return NULL;
  }
  file_size = filestatus.st_size;
  file_contents = (char*)malloc(filestatus.st_size);
  if ( file_contents == NULL) {
    fprintf(stderr, "Memory error: unable to allocate %d bytes\n", file_size);
    return NULL;
  }

  fd = fopen(filename, "rb");
  if (!fd) return NULL;
  
  /* allocate memory for main QM data structure */
  data = init_qmdata();
  if (!data) return NULL;

  data->file = fd;

  /* allocate JSON specific data */
  jsondata_t * jsondata = (jsondata_t *)calloc(1, sizeof(jsondata_t));
  if (!jsondata) return NULL;
  
  if ( fread(file_contents, file_size, 1, fd) != 1 ) {
    fprintf(stderr, "Unable to read content of %s\n", filename);
    fclose(fd);
    free(file_contents);
    return NULL;
  }

  fclose(fd);
  json = (json_char*)file_contents;
  
  value = json_parse(json,file_size);

  int i, j;
  if (value == NULL) {
    return NULL;
  }
  
  _Bool check_schema_done = FALSE;
  _Bool molecule_done = FALSE;
  _Bool driver_done = FALSE;
  _Bool model_done = FALSE;
  _Bool extras_done = FALSE;
  _Bool properties_done = FALSE;
  _Bool return_result_done = FALSE;

  /* loop over main objects */
  for (i = 0; i < value->u.object.length; i++) {
    /*printf("object[%d].name = %s\n", i, value->u.object.values[i].name);*/
    
    /* Check if the file is QCSCHEMA format */
    if (!check_schema_done && !strcmp(value->u.object.values[i].name, "schema_name")) {
      /* we have a QCSchema file */
      check_schema_done = TRUE;
    }

    else if (!molecule_done && !strcmp(value->u.object.values[i].name, "molecule")) {
      /* read molecule part */
      aux_value = value->u.object.values[i].value;
      /* in molecule, the order of 'geometry' and 'symbols' is not fixed */
      for (j = 0; j < aux_value->u.object.length; j++) {
        if (!strcmp(aux_value->u.object.values[j].name, "geometry")) {
          /* number of atoms must be array size / 3 */
          data->numatoms = aux_value->u.object.values[j].value->u.array.length / 3;
          jsondata->geom_array = aux_value->u.object.values[j].value;
        } 
        else if (!strcmp(aux_value->u.object.values[j].name, "symbols")) {  
          /* get numatoms from size */
          data->numatoms = aux_value->u.object.values[j].value->u.array.length;
          jsondata->symbol_array = aux_value->u.object.values[j].value;
        }
        else if (!strcmp(aux_value->u.object.values[j].name, "molecular_charge")) {  
          data->totalcharge = aux_value->u.object.values[j].value;
        }
        else if (!strcmp(aux_value->u.object.values[j].name, "molecular_multiplicity")) {  
          data->multiplicity = aux_value->u.object.values[j].value;
        }
      }
      molecule_done = TRUE;
    }

    else if (!driver_done && !strcmp(value->u.object.values[i].name, "driver")) {
      aux_value = value->u.object.values[i].value;
      /* there are only 4 types supported: energy, gradient, hessian, properties */
      if (!strcmp(aux_value->u.string.ptr, "energy")) {
        data->runtype = MOLFILE_RUNTYPE_ENERGY;
      }
      else if (!strcmp(aux_value->u.string.ptr, "gradient")) {
        data->runtype = MOLFILE_RUNTYPE_GRADIENT;
      }
      else if (!strcmp(aux_value->u.string.ptr, "hessian")) {
        data->runtype = MOLFILE_RUNTYPE_HESSIAN;
      }
      driver_done = TRUE;
    }

    else if (!model_done && !strcmp(value->u.object.values[i].name, "model")) {
      jsondata->model = value->u.object.values[i].value;
      model_done = TRUE;
    }

    else if (!extras_done && !strcmp(value->u.object.values[i].name, "extras")) {
      jsondata->extras = value->u.object.values[i].value;
      extras_done = TRUE;
    }

    else if (!properties_done && !strcmp(value->u.object.values[i].name, "properties")) {
      jsondata->properties = value->u.object.values[i].value;
      properties_done = TRUE;
    }

    else if (!return_result_done && !strcmp(value->u.object.values[i].name, "return_result")) {
      jsondata->return_result = value->u.object.values[i].value;
      return_result_done = TRUE;
    }
  }

  if (!check_schema_done) {
    printf("qcschemaplugin) The file is not in JSON/QCSCHEMA format!\n");
    return NULL;
  }
  

  /* allocate JSON specific data */
  data->format_specific_data = jsondata;
  *natoms = data->numatoms;
  data->num_frames = 1;
  return data;
}


/**********************************************************
 *
 * Read geometry from json_value object
 *
 *********************************************************/
static int read_qcschema_structure(void *mydata, int *optflags, 
                                 molfile_atom_t *atoms) 
{
  int i;
  char buffer[1024];
  char atname[1024];
  int num, atomicnum;
  molfile_atom_t *atom;
  qmdata_t *data = (qmdata_t *)mydata;
  jsondata_t* jsondata = (jsondata_t*)data->format_specific_data;

  ALLOCATE(data->atoms, qm_atom_t, data->numatoms);

  /* atomic number is provided by plugin.
   * (This is required for QM plugins!) */
  *optflags = MOLFILE_ATOMICNUMBER;

  float unitfac = 1.f; /* to use in case of bohr units */

  for (i=0; i<data->numatoms; i++) {
    atom = atoms+i;

    strncpy(atname,jsondata->symbol_array->u.array.values[i]->u.string.ptr,sizeof(atname));

    strncpy(atom->name, atname, sizeof(atom->name)); 
    strncpy(atom->type, atom->name, sizeof(atom->type));
    atom->atomicnumber = get_pte_idx_from_string(atname);
    strncpy(atom->resname,"MOL", sizeof(atom->resname));
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
    data->atoms[i].atomicnum = atom->atomicnumber;
    data->num_frames_read = 0;

    /* keep local copy */
    strncpy(data->atoms[i].type, atname, sizeof(data->atoms[i].type));
    data->atoms[i].atomicnum = atomicnum;
    data->atoms[i].x = jsondata->geom_array->u.array.values[i*3]->u.dbl*unitfac;
    data->atoms[i].y = jsondata->geom_array->u.array.values[i*3+1]->u.dbl*unitfac;
    data->atoms[i].z = jsondata->geom_array->u.array.values[i*3+2]->u.dbl*unitfac;
    
  }

  return MOLFILE_SUCCESS;

}


/***********************************************************
 *
 * Provide non-QM metadata for next timestep. 
 * Required by the plugin interface.
 *
 ***********************************************************/
static int read_timestep_metadata(void *mydata,
                                  molfile_timestep_metadata_t *meta) {
  
  meta->count = -1;
  meta->has_velocities = 0;

  return MOLFILE_SUCCESS;
}


/***********************************************************
 *
 * We are not reading the coefficients themselves,
 * because that could require a large amount of memory.
 *
 ***********************************************************/
static int read_qm_timestep_metadata(void *mydata,
                                    molfile_qm_timestep_metadata_t *meta) {
  qmdata_t *data = (qmdata_t *)mydata;
  jsondata_t* jsondata = (jsondata_t *)data->format_specific_data;
 
  if (data->num_frames_sent >= data->num_frames) {
    /* All frames were sent. */
    return MOLFILE_ERROR;
  }

  /* Count the number of cartesian basis functions in 
     the basis set */
  if (data->num_frames_sent == data->num_frames-1) {
    int i;
    qm_timestep_t *cur_ts;

    if (!count_orbitals(data)) return MOLFILE_ERROR;

    /* get a pointer to the current qm timestep */
    cur_ts = data->qm_timestep;
    
    for (i=0; (i<MOLFILE_MAXWAVEPERTS && i<cur_ts->numwave); i++) {
      meta->num_orbitals_per_wavef[i] = cur_ts->wave[i].num_orbitals;
      meta->has_occup_per_wavef[i]    = cur_ts->wave[i].has_occup;
      meta->has_orben_per_wavef[i]    = cur_ts->wave[i].has_orben;
    }
    meta->wavef_size   = data->wavef_size;
    meta->num_wavef    = cur_ts->numwave;
    meta->num_scfiter  = cur_ts->num_scfiter;
    meta->has_gradient = FALSE;
    meta->num_charge_sets = 0;
  }
  return MOLFILE_SUCCESS;
}



/***********************************************************
 *
 * Provides VMD with the data of the next timestep.
 *
 ***********************************************************/
static int read_timestep(void *mydata, int natoms, 
       molfile_timestep_t *ts, molfile_qm_metadata_t *qm_metadata,
                         molfile_qm_timestep_t *qm_ts) {
  int i;
  qmdata_t *data = (qmdata_t *)mydata;
  qm_timestep_t *cur_ts;

  if (data->num_frames_sent >= data->num_frames) {
    /* All frames were sent. */
    return MOLFILE_ERROR;
  }

  if (data->num_frames_sent == data->num_frames_read) {
    /* Read next coordinate block from file */
    /*fseek(data->file, data->filepos_array[data->num_frames_read], SEEK_SET);*/
    /*read_geom_block(data);*/

    printf("qcschemaplugin) Read frame %d\n", data->num_frames_read);
    data->num_frames_read++;
  }


  /* Copy the coordinates */
  for (i=0; i<natoms; i++) {
    ts->coords[3*i  ] = data->atoms[i].x;
    ts->coords[3*i+1] = data->atoms[i].y;
    ts->coords[3*i+2] = data->atoms[i].z; 
  }
  
  /*printf("qcschemaplugin) Sent frame %d\n", data->num_frames_sent); */
  data->num_frames_sent++;

  /* In MOLDEN the MOs are listed only for the last frame */
  if (data->num_frames_sent == data->num_frames) {
    
    /* get a convenient pointer to the current qm timestep */
    cur_ts = data->qm_timestep;

    read_molecular_orbitals(data);

    /* store the wave function and orbital energies */
    if (cur_ts != NULL && cur_ts->wave != NULL) {
      for (i=0; i<cur_ts->numwave; i++) {
        qm_wavefunction_t *wave = &cur_ts->wave[i];
        qm_ts->wave[i].type         = wave->type;
        qm_ts->wave[i].spin         = wave->spin;
        qm_ts->wave[i].excitation   = wave->exci;
        qm_ts->wave[i].multiplicity = wave->mult;
        qm_ts->wave[i].energy       = wave->energy;
        strncpy(qm_ts->wave[i].info, wave->info, MOLFILE_BUFSIZ);
        
        if (wave->wave_coeffs) {
          memcpy(qm_ts->wave[i].wave_coeffs, wave->wave_coeffs,
                 wave->num_orbitals*data->wavef_size*sizeof(float));
        }
        if (wave->orb_energies) {
          memcpy(qm_ts->wave[i].orbital_energies, wave->orb_energies,
                 wave->num_orbitals*sizeof(float));
        }
        if (wave->has_occup) {
          memcpy(qm_ts->wave[i].occupancies, wave->orb_occupancies,
                 wave->num_orbitals*sizeof(float));
        }
      }
    }

  }
  
  return MOLFILE_SUCCESS;
}
  

/*****************************************************
 *
 * Provide VMD with the sizes of the QM related
 * data structure arrays that need to be made
 * available.
 * Since we cannot determine the basis set meta data
 * without parsing the whole basis set section, we
 * read all basis set data here. The data is stored
 * in the qmdata_t structure for later use in
 * read_molden_rundata().
 *
 *****************************************************/
static int read_qcschema_metadata(void *mydata, 
    molfile_qm_metadata_t *metadata) {

  qmdata_t *data;
  jsondata_t* jsondata = (jsondata_t *)data->format_specific_data;
  data = (qmdata_t *)mydata;
  
  metadata->ncart = 0;
  metadata->nimag = 0;
  metadata->nintcoords = 0;

  metadata->have_sysinfo = 0;
  metadata->have_carthessian = 0;
  metadata->have_inthessian = 0;
  metadata->have_normalmodes = 0;

  metadata->num_basis_funcs = 0;
  metadata->num_basis_atoms = 0;
  metadata->num_shells = 0;
  metadata->wavef_size = 0;

  jsondata->coordsonly = 1;

  if (!jsondata->coordsonly) {
    /* Read the basis set */
    if (!get_basis(data)) return MOLFILE_ERROR; 

    /* orbital + basis set data */
    metadata->num_basis_funcs = data->num_basis_funcs;
    metadata->num_basis_atoms = data->num_basis_atoms;
    metadata->num_shells      = data->num_shells;
    metadata->wavef_size      = data->wavef_size;  
  }

  return MOLFILE_SUCCESS;
}


/******************************************************
 * 
 * Provide VMD with the static (i.e. non-trajectory)
 * data. That means we are filling the molfile_plugin
 * data structures.
 *
 ******************************************************/
static int read_qcschema_rundata(void *mydata, 
                               molfile_qm_t *qm_data) {
  qmdata_t *data = (qmdata_t *)mydata;
  int i;
  molfile_qm_hessian_t *hessian_data;
  molfile_qm_basis_t   *basis_data;
  molfile_qm_sysinfo_t *sys_data;

  if (!qm_data) return MOLFILE_ERROR;


  hessian_data = &qm_data->hess;
  basis_data   = &qm_data->basis;
  sys_data     = &qm_data->run;

  sys_data->num_electrons = data->num_electrons;
  sys_data->totalcharge = data->totalcharge;
  sys_data->runtype = data->runtype;


  /* Populate basis set data */
  if (data->num_basis_funcs) {
    for (i=0; i<data->num_basis_atoms; i++) {
      basis_data->num_shells_per_atom[i] = data->num_shells_per_atom[i];
      basis_data->atomic_number[i] = data->atomicnum_per_basisatom[i];
    }
    
    for (i=0; i<data->num_shells; i++) {
      basis_data->num_prim_per_shell[i] = data->num_prim_per_shell[i];
      basis_data->shell_types[i] = data->shell_types[i];
    }
    
    for (i=0; i<2*data->num_basis_funcs; i++) {
      basis_data->basis[i] = data->basis[i];
    }

    /* If we have MOs in the file we must provide the 
     * angular momentum exponents */
    if (data->angular_momentum) {
      for (i=0; i<3*data->wavef_size; i++) {
        basis_data->angular_momentum[i] = data->angular_momentum[i];
      }
    }
  }

  /* fill in molfile_qm_sysinfo_t */
  /*sys_data->runtype = data->runtype;
  sys_data->scftype = data->scftype;
  sys_data->nproc   = data->nproc;
  sys_data->num_electrons  = data->num_electrons;
  sys_data->totalcharge    = data->totalcharge;
  sys_data->num_occupied_A = data->num_occupied_A;
  sys_data->num_occupied_B = data->num_occupied_B;
  sys_data->status         = data->opt_status;
  */
  return MOLFILE_SUCCESS;
}


/**********************************************************
 *
 * close file and free memory
 *
 **********************************************************/
static void close_qcschema_read(void *mydata) {
  int i, j;
  qmdata_t *data = (qmdata_t *)mydata;
   
  /*fclose(data->file);*/
  free(data->atoms);
  free(data->basis);
  free(data->shell_types);
  free(data->atomicnum_per_basisatom);
  free(data->num_shells_per_atom);
  free(data->num_prim_per_shell);
  free(data->angular_momentum);

  if (data->basis_set) {
    for(i=0; i<data->num_basis_atoms; i++) {
      for (j=0; j<data->basis_set[i].numshells; j++) {
        free(data->basis_set[i].shell[j].prim);
      }
      free(data->basis_set[i].shell);
    } 
    free(data->basis_set);
  }

  free(data->format_specific_data);
  free(data->filepos_array);

  if (data->qm_timestep != NULL) {
    for (j=0; j<data->qm_timestep[0].numwave; j++) {
      free(data->qm_timestep[0].wave[j].wave_coeffs);
      free(data->qm_timestep[0].wave[j].orb_energies);
      free(data->qm_timestep[0].wave[j].orb_occupancies);
    }
    free(data->qm_timestep[0].wave);
    free(data->qm_timestep);
  } else {
    printf("close_qcschema_read(): NULL qm_timestep!\n");
  }

  free(data);
}


/* ####################################################### */
/*             End of API functions                        */
/* The following functions actually do the file parsing.   */
/* ####################################################### */


static int get_basis(qmdata_t *data) {

  return TRUE;
}


/******************************************************
 *
 * Populate the flat arrays containing the basis
 * set data.
 *
 ******************************************************/
static int fill_basis_arrays(qmdata_t *data) {
  int i, j, k;
  int shellcount = 0;
  int primcount = 0;

  float *basis;
  int *num_shells_per_atom;
  int *num_prim_per_shell;
  int *shell_types;
  int *atomicnum_per_basisatom;

  /* Count the total number of primitives which
   * determines the size of the basis array. */
  for(i=0; i<data->num_basis_atoms; i++) {
    for (j=0; j<data->basis_set[i].numshells; j++) {
      primcount += data->basis_set[i].shell[j].numprims;
    }
  }
  data->num_basis_funcs = primcount;

  /* reserve space for pointer to array containing basis
   * info, i.e. contraction coeficients and expansion 
   * coefficients; need 2 entries per basis function, i.e.
   * exponent and contraction coefficient; also,
   * allocate space for the array holding the orbital symmetry
   * information per primitive Gaussian.
   * Finally, initialize the arrays holding the number of 
   * shells per atom and the number of primitives per shell*/
  ALLOCATE(basis,                   float, 2*primcount);
  ALLOCATE(shell_types,             int,   data->num_shells);
  ALLOCATE(num_shells_per_atom,     int,   data->num_basis_atoms);
  ALLOCATE(num_prim_per_shell,      int,   data->num_shells);
  ALLOCATE(atomicnum_per_basisatom, int,   data->num_basis_atoms);



  /* store pointers in struct qmdata_t */
  data->basis = basis;
  data->shell_types = shell_types;
  data->num_shells_per_atom = num_shells_per_atom;
  data->num_prim_per_shell  = num_prim_per_shell;
  data->atomicnum_per_basisatom = atomicnum_per_basisatom;

  primcount = 0;
  for (i=0; i<data->num_basis_atoms; i++) {
    /* assign atomic number */
    data->basis_set[i].atomicnum = data->atoms[i].atomicnum;
    atomicnum_per_basisatom[i]   = data->atoms[i].atomicnum;

    num_shells_per_atom[i] = data->basis_set[i].numshells;

    for (j=0; j<data->basis_set[i].numshells; j++) {
      shell_types[shellcount]        = data->basis_set[i].shell[j].type;
      num_prim_per_shell[shellcount] = data->basis_set[i].shell[j].numprims;

      for (k=0; k<data->basis_set[i].shell[j].numprims; k++) {
        basis[2*primcount  ] = data->basis_set[i].shell[j].prim[k].exponent;
        basis[2*primcount+1] = data->basis_set[i].shell[j].prim[k].contraction_coeff;
        primcount++;
      }
      shellcount++;
    }
  } 

  return TRUE;
}


/**************************************************
 *
 * Convert shell type from char to int.
 * Note that SP_P shells are assigned in get_basis()
 *
 ************************************************ */
static int shelltype_int(char *type) {
  int shelltype;
  if      (!strcasecmp(type, "sp")) shelltype = SP_SHELL;
  else if (!strcasecmp(type, "s"))  shelltype = S_SHELL;
  else if (!strcasecmp(type, "p"))  shelltype = P_SHELL;
  else if (!strcasecmp(type, "d"))  shelltype = D_SHELL;
  else if (!strcasecmp(type, "f"))  shelltype = F_SHELL;
  else if (!strcasecmp(type, "g"))  shelltype = G_SHELL;
  else shelltype = UNK_SHELL;
  
  return shelltype;
}

static int count_orbitals(qmdata_t *data) {
  int nr;
  int num_wave_coeff=0;
  float orbenergy, occu;
  char spin[1024];
  qm_wavefunction_t *wave;
  jsondata_t *jsondata = (jsondata_t *)data->format_specific_data;
  int dummy1;
  float dummy2;

  /* Allocate memory for the qm_timestep frame */
  data->qm_timestep = (qm_timestep_t *)calloc(1, sizeof(qm_timestep_t));

  return TRUE;
}

static int read_molecular_orbitals(qmdata_t *data) {
  jsondata_t *jsondata = (jsondata_t *)data->format_specific_data;
  qm_wavefunction_t *wave;
  if (!data->qm_timestep || jsondata->coordsonly) return FALSE;
  return TRUE;
}

static int read_wave_coeffs(FILE *file, qm_wavefunction_t *wave) {
  return TRUE;
}

/*************************************************************
 *
 * plugin registration 
 *
 **************************************************************/
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "qcschema";
  plugin.prettyname = "QCSchema";
  plugin.author = "Mariano Spivak";
  plugin.majorv = 0;
  plugin.minorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "json,qcs";
  plugin.open_file_read = open_qcschema_read;
  plugin.read_structure = read_qcschema_structure;

  plugin.read_timestep_metadata    = read_timestep_metadata;
  plugin.read_timestep             = read_timestep;
  plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;

  plugin.read_qm_metadata = read_qcschema_metadata;
  plugin.read_qm_rundata  = read_qcschema_rundata;

  plugin.close_file_read = close_qcschema_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

