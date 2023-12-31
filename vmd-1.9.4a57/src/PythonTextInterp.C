/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: PythonTextInterp.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.76 $       $Date: 2020/10/21 20:33:20 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  Python text interpreter
 ***************************************************************************/

#include "py_commands.h"
#include "Inform.h"
#include "PythonTextInterp.h"
#include "config.h"
#include "VMDApp.h"
#include "TextEvent.h"

#if defined(__APPLE__)
// use the Apple-provided Python framework
#include "Python/errcode.h"
#else
#include "errcode.h"
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_vmd(void);
#else
extern "C" PyMODINIT_FUNC PyInit_vmd(void);
#endif

static PyObject *cbdict = NULL;

static PyObject *add_callback(PyObject *, PyObject *args) {
  char *type;
  PyObject *temp;

  if (!PyArg_ParseTuple(args, "sO:add_callback", &type, &temp)) 
    return NULL;

  if (!PyCallable_Check(temp)) {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
    return NULL;
  }
  PyObject *cblist = PyDict_GetItemString(cbdict, type);
  if (!cblist) {
    PyErr_SetString(PyExc_KeyError, type);
    return NULL;
  }
  PyList_Append(cblist, temp);
  Py_XDECREF(temp); // PyList_Append does not steal a reference
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *del_callback(PyObject *, PyObject *args) {
  char *type;
  PyObject *temp;

  if (!PyArg_ParseTuple(args, "sO:del_callback", &type, &temp)) 
    return NULL;

  if (!PyCallable_Check(temp)) {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
    return NULL;
  }
  PyObject *cblist = PyDict_GetItemString(cbdict, type);
  if (!cblist) {
    PyErr_SetString(PyExc_KeyError, type);
    return NULL;
  }
  int ind = PySequence_Index(cblist, temp);
  if (ind >= 0) {
    PySequence_DelItem(cblist, ind);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static void call_callbacks(const char *type, PyObject *arglist) {
  PyObject *cblist = PyDict_GetItemString(cbdict, type);
  if (!cblist) {
    msgErr << "Internal error: callback list " << type << " does not exist."
           << sendmsg;
    return;
  }

  /* The GIL must be held when callbacks are invoked, but we don't have
   * the GIL when Tk (as opposed to Tkinter) events come in.  This ensures
   * that we have the GIL in either case.
   */
  PyGILState_STATE state = PyGILState_Ensure();

  for (int i=0; i<PyList_GET_SIZE(cblist); i++) {
    PyObject *obj = PyList_GET_ITEM(cblist, i);
    PyObject *result = PyEval_CallObject(obj, arglist);
    if (result == NULL) {
      PyErr_Print();
      PySequence_DelItem(cblist, i);
      i--;
    } else {
      Py_DECREF(result);
    }
  }
  Py_DECREF(arglist);

  PyGILState_Release(state);
}
  
static PyMethodDef CallbackMethods[] = {
  {"add_callback", (PyCFunction) add_callback, METH_VARARGS },
  {"del_callback", (PyCFunction) del_callback, METH_VARARGS },
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef vmdcallbacksdef = {
    PyModuleDef_HEAD_INIT,
    "vmdcallbacks",
    NULL,
    -1, // global state, no sub-interpreters
    CallbackMethods,
    NULL,
    NULL, // m_traverse gc traversal
    NULL, // m_clear gc clear
    NULL  // m_free gc free
};
#endif

PyObject* initvmdcallbacks(void) {
#if PY_MAJOR_VERSION >= 3
  PyObject *m = PyModule_Create(&vmdcallbacksdef);
#else
  PyObject *m = Py_InitModule("vmdcallbacks", CallbackMethods);
#endif
  PyObject *dict = PyDict_New();
  PyDict_SetItemString(dict, "display_update", PyList_New(0));
  PyDict_SetItemString(dict, "frame", PyList_New(0));
  PyDict_SetItemString(dict, "initialize_structure", PyList_New(0));
  PyDict_SetItemString(dict, "molecule", PyList_New(0));
  PyDict_SetItemString(dict, "pick_atom", PyList_New(0));
  PyDict_SetItemString(dict, "pick_event", PyList_New(0));
  PyDict_SetItemString(dict, "pick_value", PyList_New(0));
  PyDict_SetItemString(dict, "timestep", PyList_New(0));
  PyDict_SetItemString(dict, "trajectory", PyList_New(0));
  PyDict_SetItemString(dict, "userkey", PyList_New(0));
  PyObject_SetAttrString(m, "callbacks", dict);
  cbdict = dict;

  return m;
}

PythonTextInterp::PythonTextInterp(VMDApp *vmdapp) : app(vmdapp) {
  const char *oo_modules[] = {"Molecule", "Label", "Material", NULL};
  PyObject *vmdmodule;
  int retval, i;

  msgInfo << "Starting Python..." << sendmsg;

#if PY_MAJOR_VERSION >= 3
  // Import VMD builtin module automatically
  // Do this before Py_initialize called
  PyImport_AppendInittab("vmd", PyInit_vmd);
#endif

  // Do emit DeprecationWarnings
#if PY_MAJOR_VERSION >= 3
  PySys_AddWarnOption(L"default");
#else
  PySys_AddWarnOption((char*) "default");
#endif

#if 0 && PY_MAJOR_VERSION >= 3
  // Set program name used to find library path etc. Defaults to 'python',
  // must occur before initialization.
  Py_SetProgramName(Py_DecodeLocale(app->argv_m[0], NULL));
#endif
  Py_Initialize();

  // Some modules (like Tk) assume that os.argv has been initialized
#if PY_MAJOR_VERSION >= 3
  int argc = app->argc_m;
  wchar_t **wargv = (wchar_t **) PyMem_Malloc((argc+1) * sizeof(wchar_t*));
  memset(wargv, 0, (argc+1) * sizeof(wchar_t*));

  for (i=0; i<argc; i++) {
    wargv[i] = Py_DecodeLocale(app->argv_m[i], NULL);
  }

  PySys_SetArgv(argc, wargv);
#else
  PySys_SetArgv(app->argc_m, (char **)app->argv_m);
#endif
  set_vmdapp(app);

  // Set up the prompts
  PySys_SetObject((char*) "ps1", as_pystring(""));
  PySys_SetObject((char*) "ps2", as_pystring("... "));

  vmdmodule = PyImport_ImportModule("vmd");
  i = 0;
  while (py_initializers[i].name) {
    const char *name = py_initializers[i].name;

    PyObject *module = (*(py_initializers[i].initfunc))();
    i++;
    if (!module) {
      msgErr << "Failed to initialize builtin module " << name << sendmsg;
      PyErr_Print();
      continue;
    }
    retval = PyModule_AddObject(vmdmodule, CAST_HACK name, module);
    if (retval || PyErr_Occurred()) {
      msgErr << "Failed to import builtin module " << name << sendmsg;
      exit(1);
    }
  }

  // Now handle the three object-oriented classes
  for (const char **tmp = oo_modules; *tmp; tmp++) {
    PyObject *module = PyImport_ImportModule(*tmp);
    if (!module) {
      msgErr << "Failed to initialize object-oriented module " << *tmp << sendmsg;
      continue;
    }
    retval = PyModule_AddObject(vmdmodule, CAST_HACK *tmp, module);
    if (retval || PyErr_Occurred()) {
      msgErr << "Failed to initialize object-oriented module " << *tmp << sendmsg;
      continue;
    }
  }

  // Make all modules accessible in the default namespace
  if (!evalString("from vmd import *")) {
    msgErr << "Failed to import VMD python modules" << sendmsg;
    exit(1);
  }

  // have_tkinter and have_vmdcallback flags are set to zero if these calls
  // ever fail so that we don't fail over and over again and fill up the
  // screen with errors.
  have_tkinter = 1;
  in_tk = 0;
  needPrompt = 1;
}

PythonTextInterp::~PythonTextInterp() {
  Py_Finalize();
  msgInfo << "Done with Python." << sendmsg;
}

int PythonTextInterp::doTkUpdate() {
  // Don't recursively call into dooneevent - it makes Tkinter crash for
  // some infathomable reason.
  if (in_tk) return 0;
  if (have_tkinter) {
    in_tk = 1;
  int rc = 1;
// Python 3 just works for whatever reason when tkinter is imported.
// I'm not sure why this is...
#if PY_MAJOR_VERSION >= 3
    rc = evalString("import _tkinter\n");
#else
    rc = evalString(
      "import Tkinter\n"
      "while Tkinter.tkinter.dooneevent(Tkinter.tkinter.DONT_WAIT):\n"
      "\tpass\n"
    );
#endif
    in_tk = 0;
    if (rc) {
      return 1; // success
    }
    // give up
    have_tkinter = 0;
  }
  return 0;
}
  
void PythonTextInterp::doEvent() {
  // Call any display loop callbacks
  // abort if the call ever fails
  PyObject *arglist = Py_BuildValue("()");
  call_callbacks("display_update", arglist);

  if (needPrompt) {
    printf(">>> ");
    fflush(stdout);
    needPrompt = 0;
  }

  if (!vmd_check_stdin()) 
    return;

  int code = PyRun_InteractiveOne(stdin, "VMD");
  needPrompt = 1;
  if (code == E_EOF) {
    // Try to change to Tcl interpreter.  If that fails, UIText will
    // bounce us back to the Python interpreter again.
    app->textinterp_change("tcl");
  }
}

int PythonTextInterp::evalString(const char *s) {
  // evaluate the string in the interpreter
  // returns success.
  // XXX should print error message if there was one.
  return !PyRun_SimpleString(s);
}

int PythonTextInterp::evalFile(const char *s) {
  FILE *fid = fopen(s, "r");
  if (!fid) { 
    msgErr << "Error opening file '" << s << "'" << sendmsg;
    return FALSE;
  }
  int code = PyRun_SimpleFile(fid, "VMD");
  fclose(fid);
  return !code;
}
 
void PythonTextInterp::frame_cb(int molid, int frame) {
  PyObject *arglist = Py_BuildValue("(i,i)", molid, frame);
  call_callbacks("frame", arglist);
}

void PythonTextInterp::initialize_structure_cb(int molid, int code) {
  PyObject *arglist = Py_BuildValue("(i,i)", molid, code);
  call_callbacks("initialize_structure", arglist);
}

void PythonTextInterp::molecule_changed_cb(int molid, int code) {
  PyObject *arglist = Py_BuildValue("(i,i)", molid, code);
  call_callbacks("molecule", arglist);
}

void PythonTextInterp::pick_atom_cb(int mol, int atom, int key_shift_state, bool ispick) {
  PyObject *arglist = Py_BuildValue("(i,i,i)", mol, atom, key_shift_state);
  call_callbacks("pick_atom", arglist);
  if (ispick) {
    // if this is a user pick event, give it its own callback event
    // to discourage devs from inappropriately overloading all pick events
    PyObject *arglist = Py_BuildValue("(i)", 1);
    call_callbacks("pick_event", arglist);
  }
}

void PythonTextInterp::pick_value_cb(float val) {
  PyObject *arglist = Py_BuildValue("(f)", val);
  call_callbacks("pick_value", arglist);
}

void PythonTextInterp::timestep_cb(int id, int frame) {
  PyObject *arglist = Py_BuildValue("(i,i)", id, frame);
  call_callbacks("timestep", arglist);
}

void PythonTextInterp::trajectory_cb(int id, const char *name) {
  PyObject *arglist = Py_BuildValue("(i,s)", id, name);
  call_callbacks("trajectory", arglist);
}

void PythonTextInterp::python_cb(const char *cmd) {
  evalString(cmd);
}

void PythonTextInterp::userkey_cb(const char *keydesc) {
  PyObject *arglist = Py_BuildValue("(s)", keydesc);
  call_callbacks("userkey", arglist);
}

