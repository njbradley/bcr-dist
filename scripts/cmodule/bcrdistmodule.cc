#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "scripts/cell.cc"
#include "scripts/table.cc"
#include "scripts/data.cc"
#include "scripts/fileio.cc"


vector<bcell> cells;

static PyObject* py_load_bd_data(PyObject *self, PyObject *args) {
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) {
        return NULL;
    }
    
    load_bd_data(command, cells);
    
    Py_RETURN_NONE;
}

static PyObject* py_save_dist_matrix(PyObject *self, PyObject *args) {
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) {
        return NULL;
    }
    
    save_dist_matrix(command, cells);
    
    Py_RETURN_NONE;
}

static PyMethodDef bcrdistMethods[] = {
  {"load_bd_data", py_load_bd_data, METH_VARARGS, "Get a random number"},
  {"save_dist_matrix", py_save_dist_matrix, METH_VARARGS, "Get a random number"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef bcrdistmodule = {
    PyModuleDef_HEAD_INIT,
    "bcrdist",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    bcrdistMethods
};

PyMODINIT_FUNC PyInit_bcrdist(void) {
    load_persistant_data();
    return PyModule_Create(&bcrdistmodule);
}
