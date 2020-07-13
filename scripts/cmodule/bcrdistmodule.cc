#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <numpy/arrayobject.h>
#include <structmember.h>
#include "scripts/cell.h"
#include "scripts/table.h"
#include "scripts/data.h"
#include "scripts/fileio.h"


vector<bcell_double> cells;



typedef struct {
    PyObject_HEAD
    bcell_double* cell;
} py_bcell_double;

static void py_bcell_double_dealloc(py_bcell_double *self) {
  delete self->cell;
  Py_TYPE(self)->tp_free((PyObject *) self);
}

// static PyObject * py_bcell_double_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
//   py_bcell_double *self;
//   self = (py_bcell_double *) type->tp_alloc(type, 0);
//   if (self != nullptr) {
//     const char* id;
//     const char* cdr1;
//     const char* cdr2;
//     const char* cdr3;
//     const char* vgene;
//     if (PyArg_ParseTuple(args, "sss", &id, &vgene, &cdr3)) {
//       self->cell = new bcell_double(id, vgene, cdr3);
//     } else if (PyArg_ParseTuple(args, "ssss", &id, &cdr1, &cdr2, &cdr3)) {
//       self->cell = new bcell_double(id, cdr1, cdr2, cdr3);
//     } else {
//       Py_DECREF(self);
//       return nullptr;
//     }
//   }
//   return (PyObject *) self;
// }

static PyTypeObject py_bcell_double_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "bcrdist.dsbcell",
  .tp_basicsize = sizeof(py_bcell_double),
  .tp_itemsize = 0,
  .tp_dealloc = (destructor) py_bcell_double_dealloc,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = "Custom objects",
  .tp_new = PyType_GenericNew,
};




typedef struct {
  PyObject_HEAD
  vector<bcell_double> cells;
  string name;
  string dist_file;
} py_bcell_double_vector;

static void py_bcell_double_vector_dealloc(py_bcell_double_vector *self) {
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static int py_bcell_double_vector_type_init(py_bcell_double_vector *self, PyObject *args) {
  self->cells = vector<bcell_double>();
  const char* newname = nullptr;
  if (!PyArg_ParseTuple(args, "|s", &newname)) {
    return -1;
  }
  if (newname != nullptr) {
    self->name = newname;
  }
  return 0;
}

static PyObject* py_bcell_double_vector_load_bd_data(py_bcell_double_vector *self, PyObject *args) {
    const char* heavy;
    const char* light;
    
    
    if (!PyArg_ParseTuple(args, "ss", &heavy, &light)) {
        return NULL;
    }
    
    string heavy_str = heavy;
    string light_str = light;
    
    if (self->name == "") {
      string newname;
      for (int i = 0; i < heavy_str.length() and i < light_str.length(); i ++) {
        if (heavy_str[i] == light_str[i]) {
          newname.push_back(heavy_str[i]);
        }
      }
      self->name = newname;
    }
    
    self->cells.clear();
    
    load_bd_data(heavy_str, light_str, self->cells);
    
    Py_RETURN_NONE;
}

static PyObject* py_bcell_double_vector_generate_dist_matrix(py_bcell_double_vector *self, PyObject *args) {
    const char* new_dist_file = nullptr;
    if (!PyArg_ParseTuple(args, "|s", &new_dist_file)) {
      return nullptr;
    }
    
    if (new_dist_file != nullptr) {
      self->dist_file = new_dist_file;
    }
    
    if (self->dist_file == "") {
      self->dist_file = self->name + "_dist.csv";
    }
    
    save_dist_matrix(self->dist_file, self->cells);
    
    Py_RETURN_NONE;
}

static PyObject* py_bcell_double_vector_get_dist_matrix(py_bcell_double_vector* self, PyObject* args) {
  if (self->dist_file == "") {
    py_bcell_double_vector_generate_dist_matrix(self, PyTuple_New(0));
  }
  
  itablestream itable(self->dist_file);
  int num_cells = itable.headers.size()-1;
  PyObject* cell_ids = PyList_New(num_cells);
  for (int i = 0; i < itable.headers.size()-1; i ++) {
    PyObject* id = PyUnicode_FromString(itable.headers[i+1].c_str());
    PyList_SET_ITEM(cell_ids, i, id);
  }
  
  double* data = new double[num_cells*num_cells];
  for (int i = 0; i < num_cells*num_cells; i ++) {
    data[i] = -1;
  }
  
  for (int i = 0; i < num_cells; i ++) {
    tablerow row(&itable);
    for (int j = 0; j < num_cells-i; j ++) {
      if (num_cells-i == j+1) {
        data[(num_cells-i-1)*num_cells + j] = 0;
      } else {
        //cout << row.get(itable.headers[j+1]) << ' ' << std::atof(row.get(itable.headers[j+1]).c_str()) << endl;
        data[(num_cells-i-1)*num_cells + j] = std::atof(row.get(itable.headers[j+1]).c_str());
      }
    }
  }
  for (int i = 0; i < num_cells; i ++) {
    for (int j = i; j < num_cells; j ++) {
      data[i*num_cells + j] = data[j*num_cells + i];
    }
  }
  const long int dims[] {num_cells, num_cells};
  PyObject* np_data = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data);
  
  return PyTuple_Pack(2, np_data, cell_ids);
}

static PyObject* py_bcell_double_vector_name(py_bcell_double_vector* self) {
  return PyUnicode_FromString(self->name.c_str());
}

static PyMethodDef py_bcell_double_vector_methods[] = {
    {"generate_dist_matrix", (PyCFunction) py_bcell_double_vector_generate_dist_matrix, METH_VARARGS,
     "calculates the distance matrix and writes it to a file"},
    {"load_bd_data", (PyCFunction) py_bcell_double_vector_load_bd_data, METH_VARARGS,
     "loads data from two bd files"},
    {"dist_matrix", (PyCFunction) py_bcell_double_vector_get_dist_matrix, METH_NOARGS,
     "returns a numpy array of the distance matrix as well as all of the cell ids"},
    {"name", (PyCFunction) py_bcell_double_vector_name, METH_NOARGS,
     "returns the name"},
    {NULL}  /* Sentinel */
};

static PyTypeObject py_bcell_double_vector_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "bcrdist.dsbcellarray",
  .tp_basicsize = sizeof(py_bcell_double_vector),
  .tp_itemsize = 0,
  .tp_dealloc = (destructor) py_bcell_double_vector_dealloc,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = "An array of double stranded bcells",
  .tp_methods = py_bcell_double_vector_methods,
  .tp_init = (initproc)py_bcell_double_vector_type_init,
  .tp_new = PyType_GenericNew,
};










static PyObject* py_load_bd_data(PyObject *self, PyObject *args) {
    const char* heavy;
    const char* light;

    if (!PyArg_ParseTuple(args, "ss", &heavy, &light)) {
        return NULL;
    }
    
    py_bcell_double_vector* cell_vector = PyObject_New(py_bcell_double_vector, &py_bcell_double_vector_type);
    PyObject* ret = PyObject_Init((PyObject*) cell_vector, &py_bcell_double_vector_type);
    load_bd_data(heavy, light, cell_vector->cells);
    return (PyObject*) cell_vector;
}

static PyObject* py_save_dist_matrix(PyObject *self, PyObject *args) {
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) {
        return NULL;
    }
    
    save_dist_matrix(command, cells);
    
    Py_RETURN_NONE;
}

static PyObject* py_get_dist_matrix(PyObject* self, PyObject* args) {
  const char* filename;
  
  if (!PyArg_ParseTuple(args, "s", &filename)) {
    return nullptr;
  }
  
  
  itablestream itable(filename);
  int num_cells = itable.headers.size()-1;
  PyObject* cell_ids = PyList_New(num_cells);
  for (int i = 0; i < itable.headers.size()-1; i ++) {
    PyObject* id = PyUnicode_FromString(itable.headers[i+1].c_str());
    PyList_SET_ITEM(cell_ids, i, id);
  }
  
  double* data = new double[num_cells*num_cells];
  for (int i = 0; i < num_cells*num_cells; i ++) {
    data[i] = -1;
  }
  
  for (int i = 0; i < num_cells; i ++) {
    tablerow row(&itable);
    for (int j = 0; j < num_cells-i; j ++) {
      if (num_cells-i == j+1) {
        data[(num_cells-i-1)*num_cells + j] = 0;
      } else {
        //cout << row.get(itable.headers[j+1]) << ' ' << std::atof(row.get(itable.headers[j+1]).c_str()) << endl;
        data[(num_cells-i-1)*num_cells + j] = std::atof(row.get(itable.headers[j+1]).c_str());
      }
    }
  }
  for (int i = 0; i < num_cells; i ++) {
    for (int j = i; j < num_cells; j ++) {
      data[i*num_cells + j] = data[j*num_cells + i];
    }
  }
  const long int dims[] {num_cells, num_cells};
  PyObject* np_data = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data);
  
  return PyTuple_Pack(2, np_data, cell_ids);
}
    

static PyMethodDef bcrdistMethods[] = {
  {"load_bd_data", py_load_bd_data, METH_VARARGS, "Get a random number"},
  {"save_dist_matrix", py_save_dist_matrix, METH_VARARGS, "Get a random number"},
  {"get_dist_matrix", py_get_dist_matrix, METH_VARARGS, "Get a random number"},
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
    
    if (PyType_Ready(&py_bcell_double_vector_type) < 0) {
      return nullptr;
    }
    if (PyType_Ready(&py_bcell_double_type) < 0) {
      return nullptr;
    }
    
    PyObject* module = PyModule_Create(&bcrdistmodule);
    import_array();
    
    Py_INCREF(&py_bcell_double_vector_type);
    if (PyModule_AddObject(module, "dsbcellarray", (PyObject *) &py_bcell_double_vector_type) < 0) {
      Py_DECREF(&py_bcell_double_vector_type);
      Py_DECREF(module);
      return nullptr;
    }
    Py_INCREF(&py_bcell_double_type);
    if (PyModule_AddObject(module, "dsbcell", (PyObject *) &py_bcell_double_type) < 0) {
      Py_DECREF(&py_bcell_double_type);
      Py_DECREF(module);
      return nullptr;
    }
    
    return module;
}
