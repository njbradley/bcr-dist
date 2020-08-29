#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <numpy/arrayobject.h>
#include <structmember.h>
#include "scripts/cell.h"
#include "scripts/table.h"
#include "scripts/data.h"
#include "scripts/fileio.h"


vector<dsbcell> cells;



// typedef struct {
//     PyObject_HEAD
//     dsbcell* cell;
// } py_dsbcell;
//
// static void py_dsbcell_dealloc(py_dsbcell *self) {
//   delete self->cell;
//   Py_TYPE(self)->tp_free((PyObject *) self);
// }

// static PyObject * py_dsbcell_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
//   py_dsbcell *self;
//   self = (py_dsbcell *) type->tp_alloc(type, 0);
//   if (self != nullptr) {
//     const char* id;
//     const char* cdr1;
//     const char* cdr2;
//     const char* cdr3;
//     const char* vgene;
//     if (PyArg_ParseTuple(args, "sss", &id, &vgene, &cdr3)) {
//       self->cell = new dsbcell(id, vgene, cdr3);
//     } else if (PyArg_ParseTuple(args, "ssss", &id, &cdr1, &cdr2, &cdr3)) {
//       self->cell = new dsbcell(id, cdr1, cdr2, cdr3);
//     } else {
//       Py_DECREF(self);
//       return nullptr;
//     }
//   }
//   return (PyObject *) self;
// }
//
// static PyTypeObject py_dsbcell_type = {
//   PyVarObject_HEAD_INIT(NULL, 0)
//   .tp_name = "cbcrdist.dsbcell",
//   .tp_basicsize = sizeof(py_dsbcell),
//   .tp_itemsize = 0,
//   .tp_dealloc = (destructor) py_dsbcell_dealloc,
//   .tp_flags = Py_TPFLAGS_DEFAULT,
//   .tp_doc = "Custom objects",
//   .tp_new = PyType_GenericNew,
// };




typedef struct {
  PyObject_HEAD
  int num_strands;
  vector<dsbcell> double_cells;
  vector<ssbcell> single_cells;
  string name;
  string dist_file;
} py_bcell_vector;

static void py_bcell_vector_dealloc(py_bcell_vector *self) {
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject* py_bcell_vector_load(py_bcell_vector* self, PyObject* args) {
  // const char* filename;
  // if (!PyArg_ParseTuple(args, "s", &filename)) {
  //   return nullptr;
  // }
  
  ifstream ifile(self->name + ".bcrsavefile.tsv");
  string line;
  getline(ifile, line, ':');
  ifile >> line;
  if (self->name == "") {
    self->name = line;
  }
  getline(ifile, line, ':');
  ifile >> line;
  if (self->dist_file == "") {
    self->dist_file = line;
  }
  
  getline(ifile, line, ':');
  ifile >> self->num_strands;
  
  getline(ifile, line);
  getline(ifile, line);
  
  if (self->num_strands == 1) {
    getline(ifile, line);
    while (!ifile.eof()) {
      stringstream ss(line);
      self->single_cells.emplace_back(ss);
      getline(ifile, line);
    }
  } else if (self->num_strands == 2) {
    getline(ifile, line);
    while (!ifile.eof()) {
      stringstream ss(line);
      self->double_cells.emplace_back(ss);
      getline(ifile, line);
    }
  }
  
  Py_RETURN_NONE;
}

static int py_bcell_vector_type_init(py_bcell_vector *self, PyObject *args) {
  self->num_strands = 0;
  const char* newname = nullptr;
  if (!PyArg_ParseTuple(args, "|s", &newname)) {
    return -1;
  }
  if (newname != nullptr) {
    self->name = newname;
    
    ifstream ifile(self->name + ".bcrsavefile.tsv");
    if (ifile.good()) {
      PyObject* newargs = PyTuple_New(0);
      py_bcell_vector_load(self, newargs);
      Py_DECREF(newargs);
    }
  }
  return 0;
}

static PyObject* py_bcell_vector_load_dekosky_data(py_bcell_vector *self, PyObject *args) {
  const char* path;
  
  if (!PyArg_ParseTuple(args, "s", &path)) {
    return nullptr;
  }
  
  self->num_strands = 2;
  if (self->name == "") {
    self->name = path;
  }
  load_dekosky_data(path, self->double_cells);
  
  Py_RETURN_NONE;
}

static PyObject* py_bcell_vector_load_bd_data(py_bcell_vector *self, PyObject *args) {
    const char* heavy;
    const char* light = nullptr;
    
    
    if (!PyArg_ParseTuple(args, "s|s", &heavy, &light)) {
        return NULL;
    }
    
    string heavy_str = heavy;
    
    if (light != nullptr) {
      self->num_strands = 2;
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
      
      load_bd_data(heavy_str, light_str, self->double_cells);
    } else {
      self->num_strands = 1;
      if (self->name == "") {
        self->name = heavy_str;
      }
      
      load_bd_data(heavy_str, self->single_cells);
    }
    
    Py_RETURN_NONE;
}

static PyObject* py_bcell_vector_load_10x_data(py_bcell_vector *self, PyObject *args) {
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }
    
    if (self->name == "") {
      self->name = filename;
    }
    
    self->num_strands = 2;
    load_10x_data(filename, self->double_cells);
    
    Py_RETURN_NONE;
}

static PyObject* py_bcell_vector_generate_dist_matrix(py_bcell_vector *self, PyObject *args) {
    const char* new_dist_file = nullptr;
    if (!PyArg_ParseTuple(args, "|s", &new_dist_file)) {
      return nullptr;
    }
    
    if (new_dist_file != nullptr) {
      self->dist_file = new_dist_file;
    }
    
    if (self->dist_file == "") {
      self->dist_file = ".bcrdistmatrix.tsv";
    }
    
    if (self->num_strands == 2) {
      save_dist_matrix(self->name + self->dist_file, self->double_cells);
    } else if (self->num_strands == 1) {
      save_dist_matrix(self->name + self->dist_file, self->single_cells);
    }
    
    Py_RETURN_NONE;
}

static PyObject* py_bcell_vector_get_dist_matrix(py_bcell_vector* self, PyObject* args) {
  if (self->dist_file == "") {
    PyObject* newargs = PyTuple_New(0);
    py_bcell_vector_generate_dist_matrix(self, newargs);
    Py_DECREF(newargs);
  }
  
  itablestream itable(self->name + self->dist_file);
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

static PyObject* py_bcell_vector_name(py_bcell_vector* self) {
  return PyUnicode_FromString(self->name.c_str());
}

static PyObject* py_bcell_vector_summary(py_bcell_vector* self) {
  std::stringstream message;
  //message << " Summary of bcell array object:" << endl;
  if (self->num_strands == 0) {
    message << "  Object uninitialized..." << endl;
  } else if (self->num_strands == 1) {
    message << "  Single stranded cell array" << endl;
    message << "  " << self->single_cells.size() << " cells" << endl;
  } else if (self->num_strands == 2) {
    message << "  Double stranded cell array" << endl;
    message << "  " << self->double_cells.size() << " cells" << endl;
  }
  
  return PyUnicode_FromString(message.str().c_str());
}

static PyObject* py_bcell_vector_repr(py_bcell_vector* self) {
  string message = "<cbcrdist.bcellarray built in object\n";
  PyObject* summary_pystr = py_bcell_vector_summary(self);
  const char* summary = PyUnicode_AsUTF8(summary_pystr);
  Py_DECREF(summary_pystr);
  message += summary;
  message += '>';
  return PyUnicode_FromString(message.c_str());
}

static PyObject* py_bcell_vector_save(py_bcell_vector* self, PyObject* args) {
  // const char* filename;
  //
  // if (!PyArg_ParseTuple(args, "s", &filename)) {
  //   return nullptr;
  // }
  
  ofstream ofile(self->name + ".bcrsavefile.tsv");
  ofile << "# name:" << self->name << endl;
  ofile << "# dist_file:" << self->dist_file << endl;
  ofile << "# num_strands:" << self->num_strands << endl;
  if (self->num_strands == 1) {
    ofile << "id\tcdr1\tcdr2\tcdr3" << endl;
    for (ssbcell& cell : self->single_cells) {
      cell.to_file(ofile);
      ofile << endl;
    }
  } else if (self->num_strands == 2) {
    ofile << "id\thcdr1\thcdr2\thcdr3\tlcdr1\tlcdr2\tlcdr3" << endl;
    for (dsbcell& cell : self->double_cells) {
      cell.to_file(ofile);
      ofile << endl;
    }
  }
  
  Py_RETURN_NONE;
}

static PyMethodDef py_bcell_vector_methods[] = {
    {"generate_dist_matrix", (PyCFunction) py_bcell_vector_generate_dist_matrix, METH_VARARGS,
     "calculates the distance matrix and writes it to a file"},
    {"loadBD", (PyCFunction) py_bcell_vector_load_bd_data, METH_VARARGS,
     "loads data from two bd files"},
    {"load10x", (PyCFunction) py_bcell_vector_load_10x_data, METH_VARARGS,
     "loads data from one 10x files"},
    {"distmatrix", (PyCFunction) py_bcell_vector_get_dist_matrix, METH_NOARGS,
     "returns a numpy array of the distance matrix as well as all of the cell ids"},
    {"name", (PyCFunction) py_bcell_vector_name, METH_NOARGS,
     "returns the name"},
    {"summary", (PyCFunction) py_bcell_vector_summary, METH_NOARGS,
     "returns a string summary of the cell array"},
    {"save", (PyCFunction) py_bcell_vector_save, METH_NOARGS,
     "returns a string summary of the cell array"},
    {"load", (PyCFunction) py_bcell_vector_load, METH_NOARGS,
     "returns a string summary of the cell array"},
    {"loaddekosky", (PyCFunction) py_bcell_vector_load_dekosky_data, METH_VARARGS,
     "loads dekosky data" },
    {NULL}  /* Sentinel */
};

static PyTypeObject py_bcell_vector_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "cbcrdist.bcellarray",
  .tp_basicsize = sizeof(py_bcell_vector),
  .tp_itemsize = 0,
  .tp_dealloc = (destructor) py_bcell_vector_dealloc,
  .tp_repr = (reprfunc) py_bcell_vector_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = "An array of double stranded bcells",
  .tp_methods = py_bcell_vector_methods,
  .tp_init = (initproc)py_bcell_vector_type_init,
  .tp_new = PyType_GenericNew,
};










// static PyObject* py_load_bd_data(PyObject *self, PyObject *args) {
//     // const char* heavy;
//     // const char* light;
//     //
//     // if (!PyArg_ParseTuple(args, "ss", &heavy, &light)) {
//     //     return NULL;
//     // }
//     //
//     PyObject* argList = PyTuple_New(0);
//     PyObject* obj = PyObject_CallObject((PyObject *) &py_bcell_vector_type, argList);
//     Py_DECREF(argList);
//
//     PyObject* result = PyObject_CallMethod(obj, "load_bd_data", "ss", args);
//     if (result == nullptr) {
//       return nullptr;
//     }
//
//     return obj;
// }


string get_path() {
  PyObject* mainmodule = PyImport_AddModule("__main__");
  PyRun_SimpleString("import os; basepath = os.path.dirname(os.path.abspath(__file__))");
  PyObject* pres = PyObject_GetAttrString(mainmodule, "basepath");
  const char* path = PyUnicode_AsUTF8(pres);
  string result = path;
  result = result + '/';
  Py_DECREF(pres);
  Py_DECREF(mainmodule);
  return result;
}


PyObject* py_init(PyObject* self, PyObject* args) {
  const char* base_path;
  if(!PyArg_ParseTuple(args, "s", &base_path)) {
    return nullptr;
  }
  
  load_persistant_data(base_path);
  
  Py_RETURN_NONE;
}

static PyMethodDef cbcrdistMethods[] = {
  // {"load_bd_data", py_load_bd_data, METH_VARARGS, "Get a random number"},
  {"init", py_init, METH_VARARGS, "loads in persistant data"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef cbcrdistmodule = {
    PyModuleDef_HEAD_INIT,
    "cbcrdist",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    cbcrdistMethods
};

PyMODINIT_FUNC PyInit_cbcrdist(void) {
    
    //string path = get_path();
    //load_persistant_data();
    
    
    if (PyType_Ready(&py_bcell_vector_type) < 0) {
      return nullptr;
    }
    // if (PyType_Ready(&py_bcell_type) < 0) {
    //   return nullptr;
    // }
    
    PyObject* module = PyModule_Create(&cbcrdistmodule);
    import_array();
    
    Py_INCREF(&py_bcell_vector_type);
    if (PyModule_AddObject(module, "bcellarray", (PyObject *) &py_bcell_vector_type) < 0) {
      Py_DECREF(&py_bcell_vector_type);
      Py_DECREF(module);
      return nullptr;
    }
    // Py_INCREF(&py_bcell_type);
    // if (PyModule_AddObject(module, "dsbcell", (PyObject *) &py_bcell_type) < 0) {
    //   Py_DECREF(&py_bcell_type);
    //   Py_DECREF(module);
    //   return nullptr;
    // }
    
    return module;
}
