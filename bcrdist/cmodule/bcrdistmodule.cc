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
//   PyObject_HEAD
//   int num_strands;
//   vector<dsbcell> double_cells;
//   vector<ssbcell> single_cells;
//   vector<bcell*> ref_cells;
//   string name;
//   string dist_file;
// } py_bcell_vector;

template <typename celltype>
struct py_bcell_vector {
  static_assert(std::is_base_of<bcell,celltype>::value,
    "py_bcell_vector can only be used with a subclass of bcell"
  );
  PyObject_HEAD
  vector<celltype> cells;
  string name;
  string dist_file;
};

template <typename celltype>
static void py_bcell_vector_dealloc(py_bcell_vector<celltype> *self) {
  Py_TYPE(self)->tp_free((PyObject *) self);
}

template <typename celltype>
static PyObject* py_bcell_vector_load(py_bcell_vector<celltype>* self, PyObject* args) {
  // const char* filename;
  // if (!PyArg_ParseTuple(args, "s", &filename)) {
  //   return nullptr;
  // }
  
  ifstream ifile(self->name + ".bcrsavefile.tsv");
  if (!ifile.good()) {
    PyErr_SetString(PyExc_RuntimeError, "There is no save file");
    return nullptr;
  }
  string line;
  getline(ifile, line, ':');
  ifile >> line;
  if (self->name == "") {
    self->name = line;
  }
  getline(ifile, line, ':');
  ifile >> line;
  if (self->dist_file == "" and line != "<UNSET>") {
    self->dist_file = line;
  }
  
  getline(ifile, line);
  getline(ifile, line);
  
  getline(ifile, line);
  while (!ifile.eof()) {
    stringstream ss(line);
    self->cells.emplace_back(ss);
    getline(ifile, line);
  }
  
  Py_RETURN_NONE;
}

template <typename celltype>
PyObject* py_bcell_vector_append(py_bcell_vector<celltype>* self, PyObject* args) {
  PyErr_SetString(PyExc_RuntimeError, "This function is not supported for this type");
  return nullptr;
}

template <>
PyObject* py_bcell_vector_append<dsbcell>(py_bcell_vector<dsbcell>* self, PyObject* args) {
  PyObject* obj;
  if (!PyArg_ParseTuple(args, "O", &obj)) {
    obj = args;
    PyErr_Clear();
  }
  const char* strings[8];
  if (PyArg_ParseTuple(obj, "ssssssss", &strings[0], &strings[1], &strings[2],
  &strings[3], &strings[4], &strings[5], &strings[6], &strings[7])) {
    self->cells.emplace_back(
      strings[0],
      bcell_chain(strings[2], strings[3], strings[4]),
      bcell_chain(strings[5], strings[6], strings[7]),
      strings[1]
    );
  } else if (PyArg_ParseTuple(obj, "ssssss", &strings[0], &strings[1], &strings[2],
  &strings[3], &strings[4], &strings[5])) {
    PyErr_Clear();
    self->cells.emplace_back(
      strings[0],
      bcell_chain(strings[2], strings[3]),
      bcell_chain(strings[4], strings[5]),
      strings[1]
    );
  } else {
    return nullptr;
  }
  
  Py_RETURN_NONE;
}

template <typename celltype>
static int py_bcell_vector_type_init(py_bcell_vector<celltype> *self, PyObject *args) {
  const char* newname = nullptr;
  PyObject* obj;
  if (PyArg_ParseTuple(args, "|s", &newname)) {
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
  } else if (PyArg_ParseTuple(args, "O", &obj)) {
    PyErr_Clear();
    PyObject* iterator = PyObject_GetIter(obj);
    PyObject* item;
    
    if (iterator == nullptr) {
      PyErr_SetString(PyExc_RuntimeError, "Non iterable passed to constructor");
      return -1;
    }
    
    while ((item = PyIter_Next(iterator))) {
      PyObject* params = Py_BuildValue("O", item);
      py_bcell_vector_append<celltype>(self, params);
      Py_DECREF(params);
      Py_DECREF(item);
    }
    
    Py_DECREF(iterator);
    return 0;
  }
  
  return -1;
}

template <typename celltype>
static PyObject* py_bcell_vector_load_bd_data(py_bcell_vector<celltype> *self, PyObject *args) {
    const char* heavy;
    const char* light = nullptr;
    
    
    if (!PyArg_ParseTuple(args, "s|s", &heavy, &light)) {
        return NULL;
    }
    
    string heavy_str = heavy;
    
    if (light != nullptr) {
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
      
      if (!load_bd_data(heavy_str, light_str, self->cells)) {
        PyErr_SetString(PyExc_RuntimeError, "File not found");
        return nullptr;
      }
    } else {
      if (self->name == "") {
        self->name = heavy_str;
      }
      
      if (!load_bd_data(heavy_str, self->cells)) {
        PyErr_SetString(PyExc_RuntimeError, "File not found");
        return nullptr;
      }
    }
    
    Py_RETURN_NONE;
}

template <typename celltype>
static PyObject* py_bcell_vector_load_10x_data(py_bcell_vector<celltype> *self, PyObject *args) {
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }
    
    if (self->name == "") {
      self->name = filename;
    }
    
    if (!load_10x_data(filename, self->cells)) {
      PyErr_SetString(PyExc_RuntimeError, "File not found");
      return nullptr;
    }
    
    Py_RETURN_NONE;
}

template <typename celltype>
PyObject* py_bcell_vector_tolist(py_bcell_vector<celltype> *self) {
  PyErr_SetString(PyExc_RuntimeError, "This type of cell vector is not compatable");
  return NULL;
}

template <>
PyObject* py_bcell_vector_tolist<ssbcell>(py_bcell_vector<ssbcell>* self) {
  PyObject* list = PyList_New(self->cells.size());
  int i = 0;
  for (ssbcell& cell : self->cells) {
    PyObject* cell_list = PyList_New(5);
    string vals[] = {cell.id, cell.clonotype, cell.chain.cdr1, cell.chain.cdr2, cell.chain.cdr3 };
    int j = 0;
    for (string val : vals) {
      PyList_SET_ITEM(cell_list, j, PyUnicode_FromString(val.c_str()));
      j ++;
    }
    PyList_SET_ITEM(list, i, cell_list);
    i++;
  }
  return list;
}

template <>
PyObject* py_bcell_vector_tolist<dsbcell>(py_bcell_vector<dsbcell>* self) {
  PyObject* list = PyList_New(self->cells.size());
  // cout << list<< endl;
  int i = 0;
  for (dsbcell& cell : self->cells) {
    // cout << i << endl;
    // cout << &cell << endl;
    PyObject* cell_list = Py_BuildValue("ssssssss",
      cell.id.c_str(), cell.clonotype.c_str(),
      cell.heavy.cdr1.c_str(), cell.heavy.cdr2.c_str(), cell.heavy.cdr3.c_str(),
      cell.light.cdr1.c_str(), cell.light.cdr2.c_str(), cell.light.cdr3.c_str()
    );
    PyList_SET_ITEM(list, i, cell_list);
    i++;
  }
  return list;
}
      
template <typename celltype>
static PyObject* py_bcell_vector_generate_dist_matrix(py_bcell_vector<celltype> *self, PyObject *args) {
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
    
    save_dist_matrix(self->name + self->dist_file, self->cells);
    
    Py_RETURN_NONE;
}

template <typename celltype>
static PyObject* py_bcell_vector_get_dist_matrix(py_bcell_vector<celltype>* self, PyObject* args) {
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
  
  float* data = new float[num_cells*num_cells];
  for (int i = 0; i < num_cells*num_cells; i ++) {
    data[i] = -1;
  }
  
  read_dist_matrix(self->name + self->dist_file, data, num_cells);
  
  // for (int i = 0; i < num_cells; i ++) {
  //   tablerow row(&itable);
  //   for (int j = 0; j < num_cells-i; j ++) {
  //     if (num_cells-i == j+1) {
  //       data[(num_cells-i-1)*num_cells + j] = 0;
  //     } else {
  //       //cout << row.get(itable.headers[j+1]) << ' ' << std::atof(row.get(itable.headers[j+1]).c_str()) << endl;
  //       data[(num_cells-i-1)*num_cells + j] = std::atof(row.get(itable.headers[j+1]).c_str());
  //     }
  //   }
  // }
  // for (int i = 0; i < num_cells; i ++) {
  //   for (int j = i; j < num_cells; j ++) {
  //     data[i*num_cells + j] = data[j*num_cells + i];
  //   }
  // }
  
  
  
  const long int dims[] {num_cells, num_cells};
  PyObject* np_data = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, data);
  
  return PyTuple_Pack(2, np_data, cell_ids);
}

template <typename celltype>
static PyObject* py_bcell_vector_name(py_bcell_vector<celltype>* self) {
  return PyUnicode_FromString(self->name.c_str());
}

template <typename celltype>
static PyObject* py_bcell_vector_summary(py_bcell_vector<celltype>* self) {
  std::stringstream message;
  message << "  cell array: '" << self->name << "'" << endl;
  if (self->dist_file.length() > 0) {
    message << "  distance matrix at: '" << self->name << self->dist_file << "'" << endl;
  }
  ifstream ftest(self->name + ".bcrsavefile.tsv");
  if (ftest.good()) {
    message << "  save file at: '" << self->name << ".bcrsavefile.tsv'" << endl;
  }
  message << "  " << self->cells.size() << " cells" << endl;
  
  return PyUnicode_FromString(message.str().c_str());
}

template <typename celltype>
static PyObject* py_bcell_vector_repr(py_bcell_vector<celltype>* self) {
  string message = "<cbcrdist.bcellarray built in object\n";
  PyObject* summary_pystr = py_bcell_vector_summary(self);
  const char* summary = PyUnicode_AsUTF8(summary_pystr);
  Py_DECREF(summary_pystr);
  message += summary;
  message += '>';
  return PyUnicode_FromString(message.c_str());
}

template <typename celltype>
static PyObject* py_bcell_vector_save(py_bcell_vector<celltype>* self, PyObject* args) {
  // const char* filename;
  //
  // if (!PyArg_ParseTuple(args, "s", &filename)) {
  //   return nullptr;
  // }
  
  ofstream ofile(self->name + ".bcrsavefile.tsv");
  ofile << "# name:" << self->name << endl;
  if (self->dist_file == "") {
    ofile << "# dist_file:<UNSET>" << endl;
  } else {
    ofile << "# dist_file:" << self->dist_file << endl;
  }
  if constexpr(std::is_same<celltype,ssbcell>::value) {
    ofile << "id\tcdr1\tcdr2\tcdr3" << endl;
  } else if constexpr(std::is_same<celltype,dsbcell>::value) {
    ofile << "id\thcdr1\thcdr2\thcdr3\tlcdr1\tlcdr2\tlcdr3" << endl;
  }
  
  for (celltype& cell : self->cells) {
    cell.to_file(ofile);
    ofile << endl;
  }
  
  Py_RETURN_NONE;
}

template <typename celltype>
static PyObject* py_bcell_vector_size(py_bcell_vector<celltype>* self) {
  return PyLong_FromLong(self->cells.size());
};

template <typename celltype>
static PyMethodDef py_bcell_vector_methods[] = {
    {"gen_distmatrix", (PyCFunction) py_bcell_vector_generate_dist_matrix<celltype>, METH_VARARGS,
     "calculates the distance matrix and writes it to a file"},
    {"append", (PyCFunction) py_bcell_vector_append<celltype>, METH_VARARGS,
     "appends a cell in the format of a tuple"},
    {"loadBD", (PyCFunction) py_bcell_vector_load_bd_data<celltype>, METH_VARARGS,
     "loads data from two bd files"},
    {"load10x", (PyCFunction) py_bcell_vector_load_10x_data<celltype>, METH_VARARGS,
     "loads data from one 10x files"},
    {"distmatrix", (PyCFunction) py_bcell_vector_get_dist_matrix<celltype>, METH_NOARGS,
     "returns a numpy array of the distance matrix as well as all of the cell ids"},
    {"name", (PyCFunction) py_bcell_vector_name<celltype>, METH_NOARGS,
     "returns the name"},
    {"summary", (PyCFunction) py_bcell_vector_summary<celltype>, METH_NOARGS,
     "returns a string summary of the cell array"},
    {"save", (PyCFunction) py_bcell_vector_save<celltype>, METH_NOARGS,
     "returns a string summary of the cell array"},
    {"load", (PyCFunction) py_bcell_vector_load<celltype>, METH_NOARGS,
     "returns a string summary of the cell array"},
    {"tolist", (PyCFunction) py_bcell_vector_tolist<celltype>, METH_NOARGS,
     "outputs a list of the data"},
    {"size", (PyCFunction) py_bcell_vector_size<celltype>, METH_NOARGS,
     "outputs a list of the data"},
    {NULL}  /* Sentinel */
};

template <typename celltype>
static PyTypeObject py_bcell_vector_type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "cbcrdist.bcellarray",
  sizeof(py_bcell_vector<celltype>),
  0,
  (destructor) py_bcell_vector_dealloc<celltype>,
  0,
  0,
  0,
  0,
  (reprfunc) py_bcell_vector_repr<celltype>,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  "An array of double stranded bcells",
  0,
  0,
  0,
  0,
  0,
  0,
  py_bcell_vector_methods<celltype>,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  (initproc)py_bcell_vector_type_init<celltype>,
  0,
  PyType_GenericNew,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0
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
    
    
    if (PyType_Ready(&py_bcell_vector_type<dsbcell>) < 0) {
      return nullptr;
    }
    // if (PyType_Ready(&py_bcell_type) < 0) {
    //   return nullptr;
    // }
    
    PyObject* module = PyModule_Create(&cbcrdistmodule);
    import_array();
    
    Py_INCREF(&py_bcell_vector_type<dsbcell>);
    if (PyModule_AddObject(module, "bcellarray", (PyObject *) &py_bcell_vector_type<dsbcell>) < 0) {
      Py_DECREF(&py_bcell_vector_type<dsbcell>);
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
