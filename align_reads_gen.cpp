#include <Python.h>
#include <iostream>

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#include "include/align_reads/aligner.hpp"

static align_reads::Aligner *gen = nullptr;

#define THREADS 20

static PyObject* initialize_cpp(PyObject *self, PyObject *args) {
    if (gen == nullptr) {
        PyObject *read_paths_list,  *item;
        const char* reads_path[2] = {nullptr, nullptr};
        
        if (!PyArg_ParseTuple(args, "O", &read_paths_list)) return NULL;
       
        if (PyList_Check(read_paths_list) && PyList_Size(read_paths_list) >= 1) {
            for (Py_ssize_t i = 0; i < PyList_Size(read_paths_list); i++) {
                item = PyList_GetItem(read_paths_list, i);
                if (!item) return NULL;
                item = PyUnicode_AsEncodedString(item, "UTF-8", "strict");
                if (!item) return NULL;
                reads_path[i] = PyBytes_AsString(item);
            }                                                   
        } else {
            return NULL;
        }
	
	    std::shared_ptr<thread_pool::ThreadPool> pool = std::make_shared<thread_pool::ThreadPool>(THREADS);         
        gen = new align_reads::Aligner(reads_path, pool, 15, 5, 0.005);           
    }
    Py_RETURN_NONE;    
}

// Module method definitions
static PyObject* generate_features_cpp(PyObject *self, PyObject *args) {
   
    gen->run();
	
    Py_RETURN_NONE;
}


static PyObject* cleanup_cpp(PyObject *self, PyObject *args) {
    delete gen;
    gen = nullptr;
    Py_RETURN_NONE;
}


static PyMethodDef align_reads_gen_methods[] = {
        {"generate_features", generate_features_cpp, METH_VARARGS,"Generate features for reads correction."},
        {"initialize", initialize_cpp, METH_VARARGS, "Initialize generator."}, 
        {"cleanup", cleanup_cpp, METH_VARARGS, "Clean up resources."},
        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef align_reads_gen_definition = {
        PyModuleDef_HEAD_INIT,
        "align_reads_gen",
        "Feature generation for correcting reads.",
        -1,
        align_reads_gen_methods
};


PyMODINIT_FUNC PyInit_align_reads_gen(void) {
    Py_Initialize();
    import_array();
    return PyModule_Create(&align_reads_gen_definition);
}
