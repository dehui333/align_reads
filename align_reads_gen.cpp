#include <Python.h>
#include <iostream>

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#include "include/align_reads/aligner.hpp"

static align_reads::Aligner *gen = nullptr;

// Module method definitions
static PyObject* generate_features_cpp(PyObject *self, PyObject *args) {
    
    align_reads::Data data = gen->next();
	
	// Creates tuple to return to python level
	PyObject* return_tuple = PyTuple_New(2);
	
	// Create lists of matrices 
	PyObject* X_list = PyList_New(data.X.size());
	PyObject* Y_list = PyList_New(data.Y.size());
    
    
	// fill in lists
	for (std::uint32_t i = 0; i < data.X.size(); i++) {
		PyList_SetItem(X_list, i, data.X[i]);				
	}
	for (std::uint32_t i = 0; i < data.Y.size(); i++) {
		PyList_SetItem(Y_list, i, data.Y[i]);				
	}
	
	PyTuple_SetItem(return_tuple, 0, X_list);
    PyTuple_SetItem(return_tuple, 1, Y_list);
	
    return return_tuple;
}

static PyObject* initialize_cpp(PyObject *self, PyObject *args) {
    if (gen == nullptr) {
        PyObject *read_paths_list, *hap_paths_list, *item;
        const char* reads_path[2] = {nullptr, nullptr};
        const char* haplotypes_path[2] = {nullptr, nullptr};
        
        if (!PyArg_ParseTuple(args, "OO", &read_paths_list, &hap_paths_list)) return NULL;
        
        if (PyList_Check(read_paths_list) && PyList_Size(read_paths_list) == 2 && PyList_Check(hap_paths_list)) {
            for (Py_ssize_t i = 0; i < 2; i++) {
                item = PyList_GetItem(read_paths_list, i);
                if (!item) return NULL;
                item = PyUnicode_AsEncodedString(item, "UTF-8", "strict");
                if (!item) return NULL;
                reads_path[i] = PyBytes_AsString(item);
            }
            
            for (Py_ssize_t i = 0; i < PyList_Size(hap_paths_list); i++) {
                item = PyList_GetItem(hap_paths_list, i);
                if (!item) return NULL;
                item = PyUnicode_AsEncodedString(item, "UTF-8", "strict");
                if (!item) return NULL;
                haplotypes_path[i] = PyBytes_AsString(item);           
            }                 
            
        } else {
            return NULL;
        }
        //align_reads::Aligner gen {reads_path, 3, 15, 5, 0.001, haplotypes_path};
        gen = new align_reads::Aligner(reads_path, 3, 15, 5, 0.001, haplotypes_path);
    }
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