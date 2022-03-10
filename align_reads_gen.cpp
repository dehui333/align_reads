#include <Python.h>
#include <iostream>

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#include "include/align_reads/aligner.hpp"

// Module method definitions
static PyObject* generate_features_cpp(PyObject *self, PyObject *args) {
    const char* reads_path[2] = {"test_data/fake_reads0.fasta", "test_data/fake_reads1.fasta"};
    const char* haplotypes_path[2] = {"test_data/fake_haplotype0.fasta", "test_data/fake_haplotype1.fasta"};
    align_reads::Aligner gen {reads_path, 3, 15, 5, 0.001, haplotypes_path};
    std::cout <<"HI" << std::endl;
    return gen.test();
}

static PyMethodDef align_reads_gen_methods[] = {
        {
                "generate_features", generate_features_cpp, METH_VARARGS,
                "Generate features for reads correction."
        },
        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef align_reads_gen_definition = {
        PyModuleDef_HEAD_INIT,
        "align_gen",
        "Feature generation for correcting reads.",
        -1,
        align_reads_gen_methods
};


PyMODINIT_FUNC PyInit_align_reads_gen(void) {
    Py_Initialize();
    import_array();
    return PyModule_Create(&align_reads_gen_definition);
}