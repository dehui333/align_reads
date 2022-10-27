#include <Python.h>
#include <iostream>
#include <string.h>

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#include "align_reads/Generator.hpp"

using namespace align_reads;

// -Need to persist over multiple python level calls
// -Too large for multiprocessing, do parallelizing in C++ layer
static align_reads::Generator *gen = nullptr;

inline static std::string convert_string(PyObject *string)
{
    string = PyUnicode_AsEncodedString(string, "UTF-8", "strict");
    return {PyBytes_AsString(string)};
}

// convert a python list of strings into a c++ vector of strings
inline static std::vector<std::string> convert_strings(PyObject *strings_list)
{
    std::vector<std::string> result;
    auto len = PyList_Size(strings_list);
    PyObject *item;
    for (Py_ssize_t i = 0; i < len; i++)
    {
        item = PyList_GetItem(strings_list, i);
        item = PyUnicode_AsEncodedString(item, "UTF-8", "strict");
        result.emplace_back(PyBytes_AsString(item));
    }

    return result;
}

// convert a list of strings at position i of a python list
inline static std::vector<std::string> convert_strings_in_pylist(PyObject *list, Py_ssize_t i)
{
    PyObject *item;
    item = PyList_GetItem(list, i);
    return convert_strings(item);
}

// convert a string at position i of a python list
inline static std::string convert_string_in_pylist(PyObject *list, Py_ssize_t i)
{
    PyObject *item;
    item = PyList_GetItem(list, i);
    return convert_string(item);
}

/* Take in inputs and initialize the Generator object
Expected input format on python level:
Currently, just let it be a list of arbitrary length

*/
static PyObject *initialize_cpp(PyObject *self, PyObject *args)
{

    if (gen == nullptr)
    {
        // will be a list of list of strings
        PyObject *paths_list;
        int num_threads;
        int window_length;
        int debug_print;
        const char *truth_path;
        // Arguments: list, number threads, window_length, debug_print(bool)
        if (!PyArg_ParseTuple(args, "Osiii", &paths_list, &truth_path, &num_threads, &window_length, &debug_print))
        {
            return NULL;
        }
        // confirm is list
        if (!PyList_Check(paths_list))
        {
            return NULL;
        }

        // Get length
        auto len = PyList_Size(paths_list);
        if (len == 0)
        {
            return NULL;
        }

        std::shared_ptr<thread_pool::ThreadPool> pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

        if (len == 1)
        {
            // Not split into haplotypes
            auto paths = convert_strings_in_pylist(paths_list, 0);
            gen = new Generator(paths, pool, window_length, debug_print);
        }
        else
        {
        }
        if (strcmp(truth_path, "") != 0)
        {
            std::string truth_path_string{truth_path};
            gen->index_truth(truth_path_string);
        }
    }

    Py_RETURN_NONE;
}

/*
static PyObject* initialize_cpp(PyObject *self, PyObject *args) {
    if (gen == nullptr) {
        PyObject *read_paths_list, *hap_paths_list, *item;
        const char* reads_path[2] = {nullptr, nullptr};
        const char* haplotypes_path[2] = {nullptr, nullptr};

        if (!PyArg_ParseTuple(args, "OO", &read_paths_list, &hap_paths_list)) return NULL;

        if (PyList_Check(read_paths_list) && PyList_Size(read_paths_list) >= 1 && PyList_Check(hap_paths_list)) {
            for (Py_ssize_t i = 0; i < PyList_Size(read_paths_list); i++) {
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

    std::shared_ptr<thread_pool::ThreadPool> pool = std::make_shared<thread_pool::ThreadPool>(THREADS);
        if (haplotypes_path[0] == nullptr || haplotypes_path[1] == nullptr) {
            gen = new align_reads::Aligner(reads_path, pool, 15, 5, 0.001, nullptr);
        } else {
            gen = new align_reads::Aligner(reads_path, pool, 15, 5, 0.001, haplotypes_path);

        }
    }
    Py_RETURN_NONE;
}
*/
static PyObject *print_window_cpp(PyObject *self, PyObject *args)
{
    int width_idx;
    int len;
    if (!PyArg_ParseTuple(args, "ii", &width_idx, &len))
    {
        return NULL;
    }
    gen->print_window(width_idx, len);
    Py_RETURN_NONE;
}

// return features. Maybe modified to fit needs
static PyObject *generate_features_cpp(PyObject *self, PyObject *args)
{
    /*align_reads::Data data = gen->next();
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

    return return_tuple;*/
    // Get a batch of feature data

    std::vector<std::vector<PyObject *>> data = gen->produce_data();

    // Creates tuple to return to python level
    PyObject *return_tuple = PyTuple_New(data.size());
    for (std::uint32_t i = 0; i < data.size(); i++)
    {
        auto &pyObjects = data[i];
        PyObject *python_list = PyList_New(pyObjects.size());
        for (std::uint32_t j = 0; j < pyObjects.size(); j++)
        {
            PyList_SetItem(python_list, j, pyObjects[j]);
        }
        PyTuple_SetItem(return_tuple, i, python_list);
    }

    return return_tuple;
}

// free memory.
static PyObject *cleanup_cpp(PyObject *self, PyObject *args)
{
    delete gen;
    gen = nullptr;
    Py_RETURN_NONE;
}

static PyObject *has_next_cpp(PyObject *self, PyObject *args)
{
    if (gen->has_next())
    {
        Py_RETURN_TRUE;
    }
    else
    {
        Py_RETURN_FALSE;
    }
}

static PyMethodDef align_reads_gen_methods[] = {
    {"generate_features", generate_features_cpp, METH_VARARGS, "Generate features for reads correction."},
    {"initialize", initialize_cpp, METH_VARARGS, "Initialize generator."},
    {"cleanup", cleanup_cpp, METH_VARARGS, "Clean up resources."},
    {"print_window", print_window_cpp, METH_VARARGS, "debug printting."},
    {"has_next", has_next_cpp, METH_VARARGS, "check if there's another target read."},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef align_reads_gen_definition = {
    PyModuleDef_HEAD_INIT,
    "align_reads_gen",
    "Feature generation for correcting reads.",
    -1,
    align_reads_gen_methods};

PyMODINIT_FUNC PyInit_align_reads_gen(void)
{
    Py_Initialize();
    import_array();
    return PyModule_Create(&align_reads_gen_definition);
}
