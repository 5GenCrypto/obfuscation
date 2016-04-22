#include <Python.h>
#include "obfuscator.h"
#include "pyutils.h"

static bool g_verbose;

static void
obf_clear_wrapper(PyObject *self)
{
    obf_clear((obf_state_t *) PyCapsule_GetPointer(self, NULL));
}

static PyObject *
obf_init_wrapper(PyObject *self, PyObject *args)
{
    long secparam, kappa, nzs, nthreads, ncores;
    enum mmap_e type;
    char *dir;
    obf_state_t *s = NULL;

    if (!PyArg_ParseTuple(args, "sllllll", &dir, &type, &secparam, &kappa,
                          &nzs, &nthreads, &ncores)) {
        PyErr_SetString(PyExc_RuntimeError, "unable to parse input");
        return NULL;
    }
    if (secparam <= 0 || kappa <= 0 || nzs <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "invalid input");
        return NULL;
    }

    s = obf_init(type, dir, secparam, kappa, nzs, nthreads, ncores, g_verbose);
    if (s == NULL)
        return NULL;

    return PyCapsule_New((void *) s, NULL, obf_clear_wrapper);
}

static PyObject *
obf_encode_layer_wrapper(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_zero_pows, *py_one_pows, *py_zero_ms, *py_one_ms;
    long inp, idx, nrows, ncols, rflag;
    ssize_t length;
    int *zero_pows, *one_pows;
    fmpz_mat_t zero, one;
    obf_state_t *s;

    if (!PyArg_ParseTuple(args, "OlllllOOOO", &py_state, &idx, &nrows, &ncols,
                          &inp, &rflag, &py_zero_pows, &py_one_pows, &py_zero_ms,
                          &py_one_ms))
        return NULL;

    s = (obf_state_t *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    length = PyList_Size(py_zero_pows);
    if (PyList_Size(py_one_pows) != length)
        return NULL;

    zero_pows = (int *) calloc(length, sizeof(int));
    one_pows = (int *) calloc(length, sizeof(int));

    for (ssize_t i = 0; i < length; ++i) {
        zero_pows[i] = PyLong_AsLong(PyList_GetItem(py_zero_pows, i));
        one_pows[i] = PyLong_AsLong(PyList_GetItem(py_one_pows, i));
    }

    fmpz_mat_init(zero, nrows, ncols);
    fmpz_mat_init(one, nrows, ncols);

    for (long i = 0; i < nrows; ++i) {
        for (long j = 0; j < ncols; ++j) {
            py_to_fmpz(fmpz_mat_entry(zero, i, j),
                       PyList_GetItem(PyList_GetItem(py_zero_ms, i), j));
            py_to_fmpz(fmpz_mat_entry(one, i, j),
                       PyList_GetItem(PyList_GetItem(py_one_ms, i), j));
        }
    }

    obf_encode_layer(s, idx, inp, nrows, ncols,
                     (encode_layer_randomization_flag_t) rflag, zero_pows,
                     one_pows, zero, one);

    fmpz_mat_clear(zero);
    fmpz_mat_clear(one);

    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate_wrapper(PyObject *self, PyObject *args)
{
    char *dir = NULL, *input = NULL;
    int iszero = -1;
    enum mmap_e type;
    uint64_t bplen, ncores;

    if (!PyArg_ParseTuple(args, "sslll", &dir, &input, &bplen, &type, &ncores)) {
        PyErr_SetString(PyExc_RuntimeError, "error parsing arguments");
        return NULL;
    }

    iszero = obf_evaluate(type, dir, input, bplen, ncores, g_verbose);
    if (iszero == -1) {
        PyErr_SetString(PyExc_RuntimeError, "zero test failed");
        return NULL;
    } else {
        return Py_BuildValue("i", iszero ? 0 : 1);
    }
}

static PyObject *
obf_wait_wrapper(PyObject *self, PyObject *args)
{
    PyObject *py_state;
    obf_state_t *s;

    if (!PyArg_ParseTuple(args, "O", &py_state))
        return NULL;

    s = (obf_state_t *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    obf_wait(s);

    Py_RETURN_NONE;
}

static PyObject *
obf_verbose(PyObject *self, PyObject *args)
{
    PyObject *py_verbose;

    if (!PyArg_ParseTuple(args, "O", &py_verbose))
        return NULL;

    g_verbose = PyObject_IsTrue(py_verbose);

    Py_RETURN_NONE;
}

static PyMethodDef
ObfMethods[] = {
    {"verbose", obf_verbose, METH_VARARGS,
     "Set verbosity."},
    {"init", obf_init_wrapper, METH_VARARGS,
     "Set up obfuscator."},
    {"encode_layer", obf_encode_layer_wrapper, METH_VARARGS,
     "Encode a branching program layer in each slot."},
    {"max_mem_usage", obf_max_mem_usage, METH_VARARGS,
     "Print out the maximum memory usage."},
    {"evaluate", obf_evaluate_wrapper, METH_VARARGS,
     "Evaluate the obfuscation."},
    {"wait", obf_wait_wrapper, METH_VARARGS,
     "Wait for threadpool to empty."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
