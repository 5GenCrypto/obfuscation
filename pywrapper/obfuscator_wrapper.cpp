#include <Python.h>
#include "obfuscator.h"
#include "pyutils.h"

static void
obf_clear_wrapper(PyObject *self)
{
    obf_clear((obf_state_t *) PyCapsule_GetPointer(self, NULL));
}

static PyObject *
obf_init_wrapper(PyObject *self, PyObject *args)
{
    long type_ = 0, secparam = 0, kappa = 0, nzs = 0, nthreads = 0,
        ncores = 0, flags = 0;
    enum mmap_e type;
    char *dir = NULL;
    obf_state_t *s = NULL;

    if (!PyArg_ParseTuple(args, "slllllll", &dir, &type_, &secparam, &kappa,
                          &nzs, &nthreads, &ncores, &flags)) {
        PyErr_SetString(PyExc_RuntimeError, "unable to parse input");
        return NULL;
    }
    if (secparam <= 0 || kappa <= 0 || nzs <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "invalid input");
        return NULL;
    }

    switch (type_) {
    case 0:
        type = MMAP_CLT;
        break;
    case 1:
        type = MMAP_GGHLITE;
        break;
    case 2:
        type = MMAP_DUMMY;
        break;
    default:
        PyErr_SetString(PyExc_RuntimeError, "invalid mmap type");
        return NULL;
    }

    s = obf_init(type, dir, secparam, kappa, nzs, nthreads, ncores, flags);
    if (s == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "initialization failed");
        return NULL;
    }

    return PyCapsule_New((void *) s, NULL, obf_clear_wrapper);
}

static PyObject *
obf_encode_layer_wrapper(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_pows, *py_mats;
    long n, idx, nrows, ncols, inp, rflag;
    ssize_t length;
    int **pows;
    fmpz_mat_t *mats;
    obf_state_t *s;

    // TODO: can probably get nrows, ncols length from matrices
    if (!PyArg_ParseTuple(args, "OlOOlllll", &py_state, &n, &py_pows, &py_mats,
                          &idx, &nrows, &ncols, &inp, &rflag))
        return NULL;

    s = (obf_state_t *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "unable to extract obf state");
        return NULL;
    }

    length = PyList_Size(PyList_GetItem(py_pows, 0));
    // TODO: check pow lengths

    pows = (int **) calloc(n, sizeof(int *));
    for (long c = 0; c < n; ++c) {
        pows[c] = (int *) calloc(length, sizeof(int));
        for (ssize_t i = 0; i < length; ++i) {
            pows[c][i] = PyLong_AsLong(
                PyList_GetItem(
                    PyList_GetItem(py_pows, c), i));
        }
    }

    mats = (fmpz_mat_t *) calloc(n, sizeof(fmpz_mat_t));
    for (long c = 0; c < n; ++c) {
        fmpz_mat_init(mats[c], nrows, ncols);
        for (long i = 0; i < nrows; ++i) {
            for (long j = 0; j < ncols; ++j) {
                py_to_fmpz(fmpz_mat_entry(mats[c], i, j),
                           PyList_GetItem(
                               PyList_GetItem(
                                   PyList_GetItem(py_mats, c), i), j));
            }
        }
    }

    obf_encode_layer(s, n, pows, mats, idx, inp,
                     (encode_layer_randomization_flag_t) rflag);

    for (long c = 0; c < n; ++c) {
        fmpz_mat_clear(mats[c]);
    }
    free(mats);
    // TODO: make sure that pows gets cleared elsewhere

    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate_wrapper(PyObject *self, PyObject *args)
{
    PyObject *py_input;
    uint64_t *input;
    char *dir = NULL;
    long iszero, type_, flags;
    enum mmap_e type;
    uint64_t len = 0, bplen = 0, ncores = 0;

    if (!PyArg_ParseTuple(args, "zOllll", &dir, &py_input, &type_, &bplen,
                          &ncores, &flags)) {
        PyErr_SetString(PyExc_RuntimeError, "error parsing arguments");
        return NULL;
    }

    switch (type_) {
    case 0:
        type = MMAP_CLT;
        break;
    case 1:
        type = MMAP_GGHLITE;
        break;
    case 2:
        type = MMAP_DUMMY;
        break;
    default:
        PyErr_SetString(PyExc_RuntimeError, "invalid mmap type");
        return NULL;
    }

    len = PyList_Size(py_input);
    input = (uint64_t *) calloc(len, sizeof(uint64_t));
    for (int i = 0; i < len; ++i) {
        input[i] = PyLong_AsLong(PyList_GetItem(py_input, i));
    }

    iszero = obf_evaluate(type, dir, len, input, bplen, ncores, flags);
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

static PyMethodDef
ObfMethods[] = {
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


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef obfmodule = {
    PyModuleDef_HEAD_INIT,
    "_obfuscator",
    NULL,
    -1,
    ObfMethods
};

PyMODINIT_FUNC
PyInit__obfuscator(void)
{
    return PyModule_Create(&obfmodule);
}
#else
PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
#endif
