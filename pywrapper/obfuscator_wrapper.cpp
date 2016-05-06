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
    PyObject *py_state, *py_zero_pows, *py_one_pows, *py_zero_ms, *py_one_ms;
    long idx, nrows, ncols, inp, rflag;
    ssize_t length;
    int *zero_pows, *one_pows;
    fmpz_mat_t zero, one;
    obf_state_t *s;

    // TODO: can probably get nrows, ncols length from matrices
    if (!PyArg_ParseTuple(args, "OlllllOOOO", &py_state,
                          &idx, &nrows, &ncols, &inp, &rflag,
                          &py_zero_pows, &py_one_pows,
                          &py_zero_ms, &py_one_ms))
        return NULL;

    s = (obf_state_t *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "unable to extract obf state");
        return NULL;
    }

    length = PyList_Size(py_zero_pows);
    if (PyList_Size(py_one_pows) != length) {
        PyErr_SetString(PyExc_RuntimeError, "pow lengths unequal");
        return NULL;
    }

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

    obf_encode_layer(s, idx, inp, (encode_layer_randomization_flag_t) rflag,
                     zero_pows, one_pows, zero, one);

    fmpz_mat_clear(zero);
    fmpz_mat_clear(one);

    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate_wrapper(PyObject *self, PyObject *args)
{
    char *dir = NULL, *input = NULL;
    long iszero, type_, flags;
    enum mmap_e type;
    uint64_t bplen = 0, ncores = 0;

    if (!PyArg_ParseTuple(args, "zzllll", &dir, &input, &type_, &bplen,
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
    default:
        PyErr_SetString(PyExc_RuntimeError, "invalid mmap type");
        return NULL;
    }

    iszero = obf_evaluate(type, dir, input, bplen, ncores, flags);
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
