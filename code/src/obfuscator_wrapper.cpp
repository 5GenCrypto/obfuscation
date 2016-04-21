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
    long secparam, kappa, nzs, nthreads, ncores;
    enum mmap_e type;
    char *dir;
    obf_state_t *s;
    fmpz_t field;
    PyObject *py_primes, *py_state;

    if (!PyArg_ParseTuple(args, "sllllll", &dir, &type, &secparam, &kappa,
                          &nzs, &nthreads, &ncores)) {
        PyErr_SetString(PyExc_RuntimeError, "unable to parse input");
        return NULL;
    }
    if (secparam <= 0 || kappa <= 0 || nzs <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "invalid input");
        return NULL;
    }

    s = obf_init(type, dir, secparam, kappa, nzs, nthreads, ncores);
    if (s == NULL)
        return NULL;

    fmpz_init(field);
    obf_get_field(s, &field);
    py_primes = PyList_New(1);
    fprintf(stderr, "FIELD: ");
    fmpz_fprint(stderr, field);
    fprintf(stderr, "\n");
    PyList_SetItem(py_primes, 0, fmpz_to_py(field));
    fmpz_clear(field);

    py_state = PyCapsule_New((void *) s, NULL, obf_clear_wrapper);
    return PyTuple_Pack(2, py_state, py_primes);
}

static PyObject *
obf_encode_layer_wrapper(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_zero_ms, *py_one_ms;
    long inp, idx, nrows, ncols;
    fmpz_mat_t zero, one;
    obf_state_t *s;

    if (!PyArg_ParseTuple(args, "OllllOO", &py_state, &idx, &nrows, &ncols,
                          &inp, &py_zero_ms, &py_one_ms))
        return NULL;

    s = (obf_state_t *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    fmpz_mat_init(zero, nrows, ncols);
    fmpz_mat_init(one, nrows, ncols);

    for (long i = 0; i < nrows; ++i) {
        for (long j = 0; j < ncols; ++j) {
            py_to_fmpz(fmpz_mat_entry(zero, i, j),
                       PyList_GET_ITEM(PyList_GET_ITEM(py_zero_ms, 0), i * nrows + j));
            py_to_fmpz(fmpz_mat_entry(one, i, j),
                       PyList_GET_ITEM(PyList_GET_ITEM(py_one_ms, 0), i * nrows + j));
        }
    }

    obf_encode_layer(s, idx, inp, nrows, ncols, zero, one);

    fmpz_mat_clear(zero);
    fmpz_mat_clear(one);

    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate_wrapper(PyObject *self, PyObject *args)
{
    char *dir, *input;
    int iszero = -1;
    enum mmap_e type;
    long bplen, ncores;

    if (!PyArg_ParseTuple(args, "sslll", &dir, &input, &bplen, &type, &ncores)) {
        PyErr_SetString(PyExc_RuntimeError, "error parsing arguments");
        return NULL;
    }

    iszero = obf_evaluate(type, dir, input, bplen, ncores);
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
