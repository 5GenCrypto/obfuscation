#include <Python.h>
#include "obfuscator.h"
#include "pyutils.h"

#include <mife/mife_internals.h>
#include <mmap/mmap_clt.h>
#include <mmap/mmap_gghlite.h>
#include <omp.h>

static void
obf_clear_wrapper(PyObject *self)
{
    obf_clear((struct state *) PyCapsule_GetPointer(self, NULL));
}

//
//
// Python functions
//
//


static PyObject *
obf_init_wrapper(PyObject *self, PyObject *args)
{
    long secparam, kappa, nzs, nthreads, ncores;
    enum mmap_e type;
    char *dir;
    struct state *s = NULL;
    PyObject *py_primes, *py_state;

    s = (struct state *) calloc(1, sizeof(struct state));
    if (s == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "memory allocation failed");
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "sllllll", &dir, &type, &secparam, &kappa,
                          &nzs, &nthreads, &ncores)) {
        PyErr_SetString(PyExc_RuntimeError, "unable to parse input");
        free(s);
        return NULL;
    }
    if (secparam <= 0 || kappa <= 0 || nzs <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "invalid input");
        free(s);
        return NULL;
    }

    (void) obf_init(s, type, dir, secparam, kappa, nzs, nthreads, ncores);

    py_primes = PyList_New(1);
    fprintf(stderr, "FIELD: ");
    fmpz_fprint(stderr, s->field);
    fprintf(stderr, "\n");
    PyList_SetItem(py_primes, 0, fmpz_to_py(s->field));

    py_state = PyCapsule_New((void *) s, NULL, obf_clear_wrapper);
    return PyTuple_Pack(2, py_state, py_primes);
}

static PyObject *
obf_encode_layer_wrapper(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_zero_ms, *py_one_ms;
    long inp, idx, nrows, ncols;
    fmpz_mat_t zero, one;
    fmpz_t rand;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OllllOO", &py_state, &idx, &nrows, &ncols,
                          &inp, &py_zero_ms, &py_one_ms))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    fmpz_init(rand);
    fmpz_randm_aes(rand, s->rand, s->field);
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

    // if (s->randomizer == NULL) {
    //     s->randomizer = (fmpz_mat_t *) malloc(sizeof(fmpz_mat_t));
    //     fmpz_mat_init(*s->randomizer, ncols, ncols);
    //     // fmpz_mat_one(*s->randomizer);
    //     for (int i = 0; i < ncols; i++)
    //         for(int j = 0; j < ncols; j++)
    //             fmpz_randm_aes(fmpz_mat_entry(*s->randomizer, i, j), s->rand, s->field);

    //     fmpz_mat_fprint_pretty(stderr, *s->randomizer);
    //     fprintf(stderr, "\n");

    //     fmpz_mat_mul(zero, zero, *s->randomizer);
    //     fmpz_mat_scalar_mod_fmpz(zero, zero, s->field);
    //     fmpz_mat_mul(one, one, *s->randomizer);
    //     fmpz_mat_scalar_mod_fmpz(one, one, s->field);

    // } else {
    //     fmpz_modp_matrix_inverse(*s->randomizer, *s->randomizer, nrows, s->field);
    //     fmpz_mat_fprint_pretty(stderr, *s->randomizer);
    //     fprintf(stderr, "\n");

    //     fmpz_mat_mul(zero, *s->randomizer, zero);
    //     fmpz_mat_scalar_mod_fmpz(zero, zero, s->field);
    //     fmpz_mat_mul(one, *s->randomizer, one);
    //     fmpz_mat_scalar_mod_fmpz(one, one, s->field);
    // }

    obf_encode_layer(s, idx, inp, nrows, ncols, zero, one);

    fmpz_clear(rand);
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
obf_wait(PyObject *self, PyObject *args)
{
    PyObject *py_state;
    struct state *s;

    if (!PyArg_ParseTuple(args, "O", &py_state))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    thpool_wait(s->thpool);

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
    {"wait", obf_wait, METH_VARARGS,
     "Wait for threadpool to empty."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
