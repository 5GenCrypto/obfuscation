#include "zobfuscator.h"
#include "circuit.h"
#include "utils.h"
#include "pyutils.h"
#include "thpool.h"
#include "thpool_fns.h"

#include <aesrand.h>
#include <mmap/mmap.h>
#include <mmap/mmap_clt.h>
#include <omp.h>

static void
zobf_state_destructor(PyObject *self)
{
    zobf_state_t *s;

    s = (zobf_state_t *) PyCapsule_GetPointer(self, NULL);
    zobf_clear(s);
}

static PyObject *
zobf_init_wrapper(PyObject *self, PyObject *args)
{
    long secparam, kappa, nzs, nthreads, ncores, flags;
    char *dir;
    PyObject *py_pows;
    int *pows;
    zobf_state_t *s;

    if (!PyArg_ParseTuple(args, "lllOslll", &secparam, &kappa, &nzs, &py_pows,
                          &dir, &nthreads, &ncores, &flags)) {
        return NULL;
    }

    pows = (int *) malloc(sizeof(int) * nzs);
    for (long i = 0; i < nzs; ++i) {
        pows[i] = (int) PyLong_AsLong(PyList_GET_ITEM(py_pows, i));
    }

    s = zobf_init(dir, secparam, kappa, nzs, pows, nthreads, ncores, flags);

    return PyCapsule_New((void *) s, NULL, zobf_state_destructor);
}

static PyObject *
zobf_encode_circuit_wrapper(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_ys, *py_xdegs;
    mpz_t *xs, *ys;
    int *xdegs;
    int n, m, ydeg;
    char *circuit;
    zobf_state_t *s;

    if (!PyArg_ParseTuple(args, "OsOOiii", &py_state, &circuit, &py_ys,
                          &py_xdegs, &ydeg, &n, &m))
        return NULL;
    s = (zobf_state_t *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    xs = (mpz_t *) calloc(n, sizeof(mpz_t));
    for (int i = 0; i < n; ++i) {
        mpz_init_set_ui(xs[i], 1);
    }
    ys = (mpz_t *) calloc(m, sizeof(mpz_t));
    for (int i = 0; i < m; ++i) {
        mpz_init(ys[i]);
        py_to_mpz(ys[i], PyList_GET_ITEM(py_ys, i));
    }

    xdegs = (int *) calloc(n, sizeof(int));
    for (int i = 0; i < n; ++i) {
        xdegs[i] = PyLong_AsLong(PyList_GET_ITEM(py_xdegs, i));
    }

    zobf_encode_circuit(s, circuit, xs, ys, xdegs, ydeg, n, m);

    free(xdegs);
    for (int i = 0; i < n; ++i) {
        mpz_clear(xs[i]);
    }
    free(xs);
    for (int i = 0; i < m; ++i) {
        mpz_clear(ys[i]);
    }
    free(ys);
    
    Py_RETURN_NONE;
}

static PyObject *
zobf_evaluate_wrapper(PyObject *self, PyObject *args)
{
    char *circuit, *dir, *input;
    long n, m, nthreads, flags;
    int iszero;

    if (!PyArg_ParseTuple(args, "sssllll", &dir, &circuit, &input, &n, &m,
                          &nthreads, &flags))
        return NULL;

    iszero = zobf_evaluate(dir, circuit, input, n, m, nthreads, flags);
    return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyMethodDef
ObfMethods[] = {
    {"init", zobf_init_wrapper, METH_VARARGS,
     "Set up obfuscator."},
    {"encode_circuit", zobf_encode_circuit_wrapper, METH_VARARGS,
     "Encode circuit."},
    {"evaluate", zobf_evaluate_wrapper, METH_VARARGS,
     "Evaluate circuit."},
    {"max_mem_usage", obf_max_mem_usage, METH_VARARGS,
     "Compute the maximum memory usage."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_zobfuscator(void)
{
    (void) Py_InitModule("_zobfuscator", ObfMethods);
}
