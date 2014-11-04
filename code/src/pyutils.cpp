#include "pyutils.h"
#include "utils.h"
#include "mpz_pylong.h"

#include <sys/resource.h>

void
state_destructor(PyObject *self)
{
    struct state *s;

    s = (struct state *) PyCapsule_GetPointer(self, NULL);
    if (s) {
        if (s->gs) {
            for (unsigned long i = 0; i < s->n; ++i) {
                mpz_clear(s->gs[i]);
            }
            free(s->gs);
        }
        if (s->crt_coeffs) {
            for (unsigned long i = 0; i < s->n; ++i) {
                mpz_clear(s->crt_coeffs[i]);
            }
            free(s->crt_coeffs);
        }
        if (s->zinvs) {
            for (unsigned long i = 0; i < s->nzs; ++i) {
                mpz_clear(s->zinvs[i]);
            }
            free(s->zinvs);
        }
        gmp_randclear(s->rng);
        mpz_clears(s->q, s->pzt, NULL);
    }
}

void *
pymalloc(const size_t size)
{
    void * r;
    if ((r = malloc(size)) == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory");
    }
    return r;
}

PyObject *
mpz_to_py(const mpz_t in)
{
    PyObject *outs, *out;
    char *buffer;

    buffer = mpz_get_str(NULL, 10, in);
    outs = PyString_FromString(buffer);
    out = PyNumber_Long(outs);
    free(buffer);
    return out;
}

void
py_to_mpz(mpz_t out, PyObject *in)
{
    (void) mpz_set_pylong(out, in);
}

PyObject *
obf_verbose(PyObject *self, PyObject *args)
{
    PyObject *py_verbose;

    if (!PyArg_ParseTuple(args, "O", &py_verbose))
        return NULL;

    g_verbose = PyObject_IsTrue(py_verbose);
    if (g_verbose == -1)
        return NULL;

    Py_RETURN_NONE;
}

PyObject *
obf_max_mem_usage(PyObject *self, PyObject *args)
{
    struct rusage usage;

    (void) getrusage(RUSAGE_SELF, &usage);
    (void) fprintf(stderr, "Max memory usage: %ld\n", usage.ru_maxrss);

    Py_RETURN_NONE;
}
