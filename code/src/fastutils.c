#include <Python.h>
#include <gmp.h>
#include <omp.h>

static gmp_randstate_t rng;

static PyObject *
fastutils_genprimes(PyObject *self, PyObject *args)
{
    const long int num;
    const long int bitlength;
    long int i;
    PyObject *py_primes;

    if (!PyArg_ParseTuple(args, "ll", &num, &bitlength))
        return NULL;

    py_primes = PyList_New(num);
    if (py_primes == NULL)
        return NULL;

#pragma omp parallel for private(i)
    for (i = 0; i < num; ++i) {
        char *buffer;
        PyObject *ps_result, *pi_result;
        mpz_t p_tmp, p_unif;

        mpz_init(p_tmp);
        mpz_init(p_unif);

        mpz_urandomb(p_unif, rng, bitlength);
        mpz_nextprime(p_tmp, p_unif);
        buffer = mpz_get_str(NULL, 10, p_tmp);
        ps_result = PyString_FromString(buffer);

        free(buffer);
        mpz_clear(p_tmp);
        mpz_clear(p_unif);

        pi_result = PyNumber_Long(ps_result);
#pragma omp critical
        {
            PyList_SET_ITEM(py_primes, i, pi_result);
        }
    }

    return py_primes;
}

static PyMethodDef
FastutilsMethods[] = {
    {"genprimes", fastutils_genprimes, METH_VARARGS,
     "Generate random primes."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initfastutils(void)
{
    (void) Py_InitModule("fastutils", FastutilsMethods);

    gmp_randinit_default(rng);
}
