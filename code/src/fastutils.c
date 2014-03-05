#include <Python.h>
#include <gmp.h>

static gmp_randstate_t rng;

static PyObject *
fastutils_genprimes(PyObject *self, PyObject *args)
{
    const long int num;
    const long int bitlength;
    long int i;
    mpz_t p_unif, p_tmp;
    PyObject *py_primes;

    if (!PyArg_ParseTuple(args, "ll", &num, &bitlength))
        return NULL;

    mpz_init(p_tmp);
    mpz_init(p_unif);
    
    py_primes = PyList_New(num);

    for (i = 0; i < num; ++i) {
        char *buffer;
        PyObject *ps_result;
        PyObject *pi_result;

        mpz_urandomb(p_unif, rng, bitlength);
        mpz_nextprime(p_tmp, p_unif);
        buffer =  mpz_get_str(NULL, 10, p_tmp);
        ps_result = PyString_FromString(buffer);
        pi_result = PyNumber_Long(ps_result);
        PyList_SET_ITEM(py_primes, i, pi_result);

        free(buffer);
        /* Py_DECREF(ps_result); */
        /* Py_DECREF(pi_result); */
    }

    mpz_clear(p_tmp);
    mpz_clear(p_unif);

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
