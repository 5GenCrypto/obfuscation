#include "pyutils.h"
#include "utils.h"
#include "mpz_pylong.h"

#include <sys/resource.h>

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

int
py_to_mpz(mpz_t out, PyObject *in)
{
    if (mpz_set_pylong(out, in) == -1) {
        fprintf(stderr, "error converting python object to mpz\n");
        return 1;
    }
    return 0;
}

PyObject *
fmpz_to_py(const fmpz_t in)
{
    PyObject *outs, *out;
    char *buffer;

    buffer = fmpz_get_str(NULL, 10, in);
    outs = PyString_FromString(buffer);
    out = PyNumber_Long(outs);
    free(buffer);
    return out;
}

int
py_to_fmpz(fmpz_t out, PyObject *in)
{
    mpz_t tmp;
    mpz_init(tmp);
    if (py_to_mpz(tmp, in) == 1)
        return 1;
    fmpz_set_mpz(out, tmp);
    mpz_clear(tmp);
    return 0;
}

PyObject *
obf_max_mem_usage(PyObject *self, PyObject *args)
{
    struct rusage usage;

    (void) getrusage(RUSAGE_SELF, &usage);
    (void) fprintf(stderr, "Max memory usage: %ld\n", usage.ru_maxrss);

    Py_RETURN_NONE;
}
