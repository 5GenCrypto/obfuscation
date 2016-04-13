#ifndef __OBFUSCATION__PYUTILS_H__
#define __OBFUSCATION__PYUTILS_H__

#include <python2.7/Python.h>
#include <gmp.h>
#include <flint/fmpz.h>

PyObject *
mpz_to_py(const mpz_t in);

int
py_to_mpz(mpz_t out, PyObject *in);

PyObject *
fmpz_to_py(const fmpz_t in);

int
py_to_fmpz(fmpz_t out, PyObject *in);

PyObject *
obf_verbose(PyObject *self, PyObject *args);

PyObject *
obf_max_mem_usage(PyObject *self, PyObject *args);

#endif
