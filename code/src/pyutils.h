#ifndef __IND_OBFUSCATION__PYUTILS_H__
#define __IND_OBFUSCATION__PYUTILS_H__

#include <Python.h>
#include <gmp.h>

void *
pymalloc(const size_t size);

PyObject *
mpz_to_py(const mpz_t in);

void
py_to_mpz(mpz_t out, PyObject *in);

#endif
