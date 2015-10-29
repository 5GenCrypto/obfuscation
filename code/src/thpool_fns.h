#ifndef THPOOL_FNS
#define THPOOL_FNS

#include <Python.h>
#include <gmp.h>

#include "clt_mlm.h"

struct mlm_encode_vector_elem_state {
    struct clt_mlm_state *mlm;
    PyObject *py_vectors;
    int *indices;
    int *pows;
};

struct mlm_encode_vector_elem_args {
    int i;
    mpz_t *elem;
    struct mlm_encode_vector_elem_state *s;
};

void *thpool_encode_vector_elem(void *voidargs);

#endif
