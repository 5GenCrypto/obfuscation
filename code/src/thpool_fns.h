#ifndef THPOOL_FNS
#define THPOOL_FNS

#include <gmp.h>

#include "clt_mlm.h"

struct mlm_encode_vector_elem_state {
    struct clt_mlm_state *mlm;
    int *indices;
    int *pows;
};

struct mlm_encode_vector_elem_s {
    mpz_t *out;
    mpz_t *elems;
    struct mlm_encode_vector_elem_state *s;
};

void *
thpool_encode_vector_elem(void *vargs);

struct write_vector_s {
    char *dir;
    mpz_t *vector;
    unsigned long length;
    char *name;
    double start;
    struct mlm_encode_vector_elem_state *state;
};

void *
thpool_write_vector(void *vargs);


#endif
