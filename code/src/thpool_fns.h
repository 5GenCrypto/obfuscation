#ifndef THPOOL_FNS
#define THPOOL_FNS

#include <gmp.h>
#include <clt13.h>

struct mlm_encode_elem_s {
    clt_state *mlm;
    aes_randstate_t *rand;
    mpz_t *out;
    size_t nins;
    mpz_t *ins;
    int *pows;
};

void *
thpool_encode_elem(void *vargs);

struct write_vector_s {
    char *dir;
    mpz_t *vector;
    unsigned long length;
    char *name;
    double start;
};

void *
thpool_write_vector(void *vargs);

struct write_layer_s {
    char *dir;
    mpz_t *zero;
    mpz_t *one;
    long inp;
    long idx;
    long nrows;
    long ncols;
    double start;
};

void *
thpool_write_layer(void *vargs);

struct write_element_s {
    char *dir;
    mpz_t *elem;
    char *name;
};

void *
thpool_write_element(void *vargs);

#endif
