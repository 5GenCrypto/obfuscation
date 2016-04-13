#ifndef THPOOL_FNS
#define THPOOL_FNS

#include "utils.h"

#include <aesrand.h>
#include <gmp.h>
#include <clt13.h>
#include <gghlite.h>

struct encode_elem_s {
    enum mmap_e mmap;
    void *mlm;
    aes_randstate_t *rand;
    void *out;
    size_t nins;
    void *ins;
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
    enum mmap_e mmap;
    char *dir;
    void *zero;
    void *one;
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
    clt_elem_t *elem;
    char *name;
};

void *
thpool_write_element(void *vargs);

#endif
