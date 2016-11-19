#ifndef THPOOL_FNS
#define THPOOL_FNS

#include "utils.h"

#include <aesrand.h>
#include <mmap/mmap.h>
#include <stdint.h>

struct encode_elem_s {
    const mmap_vtable *vtable;
    mmap_ro_sk sk;
    int n;
    fmpz_t *plaintext;
    int *group;
    mmap_enc *enc;
    aes_randstate_t *rand;
};

void *
thpool_encode_elem(void *vargs);

struct write_layer_s {
    const mmap_vtable *vtable;
    const char *dir;
    uint64_t n;
    mmap_enc_mat_t **enc_mats;
    char **names;
    long inp;
    long idx;
    long nrows;
    long ncols;
    double start;
    bool verbose;
};

void *
thpool_write_layer(void *vargs);

struct write_element_s {
    const char *dir;
    mmap_enc *elem;
    char *name;
};

#endif
