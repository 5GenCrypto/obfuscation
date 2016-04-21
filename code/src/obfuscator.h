#ifndef OBFUSCATOR_H
#define OBFUSCATOR_H

#include "thpool.h"
#include <mmap/mmap.h>
#include <aesrand.h>
#include "utils.h"

struct state {
    threadpool thpool;
    unsigned long secparam;
    enum mmap_e type;
    mmap_sk mmap;
    const mmap_vtable *vtable;
    aes_randstate_t rand;
    const char *dir;
    long nzs;
    /* fmpz_mat_t *randomizer; */
    fmpz_t field;
};

int
obf_init(struct state *s, enum mmap_e type, const char *dir,
         unsigned long secparam, unsigned long kappa, unsigned long nzs,
         unsigned long nthreads, unsigned long ncores);

void
obf_clear(struct state *s);

void
obf_encode_layer(struct state *s, long idx, long inp, long nrows, long ncols,
                 fmpz_mat_t zero, fmpz_mat_t one);

int
obf_evaluate(enum mmap_e type, char *dir, char *input, unsigned long bplen,
             unsigned long ncores);

#endif

