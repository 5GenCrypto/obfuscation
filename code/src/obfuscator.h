#ifndef OBFUSCATOR_H
#define OBFUSCATOR_H

#include "thpool.h"
#include <mmap/mmap.h>
#include <aesrand.h>
#include "utils.h"

typedef struct obf_state_s obf_state_t;

obf_state_t *
obf_init(enum mmap_e type, const char *dir, unsigned long secparam,
         unsigned long kappa, unsigned long nzs, unsigned long nthreads,
         unsigned long ncores);

void
obf_clear(obf_state_t *s);

void
obf_get_field(obf_state_t *s, fmpz_t *field);

void
obf_randomize_layer(obf_state_t *s, long nrows, long ncols, fmpz_mat_t zero,
                    fmpz_mat_t one);

void
obf_encode_layer(obf_state_t *s, long idx, long inp, long nrows, long ncols,
                 fmpz_mat_t zero, fmpz_mat_t one);

int
obf_evaluate(enum mmap_e type, char *dir, char *input, unsigned long bplen,
             unsigned long ncores);

void
obf_wait(obf_state_t *s);

#endif

