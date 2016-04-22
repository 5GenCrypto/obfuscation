#ifndef ZOBFUSCATOR_H
#define ZOBFUSCATOR_H

#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct zobf_state_s zobf_state_t;

zobf_state_t *
zobf_init(const char *dir, uint64_t secparam, uint64_t kappa, uint64_t nzs,
          int *pows, uint64_t nthreads, uint64_t ncores, bool verbose);

void
zobf_clear(zobf_state_t *s);

void
zobf_encode_circuit(zobf_state_t *s, const char *circuit, const mpz_t *ys,
                    const int *xdegs, int ydeg, int n, int m);

int
zobf_evaluate(const char *dir, const char *circuit, const char *input,
              int n, int m, uint64_t nthreads);

#ifdef __cplusplus
}
#endif

#endif
