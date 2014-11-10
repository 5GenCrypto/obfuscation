#ifndef __IND_OBFUSCATION__CIRCUIT_H__
#define __IND_OBFUSCATION__CIRCUIT_H__

#include <gmp.h>

struct circuit;

struct circuit *
circ_parse(const char *circname);

int
circ_copy_circuit(const char * const circname, const char * const newcirc);

int
circ_evaluate(const struct circuit *circ, const mpz_t *alphas,
              const mpz_t *betas, mpz_t out, const mpz_t q);

int
circ_evaluate_encoding(const struct circuit *circ, const mpz_t *xs,
                       const mpz_t *xones, const mpz_t *ys, const mpz_t *yones,
                       mpz_t out, const mpz_t q);

void
circ_cleanup(struct circuit *circ);

#endif
