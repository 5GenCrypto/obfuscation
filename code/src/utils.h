#ifndef __OBFUSCATION__UTILS_H__
#define __OBFUSCATION__UTILS_H__

#include <gmp.h>

#define OBFUSCATOR_OK 0
#define OBFUSCATOR_ERR (-1)

extern int g_verbose;

double
current_time(void);

int
seed_rng(gmp_randstate_t *rng);

int
load_mpz_scalar(const char *fname, mpz_t x);

int
save_mpz_scalar(const char *fname, const mpz_t x);

int
load_mpz_vector(const char *fname, mpz_t *m, const int len);

int
save_mpz_vector(const char *fname, const mpz_t *m, const int len);

void
mult_mats(mpz_t *result, const mpz_t *left, const mpz_t *right, const mpz_t q,
		  long m, long n, long p);

void
mult_vect_by_mat(mpz_t *v, const mpz_t *m, mpz_t q, int size, mpz_t *tmparray);

void
mult_vect_by_vect(mpz_t out, const mpz_t *m, const mpz_t *v, mpz_t q, int size);

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits);

#endif
