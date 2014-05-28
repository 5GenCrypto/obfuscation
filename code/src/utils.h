#ifndef __IND_OBFUSCATION__UTILS_H__
#define __IND_OBFUSCATION__UTILS_H__

#include <gmp.h>

#define SUCCESS 1
#define FAILURE 0

extern int g_verbose;

double
current_time(void);

int
load_mpz_scalar(const char *fname, mpz_t x);

int
save_mpz_scalar(const char *fname, const mpz_t x);

int
load_mpz_vector(const char *fname, mpz_t *m, const int len);

int
save_mpz_vector(const char *fname, const mpz_t *m, const int len);

void
mpz_mod_near(mpz_t out, const mpz_t a, const mpz_t b);

void
mat_mult(mpz_t *a, const mpz_t *b, int size);

void
mat_mult_by_vects(mpz_t out, const mpz_t *s, const mpz_t *m, const mpz_t *t,
                      int size);

int
is_zero(mpz_t c, mpz_t pzt, mpz_t x0, long nu);

#endif
