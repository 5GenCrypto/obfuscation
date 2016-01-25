#ifndef __CLT_MLM_H__
#define __CLT_MLM_H__

#include <gmp.h>

struct clt_mlm_state {
	gmp_randstate_t rng;
    unsigned long secparam;
    unsigned long n;
    unsigned long nzs;
    unsigned long rho;
    unsigned long nu;
    mpz_t q;
    mpz_t pzt;
    mpz_t *gs;
    mpz_t *crt_coeffs;
    mpz_t *zinvs;
};


int
clt_mlm_setup(struct clt_mlm_state *s, const char *dir, const long *pows,
              long kappa, long size, int verbose);

void
clt_mlm_cleanup(struct clt_mlm_state *s);

void
clt_mlm_encode(struct clt_mlm_state *s, mpz_t out, size_t nins,
			   const mpz_t *ins, unsigned int nzs, const int *indices,
			   const int *pows);

int
clt_mlm_is_zero(const mpz_t c, const mpz_t pzt, const mpz_t q, long nu);

#endif
