#ifndef __IND_OBFUSCATION__UTILS_H__
#define __IND_OBFUSCATION__UTILS_H__

#include <python2.7/Python.h>
#include <gmp.h>

#define SUCCESS 1
#define FAILURE 0

extern int g_verbose;

struct state {
    gmp_randstate_t rng;
    unsigned long secparam;
    long size;
    unsigned long n;
    unsigned long nzs;
    unsigned long rho;
    mpz_t q;
    mpz_t pzt;
    mpz_t *gs;
    mpz_t *crt_coeffs;
    mpz_t *zinvs;
    char *dir;
};

void
state_destructor(PyObject *self);

void *
pymalloc(const size_t size);

PyObject *
mpz_to_py(const mpz_t in);

void
py_to_mpz(mpz_t out, PyObject *in);

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
mpz_mod_near(mpz_t out, const mpz_t a, const mpz_t b);

void
mat_mult(mpz_t *a, const mpz_t *b, int size);

void
mat_mult_by_vects(mpz_t out, const mpz_t *s, const mpz_t *m, const mpz_t *t,
				  int size);

int
is_zero(mpz_t c, mpz_t pzt, mpz_t x0, long nu);

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits);

int
write_scalar(struct state *s, mpz_t val, char *name);

int
write_setup_params(struct state *s, long nu);

PyObject *
obf_verbose(PyObject *self, PyObject *args);

PyObject *
obf_setup(PyObject *self, PyObject *args);

PyObject *
obf_max_mem_usage(PyObject *self, PyObject *args);

#endif
