#include "utils.h"

#include <fcntl.h>
#include <gmp.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

int g_verbose;

// XXX: The use of /dev/urandom is not secure; however, the supercomputer we run
// on doesn't appear to have enough entropy, and blocks for long periods of
// time.  Thus, we use /dev/urandom instead.
#ifndef RANDFILE
#  define RANDFILE "/dev/urandom"
#endif

double
current_time(void)
{
    struct timeval t;
    (void) gettimeofday(&t, NULL);
    return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0));
}

int
seed_rng(gmp_randstate_t *rng)
{
    int file;
    if ((file = open(RANDFILE, O_RDONLY)) == -1) {
        (void) fprintf(stderr, "Error opening %s\n", RANDFILE);
        return 1;
    } else {
        unsigned long seed;
        if (read(file, &seed, sizeof seed) == -1) {
            (void) fprintf(stderr, "Error reading from %s\n", RANDFILE);
            (void) close(file);
            return 1;
        } else {
            if (g_verbose)
                (void) fprintf(stderr, "  Seed: %lu\n", seed);

            gmp_randinit_default(*rng);
            gmp_randseed_ui(*rng, seed);
        }
    }
    if (file != -1)
        (void) close(file);
    return 0;
}

int
fmpz_mod_poly_fread_raw(FILE * f, fmpz_mod_poly_t poly)
{
    slong i, length;
    fmpz_t coeff;
    ulong res;

    fmpz_init(coeff);
    if (flint_fscanf(f, "%wd ", &length) != 1) {
        fmpz_clear(coeff);
        return 0;
    }

    fmpz_inp_raw(coeff,f);
    fmpz_mod_poly_init(poly, coeff);
    fmpz_mod_poly_fit_length(poly, length);

    poly->length = length;
    flint_fscanf(f, " ");

    for (i = 0; i < length; i++)
    {
        flint_fscanf(f, " ");
        res = fmpz_inp_raw(coeff, f);

        fmpz_mod_poly_set_coeff_fmpz(poly,i,coeff);

        if (!res)
        {
            poly->length = i;
            fmpz_clear(coeff);
            return 0;
        }
    }

    fmpz_clear(coeff);
    _fmpz_mod_poly_normalise(poly);

    return 1;
}


static int
_fmpz_mod_poly_fprint_raw(FILE * file, const fmpz *poly, slong len, const fmpz_t p)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%wd ", len);
    if (r <= 0)
        return r;

    r = fmpz_out_raw(file, p);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    r = flint_fprintf(file, " ");
    if (r <= 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = flint_fprintf(file, " ");
        if (r <= 0)
            return r;
        r = fmpz_out_raw(file, poly + i);
        if (r <= 0)
            return r;
    }

    return r;
}

int
fmpz_mod_poly_fprint_raw(FILE * file, const fmpz_mod_poly_t poly)
{
    return _fmpz_mod_poly_fprint_raw(file, poly->coeffs, poly->length,
        &(poly->p));
}


int
load_mpz_scalar(const char *fname, mpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    (void) mpz_inp_raw(x, f);
    (void) fclose(f);
    return 0;
}

int
save_mpz_scalar(const char *fname, const mpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    if (mpz_out_raw(f, x) == 0) {
        (void) fclose(f);
        return 1;
    }
    (void) fclose(f);
    return 0;
}

int
save_fmpz_scalar(const char *fname, const fmpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    if (fmpz_out_raw(f, x) == 0) {
        (void) fclose(f);
        return 1;
    }
    (void) fclose(f);
    return 0;
}

int
load_mpz_vector(const char *fname, mpz_t *m, const int len)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        (void) mpz_inp_raw(m[i], f);
    }
    (void) fclose(f);
    return 0;
}

int
save_mpz_vector(const char *fname, const mpz_t *m, const int len)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        if (mpz_out_raw(f, m[i]) == 0) {
            (void) fclose(f);
            return 1;
        }
    }
    (void) fclose(f);
    return 0;
}

int
load_gghlite_enc_vector(const char *fname, gghlite_enc_t *m, const int len)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        (void) fmpz_mod_poly_fread_raw(f, m[i]);
    }
    (void) fclose(f);
    return 0;
    
}
int
save_gghlite_enc_vector(const char *fname, const gghlite_enc_t *m, const int len)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        if (fmpz_mod_poly_fprint_raw(f, m[i]) == 0) {
            (void) fclose(f);
            return 1;
        }
    }
    (void) fclose(f);
    return 0;
}

void
mult_mpz_matrices(mpz_t *result, const mpz_t *left, const mpz_t *right,
                  const mpz_t q, long m, long n, long p)
{
    mpz_t *tmparray;
    double start, end;

    start = current_time();
    tmparray = (mpz_t *) malloc(sizeof(mpz_t) * m * p);
    for (int i = 0; i < m * p; ++i) {
        mpz_init(tmparray[i]);
    }
#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            mpz_t tmp, sum;
            mpz_inits(tmp, sum, NULL);
            for (int k = 0; k < n; ++k) {
                mpz_mul(tmp,
                        left[k * m + (i * m + j) % m],
                        right[k + n * ((i * m + j) / m)]);
                mpz_add(sum, sum, tmp);
                mpz_mod(sum, sum, q);
            }
            mpz_set(tmparray[i * n + j], sum);
            mpz_clears(tmp, sum, NULL);
        }
    }
    for (int i = 0; i < m * p; ++i) {
        mpz_swap(result[i], tmparray[i]);
        mpz_clear(tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, " Multiplying took: %f\n", end - start);
}

void
mult_gghlite_enc_matrices(gghlite_enc_t *result, const gghlite_params_t pp,
                          const gghlite_enc_t *left, const gghlite_enc_t *right,
                          const fmpz_t q, long m, long n, long p)
{
    gghlite_enc_t *tmparray;
    double start, end;

    start = current_time();
    tmparray = (gghlite_enc_t *) malloc(sizeof(gghlite_enc_t) * m * p);
    for (int i = 0; i < m * p; ++i) {
        gghlite_enc_init(tmparray[i], pp);
    }
#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            gghlite_enc_t tmp, sum;
            gghlite_enc_init(tmp, pp);
            gghlite_enc_init(sum, pp);
            for (int k = 0; k < n; ++k) {
                gghlite_enc_mul(tmp, pp,
                                left[k * m + (i * m + j) % m],
                                right[k + n * ((i * m + j) / m)]);
                gghlite_enc_add(sum, pp, sum, tmp);
                // mpz_mod(sum, sum, q);
            }
            gghlite_enc_set(tmparray[i * n + j], sum);
            gghlite_enc_clear(tmp);
            gghlite_enc_clear(sum);
        }
    }
    for (int i = 0; i < m * p; ++i) {
        gghlite_enc_set(result[i], tmparray[i]);
        gghlite_enc_clear(tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, " Multiplying took: %f\n", end - start);
}


void
mult_vect_by_mat(mpz_t *v, const mpz_t *m, mpz_t q, int size, mpz_t *tmparray)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        mpz_t tmp, sum;
        mpz_inits(tmp, sum, NULL);
        for (int j = 0; j < size; ++j) {
            mpz_mul(tmp, v[j], m[i * size + j]);
            mpz_add(sum, sum, tmp);
            mpz_mod(sum, sum, q);
        }
        mpz_set(tmparray[i], sum);
        mpz_clears(tmp, sum, NULL);
    }
    for (int i = 0; i < size; ++i) {
        mpz_swap(v[i], tmparray[i]);
    }
}

void
mult_vect_by_vect(mpz_t out, const mpz_t *v, const mpz_t *u, mpz_t q, int size)
{
    mpz_set_ui(out, 0);
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mul(tmp, v[i], u[i]);
#pragma omp critical
        {
            mpz_add(out, out, tmp);
            mpz_mod(out, out, q);
        }
        mpz_clears(tmp, NULL);
    }
}

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits)
{
    mpz_t one;
    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, *rng, nbits);
    mpz_clear(one);
}
