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

void
mult_mats(mpz_t *result, const mpz_t *left, const mpz_t *right, const mpz_t q,
          long m, long n, long p)
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
