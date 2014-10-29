#include "utils.h"

#include <fcntl.h>
#include <gmp.h>
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
    int ret = SUCCESS;
    FILE *f;

    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return FAILURE;
    }
    (void) mpz_inp_raw(x, f);

    (void) fclose(f);
    return ret;
}

int
save_mpz_scalar(const char *fname, const mpz_t x)
{
    FILE *f;

    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return FAILURE;
    }
    if (mpz_out_raw(f, x) == 0) {
        (void) fprintf(stderr, "ERROR: saving value failed!\n");
        (void) fclose(f);
        return FAILURE;
    }
    (void) fclose(f);
    return SUCCESS;
}

int
load_mpz_vector(const char *fname, mpz_t *m, const int len)
{
    int ret = SUCCESS;
    FILE *f;
    
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return FAILURE;
    }

    for (int i = 0; i < len; ++i) {
        (void) mpz_inp_raw(m[i], f);
    }

    (void) fclose(f);
    return ret;
}

int
save_mpz_vector(const char *fname, const mpz_t *m, const int len)
{
    int ret = SUCCESS;
    FILE *f;

    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return FAILURE;
    }

    for (int i = 0; i < len; ++i) {
        if (mpz_out_raw(f, m[i]) == 0) {
            fprintf(stderr, "ERROR: saving value failed!\n");
            (void) fclose(f);
            return FAILURE;
        }
    }

    (void) fclose(f);
    return ret;
}

void
mpz_mod_near(mpz_t out, const mpz_t a, const mpz_t b)
{
    mpz_t res, shift;

    mpz_inits(res, shift, NULL);
    mpz_mod(res, a, b);
    mpz_tdiv_q_2exp(shift, b, 1);
    if (mpz_cmp(res, shift) > 0) {
        mpz_sub(res, res, b);
    }
    mpz_set(out, res);
    mpz_clears(res, shift, NULL);
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
						right[k + p * ((i * m + j) / m)]);
				mpz_add(sum, sum, tmp);
				mpz_mod_near(sum, sum, q);
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
mult_vect_by_mat(mpz_t *v, const mpz_t *m, mpz_t q, int size)
{
	mpz_t *tmparray;
	double start, end;

	start = current_time();
	tmparray = (mpz_t *) malloc(sizeof(mpz_t) * size);
	for (int i = 0; i < size; ++i) {
        mpz_init(tmparray[i]);
    }
#pragma omp parallel for
	for (int i = 0; i < size; ++i) {
		mpz_t tmp, sum;
		mpz_inits(tmp, sum, NULL);
		for (int j = 0; j < size; ++j) {
			mpz_mul(tmp, v[j], m[i * size + j]);
			mpz_add(sum, sum, tmp);
		}
		mpz_mod_near(sum, sum, q);
		mpz_set(tmparray[i], sum);
		mpz_clears(tmp, sum, NULL);
	}
	for (int i = 0; i < size; ++i) {
        mpz_swap(v[i], tmparray[i]);
        mpz_clear(tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "    Multiplying took: %f\n", end - start);
}

void
mult_vect_by_vect(mpz_t out, const mpz_t *v, const mpz_t *u, mpz_t q, int size)
{
	double start, end;

	start = current_time();
	mpz_set_ui(out, 0);
#pragma omp parallel for
	for (int i = 0; i < size; ++i) {
		mpz_t tmp;
		mpz_init(tmp);
		mpz_mul(tmp, v[i], u[i]);
		mpz_mod_near(tmp, tmp, q);
#pragma omp critical
		{
			mpz_add(out, out, tmp);
		}
		mpz_clears(tmp, NULL);
	}
	end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Multiplying took: %f\n",
                       end - start);
}

int
is_zero(mpz_t c, mpz_t pzt, mpz_t q, long nu)
{
    mpz_t tmp;
    int ret;

    mpz_init(tmp);
    mpz_mul(tmp, c, pzt);
    mpz_mod_near(tmp, tmp, q);
    ret = (mpz_sizeinbase(tmp, 2) < (mpz_sizeinbase(q, 2) - nu)) ? 1 : 0;
    mpz_clear(tmp);

    return ret;
}

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits)
{
    mpz_t one;
    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, *rng, nbits);
    mpz_clear(one);
}

int
write_scalar(struct state *s, mpz_t val, char *name)
{
    char *fname;
    int fnamelen;
    double start, end;

    start = current_time();
    fnamelen = strlen(s->dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return FAILURE;
    (void) snprintf(fname, fnamelen, "%s/%s", s->dir, name);
    (void) save_mpz_scalar(fname, val);
    free(fname);
    end = current_time();

    if (g_verbose)
        (void) fprintf(stderr, "  Saving to file: %f\n", end - start);
    return SUCCESS;
}

int
write_setup_params(struct state *s, long nu, long size)
{
    char *fname;
    int len;
    mpz_t tmp;
    double start, end;
    start = current_time();

    len = strlen(s->dir) + 10;

    fname = (char *) malloc(sizeof(char) * len);
    if (fname == NULL)
        return FAILURE;

    mpz_init(tmp);

    // save size
	if (size > 0) {
		mpz_set_ui(tmp, size);
		(void) snprintf(fname, len, "%s/size", s->dir);
		(void) save_mpz_scalar(fname, tmp);
	}
    // save nu
    mpz_set_ui(tmp, nu);
    (void) snprintf(fname, len, "%s/nu", s->dir);
    (void) save_mpz_scalar(fname, tmp);
    // save q
    (void) snprintf(fname, len, "%s/q", s->dir);
    (void) save_mpz_scalar(fname, s->q);
    // save pzt
    (void) snprintf(fname, len, "%s/pzt", s->dir);
    (void) save_mpz_scalar(fname, s->pzt);

    mpz_clear(tmp);

    free(fname);

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Saving to file: %f\n", end - start);

    return SUCCESS;
}

