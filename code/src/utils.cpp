#include "utils.h"

#include <omp.h>
#include <stdlib.h>
#include <sys/time.h>

int g_verbose;

double
current_time(void)
{
    struct timeval t;
    (void) gettimeofday(&t, NULL);
    return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0));
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
mat_mult(mpz_t *a, const mpz_t *b, int size)
{
    mpz_t *tmparray;
    double start, end;

    start = current_time();
    tmparray = (mpz_t *) malloc(sizeof(mpz_t) * size * size);
    for (int i = 0; i < size * size; ++i) {
        mpz_init(tmparray[i]);
    }
#pragma omp parallel for
    for (int ctr = 0; ctr < size * size; ++ctr) {
        mpz_t tmp, sum;
        mpz_inits(tmp, sum, NULL);
        for (int i = 0; i < size; ++i) {
            mpz_mul(tmp,
                    a[i * size + ctr % size],
                    b[i + size * (ctr / size)]);
            mpz_add(sum, sum, tmp);
        }
        mpz_set(tmparray[ctr], sum);
        mpz_clears(tmp, sum, NULL);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "    Multiplying took: %f\n", end - start);

    start = current_time();
    for (int i = 0; i < size * size; ++i) {
        mpz_swap(a[i], tmparray[i]);
        mpz_clear(tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "    Swapping took: %f\n", end - start);
}

void
mat_mult_by_vects(mpz_t out, const mpz_t *s, const mpz_t *m, const mpz_t *t,
                  int size)
{
    double start, end;

    mpz_set_ui(out, 0);

    start = current_time();
#pragma omp parallel for
    for (int col = 0; col < size; ++col) {
        mpz_t tmp;
        mpz_t sum;
        mpz_inits(tmp, sum, NULL);
        for (int row = 0; row < size; ++row) {
            int elem = col * size + row;
            mpz_mul(tmp, s[row], m[elem]);
            mpz_add(sum, sum, tmp);
        }
        mpz_mul(tmp, sum, t[col]);
#pragma omp critical
        {
            mpz_add(out, out, tmp);
        }
        mpz_clears(tmp, sum, NULL);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Multiplying by vectors took: %f\n",
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

