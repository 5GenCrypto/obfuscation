#include "thpool_fns.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>

void *
thpool_encode_elem(void *vargs)
{
    struct mlm_encode_elem_s *args =
        (struct mlm_encode_elem_s *) vargs;

    clt_encode(*args->out, args->mlm, args->nins, args->ins, args->pows);
    for (unsigned long j = 0; j < args->nins; ++j) {
        mpz_clear(args->ins[j]);
    }
    free(args->pows);
    free(args->ins);
    free(args);

    return NULL;
}

static int
write_vector(const char *dir, mpz_t *vector, long size, char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, name);
    (void) save_mpz_vector(fname, vector, size);
    free(fname);
    return 0;
}

void *
thpool_write_vector(void *vargs)
{
    double end;
    struct write_vector_s *args = (struct write_vector_s *) vargs;

    (void) write_vector(args->dir, args->vector, args->length, args->name);
    for (unsigned long i = 0; i < args->length; ++i) {
        mpz_clear(args->vector[i]);
    }

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       args->length, end - args->start);

    free(args->dir);
    free(args->name);
    free(args->vector);
    free(args);

    return NULL;
}

static int
write_layer(const char *dir, int inp, long idx, mpz_t *zero, mpz_t *one,
            long nrows, long ncols)
{
    mpz_t tmp;
    char *fname;
    int fnamelen;

    if (idx < 0)
        return 1;
    fnamelen = strlen(dir) + sizeof idx + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    mpz_init_set_ui(tmp, inp);
    (void) snprintf(fname, fnamelen, "%s/%ld.input", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    (void) snprintf(fname, fnamelen, "%s/%ld.zero", dir, idx);
    (void) save_mpz_vector(fname, zero, nrows * ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.one", dir, idx);
    (void) save_mpz_vector(fname, one, nrows * ncols);
    mpz_set_ui(tmp, nrows);
    (void) snprintf(fname, fnamelen, "%s/%ld.nrows", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    mpz_set_ui(tmp, ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.ncols", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    free(fname);
    mpz_clear(tmp);

    return 0;
}

void *
thpool_write_layer(void *vargs)
{
    double end;
    struct write_layer_s *args = (struct write_layer_s *) vargs;

    (void) write_layer(args->dir, args->inp, args->idx, args->zero, args->one,
                       args->nrows, args->ncols);
    for (int i = 0; i < args->nrows * args->ncols; ++i) {
        mpz_clears(args->zero[i], args->one[i], NULL);
    }
    free(args->zero);
    free(args->one);

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       2 * args->nrows * args->ncols, end - args->start);


    return NULL;
}

static int
write_element(const char *dir, mpz_t elem, const char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, name);
    (void) save_mpz_scalar(fname, elem);
    free(fname);
    return 0;
}

void *
thpool_write_element(void *vargs)
{
    struct write_element_s *args = (struct write_element_s *) vargs;

    (void) write_element(args->dir, *args->elem, args->name);

    return NULL;
}
