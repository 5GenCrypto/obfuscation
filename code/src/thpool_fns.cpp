#include "thpool_fns.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>

void *
thpool_encode_elem(void *vargs)
{
    struct encode_elem_s *args = (struct encode_elem_s *) vargs;
    aes_randstate_t rand;

    aes_randinit(rand);
    switch (args->mmap) {
    case MMAP_CLT:
        clt_encode(*((clt_elem_t *) args->out), (clt_state *) args->mlm,
                   args->nins, (clt_elem_t *) args->ins, args->pows, rand);

        for (unsigned long j = 0; j < args->nins; ++j) {
            mpz_clear(((clt_elem_t *) args->ins)[j]);
        }
        break;
    case MMAP_GGHLITE:
    {
        gghlite_clr_t e;
        gghlite_clr_init(e);
        fmpz_poly_set_coeff_fmpz(e, 0, ((fmpz_t *) args->ins)[0]);
        gghlite_enc_set_gghlite_clr(*((gghlite_enc_t *) args->out),
                                    *((gghlite_sk_t *) args->mlm),
                                    e, 1, args->pows, 0, rand);
        gghlite_clr_clear(e);
        fmpz_clear(((fmpz_t *) args->ins)[0]);
        break;
    }
    }
    aes_randclear(rand);

    if (args->pows)
        free(args->pows);
    if (args->ins)
        free(args->ins);
    free(args);

    return NULL;
}

static int
write_layer(const char *dir, enum mmap_e mmap, int inp, long idx, void *zero,
            void *one, long nrows, long ncols)
{
    mpz_t tmp;
    char *fname;
    int fnamelen;
    int ret = 1;

    if (idx < 0)
        return 1;

    mpz_init_set_ui(tmp, inp);
    fnamelen = strlen(dir) + sizeof idx + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        goto cleanup;

    (void) snprintf(fname, fnamelen, "%s/%ld.input", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    (void) snprintf(fname, fnamelen, "%s/%ld.zero", dir, idx);
    switch (mmap) {
    case MMAP_CLT:
        (void) clt_vector_save(fname, (clt_elem_t *) zero, nrows * ncols);
        break;
    case MMAP_GGHLITE:
        (void) save_gghlite_enc_vector(fname, (gghlite_enc_t *) zero,
                                       nrows * ncols);
        break;
    }
    (void) snprintf(fname, fnamelen, "%s/%ld.one", dir, idx);
    switch (mmap) {
    case MMAP_CLT:
        (void) clt_vector_save(fname, (clt_elem_t *) one, nrows * ncols);
        break;
    case MMAP_GGHLITE:
        (void) save_gghlite_enc_vector(fname, (gghlite_enc_t *) one,
                                       nrows * ncols);
        break;
    }
    mpz_set_ui(tmp, nrows);
    (void) snprintf(fname, fnamelen, "%s/%ld.nrows", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    mpz_set_ui(tmp, ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.ncols", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    ret = 0;

cleanup:
    if (fname)
        free(fname);
    mpz_clear(tmp);

    return ret;
}

void *
thpool_write_layer(void *vargs)
{
    double end;
    struct write_layer_s *args = (struct write_layer_s *) vargs;

    (void) write_layer(args->dir, args->mmap, args->inp, args->idx, args->zero,
                       args->one, args->nrows, args->ncols);
    switch (args->mmap) {
    case MMAP_CLT:
        for (int i = 0; i < args->nrows * args->ncols; ++i) {
            mpz_clears(((mpz_t *) args->zero)[i], ((mpz_t *) args->one)[i], NULL);
        }
        break;
    case MMAP_GGHLITE:
        for (int i = 0; i < args->nrows * args->ncols; ++i) {
            gghlite_enc_clear(((gghlite_enc_t *) args->zero)[i]);
            gghlite_enc_clear(((gghlite_enc_t *) args->one)[i]);
        }
        break;
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
write_element(const char *dir, clt_elem_t elem, const char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, name);
    (void) clt_elem_save(fname, elem);
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
