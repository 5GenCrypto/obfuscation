#include "obfuscator.h"
#include "thpool_fns.h"

#include <mmap/mmap_clt.h>
#include <mmap/mmap_gghlite.h>
#include <omp.h>

static void
mult_matrices(const mmap_vtable *vtable, const mmap_pp *pp, mmap_enc *result,
              const mmap_enc *left, const mmap_enc *right, long m, long n,
              long p)
{
    mmap_enc *tmparray;
    double start, end;

    start = current_time();
    tmparray = (mmap_enc *) malloc(sizeof(mmap_enc) * m * p);
    for (int i = 0; i < m * p; ++i) {
        vtable->enc->init(&tmparray[i], pp);
    }
#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            mmap_enc tmp, sum;
            vtable->enc->init(&tmp, pp);
            vtable->enc->init(&sum, pp);
            for (int k = 0; k < n; ++k) {
                vtable->enc->mul(&tmp, pp,
                                 &left[k * m + (i * m + j) % m],
                                 &right[k + n * ((i * m + j) / m)]);
                vtable->enc->add(&sum, pp, &sum, &tmp);
            }
            vtable->enc->set(&tmparray[i * n + j], &sum);
            vtable->enc->clear(&tmp);
            vtable->enc->clear(&sum);
        }
    }
    for (int i = 0; i < m * p; ++i) {
        vtable->enc->set(&result[i], &tmparray[i]);
        vtable->enc->clear(&tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, " Multiplying took: %f\n", end - start);
}

int
obf_init(struct state *s, enum mmap_e type, const char *dir,
         unsigned long secparam, unsigned long kappa, unsigned long nzs,
         unsigned long nthreads, unsigned long ncores)
{
    s->secparam = secparam;
    s->type = type;
    s->dir = dir;
    s->nzs = nzs;

    switch (s->type) {
    case MMAP_CLT:
        s->vtable = &clt13_vtable;
        break;
    case MMAP_GGHLITE:
        s->vtable = &gghlite_vtable;
        break;
    default:
        return 1;
    }

    fmpz_init(s->field);
    (void) aes_randinit(s->rand);
    s->thpool = thpool_init(nthreads);
    (void) omp_set_num_threads(ncores);

    if (g_verbose) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
    }

    s->vtable->sk->init(&s->mmap, secparam, kappa, nzs, s->rand);
    {
        char fname[100];
        FILE *fp;

        snprintf(fname, 100, "%s/params", s->dir);
        fp = fopen(fname, "w+b");
        s->vtable->pp->fwrite(s->vtable->sk->pp(&s->mmap), fp);
        fclose(fp);
    }

    s->vtable->sk->plaintext_field(&s->mmap, s->field);
    return 0;
}

void
obf_clear(struct state *s)
{
    if (s) {
        s->vtable->sk->clear(&s->mmap);
        aes_randclear(s->rand);
        thpool_destroy(s->thpool);
    }
    free(s);
}

void
obf_encode_layer(struct state *s, long idx, long inp, long nrows, long ncols,
                 fmpz_mat_t zero, fmpz_mat_t one)
{
    struct write_layer_s *wl_s;
    mmap_enc *zero_enc, *one_enc;
    double start;
    char idx_s[10];
    const mmap_pp *pp = s->vtable->sk->pp(&s->mmap);

    start = current_time();

    (void) snprintf(idx_s, 10, "%ld", idx);

    zero_enc = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
    one_enc = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
    for (ssize_t i = 0; i < nrows * ncols; ++i) {
        s->vtable->enc->init(&zero_enc[i], pp);
        s->vtable->enc->init(&one_enc[i], pp);
    }

    wl_s = (struct write_layer_s *) malloc(sizeof(write_layer_s));
    wl_s->vtable = s->vtable;
    wl_s->dir = s->dir;
    wl_s->zero_enc = zero_enc;
    wl_s->one_enc = one_enc;
    wl_s->inp = inp;
    wl_s->idx = idx;
    wl_s->nrows = nrows;
    wl_s->ncols = ncols;
    wl_s->start = start;

    (void) thpool_add_tag(s->thpool, idx_s, 2 * nrows * ncols,
                          thpool_write_layer, wl_s);

    for (long c = 0; c < 2; ++c) {
        for (long i = 0; i < nrows; ++i) {
            for (long j = 0; j < ncols; ++j) {
                fmpz_t *plaintext;
                mmap_enc *enc;
                struct encode_elem_s *args;

                plaintext = (fmpz_t *) malloc(sizeof(fmpz_t));
                if (c == 0) {
                    fmpz_init_set(*plaintext, fmpz_mat_entry(zero, i, j));
                    enc = &zero_enc[i * nrows + j];
                } else {
                    fmpz_init_set(*plaintext, fmpz_mat_entry(one, i, j));
                    enc = &one_enc[i * nrows + j];
                }
                args = (struct encode_elem_s *) malloc(sizeof(struct encode_elem_s));
                args->vtable = s->vtable;
                args->sk = &s->mmap;
                args->enc = enc;
                args->n = 1;
                args->plaintext = plaintext;
                args->group = (int *) calloc(s->nzs , sizeof(int));
                if (s->type == MMAP_CLT)
                    args->group[idx] = 1;
                thpool_add_work(s->thpool, thpool_encode_elem, (void *) args, idx_s);
            }
        }
    }
}

int
obf_evaluate(enum mmap_e type, char *dir, char *input, unsigned long bplen,
             unsigned long ncores)
{
    const mmap_vtable *vtable;
    mmap_pp pp;
    char fname[100];
    FILE *fp;
    mmap_enc *result = NULL;
    long nrows, ncols, nrows_prev;
    int err = 0, iszero = -1;
    double start, end;

    (void) omp_set_num_threads(ncores);

    switch (type) {
    case MMAP_CLT:
        vtable = &clt13_vtable;
        break;
    case MMAP_GGHLITE:
        vtable = &gghlite_vtable;
        break;
    default:
        return 1;
    }

    (void) snprintf(fname, 100, "%s/params", dir);
    fp = fopen(fname, "r+b");
    vtable->pp->fread(&pp, fp);
    fclose(fp);

    for (unsigned long layer = 0; layer < bplen; ++layer) {
        unsigned int inp;
        mmap_enc *left, *right;

        start = current_time();

        // determine the size of the matrix
        (void) snprintf(fname, 100, "%s/%ld.nrows", dir, layer);
        fp = fopen(fname, "r+b");
        fread(&nrows, sizeof nrows, 1, fp);
        fclose(fp);

        (void) snprintf(fname, 100, "%s/%ld.ncols", dir, layer);
        fp = fopen(fname, "r+b");
        fread(&ncols, sizeof ncols, 1, fp);
        fclose(fp);

        // find out the input bit for the given layer
        (void) snprintf(fname, 100, "%s/%ld.input", dir, layer);
        fp = fopen(fname, "r+b");
        fread(&inp, sizeof inp, 1, fp);
        fclose(fp);

        if (inp >= strlen(input)) {
            fprintf(stderr, "Error: invalid input: %d >= %ld\n", inp, strlen(input));
            err = 1;
            break;
        }
        if (input[inp] != '0' && input[inp] != '1') {
            fprintf(stderr, "Error: input must be 0 or 1, got %d\n", input[inp]);
            err = 1;
            break;
        }
        // load in appropriate matrix for the given input value
        if (input[inp] == '0') {
            (void) snprintf(fname, 100, "%s/%ld.zero", dir, layer);
        } else {
            (void) snprintf(fname, 100, "%s/%ld.one", dir, layer);
        }

        if (layer == 0) {
            result = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
            fp = fopen(fname, "r+b");
            for (int i = 0; i < nrows * ncols; ++i) {
                vtable->enc->init(&result[i], &pp);
                vtable->enc->fread(&result[i], fp);
            }
            fclose(fp);
            nrows_prev = nrows;
        } else {
            left = result;
            right = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
            fp = fopen(fname, "r+b");
            for (int i = 0; i < nrows * ncols; ++i) {
                vtable->enc->init(&right[i], &pp);
                vtable->enc->fread(&right[i], fp);
            }
            fclose(fp);
            result = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows_prev * ncols);
            for (int i = 0; i < nrows_prev * ncols; ++i) {
                vtable->enc->init(&result[i], &pp);
            }
            mult_matrices(vtable, &pp, result, left, right, nrows_prev, nrows, ncols);
            for (int i = 0; i < nrows_prev * nrows; ++i) {
                vtable->enc->clear(&left[i]);
            }
            for (int i = 0; i < nrows * ncols; ++i) {
                vtable->enc->clear(&right[i]);
            }
            free(left);
            free(right);
        }
        end = current_time();

        if (g_verbose)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n", end - start);
    }

    if (!err) {
        start = current_time();
        iszero = vtable->enc->is_zero(&result[1], &pp);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Zero test: %f\n", end - start);
    }

    if (result) {
        for (int i = 0; i < nrows_prev * ncols; ++i) {
            vtable->enc->clear(&result[i]);
        }
        free(result);
    }

    return iszero;
}
