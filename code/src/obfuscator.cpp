#include "obfuscator.h"
#include "thpool_fns.h"

#include <mife/mife_internals.h>
#include <mmap/mmap_clt.h>
#include <mmap/mmap_gghlite.h>
#include <omp.h>

typedef struct obf_state_s {
    threadpool thpool;
    unsigned long secparam;
    enum mmap_e type;
    mmap_sk mmap;
    const mmap_vtable *vtable;
    aes_randstate_t rand;
    const char *dir;
    long nzs;
    fmpz_mat_t *randomizer;
} obf_state_t;

fmpz_mat_t g_first;
fmpz_mat_t g_second;

obf_state_t *
obf_init(enum mmap_e type, const char *dir,
         unsigned long secparam, unsigned long kappa, unsigned long nzs,
         unsigned long nthreads, unsigned long ncores)
{
    obf_state_t *s = NULL;

    s = (obf_state_t *) calloc(1, sizeof(obf_state_t));
    if (s == NULL)
        return NULL;

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
        return NULL;
    }

    nthreads = 1;               // FIXME: remove

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

    return s;
}

void
obf_clear(obf_state_t *s)
{
    if (s) {
        s->vtable->sk->clear(&s->mmap);
        aes_randclear(s->rand);
        thpool_destroy(s->thpool);
    }
    free(s);
}

void
obf_get_field(obf_state_t *s, fmpz_t *field)
{
    s->vtable->sk->plaintext_field(&s->mmap, *field);
}

static void
_fmpz_mat_init_rand(fmpz_mat_t m, long n, aes_randstate_t rand, fmpz_t field)
{
    fmpz_mat_t inverse;

    fmpz_mat_init(m, n, n);
    fmpz_mat_init(inverse, n, n);
    while (true) {
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                fmpz_randm_aes(fmpz_mat_entry(m, i, j), rand, field);
            }
        }
        fmpz_modp_matrix_inverse(inverse, m, n, field);
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if (!fmpz_is_zero(fmpz_mat_entry(inverse, i, j)))
                    goto done;
            }
        }
    }
done:
    fmpz_mat_clear(inverse);
}

void
obf_randomize_layer(obf_state_t *s, long nrows, long ncols, int type,
                    fmpz_mat_t zero, fmpz_mat_t one)
{
    fmpz_t rand, field;

    fmpz_init(rand);
    fmpz_init(field);

    s->vtable->sk->plaintext_field(&s->mmap, field);
    fmpz_randm_aes(rand, s->rand, field);

    if (type & ENCODE_LAYER_TYPE_FIRST && type & ENCODE_LAYER_TYPE_LAST) {

    } else if (type & ENCODE_LAYER_TYPE_FIRST) {
        s->randomizer = (fmpz_mat_t *) malloc(sizeof(fmpz_mat_t));
        _fmpz_mat_init_rand(*s->randomizer, ncols, s->rand, field);
        fmpz_mat_mul(zero, zero, *s->randomizer);
        fmpz_mat_scalar_mod_fmpz(zero, zero, field);
        fmpz_mat_mul(one, one, *s->randomizer);
        fmpz_mat_scalar_mod_fmpz(one, one, field);

        // fmpz_mat_init(g_first, nrows, ncols);
        // fmpz_mat_set(g_first, one);
    } else if (type & ENCODE_LAYER_TYPE_MIDDLE) {
        fmpz_modp_matrix_inverse(*s->randomizer, *s->randomizer, nrows, field);
        fmpz_mat_mul(zero, *s->randomizer, zero);
        fmpz_mat_scalar_mod_fmpz(zero, zero, field);
        fmpz_mat_mul(one, *s->randomizer, one);
        fmpz_mat_scalar_mod_fmpz(one, one, field);
        fmpz_mat_clear(*s->randomizer);
        free(s->randomizer);
        s->randomizer = (fmpz_mat_t *) malloc(sizeof(fmpz_mat_t));
        _fmpz_mat_init_rand(*s->randomizer, ncols, s->rand, field);
        fmpz_mat_mul(zero, zero, *s->randomizer);
        fmpz_mat_scalar_mod_fmpz(zero, zero, field);
        fmpz_mat_mul(one, one, *s->randomizer);
        fmpz_mat_scalar_mod_fmpz(one, one, field);
    } else if (type & ENCODE_LAYER_TYPE_LAST) {
        fmpz_modp_matrix_inverse(*s->randomizer, *s->randomizer, nrows, field);

        fmpz_mat_mul(zero, *s->randomizer, zero);
        fmpz_mat_scalar_mod_fmpz(zero, zero, field);
        fmpz_mat_mul(one, *s->randomizer, one);
        fmpz_mat_scalar_mod_fmpz(one, one, field);

        // fprintf(stderr, "RESULT:\n");
        // fmpz_mat_init(g_second, nrows, ncols);
        // fmpz_mat_set(g_second, one);
        // fmpz_mat_mul(g_first, g_first, g_second);
        // fmpz_mat_scalar_mod_fmpz(g_first, g_first, field);
        // fmpz_mat_fprint_pretty(stderr, g_first);
        // fprintf(stderr, "\n\n");
    }

    fmpz_clear(rand);
    fmpz_clear(field);
}

void
obf_encode_layer(obf_state_t *s, long idx, long inp, long nrows, long ncols,
                 int type, fmpz_mat_t zero, fmpz_mat_t one)
{
    struct write_layer_s *wl_s;
    mmap_enc_mat_t *zero_enc, *one_enc;
    double start;
    char idx_s[10];
    const mmap_pp *pp = s->vtable->sk->pp(&s->mmap);

    start = current_time();

    (void) snprintf(idx_s, 10, "%ld", idx);

    zero_enc = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
    one_enc = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
    mmap_enc_mat_init(s->vtable, pp, *zero_enc, nrows, ncols);
    mmap_enc_mat_init(s->vtable, pp, *one_enc, nrows, ncols);

    obf_randomize_layer(s, nrows, ncols, type, zero, one);

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
                    enc = zero_enc[0]->m[i][j];
                } else {
                    fmpz_init_set(*plaintext, fmpz_mat_entry(one, i, j));
                    enc = one_enc[0]->m[i][j];
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
    mmap_enc_mat_t *result = NULL;
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
        mmap_enc_mat_t *left, *right;

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
            result = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
            mmap_enc_mat_init(vtable, &pp, *result, nrows, ncols);
            fp = fopen(fname, "r+b");
            for (int i = 0; i < nrows; ++i) {
                for (int j = 0; j < ncols; ++j) {
                    vtable->enc->fread(result[0]->m[i][j], fp);
                }
            }
            fclose(fp);
            nrows_prev = nrows;
        } else {
            left = result;
            right = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
            mmap_enc_mat_init(vtable, &pp, *right, nrows, ncols);
            fp = fopen(fname, "r+b");
            for (int i = 0; i < nrows; ++i) {
                for (int j = 0; j < ncols; ++j) {
                    vtable->enc->fread(right[0]->m[i][j], fp);
                }
            }
            fclose(fp);

            result = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
            mmap_enc_mat_init(vtable, &pp, *result, nrows_prev, ncols);
            mmap_enc_mat_mul(vtable, &pp, *result, *left, *right);
            mmap_enc_mat_clear(vtable, *left);
            mmap_enc_mat_clear(vtable, *right);
            free(left);
            free(right);
        }

        end = current_time();

        if (layer && g_verbose)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n", end - start);
    }

    if (!err) {
        start = current_time();
        iszero = vtable->enc->is_zero(result[0]->m[0][1], &pp);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Zero test: %f\n", end - start);
    }

    if (result) {
        mmap_enc_mat_clear(vtable, *result);
        free(result);
    }

    return iszero;
}

void
obf_wait(obf_state_t *s)
{
    thpool_wait(s->thpool);
}
