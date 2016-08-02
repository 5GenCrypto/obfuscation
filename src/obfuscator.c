#include "obfuscator.h"
#include "thpool.h"
#include "thpool_fns.h"

#include <mmap/mmap_clt.h>
#include <mmap/mmap_gghlite.h>
#include <mmap/mmap_dummy.h>
#include <omp.h>

typedef struct obf_state_s {
    threadpool thpool;
    uint64_t secparam;
    enum mmap_e type;
    mmap_sk mmap;
    const mmap_vtable *vtable;
    aes_randstate_t rand;
    const char *dir;
    uint64_t nzs;
    fmpz_mat_t *randomizer;
    fmpz_mat_t *inverse;
    uint64_t flags;
} obf_state_t;


static FILE *
open_indexed_file(const char *dir, const char *file, uint64_t index,
                  const char *mode)
{
    char fname[50];

    (void) snprintf(fname, 50, "%lu.%s", index, file);
    return open_file(dir, fname, mode);
}

obf_state_t *
obf_init(enum mmap_e type, const char *dir, uint64_t secparam, uint64_t kappa,
         uint64_t nzs, uint64_t nthreads, uint64_t ncores, uint64_t flags)
{
    obf_state_t *s = NULL;

    if (secparam == 0 || kappa == 0 || nzs == 0)
        return NULL;

    s = calloc(1, sizeof(obf_state_t));
    if (s == NULL)
        return NULL;

    s->secparam = secparam;
    s->type = type;
    s->dir = dir;
    s->nzs = nzs;
    s->flags = flags;
    s->randomizer = malloc(sizeof(fmpz_mat_t));
    s->inverse = malloc(sizeof(fmpz_mat_t));

    switch (s->type) {
    case MMAP_DUMMY:
        s->vtable = &dummy_vtable;
        break;
    case MMAP_CLT:
        s->vtable = &clt_vtable;
        break;
    case MMAP_GGHLITE:
        s->vtable = &gghlite_vtable;
        break;
    default:
        return NULL;
    }

    (void) aes_randinit(s->rand);
    if (nthreads == 0)
        nthreads = ncores;
    s->thpool = thpool_init(nthreads);
    (void) omp_set_num_threads(ncores);

    if (s->flags & OBFUSCATOR_FLAG_VERBOSE) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
        if (s->flags & OBFUSCATOR_FLAG_DUAL_INPUT_BP)
            fprintf(stderr, "  Using dual input branching programs\n");
        if (s->flags & OBFUSCATOR_FLAG_NO_RANDOMIZATION)
            fprintf(stderr, "  Not randomizing branching programs\n");
    }

    s->vtable->sk->init(&s->mmap, secparam, kappa, nzs, s->rand,
                        s->flags & OBFUSCATOR_FLAG_VERBOSE);
    {
        FILE *fp = open_file(dir, "params", "w+b");
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
        free(s->randomizer);
        free(s->inverse);
        thpool_destroy(s->thpool);
    }
    free(s);
}

static void
_fmpz_mat_init_square_rand(obf_state_t *s, fmpz_mat_t mat, fmpz_mat_t inverse,
                           long n, aes_randstate_t rand, fmpz_t field)
{
    int nzeros;
    do {
        nzeros = 0;
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                fmpz_randm_aes(fmpz_mat_entry(mat, i, j), rand, field);
            }
        }
        if (s->flags & OBFUSCATOR_FLAG_VERBOSE)
            fprintf(stderr, "    Finding inverse...\n");
        fmpz_modp_matrix_inverse(inverse, mat, n, field);
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if (fmpz_is_zero(fmpz_mat_entry(inverse, i, j)))
                    nzeros++;
            }
        }
    } while (nzeros == n * n);
}

static void
_fmpz_mat_init_diagonal_rand(fmpz_mat_t mat, long n, aes_randstate_t rand,
                             fmpz_t field)
{
    fmpz_mat_init(mat, n, n);
    fmpz_mat_one(mat);
    for (int i = 0; i < n; i++) {
        do {
            fmpz_randm_aes(fmpz_mat_entry(mat, i, i), rand, field);
        } while (fmpz_cmp_ui(fmpz_mat_entry(mat, i, i), 0) == 0);
    }
}

static inline void
fmpz_mat_mul_mod(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, fmpz_t p)
{
    fmpz_mat_mul(a, b, c);
    fmpz_mat_scalar_mod_fmpz(a, a, p);
}

static inline void
fmpz_layer_mul_left(uint64_t n, fmpz_mat_t *mats, fmpz_mat_t m, fmpz_t p)
{
    for (uint64_t i = 0; i < n; ++i) {
        fmpz_mat_mul_mod(mats[i], m, mats[i], p);
    }
}

static void
fmpz_layer_mul_right(uint64_t n, fmpz_mat_t *mats, fmpz_mat_t m, fmpz_t p)
{
    for (uint64_t i = 0; i < n; ++i) {
        fmpz_mat_mul_mod(mats[i], mats[i], m, p);
    }
}

static void
obf_randomize_layer(obf_state_t *s, long nrows, long ncols,
                    encode_layer_randomization_flag_t rflag,
                    uint64_t n, fmpz_mat_t *mats)
{
    fmpz_t field;

    fmpz_init(field);

    s->vtable->sk->plaintext_field(&s->mmap, field);

    if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST) {
        fmpz_mat_t first;
        _fmpz_mat_init_diagonal_rand(first, nrows, s->rand, field);
        fmpz_layer_mul_left(n, mats, first, field);
        fmpz_mat_clear(first);
    }
    if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_LAST) {
        fmpz_mat_t last;
        _fmpz_mat_init_diagonal_rand(last, ncols, s->rand, field);
        fmpz_layer_mul_right(n, mats, last, field);
        fmpz_mat_clear(last);
    }

    if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST
        && rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_LAST) {
    } else if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST) {
        fmpz_mat_init(*s->randomizer, ncols, ncols);
        fmpz_mat_init(*s->inverse, ncols, ncols);
        _fmpz_mat_init_square_rand(s, *s->randomizer, *s->inverse, ncols,
                                   s->rand, field);
        fmpz_layer_mul_right(n, mats, *s->randomizer, field);
    } else if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_MIDDLE) {
        fmpz_layer_mul_left(n, mats, *s->inverse, field);
        fmpz_mat_clear(*s->randomizer);
        fmpz_mat_clear(*s->inverse);

        fmpz_mat_init(*s->randomizer, ncols, ncols);
        fmpz_mat_init(*s->inverse, ncols, ncols);
        _fmpz_mat_init_square_rand(s, *s->randomizer, *s->inverse, ncols,
                                   s->rand, field);
        fmpz_layer_mul_right(n, mats, *s->randomizer, field);
    } else if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_LAST) {
        fmpz_layer_mul_left(n, mats, *s->inverse, field);
        fmpz_mat_clear(*s->randomizer);
        fmpz_mat_clear(*s->inverse);
    }

    {
        fmpz_t alpha;
        fmpz_init(alpha);
        for (uint64_t i = 0; i < n; ++i) {
            do {
                fmpz_randm_aes(alpha, s->rand, field);
            } while (fmpz_cmp_ui(alpha, 0) == 0);
            fmpz_mat_scalar_mul_fmpz(mats[i], mats[i], alpha);
        }
        fmpz_clear(alpha);
    }

    fmpz_clear(field);
}

static int
add_work_write_layer(obf_state_t *s, uint64_t n, long inp, long idx,
                     long nrows, long ncols, char **names, char *tag,
                     mmap_enc_mat_t **enc_mats)
{
    struct write_layer_s *wl_s;
    wl_s = malloc(sizeof(struct write_layer_s));
    wl_s->vtable = s->vtable;
    wl_s->dir = s->dir;
    wl_s->n = n;
    wl_s->enc_mats = enc_mats;
    wl_s->names = names;
    wl_s->inp = inp;
    wl_s->idx = idx;
    wl_s->nrows = nrows;
    wl_s->ncols = ncols;
    wl_s->start = current_time();
    wl_s->verbose = s->flags & OBFUSCATOR_FLAG_VERBOSE;
    if (thpool_add_tag(s->thpool, tag, n * nrows * ncols,
                       thpool_write_layer, wl_s) == OBFUSCATOR_ERR) {
        free(wl_s);
        return OBFUSCATOR_ERR;
    }
    return OBFUSCATOR_OK;
}

static void
add_work(obf_state_t *s, fmpz_mat_t *mats, mmap_enc_mat_t **enc_mats,
         int **pows, long c, long i, long j, char *tag)
{
    fmpz_t *plaintext;
    mmap_enc *enc;
    struct encode_elem_s *args;

    plaintext = malloc(sizeof(fmpz_t));
    args = malloc(sizeof(struct encode_elem_s));
    fmpz_init_set(*plaintext, fmpz_mat_entry(mats[c], i, j));
    enc = enc_mats[c][0]->m[i][j];
    args->group = pows[c];
    args->n = 1;
    args->plaintext = plaintext;
    args->vtable = s->vtable;
    args->sk = &s->mmap;
    args->enc = enc;

    thpool_add_work(s->thpool, thpool_encode_elem, (void *) args, tag);
}

int
obf_encode_layer(obf_state_t *s, uint64_t n, int **pows, fmpz_mat_t *mats,
                 long idx, long inp, encode_layer_randomization_flag_t rflag)
{
    char tag[10];
    mmap_enc_mat_t **enc_mats;
    char **names;
    const mmap_pp *pp = s->vtable->sk->pp(&s->mmap);
    long nrows, ncols;

    /* TODO: check for mismatched matrices */

    nrows = mats[0]->r;
    ncols = mats[0]->c;

    (void) snprintf(tag, 10, "%ld", idx);

    if (!(s->flags & OBFUSCATOR_FLAG_NO_RANDOMIZATION)) {
        double start, end;
        start = current_time();
        obf_randomize_layer(s, nrows, ncols, rflag, n, mats);
        end = current_time();
        if (s->flags & OBFUSCATOR_FLAG_VERBOSE)
            (void) fprintf(stderr, "  Randomizing matrix: %f\n", end - start);
    }

    enc_mats = calloc(n, sizeof(mmap_enc_mat_t *));
    for (uint64_t c = 0; c < n; ++c) {
        enc_mats[c] = malloc(sizeof(mmap_enc_mat_t));
        mmap_enc_mat_init(s->vtable, pp, *enc_mats[c], nrows, ncols);
    }
    names = calloc(n, sizeof(char *));
    for (uint64_t c = 0; c < n; ++c) {
        names[c] = calloc(10, sizeof(char));
        (void) snprintf(names[c], 10, "%lu", c);
    }

    if (add_work_write_layer(s, n, inp, idx, nrows, ncols, names, tag, enc_mats)
        == OBFUSCATOR_ERR)
        return OBFUSCATOR_ERR;

    for (uint64_t c = 0; c < n; ++c) {
        for (long i = 0; i < nrows; ++i) {
            for (long j = 0; j < ncols; ++j) {
                add_work(s, mats, enc_mats, pows, c, i, j, tag);
            }
        }
    }

    return OBFUSCATOR_OK;
}

int
obf_evaluate(enum mmap_e type, char *dir, uint64_t len, uint64_t *input,
             uint64_t bplen, uint64_t ncores, bool verbose)
{
    const mmap_vtable *vtable;
    mmap_pp pp;
    FILE *fp;
    mmap_enc_mat_t *result = NULL;
    uint64_t nrows, ncols, nrows_prev = 0;
    int err = 1, iszero = -1;
    double start, end;

    switch (type) {
    case MMAP_DUMMY:
        vtable = &dummy_vtable;
        break;
    case MMAP_CLT:
        vtable = &clt_vtable;
        break;
    case MMAP_GGHLITE:
        vtable = &gghlite_vtable;
        break;
    default:
        return err;
    }

    (void) omp_set_num_threads(ncores);

    if ((fp = open_file(dir, "params", "r+b")) == NULL)
        goto done;
    vtable->pp->fread(&pp, fp);
    fclose(fp);

    for (uint64_t layer = 0; layer < bplen; ++layer) {
        uint64_t inp;
        char str[10];
        mmap_enc_mat_t *left, *right;

        start = current_time();

        // determine the size of the matrix
        if ((fp = open_indexed_file(dir, "nrows", layer, "r+b")) == NULL)
            goto done;
        fread(&nrows, sizeof nrows, 1, fp);
        fclose(fp);
        if ((fp = open_indexed_file(dir, "ncols", layer, "r+b")) == NULL)
            goto done;
        fread(&ncols, sizeof ncols, 1, fp);
        fclose(fp);

        // find out the input bit for the given layer
        if ((fp = open_indexed_file(dir, "input", layer, "r+b")) == NULL)
            goto done;
        fread(&inp, sizeof inp, 1, fp);
        fclose(fp);

        if (inp > len) {
            fprintf(stderr, "invalid input: %lu > %ld\n", inp, len);
            goto done;
        }
        // load in appropriate matrix for the given input value
        (void) snprintf(str, 10, "%lu", input[inp]);
        if ((fp = open_indexed_file(dir, str, layer, "r+b")) == NULL)
            goto done;

        if (layer == 0) {
            result = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
            mmap_enc_mat_init(vtable, &pp, *result, nrows, ncols);
            for (uint64_t i = 0; i < nrows; ++i) {
                for (uint64_t j = 0; j < ncols; ++j) {
                    vtable->enc->fread(result[0]->m[i][j], fp);
                }
            }
            nrows_prev = nrows;
        } else {
            left = result;
            right = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
            mmap_enc_mat_init(vtable, &pp, *right, nrows, ncols);
            for (uint64_t i = 0; i < nrows; ++i) {
                for (uint64_t j = 0; j < ncols; ++j) {
                    vtable->enc->fread(right[0]->m[i][j], fp);
                }
            }

            result = (mmap_enc_mat_t *) malloc(sizeof(mmap_enc_mat_t));
            mmap_enc_mat_init(vtable, &pp, *result, nrows_prev, ncols);
            mmap_enc_mat_mul_par(vtable, &pp, *result, *left, *right);
            mmap_enc_mat_clear(vtable, *left);
            mmap_enc_mat_clear(vtable, *right);
            free(left);
            free(right);
        }

        fclose(fp);

        end = current_time();

        if (verbose && layer != 0)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n", end - start);
    }
    err = 0;

done:
    if (!err) {
        start = current_time();
        if (result[0]->nrows == 1 && result[0]->ncols == 1)
            iszero = vtable->enc->is_zero(result[0]->m[0][0], &pp);
        else
            iszero = vtable->enc->is_zero(result[0]->m[0][1], &pp);
        end = current_time();
        if (verbose)
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
