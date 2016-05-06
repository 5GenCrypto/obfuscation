#include "obfuscator.h"
#include "thpool.h"
#include "thpool_fns.h"

#include <mmap/mmap_clt.h>
#include <mmap/mmap_gghlite.h>
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

    s = (obf_state_t *) calloc(1, sizeof(obf_state_t));
    if (s == NULL)
        return NULL;

    s->secparam = secparam;
    s->type = type;
    s->dir = dir;
    s->nzs = nzs;
    s->flags = flags;

    switch (s->type) {
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
    s->thpool = thpool_init(nthreads);
    (void) omp_set_num_threads(ncores);

    if (s->flags & OBFUSCATOR_FLAG_VERBOSE) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
        if (s->flags & OBFUSCATOR_FLAG_DUAL_INPUT_BP)
            fprintf(stderr, "  Using dual input branching programs\n");
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
        if (s->randomizer)
            free(s->randomizer);
        thpool_destroy(s->thpool);
    }
    free(s);
}

static void
_fmpz_mat_init_square_rand(fmpz_mat_t m, long n, aes_randstate_t rand,
                           fmpz_t field)
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

static void
_fmpz_mat_init_rand(fmpz_mat_t mat, long m, long n, aes_randstate_t rand,
                    fmpz_t field)
{
    fmpz_mat_init(mat, m, n);
    for (int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            fmpz_randm_aes(fmpz_mat_entry(mat, i, j), rand, field);
        }
    }
}

static void
_fmpz_mat_init_diagonal_rand(fmpz_mat_t mat, long n, aes_randstate_t rand,
                             fmpz_t field)
{
    fmpz_mat_init(mat, n, n);
    fmpz_mat_one(mat);
    for (int i = 0; i < n; i++) {
        fmpz_randm_aes(fmpz_mat_entry(mat, i, i), rand, field);
    }
}

static inline void
fmpz_mat_mul_mod(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, fmpz_t p)
{
    fmpz_mat_mul(a, b, c);
    fmpz_mat_scalar_mod_fmpz(a, a, p);
}

static inline void
fmpz_layer_mul_left(fmpz_mat_t zero, fmpz_mat_t one, fmpz_mat_t m, fmpz_t p)
{
    fmpz_mat_mul_mod(zero, m, zero, p);
    fmpz_mat_mul_mod(one, m, one, p);
}

static void
fmpz_layer_mul_right(fmpz_mat_t zero, fmpz_mat_t one, fmpz_mat_t m, fmpz_t p)
{
    fmpz_mat_mul_mod(zero, zero, m, p);
    fmpz_mat_mul_mod(one, one, m, p);
}

void
obf_randomize_layer(obf_state_t *s, long nrows, long ncols,
                    encode_layer_randomization_flag_t rflag,
                    fmpz_mat_t zero, fmpz_mat_t one)
{
    fmpz_t field;

    fmpz_init(field);

    s->vtable->sk->plaintext_field(&s->mmap, field);

    if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST) {
        fmpz_mat_t first;
        _fmpz_mat_init_diagonal_rand(first, nrows, s->rand, field);
        fmpz_layer_mul_left(zero, one, first, field);
        fmpz_mat_clear(first);
    }
    if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_LAST) {
        fmpz_mat_t last;
        _fmpz_mat_init_diagonal_rand(last, ncols, s->rand, field);
        fmpz_layer_mul_right(zero, one, last, field);
        fmpz_mat_clear(last);
    }

    if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST
        && rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_LAST) {
    } else if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST) {
        s->randomizer = (fmpz_mat_t *) malloc(sizeof(fmpz_mat_t));
        _fmpz_mat_init_square_rand(*s->randomizer, ncols, s->rand, field);
        fmpz_layer_mul_right(zero, one, *s->randomizer, field);
    } else if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_MIDDLE) {
        fmpz_modp_matrix_inverse(*s->randomizer, *s->randomizer, nrows, field);
        fmpz_layer_mul_left(zero, one, *s->randomizer, field);
        fmpz_mat_clear(*s->randomizer);
        free(s->randomizer);
        s->randomizer = (fmpz_mat_t *) malloc(sizeof(fmpz_mat_t));
        _fmpz_mat_init_square_rand(*s->randomizer, ncols, s->rand, field);
        fmpz_layer_mul_right(zero, one, *s->randomizer, field);
    } else if (rflag & ENCODE_LAYER_RANDOMIZATION_TYPE_LAST) {
        fmpz_modp_matrix_inverse(*s->randomizer, *s->randomizer, nrows, field);
        fmpz_layer_mul_left(zero, one, *s->randomizer, field);
        fmpz_mat_clear(*s->randomizer);
        free(s->randomizer);
        s->randomizer = NULL;
    }

    {
        fmpz_t alpha;
        fmpz_init(alpha);
        fmpz_randm_aes(alpha, s->rand, field);
        fmpz_mat_scalar_mul_fmpz(zero, zero, alpha);
        fmpz_randm_aes(alpha, s->rand, field);
        fmpz_mat_scalar_mul_fmpz(one, one, alpha);
        fmpz_clear(alpha);
    }

    fmpz_clear(field);
}

static void
add_work(obf_state_t *s, fmpz_mat_t m0, fmpz_mat_t m1, mmap_enc_mat_t *m0_enc,
         mmap_enc_mat_t *m1_enc, int *zero_pows, int *one_pows, long i, long j,
         long c, char *tag)
{
    fmpz_t *plaintext;
    mmap_enc *enc;
    struct encode_elem_s *args;

    plaintext = malloc(sizeof(fmpz_t));
    args = malloc(sizeof(struct encode_elem_s));
    if (c == 0) {
        fmpz_init_set(*plaintext, fmpz_mat_entry(m0, i, j));
        enc = m0_enc[0]->m[i][j];
        args->group = zero_pows;
    } else {
        fmpz_init_set(*plaintext, fmpz_mat_entry(m1, i, j));
        enc = m1_enc[0]->m[i][j];
        args->group = one_pows;
    }
    args->n = 1;
    args->plaintext = plaintext;
    args->vtable = s->vtable;
    args->sk = &s->mmap;
    args->enc = enc;

    thpool_add_work(s->thpool, thpool_encode_elem, (void *) args, tag);
}

int
obf_encode_layer(obf_state_t *s, long idx, long inp, long nrows, long ncols,
                 encode_layer_randomization_flag_t rflag, int *zero_pows,
                 int *one_pows, fmpz_mat_t zero, fmpz_mat_t one)
{
    char tag[10], di_tag[10];
    mmap_enc_mat_t *zero_enc, *one_enc;
    double start, end;
    const mmap_pp *pp = s->vtable->sk->pp(&s->mmap);
    /* For dual input BPs */
    fmpz_mat_t zero_one, one_zero;
    mmap_enc_mat_t *zero_one_enc = NULL, *one_zero_enc = NULL;

    (void) snprintf(tag, 10, "%ld", idx);
    (void) snprintf(di_tag, 10, "%ld_di", idx);

    if (s->flags & OBFUSCATOR_FLAG_DUAL_INPUT_BP) {
        fmpz_t field;
        fmpz_init(field);
        s->vtable->sk->plaintext_field(&s->mmap, field);
        _fmpz_mat_init_rand(zero_one, nrows, ncols, s->rand, field);
        _fmpz_mat_init_rand(one_zero, nrows, ncols, s->rand, field);
        fmpz_clear(field);
    }

    if (!(s->flags & OBFUSCATOR_FLAG_NO_RANDOMIZATION)) {
        start = current_time();
        obf_randomize_layer(s, nrows, ncols, rflag, zero, one);
        end = current_time();
        if (s->flags & OBFUSCATOR_FLAG_VERBOSE)
            (void) fprintf(stderr, "  Randomizing matrix: %f\n", end - start);
    }

    start = current_time();

    zero_enc = malloc(sizeof(mmap_enc_mat_t));
    one_enc = malloc(sizeof(mmap_enc_mat_t));
    mmap_enc_mat_init(s->vtable, pp, *zero_enc, nrows, ncols);
    mmap_enc_mat_init(s->vtable, pp, *one_enc, nrows, ncols);
    if (s->flags & OBFUSCATOR_FLAG_DUAL_INPUT_BP) {
        zero_one_enc = malloc(sizeof(mmap_enc_mat_t));
        one_zero_enc = malloc(sizeof(mmap_enc_mat_t));
        mmap_enc_mat_init(s->vtable, pp, *zero_one_enc, nrows, ncols);
        mmap_enc_mat_init(s->vtable, pp, *one_zero_enc, nrows, ncols);
    }

    {
        struct write_layer_s *wl_s;
        wl_s = malloc(sizeof(struct write_layer_s));
        wl_s->vtable = s->vtable;
        wl_s->dir = s->dir;
        wl_s->zero_enc = zero_enc;
        wl_s->one_enc = one_enc;
        wl_s->inp = inp;
        wl_s->idx = idx;
        wl_s->nrows = nrows;
        wl_s->ncols = ncols;
        wl_s->zero = "zero";
        wl_s->one = "one";
        wl_s->start = start;
        wl_s->verbose = s->flags & OBFUSCATOR_FLAG_VERBOSE;
        if (thpool_add_tag(s->thpool, tag, 2 * nrows * ncols,
                           thpool_write_layer, wl_s) == -1) {
            mmap_enc_mat_clear(s->vtable, *zero_enc);
            mmap_enc_mat_clear(s->vtable, *one_enc);
            free(zero_enc);
            free(one_enc);
            free(wl_s);
            return OBFUSCATOR_ERR;
        }
    }

    if (s->flags & OBFUSCATOR_FLAG_DUAL_INPUT_BP) {
        struct write_layer_s *wl_s;
        wl_s = malloc(sizeof(struct write_layer_s));
        wl_s->vtable = s->vtable;
        wl_s->dir = s->dir;
        wl_s->zero_enc = zero_one_enc;
        wl_s->one_enc = one_zero_enc;
        wl_s->inp = inp;
        wl_s->idx = idx;
        wl_s->nrows = nrows;
        wl_s->ncols = ncols;
        wl_s->zero = "zero_one";
        wl_s->one = "one_zero";
        wl_s->start = start;
        wl_s->verbose = s->flags & OBFUSCATOR_FLAG_VERBOSE;
        if (thpool_add_tag(s->thpool, di_tag, 2 * nrows * ncols,
                           thpool_write_layer, wl_s) == -1) {
            mmap_enc_mat_clear(s->vtable, *zero_one_enc);
            mmap_enc_mat_clear(s->vtable, *one_zero_enc);
            free(zero_one_enc);
            free(one_zero_enc);
            free(wl_s);
            return OBFUSCATOR_ERR;
        }
    }

    for (long c = 0; c < 2; ++c) {
        for (long i = 0; i < nrows; ++i) {
            for (long j = 0; j < ncols; ++j) {
                add_work(s, zero, one, zero_enc, one_enc, zero_pows, one_pows,
                         i, j, c, tag);
                if (s->flags & OBFUSCATOR_FLAG_DUAL_INPUT_BP) {
                    add_work(s, zero_one, one_zero, zero_one_enc, one_zero_enc,
                             zero_pows, one_pows, i, j, c, di_tag);
                }
            }
        }
    }

    return OBFUSCATOR_OK;
}

int
obf_evaluate(enum mmap_e type, char *dir, char *input, uint64_t bplen,
             uint64_t ncores, bool verbose)
{
    const mmap_vtable *vtable;
    mmap_pp pp;
    FILE *fp;
    mmap_enc_mat_t *result = NULL;
    uint64_t nrows, ncols, nrows_prev = 0;
    int err = 1, iszero = -1;
    double start, end;

    switch (type) {
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

        if (inp >= strlen(input)) {
            fprintf(stderr, "invalid input: %lu >= %ld\n", inp,
                    strlen(input));
            goto done;
        }
        // load in appropriate matrix for the given input value
        switch(input[inp]) {
        case '0':
            if ((fp = open_indexed_file(dir, "zero", layer, "r+b")) == NULL)
                goto done;
            break;
        case '1':
            if ((fp = open_indexed_file(dir, "one", layer, "r+b")) == NULL)
                goto done;
            break;
        default:
            fprintf(stderr, "input must be 0 or 1, got %d\n", input[inp]);
            goto done;
        }

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
            mmap_enc_mat_mul(vtable, &pp, *result, *left, *right);
            mmap_enc_mat_clear(vtable, *left);
            mmap_enc_mat_clear(vtable, *right);
            free(left);
            free(right);
        }

        fclose(fp);

        end = current_time();

        if (verbose)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n", end - start);
    }
    err = 0;

done:
    if (!err) {
        start = current_time();
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
