#include "zobfuscator.h"
#include "thpool.h"
#include "thpool_fns.h"
#include "circuit.h"

#include <clt13.h>
#include <mmap/mmap_clt.h>
#include <gmp.h>
#include <omp.h>

typedef struct zobf_state_s {
    threadpool thpool;
    unsigned long secparam;
    clt_state mmap;
    aes_randstate_t rand;
    const char *dir;
    mpz_t nchk;
    mpz_t nev;
    uint64_t flags;
} zobf_state_t;

zobf_state_t *
zobf_init(const char *dir, uint64_t secparam, uint64_t kappa, uint64_t nzs,
         int *pows, uint64_t nthreads, uint64_t ncores, uint64_t flags)
{
    zobf_state_t *s;
    int clt_flags = CLT_FLAG_DEFAULT | CLT_FLAG_OPT_PARALLEL_ENCODE;

    s = (zobf_state_t *) malloc(sizeof(zobf_state_t));
    if (s == NULL)
        return NULL;

    s->secparam = secparam;
    s->dir = dir;
    (void) aes_randinit(s->rand);
    s->thpool = thpool_init(nthreads);
    (void) omp_set_num_threads(ncores);
    s->flags = flags;

    if (s->flags & ZOBFUSCATOR_FLAG_VERBOSE) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
        clt_flags |= CLT_FLAG_VERBOSE;
    }

    clt_state_init(&s->mmap, kappa, secparam, nzs, pows, clt_flags, s->rand);
    {
        clt_pp pp;
        FILE *fp;

        clt_pp_init(&pp, &s->mmap);
        fp = open_file(s->dir, "params", "w+b");
        clt_pp_fsave(fp, &pp);
        fclose(fp);
        clt_pp_clear(&pp);
    }

    mpz_init_set(s->nev, s->mmap.gs[0]);
    mpz_init_set(s->nchk, s->mmap.gs[1]);

    return s;
}

void
zobf_clear(zobf_state_t *s)
{
    if (s) {
        aes_randclear(s->rand);
        clt_state_clear(&s->mmap);
        thpool_destroy(s->thpool);
        mpz_clears(s->nev, s->nchk, NULL);
    }
    free(s);
}

static void
create_work(zobf_state_t *s, const char *str, int i,
            const mpz_t elem0, const mpz_t elem1, unsigned int num, ...)
{
    va_list vals;

    struct write_element_s *we_s;
    struct encode_elem_s *args;
    mmap_enc *out;
    fmpz_t *plaintext;
    const mmap_pp *pp = clt_vtable.sk->pp((mmap_sk *) &s->mmap);

    out = (mmap_enc *) malloc(sizeof(mmap_enc));
    plaintext = (fmpz_t *) calloc(2, sizeof(fmpz_t));
    clt_vtable.enc->init(out, pp);
    fmpz_init(plaintext[0]);
    fmpz_init(plaintext[1]);
    fmpz_set_mpz(plaintext[0], elem0);
    fmpz_set_mpz(plaintext[1], elem1);

    we_s = (struct write_element_s *) malloc(sizeof(struct write_element_s));
    we_s->dir = s->dir;
    we_s->elem = out;
    we_s->name = (char *) calloc(sizeof(int) + 10, sizeof(char));
    if (i >= 0) {
        (void) sprintf(we_s->name, str, i);
    } else {
        (void) strcpy(we_s->name, str);
    }
    (void) thpool_add_tag(s->thpool, we_s->name, 1, thpool_write_element, we_s);

    args = (struct encode_elem_s *) malloc(sizeof(struct encode_elem_s));
    args->vtable = &clt_vtable;
    args->sk = (mmap_sk *) &s->mmap;
    args->enc = out;
    args->n = 2;
    args->plaintext = plaintext;
    args->group = (int *) calloc(s->mmap.nzs, sizeof(int));
    va_start(vals, num);
    for (unsigned int j = 0; j < num; ++j) {
        int idx = va_arg(vals, int);
        args->group[idx] = va_arg(vals, int);
    }

    (void) thpool_add_work(s->thpool, thpool_encode_elem, (void *) args,
                           we_s->name);
}

static int
write_element(const char *dir, mmap_enc *elem, const char *name)
{
    FILE *fp;

    fp = open_file(dir, name, "w+b");
    clt_vtable.enc->fwrite(elem, fp);
    fclose(fp);
    return 0;
}

void
zobf_encode_circuit(zobf_state_t *s, const char *circuit, const mpz_t *xs,
                    const mpz_t *ys, const int *xdegs, int ydeg, int n, int m)
{
    mpz_t tmp, c_star;
    mpz_t *alphas, *betas;
    mpz_t zero, one;

    mpz_inits(c_star, tmp, zero, one, NULL);
    mpz_set_ui(zero, 0);
    mpz_set_ui(one, 1);

    alphas = (mpz_t *) malloc(sizeof(mpz_t) * n);
    for (int i = 0; i < n; ++i) {
        mpz_init(alphas[i]);
        mpz_urandomm_aes(alphas[i], s->rand, s->nchk);
    }
    betas = (mpz_t *) malloc(sizeof(mpz_t) * m);
    for (int i = 0; i < m; ++i) {
        mpz_init(betas[i]);
        mpz_urandomm_aes(betas[i], s->rand, s->nchk);
    }

    for (int i = 0; i < n; ++i) {
        mpz_t elems[2];

        mpz_inits(elems[0], elems[1], NULL);

        create_work(s, "x_\%d_0", i, zero, alphas[i], 1, 2 * i, 1);
        create_work(s, "x_\%d_1", i, xs[i], alphas[i], 1, 2 * i + 1, 1);
        create_work(s, "u_\%d_0", i, one, one, 1, 2 * i, 1);
        create_work(s, "u_\%d_1", i, one, one, 1, 2 * i + 1, 1);

        mpz_urandomm_aes(elems[0], s->rand, s->nev);
        mpz_urandomm_aes(elems[1], s->rand, s->nchk);

        create_work(s, "z_\%d_0", i, elems[0], elems[1], 3, 2 * i + 1, xdegs[i],
                    2 * n + i, 1, 3 * n + i, 1);
        create_work(s, "w_\%d_0", i, zero, elems[1], 1, 3 * n + i, 1);

        mpz_urandomm_aes(elems[0], s->rand, s->nev);
        mpz_urandomm_aes(elems[1], s->rand, s->nchk);

        create_work(s, "z_\%d_1", i, elems[0], elems[1], 3, 2 * i, xdegs[i],
                    2 * n + i, 1, 3 * n + i, 1);
        create_work(s, "w_\%d_1", i, zero, elems[1], 1, 3 * n + i, 1);

        mpz_clears(elems[0], elems[1], NULL);
    }

    for (int i = 0; i < m; ++i) {
        mpz_t elems[2];
        mpz_inits(elems[0], elems[1], NULL);

        mpz_set(elems[0], ys[i]);
        if (mpz_sgn(elems[0]) == -1) {
            mpz_mod(elems[0], elems[0], s->nev);
        }
        mpz_set(elems[1], betas[i]);
        create_work(s, "y_\%d", i, elems[0], elems[1], 1, 4 * n, 1);

        mpz_clears(elems[0], elems[1], NULL);
    }

    create_work(s, "v", -1, one, one, 1, 4 * n, 1);

    {
        struct circuit *c;

        c = circ_parse(circuit);
        {
            int circnamelen;
            char *circname;
            circnamelen = strlen(s->dir) + strlen("/circuit") + 2;
            circname = (char *) malloc(sizeof(char) * circnamelen);
            (void) snprintf(circname, circnamelen, "%s/circuit", s->dir);
            (void) circ_copy_circuit(circuit, circname);
            free(circname);
        }
        (void) circ_evaluate(c, alphas, betas, c_star, s->mmap.x0);
        circ_cleanup(c);
    }

    // Last element to encode, so don't need to parallelize
    {
        mpz_t elems[2];
        int *pows;

        mpz_init_set_ui(elems[0], 0);
        mpz_init_set(elems[1], c_star);

        // The index set is laid out as follows:
        //   - The first 2 * n entries contain X_i,0 and X_i,1
        //   - The next n entries contain Z_i
        //   - The next n entries contain W_i
        //   - The final entry contains Y
        pows = (int *) calloc(s->mmap.nzs, sizeof(int));
        
        // The C* encoding contains everything but the W_i symbols.
        // Here we calculate the appropriate indices and powers.
        for (int i = 0; i < n; ++i) {
            // X_i,0^deg(x_i)
            pows[2 * i] = xdegs[i];
            // X_i,1^deg(x_i)
            pows[2 * i + 1] = xdegs[i];
            // Z_i
            pows[2 * n + i] = 1;
        }
        // Y^deg(y)
        pows[4 * n] = ydeg;
        // Encode against these indices/powers

        clt_encode(tmp, &s->mmap, 2, elems, pows, s->rand);

        mpz_clears(elems[0], elems[1], NULL);
        free(pows);
    }
    (void) write_element(s->dir, (mmap_enc *) &tmp, "c_star");

    mpz_clears(c_star, tmp, NULL);

    thpool_wait(s->thpool);
}

int
zobf_evaluate(const char *dir, const char *circuit, const char *input,
              int n, int m, uint64_t nthreads, uint64_t flags)
{
    char fname[100];
    int fnamelen = 100;
    int iszero;
    clt_pp pp;
    mpz_t c_1, c_2, z, w;
    mpz_t *xs, *xones, *ys, *yones;

    mpz_inits(c_1, c_2, z, w, NULL);

    {
        FILE *fp = open_file(dir, "params", "r+b");
        clt_pp_fread(fp, &pp);
        fclose(fp);

    }

    xs = (mpz_t *) malloc(sizeof(mpz_t) * n);
    xones = (mpz_t *) malloc(sizeof(mpz_t) * n);
    for (int i = 0; i < n; ++i) {
        mpz_inits(xs[i], xones[i], NULL);
    }
    ys = (mpz_t *) malloc(sizeof(mpz_t) * m);
    yones = (mpz_t *) malloc(sizeof(mpz_t) * m);
    for (int i = 0; i < m; ++i) {
        mpz_inits(ys[i], yones[i], NULL);
    }

    // Check that all input choices are bits
    for (int i = 0; i < n; ++i) {
        if (input[i] != '0' && input[i] != '1') {
            return -1;
        }
    }

    (void) omp_set_num_threads(nthreads);

    // Load in appropriate input
    for (int i = 0; i < n; ++i) {
        (void) snprintf(fname, fnamelen, "%s/x_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, xs[i]);
        (void) snprintf(fname, fnamelen, "%s/u_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, xones[i]);
    }

    // Load in secret input
    for (int i = 0; i < m; ++i) {
        (void) snprintf(fname, fnamelen, "%s/y_%d", dir, i);
        (void) load_mpz_scalar(fname, ys[i]);
        (void) snprintf(fname, fnamelen, "%s/v", dir);
        (void) load_mpz_scalar(fname, yones[i]);
    }

    // Evaluate the circuit on x_1, ..., x_n
    {
        struct circuit *c;

        c = circ_parse(circuit);
        circ_evaluate_encoding(c, xs, xones, ys, yones, c_1, pp.x0);
        circ_cleanup(c);
    }

    // Load in c_2
    (void) snprintf(fname, fnamelen, "%s/c_star", dir);
    (void) load_mpz_scalar(fname, c_2);

    // Compute c_1 * \Prod z_{i,x_i} and c_2 * \Prod w_{i,x_i}
    for (int i = 0; i < n; ++i) {
        (void) snprintf(fname, fnamelen, "%s/z_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, z);
        mpz_mul(c_1, c_1, z);
        mpz_mod(c_1, c_1, pp.x0);
        (void) snprintf(fname, fnamelen, "%s/w_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, w);
        mpz_mul(c_2, c_2, w);
        mpz_mod(c_2, c_2, pp.x0);
    }

    // Compute c_1 - c_2 and zero test
    {
        mpz_t tmp;
        mpz_inits(tmp, NULL);
        mpz_sub(tmp, c_1, c_2);
        iszero = clt_is_zero(&pp, tmp);
        mpz_clears(tmp, NULL);
    }

    for (int i = 0; i < m; ++i) {
        mpz_clears(ys[i], yones[i], NULL);
    }
    free(ys);
    free(yones);
    for (int i = 0; i < n; ++i) {
        mpz_clears(xs[i], xones[i], NULL);
    }
    free(xs);
    free(xones);
    mpz_clears(c_1, c_2, z, w, NULL);

    return iszero;
}
