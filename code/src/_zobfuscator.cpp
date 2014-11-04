#include "utils.h"
#include "pyutils.h"

struct zstate {
    struct state s;
    mpz_t nchk;
    mpz_t nev;
};

static void
zstate_destructor(PyObject *self)
{
    struct zstate *s;

    s = (struct zstate *) PyCapsule_GetPointer(self, NULL);
    if (s) {
        if (s->s.gs) {
            for (unsigned long i = 0; i < s->s.n; ++i) {
                mpz_clear(s->s.gs[i]);
            }
            free(s->s.gs);
        }
        if (s->s.crt_coeffs) {
            for (unsigned long i = 0; i < s->s.n; ++i) {
                mpz_clear(s->s.crt_coeffs[i]);
            }
            free(s->s.crt_coeffs);
        }
        if (s->s.zinvs) {
            for (unsigned long i = 0; i < s->s.nzs; ++i) {
                mpz_clear(s->s.zinvs[i]);
            }
            free(s->s.zinvs);
        }
        gmp_randclear(s->s.rng);
        mpz_clears(s->s.q, s->s.pzt, s->nev, s->nchk, NULL);
    }
}

static int
write_element(const struct zstate *s, mpz_t elem, const char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(s->s.dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", s->s.dir, name);
    (void) save_mpz_scalar(fname, elem);
    free(fname);
    return 0;
}

static void
encode(struct zstate *s, mpz_t out, mpz_t in1, mpz_t in2, unsigned int num, ...)
{
    va_list indices;
    mpz_t r, tmp;

    va_start(indices, num);
    mpz_inits(r, tmp, NULL);
    mpz_set_ui(out, 0);
    for (unsigned long i = 0; i < s->s.n; ++i) {
        mpz_genrandom(r, &s->s.rng, s->s.rho);
        mpz_mul(tmp, r, s->s.gs[i]);
        if (i < s->s.n / 2) {
            mpz_add(tmp, tmp, in1);
        } else {
            mpz_add(tmp, tmp, in2);
        }
        mpz_mul(tmp, tmp, s->s.crt_coeffs[i]);
        mpz_add(out, out, tmp);
    }
    mpz_mod(out, out, s->s.q);
    for (unsigned long i = 0; i < num; ++i) {
        mpz_powm_ui(tmp, s->s.zinvs[va_arg(indices, int)],
                    va_arg(indices, int), s->s.q);
        mpz_mul(out, out, tmp);
        mpz_mod(out, out, s->s.q);
    }
    mpz_clears(r, tmp, NULL);
}

PyObject *
zobf_setup(PyObject *self, PyObject *args)
{
    long alpha, beta, eta, nu, kappa, rho_f;
    mpz_t *ps, *zs;
    PyObject *py_s;
    double start, end;
    struct zstate *s;

    s = (struct zstate *) malloc(sizeof(struct zstate));
    if (s == NULL)
        return NULL;
    if (!PyArg_ParseTuple(args, "llls", &s->s.secparam, &kappa, &s->s.nzs,
                          &s->s.dir))
        return NULL;
    py_s = PyCapsule_New((void *) s, NULL, zstate_destructor);
    if (py_s == NULL)
        return NULL;

    /* Calculate CLT parameters */
    alpha = s->s.secparam;
    beta = s->s.secparam;
    s->s.rho = s->s.secparam;
    rho_f = kappa * (s->s.rho + alpha + 2);
    eta = rho_f + alpha + 2 * beta + s->s.secparam + 8;
    nu = eta - beta - rho_f - s->s.secparam + 3;
    s->s.n = (int) (eta * log2((float) s->s.secparam));

    if (g_verbose) {
        fprintf(stderr, "  Security Parameter: %ld\n", s->s.secparam);
        fprintf(stderr, "  Kappa: %ld\n", kappa);
        fprintf(stderr, "  Alpha: %ld\n", alpha);
        fprintf(stderr, "  Beta: %ld\n", beta);
        fprintf(stderr, "  Eta: %ld\n", eta);
        fprintf(stderr, "  Nu: %ld\n", nu);
        fprintf(stderr, "  Rho: %ld\n", s->s.rho);
        fprintf(stderr, "  Rho_f: %ld\n", rho_f);
        fprintf(stderr, "  N: %ld\n", s->s.n);
        fprintf(stderr, "  Number of Zs: %ld\n", s->s.nzs);
    }

    ps = (mpz_t *) malloc(sizeof(mpz_t) * s->s.n);
    s->s.gs = (mpz_t *) malloc(sizeof(mpz_t) * s->s.n);
    s->s.crt_coeffs = (mpz_t *) malloc(sizeof(mpz_t) * s->s.n);
    zs = (mpz_t *) malloc(sizeof(mpz_t) * s->s.nzs);
    s->s.zinvs = (mpz_t *) malloc(sizeof(mpz_t) * s->s.nzs);

    seed_rng(&s->s.rng);

    /* initialize gmp variables */
    mpz_init_set_ui(s->s.q, 1);
    mpz_init_set_ui(s->s.pzt, 0);
    for (unsigned long i = 0; i < s->s.n; ++i) {
        mpz_init_set_ui(ps[i], 1);
        mpz_inits(s->s.gs[i], s->s.crt_coeffs[i], NULL);
    }
    for (unsigned long i = 0; i < s->s.nzs; ++i) {
        mpz_inits(zs[i], s->s.zinvs[i], NULL);
    }

    /* Generate p_i's and g_i's, as well as q = \prod p_i */
    start = current_time();
    mpz_init_set_ui(s->nev, 1);
    mpz_init_set_ui(s->nchk, 1);
#pragma omp parallel for
    for (unsigned long i = 0; i < s->s.n; ++i) {
        mpz_t p_unif;
        mpz_init(p_unif);
        // XXX: the primes generated here aren't officially uniform
        mpz_urandomb(p_unif, s->s.rng, eta);
        mpz_nextprime(ps[i], p_unif);
        mpz_urandomb(p_unif, s->s.rng, alpha);
        mpz_nextprime(s->s.gs[i], p_unif);
#pragma omp critical
        {
            //
            // This step is very expensive, and unfortunately it blocks the
            // parallelism of generating the primes.
            //
            if (i < s->s.n / 2)
                mpz_mul(s->nev, s->nev, ps[i]);
            else
                mpz_mul(s->nchk, s->nchk, ps[i]);
            mpz_mul(s->s.q, s->s.q, ps[i]);
        }
        mpz_clear(p_unif);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Generating p_i's and g_i's: %f\n",
                       end - start);

    /* Compute CRT coefficients */
    start = current_time();
    //
    // This step is needed for making encoding efficient.  Unfortunately, the
    // CRT coefficients take an enormous amount of memory to compute / store.
    //
#pragma omp parallel for
    for (unsigned long i = 0; i < s->s.n; ++i) {
        mpz_t q;
        mpz_init(q);
        mpz_tdiv_q(q, s->s.q, ps[i]);
        mpz_invert(s->s.crt_coeffs[i], q, ps[i]);
        mpz_mul(s->s.crt_coeffs[i], s->s.crt_coeffs[i], q);
        mpz_clear(q);
    }

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Generating CRT coefficients: %f\n",
                       end - start);

    /* Compute z_i's */
    start = current_time();
#pragma omp parallel for
    for (unsigned long i = 0; i < s->s.nzs; ++i) {
        do {
            mpz_urandomm(zs[i], s->s.rng, s->s.q);
        } while (mpz_invert(s->s.zinvs[i], zs[i], s->s.q) == 0);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Generating z_i's: %f\n", end - start);

    /* Compute pzt */
    start = current_time();
    {
        mpz_t zk;
        mpz_init_set_ui(zk, 1);
        // compute z^k mod q
        for (unsigned long i = 0; i < s->s.nzs; ++i) {
            mpz_mul(zk, zk, zs[i]);
            mpz_mod(zk, zk, s->s.q);
        }
#pragma omp parallel for
        for (unsigned long i = 0; i < s->s.n; ++i) {
            mpz_t tmp, qpi, rnd;
            mpz_inits(tmp, qpi, rnd, NULL);
            // compute (((g_i)^{-1} mod p_i) * z^k mod p_i) * r_i * (q / p_i)
            mpz_invert(tmp, s->s.gs[i], ps[i]);
            mpz_mul(tmp, tmp, zk);
            mpz_mod(tmp, tmp, ps[i]);
            mpz_genrandom(rnd, &s->s.rng, beta);
            mpz_mul(tmp, tmp, rnd);
            mpz_div(qpi, s->s.q, ps[i]);
            mpz_mul(tmp, tmp, qpi);
#pragma omp critical
            {
                mpz_add(s->s.pzt, s->s.pzt, tmp);
            }
            mpz_clears(tmp, qpi, rnd, NULL);
        }
        mpz_mod(s->s.pzt, s->s.pzt, s->s.q);
        mpz_clear(zk);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Generating pzt: %f\n", end - start);

    (void) write_setup_params(&s->s, nu, 0);

    for (unsigned long i = 0; i < s->s.n; ++i) {
        mpz_clear(ps[i]);
    }
    free(ps);
    for (unsigned long i = 0; i < s->s.nzs; ++i) {
        mpz_clear(zs[i]);
    }
    free(zs);

    return py_s;
}

static PyObject *
zobf_encode_circuit(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_y;
    mpz_t zero, one, tmp;
    int n, m;
    char *fname;
    int fnamelen = 10;
    struct zstate *s;
    double start, end;

    if (!PyArg_ParseTuple(args, "OOii", &py_state, &py_y, &n, &m))
        return NULL;
    s = (struct zstate *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    fname = (char *) malloc(sizeof(char) * fnamelen);

    mpz_init(tmp);
    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(one, 1);

    for (int i = 0; i < n; ++i) {
        mpz_t alpha, delta, gamma, out;
        mpz_inits(alpha, delta, gamma, out, NULL);
        mpz_urandomm(alpha, s->s.rng, s->nchk);

        encode(s, out, zero, alpha, 1, 2 * i, 1);
        (void) snprintf(fname, fnamelen, "x_%d_0", i);
        (void) write_element(s, out, fname);

        encode(s, out, one, one, 1, 2 * i, 1);
        (void) snprintf(fname, fnamelen, "u_%d_0", i);
        (void) write_element(s, out, fname);

        encode(s, out, one, alpha, 1, 2 * i + 1, 1);
        (void) snprintf(fname, fnamelen, "x_%d_1", i);
        (void) write_element(s, out, fname);

        encode(s, out, one, one, 1, 2 * i + 1, 1);
        (void) snprintf(fname, fnamelen, "u_%d_1", i);
        (void) write_element(s, out, fname);

        mpz_urandomm(delta, s->s.rng, s->nev);
        mpz_urandomm(gamma, s->s.rng, s->nchk);

        encode(s, out, delta, gamma, 3, 2 * n + 2 * i, 1, 5 * n + i, 1, 6 * n + i, 1);
        (void) snprintf(fname, fnamelen, "z_%d_0", i);
        (void) write_element(s, out, fname);

        encode(s, out, zero, gamma, 1, 6 * n + i);
        (void) snprintf(fname, fnamelen, "w_%d_0", i);
        (void) write_element(s, out, fname);

        mpz_urandomm(delta, s->s.rng, s->nev);
        mpz_urandomm(gamma, s->s.rng, s->nchk);

        encode(s, out, delta, gamma, 3, 2 * n + 2 * i, 1, 5 * n + i, 1, 6 * n + i, 1);
        (void) snprintf(fname, fnamelen, "z_%d_1", i);
        (void) write_element(s, out, fname);

        encode(s, out, zero, gamma, 1, 6 * n + i);
        (void) snprintf(fname, fnamelen, "w_%d_0", i);
        (void) write_element(s, out, fname);

        mpz_clears(alpha, delta, gamma, out, NULL);
    }
    for (int i = 0; i < m; ++i) {
        mpz_t beta, out, y;
        mpz_inits(beta, out, y, NULL);
        mpz_urandomm(beta, s->s.rng, s->nchk);
        py_to_mpz(y, PyList_GET_ITEM(py_y, i));
        encode(s, out, y, beta, 1, 7 * n, 1);
        (void) snprintf(fname, fnamelen, "y_%d", i);
        (void) write_element(s, out, fname);
        mpz_clears(beta, out, y, NULL);
    }
    encode(s, tmp, one, one, 1, 7 * n);
    (void) write_element(s, tmp, "v");

    // TODO: compute and encode C*.

    Py_RETURN_NONE;
}

static PyMethodDef
ObfMethods[] = {
    {"verbose", obf_verbose, METH_VARARGS,
     "Set verbosity."},
    {"setup", zobf_setup, METH_VARARGS,
     "Set up obfuscator."},
    {"encode_circuit", zobf_encode_circuit, METH_VARARGS,
     "Encode [a, b] under given index sets."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_zobfuscator(void)
{
    (void) Py_InitModule("_zobfuscator", ObfMethods);
}
