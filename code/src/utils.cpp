#include "utils.h"
#include "mpz_pylong.h"

#include <fcntl.h>
#include <gmp.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

int g_verbose;

// XXX: The use of /dev/urandom is not secure; however, the supercomputer we run
// on doesn't appear to have enough entropy, and blocks for long periods of
// time.  Thus, we use /dev/urandom instead.
#ifndef RANDFILE
#  define RANDFILE "/dev/urandom"
#endif

void
state_destructor(PyObject *self)
{
    struct state *s;

    s = (struct state *) PyCapsule_GetPointer(self, NULL);
    if (s) {
        if (s->gs) {
            for (unsigned long i = 0; i < s->n; ++i) {
                mpz_clear(s->gs[i]);
            }
            free(s->gs);
        }
        if (s->crt_coeffs) {
            for (unsigned long i = 0; i < s->n; ++i) {
                mpz_clear(s->crt_coeffs[i]);
            }
            free(s->crt_coeffs);
        }
        if (s->zinvs) {
            for (unsigned long i = 0; i < s->nzs; ++i) {
                mpz_clear(s->zinvs[i]);
            }
            free(s->zinvs);
        }
        gmp_randclear(s->rng);
        mpz_clears(s->q, s->pzt, NULL);
    }
}

void *
pymalloc(const size_t size)
{
    void * r;
    if ((r = malloc(size)) == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory");
    }
    return r;
}

PyObject *
mpz_to_py(const mpz_t in)
{
    PyObject *outs, *out;
    char *buffer;

    buffer = mpz_get_str(NULL, 10, in);
    outs = PyString_FromString(buffer);
    out = PyNumber_Long(outs);
    free(buffer);
    return out;
}

void
py_to_mpz(mpz_t out, PyObject *in)
{
    (void) mpz_set_pylong(out, in);
}

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
    fname = (char *) pymalloc(sizeof(char) * fnamelen);
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
write_setup_params(struct state *s, long nu)
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
    mpz_set_ui(tmp, s->size);
    (void) snprintf(fname, len, "%s/size", s->dir);
    (void) save_mpz_scalar(fname, tmp);
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

PyObject *
obf_verbose(PyObject *self, PyObject *args)
{
    PyObject *py_verbose;

    if (!PyArg_ParseTuple(args, "O", &py_verbose))
        return NULL;

    g_verbose = PyObject_IsTrue(py_verbose);
    if (g_verbose == -1)
        return NULL;

    Py_RETURN_NONE;
}

PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long alpha, beta, eta, nu, kappa, rho_f;
    mpz_t *ps, *zs;
    PyObject *py_s, *py_gs;
    double start, end;
    struct state *s;

    s = (struct state *) pymalloc(sizeof(struct state));
    if (s == NULL)
        goto error;
    py_s = PyCapsule_New((void *) s, NULL, state_destructor);
    if (py_s == NULL)
        goto error;
    
    if (!PyArg_ParseTuple(args, "lllls", &s->secparam, &kappa, &s->size,
                          &s->nzs, &s->dir))
        goto error;

    /* Calculate CLT parameters */
    alpha = s->secparam;
    beta = s->secparam;
    s->rho = s->secparam;
    rho_f = kappa * (s->rho + alpha + 2);
    eta = rho_f + alpha + 2 * beta + s->secparam + 8;
    nu = eta - beta - rho_f - s->secparam + 3;
    s->n = (int) (eta * log2((float) s->secparam));

    if (g_verbose) {
        fprintf(stderr, "  Security Parameter: %ld\n", s->secparam);
        fprintf(stderr, "  Kappa: %ld\n", kappa);
        fprintf(stderr, "  Alpha: %ld\n", alpha);
        fprintf(stderr, "  Beta: %ld\n", beta);
        fprintf(stderr, "  Eta: %ld\n", eta);
        fprintf(stderr, "  Nu: %ld\n", nu);
        fprintf(stderr, "  Rho: %ld\n", s->rho);
        fprintf(stderr, "  Rho_f: %ld\n", rho_f);
        fprintf(stderr, "  N: %ld\n", s->n);
        fprintf(stderr, "  Number of Zs: %ld\n", s->nzs);
        fprintf(stderr, "  Size: %ld\n", s->size);
    }

    ps = (mpz_t *) pymalloc(sizeof(mpz_t) * s->n);
    s->gs = (mpz_t *) pymalloc(sizeof(mpz_t) * s->n);
    s->crt_coeffs = (mpz_t *) pymalloc(sizeof(mpz_t) * s->n);
    zs = (mpz_t *) pymalloc(sizeof(mpz_t) * s->nzs);
    s->zinvs = (mpz_t *) pymalloc(sizeof(mpz_t) * s->nzs);
    if (!ps || !s->gs || !s->crt_coeffs || !zs || !s->zinvs)
        goto error;

    py_gs = PyList_New(s->secparam);
    if (!py_gs)
        goto error;

    seed_rng(&s->rng);

    /* initialize gmp variables */
    mpz_init_set_ui(s->q, 1);
    mpz_init_set_ui(s->pzt, 0);
    for (unsigned long i = 0; i < s->n; ++i) {
        mpz_init_set_ui(ps[i], 1);
        mpz_inits(s->gs[i], s->crt_coeffs[i], NULL);
    }
    for (unsigned long i = 0; i < s->nzs; ++i) {
        mpz_inits(zs[i], s->zinvs[i], NULL);
    }

    /* Generate p_i's and g_i's, as well as q = \prod p_i */
    start = current_time();
#pragma omp parallel for
    for (unsigned long i = 0; i < s->n; ++i) {
        mpz_t p_unif;
        mpz_init(p_unif);
        // XXX: the primes generated here aren't officially uniform
        mpz_urandomb(p_unif, s->rng, eta);
        mpz_nextprime(ps[i], p_unif);
        mpz_urandomb(p_unif, s->rng, alpha);
        mpz_nextprime(s->gs[i], p_unif);
#pragma omp critical
        {
            //
            // This step is very expensive, and unfortunately it blocks the
            // parallelism of generating the primes.
            //
            mpz_mul(s->q, s->q, ps[i]);
        }
        mpz_clear(p_unif);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Generating p_i's and g_i's: %f\n",
                       end - start);

#ifdef ATTACK
    (void) gmp_fprintf(stderr, "g_1 in use: %Zd\n", s->gs[0]);
#endif

    /* Convert g_i values to python objects */
    start = current_time();
    //
    // Only convert the first secparam g_i values since we only need to fill in
    // the first secparam slots of the plaintext space.
    //
    for (unsigned long i = 0; i < s->secparam; ++i) {
        PyList_SetItem(py_gs, i, mpz_to_py(s->gs[i]));
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Converting g_i's to python objects: %f\n",
                       end - start);

    /* Compute CRT coefficients */
    start = current_time();
    //
    // This step is needed for making encoding efficient.  Unfortunately, the
    // CRT coefficients take an enormous amount of memory to compute / store.
    //
#pragma omp parallel for
    for (unsigned long i = 0; i < s->n; ++i) {
        mpz_t q;
        mpz_init(q);
        mpz_tdiv_q(q, s->q, ps[i]);
        mpz_invert(s->crt_coeffs[i], q, ps[i]);
        mpz_mul(s->crt_coeffs[i], s->crt_coeffs[i], q);
        mpz_clear(q);
    }

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Generating CRT coefficients: %f\n",
                       end - start);

    /* Compute z_i's */
    start = current_time();
#pragma omp parallel for
    for (unsigned long i = 0; i < s->nzs; ++i) {
        do {
            mpz_urandomm(zs[i], s->rng, s->q);
        } while (mpz_invert(s->zinvs[i], zs[i], s->q) == 0);
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
        for (unsigned long i = 0; i < s->nzs; ++i) {
            mpz_mul(zk, zk, zs[i]);
            mpz_mod(zk, zk, s->q);
        }
#pragma omp parallel for
        for (unsigned long i = 0; i < s->n; ++i) {
            mpz_t tmp, qpi, rnd;
            mpz_inits(tmp, qpi, rnd, NULL);
            // compute (((g_i)^{-1} mod p_i) * z^k mod p_i) * r_i * (q / p_i)
            mpz_invert(tmp, s->gs[i], ps[i]);
            mpz_mul(tmp, tmp, zk);
            mpz_mod(tmp, tmp, ps[i]);
            mpz_genrandom(rnd, &s->rng, beta);
            mpz_mul(tmp, tmp, rnd);
            mpz_div(qpi, s->q, ps[i]);
            mpz_mul(tmp, tmp, qpi);
#pragma omp critical
            {
                mpz_add(s->pzt, s->pzt, tmp);
            }
            mpz_clears(tmp, qpi, rnd, NULL);
        }
        mpz_mod(s->pzt, s->pzt, s->q);
        mpz_clear(zk);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Generating pzt: %f\n", end - start);

    (void) write_setup_params(s, nu);

    for (unsigned long i = 0; i < s->n; ++i) {
        mpz_clear(ps[i]);
    }
    free(ps);
    for (unsigned long i = 0; i < s->nzs; ++i) {
        mpz_clear(zs[i]);
    }
    free(zs);

    return PyTuple_Pack(2, py_s, py_gs);

error:
    if (ps)
        free(ps);
    if (zs)
        free(zs);
    if (s) {
        if (s->gs)
            free(s->gs);
        if (s->crt_coeffs)
            free(s->crt_coeffs);
        if (s->zinvs)
            free(s->zinvs);
    }
    return NULL;
}

PyObject *
obf_max_mem_usage(PyObject *self, PyObject *args)
{
    struct rusage usage;

    (void) getrusage(RUSAGE_SELF, &usage);
    (void) fprintf(stderr, "Max memory usage: %ld\n", usage.ru_maxrss);

    Py_RETURN_NONE;
}

