#include <Python.h>

#include <gmp.h>
#include <omp.h>
#include <sys/resource.h>

#ifdef ATTACK
#include <fplll.h>
#endif

#include "mpz_pylong.h"
#include "utils.h"

struct state {
    gmp_randstate_t rng;
    long secparam;
    long size;
    long n;
    long nzs;
    long rho;
    mpz_t q;
    mpz_t pzt;
    mpz_t *gs;
    mpz_t *crt_coeffs;
    mpz_t *zinvs;
    char *dir;
};

static void
state_destructor(PyObject *self)
{
    struct state *s;

    s = (struct state *) PyCapsule_GetPointer(self, NULL);
    if (s) {
        if (s->gs) {
            for (int i = 0; i < s->n; ++i) {
                mpz_clear(s->gs[i]);
            }
            free(s->gs);
        }
        if (s->crt_coeffs) {
            for (int i = 0; i < s->n; ++i) {
                mpz_clear(s->crt_coeffs[i]);
            }
            free(s->crt_coeffs);
        }
        if (s->zinvs) {
            for (int i = 0; i < s->nzs; ++i) {
                mpz_clear(s->zinvs[i]);
            }
            free(s->zinvs);
        }
        gmp_randclear(s->rng);
        mpz_clears(s->q, s->pzt, NULL);
    }
}

static void *
pymalloc(const size_t size)
{
    void * r;
    if ((r = malloc(size)) == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory");
    }
    return r;
}

static int
extract_indices(PyObject *py_list, int *idx1, int *idx2)
{
    *idx1 = -1;
    *idx2 = -1;
    switch (PyList_GET_SIZE(py_list)) {
    case 2:
        *idx2 = PyLong_AsLong(PyList_GET_ITEM(py_list, 1));
        /* fallthrough */
    case 1:
        *idx1 = PyLong_AsLong(PyList_GET_ITEM(py_list, 0));
        break;
    default:
        return FAILURE;
    }
    return SUCCESS;
}

static PyObject *
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

static void
py_to_mpz(mpz_t out, PyObject *in)
{
    (void) mpz_set_pylong(out, in);
}

static void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits)
{
    mpz_t one;
    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, *rng, nbits);
    mpz_clear(one);
}

static int
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

static int
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

static int
write_vector(struct state *s, mpz_t *vector, long size, char *name)
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
    (void) save_mpz_vector(fname, vector, size);
    free(fname);

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Saving to file: %f\n", end - start);
    return SUCCESS;
}

static int
write_layer(struct state *s, int inp, long idx, mpz_t *zero, mpz_t *one, long size)
{
    mpz_t z;
    char *fname;
    int fnamelen;
    double start, end;

    start = current_time();
    fnamelen = strlen(s->dir) + 10;
    fname = (char *) pymalloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return FAILURE;
    mpz_init_set_ui(z, inp);
    (void) snprintf(fname, fnamelen, "%s/%ld.input", s->dir, idx);
    (void) save_mpz_scalar(fname, z);
    (void) snprintf(fname, fnamelen, "%s/%ld.zero", s->dir, idx);
    (void) save_mpz_vector(fname, zero, size);
    (void) snprintf(fname, fnamelen, "%s/%ld.one", s->dir, idx);
    (void) save_mpz_vector(fname, one, size);
    free(fname);
    mpz_clear(z);
    end = current_time();

    if (g_verbose)
        (void) fprintf(stderr, "  Saving to file: %f\n", end - start);
    return SUCCESS;
}

static void
encode(struct state *st, mpz_t out, const PyObject *in, const long row,
       const long idx1, const long idx2)
{
    mpz_t r, tmp;

    mpz_inits(r, tmp, NULL);

    mpz_set_ui(out, 0);

    for (long i = 0; i < st->n; ++i) {
        mpz_genrandom(r, &st->rng, st->rho);
        mpz_mul(tmp, r, st->gs[i]);

        if (i < st->secparam) {
            py_to_mpz(r, PyList_GET_ITEM(PyList_GET_ITEM(in, i), row));
            mpz_add(tmp, tmp, r);
        }
        mpz_mul(tmp, tmp, st->crt_coeffs[i]);
        mpz_add(out, out, tmp);

    }
    mpz_mod(out, out, st->q);
    if (idx1 >= 0) {
        mpz_mul(out, out, st->zinvs[idx1]);
        mpz_mod(out, out, st->q);
    }
    if (idx2 >= 0) {
        mpz_mul(out, out, st->zinvs[idx2]);
        mpz_mod(out, out, st->q);
    }

    mpz_clears(r, tmp, NULL);
}

//
//
// Python functions
//
//

static PyObject *
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

static PyObject *
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
    for (int i = 0; i < s->n; ++i) {
        mpz_init_set_ui(ps[i], 1);
        mpz_inits(s->gs[i], s->crt_coeffs[i], NULL);
    }
    for (int i = 0; i < s->nzs; ++i) {
        mpz_inits(zs[i], s->zinvs[i], NULL);
    }

    /* Generate p_i's and g_i's, as well as q = \prod p_i */
    start = current_time();
#pragma omp parallel for
    for (int i = 0; i < s->n; ++i) {
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
    (void) gmp_fprintf(stderr, "g_1 = %Zd\n", s->gs[0]);
#endif

    /* Convert g_i values to python objects */
    start = current_time();
    //
    // Only convert the first secparam g_i values since we only need to fill in
    // the first secparam slots of the plaintext space.
    //
    for (int i = 0; i < s->secparam; ++i) {
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
    for (int i = 0; i < s->n; ++i) {
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
    for (int i = 0; i < s->nzs; ++i) {
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
        for (int i = 0; i < s->nzs; ++i) {
            mpz_mul(zk, zk, zs[i]);
            mpz_mod(zk, zk, s->q);
        }
#pragma omp parallel for
        for (int i = 0; i < s->n; ++i) {
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

    for (int i = 0; i < s->n; ++i) {
        mpz_clear(ps[i]);
    }
    free(ps);
    for (int i = 0; i < s->nzs; ++i) {
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

//
// Encode N scalars across all slots of the MLM
//
static PyObject *
obf_encode_scalars(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_scalars, *py_list;
    int idx1, idx2;
    char *name;
    mpz_t val;
    int err = 0;
    double start, end;
    struct state *st;

    if (!PyArg_ParseTuple(args, "OOOs", &py_state, &py_scalars, &py_list,
                          &name))
        return NULL;

    st = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (st == NULL)
        return NULL;
    
    (void) extract_indices(py_list, &idx1, &idx2);
    
    mpz_init(val);

    start = current_time();
    encode(st, val, py_scalars, 0, idx1, idx2);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding one element: %f\n", end - start);

    (void) write_scalar(st, val, name);
    mpz_clear(val);

    if (err) {
        return NULL;
    } else {
        Py_RETURN_NONE;
    }
}

//
// Encode N vectors across all slots of the MLM
//
static PyObject *
obf_encode_vectors(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_vectors, *py_list;
    char *name;
    int idx1, idx2;
    mpz_t *vector;
    Py_ssize_t length;
    double start, end;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OOOs", &py_state, &py_vectors, &py_list, &name))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;
    
    (void) extract_indices(py_list, &idx1, &idx2);

    // We assume that all vectors have the same length, and thus just grab the
    // length of the first vector
    length = PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0));
    vector = (mpz_t *) pymalloc(sizeof(mpz_t) * length);
    if (vector == NULL)
        return NULL;

    start = current_time();
#pragma omp parallel for
    for (Py_ssize_t i = 0; i < length; ++i) {
        mpz_init(vector[i]);
        encode(s, vector[i], py_vectors, i, idx1, idx2);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       length, end - start);

    (void) write_vector(s, vector, length, name);

    for (Py_ssize_t i = 0; i < length; ++i) {
        mpz_clear(vector[i]);
    }
    free(vector);

    Py_RETURN_NONE;
}

//
// Encode N layers across all slots of the MLM
//
static PyObject *
obf_encode_layers(PyObject *self, PyObject *args)
{
    PyObject *py_zero_ms, *py_one_ms;
    PyObject *py_zero_set, *py_one_set;
    PyObject *py_state;
    int zeroidx1, zeroidx2, oneidx1, oneidx2;
    int err = 0;
    long inp, idx;
    Py_ssize_t size;
    mpz_t *zero, *one;
    double start, end;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OllOOOO", &py_state, &idx, &inp, &py_zero_ms,
                          &py_one_ms, &py_zero_set, &py_one_set))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;
    
    (void) extract_indices(py_zero_set, &zeroidx1, &zeroidx2);
    (void) extract_indices(py_one_set, &oneidx1, &oneidx2);

    if (zeroidx1 < 0 && zeroidx2 < 0)
        return NULL;
    if (oneidx1 < 0 && oneidx2 < 0)
        return NULL;

    size = PyList_GET_SIZE(PyList_GET_ITEM(py_zero_ms, 0));
    zero = (mpz_t *) pymalloc(sizeof(mpz_t) * size);
    one = (mpz_t *) pymalloc(sizeof(mpz_t) * size);
    if (!zero || !one)
        return NULL;

    start = current_time();
#pragma omp parallel for
    for (Py_ssize_t ctr = 0; ctr < 2 * size; ++ctr) {
        PyObject *py_array;
        int idx1, idx2;
        mpz_t *val;
        size_t i;

        if (ctr < size) {
            i = ctr;
            val = &zero[i];
            py_array = py_zero_ms;
            idx1 = zeroidx1;
            idx2 = zeroidx2;
        } else {
            i = ctr - size;
            val = &one[i];
            py_array = py_one_ms;
            idx1 = oneidx1;
            idx2 = oneidx2;
        }

        mpz_init(*val);
        encode(s, *val, py_array, i, idx1, idx2);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       2 * size, end - start);

    (void) write_layer(s, inp, idx, zero, one, size);

    for (int i = 0; i < size; ++i) {
        mpz_clears(zero[i], one[i], NULL);
    }
    free(zero);
    free(one);

    if (err)
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
}

static PyObject *
obf_evaluate(PyObject *self, PyObject *args)
{
    char *dir = NULL;
    char *input = NULL;
    char *fname = NULL;
    int fnamelen;
    int iszero = -1;
    mpz_t *comp, *s, *t;
    mpz_t p, tmp;
    long bplen, size;
    int err = 0;
    int islayered;

    if (!PyArg_ParseTuple(args, "ssli", &dir, &input, &bplen, &islayered))
        return NULL;

    fnamelen = strlen(dir) + 20; // XXX: should include bplen somewhere

    fname = (char *) pymalloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(tmp, p, NULL);

    // Get the size of the matrices
    (void) snprintf(fname, fnamelen, "%s/size", dir);
    (void) load_mpz_scalar(fname, tmp);
    size = mpz_get_ui(tmp);

    comp = (mpz_t *) pymalloc(sizeof(mpz_t) * size * size);
    s = (mpz_t *) pymalloc(sizeof(mpz_t) * size);
    t = (mpz_t *) pymalloc(sizeof(mpz_t) * size);
    if (!comp || !s || !t) {
        err = 1;
        goto cleanup;
    }

    for (int i = 0; i < size; ++i) {
        mpz_inits(s[i], t[i], NULL);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_init(comp[i]);
    }

    if (!islayered) {
        (void) snprintf(fname, fnamelen, "%s/p_enc", dir);
        (void) load_mpz_scalar(fname, p);
    }

    for (int layer = 0; layer < bplen; ++layer) {
        unsigned int input_idx;
        double start, end;

        start = current_time();

        // find out the input bit for the given layer
        (void) snprintf(fname, fnamelen, "%s/%d.input", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        input_idx = mpz_get_ui(tmp);
        if (input_idx < 0 || input_idx >= strlen(input)) {
            PyErr_SetString(PyExc_RuntimeError, "invalid input");
            err = 1;
            break;
        }
        if (input[input_idx] != '0' && input[input_idx] != '1') {
            PyErr_SetString(PyExc_RuntimeError, "input must be 0 or 1");
            err = 1;
            break;
        }
        // load in appropriate matrix for the given input value
        if (input[input_idx] == '0') {
            (void) snprintf(fname, fnamelen, "%s/%d.zero", dir, layer);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", dir, layer);
        }

        if (layer == 0) {
            (void) load_mpz_vector(fname, comp, size * size);
        } else {
            mpz_t *tmparray;

            tmparray = (mpz_t *) malloc(sizeof(mpz_t) * size * size);
            for (int i = 0; i < size * size; ++i) {
                mpz_init(tmparray[i]);
            }
            (void) load_mpz_vector(fname, tmparray, size * size);
            mat_mult(comp, tmparray, size);
            for (int i = 0; i < size * size; ++i) {
                mpz_clear(tmparray[i]);
            }
            free(tmparray);
        }
        if (!islayered) {
            if (input[input_idx] == '0') {
                (void) snprintf(fname, fnamelen, "%s/%d.a0_enc", dir, layer);
            } else {
                (void) snprintf(fname, fnamelen, "%s/%d.a1_enc", dir, layer);
            }
            (void) load_mpz_scalar(fname, tmp);
            mpz_mul(p, p, tmp);
        }

        end = current_time();

        if (g_verbose)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n",
                           end - start);
    }

    if (!err) {
        double start, end;
        mpz_t pzt, q, nu;

        start = current_time();

        mpz_inits(pzt, q, nu, NULL);
        (void) snprintf(fname, fnamelen, "%s/s_enc", dir);
        (void) load_mpz_vector(fname, s, size);
        (void) snprintf(fname, fnamelen, "%s/t_enc", dir);
        (void) load_mpz_vector(fname, t, size);
        mat_mult_by_vects(tmp, s, comp, t, size);
        (void) snprintf(fname, fnamelen, "%s/pzt", dir);
        (void) load_mpz_scalar(fname, pzt);
        (void) snprintf(fname, fnamelen, "%s/q", dir);
        (void) load_mpz_scalar(fname, q);
        (void) snprintf(fname, fnamelen, "%s/nu", dir);
        (void) load_mpz_scalar(fname, nu);
        if (!islayered) {
            mpz_sub(tmp, tmp, p);
        }
        iszero = is_zero(tmp, pzt, q, mpz_get_ui(nu));

        mpz_clears(pzt, q, nu, NULL);

        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Zero test: %f\n", end - start);
    }

    for (int i = 0; i < size; ++i) {
        mpz_clears(s[i], t[i], NULL);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_clear(comp[i]);
    }

 cleanup:
    mpz_clears(tmp, p, NULL);

    if (comp)
        free(comp);
    if (s)
        free(s);
    if (t)
        free(t);

    if (fname)
        free(fname);

    if (err)
        return NULL;
    else
        return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyObject *
obf_max_mem_usage(PyObject *self, PyObject *args)
{
    struct rusage usage;

    (void) getrusage(RUSAGE_SELF, &usage);
    (void) fprintf(stderr, "Max memory usage: %ld\n", usage.ru_maxrss);

    Py_RETURN_NONE;
}

static PyObject *
obf_cleanup(PyObject *self, PyObject *args)
{
    PyObject *py_state;
    struct state *s;

    if (!PyArg_ParseTuple(args, "O", &py_state))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    gmp_randclear(s->rng);
    mpz_clears(s->q, s->pzt, NULL);
    for (int i = 0; i < s->n; ++i) {
        mpz_clears(s->gs[i], s->crt_coeffs[i], NULL);
    }
    free(s->gs);
    free(s->crt_coeffs);
    for (int i = 0; i < s->nzs; ++i) {
        mpz_clear(s->zinvs[i]);
    }
    free(s->zinvs);

    Py_RETURN_NONE;
}

#ifdef ATTACK
static PyObject *
obf_attack(PyObject *self, PyObject *args)
{
    mpz_t b, out, omega, Omega, q;
    double start, end;
    struct state *st;
    long kappa, bplen;
    int islayered, nslots;
    PyObject *py_out;

    st = (struct state *) pymalloc(sizeof(struct state));

    if (!PyArg_ParseTuple(args, "slilli", &st->dir, &bplen,
                          &islayered, &st->secparam, &kappa, &nslots))
        return NULL;

    if (!islayered) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Attack not implemented for Barrington's approach");
        return NULL;
    }

    mpz_inits(b, out, omega, Omega, q, NULL);

    /* Compute omega = pzt * u for some encoded element u */

    if (g_verbose)
        (void) fprintf(stderr, "Computing omega...\n");
    start = current_time();
    {
        char *fname;
        int fnamelen;
        mpz_t *comp;
        mpz_t tmp, pzt, nu;

        fnamelen = strlen(st->dir) + 20; // XXX: should include bplen somewhere
        fname = (char *) pymalloc(sizeof(char) * fnamelen);

        mpz_inits(tmp, pzt, nu, NULL);

        // Get the size of the matrices
        (void) snprintf(fname, fnamelen, "%s/size", st->dir);
        (void) load_mpz_scalar(fname, tmp);
        st->size = mpz_get_ui(tmp);

        comp = (mpz_t *) pymalloc(sizeof(mpz_t) * st->size * st->size);
        for (int i = 0; i < st->size * st->size; ++i) {
            mpz_init(comp[i]);
        }

        for (int layer = 0; layer < bplen; ++layer) {
            (void) snprintf(fname, fnamelen, "%s/%d.zero", st->dir, layer);
            (void) load_mpz_vector(fname, comp, st->size * st->size);
            if (layer == 0) {
                mpz_set(tmp, comp[0]);
            } else {
                mpz_mul(tmp, tmp, comp[0]);
            }
        }

        (void) snprintf(fname, fnamelen, "%s/s_enc", st->dir);
        (void) load_mpz_vector(fname, comp, st->size);
        mpz_mul(tmp, tmp, comp[0]);
        (void) snprintf(fname, fnamelen, "%s/t_enc", st->dir);
        (void) load_mpz_vector(fname, comp, st->size);
        mpz_mul(tmp, tmp, comp[st->size - 1]);

        (void) snprintf(fname, fnamelen, "%s/pzt", st->dir);
        (void) load_mpz_scalar(fname, pzt);
        (void) snprintf(fname, fnamelen, "%s/q", st->dir);
        (void) load_mpz_scalar(fname, q);
        // (void) snprintf(fname, fnamelen, "%s/nu", st->dir);
        // (void) load_mpz_scalar(fname, nu);

        // {
        //     int iszero = is_zero(tmp, pzt, q, mpz_get_ui(nu));
        //     fprintf(stderr, "iszero = %d\n", iszero);
        // }

        // Compute omega
        mpz_mul(tmp, tmp, pzt);
        mpz_mod_near(omega, tmp, q);

        mpz_clears(tmp, pzt, nu, NULL);

        for (int i = 0; i < st->size * st->size; ++i) {
            mpz_clear(comp[i]);
        }
        if (comp)
            free(comp);
        if (fname)
            free(fname);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Took: %f\n", end - start);

    /* We now have omega and q.  We now need to sample everything fresh to
       compute Omega. */

    {
        long alpha, beta, eta, nu, rho_f;
        mpz_t *hs, *ms, *ps, *rs;

        /* Calculate CLT parameters */
        alpha = st->secparam;
        beta = st->secparam;
        st->rho = st->secparam;
        rho_f = kappa * (st->rho + alpha + 2);
        eta = rho_f + alpha + 2 * beta + st->secparam + 8;
        nu = eta - beta - rho_f - st->secparam + 3;
        st->n = (int) (eta * log2((float) st->secparam));
        if (g_verbose) {
            (void) fprintf(stderr, "Parameters:\n");
            (void) fprintf(stderr, "  Security Parameter: %ld\n", st->secparam);
            (void) fprintf(stderr, "  Kappa: %ld\n", kappa);
            (void) fprintf(stderr, "  Alpha: %ld\n", alpha);
            (void) fprintf(stderr, "  Beta: %ld\n", beta);
            (void) fprintf(stderr, "  Eta: %ld\n", eta);
            (void) fprintf(stderr, "  Nu: %ld\n", nu);
            (void) fprintf(stderr, "  Rho: %ld\n", st->rho);
            (void) fprintf(stderr, "  Rho_f: %ld\n", rho_f);
            (void) fprintf(stderr, "  N: %ld\n", st->n);
            (void) fprintf(stderr, "  Size: %ld\n", st->size);
        }

        hs = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);
        ms = (mpz_t *) pymalloc(sizeof(mpz_t) * nslots);
        ps = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);
        rs = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);
        st->gs = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);

        seed_rng(&st->rng);

        /* initialize gmp variables */
        mpz_init_set_ui(st->pzt, 0);
        for (int i = 0; i < nslots; ++i) {
            mpz_inits(ms[i], NULL);
        }
        for (int i = 0; i < st->n; ++i) {
            mpz_inits(hs[i], ps[i], rs[i], st->gs[i], NULL);
        }

        /* Generate p_i's and g_i's */
        start = current_time();
#pragma omp parallel for
        for (int i = 0; i < st->n; ++i) {
            mpz_t p_unif;
            mpz_init(p_unif);
            mpz_urandomb(p_unif, st->rng, eta);
            mpz_nextprime(ps[i], p_unif);
            mpz_urandomb(p_unif, st->rng, alpha);
            mpz_nextprime(st->gs[i], p_unif);
            mpz_clear(p_unif);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating p_i's and g_i's: %f\n",
                           end - start);

        /* Compute h_i's */
        start = current_time();
#pragma omp parallel for
        for (int i = 0; i < st->n; ++i) {
            mpz_genrandom(hs[i], &st->rng, beta);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating h_i's: %f\n", end - start);

        /* Compute r_i's */
        start = current_time();
#pragma omp parallel for
        for (int i = 0; i < st->n; ++i) {
            mpz_genrandom(rs[i], &st->rng, rho_f);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating r_i's: %f\n", end - start);

        /* Sample m_i's */
        start = current_time();
#pragma omp parallel for
        for (int i = 0; i < nslots; ++i) {
            mpz_urandomb(ms[i], st->rng, alpha);
            mpz_mod(ms[i], ms[i], st->gs[i]);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating m_i's: %f\n", end - start);

        /* Compute b = \Omega \cdot g_1 */

        start = current_time();
        mpz_set_ui(b, 0L);
#pragma omp parallel for
        for (int i = 0; i < nslots; ++i) {
            mpz_t tmp, qpi;
            mpz_inits(tmp, qpi, NULL);
            mpz_mul(tmp, hs[i], ms[i]);
            mpz_div(qpi, q, ps[i]);
            mpz_mul(tmp, tmp, qpi);
            for (int j = 0; j < nslots; ++j) {
                if (j != i)
                    mpz_mul(tmp, tmp, st->gs[j]);
            }
#pragma omp critical
            {
                mpz_add(b, b, tmp);
            }
            mpz_clears(tmp, qpi, NULL);
        }
#pragma omp parallel for
        for (int i = 0; i < st->n; ++i) {
            mpz_t tmp, qpi;
            mpz_inits(tmp, qpi, NULL);
            mpz_mul(tmp, hs[i], rs[i]);
            mpz_div(qpi, q, ps[i]);
            mpz_mul(tmp, tmp, qpi);
            for (int j = 0; j < nslots; ++j) {
                mpz_mul(tmp, tmp, st->gs[j]);
            }
#pragma omp critical
            {
                mpz_add(b, b, tmp);
            }
            mpz_clears(tmp, qpi, NULL);
        }
        mpz_mod(b, b, q);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Computing b: %f\n", end - start);

        for (int i = 0; i < nslots; ++i) {
            mpz_clears(ms[i], NULL);
        }
        for (int i = 0; i < st->n; ++i) {
            mpz_clears(hs[i], ps[i], rs[i], NULL);
        }
        free(hs);
        free(ms);
        free(ps);
        free(rs);
    }

    /* Compute Omega = b / g_1 */
    mpz_set(Omega, b);
    for (int i = 0; i < nslots; ++i) {
        mpz_div(Omega, Omega, st->gs[i]);
    }

    {
        int r;
        mpz_t tmp;
        ZZ_mat<mpz_t> M(2, 2);

        // mpz_init(tmp);
        // mpz_set_ui(tmp, 4L);
        // mpz_mul_ui(tmp, tmp, G_1);
        // mpz_mul_ui(tmp, tmp, G_1);
        // mpz_mul_ui(tmp, tmp, G_2);
        // mpz_mul_ui(tmp, tmp, G_2);
        // mpz_div(tmp, q, tmp);
        // r = mpz_cmp(Omega, tmp);
        // mpz_clear(tmp);
        // fprintf(stderr, "Omega = %u | q / 4g_1^2g_2^2 = %u\n",
        //         mpz_sizeinbase(Omega, 2), mpz_sizeinbase(tmp, 2));
        // switch (r) {
        // case -1:
        //     (void) fprintf(stderr, "Omega < q / (4 g_1^2 g_2^2)\n");
        //     break;
        // case 1:
        //     fprintf(stderr, "Omega > q / (4 g_1^2 g_2^2)\n");
        //     break;
        // default:
        //     break;
        // }

        /* Apply LLL to lattice basis {(Omega, omega), (0, q)} */
        if (g_verbose)
            (void) fprintf(stderr, "Applying LLL...\n");
        start = current_time();
        M(0,0).set(Omega);
        M(0,1).set(omega);
        M(1,0) = 0L;
        M(1,1).set(q);
        fplll::lllReduction(M, 0.99, 0.51, LM_WRAPPER);
        /* Divide out Omega from output to get g_1 */
        mpz_div(out, M(0,0).getData(), Omega);
        mpz_abs(out, out);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Took: %f\n", end - start);
    }

    py_out = mpz_to_py(out);

    mpz_clears(b, out, omega, Omega, NULL);

    if (st) {
        if (st->gs)
            free(st->gs);
    }

    return py_out;
}
#endif


static PyMethodDef
ObfMethods[] = {
    {"verbose", obf_verbose, METH_VARARGS,
     "Set verbosity."},
    {"setup", obf_setup, METH_VARARGS,
     "Set up obfuscator."},
    {"encode_scalars", obf_encode_scalars, METH_VARARGS,
     "Encode a scalar in each slot."},
    {"encode_vectors", obf_encode_vectors, METH_VARARGS,
     "Encode a vector in each slot."},
    {"encode_layers", obf_encode_layers, METH_VARARGS,
     "Encode a branching program layer in each slot."},
    {"max_mem_usage", obf_max_mem_usage, METH_VARARGS,
     "Print out the maximum memory usage."},
    {"cleanup", obf_cleanup, METH_VARARGS,
     "Clean up objects created during setup."},
    {"evaluate", obf_evaluate, METH_VARARGS,
     "evaluate the obfuscation."},
#ifdef ATTACK
    {"attack", obf_attack, METH_VARARGS,
     "implementation of attack from paper (Section 3.1.2)."},
#endif
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
