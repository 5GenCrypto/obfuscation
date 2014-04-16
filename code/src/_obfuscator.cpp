#include <Python.h>

#include <assert.h>
#include <fcntl.h>
#include <gmp.h>
#include <omp.h>
#include <sys/time.h>

#include "mpz_pylong.h"

#define SUCCESS 1
#define FAILURE 0

// XXX: these are never cleaned up
static gmp_randstate_t g_rng;
static long g_n;
static long g_nzs;
static long g_rho;
static long g_nu;
static mpz_t g_x0;
static mpz_t g_pzt;
static mpz_t *g_gs;
static mpz_t *g_crt_coeffs;
static mpz_t *g_zinvs;
static char *g_dir;

#include "utils.cpp"

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long alpha, beta, eta, size;
    mpz_t *ps, *zs;
    PyObject *py_gs;

    if (!PyArg_ParseTuple(args, "lllllllls", &size, &g_n, &alpha,
                          &beta, &eta, &g_nu, &g_rho, &g_nzs, &g_dir))
        return NULL;

    // initialization

    ps = (mpz_t *) mymalloc(sizeof(mpz_t) * g_n);
    g_gs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_n);
    g_crt_coeffs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_n);
    zs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_nzs);
    g_zinvs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_nzs);
    if (!ps || !g_gs || !g_crt_coeffs || !zs || !g_zinvs)
        goto error;

    py_gs = PyList_New(g_n);
    if (!py_gs)
        goto error;

    mpz_init_set_ui(g_x0, 1);
    mpz_init_set_ui(g_pzt, 0);

    for (int i = 0; i < g_n; ++i) {
        mpz_init(ps[i]);
        mpz_init(g_gs[i]);
        mpz_init(g_crt_coeffs[i]);
    }
    for (int i = 0; i < g_nzs; ++i) {
        mpz_init(zs[i]);
        mpz_init(g_zinvs[i]);
    }

    // generate p_i's and g_i's, as well as compute x0

#pragma omp parallel for
    for (int i = 0; i < g_n; ++i) {
        mpz_t p_unif;
        mpz_init(p_unif);
        // XXX: not uniform primes
        mpz_urandomb(p_unif, g_rng, eta);
        mpz_nextprime(ps[i], p_unif);
        // XXX: not uniform primes
        mpz_urandomb(p_unif, g_rng, alpha);
        mpz_nextprime(g_gs[i], p_unif);
        PyList_SetItem(py_gs, i, mpz_to_py(g_gs[i]));
#pragma omp critical
        {
            mpz_mul(g_x0, g_x0, ps[i]);
        }
        mpz_clear(p_unif);
    }

    // generate CRT coefficients

#pragma omp parallel for
    for (int i = 0; i < g_n; ++i) {
        mpz_t q;
        mpz_init(q);
        mpz_tdiv_q(q, g_x0, ps[i]);
        mpz_invert(g_crt_coeffs[i], q, ps[i]);
        mpz_mul(g_crt_coeffs[i], g_crt_coeffs[i], q);
        mpz_clear(q);
    }

    // generate z_i's

#pragma omp parallel for
    for (int i = 0; i < g_nzs; ++i) {
        do {
            mpz_urandomm(zs[i], g_rng, g_x0);
        } while (mpz_invert(g_zinvs[i], zs[i], g_x0) == 0);
    }

    // generate pzt

    {
        mpz_t zk;
        mpz_init_set_ui(zk, 1);
        // compute z^k
        for (int i = 0; i < g_nzs; ++i) {
            mpz_mul(zk, zk, zs[i]);
            mpz_mod(zk, zk, g_x0);
        }
#pragma omp parallel for
        for (int i = 0; i < g_n; ++i) {
            mpz_t tmp, x0pi, rnd;
            mpz_init(tmp);
            mpz_init(x0pi);
            mpz_init(rnd);
            // compute (((g_i)^{-1} mod p_i) * z^k mod p_i) * r_i * (x_0 / p_i)
            mpz_invert(tmp, g_gs[i], ps[i]);
            mpz_mul(tmp, tmp, zk);
            mpz_mod(tmp, tmp, ps[i]);
            mpz_genrandom(rnd, beta);
            mpz_mul(tmp, tmp, rnd);
            mpz_div(x0pi, g_x0, ps[i]);
            mpz_mul(tmp, tmp, x0pi);
#pragma omp critical
            {
                mpz_add(g_pzt, g_pzt, tmp);
            }
            mpz_clear(tmp);
            mpz_clear(x0pi);
            mpz_clear(rnd);
        }
        mpz_mod(g_pzt, g_pzt, g_x0);
        mpz_clear(zk);
    }

    // save length, nu, x0, and pzt to file

    {
        char *fname;
        int len;
        mpz_t tmp;

        len = strlen(g_dir) + 10;

        fname = (char *) malloc(sizeof(char) * len);
        if (fname == NULL)
            return NULL;

        mpz_init(tmp);

        // save length
        mpz_set_ui(tmp, size);
        (void) snprintf(fname, len, "%s/size", g_dir);
        (void) save_mpz_scalar(fname, tmp);
        // save nu
        mpz_set_ui(tmp, g_nu);
        (void) snprintf(fname, len, "%s/nu", g_dir);
        (void) save_mpz_scalar(fname, tmp);
        // save x0
        (void) snprintf(fname, len, "%s/x0", g_dir);
        (void) save_mpz_scalar(fname, g_x0);
        // save pzt
        (void) snprintf(fname, len, "%s/pzt", g_dir);
        (void) save_mpz_scalar(fname, g_pzt);

        mpz_clear(tmp);

        free(fname);
    }

    // cleanup

    for (int i = 0; i < g_n; ++i) {
        mpz_clear(ps[i]);
    }
    free(ps);
    for (int i = 0; i < g_nzs; ++i) {
        mpz_clear(zs[i]);
    }
    free(zs);

    return py_gs;

error:
    if (ps)
        free(ps);
    if (g_gs)
        free(g_gs);
    if (g_crt_coeffs)
        free(g_crt_coeffs);
    if (zs)
        free(zs);
    if (g_zinvs)
        free(g_zinvs);
    return NULL;
}

static PyObject *
obf_encode_scalars(PyObject *self, PyObject *args)
{
    PyObject *py_scalars = NULL, *py_list = NULL;
    int idx1, idx2;
    char *name;
    mpz_t val;
    int err = 0;

    if (!PyArg_ParseTuple(args, "OOs", &py_scalars, &py_list, &name))
        return NULL;
    (void) extract_indices(py_list, &idx1, &idx2);

    mpz_init(val);

    encode(val, py_scalars, 0, idx1, idx2);

    {
        int fnamelen = strlen(g_dir) + strlen(name) + 2;
        char *fname = (char *) mymalloc(sizeof(char) * fnamelen);
        if (fname == NULL) {
            err = 1;
        } else {
            (void) snprintf(fname, fnamelen, "%s/%s", g_dir, name);
            (void) save_mpz_scalar(fname, val);
            free(fname);
        }
    }

    mpz_clear(val);

    if (err) {
        return NULL;
    } else {
        Py_RETURN_NONE;
    }
}

static PyObject *
obf_encode_vectors(PyObject *self, PyObject *args)
{
    PyObject *py_vectors = NULL, *py_list = NULL;
    int idx1, idx2, err = 0;
    mpz_t *vector;
    Py_ssize_t size;
    char *name;

    if (!PyArg_ParseTuple(args, "OOs", &py_vectors, &py_list, &name))
        return NULL;
    (void) extract_indices(py_list, &idx1, &idx2);

    size = PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0));
    vector = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    if (vector == NULL)
        return NULL;

#pragma omp parallel for
    for (Py_ssize_t i = 0; i < size; ++i) {
        mpz_init(vector[i]);
        (void) encode(vector[i], py_vectors, i, idx1, idx2);
    }

    if (!err) {
        int fnamelen = strlen(g_dir) + strlen(name) + 2;
        char *fname = (char *) mymalloc(sizeof(char) * fnamelen);
        if (fname == NULL) {
            err = 1;
        } else {
            (void) snprintf(fname, fnamelen, "%s/%s", g_dir, name);
            (void) save_mpz_vector(fname, vector, size);
            free(fname);
        }
    }

    for (Py_ssize_t i = 0; i < size; ++i) {
        mpz_clear(vector[i]);
    }
    free(vector);

    if (err)
        return NULL;
    else
        Py_RETURN_NONE;
}

static PyObject *
obf_encode_layers(PyObject *self, PyObject *args)
{
    PyObject *py_zero_ms, *py_one_ms;
    PyObject *py_zero_set, *py_one_set;
    int zeroidx1, zeroidx2, oneidx1, oneidx2;
    char *fname = NULL;
    int fnamelen;
    int err = 0;
    long inp, idx;
    Py_ssize_t size;
    mpz_t *zero;
    mpz_t *one;
    mpz_t z_inp;

    if (!PyArg_ParseTuple(args, "llOOOO", &idx, &inp, &py_zero_ms, &py_one_ms,
                          &py_zero_set, &py_one_set))
        return NULL;
    (void) extract_indices(py_zero_set, &zeroidx1, &zeroidx2);
    (void) extract_indices(py_one_set, &oneidx1, &oneidx2);

    // size = PyList_GET_SIZE(py_zero_m);
    size = PyList_GET_SIZE(PyList_GET_ITEM(py_zero_ms, 0));
    zero = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    one = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    if (zero == NULL || one == NULL)
        return NULL;

    mpz_init_set_ui(z_inp, inp);

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
        (void) encode(*val, py_array, i, idx1, idx2);
    }

    if (!err) {
        fnamelen = strlen(g_dir) + 10; // XXX: needs to include length of idx!
        fname = (char *) mymalloc(sizeof(char) * fnamelen);
        if (fname == NULL) {
            err = 1;
        } else {
            (void) snprintf(fname, fnamelen, "%s/%ld.input", g_dir, idx);
            (void) save_mpz_scalar(fname, z_inp);
            (void) snprintf(fname, fnamelen, "%s/%ld.zero", g_dir, idx);
            (void) save_mpz_vector(fname, zero, size);
            (void) snprintf(fname, fnamelen, "%s/%ld.one", g_dir, idx);
            (void) save_mpz_vector(fname, one, size);
        }
    }

    mpz_clear(z_inp);
    for (int i = 0; i < size; ++i) {
        mpz_clear(zero[i]);
        mpz_clear(one[i]);
    }
    free(zero);
    free(one);

    if (fname)
        free(fname);

    if (err)
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
}

static PyObject *
obf_evaluate(PyObject *self, PyObject *args)
{
    char *input = NULL;
    char *fname = NULL;
    int fnamelen;
    int iszero = -1;
    mpz_t *comp;
    mpz_t *tmp1;
    mpz_t *tmp2;
    mpz_t *s;
    mpz_t *t;
    mpz_t p1, p2, tmp;
    long bplen, size;
    int err = 0;
    int islayered;

    if (!PyArg_ParseTuple(args, "ssli", &g_dir, &input, &bplen, &islayered))
        return NULL;

    fnamelen = strlen(g_dir) + 10; // XXX: should include bplen somewhere

    fname = (char *) mymalloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_init(tmp);

    (void) snprintf(fname, fnamelen, "%s/size", g_dir);
    (void) load_mpz_scalar(fname, tmp);
    size = mpz_get_ui(tmp);

    comp = (mpz_t *) mymalloc(sizeof(mpz_t) * size * size);
    tmp1 = (mpz_t *) mymalloc(sizeof(mpz_t) * size * size);
    tmp2 = (mpz_t *) mymalloc(sizeof(mpz_t) * size * size);
    // XXX: memory leak

    s = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    t = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    // XXX: memory leak
    
    mpz_init(p1);
    mpz_init(p2);
    for (int i = 0; i < size; ++i) {
        mpz_init(s[i]);
        mpz_init(t[i]);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_init(comp[i]);
        mpz_init(tmp1[i]);
        mpz_init(tmp2[i]);
    }

    if (!islayered) {
        (void) snprintf(fname, fnamelen, "%s/p_enc", g_dir);
        (void) load_mpz_scalar(fname, p2);
    }

    for (int layer = 0; layer < bplen; ++layer) {
        unsigned int input_idx;

        // find out the input bit for the given layer
        (void) snprintf(fname, fnamelen, "%s/%d.input", g_dir, layer);
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
            (void) snprintf(fname, fnamelen, "%s/%d.zero", g_dir, layer);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", g_dir, layer);
        }

        if (layer == 0) {
            (void) load_mpz_vector(fname, comp, size * size);
        } else {
            (void) load_mpz_vector(fname, tmp1, size * size);
            mat_mult(tmp2, comp, tmp1, size);
            for (int ctr = 0; ctr < size * size; ++ctr) {
                mpz_set(comp[ctr], tmp2[ctr]);
            }
        }
        if (!islayered) {
            if (input[input_idx] == '0') {
                (void) snprintf(fname, fnamelen, "%s/%d.a0_enc", g_dir, layer);
            } else {
                (void) snprintf(fname, fnamelen, "%s/%d.a1_enc", g_dir, layer);
            }
            load_mpz_scalar(fname, tmp);
            mpz_mul(p2, p2, tmp);
        }
    }

    if (!err) {
        (void) snprintf(fname, fnamelen, "%s/s_enc", g_dir);
        (void) load_mpz_vector(fname, s, size);
        (void) snprintf(fname, fnamelen, "%s/t_enc", g_dir);
        (void) load_mpz_vector(fname, t, size);
        mat_mult_by_vects(p1, s, comp, t, size);
        if (islayered) {
            iszero = is_zero(p1);
        } else {
            mpz_sub(tmp, p1, p2);
            iszero = is_zero(tmp);
        }
    }

    mpz_clear(tmp);
    mpz_clear(p1);
    mpz_clear(p2);
    for (int i = 0; i < size; ++i) {
        mpz_clear(s[i]);
        mpz_clear(t[i]);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_clear(comp[i]);
        mpz_clear(tmp1[i]);
        mpz_clear(tmp2[i]);
    }

    if (comp)
        free(comp);
    if (tmp1)
        free(tmp1);
    if (tmp2)
        free(tmp2);
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

static PyMethodDef
ObfMethods[] = {
    {"setup", obf_setup, METH_VARARGS,
     "Set up obfuscator."},
    {"encode_scalars", obf_encode_scalars, METH_VARARGS,
     "Encode a scalar in each slot."},
    {"encode_vectors", obf_encode_vectors, METH_VARARGS,
     "Encode a vector in each slot."},
    {"encode_layers", obf_encode_layers, METH_VARARGS,
     "Encode a branching program layer in each slot."},
    {"evaluate", obf_evaluate, METH_VARARGS,
     "evaluate the obfuscation."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    int file;
    unsigned long seed;

    (void) Py_InitModule("_obfuscator", ObfMethods);

    if ((file = open("/dev/random", O_RDONLY)) == -1) {
        (void) fprintf(stderr, "Error opening /dev/random\n");
        goto cleanup;
    }
    if (read(file, &seed, sizeof seed) == -1) {
        (void) fprintf(stderr, "Error reading from /dev/random\n");
        goto cleanup;
    }

    fprintf(stderr, "SEED = %lu\n", seed);

    gmp_randinit_default(g_rng);
    gmp_randseed_ui(g_rng, seed);

 cleanup:
    if (file != -1)
        (void) close(file);
}
