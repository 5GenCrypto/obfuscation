#include <Python.h>

#include <assert.h>
#include <fcntl.h>
#include <gmp.h>
#include <omp.h>
#include <sys/time.h>

#include "mpz_pylong.h"

#define SUCCESS 1
#define FAILURE 0

#define GROUPLENGTH_SQ GROUPLENGTH * GROUPLENGTH

// XXX: these are never cleaned up
static gmp_randstate_t g_rng;
static long g_n;
static mpz_t g_x0;
static mpz_t *g_ps;
static mpz_t *g_gs;
static mpz_t g_pzt;
static mpz_t *g_crt_coeffs;
static mpz_t *g_zinvs;
static long g_nzs;
static long g_rho;
static long g_nu;
static char *g_dir;

#include "utils.cpp"

static PyObject *
obf_genprime(PyObject *self, PyObject *args)
{
    long bitlength;
    PyObject *py_out;
    mpz_t out;

    if (!PyArg_ParseTuple(args, "l", &bitlength))
        return NULL;

    mpz_init(out);
    mpz_genrandom(out, bitlength);
    mpz_nextprime(out, out);
    py_out = mpz_to_py(out);
    mpz_clear(out);
    return py_out;
}

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long alpha, beta, eta;
    PyObject *py_gs;
    Py_ssize_t gs_len;
    mpz_t *zs;

    if (!PyArg_ParseTuple(args, "lllllllOs", &g_n, &alpha, &beta, &eta,
                          &g_nu, &g_rho, &g_nzs, &py_gs, &g_dir))
        return NULL;
    if (check_pylist(py_gs) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "eighth argument must be a list of longs");
        return NULL;
    }

    // initialization

    g_ps = (mpz_t *) mymalloc(sizeof(mpz_t) * g_n);
    if (g_ps == NULL)
        goto error;
    g_gs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_n);
    if (g_ps == NULL)
        goto error;
    g_crt_coeffs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_n);
    if (g_crt_coeffs == NULL)
        goto error;

    zs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_nzs);
    if (zs == NULL)
        goto error;
    g_zinvs = (mpz_t *) mymalloc(sizeof(mpz_t) * g_nzs);
    if (g_zinvs == NULL)
        goto error;

    mpz_init_set_ui(g_x0, 1);
    mpz_init_set_ui(g_pzt, 0);

    for (int i = 0; i < g_n; ++i) {
        mpz_init(g_ps[i]);
        mpz_init(g_gs[i]);
        mpz_init(g_crt_coeffs[i]);
    }

    for (int i = 0; i < g_nzs; ++i) {
        mpz_init(zs[i]);
        mpz_init(g_zinvs[i]);
    }

    gs_len = PyList_GET_SIZE(py_gs);

    // generate p_i's and g_'s, as well as compute x0

#pragma omp parallel for
    for (int i = 0; i < g_n; ++i) {
        mpz_t p_unif;
        mpz_init(p_unif);
        // XXX: not uniform primes
        mpz_urandomb(p_unif, g_rng, eta);
        mpz_nextprime(g_ps[i], p_unif);
        if (i < gs_len) {
            py_to_mpz(g_gs[i], PyList_GET_ITEM(py_gs, i));
        } else {
            // XXX: not uniform primes
            mpz_urandomb(p_unif, g_rng, alpha);
            mpz_nextprime(g_gs[i], p_unif);
        }
#pragma omp critical
        {
            mpz_mul(g_x0, g_x0, g_ps[i]);
        }
        mpz_clear(p_unif);
    }

    // generate CRT coefficients

#pragma omp parallel for
    for (int i = 0; i < g_n; ++i) {
        mpz_t q;
        mpz_init(q);
        mpz_tdiv_q(q, g_x0, g_ps[i]);
        mpz_invert(g_crt_coeffs[i], q, g_ps[i]);
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
            mpz_invert(tmp, g_gs[i], g_ps[i]);
            mpz_mul(tmp, tmp, zk);
            mpz_mod(tmp, tmp, g_ps[i]);
            mpz_genrandom(rnd, beta);
            mpz_mul(tmp, tmp, rnd);
            mpz_div(x0pi, g_x0, g_ps[i]);
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

    // cleanup

    for (int i = 0; i < g_nzs; ++i) {
        mpz_clear(zs[i]);
    }
    free(zs);

    Py_RETURN_NONE;

 error:
    if (g_ps)
        free(g_ps);
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
obf_encode_scalar(PyObject *self, PyObject *args)
{
    PyObject *py_val = NULL, *py_list = NULL;
    int idx1, idx2;
    char *name;
    mpz_t val;
    int err = 0;

    if (!PyArg_ParseTuple(args, "OOs", &py_val, &py_list, &name)) {
        return NULL;
    }
    if (!PyLong_Check(py_val)) {
        PyErr_SetString(PyExc_RuntimeError, "first argument must be a long");
        return NULL;
    }
    if (!PyList_Check(py_list)) {
        PyErr_SetString(PyExc_RuntimeError, "third argument must be a list");
        return NULL;
    }

    if (extract_indices(py_list, &idx1, &idx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "third argument must be a list of length 1 or 2");
        return NULL;
    }

    mpz_init(val);

    py_to_mpz(val, py_val);
    if (encode(val, val, idx1, idx2, 0) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError, "encoding failed");
        err = 1;
        goto cleanup;
    }

    // {   // XXX: TEMPORARY
    //     mpz_t tmp;
    //     mpz_init(tmp);
    //     for (int i = 0; i < hack; ++i) {
    //         mpz_set_ui(tmp, 1);
    //         encode(tmp, tmp, i, -1, 0);
    //         mpz_mul(val, val, tmp);
    //     }
    //     mpz_clear(tmp);
    // }

    {
        int fnamelen = strlen(g_dir) + strlen(name) + 2;
        char *fname = (char *) malloc(sizeof(char) * fnamelen);
        if (fname == NULL) {
            err = 1;
            goto cleanup;
        }
        (void) snprintf(fname, fnamelen, "%s/%s", g_dir, name);
        (void) save_mpz_scalar(fname, val);
        free(fname);
    }

 cleanup:
    mpz_clear(val);

    if (err) {
        return NULL;
    } else {
        Py_RETURN_NONE;
    }
}

static PyObject *
obf_encode_vector(PyObject *self, PyObject *args)
{
    PyObject *py_vals = NULL, *py_list = NULL;
    int idx1, idx2, err = 0;
    mpz_t vector[GROUPLENGTH];
    char *name;

    if (!PyArg_ParseTuple(args, "OOs", &py_vals, &py_list, &name)) {
        return NULL;
    }
    if (check_pylist(py_vals) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "first argument must be a list of longs");
        return NULL;
    }
    if (!PyList_Check(py_list)) {
        PyErr_SetString(PyExc_RuntimeError, "third argument must be a list");
        return NULL;
    }

    if (extract_indices(py_list, &idx1, &idx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "third argument must be a list of length 1 or 2");
        return NULL;
    }

// #pragma omp parallel for
    for (size_t i = 0; i < GROUPLENGTH; ++i) {
        mpz_init(vector[i]);
        py_to_mpz(vector[i], PyList_GET_ITEM(py_vals, i));
        if (encode(vector[i], vector[i], idx1, idx2, 0) == FAILURE) {
            err = 1;
        }
    }

    if (err) {
        PyErr_SetString(PyExc_RuntimeError, "encoding failed");
    } else {
        int fnamelen = strlen(g_dir) + strlen(name) + 2;
        char *fname = (char *) malloc(sizeof(char) * fnamelen);
        if (fname == NULL) {
            err = 1;
            goto cleanup;
        }
        (void) snprintf(fname, fnamelen, "%s/%s", g_dir, name);
        (void) save_mpz_vector(fname, vector, GROUPLENGTH);
        free(fname);
    }

    for (size_t i = 0; i < GROUPLENGTH; ++i) {
        mpz_clear(vector[i]);
    }

 cleanup:
    if (err) {
        return NULL;
    } else {
        Py_RETURN_NONE;
    }

}

static PyObject *
obf_encode_level(PyObject *self, PyObject *args)
{
    PyObject *py_zero_m, *py_one_m;
    PyObject *py_zero_set, *py_one_set;
    int zeroidx1, zeroidx2, oneidx1, oneidx2;
    char *fname = NULL;
    int fnamelen;
    int err = 0;
    long inp, idx;
    mpz_t zero[GROUPLENGTH_SQ];
    mpz_t one[GROUPLENGTH_SQ];
    mpz_t z_inp;

    if (!PyArg_ParseTuple(args, "llOOOO", &idx, &inp, &py_zero_m, &py_one_m,
                          &py_zero_set, &py_one_set)) {
        return NULL;
    }
    if (extract_indices(py_zero_set, &zeroidx1, &zeroidx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "second argument must be a list of length 1 or 2");
        return NULL;
    }
    if (extract_indices(py_one_set, &oneidx1, &oneidx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "third argument must be a list of length 1 or 2");
        return NULL;
    }

    mpz_init_set_ui(z_inp, inp);

#pragma omp parallel for
    for (size_t ctr = 0; ctr < 2 * GROUPLENGTH_SQ; ++ctr) {
        PyObject *py_array;
        int sidx1, sidx2;
        mpz_t *val;
        size_t i;

        if (ctr < GROUPLENGTH_SQ) {
            i = ctr;
            val = &zero[i];
            py_array = py_zero_m;
            sidx1 = zeroidx1;
            sidx2 = zeroidx2;
        } else {
            i = ctr - GROUPLENGTH_SQ;
            val = &one[i];
            py_array = py_one_m;
            sidx1 = oneidx1;
            sidx2 = oneidx2;
        }

        mpz_init(*val);
        py_to_mpz(*val, PyList_GET_ITEM(py_array, i));
        if (encode(*val, *val, sidx1, sidx2, 0) == FAILURE) {
            fprintf(stderr, "encoding failed!\n");
            err = 1;
        }
    }

    fnamelen = strlen(g_dir) + 10;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL) {
        err = 1;
        goto exit;
    }

    (void) snprintf(fname, fnamelen, "%s/%ld.input", g_dir, idx);
    (void) save_mpz_scalar(fname, z_inp);
    (void) snprintf(fname, fnamelen, "%s/%ld.zero", g_dir, idx);
    (void) save_mpz_vector(fname, zero, GROUPLENGTH_SQ);
    (void) snprintf(fname, fnamelen, "%s/%ld.one", g_dir, idx);
    (void) save_mpz_vector(fname, one, GROUPLENGTH_SQ);

    for (int i = 0; i < GROUPLENGTH_SQ; ++i) {
        mpz_clear(zero[i]);
        mpz_clear(one[i]);
    }
exit:
    if (fname)
        free(fname);

    mpz_clear(z_inp);

    if (err)
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
}

static PyObject *
obf_save(PyObject *self, PyObject *args)
{
    char *fname;
    int len;
    mpz_t nu;

    len = strlen(g_dir) + 5;

    fname = (char *) malloc(sizeof(char) * len);
    if (fname == NULL)
        return NULL;

    mpz_init_set_ui(nu, g_nu);

    // save nu
    (void) snprintf(fname, len, "%s/nu", g_dir);
    (void) save_mpz_scalar(fname, nu);
    // save x0
    (void) snprintf(fname, len, "%s/x0", g_dir);
    (void) save_mpz_scalar(fname, g_x0);
    // save pzt
    (void) snprintf(fname, len, "%s/pzt", g_dir);
    (void) save_mpz_scalar(fname, g_pzt);

    mpz_clear(nu);

    free(fname);

    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate(PyObject *self, PyObject *args)
{
    char *input = NULL;
    char *fname = NULL;
    int fnamelen;
    int iszero;
    mpz_t comp[GROUPLENGTH_SQ];
    mpz_t tmp1[GROUPLENGTH_SQ];
    mpz_t tmp2[GROUPLENGTH_SQ];
    mpz_t s[GROUPLENGTH];
    mpz_t t[GROUPLENGTH];
    mpz_t p1, p2, tmp;
    long bplen;
    
    if (!PyArg_ParseTuple(args, "ssl", &g_dir, &input, &bplen))
        return NULL;

    fnamelen = strlen(g_dir) + 10;

    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_init(tmp);
    mpz_init(p1);
    mpz_init(p2);
    for (int i = 0; i < GROUPLENGTH; ++i) {
        mpz_init(s[i]);
        mpz_init(t[i]);
    }
    for (int i = 0; i < GROUPLENGTH_SQ; ++i) {
        mpz_init(comp[i]);
        mpz_init(tmp1[i]);
        mpz_init(tmp2[i]);
    }

    (void) snprintf(fname, fnamelen, "%s/s_enc", g_dir);
    (void) load_mpz_vector(fname, s, GROUPLENGTH);
    (void) snprintf(fname, fnamelen, "%s/t_enc", g_dir);
    (void) load_mpz_vector(fname, t, GROUPLENGTH);
    (void) snprintf(fname, fnamelen, "%s/p_enc", g_dir);
    (void) load_mpz_scalar(fname, p2);

    for (int level = 0; level < bplen; ++level) {
        unsigned int input_idx;

        // find out the input bit for the given level
        (void) snprintf(fname, fnamelen, "%s/%d.input", g_dir, level);
        (void) load_mpz_scalar(fname, tmp);
        input_idx = mpz_get_ui(tmp);
        if (input_idx < 0 || input_idx >= strlen(input)) {
            return NULL;        // XXX: memory leak
        }
        assert(input[input_idx] == '0' || input[input_idx] == '1');
        // load in appropriate matrix for the given input value
        if (input[input_idx] == '0') {
            (void) snprintf(fname, fnamelen, "%s/%d.zero", g_dir, level);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", g_dir, level);
        }
        if (level == 0) {
            (void) load_mpz_vector(fname, comp, GROUPLENGTH_SQ);
        } else {
            (void) load_mpz_vector(fname, tmp1, GROUPLENGTH_SQ);
            // XXX: make in-place matrix multiplication to save the additional
            // copying needed
            mat_mult(tmp2, comp, tmp1);
            for (int ctr = 0; ctr < GROUPLENGTH_SQ; ++ctr) {
                mpz_set(comp[ctr], tmp2[ctr]);
            }
        }
        if (input[input_idx] == '0') {
            (void) snprintf(fname, fnamelen, "%s/%d.a0_enc", g_dir, level);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.a1_enc", g_dir, level);
        }
        load_mpz_scalar(fname, tmp);
        mpz_mul(p2, p2, tmp);
    }

    mat_mult_by_vects(p1, s, comp, t);

    mpz_sub(tmp, p1, p2);
    iszero = is_zero(tmp);
    // iszero = is_zero(comp[1]);

    mpz_clear(tmp);
    mpz_clear(p1);
    mpz_clear(p2);
    for (int i = 0; i < GROUPLENGTH; ++i) {
        mpz_clear(s[i]);
        mpz_clear(t[i]);
    }
    for (int i = 0; i < GROUPLENGTH_SQ; ++i) {
        mpz_clear(comp[i]);
        mpz_clear(tmp1[i]);
        mpz_clear(tmp2[i]);
    }

    if (fname)
        free(fname);

    return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyMethodDef
ObfMethods[] = {
    {"setup", obf_setup, METH_VARARGS,
     "Set up obfuscator."},
    {"genprime", obf_genprime, METH_VARARGS,
     "Generate a random prime."},
    {"encode_scalar", obf_encode_scalar, METH_VARARGS,
     "Encode a scalar."},
    {"encode_vector", obf_encode_vector, METH_VARARGS,
     "Encode a vector."},
    {"encode_level", obf_encode_level, METH_VARARGS,
     "Encode a branching program level."},
    {"save", obf_save, METH_VARARGS,
     "Save obfuscation values."},
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

    // seed = 16360538864245260915;
    fprintf(stderr, "SEED = %lu\n", seed);

    gmp_randinit_default(g_rng);
    gmp_randseed_ui(g_rng, seed);

 cleanup:
    if (file != -1)
        (void) close(file);
}
