#include <Python.h>

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
static mpz_t g_x0;
static mpz_t *g_ps;
static mpz_t *g_gs;
static mpz_t g_pzt;
static mpz_t *g_crt_coeffs;
static mpz_t *g_zinvs;
static long g_nzs;
static long g_rho;

/* static double */
/* current_time(void) */
/* { */
/*     struct timeval t; */
/*     gettimeofday(&t, NULL); */
/*     return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0)); */
/* } */

inline static PyObject *
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

inline static void
py_to_mpz(mpz_t out, PyObject *in)
{
    (void) mpz_set_pylong(out, in);
}

inline static int
check_pylist(PyObject *list)
{
    PyObject *item;
    Py_ssize_t i, len;

    if (!PyList_Check(list)) {
        return FAILURE;
    }

    len = PyList_GET_SIZE(list);
    for (i = 0; i < len; ++i) {
        item = PyList_GET_ITEM(list, i);
        if (!PyLong_Check(item)) {
            return FAILURE;
        }
    }
    return SUCCESS;
}

inline static int
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
        break;
    }
    return SUCCESS;
}

static void
mpz_genrandom(mpz_t rnd, const long nbits)
{
    mpz_t one;

    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, g_rng, (mp_bitcnt_t) nbits);
    /* mpz_sub(rnd, rnd, one); */

    mpz_clear(one);
}

static void
mpz_mod_near(mpz_t out, const mpz_t a, const mpz_t b)
{
    mpz_t res, shift;

    mpz_init(res);
    mpz_init(shift);

    mpz_mod(res, a, b);
    mpz_tdiv_q_2exp(shift, b, 1);
    if (mpz_cmp(res, shift) > 0) {
        mpz_sub(res, res, b);
    }

    mpz_set(out, res);

    mpz_clear(res);
    mpz_clear(shift);
}

static PyObject *
fastutils_genprime(PyObject *self, PyObject *args)
{
    const long bitlength;
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
fastutils_genparams(PyObject *self, PyObject *args)
{
    const long alpha, beta, eta, kappa;
    PyObject *py_x0, *py_pzt, *py_gs;
    Py_ssize_t gs_len;
    mpz_t *zs;
    long i;

    if (!PyArg_ParseTuple(args, "lllllllO", &g_n, &alpha, &beta, &eta, &kappa,
                          &g_rho, &g_nzs, &py_gs))
        return NULL;
    if (check_pylist(py_gs) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "eighth argument must be a list of longs");
        return NULL;
    }

    mpz_init_set_ui(g_x0, 1);
    mpz_init_set_ui(g_pzt, 0);

    g_ps = (mpz_t *) malloc(sizeof(mpz_t) * g_n);
    g_gs = (mpz_t *) malloc(sizeof(mpz_t) * g_n);
    g_crt_coeffs = (mpz_t *) malloc(sizeof(mpz_t) * g_n);
    if (g_ps == NULL || g_gs == NULL || g_crt_coeffs == NULL)
        return NULL;

    zs = (mpz_t *) malloc(sizeof(mpz_t) * g_nzs);
    g_zinvs = (mpz_t *) malloc(sizeof(mpz_t) * g_nzs);
    if (zs == NULL || g_zinvs == NULL)
        return NULL;

    for (i = 0; i < g_n; ++i) {
        mpz_init(g_ps[i]);
        mpz_init(g_gs[i]);
        mpz_init(g_crt_coeffs[i]);
    }

    for (i = 0; i < g_nzs; ++i) {
        mpz_init(zs[i]);
        mpz_init(g_zinvs[i]);
    }

    gs_len = PyList_GET_SIZE(py_gs);

    // Generate p_i's and g_'s, as well as compute x0
    {
#pragma omp parallel for private(i)
        for (i = 0; i < g_n; ++i) {
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
        py_x0 = mpz_to_py(g_x0);
    }

    // Generate CRT coefficients
    {
#pragma omp parallel for default(shared) private(i)
        for (i = 0; i < g_n; ++i) {
            mpz_t q;
            mpz_init(q);

            mpz_tdiv_q(q, g_x0, g_ps[i]);
            mpz_invert(g_crt_coeffs[i], q, g_ps[i]);
            mpz_mul(g_crt_coeffs[i], g_crt_coeffs[i], q);

            mpz_clear(q);
        }
    }

    // Generate z_i's
#pragma omp parallel for private(i)
    for (i = 0; i < g_nzs; ++i) {
        do {
            mpz_urandomm(zs[i], g_rng, g_x0);
        } while (mpz_invert(g_zinvs[i], zs[i], g_x0) == 0);
    }

    // Generate pzt
    {
        mpz_t zk;

        mpz_init_set_ui(zk, 1);

        // Compute z^k
        for (i = 0; i < g_nzs; ++i) {
            mpz_mul(zk, zk, zs[i]);
            mpz_mod(zk, zk, g_x0);
        }

#pragma omp parallel for private(i)
        for (i = 0; i < g_n; ++i) {
            mpz_t tmp, x0pi, rnd;

            mpz_init(tmp);
            mpz_init(x0pi);
            mpz_init(rnd);

            // Compute (((g_i)^{-1} mod p_i) * z^k mod p_i) * r_i * (x_0 / p_i)
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

        py_pzt = mpz_to_py(g_pzt);

        mpz_clear(zk);
    }

    for (i = 0; i < g_nzs; ++i) {
        mpz_clear(zs[i]);
    }
    free(zs);

    return PyTuple_Pack(2, py_x0, py_pzt);
}

static PyObject *
fastutils_loadparams(PyObject *self, PyObject *args)
{
    PyObject *py_x0 = NULL, *py_pzt = NULL;

    if (!PyArg_ParseTuple(args, "OO", &py_x0, &py_pzt))
        return NULL;
    if (!PyLong_Check(py_x0)) {
        PyErr_SetString(PyExc_RuntimeError, "first argument must be a long");
        return NULL;
    }
    if (!PyLong_Check(py_pzt)) {
        PyErr_SetString(PyExc_RuntimeError, "second argument must be a long");
        return NULL;
    }

    py_to_mpz(g_x0, py_x0);
    py_to_mpz(g_pzt, py_pzt);

    Py_RETURN_NONE;
}

static int
encode(mpz_t out, const mpz_t in, const long idx1,
       const long idx2, const long slot)
{
    mpz_t res, r, tmp;
    long i;

    if (idx1 >= g_nzs || idx2 >= g_nzs)
        return FAILURE;
    if (idx1 < 0 && idx2 < 0)
        return FAILURE;
    if (idx1 == idx2)
        return FAILURE;

    mpz_init(res);
    mpz_init(r);
    mpz_init(tmp);

    mpz_set_ui(res, 0);

    for (i = 0; i < g_n; ++i) {
        mpz_genrandom(r, g_rho);
        mpz_mul(tmp, r, g_gs[i]);
        if (i == slot) {
            mpz_add(tmp, tmp, in);
        }
        mpz_mul(tmp, tmp, g_crt_coeffs[i]);
        mpz_add(res, res, tmp);
    }
    mpz_mod(res, res, g_x0);
    if (idx1 >= 0) {
        mpz_mul(res, res, g_zinvs[idx1]);
        mpz_mod(res, res, g_x0);
    }
    if (idx2 >= 0) {
        mpz_mul(res, res, g_zinvs[idx2]);
        mpz_mod(res, res, g_x0);
    }

    mpz_set(out, res);

    mpz_clear(res);
    mpz_clear(r);
    mpz_clear(tmp);

    return SUCCESS;
}

static PyObject *
fastutils_encode_scalar(PyObject *self, PyObject *args)
{
    PyObject *py_val = NULL, *py_list = NULL, *py_out = NULL;
    int idx1, idx2;
    long slot;
    mpz_t val;

    mpz_init(val);

    if (!PyArg_ParseTuple(args, "OlO", &py_val, &slot, &py_list))
        goto error;
    if (!PyLong_Check(py_val)) {
        PyErr_SetString(PyExc_RuntimeError, "first argument must be a long");
        goto error;
    }
    if (!PyList_Check(py_list)) {
        PyErr_SetString(PyExc_RuntimeError, "third argument must be a list");
        goto error;
    }

    if (extract_indices(py_list, &idx1, &idx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "third argument must be a list of length 1 or 2");
        goto error;
    }

    py_to_mpz(val, py_val);
    if (encode(val, val, idx1, idx2, slot) == SUCCESS) {
        py_out = mpz_to_py(val);
    } else {
        PyErr_SetString(PyExc_RuntimeError, "encoding failed");
        goto error;
    }

 error:
    mpz_clear(val);

    return py_out;
}

static PyObject *
fastutils_encode_vector(PyObject *self, PyObject *args)
{
    const long slot;
    PyObject *py_vals = NULL, *py_list = NULL, *py_outs = NULL;
    Py_ssize_t i, len;
    int idx1, idx2, err = 0;

    if (!PyArg_ParseTuple(args, "OlO", &py_vals, &slot, &py_list)) {
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

    len = PyList_GET_SIZE(py_vals);
    if ((py_outs = PyList_New(len)) == NULL) {
        return NULL;
    }

#pragma omp parallel for private(i)
    for (i = 0; i < len; ++i) {
        mpz_t val;

        mpz_init(val);
        py_to_mpz(val, PyList_GET_ITEM(py_vals, i));
        if (encode(val, val, idx1, idx2, slot) == FAILURE) {
            err = 1;
        }
        PyList_SET_ITEM(py_outs, i, mpz_to_py(val));
        mpz_clear(val);
    }

    if (err)
        return NULL;
    else
        return py_outs;
}

static PyObject *
fastutils_encode_layer(PyObject *self, PyObject *args)
{
    PyObject *py_vals, *py_outs, *py_zeros, *py_ones;
    Py_ssize_t i, len, half;
    int zeroidx1, zeroidx2, oneidx1, oneidx2;
    long slot;

    if (!PyArg_ParseTuple(args, "OlOO", &py_vals, &slot, &py_zeros, &py_ones))
        return NULL;
    if (check_pylist(py_vals) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "first argument must be a list of longs");
        return NULL;
    }

    if (extract_indices(py_zeros, &zeroidx1, &zeroidx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "second argument must be a list of length 1 or 2");
        return NULL;
    }
    if (extract_indices(py_ones, &oneidx1, &oneidx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "third argument must be a list of length 1 or 2");
        return NULL;
    }

    len = PyList_GET_SIZE(py_vals);
    half = len >> 1;

    if ((py_outs = PyList_New(len)) == NULL)
        return NULL;

#pragma omp parallel for private(i)
    for (i = 0; i < len; ++i) {
        mpz_t val;

        mpz_init(val);
        py_to_mpz(val, PyList_GET_ITEM(py_vals, i));
        if (i < half) {
            encode(val, val, zeroidx1, zeroidx2, slot);
        } else {
            encode(val, val, oneidx1, oneidx2, slot);
        }
        PyList_SET_ITEM(py_outs, i, mpz_to_py(val));
        mpz_clear(val);
    }

    return py_outs;
}

static PyObject *
fastutils_is_zero(PyObject *self, PyObject *args)
{
    const long nu;
    PyObject *py_c = NULL;
    mpz_t c;
    int ret;

    if (!PyArg_ParseTuple(args, "Ol", &py_c, &nu))
        return NULL;
    if (!PyLong_Check(py_c)) {
        PyErr_SetString(PyExc_RuntimeError, "first argument must be a long");
        return NULL;
    }

    mpz_init(c);

    py_to_mpz(c, py_c);

    mpz_mul(c, c, g_pzt);
    mpz_mod_near(c, c, g_x0);

    ret = (mpz_sizeinbase(c, 2) < (mpz_sizeinbase(g_x0, 2) - nu)) ? 1 : 0;

    mpz_clear(c);

    if (ret)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyMethodDef
FastutilsMethods[] = {
    {"genparams", fastutils_genparams, METH_VARARGS,
     "Generate MLM parameters."},
    {"loadparams", fastutils_loadparams, METH_VARARGS,
     "Load MLM parameters."},
    {"genprime", fastutils_genprime, METH_VARARGS,
     "Generate a random prime."},
    {"encode_scalar", fastutils_encode_scalar, METH_VARARGS,
     "Encode a scalar."},
    {"encode_vector", fastutils_encode_vector, METH_VARARGS,
     "Encode a vector."},
    {"encode_layer", fastutils_encode_layer, METH_VARARGS,
     "Encode a branching program layer."},
    {"is_zero", fastutils_is_zero, METH_VARARGS,
     "Zero test."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initfastutils(void)
{
    int file;
    unsigned long seed;

    (void) Py_InitModule("fastutils", FastutilsMethods);

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
