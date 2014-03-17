#include <Python.h>
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

static void
mpz_genrandom(mpz_t rnd, const long nbits)
{
    mpz_t one, rndtmp;

    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_init(rndtmp);

    mpz_urandomb(rndtmp, g_rng, nbits);
    mpz_ior(rnd, rndtmp, one);

    mpz_clear(one);
    mpz_clear(rndtmp);
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
fastutils_genparams(PyObject *self, PyObject *args)
{
    const long alpha, beta, eta, kappa;
    long i;
    PyObject *py_x0, *py_pzt, *py_g0;
    mpz_t *zs;

    if (!PyArg_ParseTuple(args, "llllllO", &g_n, &alpha, &beta, &eta, &kappa,
                          &g_nzs, &py_g0))
        return NULL;

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

    // Generate p_i's and g_'s, as well as compute x0
    {
#pragma omp parallel for private(i)
        for (i = 0; i < g_n; ++i) {
            mpz_t p_unif;

            mpz_init(p_unif);

            // XXX: not uniform primes
            mpz_urandomb(p_unif, g_rng, eta);
            mpz_nextprime(g_ps[i], p_unif);
            if (i == 0) {
                py_to_mpz(g_gs[0], py_g0);
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
    PyObject *py_x0, *py_pzt;

    if (!PyArg_ParseTuple(args, "OO", &py_x0, &py_pzt))
        return NULL;

    py_to_mpz(g_x0, py_x0);
    py_to_mpz(g_pzt, py_pzt);

    Py_RETURN_NONE;
}

static int
encode(mpz_t out, const mpz_t in, const long rho, const long idx1,
       const long idx2)
{
    mpz_t res, r, tmp;
    int i;

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
        mpz_genrandom(r, rho);
        mpz_mul(tmp, r, g_gs[i]);
        if (i == 0) {
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
    const long rho, idx1, idx2;
    PyObject *py_val, *py_out;
    mpz_t val;

    if (!PyArg_ParseTuple(args, "Olll", &py_val, &rho, &idx1, &idx2))
        return NULL;

    mpz_init(val);
    py_to_mpz(val, py_val);
    if (encode(val, val, rho, idx1, idx2) == SUCCESS) {
        py_out = mpz_to_py(val);
        mpz_clear(val);
        return py_out;
    } else {
        mpz_clear(val);
        return NULL;
    }
}

static PyObject *
fastutils_encode_vector(PyObject *self, PyObject *args)
{
    const long rho, idx;
    PyObject *py_vals, *py_outs;
    Py_ssize_t i, len;
    int err = 0;

    if (!PyArg_ParseTuple(args, "Oll", &py_vals, &rho, &idx))
        return NULL;
    if (!PyList_Check(py_vals))
        return NULL;

    len = PyList_GET_SIZE(py_vals);
    if ((py_outs = PyList_New(len)) == NULL)
        return NULL;

#pragma omp parallel for private(i)
    for (i = 0; i < len; ++i) {
        mpz_t val;

        mpz_init(val);
        py_to_mpz(val, PyList_GET_ITEM(py_vals, i));
        if (encode(val, val, rho, idx, -1) == FAILURE) {
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
    const long rho;
    long zeroidx1 = -1, zeroidx2 = -1, oneidx1 = -1, oneidx2 = -1;
    PyObject *py_vals, *py_outs, *py_zeros, *py_ones;
    Py_ssize_t i, len, half;

    if (!PyArg_ParseTuple(args, "OlOO", &py_vals, &rho, &py_zeros, &py_ones))
        return NULL;
    if (!PyList_Check(py_vals))
        return NULL;
    if (!PyList_Check(py_zeros))
        return NULL;
    if (!PyList_Check(py_ones))
        return NULL;

    len = PyList_GET_SIZE(py_zeros);
    if (len == 0) {
        return NULL;
    } else if (len == 1) {
        zeroidx1 = PyLong_AsLong(PyList_GET_ITEM(py_zeros, 0));
    } else if (len == 2) {
        zeroidx1 = PyLong_AsLong(PyList_GET_ITEM(py_zeros, 0));
        zeroidx2 = PyLong_AsLong(PyList_GET_ITEM(py_zeros, 1));
    } else {
        return NULL;
    }

    len = PyList_GET_SIZE(py_ones);
    if (len == 0) {
        return NULL;
    } else if (len == 1) {
        oneidx1 = PyLong_AsLong(PyList_GET_ITEM(py_ones, 0));
    } else if (len == 2) {
        oneidx1 = PyLong_AsLong(PyList_GET_ITEM(py_ones, 0));
        oneidx2 = PyLong_AsLong(PyList_GET_ITEM(py_ones, 1));
    } else {
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
            encode(val, val, rho, zeroidx1, zeroidx2);
        } else {
            encode(val, val, rho, oneidx1, oneidx2);
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
    PyObject *py_c;
    mpz_t c;
    int ret;

    if (!PyArg_ParseTuple(args, "Ol", &py_c, &nu))
        return NULL;

    mpz_init(c);

    py_to_mpz(c, py_c);

    mpz_mul(c, c, g_pzt);
    mpz_mod_near(c, c, g_x0);

    /* fprintf(stderr, "is_zero length = %ld\n", mpz_sizeinbase(c, 2)); */

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
    (void) Py_InitModule("fastutils", FastutilsMethods);

    gmp_randinit_default(g_rng);
}
