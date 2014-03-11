#include <Python.h>
#include <gmp.h>
#include <omp.h>
#include <sys/time.h>

#include "mpz_pylong.h"

static gmp_randstate_t g_rng;
static long g_n;
static mpz_t g_x0;
static mpz_t *g_ps;
static mpz_t *g_gs;
static mpz_t g_z;
static mpz_t g_zinv;
static mpz_t g_pzt;
static mpz_t *g_crt_coeffs;

/* static double */
/* current_time(void) */
/* { */
/*     struct timeval t; */
/*     gettimeofday(&t, NULL); */
/*     return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0)); */
/* } */

inline static PyObject *
mpz_to_py(const mpz_t x)
{
    PyObject *outs, *out;
    char *buffer;

    buffer = mpz_get_str(NULL, 10, x);
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
genrandom(mpz_t rnd, const long nbits)
{
    mpz_t one, rndtmp;

    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_init(rndtmp);
    mpz_urandomb(rndtmp, g_rng, nbits);
    mpz_ior(rnd, rndtmp, one);
    mpz_clear(one);
    mpz_clear(rndtmp);
}

static PyObject *
fastutils_genparams(PyObject *self, PyObject *args)
{
    const long alpha, beta, eta, kappa;
    long i;
    PyObject *py_x0, *py_pzt;

    if (!PyArg_ParseTuple(args, "lllll", &g_n, &alpha, &beta, &eta, &kappa))
        return NULL;

    mpz_init_set_ui(g_x0, 1);
    mpz_init(g_z);
    mpz_init_set_ui(g_pzt, 0);
    g_ps = (mpz_t *) malloc(sizeof(mpz_t) * g_n);
    g_gs = (mpz_t *) malloc(sizeof(mpz_t) * g_n);
    g_crt_coeffs = (mpz_t *) malloc(sizeof(mpz_t) * g_n);
    // XXX: never free'd
    if (g_ps == NULL || g_gs == NULL || g_crt_coeffs == NULL)
        return NULL;
    for (i = 0; i < g_n; ++i) {
        mpz_init(g_ps[i]);
        mpz_init(g_gs[i]);
        mpz_init(g_crt_coeffs[i]);
    }

    // Generate p_i's and g_'s, as well as compute x0
    {
#pragma omp parallel for private(i)
        for (i = 0; i < g_n; ++i) {
            mpz_t p_unif;

            mpz_init(p_unif);

            mpz_urandomb(p_unif, g_rng, eta);
            mpz_nextprime(g_ps[i], p_unif);

            mpz_urandomb(p_unif, g_rng, alpha);
            mpz_nextprime(g_gs[i], p_unif);
            
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

    // Generate z
    do {
        mpz_urandomm(g_z, g_rng, g_x0);
    } while (mpz_invert(g_zinv, g_z, g_x0) == 0);

    // Generate pzt
    {
        mpz_t zkappa;

        mpz_init_set_ui(zkappa, 1);

        for (i = 0; i < kappa; ++i) {
            mpz_mul(zkappa, zkappa, g_z);
            mpz_mod(zkappa, zkappa, g_x0);
        }

#pragma omp parallel for private(i)
        for (i = 0; i < g_n; ++i) {
            mpz_t input, tmp, rnd;

            mpz_init(input);
            mpz_init(tmp);
            mpz_init(rnd);

            mpz_invert(input, g_gs[i], g_ps[i]);
            mpz_mul(input, input, zkappa);
            mpz_mod(input, input, g_ps[i]);
            genrandom(rnd, beta);
            mpz_mul(input, input, rnd);
            mpz_div(tmp, g_x0, g_ps[i]);
            mpz_mul(input, input, tmp);
#pragma omp critical
            {
                mpz_add(g_pzt, g_pzt, input);
            }

            mpz_clear(input);
            mpz_clear(tmp);
            mpz_clear(rnd);
        }
        py_pzt = mpz_to_py(g_pzt);

        mpz_clear(zkappa);
    }

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

static void
encode(mpz_t out, const mpz_t in, const long rho)
{
    mpz_t res, r, tmp;
    int i;

    mpz_init(res);
    mpz_init(r);
    mpz_init(tmp);

    mpz_set_ui(res, 0);

    for (i = 0; i < g_n; ++i) {
        genrandom(r, rho);
        mpz_mul(tmp, r, g_gs[i]);
        if (i == 0) {
            mpz_add(tmp, tmp, in);
        }
        mpz_mul(tmp, tmp, g_crt_coeffs[i]);
        mpz_add(res, res, tmp);
    }
    mpz_mul(res, res, g_zinv);
    mpz_mod(res, res, g_x0);

    mpz_set(out, res);

    mpz_clear(res);
    mpz_clear(r);
    mpz_clear(tmp);
}

static PyObject *
fastutils_encode(PyObject *self, PyObject *args)
{
    const long rho;
    PyObject *py_msg, *py_out;
    mpz_t msg;

    if (!PyArg_ParseTuple(args, "Ol", &py_msg, &rho))
        return NULL;

    mpz_init(msg);
    py_to_mpz(msg, py_msg);
    if (mpz_cmp(msg, g_gs[0]) >= 0) {
        py_out = NULL;
        goto cleanup;
    }
    encode(msg, msg, rho);
    py_out = mpz_to_py(msg);

 cleanup:
    mpz_clear(msg);
    return py_out;
}

static PyObject *
fastutils_encode_list(PyObject *self, PyObject *args)
{
    const long rho;
    PyObject *py_vals, *py_outs;
    Py_ssize_t i, len;

    if (!PyArg_ParseTuple(args, "Ol", &py_vals, &rho))
        return NULL;

    len = PySequence_Size(py_vals);

    if ((py_outs = PyList_New(len)) == NULL)
        return NULL;

#pragma omp parallel for private(i)
    for (i = 0; i < len; ++i) {
        mpz_t val;

        mpz_init(val);
        py_to_mpz(val, PyList_GET_ITEM(py_vals, i));
        encode(val, val, rho);
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
    mpz_t c, cmp;
    int ret;

    if (!PyArg_ParseTuple(args, "Ol", &py_c, &nu))
        return NULL;

    mpz_init(c);
    mpz_init(cmp);

    py_to_mpz(c, py_c);

    mpz_tdiv_q_2exp(cmp, g_x0, nu);

    mpz_mul(c, c, g_pzt);
    mpz_mod(c, c, g_x0);
    if (mpz_cmpabs(c, cmp) < 0)
        ret = 1;
    else
        ret = 0;

    mpz_clear(c);
    mpz_clear(cmp);

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
    {"encode", fastutils_encode, METH_VARARGS,
     "Encode a value."},
    {"encode_list", fastutils_encode_list, METH_VARARGS,
     "Encode a list of values."},
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
