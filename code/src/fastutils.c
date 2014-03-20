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
encode(mpz_t out, const PyObject *ins, const PyObject *slots,
       const Py_ssize_t length, const long idx1, const long idx2)
{
    mpz_t res, r, tmp, item;
    long i, j;

    if (idx1 >= g_nzs || idx2 >= g_nzs)
        return FAILURE;
    if (idx1 < 0 && idx2 < 0)
        return FAILURE;
    if (idx1 == idx2)
        return FAILURE;

    mpz_init(res);
    mpz_init(r);
    mpz_init(tmp);
    mpz_init(item);

    mpz_set_ui(res, 0);

    for (i = 0; i < g_n; ++i) {
        mpz_genrandom(r, g_rho);
        mpz_mul(tmp, r, g_gs[i]);
        for (j = 0; j < length; ++j) {
            if (i == PyLong_AsLong(PyList_GET_ITEM(slots, j))) {
                py_to_mpz(item, PyList_GET_ITEM(ins, j));
                mpz_add(tmp, tmp, item);
            }
        }
        /* if (i == slot) { */
        /*     mpz_add(tmp, tmp, in); */
        /* } */
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
    mpz_clear(item);

    return SUCCESS;
}

static PyObject *
fastutils_encode_scalar(PyObject *self, PyObject *args)
{
    PyObject *py_vals, *py_list, *py_slots, *py_out = NULL;
    Py_ssize_t lvals, lslots;
    int idx1, idx2;
    mpz_t out;

    mpz_init(out);

    if (!PyArg_ParseTuple(args, "OOO", &py_vals, &py_slots, &py_list))
        goto error;
    if (check_pylist(py_vals) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "first argument must be a list of longs");
        goto error;
    }
    if (!PyList_Check(py_slots)) {
        PyErr_SetString(PyExc_RuntimeError, "second argument must be a list");
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

    lvals = PyList_GET_SIZE(py_vals);
    lslots = PyList_GET_SIZE(py_slots);
    if (lvals != lslots) {
        PyErr_SetString(PyExc_RuntimeError,
                        "first two arguments must be of same length");
        goto error;
    }
    
    if (encode(out, py_vals, py_slots, lvals, idx1, idx2) == SUCCESS) {
        py_out = mpz_to_py(out);
    } else {
        PyErr_SetString(PyExc_RuntimeError, "encoding failed");
        goto error;
    }

 error:
    mpz_clear(out);

    return py_out;
}

static PyObject *
fastutils_encode_vector(PyObject *self, PyObject *args)
{
    const long veclen;
    PyObject *py_vectors, *py_slots, *py_list, *py_outs;
    Py_ssize_t numvectors;
    long i;
    int idx1, idx2;
    PyObject **lists;

    if (!PyArg_ParseTuple(args, "lOOO", &veclen, &py_vectors, &py_slots, &py_list)) {
        return NULL;
    }
    if (!PyList_Check(py_vectors)) {
        PyErr_SetString(PyExc_RuntimeError, "second argument must be a list");
        return NULL;
    }
    if (!PyList_Check(py_slots)) {
        PyErr_SetString(PyExc_RuntimeError, "third argument must be a list");
        return NULL;
    }
    if (!PyList_Check(py_list)) {
        PyErr_SetString(PyExc_RuntimeError, "fourth argument must be a list");
        return NULL;
    }

    if (extract_indices(py_list, &idx1, &idx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "fourth argument must be a list of length 1 or 2");
        return NULL;
    }

    numvectors = PyList_GET_SIZE(py_vectors);
    if (numvectors != PyList_GET_SIZE(py_slots)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "second and third arguments must be of same length");
        return NULL;
    }

    if ((py_outs = PyList_New(veclen)) == NULL) {
        return NULL;
    }

    lists = (PyObject **) malloc(sizeof(PyObject *) * veclen);
    for (i = 0; i < veclen; ++i) {
        lists[i] = PyList_New(numvectors);
    }

#pragma omp parallel for
    for (i = 0; i < veclen; ++i) {
        mpz_t val;
        PyObject *py_tmp;
        Py_ssize_t j;

        mpz_init(val);
        py_tmp = lists[i];

        for (j = 0; j < numvectors; ++j) {
            PyList_SET_ITEM(py_tmp, j,
                            PyList_GET_ITEM(PyList_GET_ITEM(py_vectors, j), i));
        }
        encode(val, py_tmp, py_slots, numvectors, idx1, idx2);
        PyList_SET_ITEM(py_outs, i, mpz_to_py(val));

        mpz_clear(val);
    }

    free(lists);

    return py_outs;
}

static PyObject *
fastutils_encode_layer(PyObject *self, PyObject *args)
{
    PyObject *py_layers, *py_slots, *py_zeros, *py_ones, *py_outs;
    Py_ssize_t numlayers, half;
    int zeroidx1, zeroidx2, oneidx1, oneidx2;
    const long layerlen;
    long i;
    PyObject **lists;

    if (!PyArg_ParseTuple(args, "lOOOO", &layerlen, &py_layers, &py_slots,
                          &py_zeros, &py_ones)) {
        return NULL;
    }
    if (!PyList_Check(py_layers)) {
        PyErr_SetString(PyExc_RuntimeError, "second argument must be a list");
        return NULL;
    }
    if (!PyList_Check(py_slots)) {
        PyErr_SetString(PyExc_RuntimeError, "third argument must be a list");
        return NULL;
    }

    if (extract_indices(py_zeros, &zeroidx1, &zeroidx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "fourth argument must be a list of length 1 or 2");
        return NULL;
    }
    if (extract_indices(py_ones, &oneidx1, &oneidx2) == FAILURE) {
        PyErr_SetString(PyExc_RuntimeError,
                        "fifth argument must be a list of length 1 or 2");
        return NULL;
    }

    numlayers = PyList_GET_SIZE(py_layers);
    if (numlayers != PyList_GET_SIZE(py_slots)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "second and third arguments must be of the same length");
        return NULL;
    }

    half = layerlen >> 1;

    if ((py_outs = PyList_New(layerlen)) == NULL) {
        return NULL;
    }

    lists = (PyObject **) malloc(sizeof(PyObject *) * layerlen);
    for (i = 0; i < layerlen; ++i) {
        lists[i] = PyList_New(numlayers);
    }

#pragma omp parallel for private(i)
    for (i = 0; i < layerlen; ++i) {
        mpz_t val;
        PyObject *py_tmp;
        Py_ssize_t j;

        mpz_init(val);
        /* py_tmp = PyList_New(numlayers); */
        py_tmp = lists[i];

        for (j = 0; j < numlayers; ++j) {
            PyList_SET_ITEM(py_tmp, j,
                            PyList_GET_ITEM(PyList_GET_ITEM(py_layers, j), i));
        }
        if (i < half) {
            encode(val, py_tmp, py_slots, numlayers, zeroidx1, zeroidx2);
        } else {
            encode(val, py_tmp, py_slots, numlayers, oneidx1, oneidx2);
        }
        PyList_SET_ITEM(py_outs, i, mpz_to_py(val));

        mpz_clear(val);
    }

    free(lists);

    return py_outs;
}

static PyObject *
fastutils_is_zero(PyObject *self, PyObject *args)
{
    const long nu;
    PyObject *py_c;
    int ret;

    if (!PyArg_ParseTuple(args, "Ol", &py_c, &nu))
        return NULL;
    if (!PyLong_Check(py_c)) {
        PyErr_SetString(PyExc_RuntimeError, "first argument must be a long");
        return NULL;
    }

    {
        mpz_t c;
        mpz_init(c);

        py_to_mpz(c, py_c);

        mpz_mul(c, c, g_pzt);
        mpz_mod_near(c, c, g_x0);

        ret = (mpz_sizeinbase(c, 2) < (mpz_sizeinbase(g_x0, 2) - nu)) ? 1 : 0;

        mpz_clear(c);
    }

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
