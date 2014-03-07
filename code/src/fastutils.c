#include <Python.h>
#include <gmp.h>
#include <omp.h>

#include "mpz_pylong.h"

static gmp_randstate_t g_rng;

inline static PyObject *
mpz_to_py(mpz_t x)
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
genrandom(mpz_t rnd, long nbits)
{
    mpz_t one, rndtmp;

    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_init(rndtmp);
    mpz_urandomb(rndtmp, g_rng, nbits);
    /* mpz_sub(rnd, rndtmp, one); */
    mpz_ior(rnd, rndtmp, one);
    mpz_clear(one);
    mpz_clear(rndtmp);
}

static PyObject *
fastutils_genparams(PyObject *self, PyObject *args)
{
    const long n, alpha, beta, eta, kappa;
    long i;
    PyObject *py_ps, *py_gs, *py_x0, *py_z, *py_zinv, *py_pzt;
    mpz_t x0, z, pzt;
    
    if (!PyArg_ParseTuple(args, "lllll", &n, &alpha, &beta, &eta, &kappa))
        return NULL;

    if ((py_ps = PyList_New(n)) == NULL)
        return NULL;
    if ((py_gs = PyList_New(n)) == NULL)
        return NULL;

    mpz_init_set_ui(x0, 1);
    mpz_init(z);
    mpz_init_set_ui(pzt, 0);

    // Generate p_i's and g_'s, as well as compute x0
    {
        mpz_t x0tmp;

        mpz_init_set(x0tmp, x0);

#pragma omp parallel for private(i)
        for (i = 0; i < n; ++i) {
            PyObject *py_p, *py_g;
            mpz_t p_tmp, p_unif;

            mpz_init(p_tmp);
            mpz_init(p_unif);

            mpz_urandomb(p_unif, g_rng, alpha);
            mpz_nextprime(p_tmp, p_unif);
            py_g = mpz_to_py(p_tmp);

            mpz_urandomb(p_unif, g_rng, eta);
            mpz_nextprime(p_tmp, p_unif);
            py_p = mpz_to_py(p_tmp);

#pragma omp critical
            {
                PyList_SET_ITEM(py_ps, i, py_p);
                PyList_SET_ITEM(py_gs, i, py_g);
                mpz_mul(x0, x0tmp, p_tmp);
                mpz_set(x0tmp, x0);
            }

            mpz_clear(p_tmp);
            mpz_clear(p_unif);
        }
        py_x0 = mpz_to_py(x0);

        mpz_clear(x0tmp);
    }

    // Generate z
    {
        mpz_t zinv;
        int ret;

        mpz_init(zinv);

        do {
            mpz_urandomm(z, g_rng, x0);
            ret = mpz_invert(zinv, z, x0);
        } while (ret == 0);

        py_z = mpz_to_py(z);
        py_zinv = mpz_to_py(zinv);

        mpz_clear(zinv);
    }

    // Generate pzt
    {
        mpz_t zkappa, pzttmp, tmp;

        mpz_init(tmp);

        mpz_init_set_ui(zkappa, 1);
        mpz_init_set_ui(pzttmp, 0);

        for (i = 0; i < kappa; ++i) {
            mpz_mul(tmp, zkappa, z);
            mpz_mod(zkappa, tmp, x0);
        }

#pragma omp parallel for private(i)
        for (i = 0; i < n; ++i) {
            mpz_t input, tmp1, tmp2, rnd, g, p;

            mpz_init(input);
            mpz_init(tmp1);
            mpz_init(tmp2);
            mpz_init(rnd);
            mpz_init(g);
            mpz_init(p);

            py_to_mpz(p, PyList_GET_ITEM(py_ps, i));
            py_to_mpz(g, PyList_GET_ITEM(py_gs, i));
            mpz_invert(input, g, p);
            mpz_mul(tmp1, input, zkappa);
            mpz_mod(tmp2, tmp1, p);
            genrandom(rnd, beta);
            mpz_mul(tmp1, tmp2, rnd);
            mpz_div(tmp2, x0, p);
            mpz_mul(input, tmp1, tmp2);
#pragma omp critical
            {
                mpz_add(pzt, pzttmp, input);
                mpz_set(pzttmp, pzt);
            }

            mpz_clear(input);
            mpz_clear(tmp1);
            mpz_clear(tmp2);
            mpz_clear(rnd);
            mpz_clear(g);
            mpz_clear(p);
        }
        py_pzt = mpz_to_py(pzt);

        mpz_clear(zkappa);
        mpz_clear(pzttmp);
        mpz_clear(tmp);
    }

    mpz_clear(x0);
    mpz_clear(z);
    mpz_clear(pzt);

    return PyTuple_Pack(6, py_x0, py_ps, py_gs, py_z, py_zinv, py_pzt);
}

/* static int */
/* convert_list(int length, PyObject *in, mpz_t *out) */
/* { */
/*     int i; */
/*     PyObject *seq; */

/*     seq = PySequence_Fast(in, "not a list"); */
/*     for (i = 0; i < length; ++i) { */
/*         (void) mpz_set_pylong(out[i], PySequence_Fast_GET_ITEM(seq, i)); */
/*     } */
/*     return 0; */
/* } */

static void
crt(mpz_t out, mpz_t a, mpz_t b, mpz_t m, mpz_t n)
{
    mpz_t g, alpha, beta, q, r, tmp, tmp2;

    mpz_init(g);
    mpz_init(alpha);
    mpz_init(beta);
    mpz_init(q);
    mpz_init(r);
    mpz_init(tmp);
    mpz_init(tmp2);

    mpz_gcdext(g, alpha, beta, m, n);
    mpz_sub(tmp, b, a);
    mpz_cdiv_qr(q, r, tmp, g);
    // TODO: check if r != 0
    mpz_mul(tmp, q, alpha);
    mpz_mul(tmp2, tmp, m);
    mpz_add(tmp, a, tmp2);
    mpz_lcm(tmp2, m, n);

    mpz_mod(out, tmp, tmp2);

    mpz_clear(g);
    mpz_clear(alpha);
    mpz_clear(beta);
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(tmp);
    mpz_clear(tmp2);
}

static PyObject *
fastutils_encode(PyObject *self, PyObject *args)
{
    const long n, rho;
    PyObject *py_msgs, *py_ps, *py_gs, *py_zinv, *py_out, *py_rs;
    mpz_t zinv, x, m, xtmp, mtmp;
    int i;

    if (!PyArg_ParseTuple(args, "llOOOOO", &n, &rho, &py_msgs, &py_ps, &py_gs, &py_zinv, &py_rs))
        return NULL;

    if (n != PySequence_Size(py_msgs)
     || n != PySequence_Size(py_ps)
     || n != PySequence_Size(py_gs))
        return NULL;

    mpz_init(zinv);
    py_to_mpz(zinv, py_zinv);

    mpz_init(x);
    mpz_init(m);
    mpz_init(xtmp);
    mpz_init(mtmp);

/* #pragma omp parallel for private(i) */
    for (i = 0; i < n; ++i) {
        mpz_t g, msg, p, r, tmp1, tmp2;

        mpz_init(g);
        mpz_init(msg);
        mpz_init(p);
        mpz_init(r);
        mpz_init(tmp1);
        mpz_init(tmp2);

        py_to_mpz(g, PyList_GET_ITEM(py_gs, i));
        py_to_mpz(msg, PyList_GET_ITEM(py_msgs, i));
        py_to_mpz(p, PyList_GET_ITEM(py_ps, i));
        py_to_mpz(r, PyList_GET_ITEM(py_rs, i));
        genrandom(r, rho);
        mpz_addmul(msg, r, g);
        mpz_mul(tmp1, msg, zinv);
        mpz_mod(tmp2, tmp1, p);
/* #pragma omp critical */
        {
            if (i == 0) {
                mpz_set(x, tmp2);
                mpz_set(m, p);
            } else {
                mpz_set(xtmp, x);
                mpz_set(mtmp, m);

                crt(x, x, tmp2, mtmp, p);
                mpz_lcm(m, mtmp, p);
            }
        }

        mpz_clear(g);
        mpz_clear(msg);
        mpz_clear(p);
        mpz_clear(r);
        mpz_clear(tmp1);
        mpz_clear(tmp2);
    }
    mpz_set(xtmp, x);
    mpz_mod(x, xtmp, m);

    py_out = mpz_to_py(x);

    mpz_clear(x);
    mpz_clear(m);
    mpz_clear(xtmp);
    mpz_clear(mtmp);

    return py_out;
}

static PyMethodDef
FastutilsMethods[] = {
    {"genparams", fastutils_genparams, METH_VARARGS,
     "Generate MLM parameters."},
    {"encode", fastutils_encode, METH_VARARGS,
     "Encode vector."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initfastutils(void)
{
    (void) Py_InitModule("fastutils", FastutilsMethods);

    gmp_randinit_default(g_rng);
}
