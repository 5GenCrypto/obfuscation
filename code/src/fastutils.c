#include <Python.h>
#include <gmp.h>
#include <omp.h>

#include "mpz_pylong.h"

static gmp_randstate_t g_rng;
/* static mpz_t *g_ps; */
/* static mpz_t *g_gs; */
/* static mpz_t g_zinv; */
/* static long g_n, g_rho; */

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
py_to_mpz(PyObject *in, mpz_t out)
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

            py_to_mpz(PyList_GET_ITEM(py_ps, i), p);
            py_to_mpz(PyList_GET_ITEM(py_gs, i), g);
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

/* static PyObject * */
/* fastutils_init(PyObject *self, PyObject *args) */
/* { */
/*     const long n; */
/*     PyObject *py_ps, *py_gs, *py_zinv; */
/*     int i; */

/*     if (!PyArg_ParseTuple(args, "lOOOl", &n, &py_ps, &py_gs, &py_zinv, &g_rho)) */
/*         return NULL; */

/*     { */
/*         int ps_len, gs_len; */
/*         ps_len = PySequence_Size(py_ps); */
/*         gs_len = PySequence_Size(py_gs); */
/*         if (n != ps_len || n != gs_len) */
/*             return NULL; */
/*         g_n = n; */
/*     } */

/*     g_ps = (mpz_t *) malloc(sizeof(mpz_t) * g_n); */
/*     g_gs = (mpz_t *) malloc(sizeof(mpz_t) * g_n); */
/*     for (i = 0; i < g_n; ++i) { */
/*         mpz_init(g_ps[i]); */
/*         mpz_init(g_gs[i]); */
/*     } */
/*     mpz_init(g_zinv); */

/*     (void) convert_list(g_n, py_ps, g_ps); */
/*     (void) convert_list(g_n, py_gs, g_gs); */

/*     /\* mpz_init(result); *\/ */

/*     /\* for (i = 0; i < num; ++i) { *\/ */
/*     /\*     res += (m[i] + g[i] * rnd) * zinv % p[i]; *\/ */
/*     /\* } *\/ */

/*     /\* mpz_mod(res, res, x0); *\/ */

/*     return Py_BuildValue(""); */
/* } */

/* static PyObject * */
/* fastutils_encode(PyObject *self, PyObject *args) */
/* { */
/*     PyObject *py_ms, *py_ms_seq; */
/*     mpz_t *ms; */
/*     mpz_t result, rnd; */
/*     int i; */

/*     if (!PyArg_ParseTuple(args, "O", &py_ms)) */
/*         return NULL; */

/*     if (g_n != PySequence_Size(py_ms)) */
/*         return NULL; */

/*     mpz_init(result); */
/*     mpz_init(rnd); */

/*     ms = (mpz_t *) malloc(sizeof(mpz_t) * g_n); */
/*     py_ms_seq = PySequence_Fast(py_ms, "not a list"); */
/*     for (i = 0; i < g_n; ++i) { */
/*         mpz_init(ms[i]); */
/*         (void) mpz_set_pylong(ms[i], PySequence_Fast_GET_ITEM(py_ms_seq, i)); */
/*         genrandom(rnd, g_rho); */
/*         /\* (ms[i] + g_gs[i] * rnd) * g_zinv; *\/ */
/*     } */
    
/* } */

/* static PyObject * */
/* fastutils_clear(PyObject *self, PyObject *args) */
/* { */
/*     int i; */

/*     for (i = 0; i < g_n; ++i) { */
/*         mpz_clear(g_ps[i]); */
/*         mpz_clear(g_gs[i]); */
/*     } */
/*     mpz_clear(g_zinv); */
/* } */

static PyMethodDef
FastutilsMethods[] = {
    /* {"init", fastutils_init, METH_VARARGS, */
    /*  "TODO."}, */
    /* {"genprimes", fastutils_genprimes, METH_VARARGS, */
    /*  "Generate random primes."}, */
    {"genparams", fastutils_genparams, METH_VARARGS,
     "Generate MLM parameters."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initfastutils(void)
{
    (void) Py_InitModule("fastutils", FastutilsMethods);

    gmp_randinit_default(g_rng);
}
