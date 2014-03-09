#include <Python.h>
#include <gmp.h>
#include <omp.h>
#include <sys/time.h>

#include "mpz_pylong.h"

static gmp_randstate_t g_rng;
static long n;
static mpz_t x0;
static mpz_t *ps;
static mpz_t *gs;
static mpz_t z;
static mpz_t zinv;
static mpz_t pzt;

static double
current_time(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0));
}

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
    const long alpha, beta, eta, kappa;
    long i;
    PyObject *py_x0, *py_pzt;

    if (!PyArg_ParseTuple(args, "lllll", &n, &alpha, &beta, &eta, &kappa))
        return NULL;

    mpz_init_set_ui(x0, 1);
    mpz_init(z);
    mpz_init_set_ui(pzt, 0);
    ps = (mpz_t *) malloc(sizeof(mpz_t) * n);
    gs = (mpz_t *) malloc(sizeof(mpz_t) * n);
    // XXX: never free'd
    if (ps == NULL || gs == NULL)
        return NULL;
    for (i = 0; i < n; ++i) {
        mpz_init(ps[i]);
        mpz_init(gs[i]);
    }

    // Generate p_i's and g_'s, as well as compute x0
    {
        mpz_t x0tmp;

        mpz_init_set(x0tmp, x0);

#pragma omp parallel for private(i)
        for (i = 0; i < n; ++i) {
            mpz_t p_unif;

            mpz_init(p_unif);

            mpz_urandomb(p_unif, g_rng, alpha);
            mpz_nextprime(gs[i], p_unif);

            mpz_urandomb(p_unif, g_rng, eta);
            mpz_nextprime(ps[i], p_unif);

#pragma omp critical
            {
                mpz_mul(x0, x0tmp, ps[i]);
                mpz_set(x0tmp, x0);
            }

            mpz_clear(p_unif);
        }
        py_x0 = mpz_to_py(x0);

        mpz_clear(x0tmp);
    }

    // Generate z
    do {
        mpz_urandomm(z, g_rng, x0);
    } while (mpz_invert(zinv, z, x0) == 0);

    // Generate pzt
    {
        mpz_t zkappa, pzttmp, tmp;

        mpz_init_set_ui(zkappa, 1);
        mpz_init_set_ui(pzttmp, 0);
        mpz_init(tmp);

        for (i = 0; i < kappa; ++i) {
            mpz_mul(tmp, zkappa, z);
            mpz_mod(zkappa, tmp, x0);
        }

#pragma omp parallel for private(i)
        for (i = 0; i < n; ++i) {
            mpz_t input, tmp1, tmp2, rnd;

            mpz_init(input);
            mpz_init(tmp1);
            mpz_init(tmp2);
            mpz_init(rnd);

            mpz_invert(input, gs[i], ps[i]);
            mpz_mul(tmp1, input, zkappa);
            mpz_mod(tmp2, tmp1, ps[i]);
            genrandom(rnd, beta);
            mpz_mul(tmp1, tmp2, rnd);
            mpz_div(tmp2, x0, ps[i]);
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
        }
        py_pzt = mpz_to_py(pzt);

        mpz_clear(zkappa);
        mpz_clear(pzttmp);
        mpz_clear(tmp);
    }

    mpz_clear(x0);
    mpz_clear(pzt);

    return PyTuple_Pack(2, py_x0, py_pzt);
}

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
    const long rho;
    PyObject *py_msgs, *py_out;
    mpz_t x, m;
    double start, end;
    int i;

    if (!PyArg_ParseTuple(args, "lO", &rho, &py_msgs))
        return NULL;
    if (n != PySequence_Size(py_msgs))
        return NULL;
    // XXX: compare msgs[i] with gs[i]

    mpz_init(x);
    mpz_init(m);

    start = current_time();

/* #pragma omp parallel for private(i) */
    for (i = 0; i < n; ++i) {
        mpz_t msg, r, tmp;

        mpz_init(msg);
        mpz_init(r);
        mpz_init(tmp);

        py_to_mpz(msg, PyList_GET_ITEM(py_msgs, i));

        genrandom(r, rho);
        mpz_addmul(msg, r, gs[i]);
        mpz_mul(tmp, msg, zinv);
        mpz_mod(tmp, tmp, ps[i]);
/* #pragma omp critical */
        {
            if (i == 0) {
                mpz_set(x, tmp);
                mpz_set(m, ps[0]);
            } else {
                crt(x, x, tmp, m, ps[i]);
                mpz_lcm(m, m, ps[i]);
            }
        }

        mpz_clear(msg);
        mpz_clear(r);
        mpz_clear(tmp);
    }
    mpz_mod(x, x, m);

    end = current_time();

    fprintf(stderr, "Took: %f s\n", (end - start));

    py_out = mpz_to_py(x);

    mpz_clear(x);
    mpz_clear(m);

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
