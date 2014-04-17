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

/* static double */
/* current_time(void) */
/* { */
/*     struct timeval t; */
/*     gettimeofday(&t, NULL); */
/*     return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0)); */
/* } */

inline static void *
mymalloc(const size_t size)
{
    void * r;
    if ((r = malloc(size)) == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory");
    }
    return r;
}

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
load_mpz_scalar(const char *fname, mpz_t x)
{
    int ret = SUCCESS;
    FILE *f;

    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return FAILURE;
    }
    (void) mpz_inp_raw(x, f);

    (void) fclose(f);
    return ret;
}

inline static int
save_mpz_scalar(const char *fname, const mpz_t x)
{
    int ret = SUCCESS;
    FILE *f;

    if ((f = fopen(fname, "w+")) == NULL) {
        perror(fname);
        return FAILURE;
    }
    (void) mpz_out_raw(f, x);

    (void) fclose(f);
    return ret;
}

inline static int
load_mpz_vector(const char *fname, mpz_t *m, const int len)
{
    int ret = SUCCESS;
    FILE *f;
    
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return FAILURE;
    }

    for (int i = 0; i < len; ++i) {
        (void) mpz_inp_raw(m[i], f);
    }

    (void) fclose(f);
    return ret;
}

inline static int
save_mpz_vector(const char *fname, const mpz_t *m, const int len)
{
    int ret = SUCCESS;
    FILE *f;

    if ((f = fopen(fname, "w+")) == NULL) {
        perror(fname);
        return FAILURE;
    }

    for (int i = 0; i < len; ++i) {
        (void) mpz_out_raw(f, m[i]);
    }

    (void) fclose(f);
    return ret;
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
        return FAILURE;
    }
    return SUCCESS;
}

inline static void
mpz_genrandom(mpz_t rnd, const long nbits)
{
    mpz_t one;
    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, g_rng, (mp_bitcnt_t) nbits);
    mpz_clear(one);
}

inline static void
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

static void
encode(mpz_t out, const PyObject *in, const long row, const long idx1,
       const long idx2)
{
    mpz_t r, tmp;

    mpz_init(r);
    mpz_init(tmp);

    mpz_set_ui(out, 0);

// #pragma omp parallel for
    for (long i = 0; i < g_n; ++i) {
        // mpz_t r;
        // mpz_init(r);
        mpz_genrandom(r, g_rho);
        mpz_mul(tmp, r, g_gs[i]);
        py_to_mpz(r, PyList_GET_ITEM(PyList_GET_ITEM(in, i), row));
        mpz_add(tmp, tmp, r);
        mpz_mul(tmp, tmp, g_crt_coeffs[i]);
// #pragma omp critical
        {
            mpz_add(out, out, tmp);
        }
        // mpz_clear(r);
    }
    mpz_mod(out, out, g_x0);
    if (idx1 >= 0) {
        mpz_mul(out, out, g_zinvs[idx1]);
        mpz_mod(out, out, g_x0);
    }
    if (idx2 >= 0) {
        mpz_mul(out, out, g_zinvs[idx2]);
        mpz_mod(out, out, g_x0);
    }

    mpz_clear(r);
    mpz_clear(tmp);
}


inline static void
mat_mult(mpz_t *out, const mpz_t *a, const mpz_t *b, int size)
{
#pragma omp parallel for
    for (int ctr = 0; ctr < size * size; ++ctr) {
        mpz_t tmp, sum;
        mpz_init(tmp);
        mpz_init_set_ui(sum, 0);
        for (int i = 0; i < size; ++i) {
            mpz_mul(tmp,
                    a[i * size + ctr % size],
                    b[i + size * (ctr / size)]);
            mpz_add(sum, sum, tmp);
        }
        mpz_set(out[ctr], sum);
        mpz_clear(tmp);
        mpz_clear(sum);
    }
}

inline static void
mat_mult_by_vects(mpz_t out, const mpz_t *s, const mpz_t *m, const mpz_t *t,
                  int size)
{
    mpz_set_ui(out, 0);

#pragma omp parallel for
    for (int col = 0; col < size; ++col) {
        mpz_t tmp;
        mpz_t sum;
        mpz_init(tmp);
        mpz_init_set_ui(sum, 0);
        for (int row = 0; row < size; ++row) {
            int elem = col * size + row;
            mpz_mul(tmp, s[row], m[elem]);
            mpz_add(sum, sum, tmp);
        }
        mpz_mul(tmp, sum, t[col]);
#pragma omp critical
        {
            mpz_add(out, out, tmp);
        }
        mpz_clear(tmp);
        mpz_clear(sum);
    }
}

inline static int
is_zero(mpz_t c)
{
    mpz_t tmp;
    int ret;

    mpz_init(tmp);
    mpz_mul(tmp, c, g_pzt);
    mpz_mod_near(tmp, tmp, g_x0);
    ret = (mpz_sizeinbase(tmp, 2) < (mpz_sizeinbase(g_x0, 2) - g_nu)) ? 1 : 0;
    mpz_clear(tmp);
    return ret;
}

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

    // XXX: this appears to be the most expensive step!

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

    // save size, nu, x0, and pzt to file

    {
        char *fname;
        int len;
        mpz_t tmp;

        len = strlen(g_dir) + 10;

        fname = (char *) malloc(sizeof(char) * len);
        if (fname == NULL)
            return NULL;

        mpz_init(tmp);

        // save size
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

//
// Encode N scalars across all slots of the MLM
//
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

//
// Encode N vectors across all slots of the MLM
//
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

    // We assume that all vectors are the same size, and thus just grab the size
    // of the first vector
    size = PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0));
    vector = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    if (vector == NULL)
        return NULL;

#pragma omp parallel for
    for (Py_ssize_t i = 0; i < size; ++i) {
        mpz_init(vector[i]);
        encode(vector[i], py_vectors, i, idx1, idx2);
    }

    {
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

//
// Encode N layers across all slots of the MLM
//
static PyObject *
obf_encode_layers(PyObject *self, PyObject *args)
{
    PyObject *py_zero_ms, *py_one_ms;
    PyObject *py_zero_set, *py_one_set;
    int zeroidx1, zeroidx2, oneidx1, oneidx2;
    int err = 0;
    long inp, idx;
    Py_ssize_t size;
    mpz_t *zero;
    mpz_t *one;

    if (!PyArg_ParseTuple(args, "llOOOO", &idx, &inp, &py_zero_ms, &py_one_ms,
                          &py_zero_set, &py_one_set))
        return NULL;
    (void) extract_indices(py_zero_set, &zeroidx1, &zeroidx2);
    (void) extract_indices(py_one_set, &oneidx1, &oneidx2);

    size = PyList_GET_SIZE(PyList_GET_ITEM(py_zero_ms, 0));
    zero = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    one = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    if (!zero || !one)
        return NULL;

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
        encode(*val, py_array, i, idx1, idx2);
    }

    {
        mpz_t z;
        int fnamelen = strlen(g_dir) + 10;
        char *fname = (char *) mymalloc(sizeof(char) * fnamelen);
        if (fname == NULL) {
            err = 1;
        } else {
            mpz_init_set_ui(z, inp);
            (void) snprintf(fname, fnamelen, "%s/%ld.input", g_dir, idx);
            (void) save_mpz_scalar(fname, z);
            (void) snprintf(fname, fnamelen, "%s/%ld.zero", g_dir, idx);
            (void) save_mpz_vector(fname, zero, size);
            (void) snprintf(fname, fnamelen, "%s/%ld.one", g_dir, idx);
            (void) save_mpz_vector(fname, one, size);
            free(fname);
            mpz_clear(z);
        }
    }

    for (int i = 0; i < size; ++i) {
        mpz_clear(zero[i]);
        mpz_clear(one[i]);
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
    s = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    t = (mpz_t *) mymalloc(sizeof(mpz_t) * size);
    if (!comp || !tmp1 || !tmp2 || !s || !t) {
        err = 1;
        goto cleanup;
    }

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
            (void) load_mpz_scalar(fname, tmp);
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

 cleanup:
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
