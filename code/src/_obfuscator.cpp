#include <Python.h>
#include "pyutils.h"
#include "utils.h"
#include "thpool.h"
#include "thpool_fns.h"

#include <aesrand.h>
#include <mmap/mmap.h>
#include <mmap/mmap_clt.h>
#include <mmap/mmap_gghlite.h>
#include <omp.h>

struct state {
    threadpool thpool;
    unsigned long secparam;
    enum mmap_e type;
    mmap_sk mmap;
    const mmap_vtable *vtable;
    aes_randstate_t rand;
    char *dir;
    long nzs;
};

static void
mult_matrices(const mmap_vtable *vtable, const mmap_pp *pp, mmap_enc *result,
              const mmap_enc *left, const mmap_enc *right, long m, long n,
              long p)
{
    mmap_enc *tmparray;
    double start, end;

    start = current_time();
    tmparray = (mmap_enc *) malloc(sizeof(mmap_enc) * m * p);
    for (int i = 0; i < m * p; ++i) {
        vtable->enc->init(&tmparray[i], pp);
    }
#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            mmap_enc tmp, sum;
            vtable->enc->init(&tmp, pp);
            vtable->enc->init(&sum, pp);
            for (int k = 0; k < n; ++k) {
                vtable->enc->mul(&tmp, pp,
                                 &left[k * m + (i * m + j) % m],
                                 &right[k + n * ((i * m + j) / m)]);
                vtable->enc->add(&sum, pp, &sum, &tmp);
            }
            vtable->enc->set(&tmparray[i * n + j], &sum);
            vtable->enc->clear(&tmp);
            vtable->enc->clear(&sum);
        }
    }
    for (int i = 0; i < m * p; ++i) {
        vtable->enc->set(&result[i], &tmparray[i]);
        vtable->enc->clear(&tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, " Multiplying took: %f\n", end - start);
}

static void
state_cleanup(struct state *s)
{
    if (s) {
        s->vtable->sk->clear(&s->mmap);
        aes_randclear(s->rand);
        thpool_destroy(s->thpool);
    }
    free(s);
}

static void
state_destructor(PyObject *self)
{
    state_cleanup((struct state *) PyCapsule_GetPointer(self, NULL));
}

//
//
// Python functions
//
//

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long kappa, size, nthreads, ncores;
    struct state *s = NULL;
    PyObject *py_primes, *py_state;

    s = (struct state *) calloc(1, sizeof(struct state));
    if (s == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "memory allocation failed");
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "llllslll", &s->secparam, &kappa, &size,
                          &s->nzs, &s->dir, &s->type, &nthreads, &ncores)) {
        PyErr_SetString(PyExc_RuntimeError, "unable to parse input");
        goto error;
    }
    if (kappa <= 0 || size < 0 || s->nzs <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "invalid input");
        return NULL;
    }

    switch (s->type) {
    case MMAP_CLT:
        s->vtable = &clt13_vtable;
        break;
    case MMAP_GGHLITE:
        s->vtable = &gghlite_vtable;
        break;
    default:
        PyErr_SetString(PyExc_RuntimeError, "invalid mmap setting");
        return NULL;
    }

    (void) aes_randinit(s->rand);
    s->thpool = thpool_init(1);
    (void) omp_set_num_threads(1);
    // s->thpool = thpool_init(nthreads);
    // (void) omp_set_num_threads(ncores);

    if (g_verbose) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
    }


    s->vtable->sk->init(&s->mmap, s->secparam, kappa, s->nzs, s->rand);
    {
        char fname[100];
        FILE *fp;

        snprintf(fname, 100, "%s/params", s->dir);
        fp = fopen(fname, "w+b");
        s->vtable->pp->fwrite(s->vtable->sk->pp(&s->mmap), fp);
        fclose(fp);
    }

    py_primes = PyList_New(1);
    if (s->type == MMAP_CLT) {
        PyList_SetItem(py_primes, 0, mpz_to_py(((clt_state *) &s->mmap)->gs[0]));
    }
    py_state = PyCapsule_New((void *) s, NULL, state_destructor);
    return PyTuple_Pack(2, py_state, py_primes);

error:
    state_cleanup(s);
    return NULL;
}

static void
_obf_encode_layers(struct state *s, long idx,
                   long inp, long nrows, long ncols,
                   PyObject *py_zero_ms, PyObject *py_one_ms)
{
    struct write_layer_s *wl_s;
    mmap_enc *zero, *one;
    double start;
    char idx_s[10];
    const mmap_pp *pp = s->vtable->sk->pp(&s->mmap);

    start = current_time();

    (void) snprintf(idx_s, 10, "%ld", idx);

    zero = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
    one = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
    for (ssize_t i = 0; i < nrows * ncols; ++i) {
        s->vtable->enc->init(&zero[i], pp);
        s->vtable->enc->init(&one[i], pp);
    }

    wl_s = (struct write_layer_s *) malloc(sizeof(write_layer_s));
    wl_s->vtable = s->vtable;
    wl_s->dir = s->dir;
    wl_s->zero = zero;
    wl_s->one = one;
    wl_s->inp = inp;
    wl_s->idx = idx;
    wl_s->nrows = nrows;
    wl_s->ncols = ncols;
    wl_s->start = start;

    (void) thpool_add_tag(s->thpool, idx_s, 2 * nrows * ncols,
                          thpool_write_layer, wl_s);

    for (Py_ssize_t ctr = 0; ctr < 2 * nrows * ncols; ++ctr) {
        PyObject *py_array;
        fmpz_t *plaintext;
        mmap_enc *enc;
        size_t i;
        struct encode_elem_s *args;

        if (ctr < nrows * ncols) {
            i = ctr;
            enc = &zero[i];
            py_array = py_zero_ms;
        } else {
            i = ctr - nrows * ncols;
            enc = &one[i];
            py_array = py_one_ms;
        }

        plaintext = (fmpz_t *) malloc(sizeof(fmpz_t));
        fmpz_init(*plaintext);
        py_to_fmpz(*plaintext, PyList_GET_ITEM(PyList_GET_ITEM(py_array, 0), i));

        args = (struct encode_elem_s *) malloc(sizeof(struct encode_elem_s));
        args->vtable = s->vtable;
        args->sk = &s->mmap;
        args->enc = enc;
        args->n = 1;
        args->plaintext = plaintext;
        args->group = (int *) calloc(s->nzs , sizeof(int));
        if (s->type == MMAP_CLT)
            args->group[idx] = 1;
        thpool_add_work(s->thpool, thpool_encode_elem, (void *) args, idx_s);
    }
}

static PyObject *
obf_encode_layers(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_zero_ms, *py_one_ms;
    long inp, idx, nrows, ncols;

    struct state *s;

    if (!PyArg_ParseTuple(args, "OllllOO", &py_state, &idx, &nrows, &ncols,
                          &inp, &py_zero_ms, &py_one_ms))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    _obf_encode_layers(s, idx, inp, nrows, ncols, py_zero_ms, py_one_ms);

    Py_RETURN_NONE;
}

static int
_obf_evaluate(const mmap_vtable *vtable, mmap_pp *pp, char *dir, char *input, long bplen)
{
    char fname[100];
    FILE *fp;
    mmap_enc *result = NULL;
    long nrows, ncols, nrows_prev;
    int err = 0, iszero = -1;
    double start, end;

    for (long layer = 0; layer < bplen; ++layer) {
        unsigned int inp;
        mmap_enc *left, *right;

        start = current_time();

        // determine the size of the matrix
        (void) snprintf(fname, 100, "%s/%ld.nrows", dir, layer);
        fp = fopen(fname, "r+b");
        fread(&nrows, sizeof nrows, 1, fp);
        fclose(fp);

        (void) snprintf(fname, 100, "%s/%ld.ncols", dir, layer);
        fp = fopen(fname, "r+b");
        fread(&ncols, sizeof ncols, 1, fp);
        fclose(fp);

        // find out the input bit for the given layer
        (void) snprintf(fname, 100, "%s/%ld.input", dir, layer);
        fp = fopen(fname, "r+b");
        fread(&inp, sizeof inp, 1, fp);
        fclose(fp);

        if (inp >= strlen(input)) {
            fprintf(stderr, "Error: invalid input: %d >= %ld\n", inp, strlen(input));
            err = 1;
            break;
        }
        if (input[inp] != '0' && input[inp] != '1') {
            fprintf(stderr, "Error: input must be 0 or 1, got %d\n", input[inp]);
            err = 1;
            break;
        }
        // load in appropriate matrix for the given input value
        if (input[inp] == '0') {
            (void) snprintf(fname, 100, "%s/%ld.zero", dir, layer);
        } else {
            (void) snprintf(fname, 100, "%s/%ld.one", dir, layer);
        }

        if (layer == 0) {
            result = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
            fp = fopen(fname, "r+b");
            for (int i = 0; i < nrows * ncols; ++i) {
                vtable->enc->init(&result[i], pp);
                vtable->enc->fread(&result[i], fp);
            }
            fclose(fp);
            nrows_prev = nrows;
        } else {
            left = result;
            right = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows * ncols);
            fp = fopen(fname, "r+b");
            for (int i = 0; i < nrows * ncols; ++i) {
                vtable->enc->init(&right[i], pp);
                vtable->enc->fread(&right[i], fp);
            }
            fclose(fp);
            result = (mmap_enc *) malloc(sizeof(mmap_enc) * nrows_prev * ncols);
            for (int i = 0; i < nrows_prev * ncols; ++i) {
                vtable->enc->init(&result[i], pp);
            }
            mult_matrices(vtable, pp, result, left, right, nrows_prev, nrows, ncols);
            for (int i = 0; i < nrows_prev * nrows; ++i) {
                vtable->enc->clear(&left[i]);
            }
            for (int i = 0; i < nrows * ncols; ++i) {
                vtable->enc->clear(&right[i]);
            }
            free(left);
            free(right);
        }
        end = current_time();

        if (g_verbose)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n", end - start);
    }

    if (!err) {
        start = current_time();
        iszero = vtable->enc->is_zero(&result[1], pp);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Zero test: %f\n", end - start);
    }

    if (result) {
        for (int i = 0; i < nrows_prev * ncols; ++i) {
            vtable->enc->clear(&result[i]);
        }
        free(result);
    }

    return iszero;
}

static PyObject *
obf_evaluate(PyObject *self, PyObject *args)
{
    char *dir, *input;
    char fname[100];
    int iszero = -1;
    enum mmap_e type;
    long bplen, nthreads;
    const mmap_vtable *vtable;
    mmap_pp pp;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "sslll", &dir, &input, &bplen, &type, &nthreads)) {
        PyErr_SetString(PyExc_RuntimeError, "error parsing arguments");
        return NULL;
    }

    (void) omp_set_num_threads(1);
    // (void) omp_set_num_threads(nthreads);

    switch (type) {
    case MMAP_CLT:
        vtable = &clt13_vtable;
        break;
    case MMAP_GGHLITE:
        vtable = &gghlite_vtable;
        break;
    default:
        PyErr_SetString(PyExc_RuntimeError, "invalid mmap setting");
        return NULL;
    }

    (void) snprintf(fname, 100, "%s/params", dir);
    fp = fopen(fname, "r+b");
    vtable->pp->fread(&pp, fp);
    fclose(fp);

    iszero = _obf_evaluate(vtable, &pp, dir, input, bplen);
    if (iszero == -1) {
        PyErr_SetString(PyExc_RuntimeError, "zero test failed");
        return NULL;
    } else {
        return Py_BuildValue("i", iszero ? 0 : 1);
    }
}

static PyObject *
obf_wait(PyObject *self, PyObject *args)
{
    PyObject *py_state;
    struct state *s;

    if (!PyArg_ParseTuple(args, "O", &py_state))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    thpool_wait(s->thpool);

    Py_RETURN_NONE;
}

static PyMethodDef
ObfMethods[] = {
    {"verbose", obf_verbose, METH_VARARGS,
     "Set verbosity."},
    {"setup", obf_setup, METH_VARARGS,
     "Set up obfuscator."},
    {"encode_layers", obf_encode_layers, METH_VARARGS,
     "Encode a branching program layer in each slot."},
    {"max_mem_usage", obf_max_mem_usage, METH_VARARGS,
     "Print out the maximum memory usage."},
    {"evaluate", obf_evaluate, METH_VARARGS,
     "Evaluate the obfuscation."},
    {"wait", obf_wait, METH_VARARGS,
     "Wait for threadpool to empty."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
