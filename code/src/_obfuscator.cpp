#include "utils.h"
#include "pyutils.h"
#include "clt_mlm.h"
#include "thpool_fns.h"

#include "C-Thread-Pool/thpool.h"

#include <omp.h>

struct state {
    threadpool thpool;
    struct clt_mlm_state mlm;
    char *dir;
};

void
state_destructor(PyObject *self)
{
    struct state *s;

    s = (struct state *) PyCapsule_GetPointer(self, NULL);
    if (s) {
        clt_mlm_cleanup(&s->mlm);
        thpool_destroy(s->thpool);
    }
    free(s);
}


static int
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
        return 1;
    }
    return 0;
}

//
//
// Python functions
//
//

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long kappa, size, nthreads;
    struct state *s = NULL;
    long *pows = NULL;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    if (!PyArg_ParseTuple(args, "llllsl", &s->mlm.secparam, &kappa, &size,
                          &s->mlm.nzs, &s->dir, &nthreads)) {
        free(s);
        return NULL;
    }

    pows = (long *) malloc(sizeof(long) * s->mlm.nzs);
    for (unsigned long i = 0; i < s->mlm.nzs; ++i) {
        pows[i] = 1L;
    }

    s->thpool = thpool_init(nthreads);
    // (void) omp_set_num_threads(nthreads);

    (void) clt_mlm_setup(&s->mlm, s->dir, pows, kappa, size, g_verbose);

    /* Convert g_i values to python objects */
    {
        PyObject *py_gs, *py_state;

        py_gs = PyList_New(s->mlm.secparam);

        //
        // Only convert the first secparam g_i values since we only need to fill
        // in the first secparam slots of the plaintext space.
        //
        for (unsigned long i = 0; i < s->mlm.secparam; ++i) {
            PyList_SetItem(py_gs, i, mpz_to_py(s->mlm.gs[i]));
        }

        /* Encapsulate state as python object */
        py_state = PyCapsule_New((void *) s, NULL, state_destructor);

        return PyTuple_Pack(2, py_state, py_gs);
    }
}

//
// Encode N vectors across all slots of the MLM
//
static PyObject *
obf_encode_vectors(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_vectors, *py_list;
    char *name;
    mpz_t *vector;
    ssize_t length;
    double start;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OOOs", &py_state, &py_vectors, &py_list, &name))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    start = current_time();

    // We assume that all vectors have the same length, and thus just grab the
    // length of the first vector
    length = PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0));
    vector = (mpz_t *) calloc(length, sizeof(mpz_t));
    for (ssize_t i = 0; i < length; ++i) {
        mpz_init(vector[i]);
    }

    {
        struct write_vector_s *wv_s;

        wv_s = (struct write_vector_s *) malloc(sizeof(write_vector_s));
        wv_s->dir = (char *) calloc(strlen(s->dir) + 1, sizeof(char));
        (void) strcpy(wv_s->dir, s->dir);
        wv_s->name = (char *) calloc(strlen(name) + 1, sizeof(char));
        (void) strcpy(wv_s->name, name);
        wv_s->vector = vector;
        wv_s->length = length;
        wv_s->start = start;

        (void) thpool_add_class(s->thpool, name, length, thpool_write_vector,
                                wv_s);
        
        for (ssize_t i = 0; i < length; ++i) {
            mpz_t *elems;
            struct mlm_encode_elem_s *args;

            elems = (mpz_t *) calloc(s->mlm.secparam, sizeof(mpz_t));
            for (unsigned long j = 0; j < s->mlm.secparam; ++j) {
                mpz_init(elems[j]);
                py_to_mpz(elems[j],
                          PyList_GET_ITEM(PyList_GET_ITEM(py_vectors, j), i));
            }
            
            args = (struct mlm_encode_elem_s *)
                malloc(sizeof(struct mlm_encode_elem_s));
            args->out = &vector[i];
            args->elems = elems;
            args->mlm = &s->mlm;
            args->indices = (int *) calloc(2, sizeof(int));
            (void) extract_indices(py_list, &args->indices[0], &args->indices[1]);
            args->pows = (int *) calloc(2, sizeof(int));
            args->pows[0] = 1;
            args->pows[1] = 1;

            thpool_add_work(s->thpool, thpool_encode_elem, (void *) args, name);
        }
    }

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
    PyObject *py_state;
    int zero_indices[2], one_indices[2];
    long inp, idx, nrows, ncols;
    mpz_t *zero, *one;
    char idx_s[10];
    double start;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OllllOOOO", &py_state, &idx, &nrows, &ncols,
                          &inp, &py_zero_ms, &py_one_ms, &py_zero_set, &py_one_set))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    start = current_time();

    (void) snprintf(idx_s, 10, "%ld", idx);

    (void) extract_indices(py_zero_set, &zero_indices[0], &zero_indices[1]);
    (void) extract_indices(py_one_set, &one_indices[0], &one_indices[1]);

    zero = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
    one = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
    for (ssize_t i = 0; i < nrows * ncols; ++i) {
        mpz_inits(zero[i], one[i], NULL);
    }

    {
        struct write_layer_s *wl_s;

        wl_s = (struct write_layer_s *) malloc(sizeof(write_layer_s));
        wl_s->dir = (char *) calloc(strlen(s->dir) + 1, sizeof(char));
        (void) strcpy(wl_s->dir, s->dir);
        wl_s->zero = zero;
        wl_s->one = one;
        wl_s->inp = inp;
        wl_s->idx = idx;
        wl_s->nrows = nrows;
        wl_s->ncols = ncols;
        wl_s->start = start;

        (void) thpool_add_class(s->thpool, idx_s, 2 * nrows * ncols,
                                thpool_write_layer, wl_s);

        for (Py_ssize_t ctr = 0; ctr < 2 * nrows * ncols; ++ctr) {
            mpz_t *elems;
            PyObject *py_array;
            int *indices;
            mpz_t *val;
            size_t i;
            struct mlm_encode_elem_s *args;

            if (ctr < nrows * ncols) {
                i = ctr;
                val = &zero[i];
                py_array = py_zero_ms;
                indices = (int *) &zero_indices;
            } else {
                i = ctr - nrows * ncols;
                val = &one[i];
                py_array = py_one_ms;
                indices = (int *) &one_indices;
            }

            elems = (mpz_t *) malloc(sizeof(mpz_t) * s->mlm.secparam);
            for (unsigned long j = 0; j < s->mlm.secparam; ++j) {
                mpz_init(elems[j]);
                py_to_mpz(elems[j],
                          PyList_GET_ITEM(PyList_GET_ITEM(py_array, j), i));
            }

            args = (struct mlm_encode_elem_s *)
                malloc(sizeof(struct mlm_encode_elem_s));
            args->out = val;
            args->elems = elems;
            args->mlm = &s->mlm;
            args->indices = (int *) calloc(2, sizeof(int));
            args->indices[0] = indices[0];
            args->indices[1] = indices[1];
            args->pows = (int *) calloc(2, sizeof(int));
            args->pows[0] = 1;
            args->pows[1] = 1;

            thpool_add_work(s->thpool, thpool_encode_elem, (void *) args, idx_s);
        }
    }

    Py_RETURN_NONE;
}

static PyObject *
obf_sz_evaluate(PyObject *self, PyObject *args)
{
    char *dir = NULL;
    char *input = NULL;
    char *fname = NULL;
    int fnamelen;
    int iszero = -1;

    mpz_t tmp, q;
    mpz_t *result = NULL;
    long bplen, nrows, ncols = -1, nrows_prev = -1, nthreads;
    int err = 0;
    double start, end;

    if (!PyArg_ParseTuple(args, "ssll", &dir, &input, &bplen, &nthreads))
        return NULL;

    fnamelen = strlen(dir) + sizeof bplen + 7;

    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(tmp, q, NULL);

    // Load q
    (void) snprintf(fname, fnamelen, "%s/q", dir);
    (void) load_mpz_scalar(fname, q);

    (void) omp_set_num_threads(nthreads);

    for (int layer = 0; layer < bplen; ++layer) {
        unsigned int input_idx;
        mpz_t *left, *right;

        start = current_time();

        // determine the size of the matrix
        (void) snprintf(fname, fnamelen, "%s/%d.nrows", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        nrows = mpz_get_ui(tmp);
        (void) snprintf(fname, fnamelen, "%s/%d.ncols", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        ncols = mpz_get_ui(tmp);

        // find out the input bit for the given layer
        (void) snprintf(fname, fnamelen, "%s/%d.input", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        input_idx = mpz_get_ui(tmp);
        if (input_idx >= strlen(input)) {
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
            (void) snprintf(fname, fnamelen, "%s/%d.zero", dir, layer);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", dir, layer);
        }

        if (layer == 0) {
            result = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
            for (int i = 0; i < nrows * ncols; ++i) {
                mpz_init(result[i]);
            }
            (void) load_mpz_vector(fname, result, nrows * ncols);
            mpz_set(tmp, result[0]);
            nrows_prev = nrows;
        } else {
            left = result;
            right = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
            for (int i = 0; i < nrows * ncols; ++i) {
                mpz_init(right[i]);
            }
            (void) load_mpz_vector(fname, right, nrows * ncols);
            result = (mpz_t *) malloc(sizeof(mpz_t) * nrows_prev * ncols);
            for (int i = 0; i < nrows_prev * ncols; ++i) {
                mpz_init(result[i]);
            }
            mult_mats(result, left, right, q, nrows_prev, nrows, ncols);
            for (int i = 0; i < nrows_prev * nrows; ++i) {
                mpz_clear(left[i]);
            }
            for (int i = 0; i < nrows * ncols; ++i) {
                mpz_clear(right[i]);
            }
            free(left);
            free(right);
        }
        end = current_time();

        if (g_verbose)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n",
                           end - start);
    }

    if (!err) {
        mpz_t pzt, nu;

        start = current_time();
        mpz_inits(pzt, nu, NULL);
        (void) snprintf(fname, fnamelen, "%s/pzt", dir);
        (void) load_mpz_scalar(fname, pzt);
        (void) snprintf(fname, fnamelen, "%s/nu", dir);
        (void) load_mpz_scalar(fname, nu);
        iszero = clt_mlm_is_zero(result[1], pzt, q, mpz_get_ui(nu));
        mpz_clears(pzt, nu, NULL);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Zero test: %f\n", end - start);
    }

    for (int i = 0; i < nrows_prev * ncols; ++i) {
        mpz_clear(result[i]);
    }
    free(result);

    mpz_clears(tmp, q, NULL);

    if (fname)
        free(fname);

    if (err)
        return NULL;
    else
        return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyObject *
obf_evaluate(PyObject *self, PyObject *args)
{
    char *dir = NULL;
    char *input = NULL;
    char *fname = NULL;
    int fnamelen;
    int iszero = -1;
    mpz_t *comp, *s, *t;
    mpz_t tmp, q;
    long bplen, size, nthreads;
    int err = 0;
    double start, end;

    if (!PyArg_ParseTuple(args, "ssll", &dir, &input, &bplen, &nthreads))
        return NULL;
    fnamelen = strlen(dir) + sizeof bplen + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(tmp, q, NULL);

    // Get the size of the matrices
    (void) snprintf(fname, fnamelen, "%s/size", dir);
    (void) load_mpz_scalar(fname, tmp);
    size = mpz_get_ui(tmp);

    // Load q
    (void) snprintf(fname, fnamelen, "%s/q", dir);
    (void) load_mpz_scalar(fname, q);

    comp = (mpz_t *) malloc(sizeof(mpz_t) * size * size);
    s = (mpz_t *) malloc(sizeof(mpz_t) * size);
    t = (mpz_t *) malloc(sizeof(mpz_t) * size);
    if (!comp || !s || !t) {
        err = 1;
        goto cleanup;
    }
    for (int i = 0; i < size; ++i) {
        mpz_inits(s[i], t[i], NULL);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_init(comp[i]);
    }

    (void) omp_set_num_threads(nthreads);

    for (int layer = 0; layer < bplen; ++layer) {
        unsigned int input_idx;

        start = current_time();
        // find out the input bit for the given layer
        (void) snprintf(fname, fnamelen, "%s/%d.input", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        input_idx = mpz_get_ui(tmp);
        if (input_idx >= strlen(input)) {
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
            (void) snprintf(fname, fnamelen, "%s/%d.zero", dir, layer);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", dir, layer);
        }
        (void) load_mpz_vector(fname, comp, size * size);

        // for the first matrix, multiply 'comp' by 's' to get a vector
        if (layer == 0) {
            (void) snprintf(fname, fnamelen, "%s/s_enc", dir);
            (void) load_mpz_vector(fname, s, size);
        }
        mult_vect_by_mat(s, comp, q, size, t);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Multiplying matrices: %f\n",
                           end - start);
    }

    if (!err) {
        start = current_time();
        (void) snprintf(fname, fnamelen, "%s/t_enc", dir);
        (void) load_mpz_vector(fname, t, size);
        mult_vect_by_vect(tmp, s, t, q, size);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Multiplying vectors: %f\n",
                           end - start);

        start = current_time();
        {
            mpz_t pzt, nu;
            mpz_inits(pzt, nu, NULL);
            (void) snprintf(fname, fnamelen, "%s/pzt", dir);
            (void) load_mpz_scalar(fname, pzt);
            (void) snprintf(fname, fnamelen, "%s/nu", dir);
            (void) load_mpz_scalar(fname, nu);
            iszero = clt_mlm_is_zero(tmp, pzt, q, mpz_get_ui(nu));
            mpz_clears(pzt, nu, NULL);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Zero test: %f\n", end - start);
    }
    for (int i = 0; i < size; ++i) {
        mpz_clears(s[i], t[i], NULL);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_clear(comp[i]);
    }
cleanup:
    mpz_clears(tmp, q, NULL);
    if (comp)
        free(comp);
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
    {"encode_vectors", obf_encode_vectors, METH_VARARGS,
     "Encode a vector in each slot."},
    {"encode_layers", obf_encode_layers, METH_VARARGS,
     "Encode a branching program layer in each slot."},
    {"max_mem_usage", obf_max_mem_usage, METH_VARARGS,
     "Print out the maximum memory usage."},
    {"evaluate", obf_evaluate, METH_VARARGS,
     "Evaluate the obfuscation."},
    {"sz_evaluate", obf_sz_evaluate, METH_VARARGS,
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
