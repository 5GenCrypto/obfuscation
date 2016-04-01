#include "utils.h"
#include "pyutils.h"
#include "thpool.h"
#include "thpool_fns.h"
#include <clt13.h>
#include <omp.h>

struct state {
    threadpool thpool;
    unsigned long secparam;
    clt_state mlm;
    char *dir;
};

void
state_destructor(PyObject *self)
{
    struct state *s;

    s = (struct state *) PyCapsule_GetPointer(self, NULL);
    if (s) {
        clt_state_clear(&s->mlm);
        thpool_destroy(s->thpool);
    }
    free(s);
}

//
//
// Python functions
//
//

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long kappa, size, nzs, nthreads, ncores;
    struct state *s = NULL;
    int *pows = NULL;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    if (!PyArg_ParseTuple(args, "llllsll", &s->secparam, &kappa, &size,
                          &nzs, &s->dir, &nthreads, &ncores)) {
        free(s);
        return NULL;
    }

    if (kappa <= 0 || size < 0 || nzs <= 0) {
        Py_RETURN_NONE;
    }

    pows = (int *) calloc(nzs, sizeof(int));
    for (long i = 0; i < nzs; ++i) {
        pows[i] = 1;
    }

    s->thpool = thpool_init(nthreads);
    (void) omp_set_num_threads(ncores);

    if (g_verbose) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
    }

    clt_state_init(&s->mlm, kappa, s->secparam, nzs, pows);
    // Write public parameters to disk
    {
        clt_pp pp;
        clt_pp_init(&pp, &s->mlm);
        clt_pp_save(&pp, s->dir);
        clt_pp_clear(&pp);
        // Needed for AGIS obfuscator
        if (size > 0) {
            char *fname;
            mpz_t tmp;
            int len;

            mpz_init(tmp);
            mpz_set_ui(tmp, size);
            len = strlen(s->dir) + 6;
            fname = (char *) calloc(len, sizeof(char));
            (void) snprintf(fname, len, "%s/size", s->dir);
            (void) save_mpz_scalar(fname, tmp);
            free(fname);
            mpz_clear(tmp);
        }
    }

    /* Convert g_i values to python objects */
    {
        PyObject *py_gs, *py_state;

        py_gs = PyList_New(s->secparam);

        //
        // Only convert the first secparam g_i values since we only need to fill
        // in the first secparam slots of the plaintext space.
        //
        for (unsigned long i = 0; i < s->secparam; ++i) {
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
    PyObject *py_state, *py_vectors;
    long index;
    char *name;
    mpz_t *vector;
    ssize_t length;
    double start;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OOls", &py_state, &py_vectors, &index, &name))
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

        (void) thpool_add_tag(s->thpool, name, length, thpool_write_vector,
                              wv_s);
        
        for (ssize_t i = 0; i < length; ++i) {
            mpz_t *elems;
            struct mlm_encode_elem_s *args;

            elems = (mpz_t *) calloc(s->secparam, sizeof(mpz_t));
            for (unsigned long j = 0; j < s->secparam; ++j) {
                mpz_init(elems[j]);
                py_to_mpz(elems[j],
                          PyList_GET_ITEM(PyList_GET_ITEM(py_vectors, j), i));
            }
            
            args = (struct mlm_encode_elem_s *)
                malloc(sizeof(struct mlm_encode_elem_s));
            args->out = &vector[i];
            args->mlm = &s->mlm;
            args->nins = s->mlm.nzs;
            args->ins = elems;
            args->pows = (int *) calloc(args->nins, sizeof(int));
            args->pows[index] = 1;

            for (unsigned long i = 0; i < args->nins; ++i) {
                printf("%d ", args->pows[i]);
            }
            printf("\n");

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
    PyObject *py_state;
    long zero_index, one_index;
    long inp, idx, nrows, ncols;
    mpz_t *zero, *one;
    char idx_s[10];
    double start;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OllllOOll", &py_state, &idx, &nrows, &ncols,
                          &inp, &py_zero_ms, &py_one_ms, &zero_index, &one_index))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    start = current_time();

    (void) snprintf(idx_s, 10, "%ld", idx);

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

        (void) thpool_add_tag(s->thpool, idx_s, 2 * nrows * ncols,
                              thpool_write_layer, wl_s);

        for (Py_ssize_t ctr = 0; ctr < 2 * nrows * ncols; ++ctr) {
            mpz_t *elems;
            PyObject *py_array;
            int index;
            mpz_t *val;
            size_t i;
            struct mlm_encode_elem_s *args;

            if (ctr < nrows * ncols) {
                i = ctr;
                val = &zero[i];
                py_array = py_zero_ms;
                index = zero_index;
            } else {
                i = ctr - nrows * ncols;
                val = &one[i];
                py_array = py_one_ms;
                index = one_index;
            }

            elems = (mpz_t *) malloc(sizeof(mpz_t) * s->secparam);
            for (unsigned long j = 0; j < s->secparam; ++j) {
                mpz_init(elems[j]);
                py_to_mpz(elems[j],
                          PyList_GET_ITEM(PyList_GET_ITEM(py_array, j), i));
            }

            args = (struct mlm_encode_elem_s *)
                malloc(sizeof(struct mlm_encode_elem_s));
            args->out = val;
            args->mlm = &s->mlm;
            args->nins = s->mlm.nzs;
            args->ins = elems;
            args->pows = (int *) calloc(args->nins, sizeof(int));
            args->pows[index] = 1;

            for (unsigned long i = 0; i < args->nins; ++i) {
                printf("%d ", args->pows[i]);
            }
            printf("\n");

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
    clt_pp pp;
    int iszero = -1;

    mpz_t tmp;
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

    mpz_inits(tmp, NULL);

    clt_pp_read(&pp, dir);

    (void) omp_set_num_threads(nthreads);

    int nmults = 0;

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
            mult_mats(result, left, right, pp.x0, nrows_prev, nrows, ncols);
            nmults++;
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

    printf("%d\n", nmults);

    if (!err) {
        start = current_time();
        iszero = clt_is_zero(&pp, result[1]);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Zero test: %f\n", end - start);
    }

    for (int i = 0; i < nrows_prev * ncols; ++i) {
        mpz_clear(result[i]);
    }
    free(result);

    clt_pp_clear(&pp);
    mpz_clears(tmp, NULL);

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
    clt_pp pp;
    mpz_t *comp, *s, *t;
    mpz_t tmp;
    long bplen, size, nthreads;
    int err = 0;
    double start, end;

    if (!PyArg_ParseTuple(args, "ssll", &dir, &input, &bplen, &nthreads))
        return NULL;
    fnamelen = strlen(dir) + sizeof bplen + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(tmp, NULL);
    clt_pp_read(&pp, dir);

    // Get the size of the matrices
    (void) snprintf(fname, fnamelen, "%s/size", dir);
    (void) load_mpz_scalar(fname, tmp);
    size = mpz_get_ui(tmp);

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
        mult_vect_by_mat(s, comp, pp.x0, size, t);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Multiplying matrices: %f\n",
                           end - start);
    }

    if (!err) {
        start = current_time();
        (void) snprintf(fname, fnamelen, "%s/t_enc", dir);
        (void) load_mpz_vector(fname, t, size);
        mult_vect_by_vect(tmp, s, t, pp.x0, size);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Multiplying vectors: %f\n",
                           end - start);

        start = current_time();
        {
            iszero = clt_is_zero(&pp, tmp);
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
    clt_pp_clear(&pp);
    mpz_clears(tmp, NULL);
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
