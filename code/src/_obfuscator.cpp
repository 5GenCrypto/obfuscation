#include "utils.h"
#include "pyutils.h"
#include "clt_mlm.h"

#include <omp.h>

struct state {
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

static int
write_vector(const char *dir, mpz_t *vector, long size, char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, name);
    (void) save_mpz_vector(fname, vector, size);
    free(fname);
    return 0;
}

static int
write_layer(const char *dir, int inp, long idx, mpz_t *zero, mpz_t *one,
            long nrows, long ncols)
{
    mpz_t tmp;
    char *fname;
    int fnamelen;

    if (idx < 0)
        return 1;
    fnamelen = strlen(dir) + sizeof idx + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    mpz_init_set_ui(tmp, inp);
    (void) snprintf(fname, fnamelen, "%s/%ld.input", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    (void) snprintf(fname, fnamelen, "%s/%ld.zero", dir, idx);
    (void) save_mpz_vector(fname, zero, nrows * ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.one", dir, idx);
    (void) save_mpz_vector(fname, one, nrows * ncols);
    mpz_set_ui(tmp, nrows);
    (void) snprintf(fname, fnamelen, "%s/%ld.nrows", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    mpz_set_ui(tmp, ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.ncols", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    free(fname);
    mpz_clear(tmp);

    return 0;
}

//
//
// Python functions
//
//

PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long kappa, size;
    struct state *s = NULL;
    long *pows = NULL;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    if (!PyArg_ParseTuple(args, "lllls", &s->mlm.secparam, &kappa, &size,
                          &s->mlm.nzs, &s->dir)) {
        free(s);
        return NULL;
    }

    pows = (long *) malloc(sizeof(long) * s->mlm.nzs);
    for (unsigned long i = 0; i < s->mlm.nzs; ++i) {
        pows[i] = 1L;
    }

    (void) clt_mlm_setup(&s->mlm, s->dir, pows, kappa, size, g_verbose);

    {
        PyObject *py_gs, *py_state;

        py_gs = PyList_New(s->mlm.secparam);

        /* Convert g_i values to python objects */
        //
        // Only convert the first secparam g_i values since we only need to fill in
        // the first secparam slots of the plaintext space.
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
    int indices[2];
    int pows[] = {1, 1};
    mpz_t *vector;
    ssize_t length;
    double start, end;
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
    vector = (mpz_t *) malloc(sizeof(mpz_t) * length);

    (void) extract_indices(py_list, &indices[0], &indices[1]);

#pragma omp parallel for
    for (ssize_t i = 0; i < length; ++i) {
        mpz_t *elems;
        mpz_init(vector[i]);
        elems = (mpz_t *) malloc(sizeof(mpz_t) * s->mlm.secparam);
        for (unsigned long j = 0; j < s->mlm.secparam; ++j) {
            mpz_init(elems[j]);
            py_to_mpz(elems[j],
                      PyList_GET_ITEM(PyList_GET_ITEM(py_vectors, j), i));
        }
        clt_mlm_encode(&s->mlm, vector[i], s->mlm.secparam, elems, 2, indices, pows);
        for (unsigned long j = 0; j < s->mlm.secparam; ++j) {
            mpz_init(elems[j]);
        }
        free(elems);
    }
    (void) write_vector(s->dir, vector, length, name);
    for (ssize_t i = 0; i < length; ++i) {
        mpz_clear(vector[i]);
    }
    free(vector);

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       length, end - start);

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
    int pows[] = {1, 1};
    long inp, idx, nrows, ncols;
    mpz_t *zero, *one;
    double start, end;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OllllOOOO", &py_state, &idx, &nrows, &ncols,
                          &inp, &py_zero_ms, &py_one_ms, &py_zero_set, &py_one_set))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    start = current_time();

    (void) extract_indices(py_zero_set, &zero_indices[0], &zero_indices[1]);
    (void) extract_indices(py_one_set, &one_indices[0], &one_indices[1]);

    zero = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
    one = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);

#pragma omp parallel for
    for (Py_ssize_t ctr = 0; ctr < 2 * nrows * ncols; ++ctr) {
        mpz_t *elems;
        PyObject *py_array;
        int *indices;
        mpz_t *val;
        size_t i;

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

        mpz_init(*val);
        elems = (mpz_t *) malloc(sizeof(mpz_t) * s->mlm.secparam);
        for (unsigned long j = 0; j < s->mlm.secparam; ++j) {
            mpz_init(elems[j]);
            py_to_mpz(elems[j],
                      PyList_GET_ITEM(PyList_GET_ITEM(py_array, j), i));
        }
        clt_mlm_encode(&s->mlm, *val, s->mlm.secparam, elems, 2, indices, pows);
        for (unsigned long j = 0; j < s->mlm.secparam; ++j) {
            mpz_init(elems[j]);
        }
        free(elems);
    }
    (void) write_layer(s->dir, inp, idx, zero, one, nrows, ncols);

    for (int i = 0; i < nrows * ncols; ++i) {
        mpz_clears(zero[i], one[i], NULL);
    }
    free(zero);
    free(one);

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       2 * nrows * ncols, end - start);

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
    long bplen, nrows, ncols = -1, nrows_prev = -1;
    int err = 0;
    double start, end;

    if (!PyArg_ParseTuple(args, "ssl", &dir, &input, &bplen))
        return NULL;

    fnamelen = strlen(dir) + sizeof bplen + 7;

    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(tmp, q, NULL);

    // Load q
    (void) snprintf(fname, fnamelen, "%s/q", dir);
    (void) load_mpz_scalar(fname, q);

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
    long bplen, size;
    int err = 0;
    double start, end;

    if (!PyArg_ParseTuple(args, "ssl", &dir, &input, &bplen))
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
obf_cleanup(PyObject *self, PyObject *args)
{
    PyObject *py_state;
    struct state *s;

    if (!PyArg_ParseTuple(args, "O", &py_state))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    clt_mlm_cleanup(&s->mlm);
    free(s);

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
    {"cleanup", obf_cleanup, METH_VARARGS,
     "Clean up objects created during setup."},
    {"evaluate", obf_evaluate, METH_VARARGS,
     "evaluate the obfuscation."},
    {"sz_evaluate", obf_sz_evaluate, METH_VARARGS,
     "evaluate the obfuscation."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
