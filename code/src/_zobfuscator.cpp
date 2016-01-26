#include "circuit.h"
#include "utils.h"
#include "pyutils.h"
#include "clt_mlm.h"
#include "thpool.h"
#include "thpool_fns.h"

#include <omp.h>

struct state {
    threadpool thpool;
    struct clt_mlm_state mlm;
    char *dir;
    mpz_t nchk;
    mpz_t nev;
};

static void
state_destructor(PyObject *self)
{
    struct state *s;

    s = (struct state *) PyCapsule_GetPointer(self, NULL);
    if (s) {
        clt_mlm_cleanup(&s->mlm);
        // thpool_destroy(s->thpool);
        mpz_clears(s->nev, s->nchk, NULL);
    }
}

static int
write_element(const char *dir, mpz_t elem, const char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, name);
    (void) save_mpz_scalar(fname, elem);
    free(fname);
    return 0;
}

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long kappa, nthreads, ncores;
    PyObject *py_pows;
    int *pows;
    struct state *s;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    if (!PyArg_ParseTuple(args, "lllOsll", &s->mlm.secparam, &kappa, &s->mlm.nzs,
                          &py_pows, &s->dir, &nthreads, &ncores)) {
        free(s);
        return NULL;
    }

    pows = (int *) malloc(sizeof(int) * s->mlm.nzs);
    for (size_t i = 0; i < s->mlm.nzs; ++i) {
        pows[i] = (int)PyLong_AsLong(PyList_GET_ITEM(py_pows, i));
    }

    s->thpool = thpool_init(nthreads);
    (void) omp_set_num_threads(ncores);

    if (g_verbose) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
    }

    (void) clt_mlm_setup(&s->mlm, s->dir, pows, kappa, 0, g_verbose);
    mpz_init_set(s->nev, s->mlm.gs[0]);
    mpz_init_set(s->nchk, s->mlm.gs[1]);
    
    {
        PyObject *py_state;
        py_state = PyCapsule_New((void *) s, NULL, state_destructor);
        if (py_state == NULL)
            return NULL;
        return py_state;
    }
}

static void
create_work(struct state *s, const char *str, int i,
            const mpz_t elem0, const mpz_t elem1, unsigned int num, ...)
{
    va_list vals;

    struct write_element_s *we_s;
    struct mlm_encode_elem_s *args;
    mpz_t *out, *elems;

    out = (mpz_t *) malloc(sizeof(mpz_t));
    elems = (mpz_t *) calloc(2, sizeof(mpz_t));
    mpz_inits(*out, elems[0], elems[1], NULL);    
    mpz_set(elems[0], elem0);
    mpz_set(elems[1], elem1);

    we_s = (struct write_element_s *) malloc(sizeof(write_element_s));
    we_s->dir = (char *) calloc(strlen(s->dir) + 1, sizeof(char));
    (void) strcpy(we_s->dir, s->dir);
    we_s->elem = out;
    we_s->name = (char *) calloc(sizeof(int) + 10, sizeof(char));
    if (i >= 0) {
        (void) sprintf(we_s->name, str, i);
    } else {
        (void) strcpy(we_s->name, str);
    }
    (void) thpool_add_tag(s->thpool, we_s->name, 1, thpool_write_element, we_s);

    args = (struct mlm_encode_elem_s *)
        malloc(sizeof(struct mlm_encode_elem_s));
    args->mlm = &s->mlm;
    args->out = out;
    args->nins = 2;
    args->ins = elems;
    args->nzs = num;
    args->indices = (int *) calloc(num, sizeof(int));
    args->pows = (int *) calloc(num, sizeof(int));

    va_start(vals, num);
    for (unsigned int i = 0; i < num; ++i) {
        args->indices[i] = va_arg(vals, int);
        args->pows[i] = va_arg(vals, int);
    }

    (void) thpool_add_work(s->thpool, thpool_encode_elem, (void *) args,
                           we_s->name);
}

static PyObject *
obf_encode_circuit(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_ys, *py_xdegs;
    mpz_t tmp, c_star;
    mpz_t *alphas, *betas;
    int n, m, ydeg;
    char *circuit;
    mpz_t zero, one;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OsOOiii", &py_state, &circuit, &py_ys,
                          &py_xdegs, &ydeg, &n, &m))
        return NULL;
    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    mpz_inits(c_star, tmp, zero, one, NULL);
    mpz_set_ui(zero, 0);
    mpz_set_ui(one, 1);

    alphas = (mpz_t *) malloc(sizeof(mpz_t) * n);
    for (int i = 0; i < n; ++i) {
        mpz_init(alphas[i]);
        mpz_urandomm(alphas[i], s->mlm.rng, s->nchk);
    }
    betas = (mpz_t *) malloc(sizeof(mpz_t) * m);
    for (int i = 0; i < m; ++i) {
        mpz_init(betas[i]);
        mpz_urandomm(betas[i], s->mlm.rng, s->nchk);
    }

    for (int i = 0; i < n; ++i) {
        mpz_t elems[2];
        int deg;

        mpz_inits(elems[0], elems[1], NULL);

        deg = PyLong_AsLong(PyList_GET_ITEM(py_xdegs, i));

        create_work(s, "x_\%d_0", i, zero, alphas[i],
                    1, 2 * i, 1);
        create_work(s, "x_\%d_1", i, one, alphas[i],
                    1, 2 * i + 1, 1);
        create_work(s, "u_\%d_0", i, one, one,
                    1, 2 * i, 1);
        create_work(s, "u_\%d_1", i, one, one,
                    1, 2 * i + 1, 1);

        mpz_urandomm(elems[0], s->mlm.rng, s->nev);
        mpz_urandomm(elems[1], s->mlm.rng, s->nchk);

        create_work(s, "z_\%d_0", i, elems[0], elems[1],
                    3, 2 * i + 1, deg, 2 * n + i, 1, 3 * n + i, 1);
        create_work(s, "w_\%d_0", i, zero, elems[1],
                    1, 3 * n + i, 1);

        mpz_urandomm(elems[0], s->mlm.rng, s->nev);
        mpz_urandomm(elems[1], s->mlm.rng, s->nchk);

        create_work(s, "z_\%d_1", i, elems[0], elems[1],
                    3, 2 * i, deg, 2 * n + i, 1, 3 * n + i, 1);
        create_work(s, "w_\%d_1", i, zero, elems[1], 1, 3 * n + i, 1);

        mpz_clears(elems[0], elems[1], NULL);
    }

    for (int i = 0; i < m; ++i) {
        mpz_t elems[2];
        mpz_inits(elems[0], elems[1], NULL);

        py_to_mpz(elems[0], PyList_GET_ITEM(py_ys, i));
        if (mpz_sgn(elems[0]) == -1) {
            mpz_mod(elems[0], elems[0], s->nev);
        }
        mpz_set(elems[1], betas[i]);
        create_work(s, "y_\%d", i, elems[0], elems[1], 1, 4 * n, 1);

        mpz_clears(elems[0], elems[1], NULL);
    }

    create_work(s, "v", -1, one, one, 1, 4 * n, 1);

    {
        struct circuit *c;

        c = circ_parse(circuit);
        {
            int circnamelen;
            char *circname;
            circnamelen = strlen(s->dir) + strlen("/circuit") + 2;
            circname = (char *) malloc(sizeof(char) * circnamelen);
            (void) snprintf(circname, circnamelen, "%s/circuit", s->dir);
            (void) circ_copy_circuit(circuit, circname);
            free(circname);
        }
        (void) circ_evaluate(c, alphas, betas, c_star, s->mlm.q);
        circ_cleanup(c);
    }

    // Last element to encode, so don't need to parallelize
    {
        mpz_t elems[2];
        int idx_set_size;
        int *indices, *pows;

        mpz_init_set_ui(elems[0], 0);
        mpz_init_set(elems[1], c_star);

        // The index set is laid out as follows:
        //   - The first 2 * n entries contain X_i,0 and X_i,1
        //   - The next n entries contain Z_i
        //   - The next n entries contain W_i
        //   - The final entry contains Y
        idx_set_size = 4 * n + 1;

        indices = (int *) malloc(sizeof(int) * idx_set_size);
        pows = (int *) malloc(sizeof(int) * idx_set_size);
        
        // The C* encoding contains everything but the W_i symbols.
        // Here we calculate the appropriate indices and powers.
        for (int i = 0; i < n; ++i) {
            int deg;

            deg = PyLong_AsLong(PyList_GET_ITEM(py_xdegs, i));
            // X_i,0^deg(x_i)
            indices[2 * i] = 2 * i;
            pows[2 * i] = deg;
            // X_i,1^deg(x_i)
            indices[2 * i + 1] = 2 * i + 1;
            pows[2 * i + 1] = deg;
            // Z_i
            indices[2 * n + i] = 2 * n + i;
            pows[2 * n + i] = 1;
        }
        // Y^deg(y)
        indices[3 * n] = 4 * n;
        pows[3 * n] = ydeg;
        // Encode against these indices/powers
        clt_mlm_encode(&s->mlm, tmp, 2, elems, 3 * n + 1, indices, pows);

        mpz_clears(elems[0], elems[1], NULL);
        free(indices);
        free(pows);
    }
    (void) write_element(s->dir, tmp, "c_star");

    mpz_clears(c_star, tmp, NULL);

    thpool_wait(s->thpool);
    
    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate(PyObject *self, PyObject *args)
{
    char *circuit, *dir, *input, *fname;
    long n, m, nthreads;
    int fnamelen;
    int iszero;
    mpz_t c_1, c_2, q, z, w;
    mpz_t *xs, *xones, *ys, *yones;

    if (!PyArg_ParseTuple(args, "ssslll", &dir, &circuit, &input, &n, &m,
                          &nthreads))
        return NULL;
    fnamelen = strlen(dir) + sizeof(int) + 5;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(c_1, c_2, q, z, w, NULL);

    xs = (mpz_t *) malloc(sizeof(mpz_t) * n);
    xones = (mpz_t *) malloc(sizeof(mpz_t) * n);
    for (int i = 0; i < n; ++i) {
        mpz_inits(xs[i], xones[i], NULL);
    }
    ys = (mpz_t *) malloc(sizeof(mpz_t) * m);
    yones = (mpz_t *) malloc(sizeof(mpz_t) * m);
    for (int i = 0; i < m; ++i) {
        mpz_inits(ys[i], yones[i], NULL);
    }

    // Load q
    (void) snprintf(fname, fnamelen, "%s/q", dir);
    (void) load_mpz_scalar(fname, q);

    // Check that all input choices are bits
    for (int i = 0; i < n; ++i) {
        if (input[i] != '0' && input[i] != '1') {
            PyErr_SetString(PyExc_RuntimeError, "input must be 0 or 1");
            return NULL;
        }
    }

    (void) omp_set_num_threads(nthreads);

    // Load in appropriate input
    for (int i = 0; i < n; ++i) {
        (void) snprintf(fname, fnamelen, "%s/x_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, xs[i]);
        (void) snprintf(fname, fnamelen, "%s/u_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, xones[i]);
    }

    // Load in secret input
    for (int i = 0; i < m; ++i) {
        (void) snprintf(fname, fnamelen, "%s/y_%d", dir, i);
        (void) load_mpz_scalar(fname, ys[i]);
        (void) snprintf(fname, fnamelen, "%s/v", dir);
        (void) load_mpz_scalar(fname, yones[i]);
    }

    // Evaluate the circuit on x_1, ..., x_n
    {
        struct circuit *c;

        c = circ_parse(circuit);
        (void) circ_evaluate_encoding(c, xs, xones, ys, yones, c_1, q);
        circ_cleanup(c);
    }

    // Load in c_2
    (void) snprintf(fname, fnamelen, "%s/c_star", dir);
    (void) load_mpz_scalar(fname, c_2);

    // Compute c_1 * \Prod z_{i,x_i} and c_2 * \Prod w_{i,x_i}
    for (int i = 0; i < n; ++i) {
        (void) snprintf(fname, fnamelen, "%s/z_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, z);
        mpz_mul(c_1, c_1, z);
        mpz_mod(c_1, c_1, q);
        (void) snprintf(fname, fnamelen, "%s/w_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, w);
        mpz_mul(c_2, c_2, w);
        mpz_mod(c_2, c_2, q);
    }

    // Compute c_1 - c_2 and zero test
    {
        mpz_t pzt, nu, tmp;
        mpz_inits(tmp, pzt, nu, NULL);
        (void) snprintf(fname, fnamelen, "%s/pzt", dir);
        (void) load_mpz_scalar(fname, pzt);
        (void) snprintf(fname, fnamelen, "%s/nu", dir);
        (void) load_mpz_scalar(fname, nu);
        mpz_sub(tmp, c_1, c_2);
        iszero = clt_mlm_is_zero(tmp, pzt, q, mpz_get_ui(nu));
        mpz_clears(tmp, pzt, nu, NULL);
    }

    free(fname);
    for (int i = 0; i < m; ++i) {
        mpz_clears(ys[i], yones[i], NULL);
    }
    free(ys);
    free(yones);
    for (int i = 0; i < n; ++i) {
        mpz_clears(xs[i], xones[i], NULL);
    }
    free(xs);
    free(xones);
    mpz_clears(c_1, c_2, q, z, w, NULL);

    return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyMethodDef
ObfMethods[] = {
    {"verbose", obf_verbose, METH_VARARGS,
     "Set verbosity."},
    {"setup", obf_setup, METH_VARARGS,
     "Set up obfuscator."},
    {"encode_circuit", obf_encode_circuit, METH_VARARGS,
     "Encode circuit."},
    {"evaluate", obf_evaluate, METH_VARARGS,
     "Evaluate circuit."},
    {"max_mem_usage", obf_max_mem_usage, METH_VARARGS,
     "Compute the maximum memory usage."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_zobfuscator(void)
{
    (void) Py_InitModule("_zobfuscator", ObfMethods);
}
