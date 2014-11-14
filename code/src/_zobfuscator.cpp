#include "circuit.h"
#include "utils.h"
#include "pyutils.h"
#include "clt_mlm.h"

struct state {
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

static void
set_indices_pows(int *indices, int *pows, unsigned int num, ...)
{
    va_list elems;

    va_start(elems, num);
    for (unsigned int i = 0; i < num; ++i) {
        indices[i] = va_arg(elems, int);
        pows[i] = va_arg(elems, int);
    }
}

static void
encode(struct clt_mlm_state *s, mpz_t out, mpz_t in1, mpz_t in2, unsigned int num,
       const int *indices, const int *pows)
{
    mpz_t r, tmp;

    mpz_inits(r, tmp, NULL);
    mpz_set_ui(out, 0);
    for (unsigned long i = 0; i < s->n; ++i) {
        mpz_genrandom(r, &s->rng, s->rho);
        mpz_mul(tmp, r, s->gs[i]);
        if (i == 0) {
            mpz_add(tmp, tmp, in1);
        } else if (i == 1) {
            mpz_add(tmp, tmp, in2);
        }
        mpz_mul(tmp, tmp, s->crt_coeffs[i]);
        mpz_add(out, out, tmp);
    }
    for (unsigned long i = 0; i < num; ++i) {
        int idx, pow;

        idx = indices[i];
        pow = pows[i];
        mpz_powm_ui(tmp, s->zinvs[idx], pow, s->q);
        mpz_mul(out, out, tmp);
        mpz_mod(out, out, s->q);
    }
    mpz_clears(r, tmp, NULL);
}

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long kappa;
    PyObject *py_pows;
    long *pows;
    struct state *s;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    if (!PyArg_ParseTuple(args, "lllOs", &s->mlm.secparam, &kappa, &s->mlm.nzs,
                          &py_pows, &s->dir)) {
        free(s);
        return NULL;
    }

    pows = (long *) malloc(sizeof(long) * s->mlm.nzs);
    for (size_t i = 0; i < s->mlm.nzs; ++i) {
        pows[i] = PyLong_AsLong(PyList_GET_ITEM(py_pows, i));
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

static PyObject *
obf_encode_circuit(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_ys, *py_xdegs;
    mpz_t zero, one, tmp, c_star;
    mpz_t *alphas, *betas;
    int n, m, ydeg;
    int *indices, *pows;
    char *circuit, *fname;
    int fnamelen = sizeof(int) + 5;
    int idx_set_size;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OsOOiii", &py_state, &circuit, &py_ys,
                          &py_xdegs, &ydeg, &n, &m))
        return NULL;
    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(c_star, tmp, NULL);
    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(one, 1);

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

    // The index set is laid out as follows:
    //   - The first 2 * n entries contain X_i,0 and X_i,1
    //   - The next n entries contain Z_i
    //   - The next n entries contain W_i
    //   - The final entry contains Y
    idx_set_size = 4 * n + 1;

    indices = (int *) malloc(sizeof(int) * idx_set_size);
    pows = (int *) malloc(sizeof(int) * idx_set_size);

    for (int i = 0; i < n; ++i) {
        mpz_t delta, gamma, out;
        int deg;

        mpz_inits(delta, gamma, out, NULL);

        deg = PyLong_AsLong(PyList_GET_ITEM(py_xdegs, i));

        set_indices_pows(indices, pows, 1, 2 * i, 1);
        encode(&s->mlm, out, zero, alphas[i], 1, indices, pows);
        (void) snprintf(fname, fnamelen, "x_%d_0", i);
        (void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 2 * i, 1);
        encode(&s->mlm, out, one, one, 1, indices, pows);
        (void) snprintf(fname, fnamelen, "u_%d_0", i);
        (void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 2 * i + 1, 1);
        encode(&s->mlm, out, one, alphas[i], 1, indices, pows);
        (void) snprintf(fname, fnamelen, "x_%d_1", i);
        (void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 2 * i + 1, 1);
        encode(&s->mlm, out, one, one, 1, indices, pows);
        (void) snprintf(fname, fnamelen, "u_%d_1", i);
        (void) write_element(s->dir, out, fname);

        mpz_urandomm(delta, s->mlm.rng, s->nev);
        mpz_urandomm(gamma, s->mlm.rng, s->nchk);

        set_indices_pows(indices, pows, 3, 2 * i + 1, deg, 2 * n + i, 1, 3 * n + i, 1);
        encode(&s->mlm, out, delta, gamma, 3, indices, pows);
        (void) snprintf(fname, fnamelen, "z_%d_0", i);
        (void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 3 * n + i, 1);
        encode(&s->mlm, out, zero, gamma, 1, indices, pows);
        (void) snprintf(fname, fnamelen, "w_%d_0", i);
        (void) write_element(s->dir, out, fname);

        mpz_urandomm(delta, s->mlm.rng, s->nev);
        mpz_urandomm(gamma, s->mlm.rng, s->nchk);

        set_indices_pows(indices, pows, 3, 2 * i, deg, 2 * n + i, 1, 3 * n + i, 1);
        encode(&s->mlm, out, delta, gamma, 3, indices, pows);
        (void) snprintf(fname, fnamelen, "z_%d_1", i);
        (void) write_element(s->dir, out, fname);

        set_indices_pows(indices, pows, 1, 3 * n + i, 1);
        encode(&s->mlm, out, zero, gamma, 1, indices, pows);
        (void) snprintf(fname, fnamelen, "w_%d_1", i);
        (void) write_element(s->dir, out, fname);

        mpz_clears(delta, gamma, out, NULL);
    }
    set_indices_pows(indices, pows, 1, 4 * n, 1);
    for (int i = 0; i < m; ++i) {
        mpz_t out, y;
        mpz_inits(out, y, NULL);

        py_to_mpz(y, PyList_GET_ITEM(py_ys, i));
        if (mpz_sgn(y) == -1) {
            mpz_mod(y, y, s->nev);
        }
        encode(&s->mlm, out, y, betas[i], 1, indices, pows);
        (void) snprintf(fname, fnamelen, "y_%d", i);
        (void) write_element(s->dir, out, fname);
        mpz_clears(out, y, NULL);
    }
    encode(&s->mlm, tmp, one, one, 1, indices, pows);
    (void) write_element(s->dir, tmp, "v");

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

    {
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
        encode(&s->mlm, tmp, zero, c_star, 3 * n + 1, indices, pows);
    }
    (void) write_element(s->dir, tmp, "c_star");

    mpz_clears(c_star, tmp, zero, one, NULL);
    
    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate(PyObject *self, PyObject *args)
{
    char *circuit, *dir, *input, *fname;
    long n, m;
    int fnamelen;
    int iszero;
    mpz_t c_1, c_2, q, z, w;
    mpz_t *xs, *xones, *ys, *yones;

    if (!PyArg_ParseTuple(args, "sssll", &dir, &circuit, &input, &n, &m))
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
        iszero = is_zero(tmp, pzt, q, mpz_get_ui(nu));
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
	mpz_clears(s->nev, s->nchk, NULL);
	free(s);

    Py_RETURN_NONE;
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
    {"cleanup", obf_cleanup, METH_VARARGS,
     "Clean up objects created during setup."},

    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_zobfuscator(void)
{
    (void) Py_InitModule("_zobfuscator", ObfMethods);
}
