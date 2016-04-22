#include "circuit.h"
#include "utils.h"
#include "pyutils.h"
#include "thpool.h"
#include "thpool_fns.h"

#include <aesrand.h>
#include <mmap/mmap.h>
#include <mmap/mmap_clt.h>
#include <omp.h>

struct state {
    threadpool thpool;
    unsigned long secparam;
    clt_state mlm;
    aes_randstate_t rand;
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
        aes_randclear(s->rand);
        clt_state_clear(&s->mlm);
        // thpool_destroy(s->thpool);
        mpz_clears(s->nev, s->nchk, NULL);
    }
}

static int
write_element(const char *dir, mmap_enc *elem, const char *name)
{
    FILE *fp;
    char fname[100];

    (void) snprintf(fname, 100, "%s/%s", dir, name);
    fp = fopen(fname, "w+b");
    clt_vtable.enc->fwrite(elem, fp);
    fclose(fp);
    return 0;
}

static PyObject *
obf_setup(PyObject *self, PyObject *args)
{
    long kappa, nthreads, ncores;
    PyObject *py_pows;
    int *pows;
    struct state *s;
    int flags = CLT_FLAG_DEFAULT | CLT_FLAG_OPT_PARALLEL_ENCODE;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    if (!PyArg_ParseTuple(args, "lllOsll", &s->secparam, &kappa, &s->mlm.nzs,
                          &py_pows, &s->dir, &nthreads, &ncores)) {
        free(s);
        return NULL;
    }

    pows = (int *) malloc(sizeof(int) * s->mlm.nzs);
    for (size_t i = 0; i < s->mlm.nzs; ++i) {
        pows[i] = (int) PyLong_AsLong(PyList_GET_ITEM(py_pows, i));
    }

    (void) aes_randinit(s->rand);
    s->thpool = thpool_init(nthreads);
    (void) omp_set_num_threads(ncores);

    if (g_verbose) {
        fprintf(stderr, "  # Threads: %ld\n", nthreads);
        fprintf(stderr, "  # Cores: %ld\n", ncores);
        flags |= CLT_FLAG_VERBOSE;
    }

    clt_state_init(&s->mlm, kappa, s->secparam, s->mlm.nzs, pows, flags,
                   s->rand);
    {
        clt_pp pp;
        char fname[100];
        FILE *fp;

        clt_pp_init(&pp, &s->mlm);
        (void) snprintf(fname, 100, "%s/params", s->dir);
        fp = fopen(fname, "w+b");
        clt_pp_fsave(fp, &pp);
        fclose(fp);
        clt_pp_clear(&pp);
    }

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
    struct encode_elem_s *args;
    mmap_enc *out;
    fmpz_t *plaintext;
    const mmap_pp *pp = clt_vtable.sk->pp((mmap_sk *) &s->mlm);

    out = (mmap_enc *) malloc(sizeof(mmap_enc));
    plaintext = (fmpz_t *) calloc(2, sizeof(fmpz_t));
    clt_vtable.enc->init(out, pp);
    fmpz_init(plaintext[0]);
    fmpz_init(plaintext[1]);
    fmpz_set_mpz(plaintext[0], elem0);
    fmpz_set_mpz(plaintext[1], elem1);

    we_s = (struct write_element_s *) malloc(sizeof(write_element_s));
    we_s->dir = s->dir;
    we_s->elem = out;
    we_s->name = (char *) calloc(sizeof(int) + 10, sizeof(char));
    if (i >= 0) {
        (void) sprintf(we_s->name, str, i);
    } else {
        (void) strcpy(we_s->name, str);
    }
    (void) thpool_add_tag(s->thpool, we_s->name, 1, thpool_write_element, we_s);

    args = (struct encode_elem_s *) malloc(sizeof(struct encode_elem_s));
    args->vtable = &clt_vtable;
    args->sk = (mmap_sk *) &s->mlm;
    args->enc = out;
    args->n = 2;
    args->plaintext = plaintext;
    args->group = (int *) calloc(s->mlm.nzs, sizeof(int));
    va_start(vals, num);
    for (unsigned int i = 0; i < num; ++i) {
        int idx = va_arg(vals, int);
        args->group[idx] = va_arg(vals, int);
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
        mpz_urandomm_aes(alphas[i], s->rand, s->nchk);
    }
    betas = (mpz_t *) malloc(sizeof(mpz_t) * m);
    for (int i = 0; i < m; ++i) {
        mpz_init(betas[i]);
        mpz_urandomm_aes(betas[i], s->rand, s->nchk);
    }

    for (int i = 0; i < n; ++i) {
        mpz_t elems[2];
        int deg;

        mpz_inits(elems[0], elems[1], NULL);

        deg = PyLong_AsLong(PyList_GET_ITEM(py_xdegs, i));

        create_work(s, "x_\%d_0", i, zero, alphas[i], 1, 2 * i, 1);
        create_work(s, "x_\%d_1", i, one, alphas[i], 1, 2 * i + 1, 1);
        create_work(s, "u_\%d_0", i, one, one, 1, 2 * i, 1);
        create_work(s, "u_\%d_1", i, one, one, 1, 2 * i + 1, 1);

        mpz_urandomm_aes(elems[0], s->rand, s->nev);
        mpz_urandomm_aes(elems[1], s->rand, s->nchk);

        create_work(s, "z_\%d_0", i, elems[0], elems[1], 3, 2 * i + 1, deg,
                    2 * n + i, 1, 3 * n + i, 1);
        create_work(s, "w_\%d_0", i, zero, elems[1], 1, 3 * n + i, 1);

        mpz_urandomm_aes(elems[0], s->rand, s->nev);
        mpz_urandomm_aes(elems[1], s->rand, s->nchk);

        create_work(s, "z_\%d_1", i, elems[0], elems[1], 3, 2 * i, deg,
                    2 * n + i, 1, 3 * n + i, 1);
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
        (void) circ_evaluate(c, alphas, betas, c_star, s->mlm.x0);
        circ_cleanup(c);
    }

    // Last element to encode, so don't need to parallelize
    {
        mpz_t elems[2];
        int *pows;

        mpz_init_set_ui(elems[0], 0);
        mpz_init_set(elems[1], c_star);

        // The index set is laid out as follows:
        //   - The first 2 * n entries contain X_i,0 and X_i,1
        //   - The next n entries contain Z_i
        //   - The next n entries contain W_i
        //   - The final entry contains Y
        pows = (int *) calloc(s->mlm.nzs, sizeof(int));
        
        // The C* encoding contains everything but the W_i symbols.
        // Here we calculate the appropriate indices and powers.
        for (int i = 0; i < n; ++i) {
            int deg;

            deg = PyLong_AsLong(PyList_GET_ITEM(py_xdegs, i));
            // X_i,0^deg(x_i)
            pows[2 * i] = deg;
            // X_i,1^deg(x_i)
            pows[2 * i + 1] = deg;
            // Z_i
            pows[2 * n + i] = 1;
        }
        // Y^deg(y)
        pows[4 * n] = ydeg;
        // Encode against these indices/powers

        clt_encode(tmp, &s->mlm, 2, elems, pows, s->rand);

        mpz_clears(elems[0], elems[1], NULL);
        free(pows);
    }
    (void) write_element(s->dir, (mmap_enc *) &tmp, "c_star");

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
    clt_pp pp;
    mpz_t c_1, c_2, z, w;
    mpz_t *xs, *xones, *ys, *yones;

    if (!PyArg_ParseTuple(args, "ssslll", &dir, &circuit, &input, &n, &m,
                          &nthreads))
        return NULL;
    fnamelen = strlen(dir) + sizeof(int) + 5;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return NULL;

    mpz_inits(c_1, c_2, z, w, NULL);

    {
        char fname[100];
        FILE *fp;

        (void) snprintf(fname, 100, "%s/params", dir);
        fp = fopen(fname, "r+b");
        clt_pp_fread(fp, &pp);
        fclose(fp);

    }

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
        circ_evaluate_encoding(c, xs, xones, ys, yones, c_1, pp.x0);
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
        mpz_mod(c_1, c_1, pp.x0);
        (void) snprintf(fname, fnamelen, "%s/w_%d_%c", dir, i, input[i]);
        (void) load_mpz_scalar(fname, w);
        mpz_mul(c_2, c_2, w);
        mpz_mod(c_2, c_2, pp.x0);
    }

    // Compute c_1 - c_2 and zero test
    {
        mpz_t tmp;
        mpz_inits(tmp, NULL);
        mpz_sub(tmp, c_1, c_2);
        iszero = clt_is_zero(&pp, tmp);
        mpz_clears(tmp, NULL);
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
    mpz_clears(c_1, c_2, z, w, NULL);

    return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyObject *
obf_verbose(PyObject *self, PyObject *args)
{
    PyObject *py_verbose;

    if (!PyArg_ParseTuple(args, "O", &py_verbose))
        return NULL;

    g_verbose = PyObject_IsTrue(py_verbose);

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
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_zobfuscator(void)
{
    (void) Py_InitModule("_zobfuscator", ObfMethods);
}
