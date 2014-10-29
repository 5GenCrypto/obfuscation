#include "utils.h"
#include "pyutils.h"

#include <omp.h>
#ifdef ATTACK
#include <fplll.h>
#endif

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
        return FAILURE;
    }
    return SUCCESS;
}

static int
write_vector(struct state *s, mpz_t *vector, long size, char *name)
{
    char *fname;
    int fnamelen;
    double start, end;

    start = current_time();
    fnamelen = strlen(s->dir) + strlen(name) + 2;
    fname = (char *) pymalloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return FAILURE;
    (void) snprintf(fname, fnamelen, "%s/%s", s->dir, name);
    (void) save_mpz_vector(fname, vector, size);
    free(fname);

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Saving to file: %f\n", end - start);
    return SUCCESS;
}

static int
write_layer(struct state *s, int inp, long idx, mpz_t *zero, mpz_t *one,
			long nrows, long ncols)
{
    mpz_t tmp;
    char *fname;
    int fnamelen;
    double start, end;

    start = current_time();
    fnamelen = strlen(s->dir) + 10;
    fname = (char *) pymalloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return FAILURE;
    mpz_init_set_ui(tmp, inp);
    (void) snprintf(fname, fnamelen, "%s/%ld.input", s->dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    (void) snprintf(fname, fnamelen, "%s/%ld.zero", s->dir, idx);
    (void) save_mpz_vector(fname, zero, nrows * ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.one", s->dir, idx);
    (void) save_mpz_vector(fname, one, nrows * ncols);
	mpz_set_ui(tmp, nrows);
    (void) snprintf(fname, fnamelen, "%s/%ld.nrows", s->dir, idx);
    (void) save_mpz_scalar(fname, tmp);
	mpz_set_ui(tmp, ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.ncols", s->dir, idx);
    (void) save_mpz_scalar(fname, tmp);

    free(fname);
    mpz_clear(tmp);
    end = current_time();

    if (g_verbose)
        (void) fprintf(stderr, "  Saving to file: %f\n", end - start);
    return SUCCESS;
}

static void
encode(struct state *st, mpz_t out, const PyObject *in, const long item,
       const long idx1, const long idx2)
{
    mpz_t r, tmp;

    mpz_inits(r, tmp, NULL);

    mpz_set_ui(out, 0);

    for (unsigned long i = 0; i < st->n; ++i) {
        mpz_genrandom(r, &st->rng, st->rho);
        mpz_mul(tmp, r, st->gs[i]);

        if (i < st->secparam) {
            py_to_mpz(r, PyList_GET_ITEM(PyList_GET_ITEM(in, i), item));
            mpz_add(tmp, tmp, r);
        }
        mpz_mul(tmp, tmp, st->crt_coeffs[i]);
        mpz_add(out, out, tmp);
    }
    mpz_mod(out, out, st->q);
    if (idx1 >= 0) {
        mpz_mul(out, out, st->zinvs[idx1]);
        mpz_mod(out, out, st->q);
    }
    if (idx2 >= 0) {
        mpz_mul(out, out, st->zinvs[idx2]);
        mpz_mod(out, out, st->q);
    }

    mpz_clears(r, tmp, NULL);
}

//
//
// Python functions
//
//

//
// Encode N vectors across all slots of the MLM
//
static PyObject *
obf_encode_vectors(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_vectors, *py_list;
    char *name;
    int idx1, idx2;
    mpz_t *vector;
    Py_ssize_t length;
    double start, end;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OOOs", &py_state, &py_vectors, &py_list, &name))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;
    
    (void) extract_indices(py_list, &idx1, &idx2);

    // We assume that all vectors have the same length, and thus just grab the
    // length of the first vector
    length = PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0));
    vector = (mpz_t *) pymalloc(sizeof(mpz_t) * length);
    if (vector == NULL)
        return NULL;

    start = current_time();
#pragma omp parallel for
    for (Py_ssize_t i = 0; i < length; ++i) {
        mpz_init(vector[i]);
        encode(s, vector[i], py_vectors, i, idx1, idx2);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       length, end - start);

    (void) write_vector(s, vector, length, name);

    for (Py_ssize_t i = 0; i < length; ++i) {
        mpz_clear(vector[i]);
    }
    free(vector);

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
    int zeroidx1, zeroidx2, oneidx1, oneidx2;
    int err = 0;
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
    
    (void) extract_indices(py_zero_set, &zeroidx1, &zeroidx2);
    (void) extract_indices(py_one_set, &oneidx1, &oneidx2);

    if (zeroidx1 < 0 && zeroidx2 < 0)
        return NULL;
    if (oneidx1 < 0 && oneidx2 < 0)
        return NULL;

    // size = PyList_GET_SIZE(PyList_GET_ITEM(py_zero_ms, 0));
    zero = (mpz_t *) pymalloc(sizeof(mpz_t) * nrows * ncols);
    one = (mpz_t *) pymalloc(sizeof(mpz_t) * nrows * ncols);
    if (!zero || !one)
        return NULL;

    start = current_time();
#pragma omp parallel for
    for (Py_ssize_t ctr = 0; ctr < 2 * nrows * ncols; ++ctr) {
        PyObject *py_array;
        int idx1, idx2;
        mpz_t *val;
        size_t i;

        if (ctr < nrows * ncols) {
            i = ctr;
            val = &zero[i];
            py_array = py_zero_ms;
            idx1 = zeroidx1;
            idx2 = zeroidx2;
        } else {
            i = ctr - nrows * ncols;
            val = &one[i];
            py_array = py_one_ms;
            idx1 = oneidx1;
            idx2 = oneidx2;
        }

        mpz_init(*val);
        encode(s, *val, py_array, i, idx1, idx2);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       2 * nrows * ncols, end - start);

    (void) write_layer(s, inp, idx, zero, one, nrows, ncols);

    for (int i = 0; i < nrows * ncols; ++i) {
        mpz_clears(zero[i], one[i], NULL);
    }
    free(zero);
    free(one);

    if (err)
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
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
	mpz_t *result;
    long bplen, nrows, ncols, nrows_prev;
    int err = 0;
	double start, end;

    if (!PyArg_ParseTuple(args, "ssl", &dir, &input, &bplen))
        return NULL;

    fnamelen = strlen(dir) + 20; // XXX: should include bplen somewhere

    fname = (char *) pymalloc(sizeof(char) * fnamelen);
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
		iszero = is_zero(result[1], pzt, q, mpz_get_ui(nu));
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
	fnamelen = strlen(dir) + 20; // XXX: should include bplen somewhere
	fname = (char *) pymalloc(sizeof(char) * fnamelen);
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
	comp = (mpz_t *) pymalloc(sizeof(mpz_t) * size * size);
	s = (mpz_t *) pymalloc(sizeof(mpz_t) * size);
	t = (mpz_t *) pymalloc(sizeof(mpz_t) * size);
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
		if (layer == 0) {
			(void) load_mpz_vector(fname, comp, size * size);
			(void) snprintf(fname, fnamelen, "%s/s_enc", dir);
			(void) load_mpz_vector(fname, s, size);
			mult_vect_by_mat(s, comp, q, size);
		} else {
			(void) load_mpz_vector(fname, comp, size * size);
			mult_vect_by_mat(s, comp, q, size);
		}
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
		{
			mpz_t pzt, nu;
			mpz_inits(pzt, nu, NULL);
			(void) snprintf(fname, fnamelen, "%s/pzt", dir);
			(void) load_mpz_scalar(fname, pzt);
			(void) snprintf(fname, fnamelen, "%s/nu", dir);
			(void) load_mpz_scalar(fname, nu);
			iszero = is_zero(tmp, pzt, q, mpz_get_ui(nu));
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

    gmp_randclear(s->rng);
    mpz_clears(s->q, s->pzt, NULL);
    for (unsigned long i = 0; i < s->n; ++i) {
        mpz_clears(s->gs[i], s->crt_coeffs[i], NULL);
    }
    free(s->gs);
    free(s->crt_coeffs);
    for (unsigned long i = 0; i < s->nzs; ++i) {
        mpz_clear(s->zinvs[i]);
    }
    free(s->zinvs);

    Py_RETURN_NONE;
}

#ifdef ATTACK
static PyObject *
obf_attack(PyObject *self, PyObject *args)
{
    mpz_t b, out, omega, Omega, q;
    double start, end;
    struct state *st;
    long kappa, bplen;
    int nslots;
    PyObject *py_out;

    st = (struct state *) pymalloc(sizeof(struct state));

    if (!PyArg_ParseTuple(args, "sllli", &st->dir, &bplen,
                          &st->secparam, &kappa, &nslots))
        return NULL;

    mpz_inits(b, out, omega, Omega, q, NULL);

    /* Compute omega = pzt * u for some encoded element u */

    if (g_verbose)
        (void) fprintf(stderr, "Computing omega...\n");
    start = current_time();
    {
        char *fname;
        int fnamelen;
        mpz_t *comp;
        mpz_t tmp, pzt, nu;

        fnamelen = strlen(st->dir) + 20; // XXX: should include bplen somewhere
        fname = (char *) pymalloc(sizeof(char) * fnamelen);

        mpz_inits(tmp, pzt, nu, NULL);

        // Get the size of the matrices
        (void) snprintf(fname, fnamelen, "%s/size", st->dir);
        (void) load_mpz_scalar(fname, tmp);
        st->size = mpz_get_ui(tmp);

        comp = (mpz_t *) pymalloc(sizeof(mpz_t) * st->size * st->size);
        for (int i = 0; i < st->size * st->size; ++i) {
            mpz_init(comp[i]);
        }

        for (int layer = 0; layer < bplen; ++layer) {
            (void) snprintf(fname, fnamelen, "%s/%d.zero", st->dir, layer);
            (void) load_mpz_vector(fname, comp, st->size * st->size);
            if (layer == 0) {
                mpz_set(tmp, comp[0]);
            } else {
                mpz_mul(tmp, tmp, comp[0]);
            }
        }

        (void) snprintf(fname, fnamelen, "%s/s_enc", st->dir);
        (void) load_mpz_vector(fname, comp, st->size);
        mpz_mul(tmp, tmp, comp[0]);
        (void) snprintf(fname, fnamelen, "%s/t_enc", st->dir);
        (void) load_mpz_vector(fname, comp, st->size);
        mpz_mul(tmp, tmp, comp[st->size - 1]);

        (void) snprintf(fname, fnamelen, "%s/pzt", st->dir);
        (void) load_mpz_scalar(fname, pzt);
        (void) snprintf(fname, fnamelen, "%s/q", st->dir);
        (void) load_mpz_scalar(fname, q);
        // (void) snprintf(fname, fnamelen, "%s/nu", st->dir);
        // (void) load_mpz_scalar(fname, nu);

        // {
        //     int iszero = is_zero(tmp, pzt, q, mpz_get_ui(nu));
        //     fprintf(stderr, "iszero = %d\n", iszero);
        // }

        // Compute omega
        mpz_mul(tmp, tmp, pzt);
        mpz_mod_near(omega, tmp, q);

        mpz_clears(tmp, pzt, nu, NULL);

        for (int i = 0; i < st->size * st->size; ++i) {
            mpz_clear(comp[i]);
        }
        if (comp)
            free(comp);
        if (fname)
            free(fname);
    }
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Took: %f\n", end - start);

    /* We now have omega and q.  We now need to sample everything fresh to
       compute Omega. */

    {
        long alpha, beta, eta, nu, rho_f;
        mpz_t *hs, *ms, *ps, *rs;

        /* Calculate CLT parameters */
        alpha = st->secparam;
        beta = st->secparam;
        st->rho = st->secparam;
        rho_f = kappa * (st->rho + alpha + 2);
        eta = rho_f + alpha + 2 * beta + st->secparam + 8;
        nu = eta - beta - rho_f - st->secparam + 3;
        st->n = (int) (eta * log2((float) st->secparam));
        if (g_verbose) {
            (void) fprintf(stderr, "Parameters:\n");
            (void) fprintf(stderr, "  Security Parameter: %ld\n", st->secparam);
            (void) fprintf(stderr, "  Kappa: %ld\n", kappa);
            (void) fprintf(stderr, "  Alpha: %ld\n", alpha);
            (void) fprintf(stderr, "  Beta: %ld\n", beta);
            (void) fprintf(stderr, "  Eta: %ld\n", eta);
            (void) fprintf(stderr, "  Nu: %ld\n", nu);
            (void) fprintf(stderr, "  Rho: %ld\n", st->rho);
            (void) fprintf(stderr, "  Rho_f: %ld\n", rho_f);
            (void) fprintf(stderr, "  N: %ld\n", st->n);
            (void) fprintf(stderr, "  Size: %ld\n", st->size);
        }

        hs = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);
        ms = (mpz_t *) pymalloc(sizeof(mpz_t) * nslots);
        ps = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);
        rs = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);
        st->gs = (mpz_t *) pymalloc(sizeof(mpz_t) * st->n);

        seed_rng(&st->rng);

        /* initialize gmp variables */
        mpz_init_set_ui(st->pzt, 0);
        for (int i = 0; i < nslots; ++i) {
            mpz_inits(ms[i], NULL);
        }
        for (unsigned int i = 0; i < st->n; ++i) {
            mpz_inits(hs[i], ps[i], rs[i], st->gs[i], NULL);
        }

        /* Generate p_i's and g_i's */
        start = current_time();
#pragma omp parallel for
        for (unsigned int i = 0; i < st->n; ++i) {
            mpz_t p_unif;
            mpz_init(p_unif);
            mpz_urandomb(p_unif, st->rng, eta);
            mpz_nextprime(ps[i], p_unif);
            mpz_urandomb(p_unif, st->rng, alpha);
            mpz_nextprime(st->gs[i], p_unif);
            mpz_clear(p_unif);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating p_i's and g_i's: %f\n",
                           end - start);

        /* Compute h_i's */
        start = current_time();
#pragma omp parallel for
        for (unsigned int i = 0; i < st->n; ++i) {
            mpz_genrandom(hs[i], &st->rng, beta);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating h_i's: %f\n", end - start);

        /* Compute r_i's */
        start = current_time();
#pragma omp parallel for
        for (unsigned int i = 0; i < st->n; ++i) {
            mpz_genrandom(rs[i], &st->rng, rho_f);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating r_i's: %f\n", end - start);

        /* Sample m_i's */
        start = current_time();
#pragma omp parallel for
        for (int i = 0; i < nslots; ++i) {
            mpz_urandomb(ms[i], st->rng, alpha);
            mpz_mod(ms[i], ms[i], st->gs[i]);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "Generating m_i's: %f\n", end - start);

        /* Compute b = \Omega \cdot g_1 */

        start = current_time();
        mpz_set_ui(b, 0L);
#pragma omp parallel for
        for (int i = 0; i < nslots; ++i) {
            mpz_t tmp, qpi;
            mpz_inits(tmp, qpi, NULL);
            mpz_mul(tmp, hs[i], ms[i]);
            mpz_div(qpi, q, ps[i]);
            mpz_mul(tmp, tmp, qpi);
            for (int j = 0; j < nslots; ++j) {
                if (j != i)
                    mpz_mul(tmp, tmp, st->gs[j]);
            }
#pragma omp critical
            {
                mpz_add(b, b, tmp);
            }
            mpz_clears(tmp, qpi, NULL);
        }
#pragma omp parallel for
        for (unsigned int i = 0; i < st->n; ++i) {
            mpz_t tmp, qpi;
            mpz_inits(tmp, qpi, NULL);
            mpz_mul(tmp, hs[i], rs[i]);
            mpz_div(qpi, q, ps[i]);
            mpz_mul(tmp, tmp, qpi);
            for (int j = 0; j < nslots; ++j) {
                mpz_mul(tmp, tmp, st->gs[j]);
            }
#pragma omp critical
            {
                mpz_add(b, b, tmp);
            }
            mpz_clears(tmp, qpi, NULL);
        }
        mpz_mod(b, b, q);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Computing b: %f\n", end - start);

        for (int i = 0; i < nslots; ++i) {
            mpz_clears(ms[i], NULL);
        }
        for (unsigned int i = 0; i < st->n; ++i) {
            mpz_clears(hs[i], ps[i], rs[i], NULL);
        }
        free(hs);
        free(ms);
        free(ps);
        free(rs);
    }

    /* Compute Omega = b / g_1 */
    mpz_set(Omega, b);
    for (int i = 0; i < nslots; ++i) {
        mpz_div(Omega, Omega, st->gs[i]);
    }

    {
        ZZ_mat<mpz_t> M(2, 2);

        /* Apply LLL to lattice basis {(Omega, omega), (0, q)} */
        if (g_verbose)
            (void) fprintf(stderr, "Applying LLL...\n");
        start = current_time();
        M(0,0).set(Omega);
        M(0,1).set(omega);
        M(1,0) = 0L;
        M(1,1).set(q);
        fplll::lllReduction(M, 0.99, 0.51, LM_WRAPPER);
        /* Divide out Omega from output to get g_1 */
        mpz_div(out, M(0,0).getData(), Omega);
        mpz_abs(out, out);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Took: %f\n", end - start);
    }

    py_out = mpz_to_py(out);

    mpz_clears(b, out, omega, Omega, NULL);

    if (st) {
        if (st->gs)
            free(st->gs);
    }

    return py_out;
}
#endif


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
#ifdef ATTACK
    {"attack", obf_attack, METH_VARARGS,
     "implementation of attack from paper (Section 3.1.2)."},
#endif
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
