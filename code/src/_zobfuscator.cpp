#include "utils.h"
#include "pyutils.h"

static void
encode(struct state *s, mpz_t out, mpz_t in1, mpz_t in2,
	   const unsigned long *indices, unsigned long nindices)
{
	mpz_t r, tmp;

	mpz_inits(r, tmp, NULL);
	mpz_set_ui(out, 0);
	for (unsigned long i = 0; i < s->n; ++i) {
        mpz_genrandom(r, &s->rng, s->rho);
        mpz_mul(tmp, r, s->gs[i]);
        if (i < s->n / 2) {
            mpz_add(tmp, tmp, in1);
        } else {
			mpz_add(tmp, tmp, in2);
		}
        mpz_mul(tmp, tmp, s->crt_coeffs[i]);
        mpz_add(out, out, tmp);
	}
	mpz_mod(out, out, s->q);
	for (unsigned long i = 0; i < nindices; ++i) {
		mpz_mul(out, out, s->zinvs[indices[i]]);
        mpz_mod(out, out, s->q);
	}
	mpz_clears(r, tmp, NULL);
}

static PyObject *
obf_encode(PyObject *self, PyObject *args)
{
	PyObject *py_state, *py_in1, *py_in2, *py_indices;
	char *name;
	mpz_t val, in1, in2;
	double start, end;
	unsigned long *indices;
	unsigned long nindices;
	struct state *s;

	if (!PyArg_ParseTuple(args, "OOOOs", &py_state, &py_in1, &py_in2,
						  &py_indices, &name))
        return NULL;
    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

	nindices = PyList_GET_SIZE(py_indices);
	indices = (unsigned long *) malloc(sizeof(unsigned long) * nindices);
	if (indices == NULL)
		return NULL;
	for (unsigned long i = 0; i < nindices; ++i) {
		indices[i] = PyLong_AsUnsignedLong(PyList_GET_ITEM(py_indices, i));
	}

	mpz_inits(val, in1, in2, NULL);

	start = current_time();
	py_to_mpz(in1, py_in1);
	py_to_mpz(in2, py_in2);
	encode(s, val, in1, in2, indices, nindices);
	(void) write_scalar(s, val, name);
	end = current_time();
	if (g_verbose)
		(void) fprintf(stderr, "  Encoding element: %f\n", end - start);

	mpz_clears(val, in1, in2, NULL);

	free(indices);

	Py_RETURN_NONE;
}

static PyMethodDef
ObfMethods[] = {
	{"verbose", obf_verbose, METH_VARARGS,
     "Set verbosity."},
	{"setup", obf_setup, METH_VARARGS,
     "Set up obfuscator."},
    {"encode", obf_encode, METH_VARARGS,
     "Encode [a, b] under given index sets."},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_zobfuscator(void)
{
    (void) Py_InitModule("_zobfuscator", ObfMethods);
}

