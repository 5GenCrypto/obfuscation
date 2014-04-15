
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
load_mpz_vector(const char *fname, mpz_t *m, const int len)
{
    int ret = SUCCESS;
    FILE *f;
    
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return FAILURE;
    }

    for (int i = 0; i < len; ++i) {
        if (mpz_inp_raw(m[i], f) == 0) {
            ret = FAILURE;
            goto cleanup;
        }
    }
 cleanup:
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
        if (mpz_out_raw(f, m[i]) == 0) {
            ret = FAILURE;
            goto cleanup;
        }
    }
 cleanup:
    (void) fclose(f);
    return ret;
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

    if (mpz_inp_raw(x, f) == 0) {
        ret = FAILURE;
        goto cleanup;
    }

 cleanup:
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
    if (mpz_out_raw(f, x) == 0) {
        ret = FAILURE;
        goto cleanup;
    }

 cleanup:
    (void) fclose(f);
    return ret;
}

inline static int
check_pylist(PyObject *list)
{
    PyObject *item;
    Py_ssize_t i, len;

    if (!PyList_Check(list)) {
        return FAILURE;
    }

    len = PyList_GET_SIZE(list);
    for (i = 0; i < len; ++i) {
        item = PyList_GET_ITEM(list, i);
        if (!PyLong_Check(item)) {
            return FAILURE;
        }
    }
    return SUCCESS;
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

static void
mpz_genrandom(mpz_t rnd, const long nbits)
{
    mpz_t one;
    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, g_rng, (mp_bitcnt_t) nbits);
    mpz_clear(one);
}

static void
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

static int
encode(mpz_t out, const mpz_t in, const long idx1, const long idx2,
       const long slot)
{
    mpz_t res, r, tmp;
    long i;

    if (idx1 >= g_nzs || idx2 >= g_nzs)
        return FAILURE;
    if (idx1 < 0 && idx2 < 0)
        return FAILURE;
    if (idx1 == idx2)
        return FAILURE;

    mpz_init(res);
    mpz_init(r);
    mpz_init(tmp);

    mpz_set_ui(res, 0);

    for (i = 0; i < g_n; ++i) {
        mpz_genrandom(r, g_rho);
        mpz_mul(tmp, r, g_gs[i]);
        if (i == slot) {
            mpz_add(tmp, tmp, in);
        }
        mpz_mul(tmp, tmp, g_crt_coeffs[i]);
        mpz_add(res, res, tmp);
    }
    mpz_mod(res, res, g_x0);
    if (idx1 >= 0) {
        mpz_mul(res, res, g_zinvs[idx1]);
        mpz_mod(res, res, g_x0);
    }
    if (idx2 >= 0) {
        mpz_mul(res, res, g_zinvs[idx2]);
        mpz_mod(res, res, g_x0);
    }

    mpz_set(out, res);

    mpz_clear(res);
    mpz_clear(r);
    mpz_clear(tmp);

    return SUCCESS;
}

static void
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

static void
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
