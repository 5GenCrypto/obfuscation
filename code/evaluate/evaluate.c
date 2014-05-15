#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <gmp.h>

static void
mpz_mod_near(mpz_t out, const mpz_t a, const mpz_t b)
{
    mpz_t res, shift;

    mpz_inits(res, shift, NULL);
    mpz_mod(res, a, b);
    mpz_tdiv_q_2exp(shift, b, 1);
    if (mpz_cmp(res, shift) > 0) {
        mpz_sub(res, res, b);
    }
    mpz_set(out, res);
    mpz_clears(res, shift, NULL);
}

static void
load_mpz_scalar(const char *fname, mpz_t x)
{
    FILE *f;

    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
    } else {
        (void) mpz_inp_raw(x, f);
        (void) fclose(f);
    }
}

static void
load_mpz_vector(const char *fname, mpz_t *m, const int len)
{
    FILE *f;
    
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
    } else {
        int i;
        for (i = 0; i < len; ++i) {
            (void) mpz_inp_raw(m[i], f);
        }
        (void) fclose(f);
    }
}

static int
is_zero(const mpz_t c, const mpz_t pzt, const mpz_t x0, const long nu)
{
    mpz_t tmp;
    int ret;

    mpz_init(tmp);
    mpz_mul(tmp, c, pzt);
    mpz_mod_near(tmp, tmp, x0);
    ret = (mpz_sizeinbase(tmp, 2) < (mpz_sizeinbase(x0, 2) - nu)) ? 1 : 0;
    mpz_clear(tmp);

    return ret;
}

static void
mat_mult(mpz_t *a, mpz_t *b, const long width)
{
    int i, ctr;
    mpz_t *tmparray;

    tmparray = (mpz_t *) malloc(sizeof(mpz_t) * width * width);
    for (i = 0; i < width * width; ++i) {
        mpz_init(tmparray[i]);
    }
#pragma omp parallel for private(ctr, i)
    for (ctr = 0; ctr < width * width; ++ctr) {
        mpz_t tmp, sum;
        mpz_inits(tmp, sum, NULL);
        for (i = 0; i < width; ++i) {
            mpz_mul(tmp,
                    a[i * width + ctr % width],
                    b[i + width * (ctr / width)]);
            mpz_add(sum, sum, tmp);
        }
        mpz_set(tmparray[ctr], sum);
        mpz_clears(tmp, sum, NULL);
    }

    for (i = 0; i < width * width; ++i) {
        mpz_swap(a[i], tmparray[i]);
        mpz_clear(tmparray[i]);
    }
    free(tmparray);
}

static void
mat_mult_by_vects(mpz_t out, const mpz_t *s, const mpz_t *m, const mpz_t *t,
                  const int width)
{
    int col, row;
    mpz_set_ui(out, 0);

#pragma omp parallel for private(col)
    for (col = 0; col < width; ++col) {
        mpz_t tmp, sum;
        mpz_inits(tmp, sum, NULL);
        for (row = 0; row < width; ++row) {
            int elem = col * width + row;
            mpz_mul(tmp, s[row], m[elem]);
            mpz_add(sum, sum, tmp);
        }
        mpz_mul(tmp, sum, t[col]);
#pragma omp critical
        {
            mpz_add(out, out, tmp);
        }
        mpz_clears(tmp, sum, NULL);
    }
}

int
main(int argc, char *argv[])
{
    int fnamelen, i, layer, iszero;
    char *dir, *fname, *input;
    mpz_t *comp, *s, *t;
    mpz_t tmp, pzt, nu, x0;
    long width;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s directory input\n", argv[0]);
        return EXIT_FAILURE;
    }

    dir = argv[1];
    input = argv[2];
    printf("Input directory: %s\n", dir);
    printf("Input: %s\n", input);

    fnamelen = strlen(dir) + 20;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return EXIT_FAILURE;

    mpz_inits(tmp, pzt, nu, x0, NULL);

    (void) snprintf(fname, fnamelen, "%s/size", dir);
    (void) load_mpz_scalar(fname, tmp);
    width = mpz_get_ui(tmp);

    comp = (mpz_t *) malloc(sizeof(mpz_t) * width * width);
    if (comp == NULL)
        return EXIT_FAILURE;
    s = (mpz_t *) malloc(sizeof(mpz_t) * width);
    if (s == NULL)
        return EXIT_FAILURE;
    t = (mpz_t *) malloc(sizeof(mpz_t) * width);
    if (t == NULL)
        return EXIT_FAILURE;

    for (i = 0; i < width; ++i) {
        mpz_inits(s[i], t[i], NULL);
    }
    for (i = 0; i < width * width; ++i) {
        mpz_init(comp[i]);
    }

    for (layer = 0; /* empty */; ++layer) {
        unsigned int input_idx;
        struct stat st;

        // find out the input bit for the given layer
        (void) snprintf(fname, fnamelen, "%s/%d.input", dir, layer);
        // if file doesn't exist, break
        if (stat(fname, &st) == -1)
            break;
        (void) load_mpz_scalar(fname, tmp);
        input_idx = mpz_get_ui(tmp);
        if (input_idx < 0 || input_idx >= strlen(input)) {
            fprintf(stderr, "error: invalid input\n");
            return EXIT_FAILURE;
        }
        if (input[input_idx] != '0' && input[input_idx] != '1') {
            fprintf(stderr, "error: input must be 0 or 1\n");
            return EXIT_FAILURE;
        }
        // load in appropriate matrix for the given input value
        if (input[input_idx] == '0') {
            (void) snprintf(fname, fnamelen, "%s/%d.zero", dir, layer);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", dir, layer);
        }

        if (layer == 0) {
            (void) load_mpz_vector(fname, comp, width * width);
        } else {
            mpz_t *tmparray;

            tmparray = (mpz_t *) malloc(sizeof(mpz_t) * width * width);
            if (tmparray == NULL)
                return EXIT_FAILURE;
            for (i = 0; i < width * width; ++i) {
                mpz_init(tmparray[i]);
            }
            (void) load_mpz_vector(fname, tmparray, width * width);
            mat_mult(comp, tmparray, width);
            for (i = 0; i < width * width; ++i) {
                mpz_clear(tmparray[i]);
            }
            free(tmparray);
        }
    }

    (void) snprintf(fname, fnamelen, "%s/s_enc", dir);
    (void) load_mpz_vector(fname, s, width);
    (void) snprintf(fname, fnamelen, "%s/t_enc", dir);
    (void) load_mpz_vector(fname, t, width);
    mat_mult_by_vects(tmp, s, comp, t, width);
    (void) snprintf(fname, fnamelen, "%s/pzt", dir);
    (void) load_mpz_scalar(fname, pzt);
    (void) snprintf(fname, fnamelen, "%s/x0", dir);
    (void) load_mpz_scalar(fname, x0);
    (void) snprintf(fname, fnamelen, "%s/nu", dir);
    (void) load_mpz_scalar(fname, nu);

    iszero = is_zero(tmp, pzt, x0, mpz_get_ui(nu));

    printf("Output: %d\n", iszero ? 0 : 1);

    return EXIT_SUCCESS;
}
