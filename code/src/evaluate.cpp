#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <gmp.h>

#include "utils.h"

int
main(int argc, char *argv[])
{
    int fnamelen, i, layer, iszero;
    char *dir, *fname, *input;
    mpz_t *comp, *s, *t;
    mpz_t tmp, pzt, nu, q;
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

    mpz_inits(tmp, pzt, nu, q, NULL);

	(void) snprintf(fname, fnamelen, "%s/q", dir);
    (void) load_mpz_scalar(fname, q);
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

	(void) snprintf(fname, fnamelen, "%s/s_enc", dir);
    (void) load_mpz_vector(fname, s, width);

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
        if (input_idx >= strlen(input)) {
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
			(void) snprintf(fname, fnamelen, "%s/s_enc", dir);
			(void) load_mpz_vector(fname, s, width);
			mult_vect_by_mat(s, comp, q, width);
        } else {
			(void) load_mpz_vector(fname, comp, width * width);
			mult_vect_by_mat(s, comp, q, width);
        }
    }

    (void) snprintf(fname, fnamelen, "%s/t_enc", dir);
    (void) load_mpz_vector(fname, t, width);
    mult_vect_by_vect(tmp, s, t, q, width);
    (void) snprintf(fname, fnamelen, "%s/pzt", dir);
    (void) load_mpz_scalar(fname, pzt);
    (void) snprintf(fname, fnamelen, "%s/nu", dir);
    (void) load_mpz_scalar(fname, nu);

    iszero = is_zero(tmp, pzt, q, mpz_get_ui(nu));

    printf("Output: %d\n", iszero ? 0 : 1);

    return EXIT_SUCCESS;
}
