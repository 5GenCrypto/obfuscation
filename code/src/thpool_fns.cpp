#include "thpool_fns.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>

void *
thpool_encode_vector_elem(void *vargs)
{
    struct mlm_encode_vector_elem_s *args =
        (struct mlm_encode_vector_elem_s *) vargs;

    clt_mlm_encode(args->s->mlm, *args->out,
                   args->s->mlm->secparam, args->elems, 2, args->s->indices,
                   args->s->pows);
    for (unsigned long j = 0; j < args->s->mlm->secparam; ++j) {
        mpz_clear(args->elems[j]);
    }
    free(args->elems);
    free(args);

    return NULL;
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

void *
thpool_write_vector(void *vargs)
{
    double end;
    struct write_vector_s *args = (struct write_vector_s *) vargs;

    // fprintf(stderr, "DONE WITH %s\n", args->name);

    (void) write_vector(args->dir, args->vector, args->length, args->name);
    for (unsigned long i = 0; i < args->length; ++i) {
        mpz_clear(args->vector[i]);
    }

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       args->length, end - args->start);

    free(args->dir);
    free(args->name);
    free(args->vector);
    // free(args->state->indices);
    // free(args->state->pows);
    // free(args->state);
    free(args);

    return NULL;
}
