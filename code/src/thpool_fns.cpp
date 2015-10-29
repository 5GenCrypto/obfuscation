#include "thpool_fns.h"

#include "pyutils.h"

void *thpool_encode_vector_elem(void *vargs)
{
    mpz_t *elems;
    struct mlm_encode_vector_elem_args *args =
        (struct mlm_encode_vector_elem_args *) vargs;

    elems = (mpz_t *) malloc(sizeof(mpz_t) * args->s->mlm->secparam);
    for (unsigned long j = 0; j < args->s->mlm->secparam; ++j) {
        mpz_init(elems[j]);
        py_to_mpz(elems[j],
                  PyList_GET_ITEM(PyList_GET_ITEM(args->s->py_vectors, j), args->i));
    }
    clt_mlm_encode(args->s->mlm, *args->elem,
                   args->s->mlm->secparam, elems, 2, args->s->indices,
                   args->s->pows);
    for (unsigned long j = 0; j < args->s->mlm->secparam; ++j) {
        mpz_clear(elems[j]);
    }
    free(elems);

    free(vargs);

    return NULL;
}
