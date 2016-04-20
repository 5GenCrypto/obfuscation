#include "thpool_fns.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>

#include <mmap/mmap.h>
#include <mmap/mmap_clt.h>
#include <mmap/mmap_gghlite.h>

void *
thpool_encode_elem(void *vargs)
{
    struct encode_elem_s *args = (struct encode_elem_s *) vargs;
    aes_randstate_t rand;

    aes_randinit(rand);
    args->vtable->enc->encode(args->enc, args->sk, args->n, args->plaintext,
                              args->group, rand);
    aes_randclear(rand);

    for (int i = 0; i < args->n; ++i)
        fmpz_clear(args->plaintext[i]);
    free(args->plaintext);
    free(args->group);
    free(args);

    return NULL;
}

void *
thpool_write_layer(void *vargs)
{
    char fname[100];
    FILE *fp;
    double end;
    struct write_layer_s *args = (struct write_layer_s *) vargs;

    (void) snprintf(fname, 100, "%s/%ld.input", args->dir, args->idx);
    fp = fopen(fname, "w+b");
    fwrite(&args->inp, sizeof args->inp, 1, fp);
    fclose(fp);

    (void) snprintf(fname, 100, "%s/%ld.nrows", args->dir, args->idx);
    fp = fopen(fname, "w+b");
    fwrite(&args->nrows, sizeof args->nrows, 1, fp);
    fclose(fp);

    (void) snprintf(fname, 100, "%s/%ld.ncols", args->dir, args->idx);
    fp = fopen(fname, "w+b");
    fwrite(&args->ncols, sizeof args->ncols, 1, fp);
    fclose(fp);

    (void) snprintf(fname, 100, "%s/%ld.zero", args->dir, args->idx);
    fp = fopen(fname, "w+b");
    for (long i = 0; i < args->nrows * args->ncols; ++i) {
        args->vtable->enc->fwrite(&args->zero[i], fp);
        args->vtable->enc->clear(&args->zero[i]);
    }
    fclose(fp);

    (void) snprintf(fname, 100, "%s/%ld.one", args->dir, args->idx);
    fp = fopen(fname, "w+b");
    for (long i = 0; i < args->nrows * args->ncols; ++i) {
        args->vtable->enc->fwrite(&args->one[i], fp);
        args->vtable->enc->clear(&args->one[i]);
    }
    fclose(fp);

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       2 * args->nrows * args->ncols, end - args->start);

    free(args->zero);
    free(args->one);
    free(args);

    return NULL;
}

void *
thpool_write_element(void *vargs)
{
    FILE *fp;
    char fname[100];
    struct write_element_s *args = (struct write_element_s *) vargs;

    (void) snprintf(fname, 100, "%s/%s", args->dir, args->name);
    fp = fopen(fname, "w+b");
    clt13_vtable.enc->fwrite(args->elem, fp);
    fclose(fp);
    return NULL;
}
