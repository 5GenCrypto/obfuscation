#include "utils.h"

#include <fcntl.h>
#include <gmp.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <mmap/mmap_gghlite.h>

double
current_time(void)
{
    struct timeval t;
    (void) gettimeofday(&t, NULL);
    return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0));
}

FILE *
open_file(const char *dir, const char *file, const char *mode)
{
    FILE *fp;
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(file) + 2;
    fname = malloc(fnamelen);
    if (fname == NULL)
        return NULL;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, file);
    fp = fopen(fname, mode);
    if (fp == NULL) {
        fprintf(stderr, "unable to open '%s'\n", fname);
    }
    free(fname);
    return fp;
}

int
load_mpz_scalar(const char *fname, mpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    (void) mpz_inp_raw(x, f);
    (void) fclose(f);
    return 0;
}
