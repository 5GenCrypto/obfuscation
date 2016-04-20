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

int g_verbose;

double
current_time(void)
{
    struct timeval t;
    (void) gettimeofday(&t, NULL);
    return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0));
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

// int
// save_mpz_scalar(const char *fname, const mpz_t x)
// {
//     FILE *f;
//     if ((f = fopen(fname, "w")) == NULL) {
//         perror(fname);
//         return 1;
//     }
//     if (mpz_out_raw(f, x) == 0) {
//         (void) fclose(f);
//         return 1;
//     }
//     (void) fclose(f);
//     return 0;
// }
