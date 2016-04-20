#ifndef __OBFUSCATION__UTILS_H__
#define __OBFUSCATION__UTILS_H__

#include <gmp.h>
#include <clt13.h>
#include <gghlite.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>

extern int g_verbose;

#define OBFUSCATOR_OK 0
#define OBFUSCATOR_ERR (-1)

enum mmap_e { MMAP_CLT, MMAP_GGHLITE };

double
current_time(void);

int
load_mpz_scalar(const char *fname, mpz_t x);
/* int */
/* save_mpz_scalar(const char *fname, const mpz_t x); */

#endif
