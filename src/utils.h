#ifndef __OBFUSCATION__UTILS_H__
#define __OBFUSCATION__UTILS_H__

#include <gmp.h>

double
current_time(void);

int
load_mpz_scalar(const char *fname, mpz_t x);

#endif
