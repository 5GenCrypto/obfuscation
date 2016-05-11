#ifndef OBFUSCATOR_H
#define OBFUSCATOR_H

#include <stdint.h>
#include <stdbool.h>
#include <flint/fmpz_mat.h>

#define OBFUSCATOR_OK 0
#define OBFUSCATOR_ERR (-1)

#define OBFUSCATOR_FLAG_NO_RANDOMIZATION 0x01
#define OBFUSCATOR_FLAG_DUAL_INPUT_BP 0x02
#define OBFUSCATOR_FLAG_VERBOSE 0x04

#ifdef __cplusplus
extern "C" {
#endif

typedef struct obf_state_s obf_state_t;

enum mmap_e { MMAP_CLT, MMAP_GGHLITE };

typedef enum {
    ENCODE_LAYER_RANDOMIZATION_TYPE_NONE = 0x00,
    ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST = 0x01,
    ENCODE_LAYER_RANDOMIZATION_TYPE_MIDDLE = 0x02,
    ENCODE_LAYER_RANDOMIZATION_TYPE_LAST = 0x04,
} encode_layer_randomization_flag_t;

obf_state_t *
obf_init(enum mmap_e type, const char *dir, uint64_t secparam, uint64_t kappa,
         uint64_t nzs, uint64_t nthreads, uint64_t ncores, uint64_t flags);

void
obf_clear(obf_state_t *s);

int
obf_encode_layer(obf_state_t *s, uint64_t n, int **pows, fmpz_mat_t *mats,
                 long idx, long inp, encode_layer_randomization_flag_t rflag);

int
obf_evaluate(enum mmap_e type, char *dir, uint64_t len, uint64_t *input,
             uint64_t bplen, uint64_t ncores, bool verbose);

void
obf_wait(obf_state_t *s);

#ifdef __cplusplus
}
#endif
    
#endif
