#include "circuit.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>

enum GATE_TYPE {
    ADD,
    SUB,
    INPUT_X,
    INPUT_Y,
    MUL,
};

struct gate {
    GATE_TYPE type;
    mpz_t value;
    mpz_t one_value;
    struct gate *left;
    struct gate *right;
    struct gate *out;
};

struct circuit {
    int n_xins;
    int n_yins;
    int ngates;
    struct gate *gates;
};

static void
gate_init(struct gate *gate, GATE_TYPE type, struct gate *left,
          struct gate *right)
{
    gate->type = type;
    mpz_init(gate->value);
    mpz_init(gate->one_value);
    gate->left = left;
    gate->right = right;
    gate->out = NULL;
}

static void
gate_cleanup(struct gate *gate)
{
    mpz_clear(gate->value);
    mpz_clear(gate->one_value);
}

static void
circuit_init(struct circuit *circ, int ngates)
{
    circ->n_xins = 0;
    circ->n_yins = 0;
    circ->ngates = ngates;
    circ->gates = (struct gate *) malloc(sizeof(struct gate) * ngates);
}

static void
circuit_cleanup(struct circuit *circ)
{
    for (int i = 0; i < circ->ngates; ++i) {
        gate_cleanup(&circ->gates[i]);
    }
    free(circ->gates);
}

struct circuit *
circ_parse(const char *circname)
{
    FILE *f;
    char *line;
    size_t len;
    ssize_t read;
    struct circuit *circ;
    int ngates = 0;
    int err = 0;

    if ((f = fopen(circname, "r")) == NULL) {
        perror(circname);
        return NULL;
    }

    len = 128;                  // XXX: fixme
    line = (char *) malloc(sizeof(char) * len);

    // Count the number of gates
    while ((read = getline(&line, &len, f)) != -1) {
        char *type;

        if (line[0] == ':' || line[0] == '#' || line[0] == '\n')
            continue;
        if (strtok(line, " ") == NULL)
            continue;
        if ((type = strtok(NULL, " ")) == NULL)
            continue;
        if (strcmp(type, "input") == 0
            || strcmp(type, "gate") == 0
            || strcmp(type, "output") == 0) {
            ngates++;
        }
    }

    rewind(f);

    circ = (struct circuit *) malloc(sizeof(struct circuit));
    if (circ == NULL) {
        (void) fclose(f);
        return NULL;
    }
    circuit_init(circ, ngates);

    while ((read = getline(&line, &len, f)) != -1) {
        int num;
        char *type;

        if (line[0] == ':' || line[0] == '#' || line[0] == '\n')
            continue;

        num = atoi(strtok(line, " "));
        if ((type = strtok(NULL, " ")) == NULL) {
            (void) fprintf(stderr, "error: no type found!\n");
            err = 1;
            goto cleanup;
        }
        if (strcmp(type, "input") != 0
            && strcmp(type, "gate") != 0
            && strcmp(type, "output") != 0) {
            (void) fprintf(stderr, "warning: unknown type '%s'\n", type);
            continue;
        }

        if (strcmp(type, "input") == 0) {
            char *inp;
            if ((inp = strtok(NULL, " ")) == NULL) {
                (void) fprintf(stderr, "error: no input type found!\n");
                err = 1;
                goto cleanup;
            }
            if (inp[0] == 'x') {
                circ->n_xins++;
                gate_init(&circ->gates[num], INPUT_X, NULL, NULL);
            } else if (inp[0] == 'y') {
                circ->n_yins++;
                gate_init(&circ->gates[num], INPUT_Y, NULL, NULL);
            } else {
                (void) fprintf(stderr, "error: unknown input type '%s'\n", inp);
                err = 1;
                goto cleanup;
            }
        } else {
            char *gtype;
            int left, right;

            gtype = strtok(NULL, " ");
            left = atoi(strtok(NULL, " "));
            right = atoi(strtok(NULL, " "));
            if (strcmp(gtype, "ADD") == 0) {
                gate_init(&circ->gates[num], ADD, &circ->gates[left],
                          &circ->gates[right]);
            } else if (strcmp(gtype, "SUB") == 0) {
                gate_init(&circ->gates[num], SUB, &circ->gates[left],
                          &circ->gates[right]);
            } else if (strcmp(gtype, "MUL") == 0) {
                gate_init(&circ->gates[num], MUL, &circ->gates[left],
                          &circ->gates[right]);
            } else {
                (void) fprintf(stderr, "error: unknown gate type '%s'", gtype);
                err = 1;
                goto cleanup;
            }
        }
    }

cleanup:
    if (err) {
        circ_cleanup(circ);
        circ = NULL;
    }
    (void) fclose(f);

    return circ;
}

int
circ_copy_circuit(const char * const circname, const char * const newcirc)
{
    FILE *in, *out;
    char *line, *tmp;
    size_t len;
    ssize_t read;

    if ((in = fopen(circname, "r")) == NULL) {
        perror(circname);
        return 1;
    }

    if ((out = fopen(newcirc, "w")) == NULL) {
        perror(newcirc);
        fclose(in);
        return 1;
    }

    len = 128;                  // XXX: fixme
    line = (char *) malloc(sizeof(char) * len);
    tmp = (char *) malloc(sizeof(char) * len);

    while ((read = getline(&line, &len, in)) != -1) {
        memset(tmp, '\0', 128);
        (void) memcpy(tmp, line, 128);
        char *num = strtok(line, " ");
        if (num) {
            char *type = strtok(NULL, " ");
            if (type) {
                if (strcmp(type, "input") == 0) {
                    char *inp = strtok(NULL, " ");
                    if (inp && inp[0] == 'y') {
                        (void) snprintf(tmp, 128, "%s %s %s hidden\n", num, type, inp);
                        len = len > 128 ? 128 : len;
                        (void) fwrite(tmp, sizeof(char), strlen(tmp), out);
                        continue;
                    }
                }
            }
        }
        (void) fwrite(tmp, sizeof(char), strlen(tmp), out);
    }

    (void) fclose(in);
    (void) fclose(out);

    return 0;
}



int
circ_evaluate(const struct circuit *circ, const mpz_t *alphas,
              const mpz_t *betas, mpz_t out, const mpz_t q)
{
    int n_xins = 0;
    int n_yins = 0;

    for (int i = 0; i < circ->ngates; ++i) {
        struct gate *gate = &circ->gates[i];
        switch (gate->type) {
        case ADD:
            mpz_add(gate->value, gate->left->value, gate->right->value);
            mpz_mod(gate->value, gate->value, q);
            break;
        case SUB:
            mpz_sub(gate->value, gate->left->value, gate->right->value);
            mpz_mod(gate->value, gate->value, q);
            break;
        case INPUT_X:
            mpz_set(gate->value, alphas[n_xins++]);
            break;
        case INPUT_Y:
            mpz_set(gate->value, betas[n_yins++]);
            break;
        case MUL:
            mpz_mul(gate->value, gate->left->value, gate->right->value);
            mpz_mod(gate->value, gate->value, q);
            break;
        default:
            (void) fprintf(stderr, "error: unknown gate type!\n");
            return 1;
        }
    }
    mpz_set(out, circ->gates[circ->ngates - 1].value);
    return 0;
}

static void
add(mpz_t out, mpz_t out_one, const mpz_t x, const mpz_t x_one, const mpz_t y,
    const mpz_t y_one, const mpz_t q)
{
    mpz_t a, b;
    mpz_inits(a, b, NULL);

    mpz_mul(a, x, y_one);
    mpz_mul(b, x_one, y);
    mpz_add(out, a, b);
    mpz_mod(out, out, q);

    mpz_mul(out_one, x_one, y_one);
    mpz_mod(out_one, out_one, q);
    mpz_clears(a, b, NULL);
}

static void
sub(mpz_t out, mpz_t out_one, const mpz_t x, const mpz_t x_one, const mpz_t y,
    const mpz_t y_one, const mpz_t q)
{
    mpz_t a, b;
    mpz_inits(a, b, NULL);

    mpz_mul(a, x, y_one);
    mpz_mul(b, x_one, y);
    mpz_sub(out, a, b);
    mpz_mod(out, out, q);

    mpz_mul(out_one, x_one, y_one);
    mpz_mod(out_one, out_one, q);
    mpz_clears(a, b, NULL);
}

static void
multiply(mpz_t out, mpz_t out_one, const mpz_t x, const mpz_t x_one,
         const mpz_t y, const mpz_t y_one, const mpz_t q)
{
    mpz_mul(out, x, y);
    mpz_mod(out, out, q);
    mpz_mul(out_one, x_one, y_one);
    mpz_mod(out_one, out_one, q);
}

int
circ_evaluate_encoding(const struct circuit *circ, const mpz_t *xs,
                       const mpz_t *xones, const mpz_t *ys, const mpz_t *yones,
                       mpz_t out, const mpz_t q)
{
    int n_xins = 0;
    int n_yins = 0;

    for (int i = 0; i < circ->ngates; ++i) {
        struct gate *gate = &circ->gates[i];
        switch (gate->type) {
        case ADD:
            add(gate->value, gate->one_value, gate->left->value,
                gate->left->one_value, gate->right->value,
                gate->right->one_value, q);
            break;
        case SUB:
            sub(gate->value, gate->one_value, gate->left->value,
                gate->left->one_value, gate->right->value,
                gate->right->one_value, q);
            break;
        case INPUT_X:
            mpz_set(gate->value, xs[n_xins]);
            mpz_set(gate->one_value, xones[n_xins]);
            n_xins++;
            break;
        case INPUT_Y:
            mpz_set(gate->value, ys[n_yins]);
            mpz_set(gate->one_value, yones[n_yins]);
            n_yins++;
            break;
        case MUL:
            multiply(gate->value, gate->one_value, gate->left->value,
                     gate->left->one_value, gate->right->value,
                     gate->right->one_value, q);
            break;
        default:
            (void) fprintf(stderr, "error: unknown gate type!\n");
            return 1;
        }
    }
    mpz_set(out, circ->gates[circ->ngates - 1].value);
    return 0;
}

void
circ_cleanup(struct circuit *circ)
{
    circuit_cleanup(circ);
    free(circ);
}
