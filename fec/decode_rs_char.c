/* General purpose Reed-Solomon decoder for 8-bit symbols or less
 * Copyright 2003 Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */

#ifdef DEBUG
#include <stdio.h>
#endif

#include <string.h>
#include <stdlib.h>

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#include "char.h"
#include "rs-common.h"
#undef MODNN
#undef MM
#undef NN
#undef ALPHA_TO
#undef INDEX_OF
#undef GENPOLY
#undef NROOTS
#undef FCR
#undef PRIM
#undef IPRIM
#undef PAD
#undef A0

int decode_rs_char(void *p, data_t *data, int *eras_pos, int no_eras)
{
    struct rs *rs = (struct rs *)p;

    int nroots = rs->nroots;
    int nn = rs->nn;
    int pad = rs->pad;
    int fcr = rs->fcr;
    int prim = rs->prim;
    int iprim = rs->iprim;
    const data_t *alpha_to = rs->alpha_to;
    const data_t *index_of = rs->index_of;
    const int A0 = nn;

    int deg_lambda, el, deg_omega;
    int i, j, r, k;
    data_t u, q, tmp, num1, num2, den, discr_r;
    int syn_error, count;
    int retval;

    data_t *lambda = (data_t *)calloc(nroots + 1, sizeof(data_t));
    data_t *s      = (data_t *)calloc(nroots, sizeof(data_t));
    data_t *b      = (data_t *)calloc(nroots + 1, sizeof(data_t));
    data_t *t      = (data_t *)calloc(nroots + 1, sizeof(data_t));
    data_t *omega  = (data_t *)calloc(nroots + 1, sizeof(data_t));
    data_t *root   = (data_t *)calloc(nroots, sizeof(data_t));
    data_t *reg    = (data_t *)calloc(nroots + 1, sizeof(data_t));
    data_t *loc    = (data_t *)calloc(nroots, sizeof(data_t));

    if (!lambda || !s || !b || !t || !omega || !root || !reg || !loc) {
        retval = -1;
        goto cleanup;
    }

    /* form the syndromes; i.e., evaluate data(x) at roots of g(x) */
    for (i = 0; i < nroots; i++)
        s[i] = data[0];

    for (j = 1; j < nn - pad; j++) {
        for (i = 0; i < nroots; i++) {
            if (s[i] == 0) {
                s[i] = data[j];
            } else {
                s[i] = data[j] ^ alpha_to[modnn(rs, index_of[s[i]] + (fcr + i) * prim)];
            }
        }
    }

    /* Convert syndromes to index form, checking for nonzero condition */
    syn_error = 0;
    for (i = 0; i < nroots; i++) {
        syn_error |= s[i];
        s[i] = index_of[s[i]];
    }

    if (!syn_error) {
        /* if syndrome is zero, data[] is a codeword and there are no
         * errors to correct. So return data[] unmodified
         */
        count = 0;
        goto finish;
    }

    memset(&lambda[1], 0, nroots * sizeof(lambda[0]));
    lambda[0] = 1;

    if (no_eras > 0) {
        /* Init lambda to be the erasure locator polynomial */
        lambda[1] = alpha_to[modnn(rs, prim * (nn - 1 - eras_pos[0]))];
        for (i = 1; i < no_eras; i++) {
            u = modnn(rs, prim * (nn - 1 - eras_pos[i]));
            for (j = i + 1; j > 0; j--) {
                tmp = index_of[lambda[j - 1]];
                if (tmp != A0)
                    lambda[j] ^= alpha_to[modnn(rs, u + tmp)];
            }
        }

#if DEBUG >= 1
        /* Test code that verifies the erasure locator polynomial */
        for (i = 1; i <= no_eras; i++)
            reg[i] = index_of[lambda[i]];

        count = 0;
        for (i = 1, k = iprim - 1; i <= nn; i++, k = modnn(rs, k + iprim)) {
            q = 1;
            for (j = 1; j <= no_eras; j++)
                if (reg[j] != A0) {
                    reg[j] = modnn(rs, reg[j] + j);
                    q ^= alpha_to[reg[j]];
                }
            if (q != 0)
                continue;
            root[count] = i;
            loc[count] = k;
            count++;
        }
        if (count != no_eras) {
            printf("count = %d no_eras = %d\n lambda(x) is WRONG\n", count, no_eras);
            count = -1;
            goto finish;
        }
#if DEBUG >= 2
        printf("\n Erasure positions as determined by roots of Eras Loc Poly:\n");
        for (i = 0; i < count; i++)
            printf("%d ", loc[i]);
        printf("\n");
#endif
#endif
    }
    for (i = 0; i < nroots + 1; i++)
        b[i] = index_of[lambda[i]];

    /*
     * Begin Berlekamp-Massey algorithm to determine error+erasure
     * locator polynomial
     */
    r = no_eras;
    el = no_eras;
    while (++r <= nroots) {
        /* Compute discrepancy at the r-th step in poly-form */
        discr_r = 0;
        for (i = 0; i < r; i++) {
            if ((lambda[i] != 0) && (s[r - i - 1] != A0)) {
                discr_r ^= alpha_to[modnn(rs, index_of[lambda[i]] + s[r - i - 1])];
            }
        }
        discr_r = index_of[discr_r]; /* Index form */
        if (discr_r == A0) {
            /* 2 lines below: B(x) <-- x*B(x) */
            memmove(&b[1], b, nroots * sizeof(b[0]));
            b[0] = A0;
        } else {
            /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
            t[0] = lambda[0];
            for (i = 0; i < nroots; i++) {
                if (b[i] != A0)
                    t[i + 1] = lambda[i + 1] ^ alpha_to[modnn(rs, discr_r + b[i])];
                else
                    t[i + 1] = lambda[i + 1];
            }
            if (2 * el <= r + no_eras - 1) {
                el = r + no_eras - el;
                /*
                 * 2 lines below: B(x) <-- inv(discr_r) * lambda(x)
                 */
                for (i = 0; i <= nroots; i++)
                    b[i] = (lambda[i] == 0) ? A0 : modnn(rs, index_of[lambda[i]] - discr_r + nn);
            } else {
                /* 2 lines below: B(x) <-- x*B(x) */
                memmove(&b[1], b, nroots * sizeof(b[0]));
                b[0] = A0;
            }
            memcpy(lambda, t, (nroots + 1) * sizeof(t[0]));
        }
    }

    /* Convert lambda to index form and compute deg(lambda(x)) */
    deg_lambda = 0;
    for (i = 0; i < nroots + 1; i++) {
        lambda[i] = index_of[lambda[i]];
        if (lambda[i] != A0)
            deg_lambda = i;
    }
    /* Find roots of the error+erasure locator polynomial by Chien search */
    memcpy(&reg[1], &lambda[1], nroots * sizeof(reg[0]));
    count = 0; /* Number of roots of lambda(x) */
    for (i = 1, k = iprim - 1; i <= nn; i++, k = modnn(rs, k + iprim)) {
        q = 1; /* lambda[0] is always 0 */
        for (j = deg_lambda; j > 0; j--) {
            if (reg[j] != A0) {
                reg[j] = modnn(rs, reg[j] + j);
                q ^= alpha_to[reg[j]];
            }
        }
        if (q != 0)
            continue; /* Not a root */
        /* store root (index-form) and error location number */
#if DEBUG >= 2
        printf("count %d root %d loc %d\n", count, i, k);
#endif
        root[count] = i;
        loc[count] = k;
        /* If we've already found max possible roots,
         * abort the search to save time
         */
        if (++count == deg_lambda)
            break;
    }
    if (deg_lambda != count) {
        /*
         * deg(lambda) unequal to number of roots => uncorrectable
         * error detected
         */
        count = -1;
        goto finish;
    }
    /*
     * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
     * x**nroots). in index form. Also find deg(omega).
     */
    deg_omega = deg_lambda - 1;
    for (i = 0; i <= deg_omega; i++) {
        tmp = 0;
        for (j = i; j >= 0; j--) {
            if ((s[i - j] != A0) && (lambda[j] != A0))
                tmp ^= alpha_to[modnn(rs, s[i - j] + lambda[j])];
        }
        omega[i] = index_of[tmp];
    }

    /*
     * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
     * inv(X(l))**(fcr-1) and den = lambda_pr(inv(X(l))) all in poly-form
     */
    for (j = count - 1; j >= 0; j--) {
        num1 = 0;
        for (i = deg_omega; i >= 0; i--) {
            if (omega[i] != A0)
                num1 ^= alpha_to[modnn(rs, omega[i] + i * root[j])];
        }
        num2 = alpha_to[modnn(rs, root[j] * (fcr - 1) + nn)];
        den = 0;

        /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
        for (i = (MIN(deg_lambda, nroots - 1) & ~1); i >= 0; i -= 2) {
            if (lambda[i + 1] != A0)
                den ^= alpha_to[modnn(rs, lambda[i + 1] + i * root[j])];
        }
#if DEBUG >= 1
        if (den == 0) {
            printf("\n ERROR: denominator = 0\n");
            count = -1;
            goto finish;
        }
#endif
        /* Apply error to data */
        if (num1 != 0 && loc[j] >= pad) {
            data[loc[j] - pad] ^=
                alpha_to[modnn(rs, index_of[num1] + index_of[num2] + nn - index_of[den])];
        }
    }
finish:
    if (eras_pos != NULL) {
        for (i = 0; i < count; i++)
            eras_pos[i] = loc[i];
    }
    retval = count;

cleanup:
    free(lambda);
    free(s);
    free(b);
    free(t);
    free(omega);
    free(root);
    free(reg);
    free(loc);
    return retval;
}
