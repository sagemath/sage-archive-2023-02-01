/*
 *  C helper functions for the computation
 *  of p-adic transcendantal functions
 *
 *********************************************/

#include <math.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

/* p-adic logarithm */
void padiclog(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo) {
    /*  Compute the p-adic logarithm of a,
        which is supposed to be a unit

        Algorithm: If a = 1 (mod p), write a as a product
            1/a = (1 - a_0*p) (1 - a_1*p^2) (1 - a_2*p^4) (1 - a_3*p^8) ...
        with 0 <= a_i < p^(2^i).
        Then compute each log(1 - a_i*p^(2^i)) using Taylor expansion
        and a binary spliting strategy.
        For general a, compute log(a^(p-1)) and then divide by p-1      */

    char congruent_to_one;
    unsigned long i, N, saveN, Np, tmp, trunc, step;
    mpz_t f, arg, trunc_mod, h, hpow, mpz_tmp, d, inv;
    mpz_t *num, *denom;

    mpz_init(mpz_tmp);
    mpz_init(arg);
    mpz_set_ui(ans, 0);

    mpz_fdiv_r_ui(mpz_tmp, a, p);
    congruent_to_one = (mpz_cmp_ui(mpz_tmp, 1) == 0);
    if (congruent_to_one) {
        mpz_set(arg, a);
    } else {
        mpz_powm_ui(arg, a, p-1, modulo);
    }

    /* Where do we need to truncate the Taylor expansion */
    N = prec;
    while(1) {
        tmp = prec + (unsigned long)(log(N)/log(p));
        if (tmp == N) break;
        N = tmp;
    }
    saveN = N;

    /* We allocate memory and initialize variables */
    mpz_init(f);
    mpz_init(h); mpz_init(hpow);
    mpz_init(d); mpz_init(inv);
    num = (mpz_t*)malloc(N*sizeof(mpz_t));
    denom = (mpz_t*)malloc(N*sizeof(mpz_t));
    for (i = 0; i < N; i++) {
        mpz_init(num[i]);
        mpz_init(denom[i]);
    }

    trunc = 2;
    mpz_init_set_ui(trunc_mod, p);
    mpz_mul_ui(trunc_mod, trunc_mod, p);
    while(1) {
        /* We compute f = 1 - a_i*p^(2^i)
           trunc_mod is p^(2^(i+1)) */
        mpz_fdiv_r(f, arg, trunc_mod);
        mpz_ui_sub(f, 2, f);
        mpz_mul(arg, arg, f);

        /* We compute the Taylor expansion of log(f)
           For now, computations are carried out over the rationals */
        for (i = 0; i < N; i++) {
            mpz_set_ui(num[i], 1);
            mpz_set_ui(denom[i], i+1);
        }
        step = 1;
        mpz_ui_sub(h, 1, f);   // we write f = 1 - h, i.e. h = a_i*p^(2^i)
        mpz_set(hpow, h);
        while(step < N) {
            for (i = 0; i < N - step; i += step << 1) {
                mpz_mul(mpz_tmp, hpow, num[i+step]);
                mpz_mul(mpz_tmp, mpz_tmp, denom[i]);
                mpz_mul(num[i], num[i], denom[i+step]);
                mpz_add(num[i], num[i], mpz_tmp);
                mpz_mul(denom[i], denom[i], denom[i+step]);
            }
            step <<= 1;
            mpz_mul(hpow, hpow, hpow);
        }

        /* We compute the p-adic valuation of the denominateur (which is N!) */
        Np = N; tmp = 0;
        while(Np > 0) { Np /= p; tmp += Np; }
        mpz_ui_pow_ui(d, p, tmp);
        mpz_divexact(mpz_tmp, num[0], d);
        mpz_mul(mpz_tmp, h, mpz_tmp);

        /* We coerce the result from Q to Zp */
        mpz_divexact(denom[0], denom[0], d);
        mpz_gcdext(d, inv, NULL, denom[0], modulo);
        mpz_mul(mpz_tmp, mpz_tmp, inv);

        /* We add this contribution to log(f) */
        mpz_add(ans, ans, mpz_tmp);

        if (trunc > prec) break;

        /* We update the variables for the next step */
        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        trunc <<= 1;
        N >>= 1;
    }

    if (! congruent_to_one) {
        mpz_set_ui(mpz_tmp, p-1);
        mpz_gcdext(d, inv, NULL, mpz_tmp, modulo);
        mpz_mul(ans, ans, inv);
    }

    mpz_fdiv_r(ans, ans, modulo);

    /* We clear memory */
    mpz_clear(f);
    mpz_clear(trunc_mod);
    mpz_clear(h);
    mpz_clear(hpow);
    mpz_clear(mpz_tmp);
    mpz_clear(d);
    mpz_clear(inv);
    for (i = 0; i < saveN; i++) {
        mpz_clear(num[i]);
        mpz_clear(denom[i]);
    }

}
