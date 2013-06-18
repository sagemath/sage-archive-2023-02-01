#include <gmp.h>

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif

// these vars are all used in rational reconstruction; they're cached so we don't
// have to recreate them with every call.
EXTERN mpz_t u, v, q, u0, u1, u2, v0, v1, v2, t0, t1, t2, x, y, ssqr, m2;
EXTERN mpq_t tmp;

EXTERN mpz_t a1, a2, mod1, sage_mod2, g, s, t, xx;

EXTERN mpz_t crtrr_a, crtrr_mod;

EXTERN mpz_t rand_val, rand_n, rand_n1;

EXTERN gmp_randstate_t rand_state;

EXTERN void init_mpz_globals();
EXTERN void clear_mpz_globals();

