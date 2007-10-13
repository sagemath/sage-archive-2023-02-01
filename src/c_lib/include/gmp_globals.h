#include <gmp.h>

// these vars are all used in rational reconstruction; they're cached so we don't
// have to recreate them with every call.
extern mpz_t u, v, q, u0, u1, u2, v0, v1, v2, t0, t1, t2, x, y, sqr, m2;
extern mpq_t tmp;

extern mpz_t a1, a2, mod1, mod2, g, s, t, xx;

extern mpz_t crtrr_a, crtrr_mod;

extern mpz_t rand_val, rand_n, rand_n1;

extern gmp_randstate_t rand_state;

void init_mpz_globals();
void clear_mpz_globals();
