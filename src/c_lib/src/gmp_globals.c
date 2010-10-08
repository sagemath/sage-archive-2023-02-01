#include "gmp_globals.h"

mpz_t u, v, q, u0, u1, u2, v0, v1, v2, t0, t1, t2, x, y, ssqr, m2;
//changed sqr to ssqr due to a collision with ntl
mpq_t tmp;

mpz_t a1, a2, mod1, sage_mod2, g, s, t, xx;

mpz_t crtrr_a, crtrr_mod;

mpz_t rand_val, rand_n, rand_n1;

gmp_randstate_t rand_state;

void init_mpz_globals() {
  mpz_init(u);  mpz_init(v); mpz_init(q);
  mpz_init(u0); mpz_init(u1); mpz_init(u2);
  mpz_init(v0); mpz_init(v1); mpz_init(v2);
  mpz_init(t0); mpz_init(t1); mpz_init(t2);
  mpz_init(x);  mpz_init(y);
  mpz_init(ssqr);  mpz_init(m2);
  mpq_init(tmp);

  mpz_init(a1); mpz_init(a2); mpz_init(mod1); mpz_init(sage_mod2);
  mpz_init(g); mpz_init(s); mpz_init(t); mpz_init(xx);

  mpz_init(crtrr_a); mpz_init(crtrr_mod);

  mpz_init(rand_val); mpz_init(rand_n); mpz_init(rand_n1);

  gmp_randinit_default(rand_state);
}

void clear_mpz_globals() {
  mpz_clear(u);  mpz_clear(v); mpz_clear(q);
  mpz_clear(u0); mpz_clear(u1); mpz_clear(u2);
  mpz_clear(v0); mpz_clear(v1); mpz_clear(v2);
  mpz_clear(t0); mpz_clear(t1); mpz_clear(t2);
  mpz_clear(x);  mpz_clear(y);
  mpz_clear(ssqr);  mpz_clear(m2);
  mpq_clear(tmp);

  mpz_clear(a1); mpz_clear(a2); mpz_clear(mod1); mpz_clear(sage_mod2);
  mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(xx);

  mpz_clear(crtrr_a); mpz_clear(crtrr_mod);

  mpz_clear(rand_val); mpz_clear(rand_n); mpz_clear(rand_n1);
}

