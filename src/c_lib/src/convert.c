/*

  convert.c
  2007 Aug 19
  Author: Craig Citro

  C code behind conversions to and from pari types.

  Currently implemented:
    Integer <--> t_INT

*/

#include "convert.h"

/*

  Both t_INT_to_ZZ and ZZ_to_t_INT convert back and forth
  from mpz_t to PARI's t_INT GEN type. Nothing fancy happens
  here -- we simply use GMP's mpz_import or mpz_export to
  basically memcopy the limbs in the integer.

*/

void t_INT_to_ZZ ( mpz_t value, GEN g )
{
  long limbs = 0;

  limbs = lgefint(g) - 2;

  mpz_realloc2( value, limbs );
  mpz_import( value, limbs, -1, sizeof(long), 0, 0, int_LSW(g) );

  if ( signe(g) == -1 )
    mpz_neg( value, value );

  return;
}


/*

  Convert back and forth from mpq_t to PARI's t_FRAC GEN type. Nothing
  fancy happens here either.

  No type checking is done.

  AUTHOR: William Stein

*/

void t_FRAC_to_QQ ( mpq_t value, GEN g )
{
    t_INT_to_ZZ(mpq_numref(value), numer(g));
    t_INT_to_ZZ(mpq_denref(value), denom(g));
}
