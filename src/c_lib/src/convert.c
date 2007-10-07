/*

  convert.c
  2007 Aug 19
  Author: Craig Citro

  C code behind conversions to and from pari types.

  Currently implemented:
    Integer <--> t_INT

*/

#include "convert.h"

void t_INT_to_ZZ ( mpz_t value, GEN g )
{
  long limbs = 0;

  limbs = lgefint(g) - 2;

  mpz_realloc2( value, limbs );
  mpz_import( value, limbs, -1, 4, 0, 0, int_LSW(g) );

  if ( signe(g) == -1 )
    mpz_neg( value, value );

  return;
}

void ZZ_to_t_INT ( GEN g, mpz_t value )
{
  long limbs = 0;

  limbs = mpz_size( value );

  g = cgetg( limbs+2, t_INT );
  setlgefint( g, limbs+2 );
  setsigne( g, mpz_sgn(value) );

  mpz_export( int_LSW(g), NULL, -1, 4, 0, 0, value );

  return;
}

