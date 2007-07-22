// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_pow.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: JG Dumas
// $Id: gmp++_int_pow.C,v 1.2 2007-01-11 18:42:51 jgdumas Exp $
// ==========================================================================
// Description:

#include "gmp++/gmp++.h"
int isperfectpower(const Integer& n) {
	return mpz_perfect_power_p((mpz_ptr)&(n.gmp_rep));
}

Integer& pow(Integer& Res, const unsigned long n, const unsigned long p)
{
  mpz_ui_pow_ui( (mpz_ptr)&(Res.gmp_rep), n, p);
  return Res;
}
Integer& pow(Integer& Res, const Integer& n, const unsigned long p)
{
__gmpz_pow_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p);
  return Res;
}

Integer pow(const Integer& n, const unsigned long p)
{
  if (p == 0) return Integer::one;

  Integer Res;
  return pow(Res,n,p);
}

Integer& pow(Integer& Res, const Integer& n, const long l) {
	return pow(Res, n, (unsigned long) GMP__ABS(l) );
}
Integer pow(const Integer& n, const long l) {
  if (l < 0)  return Integer::zero;
  return pow(n, (unsigned long) GMP__ABS(l) );
}

Integer& powmod(Integer& Res, const Integer& n, const unsigned long p, const Integer& m) {
  mpz_powm_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p, (mpz_ptr)&m.gmp_rep);
  return Res;
}

Integer powmod(const Integer& n, const unsigned long p, const Integer& m)
{
  if (p == 0) return Integer::one;
  Integer Res;
  return powmod(Res,n,p,m);
}

Integer& powmod(Integer& Res, const Integer& n, const long e, const Integer& m)
{
  return powmod (Res, n, (unsigned long)GMP__ABS(e), m);
}
Integer powmod(const Integer& n, const long e, const Integer& m)
{
  if (e < 0)  return Integer::zero;
  return powmod (n, (unsigned long)GMP__ABS(e), m);
}


Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m)
{
  mpz_powm( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, (mpz_ptr)&e.gmp_rep, (mpz_ptr)&m.gmp_rep);
  return Res;
}
Integer powmod(const Integer& n, const Integer& e, const Integer& m)
{
  if (e == 0) return Integer::one;
  if (e < 0)  return Integer::zero;
  Integer Res;
  return powmod(Res, n, e, m);
}
