// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_gcd.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_gcd.C,v 1.2 2007-01-11 18:42:51 jgdumas Exp $
// ==========================================================================
// Description:
//

#include "gmp++/gmp++.h"

// ==========================================================================
// Computes and returns the lcm of the two integers a and b.
Integer lcm(const Integer& a, const Integer& b) {
  Integer Res(Integer::one);
  mpz_lcm( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep) ) ;
  if (Res.priv_sign() <0) return -Res;
  else return Res ;
}

Integer& lcm(Integer& g, const Integer& a, const Integer& b) {
  mpz_lcm( (mpz_ptr)&(g.gmp_rep), (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep) ) ;
  if (g.priv_sign() <0) return Integer::negin(g);
  else return g ;
}


// ==========================================================================
// Computes and returns the gcd of the two integers a and b.
Integer gcd(const Integer& a, const Integer& b) {
  Integer Res(Integer::one);
  mpz_gcd( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep) ) ;
  if (Res.priv_sign() <0) return -Res;
  return Res ;
}

Integer& gcd(Integer& g, const Integer& a, const Integer& b) {
  mpz_gcd( (mpz_ptr)&(g.gmp_rep), (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep) ) ;
  if (g.priv_sign() <0) return Integer::negin(g);
  return g ;
}

Integer& inv(Integer& u, const Integer& a, const Integer& b) {
  mpz_invert( (mpz_ptr)&(u.gmp_rep), (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep) ) ;
  return u ;
}

// ==========================================================================
// Computes and returns the gcd g of the two integers a and b such that
// g = a*u + b*v .
// The algorithm used is this of Gmp.
Integer  gcd (const Integer& a, const Integer& b, Integer& u, Integer& v)
{
  v = 1; // v must not be 0 to be computed.
  Integer Res(Integer::one);
  mpz_gcdext( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&(u.gmp_rep), (mpz_ptr)&(v.gmp_rep),
              (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep) ) ;
  if (Res.priv_sign() < 0) { Integer::negin(u); Integer::negin(v); return Integer::negin(Res);}
//   { u = -u ; v = -v ; return -Res;}
  return Res;
}

Integer&  gcd (Integer& g, const Integer& a, const Integer& b, Integer& u, Integer& v)
{
  v = 1; // v must not be 0 to be computed.
  mpz_gcdext( (mpz_ptr)&(g.gmp_rep), (mpz_ptr)&(u.gmp_rep), (mpz_ptr)&(v.gmp_rep),
              (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep) ) ;
  if (g.priv_sign() < 0) { Integer::negin(u); Integer::negin(v); return Integer::negin(g);}
  return g;
}


Integer pp( const Integer& P, const Integer& Q )
{
  Integer U = P ;
  Integer V = gcd(P,Q) ;
  // -- computes the prime part U of g relatively to U
  while ( V != Integer::one )
  {
    U = U / V ;
    V = gcd( U,V) ;
  }
  return U ;
}

