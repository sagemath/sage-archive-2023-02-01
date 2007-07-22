// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_mod.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_mod.C,v 1.6 2007-01-24 15:40:38 jgdumas Exp $
// ==========================================================================

#include "gmp++/gmp++.h"


//-------------------------------------------------- operator /
Integer& Integer::modin(Integer& res, const Integer& n)
{
//  if (isZero(n)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (isZero(res)) return res;
  mpz_tdiv_r( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
  return res;
}
Integer& Integer::modin(Integer& res, const long n)
{
//  if (n ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (isZero(res)) return res;
  int sgn = GMP__SGN(n);
  mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, GMP__ABS(n));
  if (sgn <0) return res = -res;
  return res;
}
Integer& Integer::modin(Integer& res, const unsigned long n)
{
//  if (n ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (isZero(res)) return res;
  mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
  return res;
}

Integer& Integer::mod(Integer& res, const Integer& n1, const Integer& n2)
{
  if (isZero(n1)) return res = Integer::zero;
//  if (isZero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  mpz_tdiv_r( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
  return res;
}
Integer& Integer::mod(Integer& res, const Integer& n1, const long n2)
{
  if (isZero(n1)) return res = Integer::zero;
//  if (isZero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  int sgn = GMP__SGN(n2);
  mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, GMP__ABS(n2));
  if (sgn <0) return res = - res;
  return res;
}
Integer& Integer::mod(Integer& res, const Integer& n1, const unsigned long n2)
{
  if (isZero(n1)) return res = Integer::zero;
//  if (isZero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
  return res;
}


Integer& Integer::operator %= (const Integer& n)
{
//  if (isZero(n)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (isZero(*this)) return *this;
  Integer res;
  mpz_tdiv_r( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  mpz_set( (mpz_ptr)&gmp_rep , (mpz_ptr)&(res.gmp_rep) );
  return *this;
}

Integer& Integer::operator %= (const unsigned long l)
{
//  if (l ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (isZero(*this)) return *this;
  mpz_tdiv_r_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  return *this;
}

Integer& Integer::operator %= (const long l)
{
//  if (l ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (isZero(*this)) return *this;
  int sgn = GMP__SGN(l);
  mpz_tdiv_r_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, GMP__ABS(l));
//   if (sgn <0) (res.gmp_rep).size = -(res.gmp_rep).size;
//   gmp_rep = (res.gmp_rep);
  if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep) );
  return *this;
}


Integer Integer::operator % (const Integer& n) const
{
//  if (isZero(n)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (isZero(*this)) return Integer::zero;
  Integer res;
  mpz_tdiv_r( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return res;
}

// long Integer::operator % (const unsigned long l) const
// {
// //  if (l ==0) {
// //    GivMathDivZero("[Integer::/]: division by zero");
// //  }
//   if (isZero(*this)) return 0L;
//   Integer Res(Integer::one);
//   mpz_tdiv_r_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&gmp_rep, l);
// //   return Integer((res.gmp_rep));
//   return long( Res );
// }
unsigned long Integer::operator % (const unsigned long l) const {
    if (isZero(*this)) return 0L;
    if (this->priv_sign()>0)
        return  mpz_tdiv_ui( (mpz_ptr)&gmp_rep, l);
    else {
        Integer Neg;
        mpz_neg( (mpz_ptr)&(Neg.gmp_rep), (mpz_ptr)&gmp_rep );
        unsigned long res = mpz_tdiv_ui( (mpz_ptr)&(Neg.gmp_rep), l);
        if (res > 0UL) return (l-res);
        else return 0UL;
    }
}

long Integer::operator % (const long l) const
{
// //  if (l ==0) {
// //    GivMathDivZero("[Integer::/]: division by zero");
//  }
    if (l>0)
        return static_cast<long>(this->operator%( static_cast<const unsigned long>(l) ) );
    else {
        long res = static_cast<long>(this->operator%( static_cast<const unsigned long>( -l ) ) );
        return -res;
    }
}

//long Integer::operator % (const long l) const
//{
//  int sgn = GMP__SGN(l);
//  long res = mpz_tdiv_ui((mpz_ptr)&gmp_rep, GMP__ABS(l));
//  if (sgn <0) return -res;
//  else return res;
//}

//Added by Dan Roche, 6-28-04
#ifdef __USE_64_bits__
unsigned long long Integer::operator % (const unsigned long long l) const
{
  if (isZero(*this)) return 0LL;
  Integer Res(Integer::one);
  Integer Divisor(l);
  mpz_tdiv_r( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&(Divisor.gmp_rep) );
  return (unsigned long long)( Res );
}

long long Integer::operator % (const long long l) const
{
  if (isZero(*this)) return 0LL;
  Integer Res(Integer::one);
  Integer Divisor(l);
  mpz_tdiv_r( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&(Divisor.gmp_rep) );
  return (long long)( Res );
}

#endif
