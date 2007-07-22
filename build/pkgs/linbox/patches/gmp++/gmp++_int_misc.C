// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_misc.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_misc.C,v 1.11 2007-01-11 18:42:51 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <math.h>
#include "gmp++/gmp++.h"

//-------------------------------------------fact (unsigned long l)
Integer fact ( unsigned long l)
{
  Integer Res ;
  mpz_fac_ui( (mpz_ptr)&(Res.gmp_rep), l ) ;
  return Res ;
}

//-------------------------------------------square root
Integer& sqrt(Integer& q, const Integer &a)
{
  mpz_sqrt( (mpz_ptr)&(q.gmp_rep),
              (mpz_ptr)&(a.gmp_rep)) ;
  return q;
}

Integer& sqrtrem(Integer& q, const Integer &a, Integer& r)
{
  mpz_sqrtrem( (mpz_ptr)&(q.gmp_rep),
              (mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(a.gmp_rep)) ;
  return q;
}

Integer sqrt(const Integer &a)
{
  Integer q;
  return sqrt(q,a);
}

Integer sqrtrem(const Integer &a, Integer& r)
{
  Integer q;
  return sqrtrem(q,a,r);
}

bool root(Integer& q, const Integer &a, unsigned int n)
{
    return (bool)mpz_root ((mpz_ptr)&(q.gmp_rep),
                           (mpz_ptr)&(a.gmp_rep),
                           n);
}

void swap(Integer& a, Integer& b) {
    return mpz_swap( (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep));
}




// base p logarithm of a
long logp(const Integer& a, const Integer& p) {
    std::list< Integer > pows;
    Integer puiss = p, sq;
    do {
        pows.push_back( puiss );
    } while ( (puiss *= puiss) <= a );
    puiss = pows.back(); pows.pop_back();
    long res = (1 << pows.size());
    while (! pows.empty() ) {
        if ((sq = puiss * pows.back()) <= a) {
            puiss = sq;
            pows.pop_back();
            res += (1 << pows.size());
        } else
            pows.pop_back();
    }
    return res;
}

// approximation of the base 2 logarithm of a
// 1/log(2) being close to 1.44269504088896341
double logtwo(const Integer& a) {
  signed long int exp;
  double d = mpz_get_d_2exp( &exp, (mpz_ptr)&(a.gmp_rep) );
  return (double)exp+log(d)*1.44269504088896341;
}

//------------------------------------------GMP isprime
//     If this function returns 0, OP is definitely not prime.  If it
//     returns 1, then OP is `probably' prime.  The probability of a
//     false positive is (1/4)^r.  A reasonable value of r is 25.

Integer& nextprime(Integer& r, const Integer &p)
{
  mpz_nextprime ((mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(p.gmp_rep)) ;
  return r;
}

// Copied and adapted from mpz/nextprime.c
Integer& prevprime(Integer& r, const Integer &p)
{
   mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(p.gmp_rep), 1L );
   while( !mpz_probab_prime_p ( (mpz_ptr)&(p.gmp_rep), 5 ) )
   mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(p.gmp_rep), 1L );
   return r;
}

int probab_prime(const Integer &p)
{
  return mpz_probab_prime_p ((mpz_ptr)&(p.gmp_rep),1) ;
}

int probab_prime(const Integer &p, int r)
{
  return mpz_probab_prime_p ((mpz_ptr)&(p.gmp_rep),r) ;
}

// ==========================================================================
// Computes and returns the Jacobi and Legendre symbols (u/v) of the integers u and v.
// The algorithm used is Gmp's.
int jacobi(const Integer& u, const Integer& v)
{
  return mpz_jacobi ((mpz_ptr)&(u.gmp_rep),(mpz_ptr)&(v.gmp_rep)) ;
}

int legendre(const Integer& u, const Integer& v)
{
  return mpz_legendre ((mpz_ptr)&(u.gmp_rep),(mpz_ptr)&(v.gmp_rep)) ;
}



//--------------------------------------------Integer::operator <<   // shift left
Integer Integer::operator << (int l) const
{ return this->operator<<( (unsigned long)l ); }
Integer Integer::operator << (unsigned int l) const
{ return this->operator<<( (unsigned long)l ); }
Integer Integer::operator << (long l) const
{ return this->operator<<( (unsigned long)l ); }

Integer Integer::operator << (unsigned long l) const
{
	Integer tmp;
	mpz_mul_2exp((mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return tmp;
}


//--------------------------------------------Integer::operator >>   // shift right
Integer Integer::operator >> (int l) const
{ return this->operator>>( (unsigned long)l ); }

Integer Integer::operator >> (long l) const
{ return this->operator>>( (unsigned long)l ); }

Integer Integer::operator >> (unsigned int l) const
{ return this->operator>>( (unsigned long)l ); }

Integer Integer::operator >> (unsigned long l) const
{
	Integer tmp;
	mpz_tdiv_q_2exp( (mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return tmp;
}

//--------------------------------------------Integer::operator <<=   // shift left
Integer& Integer::operator <<= (int l)
{ return this->operator<<= ( (unsigned long)l ); }
Integer& Integer::operator <<=  (unsigned int l)
{ return this->operator<<= ( (unsigned long)l ); }
Integer& Integer::operator <<= (long l)
{ return this->operator<<= ( (unsigned long)l ); }

Integer& Integer::operator <<= (unsigned long l)
{
	mpz_mul_2exp((mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return *this;
}


//--------------------------------------------Integer::operator >>=   // shift right
Integer& Integer::operator >>= (int l)
{ return this->operator>>= ( (unsigned long)l ); }
Integer& Integer::operator >>= (long l)
{ return this->operator>>= ( (unsigned long)l ); }
Integer& Integer::operator >>= (unsigned int l)
{ return this->operator>>= ( (unsigned long)l ); }

Integer& Integer::operator >>= (unsigned long l)
{
	mpz_tdiv_q_2exp( (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return *this;
}

//------------------------------------------- Bit logic
    Integer Integer::operator^ (const Integer& a) {   // XOR
        Integer res(*this);
        return res ^= a;
    }
    Integer Integer::operator| (const Integer& a) {   // OR
        Integer res(*this);
        return res |= a;
    }
    Integer Integer::operator& (const Integer& a) {   // AND
        Integer res(*this);
        return res &= a;
    }
    Integer Integer::operator~ () const {   // 1 complement
        Integer res;
        mpz_com( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&(gmp_rep));
        return res;
    }
    Integer& Integer::operator^= (const Integer& a) {   // XOR
        mpz_xor( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
        return *this;
    }
    Integer& Integer::operator|= (const Integer& a) {   // OR
        mpz_ior( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
        return *this;
    }
    Integer& Integer::operator&= (const Integer& a) {   // AND
        mpz_and( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
        return *this;
    }


//------------------------------------------- convert method
//------------------------------------------- casting method
Integer::operator int() const {
	return mpz_get_si ( (mpz_srcptr)&gmp_rep);
}
Integer::operator unsigned int() const {
	return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
}
Integer::operator long() const {
	return mpz_get_si ( (mpz_srcptr)&gmp_rep);
}
Integer::operator unsigned long() const {
	return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
}
#ifdef __USE_64_bits__
Integer::operator unsigned long long() const {
	unsigned long low = (unsigned long)(*this);
	Integer rem;
	short cbtuli = CHAR_BIT*sizeof(unsigned long int);
	mpz_tdiv_q_2exp( (mpz_ptr)&(rem.gmp_rep), (mpz_srcptr)&(gmp_rep), cbtuli );
	unsigned long high = (unsigned long)(rem);
	unsigned long long tmp = high;
//	tmp <<= CHAR_BIT*sizeof(unsigned long int) ;
	cbtuli >>= 1;
	tmp <<= cbtuli ;
        tmp <<= cbtuli ;

	return tmp += low;
}
Integer::operator long long() const {
	unsigned long long tmp = (unsigned long long)(*this);
	return (long long)tmp;
}
#endif

Integer::operator double() const {
	return mpz_get_d ( (mpz_srcptr)&gmp_rep);
}
Integer::operator float() const {
	return (float)mpz_get_d ( (mpz_srcptr)&gmp_rep);
}

