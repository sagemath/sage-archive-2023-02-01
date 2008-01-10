// =================================================================== //
// Givaro : RSA scheme.
// Time-stamp: <24 Mar 05 13:52:12 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //

#ifndef _GIVARO_RSA_
#define _GIVARO_RSA_

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givintfactor.h"
#include "givaro/givrandom.h"


    // k = 2^16 + 1, is prime
#define SIMPLE_EXPONENT (Element( (1<<16)+1 ))

// =================================================================== //
// RSA public-key cipher codes
// =================================================================== //

template<class RandIter = GivRandom>
class IntRSADom : public IntFactorDom<RandIter> {
public:
    // JGD 19.02.2003 : Should work nicely, but produces:
    //                  "implicit typename is deprecated"
    // using IntFactorDom<RandIter>::Element;
    // using IntFactorDom<RandIter>::random_generator;
    typedef typename IntFactorDom<RandIter>::Element Element;
    typedef typename IntFactorDom<RandIter>::random_generator random_generator;

private:
    Element _m, _k;
    Element _u;
    long _lm;

public:

// =================================================================== //
// Constructors
// =================================================================== //
    IntRSADom(bool fi = false, RandIter g = RandIter() ) : IntFactorDom<RandIter>(g), _fast_impl(fi) { keys_gen(IntFactorDom<RandIter>::_g, 257, 255, _m, _k, _u); _lm = log(_m,1<<(8*sizeof(unsigned char))); }
    IntRSADom(const long s, bool fi = false, RandIter g = RandIter() ) : IntFactorDom<RandIter>(g), _fast_impl(fi)  { keys_gen(IntFactorDom<RandIter>::_g, (s>>1)-1, (s>>1)+1, _m, _k, _u); _lm = log(_m,1<<(8*sizeof(unsigned char))); }
    IntRSADom(const long p, const long q, bool fi = false, RandIter g = RandIter() ) : IntFactorDom<RandIter>(g), _fast_impl(fi)  { keys_gen(IntFactorDom<RandIter>::_g, p, q, _m, _k, _u); _lm = log(_m,1<<(8*sizeof(unsigned char))); }
    IntRSADom(const Element& m, const Element& k, const Element& u) : _m(m), _k(k), _u(u), _lm(log(m,1<<(8*sizeof(unsigned char)))), _fast_impl( k == SIMPLE_EXPONENT )  {}
    IntRSADom(const Element& m, const Element& k) : _m(m), _k(k), _u(0), _lm(log(m,1<<(8*sizeof(unsigned char)))), _fast_impl( k == SIMPLE_EXPONENT )  {}

// =================================================================== //
// Accesses
// =================================================================== //
    const Element& getm() const { return _m; }
    const Element& getk() const { return _k; }
    const Element& getu() const { return _u; }

// =================================================================== //
// Text conversions
// =================================================================== //
    std::ostream& encipher(std::ostream&, std::istream&) const ;
    std::ostream& decipher(std::ostream&, std::istream&) ;

protected:
// =================================================================== //
// Keys generation
// public keys are m and k, the secret key is u.
// ciphering is computing       : x^k mod m
// deciphering is computing     : b^u mod m
// since for any x, x^(k.u) = x mod m
// =================================================================== //


// =================================================================== //
// [Strong Primes Are Easy to Find, J. Gordon, EUROCRYPT'84, LNCS 209
// =================================================================== //
    Element& strong_prime(random_generator& g, long psize, Element& p) const;

// =================================================================== //
// Here m = p*q
// p and q are prime numbers of respective sizes psize, qsize
// Moreover p-1 and q-1 have one prime factor of respective size 2/3
// since k.u = 1 mod (p-1)(q-1)
// =================================================================== //
    void keys_gen(random_generator& g, long psize, long qsize, Element& m, Element& k, Element& u) const ;

// =================================================================== //
// log[10]
// =================================================================== //
    long log(const Element& n, const long) const ;

// =================================================================== //
// Text conversions
// =================================================================== //
    std::ostream& ecriture_str(std::ostream&, const Element&) const ;
    std::ostream& ecriture_str_last(std::ostream&, const Element&) const ;
    std::ostream& ecriture_Int(std::ostream&, const Element&) const ;

public:
// =================================================================== //
// Breaking codes : finding u knowing only m an k ...
// =================================================================== //
    Element& point_break(Element& u) ;

// Fast implementation
// Means simple enciphering key, and deciphering via chinese remaindering
// WARNING: this means less security !
    bool _fast_impl;
};

#include "givaro/givintrsa.inl"

#endif
