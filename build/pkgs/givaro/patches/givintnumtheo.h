// =================================================================== //
// Givaro : Euler's phi function
//          Primitive roots.
//          RSA scheme.
// Time-stamp: <30 Jun 04 10:59:26 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //

#ifndef _GIVARO_NUMTHEORY_
#define _GIVARO_NUMTHEORY_

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givintfactor.h"
#include "givaro/givrandom.h"

template<class RandIter = GivRandom>
class IntNumTheoDom : public IntFactorDom<RandIter> {
public:
    typedef typename IntFactorDom<RandIter>::Rep Rep;
    IntNumTheoDom(RandIter g = RandIter())
            :  IntFactorDom<RandIter>(g) {}
// =================================================================== //
// Euler's phi function
// =================================================================== //

template <template <class, class> class Container, template<class> class Alloc>
Rep& phi(Rep& res, const Container<Rep, Alloc<Rep> >& Lf, const Rep& n) const ;

Rep& phi(Rep& r, const Rep& n) const ;

// =================================================================== //
// Primitive Root
// =================================================================== //
    Rep& prim_root(Rep&, const Rep&) const ;
    Rep& prim_root(Rep&, unsigned long&, const Rep&) const ;
    Rep& prim_root_of_prime(Rep&, const Rep&) const ;
    template<class Array> Rep& prim_root_of_prime(Rep& A, const Array& Lf, const Rep& phin, const Rep& n) const ;

//  Polynomial-time generation of primitive roots
//  L is number of loops of Pollard partial factorization of n-1
//  10,000,000 gives at least 1-2^{-40} probability of success
//  [Dubrois & Dumas, Industrial-strength primitive roots]
//  Returns the probable primitive root and the probability of error.
    Rep& probable_prim_root(Rep&, double&, const Rep& n, const unsigned long L = 10000000) const;

//  Here L is computed so that the error is close to epsilon
    Rep& probable_prim_root(Rep&, double&, const Rep& n, const double epsilon) const;

    Rep& lowest_prim_root(Rep&, const Rep&) const ;
    bool is_prim_root(const Rep&, const Rep&) const ;
    Rep& order(Rep&, const Rep&, const Rep&) const ;
    bool isorder(const Rep&, const Rep&, const Rep&) const ;

// =================================================================== //
// Generalization of primitive roots for any modulus
// Primitive means maximal order
//    Primitive Element, Primitive invertible
//    Both functions coïncides except for m=8
//
// Lambda Function : maximal orbit size
//    lambda : Order of a primitive Element
//    lambda_inv : Order of an invertible Element
//    Both functions coïncides except for m=8
// =================================================================== //
    Rep& prim_inv(Rep & , const Rep&) const ;
    Rep& prim_elem(Rep & , const Rep&) const ;
private:
    Rep& prim_base(Rep & , const Rep&) const ;
    Rep& lambda_base(Rep & , const Rep&) const ;
public:
    Rep& lambda_primpow(Rep & , const Rep&, unsigned long) const ;
    Rep& lambda_inv_primpow(Rep & , const Rep&, unsigned long) const ;
    Rep& lambda(Rep & , const Rep&) const ;
    Rep& lambda_inv(Rep & , const Rep&) const ;

// =================================================================== //
// Möbius function
// =================================================================== //

template< template<class, class> class Container, template <class> class Alloc>
short mobius(const Container<Rep, Alloc<Rep> >& lpow) const ;

short mobius(const Rep& a) const;
};


#include "givaro/givintnumtheo.inl"

#endif
