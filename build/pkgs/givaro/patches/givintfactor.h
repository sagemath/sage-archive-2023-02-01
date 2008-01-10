// =================================================================== //
// Givaro : Prime numbers
//              Factor sets :
//              Pollard's rho method for factorization
//              Elliptic curves factorization by Lenstra
// Needs Container structures : stl ones for instance
// Time-stamp: <29 Jun 05 14:12:09 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //
#ifndef _GIVARO_FACTORISATION_H_
#define _GIVARO_FACTORISATION_H_

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givrandom.h"

// #define BOUNDARY_factor TABMAX2

#define factor_first_primes(tmp,n) (tmp = isZero(mod(tmp,n,23))?23:( isZero(mod(tmp,n,19))?19:( isZero(mod(tmp,n,17))?17:  (isZero(mod(tmp,n,2))?2:( isZero(mod(tmp,n,3))?3:( isZero(mod(tmp,n,5))?5:( isZero(mod(tmp,n,7))?7: ( isZero(mod(tmp,n,11))?11:13 ))))))))

#define factor_second_primes(tmp,n) (tmp = isZero(mod(tmp,n,31))?31:( isZero(mod(tmp,n,29))?29: ( isZero(mod(tmp,n,37))?37: ( isZero(mod(tmp,n,41))?41:( isZero(mod(tmp,n,43))?43:  ( isZero(mod(tmp,n,71))?71:( isZero(mod(tmp,n,67))?67:( isZero(mod(tmp,n,61))?61:( isZero(mod(tmp,n,59))?59: ( isZero(mod(tmp,n,53))?53:( isZero(mod(tmp,n,47))?47: ( isZero(mod(tmp,n,97))?97: ( isZero(mod(tmp,n,89))?89:( isZero(mod(tmp,n,83))?83:( isZero(mod(tmp,n,79))?79:73)))))))))))))))


// =================================================================== //
// Set or Container of divisors, factors.
// =================================================================== //

template<class RandIter = GivRandom>
class IntFactorDom : public IntPrimeDom {
private:
    // 2*3*5*7*11*13*17*19*23
    const int PROD_first_primes;
    // 29*31*37*41*43*47*53*59*61*67*71*73*79*83*89*97
    const Rep PROD_second_primes;
protected:
    RandIter _g;

public:
    typedef RandIter random_generator;

    IntFactorDom(RandIter g = RandIter()) :  IntPrimeDom(),PROD_first_primes(223092870), PROD_second_primes("10334565887047481278774629361"), _g(g) {
#ifdef __GMP_PLUSPLUS__
	    seeding( _g.seed() );
#endif
    }

        //  loops defaulted to 0 forces Pollard's factorization to
        //  be complete
    Rep& factor(Rep& r, const Rep& n, unsigned long loops = 0) const {
        if (isOne(gcd(r,n,PROD_first_primes)))
            if (isOne(gcd(r,n,PROD_second_primes))) {
#ifdef GIVARO_LENSTRA
                return Lenstra((RandIter&)_g, r, n);
#else
                return Pollard((RandIter&)_g, r, n, loops);
#endif
            } else
                return factor_second_primes(r,n);
        else
            return factor_first_primes(r,n);
    }

        //  Factors are checked for primality
    Rep& iffactorprime (Rep& r, const Rep& n, unsigned long loops = 0) const {
	if (factor(r, n, loops) != 1) {
	   if (! isprime(r,_GIVARO_ISPRIMETESTS_) ) {
		Rep nn = r; factor(r,nn, loops);
	   }
	   while (! isprime(r,_GIVARO_ISPRIMETESTS_) ) {
		Rep nn = r;
        	if (isOne(gcd(r,nn,PROD_first_primes))) {
            	   if (isOne(gcd(r,nn,PROD_second_primes))) {
                	Pollard((RandIter&)_g, r, nn, loops);
            	   } else {
                	factor_second_primes(r,nn);
		   }
		} else {
            	   factor_first_primes(r,nn);
		}
	    	if (r == nn) {
			Lenstra((RandIter&)_g, r, nn) ;
			break; // In case Lenstra fails also
		}
           }
	}
	return r;
    }

    Rep& primefactor(Rep& r, const Rep& n) const {
	while ((iffactorprime(r,n,0) == 1) && (! isprime(n, _GIVARO_ISPRIMETESTS_)) ) {}
	return r;
    }


        /// Factors with primes
        //  loops defaulted to 0 forces factorization to be complete
        //  otherwise returns if factorization is complete or not
        //  Factors are checked for primality
/* #ifndef __ECC */
/*     template< template<class> class Container> bool set */
/*         ( Container<Rep>& setint, Container<unsigned long>& setpwd,  const Rep& a, unsigned long loops = 0) const ; */
/*         /// */
/*     template< template<class> class Container> void set( Container<Rep>&,  const Rep&) const ; */
/*         /// */
/*     template< template<class> class Container> void Erathostene(Container<Rep>&, const Rep&) const ; */
/*         /// */
/*     template< template<class> class Container, template<class> class Cont2> Container<Rep>& divisors(Container<Rep>& L, const Cont2<Rep>& Lf, const Cont2<unsigned long>& Le)  const; */
/*     template< template<class> class Container> Container<Rep>& divisors(Container<Rep>&, const Rep& ) const ; */
/* #else */
    template<class Container1, class Container2> bool set
        ( Container1& setint, Container2& setpwd,  const Rep& a, unsigned long loops = 0) const ;
        ///
    template<class Container> void set( Container&,  const Rep&) const ;
        ///
    template<class Container> void Erathostene(Container&, const Rep&) const ;
        ///
    template<class Container, class Cont2, class Cont3> Container& divisors(Container& L, const Cont2& Lf, const Cont3& Le)  const;
    template<class Container> Container& divisors(Container&, const Rep& ) const ;
//#endif

        /// returns a small factor
    Rep& Erathostene(Rep&,  const Rep& p ) const ;

        // Pollard with a bound on the number of loops
        // Bound 0 is without bound
    Rep& Pollard(RandIter&, Rep&, const Rep& n, unsigned long threshold = 0) const ;
        // returns a factor by Lenstra's elliptic curves method
    Rep& Lenstra(RandIter&, Rep&, const Rep& n, const Rep& B1 = 10000000, const unsigned long curves = 30) const ;

    std::ostream& write(std::ostream& o, const Rep& n) const;
    template<class Array> std::ostream& write(std::ostream& o, Array&, const Rep& n) const;


private:
// Those are parameters for Pollard's algorithms
// Pollard_fctin : must be somewhat a "random" function in Z/nZ
// Pollard_cst can be a slight alternative for the Pfct x^2+1
#ifndef Pollard_cst
#define Pollard_cst 1
#endif

    Rep& Pollard_fctin(Rep & x, const Rep& n) const {
        mulin(x,x);
        addin(x,Pollard_cst);
        return modin(x,n);
    }

};

#include "givaro/givintfactor.inl"

#endif
