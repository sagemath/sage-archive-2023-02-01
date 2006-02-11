#ifndef _MPZ_CLASS_EXTRAS_H
#define _MPZ_CLASS_EXTRAS_H
using namespace std;
#include <valarray>


// Allow mpz_class^(unsigned long)
mpz_class mpz_power(const mpz_class & m, unsigned long pow);


// Allow mpz_class^(unsigned long)
mpz_class operator^(const mpz_class & m, unsigned long pow);


// Allow mpz_class^(size_t)
mpz_class operator^(const mpz_class & m, size_t pow);





///////////////////////////////////
// Finds the GCD of two integers //
///////////////////////////////////

mpz_class GCD(mpz_class a, mpz_class b);


///////////////////////////////////
// Finds the LCM of two integers //
///////////////////////////////////

mpz_class LCM(mpz_class a, mpz_class b);



/////////////////////////////////////////////////////////////
// Computes the valuation of an integer m mod some prime p //
/////////////////////////////////////////////////////////////
// Note: Assumes the given p is a prime.

unsigned long Valuation(const mpz_class & m, const mpz_class & p);


////////////////////////////////////////////
// Returns the smallest non-residue mod p //
////////////////////////////////////////////

mpz_class NonResidue(const mpz_class & p);



////////////////////////////////////////////////////
// Computes the square-free part of the integer m //
////////////////////////////////////////////////////

mpz_class SquarefreePart(const mpz_class & m);



/////////////////////////////////////////////////////
// Computes the core discriminant of the integer d //
/////////////////////////////////////////////////////

mpz_class CoreDiscriminant(const mpz_class d);



//////////////////////////////////////////////////////
// Computes the Legendre symbol (a/p) at a prime p. //
//////////////////////////.///////////////////////////

int LegendreSymbol(const mpz_class & a, const mpz_class & p);


/////////////////////////////////////////
// Computes the Kronecker symbol (a/b) //
/////////////////////////////////////////

long KroneckerSymbol(const mpz_class & a, const mpz_class & b);



/////////////////////////////////////////////////////
// Computes the Hilbert symbol (a,b) at a prime p. //
/////////////////////////////////////////////////////

long HilbertSymbol(const mpz_class & a, const mpz_class & b, const mpz_class & p);


////////////////////////////////////////////////////////
// Checks if m is a square in the p-adic numbers Q_p. //
////////////////////////////////////////////////////////

bool IsPadicSquare(const mpz_class & m, const mpz_class & p);


//////////////////////////////////////////
// Returns a valarray of prime divisors //
//////////////////////////////////////////

valarray<mpz_class> PrimeDivisors(const mpz_class & m);

#include "valarray_extras.h"  // Needed for Vector_Trim() routine


#endif
