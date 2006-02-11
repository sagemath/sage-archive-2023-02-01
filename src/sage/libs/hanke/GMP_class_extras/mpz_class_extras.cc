#include<gmp.h>
#include<gmpxx.h>

#include "mpz_class_extras.h"

// Allow mpz_power(a^pow) to compute a^pow where a = mpz_class and pow = unsigned long
mpz_class mpz_power(const mpz_class & m, unsigned long pow) {

  mpz_t ans;
  mpz_init (ans);

  mpz_pow_ui(ans, m.get_mpz_t(), pow);


  // Put the answer in an mpz_class
  mpz_class ans2;
  ans2 = mpz_class(ans);

  // Unallocate the mpz_t
  mpz_clear(ans);


  return ans2;
}


// Allow mpz_class^(unsigned long)
mpz_class operator^(const mpz_class & m, unsigned long pow) {

  mpz_t ans;
  mpz_init (ans);

  mpz_pow_ui(ans, m.get_mpz_t(), pow);


  // Put the answer in an mpz_class
  mpz_class ans2;
  ans2 = mpz_class(ans);

  // Unallocate the mpz_t
  mpz_clear(ans);


  return ans2;
}


// Allow mpz_class^(size_t)
mpz_class operator^(const mpz_class & m, size_t pow) {

  mpz_t ans;
  mpz_init (ans);

  mpz_pow_ui(ans, m.get_mpz_t(), pow);

  // Put the answer in an mpz_class
  mpz_class ans2;
  ans2 = mpz_class(ans);

  // Unallocate the mpz_t
  mpz_clear(ans);

  return ans2;
}





///////////////////////////////////
// Finds the GCD of two integers //
///////////////////////////////////

mpz_class GCD(mpz_class a, mpz_class b) {

  mpz_t ans;
  mpz_init (ans);

  mpz_gcd(ans, a.get_mpz_t(), b.get_mpz_t());

  // Put the answer in an mpz_class
  mpz_class ans2;
  ans2 = mpz_class(ans);

  // Unallocate the mpz_t
  mpz_clear(ans);

  return ans2;
}


///////////////////////////////////
// Finds the LCM of two integers //
///////////////////////////////////

mpz_class LCM(mpz_class a, mpz_class b) {

  mpz_t ans;
  mpz_init (ans);

  mpz_lcm(ans, a.get_mpz_t(), b.get_mpz_t());

  // Put the answer in an mpz_class
  mpz_class ans2;
  ans2 = mpz_class(ans);

  // Unallocate the mpz_t
  mpz_clear(ans);

  return ans2;
}



/////////////////////////////////////////////////////////////
// Computes the valuation of an integer m mod some prime p //
/////////////////////////////////////////////////////////////
// Note: Assumes the given p is a prime.

unsigned long Valuation(const mpz_class & m, const mpz_class & p)
{
  unsigned long val = 0;
  if (m == 0) {
    printf("Error in Valuation: Zero doesn't have a valuation!");
    return (100);
  }
  else {
    mpz_class m1;
    m1 = m;
    while(m1 % p == 0) {
      val++;
      m1 = m1 / p;
    }
  }

  return val;
}



////////////////////////////////////////////
// Returns the smallest non-residue mod p //
////////////////////////////////////////////

mpz_class NonResidue(const mpz_class & p) {

  if (p>2) {
    mpz_class r;
    r = 1;
    while (r < p) {
      if (LegendreSymbol(r,p) == mpz_class(-1))
	return r;
      r++;
    }
  }
  else
    cerr << "Error in NonResidue:  p =  " << p << " < 2. " << endl;


  // Random error checking... =)
  cerr << "Error in NonResidue " << endl;
  cerr << " p = " << p << endl;
  abort();
}



////////////////////////////////////////////////////
// Computes the square-free part of the integer m //
////////////////////////////////////////////////////

mpz_class SquarefreePart(const mpz_class & m) {

  // Convert everything to mpz_t, and do it in C... =)
  mpz_t p, p2, m1;
  mpz_init (m1);
  mpz_init (p);
  mpz_init (p2);

  //  mm = m.get_mpz_t();  // This one fails. =(
  mpz_set(m1, m.get_mpz_t());
  mpz_set_ui(p, 2);
  mpz_set_ui(p2, 4);

  while (mpz_cmp(p2, m1) <= 0) {   // Check that p2 <= m
    if (mpz_divisible_p(m1, p2) == 0) {
      mpz_nextprime(p, p);   //  Set p to the next prime
      mpz_pow_ui(p2, p, 2);  //  and recompute p^2
    }
    else
      mpz_divexact(m1, m1, p2);  // Replace m1 by m1 / p^2
  }


  mpz_class m2;
  m2 = mpz_class(m1);


  // Unallocate the mpz_t variables
  mpz_clear(p);
  mpz_clear(p2);
  mpz_clear(m1);


  return m2;

}



/////////////////////////////////////////////////////
// Computes the core discriminant of the integer d //
/////////////////////////////////////////////////////

mpz_class CoreDiscriminant(const mpz_class d) {

  mpz_class t;
  t = SquarefreePart(d);

  if ((t-1) % 4 == 0)
    return t;
  else
    return 4*t;
}




//////////////////////////////////////////////////////
// Computes the Legendre symbol (a/p) at a prime p. //
//////////////////////////.///////////////////////////

int LegendreSymbol(const mpz_class & a, const mpz_class & p) {

  // To Do: Check that p is an odd prime > 0

  if (p>2)
    return mpz_legendre(a.get_mpz_t(), p.get_mpz_t());

  cout << "\n Error in LegendreSymbol: p <= 2.\n";
  abort();
}



/////////////////////////////////////////
// Computes the Kronecker symbol (a/b) //
/////////////////////////////////////////

long KroneckerSymbol(const mpz_class & a, const mpz_class & b) {
  return mpz_kronecker(a.get_mpz_t(), b.get_mpz_t());
}









/////////////////////////////////////////////////////
// Computes the Hilbert symbol (a,b) at a prime p. //
/////////////////////////////////////////////////////

long HilbertSymbol(const mpz_class & a, const mpz_class & b, const mpz_class & p) {

  unsigned long apow, bpow, apow_new, bpow_new;
  mpz_class ares, bres;
  //mpz_class neg1 = -1;

  apow = Valuation(a, p);
  bpow = Valuation(b, p);
  ares = a / (p^apow);
  bres = b / (p^bpow);

  /*
    cout << " apow = " << apow << endl;
    cout << " bpow = " << bpow << endl;
    cout << " ares = " << ares << endl;
    cout << " bres = " << bres << endl;
    cout << " LegendreSymbol(ares, p) = " << LegendreSymbol(ares, p) << endl;
    cout << " LegendreSymbol(bres, p) = " << LegendreSymbol(bres, p) << endl;
    cout << " (p-1)/2 = " << ((p-1)/2) << endl;
    //  cout << " (LegendreSymbol(ares, p)^bpow) = " << (LegendreSymbol(ares, p)^bpow) << endl;
    //  cout << " (LegendreSymbol(bres, p)^apow) = " << (LegendreSymbol(bres, p)^apow) << endl;
  */

  mpz_class parity;

  if (p != 2) {
    if (LegendreSymbol(ares, p) == 1)
      bpow_new = 0;
    else
      bpow_new = bpow;

    if (LegendreSymbol(bres, p) == 1)
      apow_new = 0;
    else
      apow_new = apow;

    parity = apow * bpow * ((p-1)/2) + bpow_new + apow_new;
  }
  else
    parity = ((ares - 1)*(bres - 1) / 4)
      + apow * ((bres*bres - 1) / 8)
      + bpow * ((ares*ares - 1) / 8);

  /*
    cout << "p = " << p << endl;
    cout << "parity = " << parity << endl;
  */

  if (parity % 2 == 0)
    return 1;
  else
    return -1;

}





////////////////////////////////////////////////////////
// Checks if m is a square in the p-adic numbers Q_p. //
////////////////////////////////////////////////////////

bool IsPadicSquare(const mpz_class & m, const mpz_class & p) {

  // Should check that p is prime...

  unsigned long v;
  mpz_class other;

  v = Valuation(m, p);
  other = m / (p^v);

  if ((v % 2) == 0) {

    if ((p == 2) && ((other % 8) == 1))
      return true;

    if ((p>2) && (KroneckerSymbol(other,p) == 1))
      return true;

  }

  return false;

  cout << "\n Error in IsPadicSquare: Should not execute this line...   p = " << p << ", m = " << m << ", and v = " << v << endl;

}



//////////////////////////////////////////
// Returns a valarray of prime divisors //
//////////////////////////////////////////

valarray<mpz_class> PrimeDivisors(const mpz_class & m) {

  valarray<mpz_class> divisors;


  mpz_t sqrt_m;
  mpz_init(sqrt_m);
  mpz_sqrt(sqrt_m, m.get_mpz_t());

  divisors.resize(m.get_ui());


  mpz_t p_mpz;
  mpz_init(p_mpz);
  mpz_set_ui(p_mpz, 2);

  mpz_class p;
  p = 2;
  size_t ptr;
  ptr = 0;

  // This is very naive... I should improve this sometime...
  while (p <= m) {
    if ((m % p) == 0) {
      divisors[ptr] = p;
      ptr++;
    }
    //       cout << " m = " << m << ", p = " << p << ", and m % p = " << (m%p) << endl;
    mpz_nextprime(p_mpz, p_mpz);
    p = mpz_class(p_mpz);
    //       cout << " Now using p = " << p << endl;
  }

  /*
  cout << " Here is the untrimmed vector: ";
  PrintV(divisors);
  */


  // Unallocate the mpz_t variables
  mpz_clear(sqrt_m);
  mpz_clear(p_mpz);


  return VectorTrim(divisors);
}
