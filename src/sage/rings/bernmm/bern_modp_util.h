/*
   bern_modp_util.h:  number-theoretic utility functions

   Copyright (C) 2008, David Harvey

   This file is part of the bernmm package (version 1.0).

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef BERNMM_BERN_MODP_UTIL_H
#define BERNMM_BERN_MODP_UTIL_H


#include <vector>
#include <cassert>
#include "pyport.h"

namespace bernmm {


/*
   Same as NTL's PowerMod, but also accepts an _ninv_ parameter, which is the
   same as the ninv parameter for NTL's MulMod routines, i.e. should have
   ninv = 1 / ((double) n).

   (Implementation is adapted from ZZ.c in NTL 5.4.1.)
*/
long PowerMod(long a, long ee, long n, double ninv);


/*
   Represents the factorisation of an integer n into distinct prime factors.

   (Very naive implementation!)
*/
class Factorisation
{
protected:
   /*
      Finds distinct prime factors of m in the range k < p <= m.
      Assumes that m does not have any prime factors p <= k.
      Appends factors found to _factors_.
   */
   void helper(long k, long m);

public:
   // the integer
   long n;

   // the distinct factors (in increasing order)
   std::vector<long> factors;

   // initialises with given integer
   Factorisation(long n);
};



class PrimeTable
{
private:
   std::vector<long> data;   // bit-vector; 0 means prime, 1 means composite

   // read bit from index i
   inline bool get(long i) const
   {
      return (data[i / LONG_BIT] >> (i % LONG_BIT)) & 1;
   }

   // set bit at index i
   inline void set(long i)
   {
      data[i / LONG_BIT] |= (1L << (i % LONG_BIT));
   }


public:
   // initialise with primes up to given bound
   PrimeTable(long bound);

   // test whether n is prime by table lookup
   inline bool is_prime(long n) const
   {
      return !get(n);
   }

   // returns smallest prime p that is larger than n
   long next_prime(long n) const
   {
      for (n++; !is_prime(n); n++);
      return n;
   }
};



/*
   Returns 1 if n is prime.
*/
int is_prime(long n);


/*
   Returns smallest prime larger than p.
*/
long next_prime(long p);


/*
   Computes order of x mod p, given the factorisation F of p-1.
*/
long order(long x, long p, double pinv, const Factorisation& F);


/*
   Finds the smallest primitive root mod p, given the factorisation F of p-1.
*/
long primitive_root(long p, double pinv, const Factorisation& F);


};    // end namespace


#endif

// end of file ================================================================
