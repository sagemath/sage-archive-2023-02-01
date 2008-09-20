/*
   bern_modp_util.cpp:  number-theoretic utility functions

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


#include <NTL/ZZ.h>
#include "bern_modp_util.h"


NTL_CLIENT;


namespace bernmm {


long PowerMod(long a, long ee, long n, double ninv)
{
   long x, y;

   unsigned long e;

   if (ee < 0)
      e = - ((unsigned long) ee);
   else
      e = ee;

   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = MulMod(x, y, n, ninv);
      y = MulMod(y, y, n, ninv);
      e = e >> 1;
   }

   if (ee < 0) x = InvMod(x, n);

   return x;
}



void Factorisation::helper(long k, long m)
{
   if (m == 1)
      return;

   for (long i = k + 1; i * i <= m; i++)
   {
      if (m % i == 0)
      {
         // found a factor
         factors.push_back(i);
         // remove that factor entirely
         for (m /= i; m % i == 0; m /= i);
         // recurse
         helper(i, m);
         return;
      }
   }

   // no more factors
   factors.push_back(m);
}


Factorisation::Factorisation(long n)
{
   this->n = n;
   helper(1, n);
}


PrimeTable::PrimeTable(long bound)
{
   long size = (bound - 1) / ULONG_BITS + 1;   // = ceil(bound / ULONG_BITS)
   data.resize(size);

   for (long i = 2; i * i < bound; i++)
      if (is_prime(i))
         for (long j = 2*i; j < bound; j += i)
            set(j);
}


long order(long x, long p, double pinv, const Factorisation& F)
{
   // in the loop below, m is always some multiple of the order of x
   long m = p - 1;

   // try to remove factors from m until we can't remove any more
   for (int i = 0; i < F.factors.size(); i++)
   {
      long q = F.factors[i];

      while (m % q == 0)
      {
         long mm = m / q;
         if (PowerMod(x, mm, p, pinv) != 1)
            break;
         m = mm;
      }
   }

   return m;
}



long primitive_root(long p, double pinv, const Factorisation& F)
{
   if (p == 2)
      return 1;

   long g = 2;
   for (; g < p; g++)
      if (order(g, p, pinv, F) == p - 1)
         return g;

   // no generator exists!?
   abort();
}



};    // end namespace


// end of file ================================================================
