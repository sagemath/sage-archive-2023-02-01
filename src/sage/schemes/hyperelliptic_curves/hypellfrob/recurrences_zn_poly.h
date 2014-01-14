/* ============================================================================

   recurrences_zn_poly.h:  header for recurrences_zn_poly.cpp

   This file is part of hypellfrob (version 2.1.1).

   Copyright (C) 2007, 2008, David Harvey

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

============================================================================ */


#include <vector>
#include <NTL/ZZ.h>
#include <zn_poly/zn_poly.h>



namespace hypellfrob {


/*
This struct stores precomputed information that can then be used to shift
evaluation values of a polynomial F(x) of degree d.

Specifically, given the values
  F(0), F(b), F(2*b), ..., F(d*b),
the shift() method computes
  F(a), F(a + b), F(a + 2*b), ..., F(a + d*b).

PRECONDITIONS:
   d >= 0
   1, 2, ..., d + 1 are invertible
   a + i*b are invertible for -d <= i <= d

*/
struct Shifter
{
   ulong d;

   // input_twist is a vector of length d + 1.
   // The i-th entry is \prod_{0 <= j <= d, j != i} (i-j)^(-1).
   // todo: this is symmetric, so can make it d/2 + 1, but need to
   // deal with even/odd cases separately?
   ulong* input_twist;

   // output_twist is a vector of length d + 1.
   // The i-th entry is b^(-d) \prod_{0 <= j <= d} (a + (i-j)*b).
   ulong* output_twist;

   // precomputed info for performing middle product against a "kernel"
   // polynomial of degree 2d. The coefficients of "kernel" are
   // (a + k*b)^(-1) for -d <= k <= d.
   zn_array_mulmid_precomp1_t kernel_precomp;

   // Scratch space for shift(), length d + 1
   ulong* scratch;

   // zn_poly modulus object
   const zn_mod_struct* mod;

   ~Shifter();

   // Constructor (performs various precomputations)
   Shifter(ulong d, ulong a, ulong b, const zn_mod_t mod);

   // Shifts evaluation values as described above.
   // Assumes both output and input have length d + 1.
   void shift(ulong* output, const ulong* input);
};


/*
   Grrr..... I would rather use vector<ulong>, but then strictly speaking
   I can't be guaranteed that the data is stored in a simple array format,
   so here I go rolling my own....
*/
struct ulong_array
{
   ulong* data;

   ulong_array()
   {
      data = NULL;
   }

   ulong_array(ulong amount)
   {
      data = (ulong*) malloc(sizeof(ulong) * amount);
   }

   ~ulong_array()
   {
      if (data)
         free(data);
   }

   void resize(ulong amount)
   {
      if (data)
         data = (ulong*) realloc(data, sizeof(ulong) * amount);
      else
         data = (ulong*) malloc(sizeof(ulong) * amount);
   }
};



int check_params(ulong k, ulong u, const zn_mod_t mod);



struct LargeEvaluator
{
   int r;
   ulong k, u, k2, odd;
   const std::vector<std::vector<ulong> >& M0;
   const std::vector<std::vector<ulong> >& M1;
   const zn_mod_t& mod;
   Shifter* shifter;
   std::vector<ulong_array> scratch;

   LargeEvaluator(int r, ulong k, ulong u,
                  const std::vector<std::vector<ulong> >& M0,
                  const std::vector<std::vector<ulong> >& M1,
                  const zn_mod_t& mod);

   ~LargeEvaluator();

   void evaluate(int half, std::vector<ulong_array>& output, ulong offset);
   void evaluate_all(std::vector<ulong_array>& output);
};



int zn_poly_interval_products(
         std::vector<std::vector<std::vector<ulong> > >& output,
         const std::vector<std::vector<ulong> >& M0,
         const std::vector<std::vector<ulong> >& M1,
         const std::vector<NTL::ZZ>& target, const zn_mod_t& mod);


};

// ----------------------- end of file
