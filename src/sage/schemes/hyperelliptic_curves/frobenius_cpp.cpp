/* ============================================================================

An implementation of the algorithm described in
    "Kedlaya's Algorithm in Larger Characteristic"

version 1.1

Copyright (C) 2007, David Harvey

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

===============================================================================

This program must be linked against Victor Shoup's NTL library. You probably
want to enable GMP in NTL for extra speed.

The main function exported from this module is frobenius(). See frobenius.h.

The matrix product evaluation algorithms implemented here follow very closely
the algorithms described in "Linear Recurrences with Polynomial Coefficients
and Application to Integer Factorization and Cartier-Manin Operator", by Alin
Bostan, Pierrick Gaudry and Eric Schost, referred to hereunder as [BGS].

*/


#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>
#include <vector>
#include <cassert>
#include "frobenius_cpp.h"


NTL_CLIENT


inline long ZZ_to_long(const ZZ& x)
{
   return trunc_long(x, NTL_BITS_PER_LONG);
}


/* ============================================================================

   Some template stuff

The matrix evaluation code is templated so that it can work over either ZZ_p
or zz_p. There are several template parameters, which can have two settings:

   SCALAR:       ZZ_p           zz_p
   POLY:         ZZ_pX          zz_pX
   VECTOR:       vec_ZZ_p       vec_zz_p
   MATRIX:       mat_ZZ_p       mat_zz_p
   POLYMODULUS:  ZZ_pXModulus   zz_pXModulus
   FFTREP:       FFTRep         fftRep

For the most part NTL uses the same function names for both columns, which
makes life easy, but there's a few I needed to define explicitly:

   to_scalar() convert a ZZ or int into a SCALAR
   forward_fft() runs a forward FFT (either ToFFTRep() or TofftRep())
   inverse_fft() runs an inverse FFT (either FromFFTRep() or FromfftRep())

============================================================================ */



// to_scalar(ZZ)

template <typename SCALAR> SCALAR to_scalar(const ZZ& input);

template<> inline ZZ_p to_scalar<ZZ_p>(const ZZ& input)
{
   return to_ZZ_p(input);
}

template<> inline zz_p to_scalar<zz_p>(const ZZ& input)
{
   return to_zz_p(input);
}


// to_scalar(int)

template <typename SCALAR> SCALAR to_scalar(int input);

template<> inline ZZ_p to_scalar<ZZ_p>(int input)
{
   return to_ZZ_p(input);
}

template<> inline zz_p to_scalar<zz_p>(int input)
{
   return to_zz_p(input);
}


// forward_fft

template <typename POLY, typename FFTREP>
void forward_fft(FFTREP& y, const POLY& x, long k, long lo, long hi);

template<> inline void
forward_fft<ZZ_pX, FFTRep>(FFTRep& y, const ZZ_pX& x, long k, long lo, long hi)
{
   ToFFTRep(y, x, k, lo, hi);
}

template<> inline void
forward_fft<zz_pX, fftRep>(fftRep& y, const zz_pX& x, long k, long lo, long hi)
{
   TofftRep(y, x, k, lo, hi);
}


// inverse_fft

template <typename POLY, typename FFTREP>
void inverse_fft(POLY& x, FFTREP& y, long lo, long hi);

template<> inline void
inverse_fft<ZZ_pX, FFTRep>(ZZ_pX& x, FFTRep& y, long lo, long hi)
{
   FromFFTRep(x, y, lo, hi);
}

template<> inline void
inverse_fft<zz_pX, fftRep>(zz_pX& x, fftRep& y, long lo, long hi)
{
   FromfftRep(x, y, lo, hi);
}



/* ============================================================================

   Dyadic evaluation stuff

This section essentially implements Theorem 8 of [BGS], for the particular
case of the parameters that are used in Theorem 15.

============================================================================ */


/*
Assume that f has degree d = 2^n and that g has degree 2d.

This function computes a polynomial h whose x^d through x^{2d} coefficients
(inclusive) are the same as those of f*g. The bottom d coefficients of h
will be junk.

The parameter g_fft should be the precomputed length 2d FFT of g.

The algorithm in this function is based on the paper "The middle product
algorithm" by Hanrot, Quercia and Zimmermann. (Many thanks to Victor Shoup
for writing such wonderfully modular FFT code. It really made my day.)

*/
template <typename SCALAR, typename POLY, typename FFTREP>
void middle_product(POLY& h, const POLY& f,
                    const POLY& g, const FFTREP& g_fft, int n)
{
   int d = 1 << n;
   h.rep.SetLength(2*d + 1);

   FFTREP f_fft(INIT_SIZE, n+1);

   // Compute length 2d cyclic convolutions of f and g, letting the top
   // third of f*g "wrap around" to the bottom half of the output.
   forward_fft<POLY, FFTREP>(f_fft, f, n+1, 0, 2*d);
   mul(f_fft, f_fft, g_fft);
   inverse_fft<POLY, FFTREP>(h, f_fft, 0, 2*d);

   // Need to correct for the x^{2d} term of g which got wrapped around to the
   // constant term.
   h.rep[d] -= g.rep[2*d] * f.rep[d];

   // Now h contains terms x^d through x^{2d-1} of f*g.
   // To finish off, we just need the x^{2d} term.
   SCALAR temp;
   SCALAR& sum = h.rep[2*d];
   sum = 0;
   for (int i = 0; i <= d; i++)
   {
      mul(temp, f.rep[i], g.rep[2*d-i]);
      add(sum, sum, temp);
   }
}



/*
This struct stores precomputed information that can then be used to shift
evaluation values of a polynomial F(x) of degree d = 2^n.

Specifically, given the values
  F(0), F(b), F(2*b), ..., F(d*b),
the shift() method computes
  F(a), F(a + b), F(a + 2*b), ..., F(a + d*b).

PRECONDITIONS:
   n >= 1
   1, 2, ..., d + 1 are invertible
   a + i*b are invertible for -d <= i <= d

*/
template <typename SCALAR, typename POLY, typename VECTOR, typename FFTREP>
struct DyadicShifter
{
   int d, n;

   // input_twist is a vector of length d/2 + 1.
   // The i-th entry is \prod_{0 <= j <= d, j != i} (i-j)^(-1).
   VECTOR input_twist;

   // output_twist is a vector of length d+1.
   // The i-th entry is b^(-d) \prod_{0 <= j <= d} (a + (i-j)*b).
   VECTOR output_twist;

   // kernel is a polynomial of degree 2d.
   // The coefficients are (a + k*b)^(-1) for -d <= k <= d.
   // We also store its length 2d FFT.
   POLY kernel;
   FFTREP kernel_fft;

   // Polynomials for scratch space in shift()
   POLY scratch, scratch2;

   // Constructor (performs various precomputations)
   DyadicShifter(int n, const SCALAR& a, const SCALAR& b)
   {
      assert(n >= 1);
      this->n = n;
      d = 1 << n;

      // ------------------------ compute input_twist -------------------------

      input_twist.SetLength(d/2 + 1);

      // prod = (d!)^(-1)
      SCALAR prod;
      prod = 1;
      for (int i = 2; i <= d; i++)
         mul(prod, prod, i);
      prod = 1 / prod;

      // input_twist[i] = ((d-i)!)^(-1)
      input_twist[0] = prod;
      for (int i = 1; i <= d/2; i++)
         mul(input_twist[i], input_twist[i-1], d-(i-1));

      // input_twist[i] = ((d-i)!*i!)^(-1)
      prod = input_twist[d/2];
      for (int i = d/2; i >= 0; i--)
      {
         mul(input_twist[i], input_twist[i], prod);
         mul(prod, prod, i);
      }

      // input_twist[i] = \prod_{0 <= j <= d, j != i} (i-j)^(-1)    :-)
      for (int i = 1; i <= d/2; i += 2)
         NTL::negate(input_twist[i], input_twist[i]);

      // ----------------- compute output_twist and kernel --------------------

      output_twist.SetLength(d+1);

      // c[i] = c_i = a + (i-d)*b     for 0 <= i <= 2d
      VECTOR c;
      c.SetLength(2*d+1);
      c[0] = a - d*b;
      for (int i = 1; i <= 2*d; i++)
         add(c[i], c[i-1], b);

      // accum[i] = c_0 * c_1 * ... * c_i    for 0 <= i <= 2d
      VECTOR accum;
      accum.SetLength(2*d+1);
      accum[0] = c[0];
      for (int i = 1; i <= 2*d; i++)
         mul(accum[i], accum[i-1], c[i]);

      // accum_inv[i] = (c_0 * c_1 * ... * c_i)^(-1)    for 0 <= i <= 2d
      VECTOR accum_inv;
      accum_inv.SetLength(2*d+1);
      accum_inv[2*d] = 1 / accum[2*d];
      for (int i = 2*d-1; i >= 0; i--)
         mul(accum_inv[i], accum_inv[i+1], c[i+1]);

      // kernel[i] = (c_i)^(-1)    for 0 <= i <= 2d
      kernel.rep.SetLength(2*d+1);
      kernel.rep[0] = accum_inv[0];
      for (int i = 1; i <= 2*d; i++)
         mul(kernel.rep[i], accum_inv[i], accum[i-1]);

      // precompute transform of kernel
      forward_fft<POLY, FFTREP>(kernel_fft, kernel, n+1, 0, 2*d);

      // output_twist[i] = b^{-d} * c_i * c_{i+1} * ... * c_{i+d}
      // for 0 <= i <= d
      SCALAR factor = power(b, -d);
      SCALAR temp;
      output_twist.SetLength(d+1);
      output_twist[0] = factor * accum[d];
      for (int i = 1; i <= d; i++)
      {
         mul(temp, factor, accum[i+d]);
         mul(output_twist[i], temp, accum_inv[i-1]);
      }
   }


   // Shifts evaluation values as described above.
   // Assumes both output and input have length d + 1.
   void shift(VECTOR& output, const VECTOR& input)
   {
      assert(input.length() == d+1);
      assert(output.length() == d+1);

      // multiply inputs pointwise by input_twist
      scratch.rep.SetLength(d+1);
      for (int i = 0; i <= d/2; i++)
         mul(scratch.rep[i], input[i], input_twist[i]);
      for (int i = 1; i <= d/2; i++)
         mul(scratch.rep[i+d/2], input[i+d/2], input_twist[d/2-i]);

      // do polynomial multiplication
      middle_product<SCALAR, POLY, FFTREP>(scratch2, scratch, kernel,
                                           kernel_fft, n);

      // multiply outputs pointwise by output_twist
      for (int i = 0; i <= d; i++)
         mul(output[i], scratch2.rep[i+d], output_twist[i]);
   }
};



/*
Let M0 and M1 be square matrices of size n*n. Let M(x) = M0 + x*M1; this is a
matrix of linear polys in x. Let P(x) = M(x+1) M(x+2) ... M(x+2^s); this is a
matrix of polynomials of degree 2^s. This function computes the values
    P(a), P(a + 2^t), P(a + 2*2^t), ..., P(a + 2^s*2^t).

The output array should have length n^2. Each entry should be a vector of
length 2^s+1, pre-initialised to all zeroes. The (y*n + x)-th vector will be
the values of the (y, x) entries of the above list of matrices. (This data
format is optimised for the case that 2^s+1 is much larger than n.)

PRECONDITIONS:
   0 <= s <= t
   2, 3, ..., 2^t + 1 must be invertible

*/
template <typename SCALAR, typename POLY, typename VECTOR,
          typename MATRIX, typename FFTREP>
void dyadic_evaluation(vector<VECTOR>& output,
                       const MATRIX& M0, const MATRIX& M1,
                       int s, int t, const SCALAR& a)
{
   int n = M0.NumRows();

   // base cases; just evaluate naively
   if (s <= 1)
   {
      MATRIX X[3];

      if (s == 0)
      {
         X[0] = M0 + (a+1) * M1;
         X[1] = M0 + (a+1 + (1 << t)) * M1;
      }
      else
      {
         for (int i = 0; i <= 2; i++)
            X[i] = (M0 + (a+1 + (i << t)) * M1) * (M0 + (a+2 + (i << t)) * M1);
      }

      for (int x = 0; x < n; x++)
         for (int y = 0; y < n; y++)
            for (int i = 0; i < output[0].length(); i++)
               output[y*n + x][i] = X[i][y][x];

      return;
   }

   // General case.
   // Let Q(x) = M(x+1) M(x+2) ... M(x+2^(s-1)).

   // Recursively compute Q(a), Q(a + 2^t), ..., Q(a + 2^(s-1)*2^t).
   vector<VECTOR> X(n*n);
   for (int i = 0; i < n*n; i++)
      X[i].SetLength((1 << (s-1)) + 1);
   dyadic_evaluation<SCALAR, POLY, VECTOR, MATRIX, FFTREP>
                      (X, M0, M1, s-1, t, a);

   // Do precomputations for shifting by 2^(s-1) and by (2^(s-1)+1)*2^t
   SCALAR c, b;
   c = 1 << (s-1);
   b = 1 << t;
   DyadicShifter<SCALAR, POLY, VECTOR, FFTREP> shifter1(s-1, c, b);
   DyadicShifter<SCALAR, POLY, VECTOR, FFTREP> shifter2(s-1, (c + 1) * b, b);

   // Shift by 2^(s-1) to obtain
   // Q(a + 2^(s-1)), Q(a + 2^t + 2^(s-1)), ..., Q(a + 2^(s-1)*2^t + 2^(s-1))
   vector<VECTOR> Y(n*n);
   for (int i = 0; i < n*n; i++)
   {
      Y[i].SetLength((1 << (s-1)) + 1);
      shifter1.shift(Y[i], X[i]);
   }

   // Multiply matrices to obtain
   // P(a), P(a + 2^t), ..., P(a + 2^(s-1)*2^t).
   SCALAR temp;
   for (int i = 0; i <= (1 << (s-1)); i++)
      for (int x = 0; x < n; x++)
         for (int y = 0; y < n; y++)
            for (int z = 0; z < n; z++)
            {
               mul(temp, X[y*n + z][i], Y[z*n + x][i]);
               output[y*n + x][i] += temp;
            }

   // Shift original sequence by (2^(s-1)+1)*2^t to obtain
   // Q(a + (2^(s-1)+1)*2^t), Q(a + (2^(s-1)+2)*2^t), ..., Q(a + (2^s+1)*2^t).
   for (int i = 0; i < n*n; i++)
      shifter2.shift(Y[i], X[i]);

   // Shift again by 2^(s-1) to obtain
   // Q(a + (2^(s-1)+1)*2^t + 2^(s-1)), Q(a + (2^(s-1)+2)*2^t + 2^(s-1)), ...,
   //                                             Q(a + (2^s+1)*2^t + 2^(s-1)).
   for (int i = 0; i < n*n; i++)
      shifter1.shift(X[i], Y[i]);

   // Multiply matrices to obtain
   // P(a + (2^(s-1)+1)*2^t), P(a + (2^(s-1)+2)*2^t), ..., P(a + (2^s+1)*2^t).
   // (we throw out the last one since it's surplus to requirements)
   for (int i = 0; i < (1 << (s-1)); i++)
      for (int x = 0; x < n; x++)
         for (int y = 0; y < n; y++)
            for (int z = 0; z < n; z++)
            {
               mul(temp, Y[y*n + z][i], X[z*n + x][i]);
               output[y*n + x][i + (1 << (s-1)) + 1] += temp;
            }
}



/* ============================================================================

   General evaluation stuff

This section essentially implements Corollary 10 of [BGS].

============================================================================ */


/*
This struct stores the product tree associated to a vector a[0], ..., a[n-1].

The top node stores the polynomial product
     (x - a[0]) ... (x - a[n-1]).
The two children nodes store
     (x - a[0]) ... (x - a[m-1])
and
     (x - a[m]) ... (x - a[n-1])
where m = floor(n/2). This continues recursively until we reach n = 1,
in which case just the polynomial x - a[0] is stored, and no children.

*/
template <typename SCALAR, typename POLY, typename VECTOR>
struct ProductTree
{
   // polynomial product stored at this node
   POLY poly;

   // children for left and right halves, if deg(poly) > 1
   ProductTree* child1;
   ProductTree* child2;

   // These are temp polys used by the Evaluator and Interpolator classes.
   // It's not very hygienic to keep them here... but it makes things more
   // efficient, because we need two temps for each node, and this prevent
   // unnecessary reallocations. (The lengths will be the same on repeated
   // calls to evaluate() and interpolate().)
   POLY scratch1, scratch2;

   // Constructs product tree for the supplied vector.
   ProductTree(const VECTOR& points)
   {
      build(points, 0, points.length());
   }

   ProductTree(const VECTOR& points, int start, int end)
   {
      build(points, start, end);
   }

   // Constructs product tree recursively for the subset [start, end) of
   // the supplied vector.
   void build(const VECTOR& points, int start, int end)
   {
      assert(end - start >= 1);
      assert(start >= 0);
      assert(end <= points.length());

      if (end - start == 1)
      {
         SetCoeff(poly, 1, 1);
         SetCoeff(poly, 0, -points[start]);
      }
      else
      {
         int m = (end - start) / 2;
         child1 = new ProductTree(points, start, start + m);
         child2 = new ProductTree(points, start + m, end);
         mul(poly, child1->poly, child2->poly);
      }
   }

   ~ProductTree()
   {
      if (deg(poly) > 1)
      {
         delete child1;
         delete child2;
      }
   }
};



/*
Given a list of evaluation points a[0], ..., a[n-1], this struct stores some
precomputed information to permit evaluating an arbitrary polynomial at those
points.
*/
template <typename SCALAR, typename POLY,
          typename POLYMODULUS, typename VECTOR>
struct Evaluator
{
   // The product tree for the evaluation points
   ProductTree<SCALAR, POLY, VECTOR>* tree;

   // A list of NTL ZZ_pXModulus/zz_pXModulus objects corresponding to the
   // polynomials in the product tree, in the order that they get used as the
   // tree is traversed in recursive_evaluate().
   vector<POLYMODULUS> moduli;

   // Constructs evaluator object for the given list of evaluation points
   Evaluator(const VECTOR& points)
   {
      assert(points.length() >= 1);
      tree = new ProductTree<SCALAR, POLY, VECTOR>(points);
      moduli.reserve(2*points.length());
      build(tree);
      assert(moduli.size() <= 2*points.length());
   }

   // Compute modulus objects for each polynomial under the supplied node of
   // the product tree; appends them in traversal order to "moduli".
   void build(const ProductTree<SCALAR, POLY, VECTOR>* node)
   {
      if (deg(node->poly) > 1)
      {
         moduli.push_back(POLYMODULUS(node->poly));
         build(node->child1);
         build(node->child2);
      }
   }

   ~Evaluator()
   {
      delete tree;
   }

   // Evaluates the input polynomial at the evaluation points, writes the
   // results to output. The output array must have the correct length.
   void evaluate(VECTOR& output, const POLY& input)
   {
      recursive_evaluate(output, input, tree, 0, 0);
   }

   // Evaluates the input polynomial at the subset [start, end) of the
   // evaluation points, which should correspond to the supplied product tree
   // node. (The length of the interval is implied by the degree of the poly
   // at that node.) Writes the output to the subset [start, end) of the
   // output array. The index parameter indicates which modulus in "moduli"
   // to use for this node of the tree. The return value is the index for
   // the modulus that should be used immediately after this call.
   int recursive_evaluate(VECTOR& output, const POLY& input,
                          ProductTree<SCALAR, POLY, VECTOR>* node,
                          int start, int index)
   {
      if (deg(node->poly) == 1)
      {
         eval(output[start], input, -coeff(node->poly, 0));
      }
      else
      {
         rem(node->scratch1, input, moduli[index++]);
         index = recursive_evaluate(output, node->scratch1, node->child1,
                                    start, index);
         index = recursive_evaluate(output, node->scratch1, node->child2,
                                    start + deg(node->child1->poly), index);
      }
      return index;
   }
};


/*
Given an integer L >= 1, this struct does some precomputations to permit
interpolating a polynomial whose values at 0, 1, ..., L are known.

PRECONDITIONS:
   1, 2, ..., L must be invertible.

*/
template <typename SCALAR, typename POLY, typename VECTOR>
struct Interpolator
{
   ProductTree<SCALAR, POLY, VECTOR>* tree;
   int L;

   // input_twist is a vector of length L+1.
   // The i-th entry is \prod_{0 <= j <= L, j != i} (i-j)^(-1).
   VECTOR input_twist;

   // vector of length L+1, used in interpolate()
   VECTOR temp;

   // Performs various precomputations for the given L.
   Interpolator(int L)
   {
      this->L = L;
      temp.SetLength(L+1);

      // Build a product tree for the evaluation points
      for (int i = 0; i <= L; i++)
         temp[i] = i;
      tree = new ProductTree<SCALAR, POLY, VECTOR>(temp);

      // prod = (L!)^(-1)
      SCALAR prod;
      prod = 1;
      for (int i = 2; i <= L; i++)
         mul(prod, prod, i);
      prod = 1 / prod;

      // input_twist[i] = (i!)^(-1),   0 <= i <= L
      input_twist.SetLength(L+1);
      input_twist[L] = prod;
      for (int i = L; i >= 1; i--)
         mul(input_twist[i-1], input_twist[i], i);

      // input_twist[i] = \prod_{0 <= j <= L, j != i} (i-j)^(-1).
      for (int i = 0; i <= L/2; i++)
      {
         mul(input_twist[i], input_twist[i], input_twist[L-i]);
         input_twist[L-i] = input_twist[i];
      }
      for (int i = L-1; i >= 0; i -= 2)
         NTL::negate(input_twist[i], input_twist[i]);
   }

   ~Interpolator()
   {
      delete tree;
   }


   // Returns the polynomial
   //    \sum_{i=start}^{end-1} values[i] * (x-start) (x-start+1) ... (x-end-1)
   // where [start, end) is the interval associated to the supplied product
   // tree node, and where the (x-i) term is omitted in each product.
   void combine(POLY& output, const VECTOR& values,
                ProductTree<SCALAR, POLY, VECTOR>* node, int start)
   {
      if (deg(node->poly) == 1)
      {
         // base case
         clear(output);
         SetCoeff(output, 0, values[start]);
      }
      else
      {
         // recursively build up from two halves
         // i.e. if f1, f2 are the results of "combine" for the two halves,
         // and if p1, p2 are the associated product tree polys, we compute
         // f1*p2 + f2*p1

         combine(node->scratch1, values, node->child1, start);
         mul(output, node->scratch1, node->child2->poly);

         combine(node->scratch1, values, node->child2,
                 start + deg(node->child1->poly));
         mul(node->scratch2, node->scratch1, node->child1->poly);

         add(output, output, node->scratch2);
      }
   }

   // Returns a polynomial F(x) of degree at most L such that F(i) = values[i]
   // for each 0 <= i <= L.
   void interpolate(POLY& output, const VECTOR& values)
   {
      assert(values.length() == L+1);

      // multiply input values pointwise by input_twist; this corrects for
      // the factor (i-0) (i-1) ... (i-L) (where the i-i factor is omitted).
      for (int i = 0; i <= L; i++)
         mul(temp[i], values[i], input_twist[i]);

      // do the interpolation
      combine(output, temp, tree, 0);
   }
};



/* ============================================================================

   Matrix products over arbitrary, relatively short intervals

This section implements something similar to steps 1, 2, ... and the final
refining step of Theorem 15 of [BGS].

============================================================================ */


/*
Let M0 and M1 be matrices of constants. This function evaluates
   M(x) = M0 + x*M1
at x = a.

The output matrix must already have the correct dimensions.

*/
template <typename SCALAR, typename MATRIX>
void eval_matrix(MATRIX& output, const MATRIX& M0, const MATRIX& M1,
                 const SCALAR& a)
{
   int n = M0.NumRows();
   for (int x = 0; x < n; x++)
      for (int y = 0; y < n; y++)
      {
         mul(output[x][y], a, M1[x][y]);
         add(output[x][y], output[x][y], M0[x][y]);
      }
}



/*
Let M(x) be the matrix M0 + x*M1; this is a matrix of linear polys in x.
Let M(a, b) = M(a + 1) M(a + 2) ... M(b). This function evaluates the products
M(a[i], b[i]) for some sequence of intervals
  a[0] < b[0] <= a[1] < b[1] <= ... <= a[n-1] < b[n-1].

The intervals are supplied in "target", simply as the list
  a[0], b[0], a[1], b[1], ...

NOTE:
   This is used as a subroutine of large_interval_products() to handle the
   smaller "refining" subintervals. Its asymptotic complexity theoretically has
   an extra logarithmic factor over that of large_interval_products().

PRECONDITIONS:
   Let d = sum of lengths of intervals. Then 2, 3, ... 1 + floor(sqrt(d)) must
   all be invertible.

*/
template <typename SCALAR, typename POLY, typename POLYMODULUS,
          typename VECTOR, typename MATRIX>
void small_interval_products(vector<MATRIX>& output,
                             const MATRIX& M0, const MATRIX& M1,
                             const vector<ZZ>& target)
{
   output.clear();

   if (target.size() == 0)
      return;

   int dim = M0.NumRows();
   int num_intervals = target.size() / 2;

   // Determine maximum target interval length
   int max_length = -1;
   for (int i = 0; i < target.size(); i += 2)
   {
      int temp = ZZ_to_long(target[i+1] - target[i]);
      if (temp > max_length)
         max_length = temp;
   }

   // Select an appropriate length for the matrix products we'll use
   int L, max_eval_points;
   if (max_length > 2*num_intervals)
   {
      // The intervals are still pretty long relative to the number of
      // intervals, so we're only going to do a single multipoint
      // evaluation.
      L = 1 + ZZ_to_long(SqrRoot(num_intervals * to_ZZ(max_length)));
      max_eval_points = L;
   }
   else
   {
      // The intervals are getting pretty short, so we probably will need
      // to do several shorter multipoint evaluations.
      L = 1 + max_length/2;
      max_eval_points = num_intervals;
   }

   // =========================================================================
   // Step 1: compute entries of M(X, X+L) as polynomials in X.

   vector<POLY> polys(dim*dim);
   {
      // left_accum[i]  = M(L-i-1, L)  for 0 <= i <= L-1
      // right_accum[i] = M(L, L+i+1)  for 0 <= i <= L-1
      vector<MATRIX> left_accum(L), right_accum(L);

      MATRIX temp;
      temp.SetDims(dim, dim);

      left_accum[0].SetDims(dim, dim);
      eval_matrix<SCALAR, MATRIX>(left_accum[0], M0, M1, to_scalar<SCALAR>(L));
      for (int i = L-1; i >= 1; i--)
      {
         eval_matrix<SCALAR, MATRIX>(temp, M0, M1, to_scalar<SCALAR>(i));
         mul(left_accum[L-i], temp, left_accum[L-i-1]);
      }

      right_accum[0].SetDims(dim, dim);
      eval_matrix<SCALAR, MATRIX>(right_accum[0], M0, M1,
                                  to_scalar<SCALAR>(L+1));
      for (int i = 1; i <= L-1; i++)
      {
         eval_matrix<SCALAR, MATRIX>(temp, M0, M1, to_scalar<SCALAR>(L+1+i));
         mul(right_accum[i], right_accum[i-1], temp);
      }

      // Use left_accum and right_accum to compute:
      // initial[i] = M(i, L+i)    for 0 <= i <= L
      // i.e. initial[i] are the values of M(X, X+L) at X = 0, 1, ..., L.
      vector<MATRIX> initial(L+1);
      initial[0] = left_accum.back();
      initial[L] = right_accum.back();
      for (int i = 1; i <= L-1; i++)
         mul(initial[i], left_accum[L-1-i], right_accum[i-1]);

      // Now interpolate entries of initial[i] to get entries of M(X, X+L)
      // as polynomials of degree L.
      Interpolator<SCALAR, POLY, VECTOR> interpolator(L);
      VECTOR values;
      values.SetLength(L+1);
      for (int x = 0; x < dim; x++)
         for (int y = 0; y < dim; y++)
         {
            for (int j = 0; j <= L; j++)
               values[j] = initial[j][y][x];
            interpolator.interpolate(polys[y*dim + x], values);
         }
   }

   // =========================================================================
   // Step 2: decompose intervals into subintervals of length L which we'll
   // attack by direct multipoint evaluation, plus leftover pieces that we'll
   // handle with a recursive call to small_interval_products().

   // eval_points holds all the values of X for which we want to
   // evaluate M(X, X+L)
   VECTOR eval_points;
   eval_points.SetMaxLength(max_eval_points);

   // leftover_target is the list of leftover intervals that we're going to
   // later do recursively
   vector<ZZ> leftover_target;
   leftover_target.reserve(target.size());

   ZZ current, next;
   for (int i = 0; i < target.size(); i += 2)
   {
      current = target[i];
      next = current + L;
      while (next <= target[i+1])
      {
         // [current, next) fits inside this interval, so peel it off into
         // eval_points
         append(eval_points, to_scalar<SCALAR>(current));
         swap(current, next);
         next = current + L;
      }
      if (current < target[i+1])
      {
         // the rest of this interval is too short to handle with M(X, X+L),
         // so put it in the leftover bin
         leftover_target.push_back(current);
         leftover_target.push_back(target[i+1]);
      }
   }

   // =========================================================================
   // Step 3: recursively handle leftover pieces

   // leftover_matrices[i] holds the matrix for leftover interval #i
   vector<MATRIX> leftover_matrices;
   small_interval_products<SCALAR, POLY, POLYMODULUS, VECTOR, MATRIX>
                    (leftover_matrices, M0, M1, leftover_target);

   // =========================================================================
   // Step 4: evaluate M(X, X+L) at each of the evaluation points. We do this
   // by breaking up the list of evaluation points into blocks of length at
   // most L+1, and using multipoint evaluation on each block.

   // main_matrices[i] will hold M(X, X+L) for the i-th evaluation point X.
   vector<MATRIX> main_matrices(eval_points.length());
   for (int i = 0; i < main_matrices.size(); i++)
      main_matrices[i].SetDims(dim, dim);

   VECTOR block, values;
   block.SetMaxLength(L+1);
   values.SetMaxLength(L+1);

   // for each block...
   for (int i = 0; i < eval_points.length(); i += (L+1))
   {
      // determine length of this block, which is at most L+1
      int length = eval_points.length() - i;
      if (length >= (L+1))
         length = (L+1);
      block.SetLength(length);

      // construct Evaluator object for evaluating at these points
      for (int j = 0; j < length; j++)
         block[j] = eval_points[i+j];
      Evaluator<SCALAR, POLY, POLYMODULUS, VECTOR> evaluator(block);

      // evaluate each entry of M(X, X+L) at those points
      for (int x = 0; x < dim; x++)
         for (int y = 0; y < dim; y++)
         {
            evaluator.evaluate(values, polys[y*dim + x]);
            for (int k = 0; k < length; k++)
               main_matrices[i+k][y][x] = values[k];
         }
   }

   // =========================================================================
   // Step 5: merge together the matrices obtained from the multipoint
   // evaluation step and the recursive leftover interval step.

   output.clear();
   output.resize(target.size() / 2);
   for (int i = 0; i < target.size()/2; i++)
      output[i].SetDims(dim, dim);

   int main_index = 0;       // index into main_matrices
   int leftover_index = 0;   // index into leftover_matrices

   MATRIX temp;
   temp.SetDims(dim, dim);

   for (int i = 0; i < target.size(); i += 2)
   {
      current = target[i];
      next = current + L;
      ident(output[i/2], dim);

      while (next <= target[i+1])
      {
         // merge in a matrix from multipoint evaluation step
         mul(temp, output[i/2], main_matrices[main_index++]);
         swap(temp, output[i/2]);
         swap(current, next);
         next = current + L;
      }
      if (current < target[i+1])
      {
         // merge in a matrix from leftover interval step
         mul(temp, output[i/2], leftover_matrices[leftover_index++]);
         swap(temp, output[i/2]);
      }
   }
}


/* ============================================================================

   Matrix products over arbitrary, long intervals

This section implements an algorithm similar to Theorem 15 of [BGS].

============================================================================ */


/*
Let M(x) be the matrix M0 + x*M1; this is a matrix of linear polys in x.
Let M(a, b) = M(a + 1) M(a + 2) ... M(b). This function evaluates the products
M(a[i], b[i]) for some sequence of intervals
  a[0] < b[0] <= a[1] < b[1] <= ... <= a[n-1] < b[n-1].

The intervals are supplied in "target", simply as the list
  a[0], b[0], a[1], b[1], ...

NOTE:
   This algorithm works best if the intervals are very long and don't have
   much space between them. The case where the gaps are relatively large is
   best handled by small_interval_products().

PRECONDITIONS:
   Let d = b[n-1] - a[0]. Then 2, 3, ... 1 + floor(sqrt(d)) must all be
   invertible.

*/
template <typename SCALAR, typename POLY, typename POLYMODULUS,
          typename VECTOR, typename MATRIX, typename FFTREP>
void large_interval_products(vector<MATRIX>& output,
                             const MATRIX& M0, const MATRIX& M1,
                             const vector<ZZ>& target)
{
   assert(target.size() % 2 == 0);
   output.resize(target.size() / 2);

   int dim = M0.NumRows();
   assert(dim == M0.NumCols());
   assert(dim == M1.NumRows());
   assert(dim == M1.NumCols());

   // =========================================================================
   // Step 0: get as many intervals as possible using dyadic_evaluation().

   // step0_matrix[i] is the transition matrix between step0_index[2*i]
   // and step0_index[2*i+1].
   vector<MATRIX> step0_matrix;
   vector<ZZ> step0_index;
   // preallocate the maximum number of matrices that could arise (plus safety)
   int reserve_size = target.size() +
                      4*NumBits(target.back() - target.front());
   step0_matrix.reserve(reserve_size);
   step0_index.reserve(2 * reserve_size);

   ZZ current_index = target.front();
   int next_target = 0;   // index into "target" array

   // This flag indicates whether the last entry of step0_matrix is
   // still accumulating matrices (in which the right endpoint of the
   // corresponding interval hasn't been written to step0_index yet).
   int active = 0;

   MATRIX temp_mat;
   temp_mat.SetDims(dim, dim);

   while (current_index < target.back() - 3)
   {
      // find largest t such that 2^t*(2^t + 1) <= remaining distance to go
      ZZ remaining = target.back() - current_index;
      int t = 0;
      while ((to_ZZ(1) << (2*t)) + (1 << t) <= remaining)
         t++;
      t--;

      // evaluate matrices for 2^t+1 intervals of length 2^t
      vector<VECTOR> dyadic_output(dim*dim);
      for (int i = 0; i < dim*dim; i++)
         dyadic_output[i].SetLength((1 << t) + 1);
      dyadic_evaluation<SCALAR, POLY, VECTOR, MATRIX, FFTREP>
               (dyadic_output, M0, M1, t, t, to_scalar<SCALAR>(current_index));

      // Walk through the intervals we just computed. Find maximal subsequences
      // of intervals none of which contain any target endpoints. Merge them
      // together (by multiplying the appropriate matrices) and store results
      // in step0_matrix, step0_index.

      SCALAR scratch;

      for (int i = 0; i <= (1 << t); i++, current_index += (1 << t))
      {
         assert(next_target == target.size() ||
                target[next_target] >= current_index);

         // Skip over target endpoints which are exactly at the beginning
         // of this interval
         while ((next_target < target.size()) &&
                (target[next_target] == current_index))
         {
            // if there's an active matrix, don't forget to close it off
            if (active)
            {
               step0_index.push_back(current_index);
               active = 0;
            }
            next_target++;
         }

         // Test if any target endpoints are strictly within this interval.
         if ((next_target == target.size()) ||
             (target[next_target] >= current_index + (1 << t)))
         {
            // There are no target endpoints in this interval.
            if (active)
            {
               // Merge this matrix with the active one
               MATRIX& active_mat = step0_matrix.back();
               for (int y = 0; y < dim; y++)
                  for (int x = 0; x < dim; x++)
                  {
                     SCALAR& accum = temp_mat[y][x];
                     accum = 0;
                     for (int z = 0; z < dim; z++)
                     {
                        mul(scratch, active_mat[y][z],
                                     dyadic_output[z*dim + x][i]);
                        add(accum, accum, scratch);
                     }
                  }

               swap(temp_mat, active_mat);
            }
            else
            {
               // Make this matrix into a new active one
               step0_index.push_back(current_index);
               step0_matrix.resize(step0_matrix.size() + 1);
               MATRIX& X = step0_matrix.back();
               X.SetDims(dim, dim);
               for (int y = 0; y < dim; y++)
                  for (int x = 0; x < dim; x++)
                     X[y][x] = dyadic_output[y*dim + x][i];
               active = 1;
            }
         }
         else
         {
            // There are target endpoints in this interval.
            if (active)
            {
               // If there is still an active matrix, close it off.
               step0_index.push_back(current_index);
               active = 0;
            }

            // skip over any other endpoints in this interval
            while ((next_target < target.size()) &&
                   (target[next_target] < current_index + (1 << t)))
            {
               next_target++;
            }
         }
      }
   }

   // If there is still an active matrix, close it off.
   if (active)
      step0_index.push_back(current_index);

   assert(step0_index.size() == 2*step0_matrix.size());

   // =========================================================================
   // Step 1: Make a list of all subintervals that we are going to need in
   // the refining steps.

   int next_step0 = 0;        // index into step0_index
   vector<ZZ> step1_index;    // list of pairs of endpoints of needed intervals
   step1_index.reserve(2*target.size());

   // add sentinel endpoints to make the next loop simpler:
   step0_index.push_back(target.back() + 10);
   step0_index.push_back(target.back() + 20);

   for (next_target = 0; next_target < target.size(); next_target += 2)
   {
      // skip dyadic intervals that come before this target interval
      while (step0_index[next_step0+1] <= target[next_target])
         next_step0 += 2;

      if (step0_index[next_step0] < target[next_target+1])
      {
         // The next dyadic interval starts before the end of this target
         // interval.
         if (step0_index[next_step0] > target[next_target])
         {
            // The next dyadic interval starts strictly within this target
            // interval, so we need to create a refining subinterval for the
            // initial segment of this target interval.
            step1_index.push_back(target[next_target]);
            step1_index.push_back(step0_index[next_step0]);
         }

         // Skip over dyadic intervals to find the last one still contained
         // within this target interval.
         while (step0_index[next_step0+3] <= target[next_target+1])
            next_step0 += 2;

         if (step0_index[next_step0+1] < target[next_target+1])
         {
            // The next dyadic interval finishes strictly within this target
            // interval, so we need to create a refining subinterval for the
            // final segment of this target interval.
            step1_index.push_back(step0_index[next_step0+1]);
            step1_index.push_back(target[next_target+1]);
         }

         // Move on to next dyadic interval
         next_step0 += 2;
      }
      else
      {
         // The next dyadic interval starts beyond (or just at the end of)
         // this target interval, so we need to create a refining subinterval
         // for this *whole* target interval.
         step1_index.push_back(target[next_target]);
         step1_index.push_back(target[next_target+1]);
      }
   }

   // remove sentinels for my sanity
   step0_index.pop_back();
   step0_index.pop_back();

   // Step 1b: Compute matrix products over those refining subintervals.
   vector<MATRIX> step1_matrix;
   small_interval_products<SCALAR, POLY, POLYMODULUS, VECTOR, MATRIX>
                                   (step1_matrix, M0, M1, step1_index);

   assert(step1_index.size() == 2*step1_matrix.size());

   // =========================================================================
   // Step 2: Merge together the dyadic intervals and refining intervals into
   // a single list, in sorted order.

   vector<MATRIX> step2_matrix(step0_matrix.size() + step1_matrix.size());
   vector<ZZ> step2_index(step0_index.size() + step1_index.size());

   // add sentinels to make the next loop simpler
   step0_index.push_back(target.back() + 10);
   step0_index.push_back(target.back() + 20);
   step1_index.push_back(target.back() + 10);
   step1_index.push_back(target.back() + 20);

   next_step0 = 0;       // index into step0_matrix
   int next_step1 = 0;   // index into step1_matrix

   for (int next_step2 = 0; next_step2 < step2_matrix.size(); next_step2++)
   {
      if (step0_index[2*next_step0] < step1_index[2*next_step1])
      {
         // grab a matrix and pair of indices from step0
         swap(step2_matrix[next_step2], step0_matrix[next_step0]);
         step2_index[2*next_step2] = step0_index[2*next_step0];
         step2_index[2*next_step2+1] = step0_index[2*next_step0+1];
         next_step0++;
      }
      else
      {
         // grab a matrix and pair of indices from step1
         swap(step2_matrix[next_step2], step1_matrix[next_step1]);
         step2_index[2*next_step2] = step1_index[2*next_step1];
         step2_index[2*next_step2+1] = step1_index[2*next_step1+1];
         next_step1++;
      }
   }

   // remove sentinels for my sanity
   step0_index.pop_back();
   step0_index.pop_back();
   step1_index.pop_back();
   step1_index.pop_back();

   assert(step2_index.size() == 2*step2_matrix.size());

   // =========================================================================
   // Step 3: Walk through target intervals, and merge together appropriate
   // intervals from step2 to get those target intervals.

   int next_step2 = 0;    // index into step2_matrix

   // add sentinels to make the next loop simpler
   step2_index.push_back(target.back() + 1);
   step2_index.push_back(target.back() + 2);

   for (int next_target = 0; next_target < target.size(); next_target += 2)
   {
      // search for step2 interval matching the start of this target interval
      while (step2_index[2*next_step2] < target[next_target])
         next_step2++;

      assert(step2_index[2*next_step2] == target[next_target]);

      // merge together matrices for step2 intervals contained in this target
      // interval
      swap(step2_matrix[next_step2++], output[next_target/2]);
      while (step2_index[2*next_step2+1] <= target[next_target+1])
      {
         mul(temp_mat, output[next_target/2], step2_matrix[next_step2]);
         swap(temp_mat, output[next_target/2]);
         next_step2++;
      }
   }
}



/*
This is a wrapper for large_interval_products(), which takes input data in ZZ_p
format, and calls either the ZZ_p or zz_p version of large_interval_products(),
depending on how big the current modulus is.
*/
void large_interval_products_wrapper(vector<mat_ZZ_p>& output,
                                     const mat_ZZ_p& M0, const mat_ZZ_p& M1,
                                     const vector<ZZ>& target)
{
   const ZZ& modulus = ZZ_p::modulus();

   if (!modulus.SinglePrecision())
   {
      // Modulus is too big... just go straight to the ZZ_p version
      large_interval_products
         <ZZ_p, ZZ_pX, ZZ_pXModulus, vec_ZZ_p, mat_ZZ_p, FFTRep>
               (output, M0, M1, target);
   }
   else
   {
      // Modulus is small enough to work at single precision.

      // Save old single-precision modulus.
      zz_pContext context;
      context.save();

      // Set new single-precision modulus
      zz_p::init(to_long(modulus));

      // Convert input data to single-precision format
      int dim = M0.NumRows();
      mat_zz_p M0_sp, M1_sp;
      M0_sp.SetDims(dim, dim);
      M1_sp.SetDims(dim, dim);
      for (int x = 0; x < dim; x++)
         for (int y = 0; y < dim; y++)
         {
            M0_sp[x][y] = to_zz_p(rep(M0[x][y]));
            M1_sp[x][y] = to_zz_p(rep(M1[x][y]));
         }

      // Run large_interval_products at single precision
      vector<mat_zz_p> output_sp;
      large_interval_products
         <zz_p, zz_pX, zz_pXModulus, vec_zz_p, mat_zz_p, fftRep>
               (output_sp, M0_sp, M1_sp, target);

      // convert output back to ZZ_p format
      output.clear();
      output.resize(output_sp.size());
      for (int i = 0; i < output_sp.size(); i++)
      {
         output[i].SetDims(dim, dim);
         for (int x = 0; x < dim; x++)
            for (int y = 0; y < dim; y++)
               output[i][x][y] = to_ZZ_p(rep(output_sp[i][x][y]));
      }

      // restore old single-precision modulus
      context.restore();
   }
}



/* ============================================================================

   Main routine for computing Frobenius matrix

============================================================================ */


/*
Changes the modulus of a vector over ZZ_p, by lifting it from input_modulus
to ZZ, then reducing it back to output_modulus.
*/
void change_vec_modulus(vec_ZZ_p& output, const ZZ_pContext& output_modulus,
                       const vec_ZZ_p& input, const ZZ_pContext& input_modulus)
{
   int n = input.length();
   output.SetLength(n);

   for (int i = 0; i < n; i++)
   {
      input_modulus.restore();
      ZZ temp = rep(input[i]);
      output_modulus.restore();
      output[i] = to_ZZ_p(temp);
   }
}


/*
Given polynomials f and g defined over Z/p^N, this routine computes polynomials
a and b such that a*f + b*g = 1 mod p^N, with deg(a) < deg(g) and
deg(b) < deg(f).

PRECONDITIONS:
   f and g must be relatively prime mod p.
   The leading coefficients of f and g must be invertible mod p.
   The current NTL modulus should be p^N, and should be the one used to create
   the polynomials f and g.

Returns 1 on success, 0 if f and g are not relatively prime mod p.

NOTE:
   NTL can almost do this, but not quite; its xgcd routine is only guaranteed
   to work if the coefficient ring is a field. If you try it over Z/p^N, it
   usually works, but sometimes you get unlucky, and it crashes and burns.

*/
int padic_xgcd(ZZ_pX& a, ZZ_pX& b, const ZZ_pX& f, const ZZ_pX& g,
               const ZZ& p, int N)
{
   ZZ_pContext modulus;
   modulus.save();

   // =================================================
   // First do it mod p using NTL's xgcd routine

   // lift to Z
   ZZX f_lift, g_lift;
   conv(f_lift, f);
   conv(g_lift, g);

   // reduce down to Z/p
   ZZ_p::init(p);
   ZZ_pX f_red, g_red;
   conv(f_red, f_lift);
   conv(g_red, g_lift);

   // do xgcd mod p
   ZZ_pX a_red, b_red, d_red;
   XGCD(d_red, a_red, b_red, f_red, g_red);

   // lift results to Z
   ZZX a_lift, b_lift, d_lift;
   conv(a_lift, a_red);
   conv(b_lift, b_red);
   conv(d_lift, d_red);

   modulus.restore();

   if (deg(d_lift) != 0)
      return 0;

   // reduce back to Z/p^N
   conv(a, a_lift);
   conv(b, b_lift);
   ZZ_pX d;
   conv(d, d_lift);

   a /= d;
   b /= d;

   // =================================================
   // Now improve the approximation until we have enough precision

   for (int prec = 1; prec < N; prec *= 2)
   {
      ZZ_pX h = a*f + b*g - 1;
      ZZ_pX a_adj = -(h * a) % g;
      ZZ_pX b_adj = -(h * b) % f;
      a += a_adj;
      b += b_adj;
   }

   return 1;
}



/*
Given a matrix A over Z/p^N that is invertible mod p, this routine computes
the inverse B = A^(-1) over Z/p^N.

PRECONDITIONS:
   The current NTL modulus should be p^N, and should be the one used to create
   the matrix A.

NOTE:
   It's not clear to me (from the code or the documentation) whether NTL's
   matrix inversion is guaranteed to work over a non-field. So I implemented
   this one just in case.

*/
void padic_invert_matrix(mat_ZZ_p& B, const mat_ZZ_p& A, const ZZ& p, int N)
{
   ZZ_pContext modulus;
   modulus.save();

   int n = A.NumRows();

   // =================================================
   // First do it mod p using NTL matrix inverse

   // lift to Z
   mat_ZZ A_lift;
   A_lift.SetDims(n, n);
   for (int y = 0; y < n; y++)
      for (int x = 0; x < n; x++)
         A_lift[y][x] = rep(A[y][x]);

   // reduce down to Z/p
   ZZ_p::init(p);
   mat_ZZ_p A_red;
   A_red.SetDims(n, n);
   for (int y = 0; y < n; y++)
      for (int x = 0; x < n; x++)
         conv(A_red[y][x], A_lift[y][x]);

   // invert matrix mod p
   mat_ZZ_p B_red;
   inv(B_red, A_red);

   // lift result to Z
   mat_ZZ B_lift;
   B_lift.SetDims(n, n);
   for (int y = 0; y < n; y++)
      for (int x = 0; x < n; x++)
         B_lift[y][x] = rep(B_red[y][x]);

   // reduce back to Z/p^N
   modulus.restore();
   B.SetDims(n, n);
   for (int y = 0; y < n; y++)
      for (int x = 0; x < n; x++)
         conv(B[y][x], B_lift[y][x]);

   // =================================================
   // Now improve the approximation until we have enough precision

   mat_ZZ_p two;
   ident(two, n);
   two *= 2;

   for (int prec = 1; prec < N; prec *= 2)
      // if BA = I + error, then ((2I - BA)B)A = (I - err)(I + err) = I - err^2
      B = (two - B*A) * B;
}



/*
The main function exported from this module. See frobenius.h for information.
*/
int frobenius(mat_ZZ& output, const ZZ& p, int N, const ZZX& Q)
{
   // check input conditions
   if (N < 1 || p < 3)
      return 0;

   if ((deg(Q) < 3) || (deg(Q) % 2 == 0) || (!IsOne(LeadCoeff(Q))))
      return 0;

   int g = (deg(Q) - 1) / 2;

   if (p <= (2*N - 1) * deg(Q))
      return 0;

   // Back up the caller's modulus
   ZZ_pContext modulus_caller;
   modulus_caller.save();

   // Create moduli for working mod p^N and mod p^(N+1)

   ZZ_pContext modulus0, modulus1;

   ZZ_p::init(power(p, N));
   modulus0.save();
   ZZ_p::init(power(p, N+1));
   modulus1.save();

   // =========================================================================
   // Compute vertical reduction matrix M_V(t) and denominator D_V(t).

   // If N > 1 we need to do this to precision p^(N+1).
   // If N = 1 we can get away with precision N, since the last reduction
   // of length (p-1)/2 doesn't involve any divisions by p.
   ZZ_pContext modulusV = (N == 1) ? modulus0 : modulus1;
   modulusV.restore();

   // To compute the vertical reduction matrices, for each 0 <= i < 2g, we need
   // to find R_i(x) and S_i(x) satisfying x^i = R_i(x) Q(x) + S_i(x) Q'(x).
   vector<ZZ_pX> R, S;
   ZZ_pX Qred, Qred_diff, S0, R0;

   conv(Qred, Q);
   Qred_diff = diff(Qred);
   int result = padic_xgcd(R0, S0, Qred, Qred_diff, p, N+1);
   if (!result)
   {
      // failure if Q(x) and Q'(x) were not relatively prime mod p
      modulus_caller.restore();
      return 0;
   }

   S.push_back(S0);
   for (int i = 1; i < 2*g; i++)
      S.push_back(LeftShift(S.back(), 1) % Qred);

   R.push_back(R0);
   for (int i = 1; i < 2*g; i++)
      R.push_back(LeftShift(R.back(), 1) % Qred_diff);

   // Then the (j, i) entry of M_V(t) is given by the x^j coefficient of
   //    (2t-1) R_i(x) + 2 S_i'(x).
   // We store the constant term of M_V(t) in MV0, and the linear term in MV1.
   mat_ZZ_p MV0, MV1;
   MV0.SetDims(2*g, 2*g);
   MV1.SetDims(2*g, 2*g);
   for (int i = 0; i < 2*g; i++)
      for (int j = 0; j < 2*g; j++)
      {
         MV0[j][i] = -coeff(R[i], j) + 2*(j+1)*coeff(S[i], j+1);
         MV1[j][i] = 2*coeff(R[i], j);
      }

   // For convenience we also store the vertical reduction denominator
   // D_V(t) = 2t - 1 as a pair of 1x1 matrices.
   mat_ZZ_p DV0, DV1;
   DV0.SetDims(1, 1);
   DV1.SetDims(1, 1);
   DV0[0][0] = -1;
   DV1[0][0] = 2;

   // =========================================================================
   // Compute horizontal reduction matrices M_H^t(s) and denominators D_H^t(s),
   // where t = (p(2j+1) - 1)/2, for j in the range 0 <= j < N.
   // This is done to precision p^N.
   modulus0.restore();

   // We have D_H^t(s) = (2g+1)(2t-1) - 2s.
   // The matrix M_H^t(s) has D_H^t(s) on the sub-diagonal entries, and the
   // coefficients of the polynomial 2sQ(x) - (2t-1)xQ'(x) in the last column.

   vector<mat_ZZ_p> MH0(N), MH1(N), DH0(N), DH1(N);

   for (int j = 0; j < N; j++)
   {
      ZZ t = (p*(2*j+1) - 1)/2;

      MH0[j].SetDims(2*g+1, 2*g+1);
      MH1[j].SetDims(2*g+1, 2*g+1);
      DH0[j].SetDims(1, 1);
      DH1[j].SetDims(1, 1);

      DH0[j][0][0] = to_ZZ_p((2*g+1)*(2*t-1));
      DH1[j][0][0] = -2;
      // subdiagonal entries:
      for (int i = 0; i < 2*g; i++)
      {
         MH0[j][i+1][i] = DH0[j][0][0];
         MH1[j][i+1][i] = DH1[j][0][0];
      }
      // last column:
      for (int i = 0; i < 2*g+1; i++)
      {
         MH0[j][i][2*g] = -to_ZZ_p(i*(2*t-1))*coeff(Qred, i);
         MH1[j][i][2*g] = 2*coeff(Qred, i);
      }
   }

   // Compute inverse of the vandermonde matrix that will be used in the
   // horizontal reduction phase
   mat_ZZ_p vand_in, vand;
   vand_in.SetDims(N, N);
   for (int y = 1; y <= N; y++)
   {
      vand_in[y-1][0] = 1;
      for (int x = 1; x < N; x++)
         vand_in[y-1][x] = vand_in[y-1][x-1] * y;
   }
   padic_invert_matrix(vand, vand_in, p, N);

   // =========================================================================
   // Compute B_{j, r} coefficients.
   // These are done to precision p^(N+1).

   modulus1.restore();

   // The relevant binomial coefficients.
   vector<vector<ZZ_p> > binomial(2*N);
   binomial[0].push_back(to_ZZ_p(1));
   for (int j = 1; j < 2*N; j++)
   {
      binomial[j].push_back(to_ZZ_p(1));
      for (int i = 1; i < j; i++)
         binomial[j].push_back(binomial[j-1][i-1] + binomial[j-1][i]);
      binomial[j].push_back(to_ZZ_p(1));
   }

   // For 0 <= j < N, compute Btemp[j] =
   //           (-1)^j \sum_{k=j}^{N-1} 4^{-k} {2k \choose k} {k \choose j}
   ZZ_p fourth = to_ZZ_p(1) / 4;
   ZZ_p fourth_pow = to_ZZ_p(1);
   vector<ZZ_p> Btemp(N);
   for (int k = 0; k < N; k++, fourth_pow = fourth_pow * fourth)
   {
      // even j
      for (int j = 0; j <= k; j += 2)
         Btemp[j] += binomial[2*k][k] * binomial[k][j] * fourth_pow;
      // odd j
      for (int j = 1; j <= k; j += 2)
         Btemp[j] -= binomial[2*k][k] * binomial[k][j] * fourth_pow;
   }

   // Compute the coefficients   B_{j, r} = p C_{j, r} (-1)^j
   //             \sum_{k=j}^{N-1} 4^{-k} {2k \choose k} {k \choose j},
   // where C_{j, r} is the coefficient of x^r in Q(x)^j.
   vector<vector<ZZ_p> > B(N);
   ZZ_pX Qpow = to_ZZ_pX(p);
   for (int j = 0; j < N; j++, Qpow *= Qred)
      for (int r = 0; r <= (2*g+1)*j; r++)
         B[j].push_back(Btemp[j] * coeff(Qpow, r));

   // =========================================================================
   // Horizontal reductions.

   // reduced[j][i] will hold the reduction to x^0 of the i-th differential
   // along row j (each entry is a vector of length 2*g)
   vector<vector<vec_ZZ_p> > reduced(N, vector<vec_ZZ_p>(2*g));

   for (int j = 0; j < N; j++)
   {
      // Compute the transition matrices, mod p^N, between the indices
      // kp-2g-2 and kp, for each k. We compute at most N matrices (the others
      // will be deduced more efficiently using the vandermonde matrix below)
      modulus0.restore();
      vector<ZZ> s;
      int L = (2*g+1)*(j+1) - 1;
      int L_effective = min(L, N);
      for (int k = 0; k < L_effective; k++)
      {
         s.push_back(k*p);
         s.push_back((k+1)*p - 2*g - 2);
      }
      vector<mat_ZZ_p> MH, DH;
      MH.reserve(L);
      DH.reserve(L);
      large_interval_products_wrapper(MH, MH0[j], MH1[j], s);
      large_interval_products_wrapper(DH, DH0[j], DH1[j], s);

      // Use the vandermonde matrix to extend all the way up to L if necessary.
      if (L > N)
      {
         // First compute X[r] = F^{(r)}(0) p^r / r! for r = 0, ..., N-1,
         // where F is the matrix corresponding to M as described in the paper;
         // and do the same for Y corresponding to the denominator D.
         vector<mat_ZZ_p> X(N);
         vector<ZZ_p> Y(N);
         for (int r = 0; r < N; r++)
         {
            X[r].SetDims(2*g+1, 2*g+1);
            for (int h = 0; h < N; h++)
            {
               ZZ_p& v = vand[r][h];

               for (int y = 0; y < 2*g+1; y++)
                  for (int x = 0; x < 2*g+1; x++)
                     X[r][y][x] += v * MH[h][y][x];

               Y[r] += v * DH[h][0][0];
            }
         }

         // Now use those taylor coefficients to get the remaining MH's
         // and DH's.
         MH.resize(L);
         DH.resize(L);
         for (int k = N; k < L; k++)
         {
            MH[k].SetDims(2*g+1, 2*g+1);
            DH[k].SetDims(1, 1);

            // this is actually a power of k+1, since the indices into
            // MH and DH are one-off from what's written in the paper
            ZZ_p k_pow;
            k_pow = 1;

            for (int h = 0; h < N; h++, k_pow *= (k+1))
            {
               for (int y = 0; y < 2*g+1; y++)
                  for (int x = 0; x < 2*g+1; x++)
                     MH[k][y][x] += k_pow * X[h][y][x];

               DH[k][0][0] += k_pow * Y[h];
            }
         }
      }

      // Divide out each MH[k] by DH[k].
      for (int k = 0; k < MH.size(); k++)
      {
         ZZ_p Dinv = 1 / DH[k][0][0];

         for (int x = 0; x < 2*g+1; x++)
            for (int y = 0; y < 2*g+1; y++)
               MH[k][y][x] *= Dinv;
      }

      // Reduce all differentials on this row into the x^0 position.
      for (int i = 0; i < 2*g; i++)
      {
         modulus1.restore();
         vec_ZZ_p sum;
         sum.SetLength(2*g+1);

         // loop over blocks in this row
         for (int k = (2*g+1)*(j+1) - 1; k >= 1; k--)
         {
            // add in new term if necessary
            if (k >= i+1 && k <= i+1 + (2*g+1)*j)
               sum[0] += B[j][k - (i+1)];

            // We do the following reductions "by hand" for a little more
            // efficiency, since the matrices are all sparse.

            // tt is 2t - 1 in the notation of the paper
            ZZ tt = (2*j+1)*p - 2;

            // Reduce from kp - 1 down to kp - 2g - 1.
            for (int u = 1; u <= 2*g; u++)
            {
               // c is the last component
               ZZ_p c = sum[2*g];

               // shuffle over all components by one step
               for (int m = 2*g; m >= 1; m--)
                  sum[m] = sum[m-1];
               sum[0] = 0;

               // add in contribution from c (these follow from the explicit
               // formulae for the horizontal reductions)
               ZZ s = k*p - u;
               c /= to_ZZ_p((2*g+1)*tt - 2*s);
               for (int m = 0; m <= 2*g; m++)
                  sum[m] += c * to_ZZ_p(2*s - tt*m) * coeff(Qred, m);
            }

            // Reduce from kp - 2g - 1 to kp - 2g - 2.
            // This is a little different to the previous block owing to the
            // denominator being divisible by p.
            {
               // c is the last component
               ZZ_p c = sum[2*g];

               assert(rep(c) % p == 0);
               c = to_ZZ_p(rep(c) / p);

               // shuffle over all components by one step
               for (int m = 2*g; m >= 1; m--)
                  sum[m] = sum[m-1];
               sum[0] = 0;

               // add in contribution from c
               ZZ s = k*p - (2*g + 1);
               ZZ denom = (2*g+1)*tt - 2*s;
               assert(denom % p == 0);
               c /= to_ZZ_p(denom / p);
               for (int m = 0; m <= 2*g; m++)
                  sum[m] += c * to_ZZ_p(2*s - tt*m) * coeff(Qred, m);
            }

            // Use big fat matrices to reduce from kp - 2g - 2 to (k-1)p.
            {
               // first drop precision to p^N
               modulus0.restore();
               vec_ZZ_p temp;
               change_vec_modulus(temp, modulus0, sum, modulus1);

               // apply the matrix
               modulus0.restore();
               temp = MH[k-1] * temp;

               // lift precision back to p^(N+1)
               change_vec_modulus(sum, modulus1, temp, modulus0);
            }

            // Reduce from (k-1)p down to (k-1)p - 1.
            {
               // c is the last component
               ZZ_p c = sum[2*g];

               // shuffle over all components by one step
               for (int m = 2*g; m >= 1; m--)
                  sum[m] = sum[m-1];
               sum[0] = 0;

               // add in contribution from c
               ZZ s = (k-1)*p;
               c /= to_ZZ_p((2*g+1)*tt - 2*s);
               for (int m = 0; m <= 2*g; m++)
                  sum[m] += c * to_ZZ_p(2*s - tt*m) * coeff(Qred, m);
            }
         }

         // Store the reduced vector. Note that currently it has length 2g+1,
         // but the lowest term is zero, and the vertical reductions will
         // operate on length 2g vectors, so we kill the lowest term.
         modulus1.restore();
         reduced[j][i].SetLength(2*g);
         for (int f = 0; f < 2*g; f++)
            reduced[j][i][f] = sum[f+1];
      }
   }

   if (N == 1)
   {
      // If N == 1, then the vertical matrices are stored at precision 1
      // (instead of at N+1), so we reduce "reduced" to precision 1 here.
      modulus0.restore();
      vector<vector<vec_ZZ_p> > reduced_temp(N, vector<vec_ZZ_p>(2*g));
      for (int j = 0; j < N; j++)
         for (int i = 0; i < 2*g; i++)
            change_vec_modulus(reduced_temp[j][i], modulus0,
                               reduced[j][i], modulus1);

      reduced.swap(reduced_temp);
   }

   // =========================================================================
   // Vertical reductions.

   // Compute the vertical reduction maps between indices
   // 0, (p-1)/2, (3p-1)/2, (5p-1)/2, ..., ((2N-1)p-1)/2.
   vector<ZZ> s;
   s.push_back(to_ZZ(0));
   s.push_back((p-1)/2);
   for (int u = 1; u < N; u++)
   {
      s.push_back(s.back());
      s.push_back(s.back() + p);
   }
   modulusV.restore();
   vector<mat_ZZ_p> MV, DV;
   large_interval_products_wrapper(MV, MV0, MV1, s);
   large_interval_products_wrapper(DV, DV0, DV1, s);

   // Divide out each MV[j] by DV[j]. Note that for 1 <= j < N, both of them
   // should be divisible by p too, so we take that into account here.
   for (int j = 0; j < N; j++)
   {
      ZZ u = rep(DV[j][0][0]);
      if (j != 0)
      {
         assert(u % p == 0);
         u /= p;
      }
      InvMod(u, u, ZZ_p::modulus());

      for (int x = 0; x < 2*g; x++)
         for (int y = 0; y < 2*g; y++)
         {
            ZZ c = rep(MV[j][y][x]);
            if (j != 0)
            {
               assert(c % p == 0);
               c /= p;
            }
            MV[j][y][x] = to_ZZ_p(c * u);
         }
   }

   output.SetDims(2*g, 2*g);
   ZZ p_to_N = power(p, N);

   // For each i, use the above reductions matrices to vertically reduce the
   // terms corresponding to x^i dx/y.
   for (int i = 0; i < 2*g; i++)
   {
      vec_ZZ_p sum;
      sum.SetLength(2*g);

      for (int j = N-1; j >= 0; j--)
      {
         // pick up a new term,
         sum += reduced[j][i];
         // and reduce:
         sum = MV[j] * sum;
      }

      // store the result
      for (int f = 0; f < 2*g; f++)
         output[f][i] = rep(sum[f]) % p_to_N;
   }

   // Restore the caller's modulus
   modulus_caller.restore();

   return 1;
}


// ----------------------- end of file
