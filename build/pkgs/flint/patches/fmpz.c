/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

   fmpz.c: "flat" integer format

   Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "fmpz.h"
#include "flint.h"
#include "memory-manager.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "F_mpn_mul-tuning.h"
#include "long_extras.h"
#include "zn_poly/src/zn_poly.h"

#define SWAP_PTRS(x_dummy_p, y_dummy_p) \
do { \
   fmpz_t swap_temp_p = x_dummy_p; \
   x_dummy_p = y_dummy_p; \
   y_dummy_p = swap_temp_p; \
} while(0);

void fmpz_check_normalisation(const fmpz_t x)
{
   if ((x[0]) && (!x[ABS(x[0])]))
   {
      printf("Error: fmpz_t not normalised!\n");
      abort();
   }
}

void mpz_to_fmpz(fmpz_t res, const mpz_t x)
{
   if (mpz_sgn(x))
   {
      size_t countp;
      mpz_export(res + 1, &countp, -1, sizeof(mp_limb_t), 0, 0, x);
      res[0] = ((long) mpz_sgn(x) > 0L) ? (long) countp : (long) -countp;
   }
   else
      res[0] = 0L;
}

void fmpz_to_mpz(mpz_t res, const fmpz_t x)
{
   long size = x[0];

   if (size == 0)
      mpz_set_ui(res, 0);
   else
   {
      mpz_import(res, ABS(size), -1, sizeof(mp_limb_t), 0, 0, x + 1);

      if (size < 0)
         mpz_neg(res, res);
   }
}

void fmpz_print(fmpz_t in)
{
   mpz_t coeff;
   mpz_init(coeff);
   fmpz_to_mpz(coeff, in);
   gmp_printf("%Zd", coeff);
   mpz_clear(coeff);
}

/*
   Generate a random fmpz_t with n limbs with longs strings of 1's and 0's
*/

void fmpz_random_limbs2(fmpz_t x, unsigned long n)
{
   if (n == 0)
   {
      x[0] = 0L;
      return;
   }
   mpn_random2(x + 1, n);
   x[0] = n;
}

/*
    Adds two fmpz's together
*/

void fmpz_add(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2)
{
   fmpz_t coeffs1 = in1;
   fmpz_t coeffs2 = in2;

   long carry;
   unsigned long size1 = ABS(coeffs1[0]);
   unsigned long size2 = ABS(coeffs2[0]);

   if (size1 < size2)
   {
      SWAP_PTRS(coeffs1, coeffs2);
      size1 = ABS(coeffs1[0]);
      size2 = ABS(coeffs2[0]);
   }

   if (!size1)
   {
      if (!size2) coeffs_out[0] = 0L;
      else
      {
         if (coeffs_out != coeffs2) F_mpn_copy(coeffs_out, coeffs2, size2+1);
      }
   } else if (!size2)
   {
      if (coeffs_out != coeffs1) F_mpn_copy(coeffs_out, coeffs1, size1+1);
   } else if ((long) (coeffs1[0] ^ coeffs2[0]) >= 0L)
   {
      coeffs_out[0] = coeffs1[0];
      carry = mpn_add(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
      if (carry)
      {
         coeffs_out[size1+1] = carry;
         if ((long) coeffs_out[0] < 0L) coeffs_out[0]--;
         else coeffs_out[0]++;
      }
   } else
   {
      carry = 0;
      if (size1 != size2) carry = 1;
      else carry = mpn_cmp(coeffs1+1, coeffs2+1, size1);

      if (carry == 0) coeffs_out[0] = 0L;
      else if (carry > 0)
      {
         mpn_sub(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
         coeffs_out[0] = coeffs1[0];
         NORM(coeffs_out);
      }
      else
      {
         mpn_sub_n(coeffs_out+1, coeffs2+1, coeffs1+1, size1);
         coeffs_out[0] = -coeffs1[0];
         NORM(coeffs_out);
      }
   }
}

/*
     Add an unsigned long to an fmpz, inplace
*/

void fmpz_add_ui_inplace(fmpz_t output, const unsigned long x)
{
   unsigned long carry;

   if (x)
   {
      if (!output[0])
      {
         output[1] = x;
         output[0] = 1;
      } else if ((long) output[0] > 0)
      {
         carry = mpn_add_1(output + 1, output + 1, output[0], x);
         if (carry)
         {
            output[output[0]+1] = carry;
            output[0]++;
         }
      } else if ((long) output[0] < -1L)
      {
         mpn_sub_1(output + 1, output + 1, ABS(output[0]), x);
         NORM(output);
      } else
      {
         if (x <= output[1])
         {
            output[1] -= x;
            if (!output[1]) output[0] = 0;
         } else
         {
            output[1] = x - output[1];
            output[0] = 1;
         }
      }
   }
}

/*
     Add an unsigned long to an fmpz, inplace
*/

void fmpz_add_ui(fmpz_t output, const fmpz_t input, const unsigned long x)
{
   unsigned long carry;

   if (x)
   {
      if (!input[0])
      {
         output[1] = x;
         output[0] = 1;
      } else if ((long) input[0] > 0)
      {
         carry = mpn_add_1(output + 1, input + 1, input[0], x);
         output[0] = input[0];
         if (carry)
         {
            output[output[0]+1] = carry;
            output[0]++;
         }
      } else if ((long) input[0] < -1L)
      {
         mpn_sub_1(output + 1, input + 1, ABS(input[0]), x);
         output[0] = input[0];
         NORM(output);
      } else
      {
         if (x <= input[1])
         {
            output[1] = input[1] - x;
            if (!output[1]) output[0] = 0L;
            else output[0] = -1L;
         } else
         {
            output[1] = x - input[1];
            output[0] = 1L;
         }
      }
   } else
   {
      fmpz_set(output, input);
   }
}

/*
   Add an unsigned long to a coefficient.
   Assumes the output coefficient is non-negative.
*/

void __fmpz_add_ui_inplace(fmpz_t output, const unsigned long x)
{
   unsigned long carry;

   if (x)
   {
      if (!output[0])
      {
         output[1] = x;
         output[0] = 1;
      } else
      {
         carry = mpn_add_1(output + 1, output + 1, output[0], x);
         if (carry)
         {
            output[output[0]+1] = carry;
            output[0]++;
         }
      }
   }
}

void fmpz_sub(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2)
{
   fmpz_t coeffs1 = in1;
   fmpz_t coeffs2 = in2;

   long carry;
   unsigned long size1 = ABS(coeffs1[0]);
   unsigned long size2 = ABS(coeffs2[0]);
   int in_order = 1;

   if (size1 < size2)
   {
      SWAP_PTRS(coeffs1, coeffs2);
      size1 = ABS(coeffs1[0]);
      size2 = ABS(coeffs2[0]);
      in_order = 0;
   }

   if (!size1)
   {
      if (!size2) coeffs_out[0] = 0L;
      else
      {
         if (coeffs2 != coeffs_out) F_mpn_copy(coeffs_out, coeffs2, size2+1);
         if (in_order) coeffs_out[0] = -coeffs_out[0];
      }
   } else if (!size2)
   {
      if (coeffs1 != coeffs_out) F_mpn_copy(coeffs_out, coeffs1, size1+1);
      if (!in_order) coeffs_out[0] = -coeffs_out[0];
   } else if ((long) (coeffs1[0] ^ coeffs2[0]) < 0)
   {
      if (in_order) coeffs_out[0] = coeffs1[0];
      else coeffs_out[0] = -coeffs1[0];
      carry = mpn_add(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
      if (carry)
      {
         coeffs_out[size1+1] = carry;
         if ((long) coeffs_out[0] < 0) coeffs_out[0]--;
         else coeffs_out[0]++;
      }
   } else
   {
      carry = 0;
      if (size1 != size2) carry = 1;
      else carry = mpn_cmp(coeffs1+1, coeffs2+1, size1);

      if (carry == 0) coeffs_out[0] = 0L;
      else if (carry > 0)
      {
         mpn_sub(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
         if (in_order) coeffs_out[0] = coeffs1[0];
         else coeffs_out[0] = -coeffs1[0];
         NORM(coeffs_out);
      }
      else
      {
         mpn_sub_n(coeffs_out+1, coeffs2+1, coeffs1+1, size1);
         if (in_order) coeffs_out[0] = -coeffs1[0];
         else coeffs_out[0] = coeffs1[0];
         NORM(coeffs_out);
      }
   }
}

void fmpz_sub_ui_inplace(fmpz_t output, const unsigned long x)
{
   unsigned long carry;

   if (x)
   {
      if (!output[0])
      {
         output[1] = x;
         output[0] = -1L;
      } else if ((long) output[0] < 0)
      {
         carry = mpn_add_1(output + 1, output + 1, ABS(output[0]), x);
         if (carry)
         {
            output[ABS(output[0])+1] = carry;
            output[0]--;
         }
      } else if ((long) output[0] > 1L)
      {
         mpn_sub_1(output + 1, output + 1, output[0], x);
         NORM(output);
      } else
      {
         if (x <= output[1])
         {
            output[1] -= x;
            if (!output[1]) output[0] = 0;
         } else
         {
            output[1] = x - output[1];
            output[0] = -1L;
         }
      }
   }
}

void fmpz_sub_ui(fmpz_t output, const fmpz_t input, const unsigned long x)
{
   unsigned long carry;

   if (x)
   {
      if (!input[0])
      {
         output[1] = x;
         output[0] = -1L;
      } else if ((long) input[0] < 0)
      {
         carry = mpn_add_1(output + 1, input + 1, ABS(input[0]), x);
         output[0] = input[0];
         if (carry)
         {
            output[ABS(output[0])+1] = carry;
            output[0]--;
         }
      } else if ((long) input[0] > 1L)
      {
         mpn_sub_1(output + 1, input + 1, input[0], x);
         output[0] = input[0];
         NORM(output);
      } else
      {
         if (x <= input[1])
         {
            output[1] = input[1] - x;
            if (!output[1]) output[0] = 0;
            else output[0] = 1L;
         } else
         {
            output[1] = x - input[1];
            output[0] = -1L;
         }
      }
   } else
   {
      fmpz_set(output, input);
   }
}

/*
   Multiplies two fmpz's
   Assumes no overlap
*/

void fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b)
{
      long a0 = a[0];
      long b0 = b[0];
      unsigned long sizea = FLINT_ABS(a0);
      unsigned long sizeb = FLINT_ABS(b0);
      while ((!a[sizea]) && (sizea)) sizea--;
      while ((!b[sizeb]) && (sizeb)) sizeb--;

      mp_limb_t mslimb;
      fmpz_t temp;

      if ((sizea == 0) || (sizeb == 0))
      {
        res[0] = 0;
      } else if (sizea + sizeb < 100)
      {
         temp = (fmpz_t) flint_stack_alloc_small(sizea + sizeb + 1);
         if (sizea > sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
         else if (sizea == sizeb)
			{
				mpn_mul_n(temp+1, a+1, b+1, sizeb);
				mslimb = temp[2*sizeb];
			}
			else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
         temp[0] = sizea + sizeb - (mslimb == 0);
         F_mpn_copy(res, temp, temp[0]+1);
         if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
         flint_stack_release_small();
      } else if (sizea + sizeb < 2*FLINT_FFT_LIMBS_CROSSOVER)
      {
         temp = (fmpz_t) flint_stack_alloc(sizea + sizeb + 1);
         if (sizea > sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
			else if (sizea == sizeb)
			{
				mpn_mul_n(temp+1, a+1, b+1, sizeb);
				mslimb = temp[2*sizeb];
			}
			else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
         temp[0] = sizea + sizeb - (mslimb == 0);
         F_mpn_copy(res, temp, temp[0]+1);
         if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
         flint_stack_release();
      } else
      {
         if (sizea >= sizeb) mslimb = F_mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = F_mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea+sizeb - (mslimb == 0);
         if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
      }
}

/*
   Used internally by fmpz_poly.
   Multiplies two fmpz's but assumes res has enough space to contain the
   number of limbs in _a_ plus the number of limbs of _b_, whenever
   this sum is less than 2*FLINT_FFT_LIMBS_CROSSOVER
   Assumes no overlap
*/

void __fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b)
{
      long a0 = a[0];
      long b0 = b[0];
      unsigned long sizea = FLINT_ABS(a0);
      unsigned long sizeb = FLINT_ABS(b0);
      while ((!a[sizea]) && (sizea)) sizea--;
      while ((!b[sizeb]) && (sizeb)) sizeb--;

      mp_limb_t mslimb;
      fmpz_t temp;

      if ((sizea == 0) || (sizeb == 0))
      {
        res[0] = 0;
      }
      else if (sizea + sizeb < 100)
      {
         if (sizea >= sizeb) mslimb = mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea + sizeb - (mslimb == 0);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
      } else
      {
         if (sizea >= sizeb) mslimb = F_mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = F_mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea+sizeb - (mslimb == 0);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
      }
}


void fmpz_mul_ui(fmpz_t output, const fmpz_t input, const unsigned long x)
{
   if (x == 0)
   {
      output[0] = 0;
      return;
   }

   mp_limb_t mslimb;

   if (output[0] = input[0]) // This isn't a typo
   {
      mslimb = mpn_mul_1(output+1, input+1, FLINT_ABS(input[0]), x);
      if (mslimb)
      {
         output[FLINT_ABS(input[0])+1] = mslimb;
         if ((long) output[0] > 0) output[0]++;
         else output[0]--;
      }
   }
}

/*
   Sets res to res+a*b
   Assumes no overlap
*/

void fmpz_addmul(fmpz_t res, const fmpz_t a, const fmpz_t b)
{
      long a0 = a[0];
      long b0 = b[0];
      unsigned long sizea = FLINT_ABS(a0);
      unsigned long sizeb = FLINT_ABS(b0);
      while ((!a[sizea]) && (sizea)) sizea--;
      while ((!b[sizeb]) && (sizeb)) sizeb--;

      fmpz_t temp;
      mp_limb_t mslimb;

      if (sizea && sizeb)
      {
         if (sizea + sizeb < 100)
         {
            temp = (fmpz_t) flint_stack_alloc_small(sizea + sizeb + 1);
            if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
            else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
            temp[0] = sizea + sizeb - (mslimb == 0);
            if ((long) (a[0] ^ b[0]) < 0) temp[0] = -temp[0];
            fmpz_add(res, res, temp);
            flint_stack_release_small();
         } else
         {
            temp = (fmpz_t) flint_stack_alloc(sizea + sizeb + 1);
            if (sizea >= sizeb) mslimb = F_mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
            else mslimb = F_mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
            temp[0] = sizea + sizeb - (mslimb == 0);
            if ((long) (a[0] ^ b[0]) < 0) temp[0] = -temp[0];
            fmpz_add(res, res, temp);
            flint_stack_release();
         }
      }
}

/*
   Sets res to a / b
   Assumes no overlap
   Rounding occurs towards zero
*/

void fmpz_tdiv(fmpz_t res, const fmpz_t a, const fmpz_t b)
{
   long a0 = a[0];
   long b0 = b[0];
   unsigned long sizea = FLINT_ABS(a0);
   unsigned long sizeb = FLINT_ABS(b0);
   while ((!a[sizea]) && (sizea)) sizea--;
   while ((!b[sizeb]) && (sizeb)) sizeb--;


   mp_limb_t mslimb;
   fmpz_t temp;

   if (sizeb == 0)
   {
      printf("Error: division by zero!\n");
      abort();
   } else if (sizea < sizeb) // Todo: make this deal with sizea == sizeb but a < b
   {
      res[0] = 0;
   } else
   {
      temp = (fmpz_t) flint_stack_alloc(sizeb);
      mpn_tdiv_qr(res+1, temp, 0, a+1, sizea, b+1, sizeb);
      res[0] = sizea - sizeb + 1;
      if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
      flint_stack_release();
   }
   NORM(res);
}

/*
   Sets res to a / b
   Assumes no overlap
   Rounding occurs towards minus infinity
*/

void fmpz_fdiv(fmpz_t res, const fmpz_t a, const fmpz_t b)
{
   long a0 = a[0];
   long b0 = b[0];
   unsigned long sizea = FLINT_ABS(a0);
   unsigned long sizeb = FLINT_ABS(b0);
   while ((!a[sizea]) && (sizea)) sizea--;
   while ((!b[sizeb]) && (sizeb)) sizeb--;


   mp_limb_t mslimb;
   fmpz_t temp;

   if (sizeb == 0)
   {
      printf("Error: division by zero!\n");
      abort();
   } else if (sizea < sizeb) // Todo: make this deal with sizea == sizeb but a < b
   {
      if (((long) (a0 ^ b0) < 0L) && (a0))
      {
         res[0] = -1L;
         res[1] = 1;
      } else res[0] = 0;
      return;
   } else
   {
      temp = (fmpz_t) flint_stack_alloc(sizeb);
      mpn_tdiv_qr(res+1, temp, 0, a+1, sizea, b+1, sizeb);
      res[0] = sizea - sizeb + 1;
      if ((long) (a0 ^ b0) < 0L) res[0] = -res[0];
      NORM(res);
      if ((long) (a0 ^ b0) < 0L)
      {
         unsigned long i = 0;
         for (; i < sizeb; i++)
         {
            if (temp[i]) break;
         }
         if (i < sizeb)
         {
            fmpz_sub_ui_inplace(res, 1UL);
         }
      }
      flint_stack_release();
   }
}

/*
   Reduce a mod b, assuming b is positive
*/

void fmpz_mod(fmpz_t res, const fmpz_t a, const fmpz_t b)
{
   long a0 = a[0];
   long b0 = b[0];
   unsigned long sizea = FLINT_ABS(a0);
   unsigned long sizeb = b0;
   while ((!a[sizea]) && (sizea)) sizea--;
   while ((!b[sizeb]) && (sizeb)) sizeb--;


   mp_limb_t mslimb;
   fmpz_t temp, temp2;

   if (sizeb == 0)
   {
      printf("Error: division by zero!\n");
      abort();
   } else if (sizea < sizeb)
   {
      if ((long) a0 < 0L)
      {
         temp = (fmpz_t) flint_stack_alloc(sizeb + 2);
         fmpz_add(temp, a, b);
         fmpz_set(res, temp);
         flint_stack_release();
      } else fmpz_set(res, a);
      return;
   } else if ((sizea == sizeb) && (fmpz_cmpabs(a, b) < 0L))
   {
      if (fmpz_sgn(a) < 0) fmpz_add(res, a, b);
      else fmpz_set(res, a);
	  return;
   } else
   {
      if (fmpz_is_one(b))
	  {
	     fmpz_set_ui(res, 0L);
		 return;
	  }
	  temp = (fmpz_t) flint_stack_alloc(sizea - sizeb + 1);
      temp2 = (fmpz_t) flint_stack_alloc(sizeb + 2);
      mpn_tdiv_qr(temp, temp2+1, 0, a+1, sizea, b+1, sizeb);
      temp2[0] = sizeb;
      NORM(temp2);
      if (a0 < 0L)
      {
         unsigned long i = 0;
         for (; i < sizeb; i++)
         {
            if (temp2[i+1]) break;
         }
         if (i < sizeb)
         {
            fmpz_sub(temp2, b, temp2);
         }
         fmpz_set(res, temp2);
      } else fmpz_set(res, temp2);
      flint_stack_release();
      flint_stack_release();
   }
}

/*
   If b divides a then return 1 and the quotient, else return 0
*/

int fmpz_divides(fmpz_t q, const fmpz_t a, const fmpz_t b)
{
   long a0 = a[0];
   long b0 = b[0];
   unsigned long sizea = FLINT_ABS(a0);
   unsigned long sizeb = FLINT_ABS(b0);

   mp_limb_t mslimb;
   fmpz_t temp;

   if (sizeb == 0)
   {
      printf("Error: division by zero!\n");
      abort();
   } else if (sizea == 0)
	{
		q[0] = 0;
		return 1;
	} else if (sizea < sizeb)
   {
      return 0;
   } else if (sizea == sizeb)
	{
		int cmp = fmpz_cmpabs(a, b);
		if (cmp < 0) return 0;
      if (cmp == 0)
		{
			if ((long) (a0 ^ b0) < 0L) q[0] = -1L;
			else q[0] = 1L;
			q[1] = 1;
		   return 1;
		}
   }

   if (fmpz_is_one(b))
	{
	   fmpz_set(q, a);
		return 1;
	}

	if (fmpz_is_m1(b))
	{
		fmpz_neg(q, a);
		return 1;
	}

	temp = (fmpz_t) flint_stack_alloc(sizeb + 2);
   mpn_tdiv_qr(q+1, temp+1, 0, a+1, sizea, b+1, sizeb);

	temp[0] = sizeb;
   NORM(temp);

	if (temp[0] != 0)
	{
		flint_stack_release();
		return 0;
	}

	q[0] = sizea - sizeb + 1;
	NORM(q);

	if ((long) (a0 ^ b0) < 0L) q[0] = -q[0];

	flint_stack_release();
	return 1;
}

void fmpz_tdiv_ui(fmpz_t output, const fmpz_t input, const unsigned long x)
{
   output[0] = input[0];
   unsigned long size = FLINT_ABS(input[0]);

   mpn_divmod_1(output+1, input+1, size, x);

   NORM(output);
}

/*
   Returns input % x. Output will be reduced mod x.
*/

unsigned long fmpz_mod_ui(const fmpz_t input, const unsigned long x)
{
   unsigned long size = FLINT_ABS(input[0]);
   unsigned long mod;

   mod = mpn_mod_1(input+1, size, x);

   if (!mod) return mod;
   else if ((long) input[0] < 0L)
   {
      return x - mod;
   } else return mod;
}

// truncated multiplication for fmpz
void fmpz_mul_trunc(fmpz_t res, fmpz_t a, fmpz_t b, unsigned long trunc)
{
    unsigned long sizea = FLINT_MIN(fmpz_size(a), trunc);
    unsigned long sizeb = FLINT_MIN(fmpz_size(b), trunc);
    while ((!a[sizea]) && (sizea)) sizea--;
    while ((!b[sizeb]) && (sizeb)) sizeb--;

    if ((sizea == 0) || (sizeb == 0)) {
        res[0] = 0;
        return;
    }

    if (trunc >= sizea + sizeb) {
        mp_limb_t mslimb;
        if (sizea >= sizeb) mslimb = F_mpn_mul(res+1, a+1, sizea, b+1, sizeb);
        else mslimb = F_mpn_mul(res+1, b+1, sizeb, a+1, sizea);
        res[0] = sizea + sizeb - (mslimb == 0);
    } else {
        mp_limb_t mslimb;
        fmpz_t temp = flint_stack_alloc(sizea + sizeb + 1);
        if (sizea >= sizeb) mslimb = F_mpn_mul_trunc(temp+1, a+1, sizea, b+1, sizeb, trunc);
        else mslimb = F_mpn_mul_trunc(temp+1, b+1, sizeb, a+1, sizea, trunc);
        temp[0] = trunc;
        if (UNLIKELY(!mslimb))
            __fmpz_normalise(temp); // normalise if most significant limb == 0
        fmpz_set(res, temp);
        flint_stack_release();
    }
    if ((long) (a[0] ^ b[0]) < 0L) res[0] = -res[0];
}

// Compute a*b mod m
void fmpz_mulmod(fmpz_t res, fmpz_t a, fmpz_t b, fmpz_t m)
{
    fmpz_t ab = fmpz_init(fmpz_size(a) + fmpz_size(b));
    fmpz_mul(ab, a, b);
    fmpz_mod(res, ab, m);

    fmpz_clear(ab);
}

// Compute a/b mod m
void fmpz_divmod(fmpz_t res, fmpz_t a, fmpz_t b, fmpz_t m)
{
    fmpz_t b_inv = fmpz_init(fmpz_size(m));
    fmpz_invert(b_inv, b, m);

    fmpz_t ab = fmpz_init(fmpz_size(a) + fmpz_size(b_inv));
    fmpz_mul(ab, a, b_inv);
    fmpz_mod(res, ab, m);

    fmpz_clear(ab);
    fmpz_clear(b_inv);
}
/*
   Raise input to the power exp
   Very simplistic at this point. It just converts to an mpz_t and uses
   GMP's mpz_pow_ui function
*/

void fmpz_pow_ui(fmpz_t output, const fmpz_t input, const unsigned long exp)
{
   mpz_t power;
   mpz_init(power);
   fmpz_to_mpz(power, input);
   mpz_pow_ui(power, power, exp);
   mpz_to_fmpz(output, power);
   mpz_clear(power);
}

unsigned long __fmpz_power_of_two(const fmpz_t x)
{
   if (x[0] == 0) return -1L;
   return mpn_scan1(x + 1, 0);
}

void fmpz_mul_2exp(fmpz_t output, fmpz_t x, unsigned long exp)
{
   unsigned long limbs = (exp >> FLINT_LG_BITS_PER_LIMB);
   unsigned long bits = (exp & (FLINT_BITS - 1));
   mp_limb_t msl = 0L;

   if (x[0] == 0)
   {
      output[0] = 0L;
      return;
   }

   if (bits)
   {
      msl = mpn_lshift(output + limbs + 1, x + 1, FLINT_ABS(x[0]), bits);
      if (msl) output[limbs + FLINT_ABS(x[0]) + 1] = msl;
   } else F_mpn_copy(output + limbs + 1, x + 1, FLINT_ABS(x[0]));
   if (limbs) F_mpn_clear(output + 1, limbs);
   if ((long) x[0] >= 0L) output[0] = x[0] + limbs + (msl != 0L);
   else output[0] = x[0] - limbs - (msl != 0L);
}

void fmpz_div_2exp(fmpz_t output, fmpz_t x, unsigned long exp)
{
   unsigned long limbs = (exp >> FLINT_LG_BITS_PER_LIMB);
   unsigned long bits = (exp & (FLINT_BITS - 1));

   if ((x[0] == 0) || (limbs >= FLINT_ABS(x[0])))
   {
      output[0] = 0L;
      return;
   }

   if (bits)
   {
      fmpz_t temp = fmpz_init(FLINT_ABS(x[0]) - limbs);
      mpn_rshift(temp + 1, x + limbs + 1, FLINT_ABS(x[0]) - limbs, bits);
      if ((long) x[0] >= 0L) temp[0] = x[0] - limbs;
      else temp[0] = limbs + x[0];
      NORM(temp);
      fmpz_set(output, temp);
      fmpz_clear(temp);
   } else
   {
      F_mpn_copy(output + 1, x + limbs + 1, FLINT_ABS(x[0]) - limbs);
      if ((long) x[0] >= 0L) output[0] = x[0] - limbs;
      else output[0] = limbs + x[0];
   }
}

/*void fmpz_gcd(fmpz_t output, fmpz_t x1, fmpz_t x2)
{
   if (x1[0] == 0)
   {
      fmpz_set(output, x2);
      if ((long) output[0] < 0L) output[0] = -output[0];
      return;
   }

   if (x2[0] == 0)
   {
      fmpz_set(output, x1);
      if ((long) output[0] < 0L) output[0] = -output[0];
      return;
   }

   unsigned long twos1 = __fmpz_power_of_two(x1);
   unsigned long twos2 = __fmpz_power_of_two(x2);
   unsigned long n1, n2;

   fmpz_t a1 = fmpz_stack_init(FLINT_ABS(x1[0]) - (twos1 >> FLINT_LG_BITS_PER_LIMB));
   fmpz_t a2 = fmpz_stack_init(FLINT_ABS(x2[0]) - (twos2 >> FLINT_LG_BITS_PER_LIMB));

   fmpz_div_2exp(a1, x1, twos1);
   if ((long) a1[0] < 0L) a1[0] = -a1[0];
   fmpz_div_2exp(a2, x2, twos2);
   if ((long) a2[0] < 0L) a2[0] = -a2[0];

   if (fmpz_is_one(a1) || fmpz_is_one(a2))
   {
      fmpz_set_ui(output, 1UL);
   } else
   {
      n1 = FLINT_ABS(a1[0]);
      n2 = FLINT_ABS(a2[0]);
      if (fmpz_bits(a1) >= fmpz_bits(a2)) output[0] = mpn_gcd(output + 1, a1 + 1, n1, a2 + 1, n2);
      else output[0] = mpn_gcd(output + 1, a2 + 1, n2, a1 + 1, n1);
   }

   unsigned long min = FLINT_MIN(twos1, twos2);
   fmpz_mul_2exp(output, output, min);

   fmpz_stack_release(); // release a2
   fmpz_stack_release(); // release a1
}*/

void fmpz_gcd(fmpz_t output, fmpz_t x1, fmpz_t x2)
{
	ulong size1 = FLINT_ABS(x1[0]);
	ulong size2 = FLINT_ABS(x2[0]);
   mpz_t m1, m2, m0, m_out;
	m1->_mp_d = x1 + 1;
   m2->_mp_d = x2 + 1;
	m1->_mp_alloc = size1;
	m1->_mp_size = size1;
	m2->_mp_alloc = size2;
	m2->_mp_size = size2;
	ulong size_out = FLINT_MAX(size1, size2)+1;
	m0->_mp_d = flint_stack_alloc(size_out);
   m0->_mp_alloc = size_out;
   m0->_mp_size = 0;
	mpz_gcd(m0, m1, m2);
	size_out = m0->_mp_size;
	F_mpn_copy(output + 1, m0->_mp_d, size_out);
	output[0] = size_out;
	flint_stack_release();
}

/*
   Sets sqrt to the square root of n and sets rem to the remainder
   No aliasing of sqrt and n
*/

void fmpz_sqrtrem(fmpz_t sqrt, fmpz_t rem, fmpz_t n)
{
   long size = n[0];

   if (size < 0L)
   {
      printf("Cannot take the square root of a negative number!\n");
      abort();
   }

   if (!size)
   {
      fmpz_set_ui(sqrt, 0L);
      fmpz_set_ui(rem, 0L);
      return;
   }

   rem[0] = mpn_sqrtrem(sqrt+1, rem+1, n+1, size);

   sqrt[0] = (size+1)/2;
}

/*
   Invert x modulo m assuming m is positive
*/

void fmpz_invert(fmpz_t res, fmpz_t x, fmpz_t m)
{
   if (m[0] == 0)
   {
       printf("Error: division by zero!\n");
       abort();
   }

   fmpz_t s0, U, V, temp;
   unsigned long size = fmpz_size(m);
   // U may use fmpz_size(m) + 1 limbs after the sum, and gmp requires +1
   U = fmpz_init(size + 2);
   V = fmpz_init(size + 2);
   s0 = fmpz_init(size + 2);
   temp = fmpz_init(size + 2);

   // U := (x%m) + abs(m)
   // V := abs(m)
   fmpz_abs(V, m);
   fmpz_mod(U, x, V);
   fmpz_add(U, U, V);

   // Compute s0 such that 1 = s0 * x  %  m
   mpn_gcdext(temp+1, s0+1, (long *) s0, U+1, fmpz_size(U), V+1, fmpz_size(V));
   fmpz_mod(res, s0, m);

   fmpz_clear(temp);
   fmpz_clear(s0);
   fmpz_clear(V);
   fmpz_clear(U);
}

void fmpz_comb_init(fmpz_comb_t comb, ulong * primes, ulong num_primes)
{
   ulong i, j, k;

	comb->primes = primes;
   comb->num_primes = num_primes;

   ulong n = 0L;
   while (num_primes > (1L<<n)) n++;
   comb->n = n;
	ulong num;


   // create zn_poly modulus information
	comb->mod = (zn_mod_t *) flint_heap_alloc_bytes(sizeof(zn_mod_t)*num_primes);
   for (ulong i = 0; i < num_primes; i++)
      zn_mod_init(comb->mod[i], primes[i]);

	if (n == 0) return; // nothing to do

	// allocate space for comb
   comb->comb = (fmpz_t **) flint_heap_alloc(n);

   j = (1L<<(n - 1));
   ulong size = 2;
   mp_limb_t * ptr;
   for (i = 0; i < n; i++)
   {
      comb->comb[i] = (fmpz_t *) flint_heap_alloc(j);

      ptr = (mp_limb_t *) flint_heap_alloc((1L<<n) + j);
      for (k = 0; k < j; k++, ptr += (size + 1))
      {
         comb->comb[i][k] = ptr;
      }

      j/=2;
      size*=2;
   }

	// allocate space for res
	comb->res = (fmpz_t **) flint_heap_alloc(n);
   j = (1L<<(n - 1));
   size = 2;
   for (i = 0; i < n; i++)
   {
      comb->res[i] = (fmpz_t *) flint_heap_alloc(j);

      ptr = (mp_limb_t *) flint_heap_alloc((1L<<n) + j);
      for (k = 0; k < j; k++, ptr += (size + 1))
      {
         comb->res[i][k] = ptr;
      }

		j/=2;
      size*=2;
   }

	// compute products of pairs of primes and place in comb
	for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
   {
      fmpz_set_ui(comb->comb[0][j], primes[i]);
      fmpz_mul_ui(comb->comb[0][j], comb->comb[0][j], primes[i+1]);
   }
   if (i < num_primes) // in case number of primes is odd
   {
      fmpz_set_ui(comb->comb[0][j], primes[i]);
	   i+=2;
	   j++;
	}
   num = (1L<<n); // set the rest of the entries on that row of the comb to 1
	for (; i < num; i += 2, j++)
   {
      fmpz_set_ui(comb->comb[0][j], 1L);
	}

	// compute rest of comb by multiplying in pairs
	ulong log_comb = 1;
	num /= 2;
   while (num >= 2)
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         fmpz_mul(comb->comb[log_comb][j], comb->comb[log_comb-1][i], comb->comb[log_comb-1][i+1]);
      }
      log_comb++;
      num /= 2;
   }

	// compute inverses from pairs of primes
	fmpz_t temp = (fmpz_t) flint_stack_alloc(2);
   fmpz_t temp2 = (fmpz_t) flint_stack_alloc(2);
   for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
   {
       fmpz_set_ui(temp, primes[i]);
       fmpz_set_ui(temp2, primes[i+1]);
       fmpz_invert(comb->res[0][j], temp, temp2);
	}
   flint_stack_release(); //temp2
   flint_stack_release(); //temp

   ulong log_res = 1;
   num = (1L<<(n - 1));

   // compute remaining inverses, each level combining pairs from the level below
	while (log_res < n)
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         fmpz_invert(comb->res[log_res][j], comb->comb[log_res-1][i], comb->comb[log_res-1][i+1]);
      }
      log_res++;
      num /= 2;
   }
}

void fmpz_comb_clear(fmpz_comb_t comb)
{
   unsigned long n = comb->n;

   for (unsigned i = 0; i < n; i++)
   {
      flint_heap_free(comb->comb[i][0]);
      flint_heap_free(comb->comb[i]);
      flint_heap_free(comb->res[i][0]);
      flint_heap_free(comb->res[i]);
   }
	if (n)
	{
		flint_heap_free(comb->res);
      flint_heap_free(comb->comb);
	}

   flint_heap_free(comb->mod);
}

unsigned long fmpz_multi_mod_ui_basecase(unsigned long * out, fmpz_t in,
                               unsigned long * primes, unsigned long num_primes)
{
   for (unsigned long i = 0; i < num_primes; i++)
   {
      out[i] = fmpz_mod_ui(in, primes[i]);
   }
}

#define FLINT_LOG_MULTI_MOD_CUTOFF 2

void fmpz_multi_mod_ui(unsigned long * out, fmpz_t in, fmpz_comb_t comb)
{
   ulong i, j, k;
   ulong n = comb->n;
   long log_comb;
	ulong size;
	mp_limb_t * ptr;
	ulong num;
	unsigned long num_primes = comb->num_primes;

   if (num_primes == 1) // we are reducing modulo a single prime which is assumed to be big enough
   {
	  if ((long)in[0] > 0L) out[0] = in[1];
	  else if ((long)in[0] < 0L) out[0] = comb->primes[0] - in[1];
	  else out[0] = 0L;
	  return;
   }

   log_comb = n - 1;

	// allocate space for temp
	fmpz_t ** temp = (fmpz_t **) flint_heap_alloc(n);
   j = (1L<<(n - 1));
   size = 2;
   for (i = 0; i < n; i++)
   {
      temp[i] = (fmpz_t *) flint_heap_alloc(j);

      ptr = (mp_limb_t *) flint_heap_alloc((1L<<n)+j);
      for (k = 0; k < j; k++, ptr += (size + 1))
      {
         temp[i][k] = ptr;
      }

      j/=2;
      size*=2;
   }

   // find level in comb with entries bigger than the input integer
	log_comb = 0;
	if ((long) in[0] < 0L)
	  while ((fmpz_bits(in) >= fmpz_bits(comb->comb[log_comb][0]) - 1) && (log_comb < comb->n - 1)) log_comb++;
   else
		while (fmpz_cmpabs(in, comb->comb[log_comb][0]) >= 0) log_comb++;
   num = (1L<<(n - log_comb - 1));

	// set each entry of this level of temp to the input integer
	for (i = 0; i < num; i++)
   {
      fmpz_set(temp[log_comb][i], in);
   }
   log_comb--;
   num *= 2;

   // fill in other entries of temp by taking entries of temp at higher level mod pairs from comb
	while (log_comb > FLINT_LOG_MULTI_MOD_CUTOFF) // keep going until we reach the basecase
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         fmpz_mod(temp[log_comb][i], temp[log_comb + 1][j], comb->comb[log_comb][i]);
         fmpz_mod(temp[log_comb][i+1], temp[log_comb + 1][j], comb->comb[log_comb][i+1]);
      }
      num *= 2;
      log_comb--;
   }

   // do basecase
	num /= 2;
   log_comb++;
   ulong stride = (1L << (log_comb + 1));
   for (i = 0, j = 0; j < num_primes; i++, j += stride)
   {
	   fmpz_multi_mod_ui_basecase(out + j, temp[log_comb][i], comb->primes + j, FLINT_MIN(stride, num_primes - j));
   }

	// free temp
	for (i = 0; i < n; i++)
   {
      flint_heap_free(temp[i][0]);
      flint_heap_free(temp[i]);
	}
	flint_heap_free(temp);
}

void fmpz_multi_CRT_ui_unsigned(fmpz_t output, unsigned long * residues, fmpz_comb_t comb)
{
   ulong i, j, k;

   ulong n = comb->n;
   ulong num;
   ulong log_res;
	mp_limb_t * ptr;
	ulong size;
   ulong num_primes = comb->num_primes;

   if (num_primes == 1) // the output is less than a single prime, so just output the result
   {
      fmpz_set_ui(output, residues[0]);
	  return;
	}

   // allocate space for comb_temp
	fmpz_t ** comb_temp = (fmpz_t **) flint_heap_alloc(n);
   j = (1L<<(n - 1));
   size = 2;
   for (i = 0; i < n; i++)
   {
      comb_temp[i] = (fmpz_t *) flint_heap_alloc(j);

      ptr = (mp_limb_t *) flint_heap_alloc((1L<<n)+j);
      for (k = 0; k < j; k++, ptr += (size + 1))
      {
         comb_temp[i][k] = ptr;
      }

      j/=2;
      size*=2;
   }

	// first layer of reconstruction
	num = (1L<<n);
   fmpz_t temp = (fmpz_t) flint_stack_alloc(3);
   fmpz_t temp2 = (fmpz_t) flint_stack_alloc(3);
   for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
   {
      fmpz_set_ui(temp, residues[i]);
      fmpz_set_ui(temp2, fmpz_mod_ui(temp, comb->primes[i+1]));
      fmpz_sub_ui_inplace(temp2, residues[i + 1]);
      temp2[0] = -temp2[0];
      fmpz_mul(temp, temp2, comb->res[0][j]);
      fmpz_set_ui(temp2, fmpz_mod_ui(temp, comb->primes[i+1]));
      fmpz_mul_ui(temp, temp2, comb->primes[i]);
      fmpz_add_ui(comb_temp[0][j], temp, residues[i]);
   }
   if (i < num_primes) fmpz_set_ui(comb_temp[0][j], residues[i]);
   flint_stack_release(); //temp2
   flint_stack_release(); //temp

   // compute other layers of reconstruction
	temp = (fmpz_t) flint_stack_alloc(2*num + 1);
   temp2 = (fmpz_t) flint_stack_alloc(2*num + 1);
   num /= 2;
   log_res = 1;
   while (log_res < n)
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         if (fmpz_is_one(comb->comb[log_res-1][i+1]))
		   {
		      if (!fmpz_is_one(comb->comb[log_res-1][i])) fmpz_set(comb_temp[log_res][j], comb_temp[log_res-1][i]);
		   } else
		   {
			   fmpz_mod(temp2, comb_temp[log_res-1][i], comb->comb[log_res-1][i+1]);
            fmpz_sub(temp, temp2, comb_temp[log_res-1][i+1]);
            temp[0] = -temp[0];
            fmpz_mul(temp2, temp, comb->res[log_res][j]);
            fmpz_mod(temp, temp2, comb->comb[log_res-1][i+1]);
            fmpz_mul(temp2, temp, comb->comb[log_res-1][i]);
            fmpz_add(comb_temp[log_res][j], temp2, comb_temp[log_res-1][i]);
		   }
      }
      log_res++;
      num /= 2;
   }

   // write out the output
	fmpz_set(output, comb_temp[log_res - 1][0]);

	// free comb_temp
	for (i = 0; i < n; i++)
   {
      flint_heap_free(comb_temp[i][0]);
      flint_heap_free(comb_temp[i]);
	}
	flint_heap_free(comb_temp);

   flint_stack_release(); //temp2
   flint_stack_release(); //temp
}

void fmpz_multi_CRT_ui(fmpz_t output, unsigned long * residues, fmpz_comb_t comb)
{
   ulong i, j, k;

   ulong n = comb->n;
   ulong num;
   ulong log_res;
	mp_limb_t * ptr;
	ulong size;
   ulong num_primes = comb->num_primes;

   if (num_primes == 1) // the output is less than a single prime, so just output the result
   {
	   unsigned long p = comb->primes[0];
	   if ((p - residues[0]) < residues[0]) fmpz_set_si(output, (long) (residues[0] - p));
	   else fmpz_set_ui(output, residues[0]);
	   return;
	}

   // allocate space for comb_temp
	fmpz_t ** comb_temp = (fmpz_t **) flint_heap_alloc(n);
   j = (1L<<(n - 1));
   size = 2;
   for (i = 0; i < n; i++)
   {
      comb_temp[i] = (fmpz_t *) flint_heap_alloc(j);

      ptr = (mp_limb_t *) flint_heap_alloc((1L<<n)+j);
      for (k = 0; k < j; k++, ptr += (size + 1))
      {
         comb_temp[i][k] = ptr;
      }

      j/=2;
      size*=2;
   }

	// first layer of reconstruction
	num = (1L<<n);
   fmpz_t temp = (fmpz_t) flint_stack_alloc(3);
   fmpz_t temp2 = (fmpz_t) flint_stack_alloc(3);
   for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
   {
      fmpz_set_ui(temp, residues[i]);
      fmpz_set_ui(temp2, fmpz_mod_ui(temp, comb->primes[i+1]));
      fmpz_sub_ui_inplace(temp2, residues[i + 1]);
      temp2[0] = -temp2[0];
      fmpz_mul(temp, temp2, comb->res[0][j]);
      fmpz_set_ui(temp2, fmpz_mod_ui(temp, comb->primes[i+1]));
      fmpz_mul_ui(temp, temp2, comb->primes[i]);
      fmpz_add_ui(comb_temp[0][j], temp, residues[i]);
   }
   if (i < num_primes) fmpz_set_ui(comb_temp[0][j], residues[i]);
   flint_stack_release(); //temp2
   flint_stack_release(); //temp

   // compute other layers of reconstruction
	temp = (fmpz_t) flint_stack_alloc(2*num + 1);
   temp2 = (fmpz_t) flint_stack_alloc(2*num + 1);
   num /= 2;
   log_res = 1;
   while (log_res < n)
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         if (fmpz_is_one(comb->comb[log_res-1][i+1]))
		   {
		      if (!fmpz_is_one(comb->comb[log_res-1][i])) fmpz_set(comb_temp[log_res][j], comb_temp[log_res-1][i]);
		   } else
		   {
			   fmpz_mod(temp2, comb_temp[log_res-1][i], comb->comb[log_res-1][i+1]);
            fmpz_sub(temp, temp2, comb_temp[log_res-1][i+1]);
            temp[0] = -temp[0];
            fmpz_mul(temp2, temp, comb->res[log_res][j]);
            fmpz_mod(temp, temp2, comb->comb[log_res-1][i+1]);
            fmpz_mul(temp2, temp, comb->comb[log_res-1][i]);
            fmpz_add(comb_temp[log_res][j], temp2, comb_temp[log_res-1][i]);
		   }
      }
      log_res++;
      num /= 2;
   }

   // write out the output
	__fmpz_multi_CRT_sign(comb_temp[log_res - 1][0], comb_temp[log_res - 1][0], comb);
	fmpz_set(output, comb_temp[log_res - 1][0]);

	// free comb_temp
	for (i = 0; i < n; i++)
   {
      flint_heap_free(comb_temp[i][0]);
      flint_heap_free(comb_temp[i]);
	}
	flint_heap_free(comb_temp);

   flint_stack_release(); //temp2
   flint_stack_release(); //temp
}

void __fmpz_multi_CRT_sign(fmpz_t output, fmpz_t input, fmpz_comb_t comb)
{
   unsigned long n = comb->n;
   if (n == 0L)
   {
      if (input[0] == 0L)
	   {
	      fmpz_set_ui(output, 0L);
	      return;
	   }

	   unsigned long p = comb->primes[0];
	   if ((p - input[1]) < input[1]) fmpz_set_si(output, (long) (input[1] - p));
	   else fmpz_set_ui(output, input[1]);
	   return;
   }

   fmpz_t temp = fmpz_init(fmpz_size(comb->comb[n-1][0]) + 1);

   fmpz_sub(temp, input, comb->comb[comb->n - 1][0]);

   if (fmpz_cmpabs(temp, input) <= 0L) fmpz_set(output, temp);
   else fmpz_set(output, input);

   fmpz_clear(temp);
   return;
}

double fmpz_get_d(fmpz_t x)
{
	mpz_t x_m;
	x_m->_mp_d = x + 1;
	x_m->_mp_size = x[0];
	x_m->_mp_alloc = FLINT_ABS(x[0]);
	return mpz_get_d(x_m);
}

#include "fmpz_montgomery.c"

// *************** end of file
