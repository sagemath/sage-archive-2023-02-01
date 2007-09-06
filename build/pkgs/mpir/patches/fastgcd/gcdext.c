/* mpn_gcdext -- Extended Greatest Common Divisor.

Copyright 1996, 1998, 2000, 2001, 2002, 2003, 2004 Free Software Foundation,
Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

You should have received a copy of the GNU General Public License
along with the GNU MP Library; see the file COPYING.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#define WANT_TRACE 0

/* Default to binary gcdext_1, since it is best on most current machines.
   We should teach tuneup to choose the right gcdext_1.  */
#define GCDEXT_1_USE_BINARY 1

#if WANT_TRACE
# include <stdio.h>
# include <stdarg.h>
#endif

#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"

#ifndef NULL
# define NULL ((void *) 0)
#endif

#if WANT_TRACE
static void
trace (const char *format, ...)
{
  va_list args;
  va_start (args, format);
  gmp_vfprintf (stderr, format, args);
  va_end (args);
}
#endif

/* Comparison of _normalized_ numbers. */

#define MPN_EQUAL_P(ap, asize, bp, bsize)			\
((asize) == (bsize) && mpn_cmp ((ap), (bp), (asize)) == 0)

#define MPN_LEQ_P(ap, asize, bp, bsize)				\
((asize) < (bsize) || ((asize) == (bsize)			\
		       && mpn_cmp ((ap), (bp), (asize)) <= 0))

/* Returns g, u and v such that g = u A - v B. There are three
   different cases for the result:

     g = u A - v B, 0 < u < b, 0 < v < a
     g = A          u = 1, v = 0
     g = B          u = B, v = A - 1

   We always return with 0 < u <= b, 0 <= v < a.
*/
#if GCDEXT_1_USE_BINARY

static mp_limb_t
gcdext_1_odd (mp_limb_t *up, mp_limb_t *vp, mp_limb_t a, mp_limb_t b)
{
  mp_limb_t u0;
  mp_limb_t v0;
  mp_limb_t v1;
  mp_limb_t u1;

  mp_limb_t B = b;
  mp_limb_t A = a;

  /* Through out this function maintain

     a = u0 A - v0 B
     b = u1 A - v1 B

     where A and B are odd. */

  u0 = 1; v0 = 0;
  u1 = b; v1 = a-1;

  if (A == 1)
    {
      *up = u0; *vp = v0;
      return 1;
    }
  else if (B == 1)
    {
      *up = u1; *vp = v1;
      return 1;
    }

  while (a != b)
    {
      mp_limb_t mask;

      ASSERT (a % 2 == 1);
      ASSERT (b % 2 == 1);

      ASSERT (0 < u0); ASSERT (u0 <= B);
      ASSERT (0 < u1); ASSERT (u1 <= B);

      ASSERT (0 <= v0); ASSERT (v0 < A);
      ASSERT (0 <= v1); ASSERT (v1 < A);

      if (a > b)
	{
	  MP_LIMB_T_SWAP (a, b);
	  MP_LIMB_T_SWAP (u0, u1);
	  MP_LIMB_T_SWAP (v0, v1);
	}

      ASSERT (a < b);

      /* Makes b even */
      b -= a;

      mask = - (mp_limb_t) (u1 < u0);
      u1 += B & mask;
      v1 += A & mask;
      u1 -= u0;
      v1 -= v0;

      ASSERT (b % 2 == 0);

      do
	{
	  /* As b = u1 A + v1 B is even, while A and B are odd,
	     either both or none of u1, v1 is even */

	  ASSERT (u1 % 2 == v1 % 2);

	  mask = -(u1 & 1);
	  u1 = u1 / 2 + ((B / 2) & mask) - mask;
	  v1 = v1 / 2 + ((A / 2) & mask) - mask;

	  b /= 2;
	}
      while (b % 2 == 0);
    }

  /* Now g = a = b */
  ASSERT (a == b);
  ASSERT (u1 <= B);
  ASSERT (v1 < A);

  ASSERT (A % a == 0);
  ASSERT (B % a == 0);
  ASSERT (u0 % (B/a) == u1 % (B/a));
  ASSERT (v0 % (A/a) == v1 % (A/a));

  *up = u0; *vp = v0;

  return a;
}

static mp_limb_t
gcdext_1 (mp_limb_t *up, mp_limb_t *vp, mp_limb_t a, mp_limb_t b)
{
  unsigned shift = 0;
  mp_limb_t g;
  mp_limb_t u;
  mp_limb_t v;

  /* We use unsigned values in the range 0, ... B - 1. As the values
     are uniquely determined only modulo B, we can add B at will, to
     get numbers in range or flip the least significant bit. */
  /* Deal with powers of two */
  while ((a | b) % 2 == 0)
    {
      a /= 2; b /= 2; shift++;
    }

  if (b % 2 == 0)
    {
      unsigned k = 0;

      do {
	b /= 2; k++;
      } while (b % 2 == 0);

      g = gcdext_1_odd (&u, &v, a, b);

      while (k--)
	{
	  /* We have g = u a + v b, and need to construct
	     g = u'a + v'(2b).

	     If v is even, we can just set u' = u, v' = v/2
	     If v is odd, we can set v' = (v + a)/2, u' = u + b
	  */

	  if (v % 2 == 0)
	    v /= 2;
	  else
	    {
	      u = u + b;
	      v = v/2 + a/2 + 1;
	    }
	  b *= 2;
	}
    }
  else if (a % 2 == 0)
    {
      unsigned k = 0;

      do {
	a /= 2; k++;
      } while (a % 2 == 0);

      g = gcdext_1_odd (&u, &v, a, b);

      while (k--)
	{
	  /* We have g = u a + v b, and need to construct
	     g = u'(2a) + v'b.

	     If u is even, we can just set u' = u/2, v' = v.
	     If u is odd, we can set u' = (u + b)/2
	  */

	  if (u % 2 == 0)
	    u /= 2;
	  else
	    {
	      u = u/2 + b/2 + 1;
	      v = v + a;
	    }
	  a *= 2;
	}
    }
  else
    /* Ok, both are odd */
    g = gcdext_1_odd (&u, &v, a, b);

  *up = u;
  *vp = v;

  return g << shift;
}

#else /* ! GCDEXT_1_USE_BINARY */
static mp_limb_t
gcdext_1_u (mp_limb_t *up, mp_limb_t a, mp_limb_t b)
{
  /* Maintain

     a =   u0 A mod B
     b = - u1 A mod B
  */
  mp_limb_t u0 = 1;
  mp_limb_t u1 = 0;
  mp_limb_t B = b;

  ASSERT (a >= b);
  ASSERT (b > 0);

  for (;;)
    {
      mp_limb_t q;

      q = a / b;
      a -= q * b;

      if (a == 0)
	{
	  *up = B - u1;
	  return b;
	}
      u0 += q * u1;

      q = b / a;
      b -= q * a;

      if (b == 0)
	{
	  *up = u0;
	  return a;
	}
      u1 += q * u0;
    }
}

static mp_limb_t
gcdext_1 (mp_limb_t *up, mp_limb_t *vp, mp_limb_t a, mp_limb_t b)
{
  /* Maintain

     a =   u0 A - v0 B
     b = - u1 A + v1 B = (B - u1) A - (A - v1) B
  */
  mp_limb_t u0 = 1;
  mp_limb_t v0 = 0;
  mp_limb_t u1 = 0;
  mp_limb_t v1 = 1;

  mp_limb_t A = a;
  mp_limb_t B = b;

  ASSERT (a >= b);
  ASSERT (b > 0);

  for (;;)
    {
      mp_limb_t q;

      q = a / b;
      a -= q * b;

      if (a == 0)
	{
	  *up = B - u1;
	  *vp = A - v1;
	  return b;
	}
      u0 += q * u1;
      v0 += q * v1;

      q = b / a;
      b -= q * a;

      if (b == 0)
	{
	  *up = u0;
	  *vp = v0;
	  return a;
	}
      u1 += q * u0;
      v1 += q * v0;
    }
}
#endif /* ! GCDEXT_1_USE_BINARY */

/* FIXME: Duplicated in gcd.c */
static mp_size_t
hgcd_tdiv (mp_ptr qp,
	   mp_ptr rp, mp_size_t *rsizep,
	   mp_srcptr ap, mp_size_t asize,
	   mp_srcptr bp, mp_size_t bsize)
{
  mp_size_t qsize;
  mp_size_t rsize;

  mpn_tdiv_qr (qp, rp, 0, ap, asize, bp, bsize);

  rsize = bsize;
  MPN_NORMALIZE (rp, rsize);
  *rsizep = rsize;

  qsize = asize - bsize + 1;
  qsize -= (qp[qsize - 1] == 0);

  if (qsize == 1 && qp[0] == 1)
    return 0;

  return qsize;
}

/* FIXME: Duplicated in hgcd.c */
static mp_limb_t
mpn_addmul2_n_1 (mp_ptr rp, mp_size_t n,
		 mp_ptr ap, mp_limb_t u,
		 mp_ptr bp, mp_limb_t v)
{
  mp_limb_t h;
  mp_limb_t cy;

  h = mpn_mul_1 (rp, ap, n, u);
  cy = mpn_addmul_1 (rp, bp, n, v);
  h += cy;
#if GMP_NAIL_BITS == 0
  rp[n] = h;
  return (h < cy);
#else /* GMP_NAIL_BITS > 0 */
  rp[n] = h & GMP_NUMB_MASK;
  return h >> GMP_NUMB_BITS;
#endif /* GMP_NAIL_BITS > 0 */
}


/* Computes u2 = u0 + q u1

   Returns new size.

   FIXME: Notation in the function not quite consistent
   FIXME: Severe code duplication with hgcd_update_uv */

static mp_size_t
hgcd_update_u (struct hgcd_row *r, mp_size_t usize,
	       mp_srcptr qp, mp_size_t qsize,
	       /* Limbs allocated for the new u, for sanity
		   checking */
	       mp_size_t alloc)
{
  mp_srcptr u0p = r[0].uvp[0];
  mp_srcptr u1p = r[1].uvp[0];
  mp_ptr u2p = r[2].uvp[0];

  ASSERT (usize < alloc);

  /* u1 = 0 is an exceptional case. Except for this, u1 should be
     normalized. */

  ASSERT ((usize == 1 && u1p[0] == 0) || u1p[usize - 1] != 0);

  /* Compute u2  = u0 + q u1 */

  if (usize == 1 && u1p[0] == 0)
    {
      /* u1 == 0 is a special case, then q might be large, but it
	 doesn't matter. Can happen only when u0 = v1 = 1, u1 = v0 =
	 0, and hence usize == 1. */
      MPN_COPY (u2p, u0p, usize);
    }
  else if (qsize == 0)
    /* Represents a unit quotient */
    {
      mp_limb_t cy = mpn_add_n (u2p, u0p, u1p, usize);
      u2p[usize] = cy;
      usize += (cy != 0);
    }
  else if (qsize == 1)
    {
      mp_limb_t cy;

      cy = mpn_mul_1 (u2p, u1p, usize, qp[0]);
      cy += mpn_add_n (u2p, u2p, u0p, usize);

      u2p[usize] = cy;
      usize += (cy != 0);
    }
  else
    {
      if (qsize <= usize)
	mpn_mul (u2p, u1p, usize, qp, qsize);
      else
	mpn_mul (u2p, qp, qsize, u1p, usize);

      ASSERT_NOCARRY (mpn_add (u2p,
			       u2p, usize + qsize,
			       u0p, usize));

      usize += qsize;
      usize -= (u2p[usize - 1] == 0);
    }
  ASSERT (mpn_cmp (r[1].uvp[0], r[2].uvp[0], usize) <= 0);
  ASSERT (r[2].uvp[0][usize - 1] != 0);

  return usize;
}


/* Computes Y = R * X. No overlap allowed. */
static mp_size_t
hgcd2_mul_vector (struct hgcd_row *Y,
		  mp_size_t alloc,
		  const struct hgcd2_row *R,
		  const struct hgcd_row *X, mp_size_t n)
{
  unsigned i;
  int grow = 0;
  mp_limb_t h = 0;

  ASSERT (n < alloc);

  for (i = 0; i < 2; i++)
    {
      /* Set Y[i] = R[i, 0] X[0] + R[i,1] X[1]
		  = u X[0] + v X[0] */
      mp_limb_t cy;

      cy = mpn_addmul2_n_1 (Y[i].uvp[0], n,
			    X[0].uvp[0], R[i].u,
			    X[1].uvp[0], R[i].v);

      if (cy)
	{
	  ASSERT (n + 2 <= alloc);
	  Y[i].uvp[0][n+1] = cy;
	  grow = 1;
	}
      else
	h |= Y[i].uvp[0][n];
    }
  if (grow)
    return n + 2;
  else
    /* Don't add redundant zeroes */
    return n + (h != 0);
}

/* Sets (a, b, c)  <--  (b, c, a) */
#define HGCD_SWAP3_LEFT(row)				\
do {							\
  struct hgcd_row __hgcd_swap4_left_tmp = row[0];	\
  row[0] = row[1];					\
  row[1] = row[2];					\
  row[2] = __hgcd_swap4_left_tmp;			\
} while (0)

/* Sets (a, b, c, d)  <--  (c, d, a, b) */
#define HGCD_SWAP4_2(row)				\
do {							\
  struct hgcd_row __hgcd_swap4_2_tmp = row[0];	\
  row[0] = row[2];					\
  row[2] = __hgcd_swap4_2_tmp;				\
  __hgcd_swap4_2_tmp = row[1];				\
  row[1] = row[3];					\
  row[3] = __hgcd_swap4_2_tmp;				\
} while (0)

static mp_size_t
gcdext_lehmer_itch (mp_size_t asize, mp_size_t bsize)
{
  mp_size_t ralloc = asize + 1;
  mp_size_t ualloc = bsize + 1;

  return 4 * ralloc + 4 * ualloc + asize;
}

static mp_size_t
gcdext_lehmer (mp_ptr gp, mp_ptr up, mp_size_t *usize,
	       mp_srcptr ap, mp_size_t asize,
	       mp_srcptr bp, mp_size_t bsize,
	       mp_ptr tp, mp_size_t talloc)
{
  struct hgcd_row r[4];
  /* Size and sign of u fields. The largest u should be normalized to
     this size, and except for the case u1 = 0, that is the latest
     u. */
  int rsize;
  int rsign;

  mp_ptr qp;
  mp_size_t qsize;
  mp_size_t ralloc = asize + 1;
  mp_size_t ualloc = bsize + 1;

  struct hgcd2 hgcd;
  int res;

  ASSERT (asize >= bsize);
  ASSERT (asize > 1);
  ASSERT (bsize > 0);

  ASSERT (MPN_LEQ_P (bp, bsize, ap, asize));

  ASSERT (4 * ralloc + 4*ualloc + asize <= talloc);

  r[0].rp = tp; tp += ralloc; talloc -= ralloc;
  r[1].rp = tp; tp += ralloc; talloc -= ralloc;
  r[2].rp = tp; tp += ralloc; talloc -= ralloc;
  r[3].rp = tp; tp += ralloc; talloc -= ralloc;

  /* Must zero out the u fields. We don't use the v fields. */
  MPN_ZERO (tp, 4 * ualloc);

  r[0].uvp[0] = tp; tp += ualloc; talloc -= ualloc;
  r[1].uvp[0] = tp; tp += ualloc; talloc -= ualloc;
  r[2].uvp[0] = tp; tp += ualloc; talloc -= ualloc;
  r[3].uvp[0] = tp; tp += ualloc; talloc -= ualloc;

  qp = tp; tp += asize; talloc -= asize;

  res = mpn_hgcd2_lehmer_step (&hgcd,
			       ap, asize,
			       bp, bsize,
			       NULL);

  if (res == 0 || (res == 2 && hgcd.row[0].v == 0))
    {
      qsize = hgcd_tdiv (qp, r[1].rp, &r[1].rsize,
			 ap, asize,
			 bp, bsize);
      MPN_COPY (r[0].rp, bp, bsize);
      r[0].rsize = bsize;

      r[0].uvp[0][0] = 0;
      r[1].uvp[0][0] = 1;
      rsign = -1;
    }
  else
    {
      const struct hgcd2_row *s = hgcd.row + (res - 2);
      rsign = hgcd.sign;
      if (res == 3)
	rsign = ~rsign;

      /* s[0] and s[1] correct. */
      r[0].rsize
	= mpn_hgcd2_fix (r[0].rp, ralloc,
			 rsign,
			 s[0].u, ap, asize,
			 s[0].v, bp, bsize);

      r[1].rsize
	= mpn_hgcd2_fix (r[1].rp, ralloc,
			 ~rsign,
			 s[1].u, ap, asize,
			 s[1].v, bp, bsize);

      r[0].uvp[0][0] = s[0].u;
      r[1].uvp[0][0] = s[1].u;
    }
  rsize = 1;

  while (r[0].rsize >= 2 && r[1].rsize > 0)
    {
      res = mpn_hgcd2_lehmer_step (&hgcd,
				   r[0].rp, r[0].rsize,
				   r[1].rp, r[1].rsize,
				   NULL);

      if (res == 0 || (res == 2 && hgcd.row[0].v == 0))
	{
	  qsize = hgcd_tdiv (qp, r[2].rp, &r[2].rsize,
			     r[0].rp, r[0].rsize,
			     r[1].rp, r[1].rsize);
	  rsize = hgcd_update_u (r, rsize, qp, qsize, ualloc);
	  HGCD_SWAP3_LEFT (r);
	  rsign = ~rsign;
	}
      else
	{
	  const struct hgcd2_row *s = hgcd.row + (res - 2);
	  int sign = hgcd.sign;
	  if (res == 3)
	    sign = ~sign;

	  /* s[0] and s[1] correct. */
	  r[2].rsize
	    = mpn_hgcd2_fix (r[2].rp, ralloc,
			     sign,
			     s[0].u, r[0].rp, r[0].rsize,
			     s[0].v, r[1].rp, r[1].rsize);

	  r[3].rsize
	    = mpn_hgcd2_fix (r[3].rp, ralloc,
			     ~sign,
			     s[1].u, r[0].rp, r[0].rsize,
			     s[1].v, r[1].rp, r[1].rsize);

	  rsize = hgcd2_mul_vector (r + 2, ralloc, s, r, rsize);
	  rsign ^= sign;
	  HGCD_SWAP4_2 (r);
	}
    }

  if (r[1].rsize == 0)
    {
      MPN_NORMALIZE (r[0].uvp[0], rsize);
      MPN_COPY (gp, r[0].rp, r[0].rsize);
      MPN_COPY (up, r[0].uvp[0], rsize);

      *usize = (rsign >= 0) ? rsize : -rsize;
      return r[0].rsize;
    }
  else
    {
      mp_limb_t cy;
      mp_limb_t u;
      mp_limb_t v;

      gp[0] = gcdext_1 (&u, &v, r[0].rp[0], r[1].rp[0]);
      cy = mpn_addmul2_n_1 (up, rsize,
			    r[0].uvp[0], u,
			    r[1].uvp[0], v);
      rsize++;
      if (cy)
	up[rsize++] = cy;
      else
	MPN_NORMALIZE (up, rsize);

      *usize = (rsign >= 0) ? rsize : -rsize;
      return 1;
    }
}

/* Computes Y = R * X. No overlap allowed.

   Temporary space is needed for two numbers smaller than the
   resulting matrix elements, i.e. bounded by 2*L <= N.

   FIXME: Severe code duplication with hgcd.c: hgcd_mul. */

static mp_size_t
hgcd_mul_vector (struct hgcd_row *Y, mp_size_t alloc,
		 const struct hgcd_row *R, mp_size_t rsize,
		 const struct hgcd_row *X, mp_size_t xsize,
		 mp_ptr tp, mp_size_t talloc)
{
  unsigned i;

  mp_size_t ysize;
  mp_limb_t h;
  int grow;

  MPN_NORMALIZE (R[1].uvp[1], rsize);
  /* u1 = 0 is an exceptional case. Except for this, u1 should be
     normalized. */
  ASSERT ((xsize == 1 && X[1].uvp[0][0] == 0)
	  || X[1].uvp[0][xsize - 1] != 0);

  if (xsize == 1 && X[1].uvp[0][0] == 0)
    {
      /* Special case. Set Y[i, 0] = R[i, 0] */
      ASSERT (X[0].uvp[0][0] == 1);

      if (rsize > 1)
	MPN_NORMALIZE (R[1].uvp[0], rsize);
      MPN_COPY (Y[0].uvp[0], R[0].uvp[0], rsize);
      MPN_COPY (Y[1].uvp[0], R[1].uvp[0], rsize);

      return rsize;
    }

  ysize = rsize + xsize;
  ASSERT (ysize <= talloc);

  h = 0; grow = 0;

  if (rsize >= xsize)
    {
      for (i = 0; i < 2; i++)
	{
	  /* Set Y[i, 0] = R[i, 0] X[0, 0] + R[i,1] X[1, 0] */
	  mp_limb_t cy;

	  mpn_mul (Y[i].uvp[0], R[i].uvp[0], rsize, X[0].uvp[0], xsize);
	  mpn_mul (tp, R[i].uvp[1], rsize, X[1].uvp[0], xsize);

	  cy = mpn_add_n (Y[i].uvp[0], Y[i].uvp[0], tp, ysize);

	  if (cy)
	    {
	      ASSERT (ysize + 1 < alloc);
	      Y[i].uvp[0][ysize] = cy;
	      grow = 1;
	    }
	  else
	    h |= Y[i].uvp[0][ysize - 1];
	}
    }
  else
    {
      for (i = 0; i < 2; i++)
	{
	  /* Set Y[i, 0] = R[i, 0] X[0, 0] + R[i,1] X[1, 0] */
	  mp_limb_t cy;

	  mpn_mul (Y[i].uvp[0], X[0].uvp[0], xsize, R[i].uvp[0], rsize);
	  mpn_mul (tp, X[1].uvp[0], xsize, R[i].uvp[1], rsize);

	  cy = mpn_add_n (Y[i].uvp[0], Y[i].uvp[0], tp, ysize);

	  if (cy)
	    {
	      ASSERT (ysize + 1 < alloc);
	      Y[i].uvp[0][ysize] = cy;
	      grow = 1;
	    }
	  else
	    h |= Y[i].uvp[0][ysize - 1];
	}
    }

  if (grow)
    ysize++;
  else
    ysize -= (h == 0);

  ASSERT ((ysize == 1 && Y[1].uvp[0][0] == 0) || Y[1].uvp[0][ysize - 1] != 0);

  return ysize;
}

#define COMPUTE_V_ITCH(asize, bsize, usize) \
  ((usize) + (asize) + 1 + (bsize))

/* Computes |v| = |(c - u a)| / b, where u may be positive or negative,
   and v is of the opposite sign. Requires that b, c, |u| <= a. */
static mp_size_t
compute_v (mp_ptr vp, mp_size_t valloc,
	   mp_srcptr ap, mp_size_t asize,
	   mp_srcptr bp, mp_size_t bsize,
	   mp_srcptr cp, mp_size_t csize,
	   mp_srcptr up, mp_size_t usize,
	   mp_ptr tp, mp_size_t talloc)
{
  mp_size_t size;
  mp_size_t vsize;
  mp_ptr rp;

  ASSERT (asize);
  ASSERT (bsize);
  ASSERT (csize);
  ASSERT (asize >= bsize);

#if 0
  trace ("compute_v: a = %Nd\n"
	"           b = %Nd\n"
	"           c = %Nd\n"
	"           u = %Nd\n",
	ap, asize, bp, bsize, cp, csize, up, usize);
#endif

  ASSERT (usize);

  size = ABS (usize);

  ASSERT (size <= asize);
  ASSERT (asize + size <= talloc);

  mpn_mul (tp, ap, asize, up, size);
  size += asize;

  ASSERT (csize <= size);

  if (usize > 0)
    {
      /* |v| = -v = (u a - c) / b */

      ASSERT_NOCARRY (mpn_sub (tp, tp, size, cp, csize));
      MPN_NORMALIZE (tp, size);
      if (size == 0)
	return 0;
    }
  else
    { /* usize < 0 */
      /* |v| = v = (c - u a) / b = (c + |u| a) / b */
      mp_limb_t cy = mpn_add (tp, tp, size, cp, csize);
      if (cy)
	{
	  ASSERT (size < talloc);
	  tp[size++] = cy;
	}
    }

  /* Now divide t / b. There must be no remainder */

  ASSERT (size >= bsize);
  ASSERT (size + bsize <= talloc);
  rp = tp + size;

  vsize = size + 1 - bsize;
  ASSERT (vsize <= valloc);

  mpn_tdiv_qr (vp, rp, 0, tp, size, bp, bsize);
  MPN_NORMALIZE (vp, vsize);

  /* Remainder must be zero */
#if WANT_ASSERT
  {
    mp_size_t i;
    for (i = 0; i < bsize; i++)
      {
	ASSERT (rp[i] == 0);
      }
  }
#endif
  return vsize;
}

static mp_size_t
gcdext_schoenhage_itch (mp_size_t asize, mp_size_t bsize)
{
  mp_size_t itch;

  mp_size_t ralloc = asize + 1;
  mp_size_t ualloc = bsize + 1;
  /* Input size for hgcd calls */
  mp_size_t halloc = (asize + 1) / 2;

  /* Storage for the rows and quotient */
  mp_size_t rstorage = 4 * ralloc + 4 * ualloc + asize;

  /* Storage for hgcd calls */
  mp_size_t tstorage = mpn_hgcd_init_itch (halloc)
    + qstack_itch (halloc)
    + mpn_hgcd_itch (halloc);

  /* Storage needed for final gcdext_lehmer */
  mp_size_t lstorage
    = gcdext_lehmer_itch (GCDEXT_SCHOENHAGE_THRESHOLD,
			  GCDEXT_SCHOENHAGE_THRESHOLD);

  /* Storage needed after final nhgcd_gcdext_lehmer */
  mp_size_t fstorage
    = COMPUTE_V_ITCH (GCDEXT_SCHOENHAGE_THRESHOLD,
		      GCDEXT_SCHOENHAGE_THRESHOLD,
		      ualloc);

  /* We need rstorage + MAX (tstorage, lstorage, fstorage) */

  itch = tstorage;
  if (lstorage > tstorage)
    itch = lstorage;
  if (fstorage > itch)
    itch = fstorage;

  return rstorage + itch;
}

#if WANT_ASSERT
static void
sanity_check_row (mp_srcptr ap, mp_size_t asize,
		  mp_srcptr bp, mp_size_t bsize,
		  int sign, mp_size_t usize,
		  const struct hgcd_row *r)
{
  /* Check that x = u * a + v * b, for some v, i.e. that
     x - u*a is divisible by b. */
  mp_srcptr up = r->uvp[0];
  mp_srcptr xp = r->rp;
  mp_size_t xsize = r->rsize;
  mp_ptr tp;
  mp_size_t tsize;
  mp_ptr qp;
  mp_size_t qsize;
  mp_ptr rp;
  mp_size_t i;
  TMP_DECL;
  TMP_MARK;

  ASSERT (asize > 0 && ap[asize - 1] != 0);
  ASSERT (bsize > 0 && bp[bsize - 1] != 0);
  ASSERT (xsize == 0 || xp[xsize - 1] != 0);
  ASSERT (MPN_LEQ_P (xp, xsize, ap, asize));
  ASSERT (MPN_LEQ_P (up, usize, bp, bsize));

  MPN_NORMALIZE (up, usize);
  if (usize == 0)
    {
      ASSERT (MPN_EQUAL_P (xp, xsize, bp, bsize));
      return;
    }

  tp = TMP_ALLOC_LIMBS (usize + asize + 1);
  qp = TMP_ALLOC_LIMBS (usize + asize + 2 - bsize);
  rp = TMP_ALLOC_LIMBS (bsize);

  mpn_mul (tp, ap, asize, up, usize);
  tsize = asize + usize;
  tsize -= (tp[tsize - 1] == 0);

  if (sign >= 0)
    {
      ASSERT_NOCARRY (mpn_sub (tp, tp, tsize, xp, xsize));
      MPN_NORMALIZE (tp, tsize);
    }
  else
    {
      mp_limb_t cy = mpn_add (tp, tp, tsize, xp, xsize);
      tp[tsize] = cy;
      tsize += (cy != 0);
    }

  if (tsize > 0)
    {
      mpn_tdiv_qr (qp, rp, 0, tp, tsize, bp, bsize);
      for (i = 0; i < bsize; i++)
	ASSERT (rp[i] == 0);
      qsize = tsize - bsize;
      qsize += (qp[qsize] != 0);
      ASSERT (MPN_LEQ_P (qp, qsize, ap, asize));
    }
  TMP_FREE;
}
# define ASSERT_ROW(ap, asize, bp, bsize, sign, usize, r) \
sanity_check_row (ap, asize, bp, bsize, sign, usize, r)

#else /* !WANT_ASSERT */
# define ASSERT_ROW(ap, asize, bp, bsize, sign, usize, r)
#endif /* !WANT_ASSERT */

static mp_size_t
gcdext_schoenhage (mp_ptr gp, mp_ptr up, mp_size_t *usizep,
		   mp_srcptr ap, mp_size_t asize,
		   mp_srcptr bp, mp_size_t bsize,
		   mp_ptr tp, mp_size_t talloc)
{
  mp_size_t scratch;
  struct hgcd hgcd;
  struct qstack quotients;
  struct hgcd_row r[4];

  /* Size and sign of u fields. The largest u should be normalized to
     this size, and except for the case u1 = 0, that is the latest
     u. */
  int rsize;
  int rsign;

  mp_ptr qp;
  mp_size_t qsize;
  mp_size_t ralloc = asize + 1;
  mp_size_t ualloc = bsize + 1;

  ASSERT (asize >= bsize);
  ASSERT (bsize > 0);

  ASSERT (MPN_LEQ_P (bp, bsize, ap, asize));

  ASSERT (4 * ralloc + 4*ualloc + asize <= talloc);

  r[0].rp = tp; tp += ralloc; talloc -= ralloc;
  r[1].rp = tp; tp += ralloc; talloc -= ralloc;
  r[2].rp = tp; tp += ralloc; talloc -= ralloc;
  r[3].rp = tp; tp += ralloc; talloc -= ralloc;

  /* Must zero out the u fields */
  MPN_ZERO (tp, 4 * ualloc);

  r[0].uvp[0] = tp; tp += ualloc; talloc -= ualloc;
  r[1].uvp[0] = tp; tp += ualloc; talloc -= ualloc;
  r[2].uvp[0] = tp; tp += ualloc; talloc -= ualloc;
  r[3].uvp[0] = tp; tp += ualloc; talloc -= ualloc;

  qp = tp; tp += asize; talloc -= asize;

  ASSERT (asize >= bsize);
  ASSERT (bsize > 0);
  MPN_COPY (r[0].rp, ap, asize); r[0].rsize = asize;
  MPN_COPY (r[1].rp, bp, bsize); r[1].rsize = bsize;

  r[0].uvp[0][0] = 1;
  r[1].uvp[0][0] = 0;

  /* We don't use the v fields. */
  rsize = 1;
  rsign = 0;

  scratch = mpn_hgcd_init_itch ((asize + 1) / 2);
  ASSERT (scratch <= talloc);
  mpn_hgcd_init (&hgcd, (asize + 1) / 2, tp);
  tp += scratch; talloc -= scratch;

  {
    mp_size_t nlimbs = qstack_itch ((asize + 1) / 2);

    ASSERT (nlimbs <= talloc);
    qstack_init (&quotients, (asize + 1) / 2, tp, nlimbs);

    tp += nlimbs;
    talloc -= nlimbs;
    scratch += nlimbs;
  }

  while (ABOVE_THRESHOLD (r[0].rsize, GCDEXT_SCHOENHAGE_THRESHOLD)
	 && r[1].rsize > 0)
    {
      mp_size_t k = r[0].rsize / 2;
      int res;

      ASSERT_ROW (ap, asize, bp, bsize, rsign, rsize, r);
      ASSERT_ROW (ap, asize, bp, bsize, ~rsign, rsize, r + 1);

      if (r[1].rsize <= k)
	goto euclid;

      qstack_reset (&quotients, r[0].rsize - k);

      res = mpn_hgcd (&hgcd,
		      r[0].rp + k, r[0].rsize - k,
		      r[1].rp + k, r[1].rsize - k,
		      &quotients,
		      tp, talloc);

      if (res == 0 || res == 1)
	{
	euclid:
	  qsize = hgcd_tdiv (qp, r[2].rp, &r[2].rsize,
			     r[0].rp, r[0].rsize,
			     r[1].rp, r[1].rsize);
	  rsize = hgcd_update_u (r, rsize, qp, qsize, ualloc);
	  ASSERT (rsize < ualloc);

	  ASSERT_ROW (ap, asize, bp, bsize, rsign, rsize, r + 2);

	  HGCD_SWAP3_LEFT (r);
	  rsign = ~rsign;
	}
      else
	{
	  const struct hgcd_row *s = hgcd.row + (res - 2);
	  int sign = hgcd.sign;
	  if (res == 3)
	    sign = ~sign;

	  /* s[0] and s[1] are correct */
	  r[2].rsize
	    = mpn_hgcd_fix (k, r[2].rp, ralloc,
			    sign, hgcd.size, s,
			    r[0].rp, r[1].rp,
			    tp, talloc);

	  r[3].rsize
	    = mpn_hgcd_fix (k, r[3].rp, ralloc,
			    ~sign, hgcd.size, s+1,
			    r[0].rp, r[1].rp,
			    tp, talloc);

	  rsize = hgcd_mul_vector (r + 2, ualloc, s, hgcd.size,
				   r, rsize, tp, talloc);
	  ASSERT (rsize < ualloc);

	  rsign ^= sign;
	  ASSERT_ROW (ap, asize, bp, bsize, rsign, rsize, r + 2);
	  ASSERT_ROW (ap, asize, bp, bsize, ~rsign, rsize, r + 3);

	  HGCD_SWAP4_2 (r);
	}
    }
  if (r[1].rsize == 0)
    {
      MPN_COPY (gp, r[0].rp, r[0].rsize);
      MPN_NORMALIZE (r[0].uvp[0], rsize);
      MPN_COPY (up, r[0].uvp[0], rsize);

      *usizep = (rsign >= 0) ? rsize : - rsize;
      return r[0].rsize;
    }
  else if (r[0].rsize == 1)
    {
      mp_limb_t u;
      mp_limb_t v;
      mp_limb_t cy;

      gp[0] = gcdext_1 (&u, &v, r[0].rp[0], r[1].rp[0]);

      /* g = u r0 + v r1 = (u u0 + v u1) a + (...) b */
      cy = mpn_addmul2_n_1 (up, rsize,
			    r[0].uvp[0], u,
			    r[1].uvp[0], v);

      rsize++;
      if (cy)
	up[rsize++] = cy;
      else
	MPN_NORMALIZE (up, rsize);

      *usizep = (rsign >= 0) ? rsize : -rsize;
      return 1;

    }
  else
    {
      /* We have r0 = u0 a + v0 b,
		 r1 = u1 a + v1 b

	 Compute g = u r0 + v r1 = (u u0 + v u1) a + (...) b
	 In the expression (u u0 + v u1), we have

	 u  <= r1,
	 u0 <= b/r0 (except if r0 = a, which should never be the case here)
	 v  <= r0
	 u1 <= b/r0
      */

      mp_size_t gsize;
      mp_size_t usize;
      mp_size_t vsize;

      /* u1 should be non-zero, and normalized */
      ASSERT (rsize);
      ASSERT (r[1].uvp[0][rsize - 1] != 0);
#if WANT_TRACE
      trace ("gcdext: \n"
	     "r0 = %Nd\n"
	     "r1 = %Nd\n"
	     "u0 = %Nd\n"
	     "u1 = %Nd\n",
	     r[0].rp, r[0].rsize, r[1].rp, r[1].rsize,
	     r[0].uvp[0], rsize, r[1].uvp[0], rsize);
#endif
      /* We don't need the space for hgcd and the quotient stack any more */
      tp -= scratch; talloc += scratch;

      /* Stores u in r[2] and v in r[3] */
      gsize = gcdext_lehmer (gp, r[2].uvp[0], &usize,
			     r[0].rp, r[0].rsize,
			     r[1].rp, r[1].rsize,
			     tp, talloc);

      if (usize == 0)
	{
	  /* u == 0  ==>  v = g / b == 1  ==> g = u1 a + (...) b */

	  MPN_NORMALIZE (r[1].uvp[0], rsize);
	  MPN_COPY (up, r[1].uvp[0], rsize);
	  *usizep = (rsign >= 0) ? - rsize : rsize;

	  return gsize;
	}

      /* Compute v = (g - s r0) / r1, storing it in r[3] */
      vsize = compute_v (r[3].uvp[0], ualloc,
			 r[0].rp, r[0].rsize, r[1].rp, r[1].rsize,
			 gp, gsize,
			 r[2].uvp[0], usize,
			 tp, talloc);

      if (usize < 0)
	{
	  usize = - usize;
	  rsign = ~rsign;
	}

      /* It's possible that u0 = 0, u1 = 1 */
      if (rsize == 1 && r[0].uvp[0][0] == 0)
	{
	  /* u0 == 0 ==> u u0 + v u1 = v */
	  MPN_COPY (up, r[3].uvp[0], vsize);
	  *usizep = (rsign >= 0) ? vsize : - vsize;

	  return gsize;
	}

      /* Ok, now u0, u1, u are non-zero. We may still have v == 0 */
      ASSERT (usize + rsize <= ualloc);
      ASSERT (vsize + rsize <= ualloc);

      /* Compute u u0 */
      if (usize <= rsize)
	/* Should be the common case */
	mpn_mul (up,
		 r[0].uvp[0], rsize,
		 r[2].uvp[0], usize);
      else
	mpn_mul (up,
		 r[2].uvp[0], usize,
		 r[0].uvp[0], rsize);

      usize += rsize;

      /* There may be more than one zero limb, if #u0 < #u1 */
      MPN_NORMALIZE (up, usize);
      ASSERT (usize < ualloc);

      if (vsize)
	{
	  mp_limb_t cy;

	  /* Overwrites old r[2].uvp[0] value */
	  if (vsize <= rsize)
	    /* Should be the common case */
	    cy = mpn_mul (r[2].uvp[0],
			  r[1].uvp[0], rsize,
			  r[3].uvp[0], vsize);
	  else
	    cy = mpn_mul (r[2].uvp[0],
			  r[3].uvp[0], vsize,
			  r[1].uvp[0], rsize);

	  vsize += rsize - (cy == 0);
	  ASSERT (vsize < ualloc);

	  if (vsize <= usize)
	    cy = mpn_add (up, up, usize, r[2].uvp[0], vsize);
	  else
	    {
	      cy = mpn_add (up, r[2].uvp[0], vsize, up, usize);
	      usize = vsize;
	    }
	  up[usize] = cy;
	  usize += (cy != 0);

	  ASSERT (usize < ualloc);
	}
      *usizep = (rsign >= 0) ? usize : -usize;

      return gsize;
    }
}

mp_size_t
mpn_gcdext (mp_ptr gp, mp_ptr up, mp_size_t *usizep,
	    mp_ptr ap, mp_size_t asize, mp_ptr bp, mp_size_t bsize)
{
  ASSERT (asize >= bsize);
  ASSERT (bsize > 0);

  if (asize == 1)
    {
#if GCDEXT_1_USE_BINARY
      mp_limb_t v;
      *gp = gcdext_1 (up, &v, ap[0], bp[0]);
#else
      *gp = gcdext_1_u (up, ap[0], bp[0]);
#endif
      *usizep = (up[0] != 0);
      ASSERT(gp[0] != 0);
      return 1;
    }
  else if (BELOW_THRESHOLD (asize, GCDEXT_SCHOENHAGE_THRESHOLD))
    {
      mp_size_t gsize;
      mp_ptr tp;
      mp_size_t talloc = gcdext_lehmer_itch (asize, bsize);
      TMP_DECL;
      TMP_MARK;

      tp = TMP_ALLOC_LIMBS (talloc);
      gsize = gcdext_lehmer (gp, up, usizep, ap, asize, bp, bsize,
			     tp, talloc);
      TMP_FREE;
      return gsize;
    }
  else
    {
      mp_size_t gsize;
      mp_ptr tp;
      mp_size_t talloc = gcdext_schoenhage_itch (asize, bsize);
      TMP_DECL;
      TMP_MARK;

      tp = TMP_ALLOC_LIMBS (talloc);
      gsize = gcdext_schoenhage (gp, up, usizep, ap, asize, bp, bsize,
				 tp, talloc);
      TMP_FREE;
      return gsize;
    }
}
