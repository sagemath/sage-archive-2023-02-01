/* hgcd2.c

   THE FUNCTIONS IN THIS FILE ARE INTERNAL WITH MUTABLE INTERFACES.  IT IS ONLY
   SAFE TO REACH THEM THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST
   GUARANTEED THAT THEY'LL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.

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

#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"

#if GMP_NAIL_BITS == 0

/* Copied from mpn/generic/gcdext.c, and modified slightly to return
   the remainder. */
/* Two-limb division optimized for small quotients.  */
static inline mp_limb_t
div2 (mp_ptr rp,
      mp_limb_t nh, mp_limb_t nl,
      mp_limb_t dh, mp_limb_t dl)
{
  mp_limb_t q = 0;

  if ((mp_limb_signed_t) nh < 0)
    {
      int cnt;
      for (cnt = 1; (mp_limb_signed_t) dh >= 0; cnt++)
	{
	  dh = (dh << 1) | (dl >> (GMP_LIMB_BITS - 1));
	  dl = dl << 1;
	}

      while (cnt)
	{
	  q <<= 1;
	  if (nh > dh || (nh == dh && nl >= dl))
	    {
	      sub_ddmmss (nh, nl, nh, nl, dh, dl);
	      q |= 1;
	    }
	  dl = (dh << (GMP_LIMB_BITS - 1)) | (dl >> 1);
	  dh = dh >> 1;
	  cnt--;
	}
    }
  else
    {
      int cnt;
      for (cnt = 0; nh > dh || (nh == dh && nl >= dl); cnt++)
	{
	  dh = (dh << 1) | (dl >> (GMP_LIMB_BITS - 1));
	  dl = dl << 1;
	}

      while (cnt)
	{
	  dl = (dh << (GMP_LIMB_BITS - 1)) | (dl >> 1);
	  dh = dh >> 1;
	  q <<= 1;
	  if (nh > dh || (nh == dh && nl >= dl))
	    {
	      sub_ddmmss (nh, nl, nh, nl, dh, dl);
	      q |= 1;
	    }
	  cnt--;
	}
    }

  rp[0] = nl;
  rp[1] = nh;

  return q;
}
#else /* GMP_NAIL_BITS != 0 */
/* Two-limb division optimized for small quotients. Input words
   include nails, which must be zero. */
static inline mp_limb_t
div2 (mp_ptr rp,
      mp_limb_t nh, mp_limb_t nl,
      mp_limb_t dh, mp_limb_t dl)
{
  mp_limb_t q = 0;
  int cnt;

  ASSERT_LIMB(nh);
  ASSERT_LIMB(nl);
  ASSERT_LIMB(dh);
  ASSERT_LIMB(dl);

  /* FIXME: Always called with nh > 0 and dh >0. Then it should be
     enough to look at the high limbs to select cnt. */
  for (cnt = 0; nh > dh || (nh == dh && nl >= dl); cnt++)
  {
    dh = (dh << 1) | (dl >> (GMP_NUMB_BITS - 1));
    dl = (dl << 1) & GMP_NUMB_MASK;
  }

  while (cnt)
    {
      dl = (dh << (GMP_NUMB_BITS - 1)) | (dl >> 1);
      dh = dh >> 1;
      dl &= GMP_NUMB_MASK;

      q <<= 1;
      if (nh > dh || (nh == dh && nl >= dl))
       {
	 /* FIXME: We could perhaps optimize this by unrolling the
	    loop 2^GMP_NUMB_BITS - 1 times? */
	 nl -= dl;
	 nh -= dh;
	 nh -= (nl >> (GMP_LIMB_BITS - 1));
	 nl &= GMP_NUMB_MASK;

	 q |= 1;
       }
      cnt--;
    }
  ASSERT (nh < dh || (nh == dh && nl < dl));
  rp[0] = nl;
  rp[1] = nh;

  return q;
}
#endif /* GMP_NAIL_BITS != 0 */

#define SUB_2(w1,w0, x1,x0, y1,y0)                      \
  do {                                                  \
    ASSERT_LIMB (x1);                                   \
    ASSERT_LIMB (x0);                                   \
    ASSERT_LIMB (y1);                                   \
    ASSERT_LIMB (y0);                                   \
                                                        \
    if (GMP_NAIL_BITS == 0)                             \
      sub_ddmmss (w1,w0, x1,x0, y1,y0);                 \
    else                                                \
      {                                                 \
        mp_limb_t   __w0, __c;                          \
        SUBC_LIMB (__c, __w0, x0, y0);                  \
        (w1) = ((x1) - (y1) - __c) & GMP_NUMB_MASK;     \
        (w0) = __w0;                                    \
      }                                                 \
  } while (0)

static inline void
qstack_push_0 (struct qstack *stack)
{
  ASSERT_QSTACK (stack);

  if (stack->size_next >= QSTACK_MAX_QUOTIENTS)
    qstack_rotate (stack, 0);

  stack->size[stack->size_next++] = 0;
}

static inline void
qstack_push_1 (struct qstack *stack, mp_limb_t q)
{
  ASSERT (q >= 2);

  ASSERT_QSTACK (stack);

  if (stack->limb_next >= stack->limb_alloc)
    qstack_rotate (stack, 1);

  else if (stack->size_next >= QSTACK_MAX_QUOTIENTS)
    qstack_rotate (stack, 0);

  stack->size[stack->size_next++] = 1;
  stack->limb[stack->limb_next++] = q;

  ASSERT_QSTACK (stack);
}

/* Produce r_k from r_i and r_j, and push the corresponding
   quotient. */
#if __GMP_HAVE_TOKEN_PASTE
#define HGCD2_STEP(i, j, k) do {			\
  SUB_2 (rh ## k, rl ## k,				\
	 rh ## i, rl ## i,				\
	 rh ## j, rl ## j);				\
							\
  /* Could check here for the special case rh3 == 0,	\
     but it's covered by the below condition as well */	\
  if (       rh ## k <  rh ## j				\
      || (   rh ## k == rh ## j				\
	  && rl ## k <  rl ## j))			\
    {							\
	  /* Unit quotient */				\
	  u ## k = u ## i + u ## j;			\
	  v ## k = v ## i + v ## j;			\
							\
	  if (quotients)				\
	    qstack_push_0 (quotients);			\
	}						\
      else						\
	{						\
	  mp_limb_t r[2];				\
	  mp_limb_t q = 1 + div2 (r, rh ## k, rl ## k,	\
				     rh ## j, rl ## j);	\
	  rl ## k = r[0]; rh ## k = r[1];		\
	  u ## k = u ## i + q * u ## j;			\
	  v ## k = v ## i + q * v ## j;			\
							\
	  if (quotients)				\
	    qstack_push_1 (quotients, q);		\
	}						\
} while (0)
#else /* ! __GMP_HAVE_TOKEN_PASTE */
#define HGCD2_STEP(i, j, k) do {			\
  SUB_2 (rh/**/k, rl/**/k,				\
	 rh/**/i, rl/**/i,				\
	 rh/**/j, rl/**/j);				\
							\
  /* Could check here for the special case rh3 == 0,	\
     but it's covered by the below condition as well */	\
  if (       rh/**/k <  rh/**/j				\
      || (   rh/**/k == rh/**/j				\
	  && rl/**/k <  rl/**/j))			\
    {							\
	  /* Unit quotient */				\
	  u/**/k = u/**/i + u/**/j;			\
	  v/**/k = v/**/i + v/**/j;			\
							\
	  if (quotients)				\
	    qstack_push_0 (quotients);			\
	}						\
      else						\
	{						\
	  mp_limb_t r[2];				\
	  mp_limb_t q = 1 + div2 (r, rh/**/k, rl/**/k,	\
				     rh/**/j, rl/**/j);	\
	  rl/**/k = r[0]; rh/**/k = r[1];		\
	  u/**/k = u/**/i + q * u/**/j;			\
	  v/**/k = v/**/i + q * v/**/j;			\
							\
	  if (quotients)				\
	    qstack_push_1 (quotients, q);		\
	}						\
} while (0)
#endif /* ! __GMP_HAVE_TOKEN_PASTE */

/* Repeatedly divides A by B, until the remainder is a single limb.
   Stores cofactors in HGCD, and pushes the quotients on STACK (unless
   STACK is NULL). On success, HGCD->row[0, 1, 2] correspond to
   remainders that are larger than one limb, while HGCD->row[3]
   correspond to a remainder that fit in a single limb.

   Return 0 on failure (if B or A mod B fits in a single limb). Return
   1 if r0 and r1 are correct, but we still make no progress because
   r0 = A, r1 = B.

   Otherwise return 2, 3 or 4 depending on how many of the r:s that
   satisfy Jebelean's criterion. */
/* FIXME: There are two more micro optimizations that could be done to
   this code:

   The div2 function starts with checking the most significant bit of
   the numerator. When we call div2, that bit is know in advance for
   all but the one or two first calls, so we could split div2 in two
   functions, and call the right one.

   We could also have two versions of this code, with and without the
   quotient argument, to avoid checking if it's NULL in the middle of
   the loop. */

int
mpn_hgcd2 (struct hgcd2 *hgcd,
	   mp_limb_t ah, mp_limb_t al,
	   mp_limb_t bh, mp_limb_t bl,
	   struct qstack *quotients)
{
  /* For all divisions, we special case q = 1, which accounts for
     approximately 41% of the quotients for random numbers (Knuth,
     TAOCP 4.5.3) */

  /* Use scalar variables */
  mp_limb_t rh1, rl1, u1, v1;
  mp_limb_t rh2, rl2, u2, v2;
  mp_limb_t rh3, rl3, u3, v3;

  ASSERT_LIMB(ah);
  ASSERT_LIMB(al);
  ASSERT_LIMB(bh);
  ASSERT_LIMB(bl);
  ASSERT (ah > bh || (ah == bh && al >= bl));

  if (bh == 0)
    return 0;

  {
    mp_limb_t rh0, rl0, u0, v0;

    /* Initialize first two rows */
    rh0 = ah; rl0 = al; u0 = 1; v0 = 0;
    rh1 = bh; rl1 = bl; u1 = 0; v1 = 1;

    SUB_2 (rh2, rl2, rh0, rl0, rh1, rl1);

    if (rh2 == 0)
      return 0;

    if (rh2 < rh1 || (rh2 == rh1 && rl2 <  rl1))
      {
	/* Unit quotient */
	v2 = 1;

	if (quotients)
	  qstack_push_0 (quotients);
      }
    else
      {
	mp_limb_t r[2];
	mp_limb_t q = 1 + div2 (r, rh2, rl2, rh1, rl1);

	rl2 = r[0]; rh2 = r[1];

	if (rh2 == 0)
	  return 0;

	v2 = q;

	if (quotients)
	  qstack_push_1 (quotients, q);
      }

    u2 = 1;

    /* The simple version of the loop is as follows:
     |
     |   hgcd->sign = 0;
     |   for (;;)
     |     {
     |       (q, rh3, rl3]) = divmod (r1, r2);
     |       u[3] = u1 + q * u2;
     |       v[3] = v1 + q * v2;
     |       qstack_push_1 (quotients, q);
     |
     |       if (rh3 == 0)
     |         break;
     |
     |       HGCD2_SHIFT4_LEFT (hgcd->row);
     |       hgcd->sign = ~hgcd->sign;
     |     }
     |
     |   But then we special case for q = 1, and unroll the loop four times
     |   to avoid data movement. */

    for (;;)
      {
	HGCD2_STEP (1, 2, 3);
	if (rh3 == 0)
	  {
	    hgcd->row[0].u = u0; hgcd->row[0].v = v0;

	    hgcd->sign = 0;

	    break;
	  }
	HGCD2_STEP (2, 3, 0);
	if (rh0 == 0)
	  {
	    hgcd->row[0].u = u1; hgcd->row[0].v = v1;

	    rh1 = rh2; rl1 = rl2; u1 = u2; v1 = v2;
	    rh2 = rh3; rl2 = rl3; u2 = u3; v2 = v3;
	    rh3 = rh0; rl3 = rl0; u3 = u0; v3 = v0;

	    hgcd->sign = -1;
	    break;
	  }

	HGCD2_STEP (3, 0, 1);
	if (rh1 == 0)
	  {
	    hgcd->row[0].u = u2; hgcd->row[0].v = v2;
	    rh2 = rh0; rl2 = rl0; u2 = u0; v2 = v0;

	    MP_LIMB_T_SWAP (rh1, rh3); MP_LIMB_T_SWAP (rl1, rl3);
	    MP_LIMB_T_SWAP ( u1,  u3); MP_LIMB_T_SWAP ( v1,  v3);

	    hgcd->sign = 0;
	    break;
	  }

	HGCD2_STEP (0, 1, 2);
	if (rh2 == 0)
	  {
	    hgcd->row[0].u = u3; hgcd->row[0].v = v3;

	    rh3 = rh2; rl3 = rl2; u3 = u2; v3 = v2;
	    rh2 = rh1; rl2 = rl1; u2 = u1; v2 = v1;
	    rh1 = rh0; rl1 = rl0; u1 = u0; v1 = v0;

	    hgcd->sign = -1;
	    break;
	  }
      }
  }

  ASSERT (rh1 != 0);
  ASSERT (rh2 != 0);
  ASSERT (rh3 == 0);
  ASSERT (rh1 > rh2 || (rh1 == rh2 && rl1 > rl2));
  ASSERT (rh2 > rh3 || (rh2 == rh3 && rl2 > rl3));

  /* Coefficients to be returned */
  hgcd->row[1].u = u1; hgcd->row[1].v = v1;
  hgcd->row[2].u = u2; hgcd->row[2].v = v2;
  hgcd->row[3].u = u3; hgcd->row[3].v = v3;

  /* Rows 1, 2 and 3 are used below, rh0, rl0, u0 and v0 are not. */
#if GMP_NAIL_BITS == 0
  {
    mp_limb_t sh;
    mp_limb_t sl;
    mp_limb_t th;
    mp_limb_t tl;

    /* Check r2 */
    /* We always have r2 > u2, v2 */

    if (hgcd->sign >= 0)
      {
	/* Check if r1 - r2 >= u2 - u1 = |u2| + |u1| */
	sl = u2 + u1;
	sh = (sl < u1);
      }
    else
      {
	/* Check if r1 - r2 >= v2 - v1 = |v2| + |v1| */
	sl = v2 + v1;
	sh = (sl < v1);
      }

    sub_ddmmss (th, tl, rh1, rl1, rh2, rl2);

    if (th < sh || (th == sh && tl < sl))
      return 2 - (hgcd->row[0].v == 0);

    /* Check r3 */

    if (hgcd->sign >= 0)
      {
	/* Check r3 >= max (-u3, -v3) = |u3| */
	if (rl3 < u3)
	  return 3;

	/* Check r3 - r2 >= v3 - v2 = |v2| + |v1|*/
	sl = v3 + v2;
	sh = (sl < v2);
      }
    else
      {
	/* Check r3 >= max (-u3, -v3) = |v3| */
	if (rl3 < v3)
	  return 3;

	/* Check r3 - r2 >= u3 - u2 = |u2| + |u1| */
	sl = u3 + u2;
	sh = (sl < u2);
      }

    sub_ddmmss (th, tl, rh2, rl2, 0, rl3);

    if (th < sh || (th == sh && tl < sl))
      return 3;

    return 4;
  }
#else /* GMP_NAIL_BITS > 0 */
  {
    mp_limb_t sl;
    mp_limb_t th;
    mp_limb_t tl;

    /* Check r2 */
    /* We always have r2 > u2, v2 */

    if (hgcd->sign >= 0)
      {
       /* Check if r1 - r2 >= u2 - u1 = |u2| + |u1| */
       sl = u2 + u1;
      }
    else
      {
       /* Check if r1 - r2 >= v2 - v1 = |v2| + |v1| */
       sl = v2 + v1;
      }

    tl = rl1 - rl2;
    th = rh1 - rh2 - (tl >> (GMP_LIMB_BITS - 1));
    ASSERT_LIMB(th);

    if (th < (CNST_LIMB(1) << GMP_NAIL_BITS)
       && ((th << GMP_NUMB_BITS) | (tl & GMP_NUMB_MASK)) < sl)
      return 2 - (hgcd->row[0].v == 0);

    /* Check r3 */

    if (hgcd->sign >= 0)
      {
       /* Check r3 >= max (-u3, -v3) = |u3| */
       if (rl3 < u3)
	 return 3;

       /* Check r3 - r2 >= v3 - v2 = |v2| + |v1|*/
       sl = v3 + v2;
      }
    else
      {
       /* Check r3 >= max (-u3, -v3) = |v3| */
       if (rl3 < v3)
	 return 3;

       /* Check r3 - r2 >= u3 - u2 = |u2| + |u1| */
       sl = u3 + u2;
      }

    tl = rl2 - rl3;
    th = rh2 - (tl >> (GMP_LIMB_BITS - 1));
    ASSERT_LIMB(th);

    if (th < (CNST_LIMB(1) << GMP_NAIL_BITS)
       && ((th << GMP_NUMB_BITS) | (tl & GMP_NUMB_MASK)) < sl)
      return 3;

    return 4;
  }
#endif /* GMP_NAIL_BITS > 0 */
}

mp_size_t
mpn_hgcd2_fix (mp_ptr rp, mp_size_t ralloc,
	       int sign,
	       mp_limb_t u, mp_srcptr ap, mp_size_t asize,
	       mp_limb_t v, mp_srcptr bp, mp_size_t bsize)
{
  mp_size_t rsize;
  mp_limb_t cy;

  ASSERT_LIMB(u);
  ASSERT_LIMB(v);

  if (sign < 0)
    {
      MP_LIMB_T_SWAP (u,v);
      MPN_SRCPTR_SWAP (ap, asize, bp, bsize);
    }

  ASSERT (u > 0);

  ASSERT (asize <= ralloc);
  rsize = asize;
  cy = mpn_mul_1 (rp, ap, asize, u);
  if (cy)
    {
      ASSERT (rsize < ralloc);
      rp[rsize++] = cy;
    }

  if (v > 0)
    {
      ASSERT (bsize <= rsize);
      cy = mpn_submul_1 (rp, bp, bsize, v);
      if (cy)
	{
	  ASSERT (bsize < rsize);
	  ASSERT_NOCARRY (mpn_sub_1 (rp + bsize,
				     rp + bsize, rsize - bsize, cy));
	}

      MPN_NORMALIZE (rp, rsize);
    }
  return rsize;
}

#undef HGCD2_STEP
