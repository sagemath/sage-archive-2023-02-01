/* qstack.c

   THE FUNCTIONS IN THIS FILE ARE INTERNAL WITH MUTABLE INTERFACES.  IT IS ONLY
   SAFE TO REACH THEM THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST
   GUARANTEED THAT THEY'LL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.

Copyright 2003 Free Software Foundation, Inc.

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

mp_size_t
qstack_itch (mp_size_t size)
{
  /* Limit on the recursion depth */

  unsigned k = mpn_hgcd_max_recursion (size);

  ASSERT (2 * k < QSTACK_MAX_QUOTIENTS);

  return 2 * size + 12 * (k+1);
}

void
qstack_reset (struct qstack *stack,
	      mp_size_t asize)
{
  /* Limit on the recursion depth */
  unsigned k = mpn_hgcd_max_recursion (asize);

  stack->size_next = 0;
  stack->limb_next= 0;
  stack->nkeep = 2 * (k + 1);
  ASSERT (stack->nkeep < QSTACK_MAX_QUOTIENTS);

  ASSERT_QSTACK (stack);
}

void
qstack_init (struct qstack *stack,
	     mp_size_t asize,
	     mp_limb_t *limbs, mp_size_t alloc)
{
  stack->limb = limbs;
  stack->limb_alloc = alloc;
  /* Should depend on input size, we need 2 * recursion depth */
  stack->nkeep = QSTACK_MAX_QUOTIENTS - 1;

  qstack_reset (stack, asize);
}

/* Drop all but the nkeep latest quotients. Drop additional quotients
   if needed to reclain at least SIZE limbs of storage. */
void
qstack_rotate (struct qstack *stack,
	       mp_size_t size)
{
  unsigned dropped_sizes;
  unsigned kept;
  unsigned i;

  mp_size_t dropped_limbs;

  ASSERT_QSTACK (stack);

  if (stack->size_next > stack->nkeep)
    dropped_sizes = stack->size_next - stack->nkeep;
  else
    dropped_sizes = 0;

  for (i = 0, dropped_limbs = 0; i < dropped_sizes; i++)
    dropped_limbs += stack->size[i];

  for (; dropped_limbs < size; dropped_sizes++)
    {
      ASSERT (dropped_sizes < stack->size_next);
      dropped_limbs += stack->size[dropped_sizes];
    }

  ASSERT (dropped_limbs <= stack->limb_next);

  kept = stack->size_next - dropped_sizes;

  if (dropped_sizes)
    /* memmove isn't portable */
    for (i = 0; i < kept; i++)
      stack->size[i] = stack->size[i + dropped_sizes];

  stack->size_next = kept;

  if (dropped_limbs)
    {
      if (dropped_limbs < stack->limb_next)
	{
	  MPN_COPY_INCR (stack->limb, stack->limb + dropped_limbs,
			stack->limb_next - dropped_limbs);
	  ASSERT (dropped_limbs <= stack->limb_next);
	  stack->limb_next -= dropped_limbs;
	}
      else
	stack->limb_next = 0;
    }
  ASSERT_QSTACK (stack);
}

#if WANT_ASSERT
void
__gmpn_qstack_sanity (struct qstack *stack)
{
  mp_size_t next;
  unsigned i;

  for (i = 0, next = 0; i < stack->size_next; i++)
    {
      mp_size_t qsize = stack->size[i];
      ASSERT (next <= stack->limb_alloc);

      ASSERT (qsize == 0 || stack->limb[next+qsize - 1] != 0);
      next += qsize;
    }

  ASSERT (next == stack->limb_next);
}
#endif
