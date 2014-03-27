#ifndef RANDOM__H
#define RANDOM__H

/* Copyright (C) 1991, 1996 Free Software Foundation, Inc.

   ----------------------------------------------------------
   Random numbers that return the same sequence on all platforms.
   Implementation taken from the GNU C Library. The same copyright
   applies.
   ----------------------------------------------------------

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

void portable_srand(unsigned int seed);
int portable_rand(void);

#endif /* RANDOM__H */


/* The largest number rand will return (same as INT_MAX).  */
#define RAND_MAX        2147483647
/* intentionally outside the #ifdef RANDOM__H block        */
/* You must overwrite the system RAND_MAX!                 */

