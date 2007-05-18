#ifndef GRAYFLEX_H
#define GRAYFLEX_H

/**
 * gray code generation used by the M4RI algorithm.
 *
 * AUTHOR: malb
 */

/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
*       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
*
*  Distributed under the terms of the GNU General Public License (GPL)
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/

struct codestruct {
  int *ord;
  int *inc;
};

typedef struct codestruct /*renamed as*/ code;

extern code **codebook;

#define TWOPOW(i) (1<<(i))

void printBitString(int number, int length);

/**
 * swaps length bits in v naively.
 *
 * WARNING: Uppper bits of return value may contain garbage after
 * operation.
 */
int swap_bits(int v,int length);

/**
 * Returns the 'number'-th gray code entry for a gray code of length
 * $2^{length}$.
 *
 * INPUT:
 *     number -- index in the gray code table
 *     length -- length of the gray code
 *
 * OUTPUT:
 *      number-th gray code entry
 *
 * AUTHOR: malb
 *
 * THANKS: Soroosh Yazdani explained the repeated sum idea to me.
 *
 */

int grayCode(int number, int length);

/**
 * Fills in 'ord' and 'inc' with gray code data for a gray code of
 * length $2^{length}$.
 *
 * INPUT:
 *    ord -- will hold gray code data, must be preallocated with correct size
 *    inc -- will hold some increment data, must be preallocated with correct size
 *
 * AUTHOR: malb
 *
 * THANKS: Robert Miller had the idea for a non-recursive implementation.
 *
 */

void buildCodeFlex(int *ord, int *inc, int length);

void buildAllCodes();

void destroyAllCodes();

#endif //GRAYFLEX_H
