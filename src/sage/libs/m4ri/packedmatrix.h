#ifndef PACKEDMATRIX_H
#define PACKEDMATRIX_H

/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
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

#include "matrix.h"
#include <stdio.h>

#define RADIX 64
#define SAFECHAR 85 /* Radix*1.25 plus a few. */
#define ONE 1ULL


typedef unsigned long long word;

struct packedmatrixstruct {
  word *values;

  int rows;
  int cols;
  int width; /* rounded as floor(rows/RADIX)+1. */

  int *rowswap;

};

typedef struct packedmatrixstruct packedmatrix;

extern word packingmask[RADIX];
extern word bytemask[RADIX/8];
extern word sixteenmask[RADIX/16];

packedmatrix *createPackedMatrix(int r, int c);

void destroyPackedMatrix( packedmatrix *condemned );

inline void rowSwapPacked( packedmatrix *m, int rowa, int rowb );

/* Internal: do not call */
void setupPackingMasks();

inline bit readPackedCell( packedmatrix *m, int row, int col );

inline void writePackedCell( packedmatrix *m, int row, int col, bit value);

/* Keep in mind that the row, col refer to a row and column (of bits), and
   you can address the block by any of the RADIX (usually 64) A_ijs there. */

inline void xorPackedBlock( packedmatrix *m, int row, int col, word value);

/* Keep in mind that the row, col refer to a row and column (of bits), and
   you can address the block by any of the RADIX (usually 64) A_ijs there. */

inline void writePackedBlock( packedmatrix *m, int row, int col, word value);

word fetchByte( word data, int which );

/* The RADIX-bit word can be divided into RADIX/8 octets or 8 bit bytes */
/* The most significant byte is numbered 0. */
word fetch16( word data, int which );

/* Warning: I assume *destination has RADIX+1 bytes available */

void wordToString( char *destination, word data);

/* Warning: I assume *destination has 9 bytes available */

void byteToString( char *destination, word data);

/* Warning: I assume *destination has RADIX*1.25 bytes available */

void wordToStringComma( char *destination, word data);

/* Important note: You can specify any of the RADIX bits (64 bits usually),
   inside of the block, and it will still return the correct entire block */
inline word readPackedBlock( packedmatrix *m, int row, int col );

void printPackedMatrix( packedmatrix *m );

void printPackedMatrixTight( packedmatrix *m );

/**********************************************************************/
/* this adds rows sourcerow and destrow and stores the total in row
   destrow, but only begins at the column coloffset */

void rowAddPackedOffset( packedmatrix *m, int sourcerow, int destrow,
			 int coloffset );

void rowClearPackedOffset(packedmatrix *m, int row, int coloffset);

void rowAddPacked( packedmatrix *m, int sourcerow, int destrow);

int gaussianPackedDelayed(packedmatrix *m, int startcol, int full);

int gaussianPacked(packedmatrix *m, int full);

inline bit dotProductPacked( word a, word b );

packedmatrix *transposePacked( packedmatrix *data );

bit bigDotProductPacked( packedmatrix *a, packedmatrix *bT, int rowofa,
			 int rowofb );

/* MEMLEAK use destroyPackedMatrix */

packedmatrix *matrixTimesMatrixTransposePacked( packedmatrix *a,
						packedmatrix *bT );

/* MEMLEAK: use destroyPackedMatrix */
/* Normally, if you will multiply several times by b, it is smarter to
  calculate bT yourself, and keep it, and then use the function called
  matrixTimesMatrixTransposePacked */
packedmatrix *matrixTimesMatrixPacked( packedmatrix *a,
				       packedmatrix *b);
word randomWord();

void fillRandomlyPacked( packedmatrix *a );

void makeIdentityPacked( packedmatrix *a );

matrix *unpackMatrix( packedmatrix *pm );

packedmatrix *packMatrix( matrix *m );

inline bit isEqualWord( word a, word b);

bit equalPackedMatrix( packedmatrix *a, packedmatrix *b );

int comparePackedMatrix(packedmatrix *a, packedmatrix *b);

/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *clonePacked( packedmatrix *p);

/* MEMLEAK: use destroyPackedMatrix */
/* This is sometimes called augment */
packedmatrix *concatPacked( packedmatrix *a, packedmatrix *b);


/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *copySubMatrixPacked( packedmatrix *a, int lowr, int lowc,
				   int highr, int highc);

/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *invertPackedGaussian(packedmatrix *target,
				   packedmatrix *identity);

/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *addPacked(packedmatrix *left, packedmatrix *right);

void lazyPrint(packedmatrix *a);


#endif //PACKEDMATRIX_H
