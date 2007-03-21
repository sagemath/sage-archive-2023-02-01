#ifndef BRILLIANTRUSSIAN_H
#define BRILLIANTRUSSIAN_H

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


#include <math.h>
#include <string.h>
#include "packedmatrix.h"

#define TWOPOW(i) (1<<(i))

/**
 * Finds a pivot row between xstart and xstop. The column where this
 * pivot is search is y. Returns YES if such a pivot row was
 * found. Also, the appropriate row is swapped to the top (== xstart).
 *
 * INPUT:
 *     m -- matrix to operate on
 *     xstart -- start row
 *     xstop  -- stop row (including)
 *     y -- column to read
 *
 * OUTPUT:
 *     YES if a pivot row was found
 */

int forceNonZero2PackedFlex(packedmatrix *m, int xstart, int xstop, int y);


/***************************************************/

/**
 * Performs Gaussian elimination on a submatrix of 3k x k starting at
 * point (homepoint, homepoint) of m.
 *
 * INPUT:
 *     m -- matrix to operate on
 *     homepoint -- row,col where to start
 *     k -- the parameter k of M4RI
 *
 * OUTPUT:
 *     returns rank of 3k x k submatrix.
 */

int prepPackedFlex(packedmatrix *m, int ai, int k);

/** Adds row1 of s1, starting with col1 to the end, to row2 of s2,
 *  starting with col2 to the end. This gets stored in dest, in row3,
 *  starting with col3
 *
 *  row3[col3:] = row1[col1:] + row2[col2:]
 *
 */

void combineFlex( packedmatrix * s1, int row1, int startblock1,
	          packedmatrix * s2, int row2, int startblock2,
	          packedmatrix * dest, int row3, int startblock3 );

/*
 * Constructs all possible 2^k row combinations using the gray code
 * table.
 *
 * INPUT:
 *     m -- matrix to operate on
 *     ai -- the starting position
 *     k -- the k parameter of M4RI
 *     tablepacked -- prealloced matrix of dimension 2^k x m->cols
 *     lookuppacked -- prealloced table of length 2^k
 *     full -- touch columns before ai?
 *
 */

void makeTablePackedFlex( packedmatrix *m, int ai, int k, packedmatrix *tablepacked, int *lookuppacked, int full);

/**
 * returns k bits/entries starting at m[x,y].
 *
 * WARNING: Precondition is k < RADIX
 */

inline int getValueFlex(packedmatrix *m, int x, int y, int k);

/**
 *
 * Adds the correct row from tablepacked to the row 'row' in 'm' starting at homecol.
 *
 * INPUT:
 *     m -- matrix to operate on
 *     row -- the row which is operated on
 *     homecol -- starting column for addition
 *     k -- M4RI parameter, used for gray table lookup
 *     tablepacked -- contains the correct row to be added
 *     lookuptable -- contains row number to be addede
 */

void processRowPackedFlex(packedmatrix *m, int row, int homecol, int k, packedmatrix *tablepacked, int *lookuppacked);

/*
 *  Iterates proccessPackedRowFlex from startrow to stoprow.
 */

void processPackedFlex(packedmatrix *m, int startrow, int stoprow, int startcol, int k, packedmatrix *tablepacked, int *lookuppacked);

/**
 * This is the actual heart of the M4RI algorithm.
 *
 */

int doAByteColumnFlex(packedmatrix *m, int full, int k, int ai, packedmatrix *tablepacked, int *lookuppacked);

/*
 *
 *
 */

int fourRussiansPackedFlex(packedmatrix *m, int full, int k, packedmatrix *tablepacked, int *lookuppacked);

/**
 * Top level function to call for M4RI
 */

int simpleFourRussiansPackedFlex(packedmatrix *m, int full, int k);


/**
 * Inverts the matrix m using the M4RI algorithm. To avoid recomputing
 * the identity matrix over and over again, I may be passed in as
 * identity parameter.
 */

packedmatrix *invertPackedFlexRussian(packedmatrix *m, packedmatrix *identity, int k);


packedmatrix *m4rmPacked(packedmatrix *A, packedmatrix *B, int k);

packedmatrix *m4rmTransposePacked(packedmatrix *A, packedmatrix *B, int k);

#endif //BRILLIANTRUSSIAN_H
