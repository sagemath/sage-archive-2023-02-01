#ifndef MATRIX_H
#define MATRIX_H

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

#include "grayflex.h"

/** code up: addToCell for XOR, like readCell and writeCell. **/

#define YES 1
#define NO 0

typedef unsigned char /* renamed as */ bit;

struct matrixstruct {
  bit *cells;
  int *rowswap;
  int *colswap;
  int rows;
  int cols;
};

typedef struct matrixstruct /* renamed as */ matrix;

extern bit *table;
extern int *lookup;
extern int numcols;

void die(char *errormessage);

int min ( int a, int b);

/* MEMLEAK, use free */
void *safeCalloc( int count, int size );

/* MEMLEAK, use free */
void *safeMalloc( int count, int size );

/* MEMLEAK (see destroyMatrix) */
matrix *createMatrix( int rows, int cols);

void destroyMatrix(matrix *condemned);

void swapRow( matrix *m, int a, int b);

void swapCol( matrix *m, int a, int b);

void addTrueCell(matrix *thematrix, int srow, int scol, int drow, int dcol);

void writeTrueCell(matrix *thematrix, int row, int col, bit value);

void addCell(matrix *thematrix, int srow, int scol, int drow, int dcol);

void writeCell(matrix *thematrix, int row, int col, bit value);

bit readTrueCell(matrix *thematrix, int row, int col);

bit readCell(matrix *thematrix, int row, int col);

void zeroMatrix(matrix *thematrix);

void printMatrix(matrix *thematrix);

bit coinFlip();

void fillRandomMatrix(matrix *thematrix);

void setIdentityMatrix(matrix *thematrix);

matrix *concatMatrix(matrix *left, matrix *right);

int forceNonZero(matrix *m, int x, int y);

void rowAdd(matrix *m, int rowa, int rowb);


int delayedGaussian(matrix *m, int firstcol, int full);

int gaussian(matrix *m, int full);

/* MEMLEAK, use destroyMatrix() */
matrix *copySubMatrix(matrix *m, int rlow, int clow, int rhigh, int chigh);

/* MEMLEAK, use destroyMatrix */
matrix *cloneMatrix( matrix *a );

/* MEMLEAK, use destroyMatrix */
matrix *matrixTimesMatrix(matrix *a, matrix *b);

int equalMatrix(matrix *a, matrix *b);

int forceNonZero2(matrix *m, int xstart, int xstop, int y);

int prep(matrix *m, int homerow, int homecol, int k);


void makeTable (matrix *m, int homerow, int homecol, int k);

void processRow(matrix *m, int row, int homecol, int k);

void process(matrix *m, int startrow, int stoprow, int homecol, int k);

int fourRussians(matrix *m, int k, int full);

#endif //MATRIX_H
