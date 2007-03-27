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

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

/***************************************************/
/***************************************************/
void die(char *errormessage) {
  /*This function prints the error message and exits
    the program.*/

  fprintf(stderr, "\a%s\n", errormessage);
  exit(1);
}

/***************************************************/
/***************************************************/
int min ( int a, int b) {
  if (a<=b) return a;
  else return b;
}

/***************************************************/
/***************************************************/

/* MEMLEAK, use free */
void *safeCalloc( int count, int size ) {
  /* this function calls calloc with the given inputs,
     but dies with an error message if a NULL is returned */

  void *newthing=calloc( count, size );
  if (newthing==NULL) {
    die("calloc returned NULL");
    return NULL; /* unreachable. */
  }
  else return newthing;
}

/***************************************************/
/***************************************************/

/* MEMLEAK, use free */
void *safeMalloc( int count, int size ) {
  /* this function calls malloc with the inputs, which are
     to be provided in calloc notation. If the result is
     NULL, the program dies with an error message.*/

  void *newthing=malloc( count*size );
  if (newthing==NULL) {
    die("malloc returned NULL");
    return NULL; /* unreachable */
  }
  else return newthing;
}

/***************************************************/
/***************************************************/

/* MEMLEAK (see destroyMatrix) */
matrix *createMatrix( int rows, int cols) {
  matrix *newmatrix=(matrix *)safeMalloc(1,sizeof(matrix));

  bit *newcells=(bit *)safeCalloc(rows*cols,sizeof(bit));
  int *newrowswap=(int *)safeMalloc(rows, sizeof(int));
  int *newcolswap=(int *)safeMalloc(cols, sizeof(int));
  int i;

  newmatrix->cells=newcells;
  newmatrix->rowswap=newrowswap;
  newmatrix->colswap=newcolswap;

  newmatrix->rows=rows;
  newmatrix->cols=cols;

  for (i=0; i<rows; i++) { newmatrix->rowswap[i]=i; }
  for (i=0; i<cols; i++) { newmatrix->colswap[i]=i; }

  return newmatrix;
}

/***************************************************/
/***************************************************/
void destroyMatrix(matrix *condemned) {
  if (condemned==NULL) {
    printf("\aTried to free a NULL data structure.\n");
    return;
  }
  free(condemned->cells);
  free(condemned->rowswap);
  free(condemned->colswap);
  free(condemned);
  return;
}

/***************************************************/
/***************************************************/
void swapRow( matrix *m, int a, int b) {
  int temp=m->rowswap[a];
  m->rowswap[a]=m->rowswap[b];
  m->rowswap[b]=temp;
}

/***************************************************/
/***************************************************/
void swapCol( matrix *m, int a, int b) {
  int temp=m->colswap[a];
  m->colswap[a]=m->colswap[b];
  m->colswap[b]=temp;
}

/***************************************************/
/***************************************************/

void addTrueCell(matrix *thematrix, int srow, int scol, int drow, int dcol) {
  int sindex=srow * (thematrix->cols) + scol;
  int dindex=drow * (thematrix->cols) + dcol;

  thematrix->cells[dindex]=
    (thematrix->cells[dindex]+thematrix->cells[sindex]) % 2;
}

/***************************************************/
/***************************************************/

void writeTrueCell(matrix *thematrix, int row, int col, bit value) {
  int index=row * (thematrix->cols) + col;

  thematrix->cells[index]=value;
}

/***************************************************/
/***************************************************/

void addCell(matrix *thematrix, int srow, int scol, int drow, int dcol) {
  int newsrow=thematrix->rowswap[srow];
  int newscol=thematrix->colswap[scol];
  int newdrow=thematrix->rowswap[drow];
  int newdcol=thematrix->colswap[dcol];

  addTrueCell(thematrix, newsrow, newscol, newdrow, newdcol);
}

/***************************************************/
/***************************************************/

void writeCell(matrix *thematrix, int row, int col, bit value) {
  int newrow=thematrix->rowswap[row];
  int newcol=thematrix->colswap[col];

  writeTrueCell(thematrix, newrow, newcol, value);
}
/***************************************************/
/***************************************************/

bit readTrueCell(matrix *thematrix, int row, int col) {
  int index=row * (thematrix->cols) + col;

  return thematrix->cells[index];
}
/***************************************************/
/***************************************************/

bit readCell(matrix *thematrix, int row, int col) {
  int newrow=thematrix->rowswap[row];
  int newcol=thematrix->colswap[col];

  return readTrueCell(thematrix, newrow, newcol);
}

/***************************************************/
/***************************************************/

void zeroMatrix(matrix *thematrix) {
  int i, j;

  for (i=0; i<thematrix->rows; i++) {
     for (j=0; j<thematrix->cols; j++) {
       writeCell(thematrix, i, j, 0);
    }
  }
}

/***************************************************/
/***************************************************/

void printMatrix(matrix *thematrix) {
  int i, j;

  for (i=0; i<thematrix->rows; i++) {
    printf("[ ");
    for (j=0; j<thematrix->cols; j++) {
      printf("%u ", readCell(thematrix, i,j) );
    }
    printf("]\n");
  }

  printf("\n");
}

/***************************************************/
/***************************************************/

bit coinFlip() {


  if (rand() < RAND_MAX/2) {
    return 0;
  }  else {
    return 1;
  }
}

/***************************************************/
/***************************************************/

void fillRandomMatrix(matrix *thematrix) {
  int i, j;

  for (i=0; i<thematrix->rows; i++) {
     for (j=0; j<thematrix->cols; j++) {
       writeCell(thematrix, i, j, coinFlip());
    }
  }
}

/***************************************************/
/***************************************************/

void setIdentityMatrix(matrix *thematrix) {
  int i;
  int limit=min(thematrix->rows, thematrix->cols);

  zeroMatrix(thematrix);

  for (i=0; i<limit; i++) {
    writeCell(thematrix, i, i, 1);
  }
}

/***************************************************/
/***************************************************/

/* MEMLEAK, use destroyMatrix */
matrix *concatMatrix(matrix *left, matrix *right) {
  int i,j,newwidth;
  matrix *newmatrix;

  if (left->rows != right->rows) {
    printf("Attempted to concat matrix of height %d with height %d.",
	   left->rows, right->rows);
    die("Terminating");
  }

  newwidth= left->cols + right->cols;

  newmatrix=createMatrix(left->rows, newwidth);

  for (i=0; i<left->rows; i++) {
    for (j=0; j<newwidth; j++) {
      if (j<left->cols) {
	writeCell(newmatrix, i, j, readCell(left, i, j));
      } else {
	writeCell(newmatrix, i, j, readCell(right, i, j-left->cols));
      }
    }
  }

  return newmatrix;
}

/***************************************************/
/***************************************************/

int forceNonZero(matrix *m, int x, int y) {
  int i;
  int j=y;

  //for (j=y; j<m->cols; j++) {
    for (i=x; i<m->rows; i++) {
      if (readCell(m, i, j)==1) {
	if (i!=x) swapRow(m, i, x);
	if (j!=y) swapCol(m, j, y);
	return YES;
      }
    }
    //}

  return NO;
}

/***************************************************/
/***************************************************/

void rowAdd(matrix *m, int rowa, int rowb) {
  int i;

  for (i=0; i<m->cols; i++) {
    /*    writeCell(m, rowb, i,
	  (readCell(m, rowa, i) + readCell(m, rowb, i)) % 2 );*/
    addCell(m, rowa, i, rowb, i);
  }
}

/***************************************************/
/***************************************************/

int delayedGaussian(matrix *m, int firstcol, int full) {
  int pc; /* pivot column */
  int tr; /* target row */
  int stop=min(m->rows, m->cols);
  int good;

  for (pc=firstcol; pc<stop; pc++) {
    /* Step one, find a pivot row in this column.*/
    good=forceNonZero(m, pc, pc);
    if (good==NO) return NO;

    for (tr=( full ? 0 : pc+1 ); tr<m->rows; tr++) {
      /* Step two, add this pivot row to other rows as needed. */
      if (tr==pc) continue;

      if (readCell(m, tr, pc)==0) continue;

      rowAdd(m, pc, tr);
    }
  }

  return YES;
}

int gaussian(matrix *m, int full) { return delayedGaussian(m,0, full); }

/***************************************************/
/***************************************************/

/* MEMLEAK, use destroyMatrix() */
matrix *copySubMatrix(matrix *m, int rlow, int clow, int rhigh, int chigh) {
  int i,j, newrows, newcols;
  matrix *newmatrix;

  if (rlow>rhigh) die("Bad args to subMatrix.\n");
  if (clow>chigh) die("Bad args to subMatrix.\n");
  if (clow>m->cols) die("Bad args to subMatrix.\n");
  if (chigh>m->cols) die("Bad args to subMatrix.\n");
  if (rlow>m->rows) die("Bad args to subMatrix.\n");
  if (rhigh>m->rows) die("Bad args to subMatrix.\n");

  newrows=rhigh-rlow+1;
  newcols=chigh-clow+1;

  newmatrix=createMatrix(newrows, newcols);

  for (i=0; i<newrows; i++) {
    for (j=0; j<newcols; j++) {
     	writeCell(newmatrix, i, j, readCell(m, i+rlow, j+clow));
    }
  }

  return newmatrix;
}

/***************************************************************/

/* MEMLEAK, use destroyMatrix */
matrix *cloneMatrix( matrix *a ) {
  return copySubMatrix( a, 0, 0, (a->rows)-1, (a->cols)-1);
}

/***************************************************************/
/***************************************************************/

/* MEMLEAK, use destroyMatrix */
matrix *matrixTimesMatrix(matrix *a, matrix *b) {
  int i,j,k, total;
  matrix *newmatrix;

  if (a->cols!=b->rows) die("Incompatible Dimensions in Matrix Multiplic.");

  newmatrix=createMatrix(a->rows, b->cols);

  for (i=0; i<a->rows; i++) {
    for (j=0; j<b->cols; j++) {
      total=0;
      for (k=0; k<a->cols; k++) {
	total+= readCell(a, i, k) * readCell(b, k, j);
      }
      writeCell(newmatrix, i, j, total % 2);
    }
  }

  return newmatrix;
}

/***************************************************************/
/***************************************************************/

int equalMatrix(matrix *a, matrix *b) {
  int i, j;

  if (a->rows!=b->rows) return NO;
  if (a->cols!=b->cols) return NO;

  for (i=0; i<a->rows; i++) {
    for (j=0; j<a->cols; j++) {
      if (readCell(a,i,j) != readCell(b,i,j)) return NO;
    }
  }

  return YES;
}

/***************************************************************/
/***************************************************************/


int forceNonZero2(matrix *m, int xstart, int xstop, int y) {
  int i;
  int j=y;

  //for (j=y; j<m->cols; j++) {
    for (i=xstart; i<=xstop; i++) {
      if (readCell(m, i, j)==1) {
	if (i!=xstart) swapRow(m, i, xstart);
	if (j!=y) swapCol(m, j, y);
	return YES;
      }
    }
    //}

  return NO;
}

/***************************************************/
/***************************************************/

int prep(matrix *m, int homerow, int homecol, int k) {
  int pc; /* pivot column */
  int tr; /* target row */
  //int stop=min(m->rows, m->cols);
  int good;

  for (pc=homecol; pc<min(homecol+k,m->cols); pc++) {
    /* Step one, find a pivot row in this column.*/
    good=forceNonZero2(m, pc, min(homerow+3*k-1,m->rows-1), pc);
    if (good==NO) return NO;

    for (tr=homerow; tr<min(homerow+3*k,m->rows); tr++) {
      /* Step two, add this pivot row to other rows as needed. */
      if (tr==pc) continue;

      if (readCell(m, tr, pc)==0) continue;

      rowAdd(m, pc, tr);
    }
  }

  return YES;
}


/*********************************/

bit *table;
int *lookup;
int numcols;

void makeTable (matrix *m, int homerow, int homecol, int k) {
  int i,j;
  int row, id;
  int twok= 1<<k;

  numcols=m->cols - homecol;
  table=(bit *)safeMalloc(twok*numcols, sizeof(bit));
  lookup=(int *)safeMalloc(twok, sizeof(int));

  for (i=0; i<numcols; i++) {
    table[i + 0*numcols]=(bit)0;
    lookup[0]=0;
  }

  for (i=1; i<twok; i++) {
    row=codebook[k]->inc[i-1]+homerow;
    id=codebook[k]->ord[i];

    lookup[id]=i;

    for  (j=0; j<numcols; j++) {
      table[j + i*numcols] = (table[j + (i-1)*numcols] +
			      readCell(m, row, j+homecol)) % 2;
    }
  }

#ifdef VERBOSE
  printf("\n Gray Table: \n");

  for (i=0; i<twok; i++) {
    printf("%2d: [ ", i);
    for (j=0; j<numcols; j++) {
      printf("%d ", table[j + i*numcols]);
    }
    printf("]\n");
  }

  printf("\nIf you want ---> Ask for\n");

  for (i=0; i<twok; i++) {
    printf("%d ---> %d\n", i, lookup[i]);
  }
#endif
}

/*********************************/
void processRow(matrix *m, int row, int homecol, int k) {
  int i, num=0;
  int tablerow;

  for (i=min(homecol+k-1, m->cols); i>=homecol; i--) {
    num*=2;
    num+=(int)readCell(m, row, i);
  }

  /* printf("I think row %d is a %d.\n", row, num);  */

  tablerow=lookup[num];

  for (i=homecol; i<m->cols; i++) {
    writeCell(m, row, i,
	    (readCell(m, row, i) + table[i-homecol + tablerow*numcols]) % 2);
  }
}

/*********************************/
void process(matrix *m, int startrow, int stoprow, int homecol, int k) {
  int i;

  /* printf("starting row work\n"); */
  for (i=startrow; i<=stoprow; i++) {
    processRow(m, i, homecol, k);
  }
}

/*********************************/

int fourRussians(matrix *m, int k, int full) {
  int ell=min(m->rows, m->cols);
  int i, fullrank;
  //int twok=powersoftwo[k];

  for (i=0; i<ell; i+=k) {
    if (i+k-1>=m->rows) {
      printf("Not enough room for the submatrix, switching to Gaussian Elimination!\n");

      return delayedGaussian(m, i, full);
    }

    fullrank=prep(m,i,i,k);

    if (fullrank==NO) {
      printf("Aborting, submatrix not full rank.\n");
      return NO;
    } else;

    makeTable(m, i, i, k);

    process(m, 3*k+i, m->rows-1, i, k);

    if (full==YES) process(m, 0, i-1, i, k);

    free(table);
    free(lookup);
  }

  return YES;
}









