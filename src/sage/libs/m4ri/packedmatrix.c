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

#include <stdlib.h>
#include <string.h>
#include "packedmatrix.h"


word packingmask[RADIX];
word bytemask[RADIX/8];
word sixteenmask[RADIX/16];

/****************************************/

/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *createPackedMatrix(int r, int c) {
  packedmatrix *newmatrix;
  int i;

  newmatrix=(packedmatrix *)calloc(1, sizeof(packedmatrix));

  if ((c % RADIX)==0) newmatrix->width=(c/RADIX);
  else newmatrix->width=(c/RADIX) + 1;

  newmatrix->cols=c;
  newmatrix->rows=r;

  newmatrix->values=(word *)safeCalloc( (newmatrix->width)*r, sizeof(word));

  newmatrix->rowswap=(int *)safeMalloc( r, sizeof(int));

  // Rowswap does not contain the rowswap index i but the correct
  // offset in the values table. Rowswap is exclusively used to access
  // elements in that table and this speeds up computation a little. (malb)

  for (i=0; i<r; i++) { newmatrix->rowswap[i]=i*(newmatrix->width); }

  return newmatrix;
}

/************************************************************/

void destroyPackedMatrix( packedmatrix *condemned) {
  free(condemned->values);
  free(condemned->rowswap);
  free(condemned);
}


/************************************************************/

/* Internal: do not call */
void setupPackingMasks() {
  int i, j;
  word x=1;

  for (i=RADIX-1; i>=0; i--) {
    packingmask[i]=x;
    x<<=1;
    /* printf("%016llx\n",packingmask[i]);*/
  }

  for (i=0; i<RADIX/8; i++) {
    x=0;
    for (j=0; j<8; j++) {
      x|=packingmask[j+i*8];
    }
    bytemask[i]=x;
    /*printf("%016llx\n",bytemask[i]);*/
  }


  for (i=0; i<RADIX/16; i++) {
    x=0;
    for (j=0; j<16; j++) {
      x|=packingmask[j+i*16];
    }
    sixteenmask[i]=x;
    /*printf("%016llx\n",bytemask[i]);*/
  }

}

/************************************************************/

inline void rowSwapPacked( packedmatrix *m, int rowa, int rowb ) {
  int temp=m->rowswap[rowa];
  m->rowswap[rowa]=m->rowswap[rowb];
  m->rowswap[rowb]=temp;
}

/************************************************************/

inline bit readPackedCell( packedmatrix *m, int row, int col ) {
  int block=col/RADIX;
  int spot=col % RADIX;
  int truerow=m->rowswap[row];

  word entry=m->values[ block + truerow ];

  word resolved=entry & ((ONE)<<(RADIX - spot - 1));

  return (resolved >> (RADIX - spot -1));
}

/************************************************************/

inline void writePackedCell( packedmatrix *m, int row, int col, bit value) {
  int block=col/RADIX;
  int spot=col % RADIX;
  int truerow=m->rowswap[row];

  if (value==0) {
    m->values[ block + truerow ] &= ~((ONE) <<(RADIX - spot - 1));
  } else {
    m->values[ block + truerow ] |= ((ONE)<<(RADIX - spot - 1));
  }
}

/**********************************************************************/

/* Keep in mind that the row, col refer to a row and column (of bits), and
   you can address the block by any of the RADIX (usually 64) A_ijs there. */
inline void xorPackedBlock( packedmatrix *m, int row, int col, word value) {
  int block=col/RADIX;
  int truerow=m->rowswap[row];

  word *entry=m->values + block + truerow;
  *entry ^= value;
}


/**********************************************************************/

/* Keep in mind that the row, col refer to a row and column (of bits), and
   you can address the block by any of the RADIX (usually 64) A_ijs there. */
inline void writePackedBlock( packedmatrix *m, int row, int col, word value) {
  int block=col/RADIX;
  int truerow=m->rowswap[row];

  m->values[ block + truerow] = value;
}

/**********************************************************************/
/* The RADIX-bit word can be divided into RADIX/8 octets or 8 bit bytes */
/* The most significant byte is numbered 0. */
word fetchByte( word data, int which ) {
  word masked=data & bytemask[which];
  int chunks=RADIX/8;
  int moved=(chunks-which-1)*8;

  masked>>=moved;

  return masked;
}

/**********************************************************************/
/* The RADIX-bit word can be divided into RADIX/8 octets or 8 bit bytes */
/* The most significant byte is numbered 0. */
word fetch16( word data, int which ) {
  word masked=data & sixteenmask[which];
  int chunks=RADIX/16;
  int moved=(chunks-which-1)*16;

  masked>>=moved;

  return masked;
}

/**********************************************************************/

/* Warning: I assume *destination has RADIX+1 bytes available */

/* Works */
void wordToString( char *destination, word data) {
  int i;

  for (i=0; i<RADIX; i++) {
    if ((data & packingmask[i]) == 0) {
      destination[i]='0';
    } else destination[i]='1';
  }

  destination[RADIX]='\0';
}

/* Warning: I assume *destination has 9 bytes available */

/* Works */
void byteToString( char *destination, word data) {
  int i;

  for (i=0; i<8; i++) {
    if ((data & packingmask[i+RADIX-8]) == 0) {
      destination[i]='0';
    } else destination[i]='1';
  }

  destination[8]='\0';
}

/**********************************************************************/

/* Warning: I assume *destination has RADIX*1.25 bytes available */

/* Works */
void wordToStringComma( char *destination, word data) {
  int i;
  int j=0;

  for (i=0; i<RADIX; i++) {
    if ((data & packingmask[i]) == 0) {
      destination[j]='0';
    } else destination[j]='1';
    j++;

    if (((i % 4)==3) && (i!=RADIX-1)) {
      destination[j]=':';
      j++;
    }
  }

  destination[(int)(RADIX*1.25)-1]='\0';
}

/************************************************************/

/* Works */
/* Important note: You can specify any of the RADIX bits (64 bits usually),
   inside of the block, and it will still return the correct entire block */
inline word readPackedBlock( packedmatrix *m, int row, int col ) {
  int block=col/RADIX;
  int truerow=m->rowswap[row];

  word entry=m->values[ block + truerow ];

  return entry;
}

/**********************************************************************/

/* Works */
void printPackedMatrix( packedmatrix *m ) {
  int i, j;
  char temp[SAFECHAR];
  word block;

  for (i=0; i< m->rows; i++ ) {
    printf("[ ");

    for (j=0; j< m->cols; j+=RADIX) {
      block=readPackedBlock(m, i, j);
      wordToStringComma(temp, block);
      printf("%s ", temp);
    }
    printf("]\n");
  }
}


/**********************************************************************/

/* Works */
void printPackedMatrixTight( packedmatrix *m ) {
  int i, j;
  char temp[SAFECHAR];
  word block;

  for (i=0; i< m->rows; i++ ) {
    printf("[");

    for (j=0; j< m->cols; j+=RADIX) {
      block=readPackedBlock(m, i, j);
      wordToString(temp, block);
      printf("%s", temp);
    }
    printf("]\n");
  }

  printf("\n\n\n");
}

/**********************************************************************/
/* this clears the row, but only begins at the column coloffset */
void rowClearPackedOffset(packedmatrix *m, int row, int coloffset) {
  int startblock= coloffset/RADIX;
  int i;
  word temp;

  // make sure to start clearing at coloffset
  if (coloffset%RADIX) {
    temp=readPackedBlock(m, row, coloffset);
    temp &=  ~(((ONE<<(RADIX-coloffset%RADIX))) - ONE);
  } else {
    temp = 0;
  }
  writePackedBlock(m, row, coloffset, temp);

  temp=0;

  for ( i=startblock+1; i < (m->width); i++ ) {
    writePackedBlock(m, row, i*RADIX, temp);
  }
}


/**********************************************************************/
/* this adds rows sourcerow and destrow and stores the total in row
   destrow, but only begins at the column coloffset */

void rowAddPackedOffset( packedmatrix *m, int sourcerow, int destrow,
		   int coloffset ) {

  int startblock= coloffset/RADIX;
  int i;
  word temp;

  // make sure to start adding at coloffset
  temp=readPackedBlock(m, sourcerow, startblock*RADIX);
  if (coloffset%RADIX)
    temp &= (ONE<<(RADIX - (coloffset%RADIX))) - ONE;
  xorPackedBlock(m, destrow, startblock*RADIX, temp);

  for ( i=startblock+1; i < (m->width); i++ ) {
    temp=readPackedBlock(m, sourcerow, i*RADIX);
    xorPackedBlock(m, destrow, i*RADIX, temp);
  }
}


void rowAddPacked( packedmatrix *m, int sourcerow, int destrow) {
  rowAddPackedOffset(m, sourcerow, destrow, 0);
}

/**********************************************************************/
/* This will do Gaussian Elimination on the matrix m but will start not
 at column 0 necc but at column "startcol". If full=NO, then it will do
 triangular style elimination, and if full=YES, it will do Gauss-Jordan style,
 or full elimination.*/

int gaussianPackedDelayed(packedmatrix *m, int startcol, int full) {
  int i,j;
  int start;

  int startrow = startcol;
  int ii;
  int pivots = 0;
  for (i=startcol ; i<m->cols ; i++) {

    for(j=startrow ; j < m->rows; j++) {
      if (readPackedCell(m,j,i)) {
	rowSwapPacked(m,startrow,j);
	pivots++;

	if (full==YES) start=0; else start=i+1;

	for(ii=start ;  ii < m->rows ; ii++) {
	  if (ii != startrow) {
	    if (readPackedCell(m, ii, i)) {
	      rowAddPackedOffset(m, startrow, ii, i);
	    }
	  }
	}
	startrow = startrow + 1;
	break;
      }
    }
  }

  return pivots;
}


/* This will do Gaussian Elimination on the matrix m.
   If full=NO, then it will do
 triangular style elimination, and if full=YES, it will do Gauss-Jordan style,
 or full elimination.*/
int gaussianPacked(packedmatrix *m, int full) {
  return gaussianPackedDelayed(m,0, full);
}


/**********************************************************************/

inline bit dotProductPacked( word a, word b ) {
  word temp=a & b;
  //int i,
  int total=0;

/*   for (i=0; i<RADIX; i++) { */
/*     if ((temp & packingmask[i])!=0) total++; */
/*   } */
/*   return (total % 2); */
  while (temp)  {
    total = !total;
    temp = temp & (temp - 1);
  }
  return total;
}

/**********************************************************************/

/* Works */
/* MEMLEAK, use destroyPackedMatrix */
/* This is not efficient, but it is quadratic time, so who cares? */
/* Efficient, would be to use the fact that: */
/* [ A B ]T    [AT CT]
   [ C D ]  =  [BT DT] and thus rearrange the blocks recursively. */
packedmatrix *transposePacked( packedmatrix *data ) {
  packedmatrix *newmatrix=createPackedMatrix( data->cols, data->rows );
  int i,j,k;
  word temp;

  for (i=0; i<newmatrix->rows; i++) {
    for (j=0; j<newmatrix->width; j++) {
      temp=(word)0;
      for (k=0; k<RADIX; k++) {
	if (  (j*RADIX+k) < data->rows ) {
	  if (readPackedCell(data, j*RADIX+k, i)==1)
	    temp=temp | packingmask[k];
	}
      }
      writePackedBlock(newmatrix, i, j*RADIX, temp);
    }
  }

  return newmatrix;
}

/********************************************************/

/* Works */
/* Internal to naive matrix mult */
bit bigDotProductPacked( packedmatrix *a, packedmatrix *bT, int rowofa,
			 int rowofb ) {
  /* ``a slot'' is a row of A, and a column of B when calcing AB */
  /* but since we use B^T so that we are working only with rows, */
  /* ``a slot'' of A is a row, ``a slot'' of B is a row of B^T */
  int total, i;

  total=0;
  for (i=0; i< a->width; i++) {
    if (  (i*RADIX) < a->rows )
    total+=dotProductPacked( readPackedBlock(a, rowofa, i*RADIX),
		       readPackedBlock(bT, rowofb, i*RADIX) );
  }

  return (bit)(total % 2);
}

/********************************************************/

/* Works */
/* MEMLEAK use destroyPackedMatrix */

packedmatrix *matrixTimesMatrixTransposePacked( packedmatrix *a,
						packedmatrix *bT ) {
  int newrows=a->rows;
  int newcols=bT->rows;
  packedmatrix *newmatrix=createPackedMatrix(newrows, newcols);
  int i, j;

  for (i=0; i<newrows; i++) {
    for (j=0; j<newcols; j++) {
      writePackedCell(newmatrix, i, j, bigDotProductPacked( a, bT, i, j ) );
    }
  }

  return newmatrix;
}

/********************************************************/

/* Works */
/* MEMLEAK: use destroyPackedMatrix */
/* Normally, if you will multiply several times by b, it is smarter to
  calculate bT yourself, and keep it, and then use the function called
  matrixTimesMatrixTransposePacked */
packedmatrix *matrixTimesMatrixPacked( packedmatrix *a,
				 packedmatrix *b) {
  packedmatrix *bT=transposePacked(b);

  packedmatrix *product=matrixTimesMatrixTransposePacked( a, bT );

  destroyPackedMatrix(bT);

  return product;
}

/********************************************************/

/* Works */
word randomWord() {
  int i;
  word temp=0;

  for (i=0; i<RADIX; i++) {
    if (coinFlip()==1) {
      temp|=packingmask[i];
    }
  }

  return temp;
}

/********************************************************/

/* Works */
void fillRandomlyPacked( packedmatrix *a ) {
  int i, j;

  for (i=0; i < (a->rows); i++) {
    for (j=0; j < (a->cols); j++) {
      writePackedCell(a, i, j, coinFlip() );
    }
  }
}

/********************************************************/

/* Works */
void makeIdentityPacked( packedmatrix *a ) {
  int i,j;
  int stop=min(a->rows, a->cols);

  for (i=0; i< (a->rows); i++) {
    for (j=0; j< (a->width); j++) {

      writePackedBlock(a, i, j*RADIX, 0);
    }
  }

  for (i=0; i<stop; i++) {
    writePackedCell(a, i, i, 1);
  }
}

/********************************************************/

/* Works */
/* MEMLEAK: use destroyMatrix */
/* I'm not sure there is any reason to use matrix vice packedmatrix, but
   perhaps there is and we don't realize just yet */
matrix *unpackMatrix( packedmatrix *pm ) {
  matrix *newmatrix=createMatrix( pm->rows, pm->cols );
  int i, j;

  for (i=0; i<pm->rows; i++) {
    for (j=0; j<pm->cols; j++) {
      writeCell(newmatrix, i, j, readPackedCell( pm, i, j ));
    }
  }

  return newmatrix;
}

/********************************************************/

/* Works */
/* MEMLEAK: use destroyPackedMatrix */

packedmatrix *packMatrix( matrix *m ) {
  packedmatrix *newmatrix=createPackedMatrix( m->rows, m->cols );
  int i, j;

  for (i=0; i<m->rows; i++) {
    for (j=0; j<m->cols; j++) {
      writePackedCell(newmatrix, i, j, readCell( m, i, j ));
    }
  }

  return newmatrix;
}

/********************************************************/

/* Works */
inline bit isEqualWord( word a, word b) {
  if (a==b) return YES;
  else return NO;
}

/********************************************************/

/* Works */
bit equalPackedMatrix( packedmatrix *a, packedmatrix *b ) {
  int i, j;
  word block1, block2;

  if (a->rows!=b->rows) return NO;
  if (a->cols!=b->cols) return NO;

  for (i=0; i< a->rows; i++) {
    for (j=0; j< a->width; j++) {
      block1=readPackedBlock(a, i, j*RADIX);
      block2=readPackedBlock(b, i, j*RADIX);
      if (isEqualWord(block1, block2)==NO) return NO;
    }
  }
  return YES;
}

int comparePackedMatrix(packedmatrix *a, packedmatrix *b) {

  int i,j;

  if(a->rows < b->rows) return -1;
  if(b->rows < a->rows) return 1;
  if(a->cols < b->cols) return -1;
  if(b->cols < a->cols) return 1;

  for(i=0; i < a->rows ; i++) {
    for(j=0 ; j< a->width ; j++) {
      if ( a->values[a->rowswap[i] + j] < b->values[b->rowswap[i] + j])
	return -1;
      else if( a->values[a->rowswap[i] + j] > b->values[b->rowswap[i] + j])
	return 1;
    }
  }
  return 0;
}

/********************************************************/

/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *clonePacked( packedmatrix *p) {
  packedmatrix *newmatrix=createPackedMatrix(p->rows, p->cols);
  int i, j;
  word entry;

  for (i=0; i<p->rows; i++) {
    for (j=0; j<p->width; j++) {

      if (  (j*RADIX) < p->cols ) {
	entry=readPackedBlock(p, i, j*RADIX);
	writePackedBlock(newmatrix, i, j*RADIX, entry);
      }
    }
  }

  return newmatrix;
}


/********************************************************/

/* MEMLEAK: use destroyPackedMatrix */
/* This is sometimes called augment */
packedmatrix *concatPacked( packedmatrix *a, packedmatrix *b) {
  packedmatrix *newmatrix;
  int i, j;
  //word entry;

  if (a->rows!=b->rows) {
    die("Bad arguments to concatPacked!\n");
  }

  newmatrix=createPackedMatrix(a->rows, a->cols + b->cols);

  for (i=0; i<a->rows; i++) {
    for (j=0; j<a->cols; j++) {
      writePackedCell(newmatrix, i, j, readPackedCell(a, i, j) );
    }
  }

  for (i=0; i<b->rows; i++) {
    for (j=0; j<b->cols; j++) {
      writePackedCell(newmatrix, i, j+(a->cols),
		      readPackedCell(b, i, j) );
    }
  }

  return newmatrix;
}


/********************************************************/
/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *copySubMatrixPacked( packedmatrix *a, int lowr, int lowc,
				   int highr, int highc) {
  packedmatrix *newmatrix;
  int i, j;

  newmatrix=createPackedMatrix(highr-lowr+1, highc-lowc+1);

  for (i=lowr; i<=highr; i++) {
    for (j=lowc; j<=highc; j++) {
      writePackedCell(newmatrix, i-lowr, j-lowc, readPackedCell(a, i, j) );
    }
  }

  return newmatrix;
}

/*********************************************************/
/* MEMLEAK: use destroyPackedMatrix */
packedmatrix *invertPackedGaussian(packedmatrix *target,
				   packedmatrix *identity) {
  packedmatrix *huge, *inverse;
  int x;

  huge=concatPacked(target, identity);

  //startWatch();
  x=gaussianPacked(huge, YES);
  //stopWatch();

  if (x==NO) { destroyPackedMatrix(huge); return NULL; }

  inverse=copySubMatrixPacked(huge, 0, target->cols, target->rows-1,
			      target->cols*2-1);

  destroyPackedMatrix(huge);
  return inverse;
}

packedmatrix *addPacked(packedmatrix *left, packedmatrix *right) {

  packedmatrix *ret;
  int i,j,width;

  if (left->rows != right->rows || left->cols != right->cols) {
    die("rows and columns must match");
  }

  width = left->width;

  ret = createPackedMatrix(left->rows, left->cols);
  memcpy(ret->values, left->values, (RADIX>>3) * width * left->rows);
  memcpy(ret->rowswap, left->rowswap, left->rows * sizeof(int));

  for(i=0; i < right->rows; i++) {
    for(j=0; j < width; j++) {
      ret->values[  ret->rowswap[i] + j ]  ^= right->values[ right->rowswap[i] + j];
    }
  }
  return ret;
}

void lazyPrint(packedmatrix *a) {
  matrix *b=unpackMatrix(a);

  printMatrix(b);

  destroyMatrix(b);
}
