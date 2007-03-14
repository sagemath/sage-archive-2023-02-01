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

#include "grayflex.h"

code **codebook;
const int MAXKAY = 16;

void printBitString(int number, int length){
  int i;
  for(i=length-1 ; i>=0; i--) {
    ((1<<i) & number) ? printf("1") : printf("0");
  }
  printf("\n");
}

int swap_bits(int v,int length) {
 unsigned int r = v; // r will be reversed bits of v; first get LSB of v
 int s = length - 1; // extra shift needed at end

 for (v >>= 1; v; v >>= 1)
    {
      r <<= 1;
      r |= v & 1;
      s--;
    }
 r <<= s;
 return r;
}

int grayCode(int number, int length) {
  int lastbit = 0;
  int res = 0;
  int i,bit;
  for(i=length-1; i>=0;  i--) {
    bit = number & (1<<i);
    res |= ((lastbit>>1) ^ bit);
    lastbit = bit;
  };
  return swap_bits(res,length) & ((1<<length)-1);
  //return res;
}

void buildCodeFlex(int *ord, int *inc, int length) {
  int i,j;

  // this one is easy.
  for(i=0 ; i < TWOPOW(length) ; i++) {
    ord[i] = grayCode(i,length);
  }

  for(i = length ; i>0 ; i--) {
    for(j=1 ; j < TWOPOW(i) + 1 ; j++) {
      inc[j *TWOPOW(length-i) -1 ] = length - i;
    }
  }
}

void buildAllCodes() {
  int k;
  codebook=calloc(MAXKAY+1, sizeof(code *));

  for(k=1 ; k<MAXKAY+1; k++) {
    codebook[k] = (code *)calloc(sizeof(code),1);
    codebook[k]->ord =(int *)calloc(TWOPOW(k),sizeof(int));
    codebook[k]->inc =(int *)calloc(TWOPOW(k),sizeof(int));
    buildCodeFlex(codebook[k]->ord, codebook[k]->inc, k);
  }
}

void destroyAllCodes() {
  int i;
  for(i=1; i<MAXKAY+1; i++) {
    free(codebook[i]->inc);
    free(codebook[i]->ord);
    free(codebook[i]);
  }
  free(codebook);
}
