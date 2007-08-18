/*============================================================================
    Copyright 2006 William Hart

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

============================================================================*/
#ifndef LANCZOS_H
#define LANCZOS_H

#ifdef __sun
#define u_int32_t unsigned int
#define u_int64_t unsigned long long
#endif


#include <stdlib.h>

typedef struct {
	unsigned long *data;		/* The list of occupied rows in this column */
	unsigned long weight;		/* Number of nonzero entries in this column */
	unsigned long orig;         /* Original relation number */
} la_col_t;

u_int64_t getNullEntry(u_int64_t *, long, long);
void reduce_matrix(unsigned long *, unsigned long *, la_col_t *);
u_int64_t * block_lanczos(unsigned long, unsigned long, unsigned long, la_col_t*);

/*==========================================================================
   insertColEntry:

   Function: insert an entry into a column of the matrix,
   reallocating the space for the column if necessary

===========================================================================*/
static inline void insertColEntry(la_col_t* colarray, unsigned long colNum, unsigned long entry)
{
   unsigned long* temp;

       if ((((colarray[colNum].weight)>>4)<<4)==colarray[colNum].weight) //need more space
       {
           temp = colarray[colNum].data;
           colarray[colNum].data = (unsigned long*)malloc((colarray[colNum].weight+16)*sizeof(unsigned long));
           for (long i = 0; i<colarray[colNum].weight; i++)
           {
               colarray[colNum].data[i] = temp[i];
           }
           if (colarray[colNum].weight!=0) free(temp);
       }

   colarray[colNum].data[colarray[colNum].weight] = entry;
   colarray[colNum].weight++;
   colarray[colNum].orig = colNum;
}

/*==========================================================================
   xorColEntry:

   Function: xor entry corresponding to a prime dividing A, which will be
   either the last entry in the column, or not there at all, so we either
   add it in or take it out

===========================================================================*/
static inline void xorColEntry(la_col_t* colarray, unsigned long colNum, unsigned long entry)
{
   for (long i = 0; i < colarray[colNum].weight; i++)
     if (colarray[colNum].data[i] == entry)
     {
        for (unsigned long j = i; j < colarray[colNum].weight - 1; j++)
          colarray[colNum].data[j] = colarray[colNum].data[j+1];
        colarray[colNum].weight--;
        return;
     }
   insertColEntry(colarray,colNum,entry);
}

/*==========================================================================
   clearCol:

   Function: clear a column

===========================================================================*/
static inline void clearCol(la_col_t* colarray, unsigned long colNum)
{
   colarray[colNum].weight =0;
}

#endif
