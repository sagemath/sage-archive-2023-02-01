/* allfaces.c: Test program to call the cdd library cddlib
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.94, August 4, 2005
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
*/

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "setoper.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

dd_boolean SetInputFile(FILE **f, dd_DataFileType fname)
{
  dd_boolean success=dd_FALSE;
  success=dd_FALSE;

  if ( ( *f = fopen(fname, "r") )!= NULL) {
    printf("input file %s is open\n", fname);
    success=dd_TRUE;
  }
  else{
    printf("The input file %s not found\n",fname);
  }
  return success;
}

dd_boolean SetWriteFile(FILE **f, dd_DataFileType fname)
{
  dd_boolean success=dd_FALSE;

  if ( (*f = fopen(fname, "w")) != NULL){
    printf("output file %s is open\n",fname);
    success=dd_TRUE;
  }
  else{
    printf("The output file %s cannot be opened\n",fname);
  }
  return success;
}

dd_boolean FaceEnum(dd_MatrixPtr M, dd_rowset R, dd_rowset S, dd_boolean rip, dd_colrange mindim)
{
  dd_ErrorType err;
  dd_rowset LLL, ImL, RR, SSS, Lbasis;
  dd_rowrange i,iprev=0;
  dd_colrange j,dim;
  dd_LPSolutionPtr lps=NULL;
  dd_boolean success=dd_FALSE;

  set_initialize(&LLL, M->rowsize);
  set_initialize(&RR, M->rowsize);
  set_initialize(&SSS, M->rowsize);
  set_copy(LLL, M->linset); /* rememer the linset. */
  set_copy(RR, R); /* copy of R. */
  set_copy(SSS, S); /* copy of S. */
  if (dd_ExistsRestrictedFace(M, R, S, &err)){
    set_uni(M->linset, M->linset, R);
	dd_FindRelativeInterior(M, &ImL, &Lbasis, &lps, &err);
	dim=M->colsize - set_card(Lbasis)-1;
    set_uni(M->linset, M->linset, ImL);
	fprintf(stdout,"%ld: ", dim); set_fwrite(stdout,M->linset);
	if (rip){
    	/* Write an interior point. */
	  printf("RIP: (");
	  for (j=1; j <(lps->d)-1; j++) {
	    dd_WriteNumber(stdout,lps->sol[j]);
	  }
	  printf(")\n");
	}
	dd_FreeLPSolution(lps);
	set_free(ImL);
	set_free(Lbasis);

    if (dim>mindim){
	  for (i=1; i<=M->rowsize; i++){
	    if (!set_member(i, M->linset) && !set_member(i, S)){
		  set_addelem(RR, i);
		  if (iprev) {
		    set_delelem(RR,iprev);
		    set_delelem(M->linset,iprev);
		    set_addelem(SSS, iprev);
		  }
		  iprev=i;
		  FaceEnum(M, RR, SSS, rip, mindim);
		}
	  }
	}
  } else if (err!=dd_NoError) goto _L99;
  success=dd_TRUE;

_L99:
  set_copy(M->linset, LLL); /* restore the linset */
  set_free(LLL);
  set_free(RR);
  set_free(SSS);
  return success;
}



int main(int argc, char *argv[])
{
  dd_MatrixPtr M=NULL;
  dd_rowrange m;
  dd_ErrorType err=dd_NoError;
  dd_rowset R, S;
  dd_DataFileType inputfile;
  FILE *reading=NULL;
  char ch;
  dd_colrange mindim;
  dd_boolean rip=dd_FALSE;


  dd_set_global_constants();  /* First, this must be called. */

  if (argc>1) strcpy(inputfile,argv[1]);
  if (argc<=1 || !SetInputFile(&reading,argv[1])){
    dd_WriteProgramDescription(stdout);
    fprintf(stdout,"\ncddlib test program to list all faces of an H-polyhedron.\n");
    dd_SetInputFile(&reading,inputfile, &err);
  }
  if (err==dd_NoError) {
    M=dd_PolyFile2Matrix(reading, &err);
  }
  else {
    fprintf(stderr,"Input file not found\n");
    goto _L99;
  }

  if (err!=dd_NoError) goto _L99;

  if (M->representation==dd_Generator){
    printf("The input is V-representation.  Consider it as H-representation (N)? ");
    ch=getchar(); getchar();
    if (ch!='y' && ch!='Y') goto _L99;
  }

  m=M->rowsize;

  set_initialize(&R, M->rowsize);
  set_initialize(&S, M->rowsize);

  printf("Output relative interior points (N)? ");
  ch=getchar();
  if (ch=='y' || ch=='Y') rip=dd_TRUE;
  printf("Minimum dimension of faces to list (0..%ld) ? ",M->colsize-1);
  scanf("%ld", &mindim);
  if (mindim>=M->colsize) mindim=M->colsize-1;
  printf("Minimum dimension is set to %ld.", mindim);

  printf("\n--- FaceEnum (dim: active set) ---\nbegin\n");
  FaceEnum(M, R, S, rip, mindim);
  fprintf(stderr,"end\nFaceEnum completed.\n");

  dd_FreeMatrix(M);
  set_free(R);
  set_free(S);
_L99:;
  if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err);
  dd_free_global_constants();  /* At the end, this should be called. */
  return 0;
}


/* end of allfaces.c */
