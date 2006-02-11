#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int main(int argc, char **argv)
{char INCURV[128]; int i=1; int GPMEM=4000000; int STDIN=1; int CLEAN=0;
 int dorank=0; int domoddeg=0; int dox0isog=0; int doanalrank=0;
 int count=0; int FIXIT=0; int doarderivs=0; int CHECKIT=0;

 ELLACC=3; TRACE=0; VERBOSE=1; *INCURV=0; PRINT=1; DISK=0; READ=0; PCOUNT=0;
 INITANALRANK=0; plimit=1000000; DOTWIST=0; LOOKUP=0; HIRANK=0; MDFORCE=0;
 DOADD=0; SPIT=0; OFILE=0;

 while(i<argc)
 {if (!strcmp(argv[i],"-all"))
  {domoddeg=1; dorank=1; dox0isog=1; doanalrank=1; i++;}
  else if (!strcmp(argv[i],"-arderivs")) {doarderivs=1; i++;}
  else if (!strcmp(argv[i],"-analrank")) {doanalrank=1; i++;}
  else if (!strcmp(argv[i],"-checkit"))
  {CHECKIT=atoi(argv[i+1]); STDIN=0; i+=2;}
  else if (!strcmp(argv[i],"-clean")) {CLEAN=1; STDIN=0; i++;}
  else if (!strcmp(argv[i],"-count")) {count=1; STDIN=0; i++;}
  else if (!strcmp(argv[i],"-curve"))
  {strcpy(INCURV,argv[i+1]); STDIN=0; i+=2;}
  else if (!strcmp(argv[i],"-doadd")) {DOADD=1; i++;}
  else if (!strcmp(argv[i],"-dodisk")) {DISK=1; i++;}
  else if (!strcmp(argv[i],"-fixit"))
  {FIXIT=1; OUTFILE=fopen(argv[i+1],"r"); STDIN=0; PRINT=0; i+=2;
   printf("FILE %s\n",argv[i-1]);
  }
  else if (!strcmp(argv[i],"-gpmem")) {GPMEM=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-hirank")) {HIRANK=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-lookup")) {LOOKUP=1; i++;}
  else if (!strcmp(argv[i],"-mdforce")) {MDFORCE=1; i++;}
  else if (!strcmp(argv[i],"-moddeg")) {domoddeg=1; dox0isog=1; i++;}
  else if (!strcmp(argv[i],"-ofile"))
  {OUTFILE=fopen(argv[i+1],"a"); OFILE=1; i+=2;}
  else if (!strcmp(argv[i],"-pcount")) {PCOUNT=1; i++;}
  else if (!strcmp(argv[i],"-plimit")) {plimit=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-primelimit")) {plimit=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-printoff")) {PRINT=0; i++;}
  else if (!strcmp(argv[i],"-quiet")) {VERBOSE=0; i++;}
  else if (!strcmp(argv[i],"-rank")) {dorank=1; i++;}
  else if (!strcmp(argv[i],"-read")) {READ=1; i++;}
  else if (!strcmp(argv[i],"-spit")) {SPIT=1; i++;}
  else if (!strcmp(argv[i],"-stdin")) {STDIN=1; i++;}
  else if (!strcmp(argv[i],"-trace")) {TRACE=1; i++;}
  else if (!strcmp(argv[i],"-twist")) {DOTWIST=1; i++;}
  else if (!strcmp(argv[i],"-verbose")) {VERBOSE=2; i++;}
  else if (!strcmp(argv[i],"-x0isog")) {dox0isog=1; i++;}
  else if (!strcmp(argv[i],"-x0opt")) {dox0isog=1; i++;}
 }

 pari_init(GPMEM,plimit); PMAX=plimit; PSIZE=plimit;
 if ((!DISK) && (!READ) && (PRINT))
   printf("ec version 1.12b  %s\n",PARIVERSION);

 if (*INCURV!=0)
   docurve(INCURV,dorank,dox0isog,domoddeg,doanalrank,doarderivs);
 if (DISK) {dodisk(); STDIN=0;} if (READ) {doread(); STDIN=0;}
 if (MDFORCE) {mdforce(); STDIN=0;} if (count) docount();
 if (FIXIT) {fixit(); printf("FIXED %s\n",argv[i-1]);}
 if (CHECKIT) checkit(CHECKIT);
 while(STDIN)
 {scanf("%s",INCURV); if (!(strcmp(INCURV,"0"))) break;
  docurve(INCURV,dorank,dox0isog,domoddeg,doanalrank,doarderivs);
 }

 return(100);
}

