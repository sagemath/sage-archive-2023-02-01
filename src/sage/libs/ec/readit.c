#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int FLAG=0;

void spitout(GEN C)
{GEN D;

 D=gabs((GEN) C[12],-1); if (!isprime(D)) return;
 if ((mpcmp(D,stoi(100000000))==-1) &&
     (mpcmp(gabs((GEN) C[5],-1),gmul2n(gun,31))==1)) output(C);
 return;
}

void readit(int i)
{int j,mds; int64 moddeg;

 X0NUM=(((int) leadbyte[i])>>6)+1;
 mds=(leadbyte[i]>>3)&7; moddeg=0; if (mds>4) FLAG=1; if (!mds) FLAG=0;
 for (j=0;j<mds;j++) {moddeg=moddeg<<8; moddeg+=(int) moddegbytes[8*i+j];}
 if (FLAG) degphi((GEN) CL[X0NUM],gun,1); else MODDEG=lltoi(moddeg);
 if (moddeg)
 {if (ggval(MODDEG,gdeux)+1<itos((GEN) matsize(CF)[1]))
   degphi((GEN) CL[X0NUM],gun,1);
 }
 if (mds<=4) FLAG=0;
}

void doread()
{char filename[32]; int i,j; GEN ib,C4,C6; pari_sp memspot=avma;

 PRINT=0; if (VERBOSE==1) VERBOSE=0; TWPROD=gun; ib=stoi(100000000);
 while (1)
 {scanf("%s",filename); if (!strcmp(filename,"END")) break; readfile(filename);
  for (i=0;i<numcurves;i++)
  {memspot=avma; C4=gzero; C6=gzero;
   for (j=0;j<c4s;j++)
   {C4=gadd(gmul2n(C4,8),stoi((int) curvedata[i*recsize+j]));}
   for (j=0;j<c6s;j++)
   {C6=gadd(gmul2n(C6,8),stoi((int) curvedata[i*recsize+c4s+j]));}
   if (c4neg) C4=gneg(C4); if (c6neg) C6=gneg(C6);
   C4=gadd(stoi(c4r),gmul(C4,stoi(576)));
   C6=gadd(stoi(c6r),gmul(C6,stoi(1728)));
   CURVE=avecfromc4c6(C4,C6);
   ROOTNO=-1; TAMA=gun; COND=gun; ISOG=0; PLACE=0; CURVES=gzero;
   ISPRIME=0; ISSQFREE=0; ISSETZER=0; MODDEG=gzero; CM=0; AISOG=0; NI=0;
   CURVE=minimalequation(); minimaltwist(); symsq(); isogenies();
   AISOG=absval(ISOG); NI=classsize();
   for (j=1;j<=NI;j++) CL[j]=(long) ellinit0((GEN) CL[j],1,ELLACC);
   ISSETZER=KnownDiffOptimal((GEN) CL[1]); X1VOL=myvol((GEN) CL[1]);
   if (anarray) free(anarray); if (antwarray) free(antwarray);
   antwsize=1; PREVD=0; ansize=1;
   anarray=malloc(u8((1+ansize)*sizeof(double))); anarray[1]=1.0;
   antwarray=malloc(u8((ansize+1)*sizeof(double))); antwarray[1]=1.0;
   readit(i); if (gequal(TWPROD,gun) && (mpcmp(ib,TWCOND)==1)) j=whichtwists();
   avma=memspot;
  }
  free(curvedata);
  free(size16bytes); free(size256bytes); free(size4096bytes);
  free(leadbyte); free(moddegbytes); free(localsizebytes);
  for (i=0;i<numcurves;i++) free(twdata[i]); free(twdata);
 }
}

