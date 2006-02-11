#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void dodisk()
{char filename[32]; int i,j; int64 modulardeg,ss; int mds,number;
 GEN C4,C6,p,l,X0CURVE,MDBOUND; int which,fersize,s,numdone;
 pari_sp memspot=avma;

 PRINT=0; VERBOSE=0; numdone=0; MDBOUND=stoi(500000); TWPROD=gun;
 PMAX=300000; PSIZE=300000;
 while (1)
 {scanf("%s",filename); if (!strcmp(filename,"END")) break; readfile(filename);
  printf("Opened %s ... %i curves here.\n",filename,numcurves);
  for (i=0;i<numcurves;i++)
  {if (!(leadbyte[i]&1))
   {if (i!=0) printf("Picking up at curve %i\n",i); i=numcurves;}}
  if (i==numcurves) printf("File already done!\n"); fflush(stdout);
  for (i=0;i<numcurves;i++)
  {memspot=avma; if (!(leadbyte[i]&1))
   {C4=gzero; C6=gzero;
    for (j=0;j<c4s;j++)
    {C4=gadd(gmul2n(C4,8),stoi((int) curvedata[i*recsize+j]));}
    for (j=0;j<c6s;j++)
    {C6=gadd(gmul2n(C6,8),stoi((int) curvedata[i*recsize+c4s+j]));}
    if (c4neg) C4=gneg(C4); if (c6neg) C6=gneg(C6);
    C4=gadd(stoi(c4r),gmul(C4,stoi(576)));
    C6=gadd(stoi(c6r),gmul(C6,stoi(1728)));
    CURVE=avecfromc4c6(C4,C6); TWCURVE=gcopy(CURVE);
    antwsize=1; PREVD=0; ansize=1;
    anarray=malloc(u8((1+ansize)*sizeof(double))); anarray[1]=1.0;
    antwarray=malloc(u8((ansize+1)*sizeof(double))); antwarray[1]=1.0;
    //printcurve(CURVE); fflush(stdout);
    ISSQFREE=1; ISPRIME=0; ISSETZER=0;
    ROOTNO=-1; TAMA=gun; COND=gun; TWCOND=gun; ISOG=0; PLACE=0; CURVES=gzero;
    CF=factor(gabs((GEN) CURVE[12],-1)); number=itos((GEN) matsize(CF)[1]);
    for(j=1;j<=number;j++)
    {p=(GEN) ((GEN) CF[1])[j];
     ROOTNO*=ellrootno(CURVE,p); l=elllocalred(CURVE,p);
     COND=gmul(COND,gpow(p,(GEN) l[1],-1)); TAMA=gmul(TAMA,(GEN) l[4]);
     ((GEN) (CF)[2])[j]=lcopy((GEN) l[1]);
     if (mpcmp((GEN) l[1],gun)==1) ISSQFREE=0;
    }
    if ((number==1) && (gequal((GEN) ((GEN) CF[2])[1],gun))) ISPRIME=1;
    TWCF=gcopy(CF); TWCOND=gcopy(COND); symsq(); isogenies();
    if (ISOG!=1)
    {if (ISPRIME)
     {which=1;
      if (carreparfait(gsub(COND,gmul2n(gun,6)))) which=2;
      if (gequal(COND,stoi(11))) which=2; if (gequal(COND,stoi(17))) which=3;
      if (gequal(COND,stoi(19))) which=2; if (gequal(COND,stoi(37))) which=2;
      X0CURVE=ellinit0((GEN) CURVES[which],1,ELLACC);
     }
     else {which=x0isogeny(); X0CURVE=ellinit0((GEN) CURVES[which],1,ELLACC);}
    }
    else {X0CURVE=ellinit0((GEN) CURVES[1],1,ELLACC); which=1;}
    TWFERRY=malloc(0); pile=malloc(0);
    fersize=whichtwists(); twdata[i]=malloc(u8(fersize));
    for (s=0;s<fersize;s++) twdata[i][s]=TWFERRY[s];
    free(TWFERRY); free(pile); pile=0; free(antwarray);
    if (mpcmp(SSCOND,MDBOUND)!=1)
    {degphi(X0CURVE,gun,0); modulardeg=gint64(MODDEG);}
    else modulardeg=0; free(anarray);

    mds=0; for (ss=1;ss<=modulardeg;ss*=256) mds++;
    leadbyte[i]=((which-1)<<6)+(mds<<3)+1;
    // WHICHX0 2 bits, MDsize 3bits, pointsearch? 1bit, done? 1bit
    if (fersize<240) localsizebytes[4*i]=fersize;
    else if (fersize<240*16)
    {localsizebytes[4*i]=240+(fersize>>8); localsizebytes[4*i+1]=fersize&255;}
    else if (fersize<240*256)
    {localsizebytes[4*i]=255; localsizebytes[4*i+1]=fersize>>8;
     localsizebytes[4*i+2]=fersize&255;
    }
    else
    {localsizebytes[4*i]=255; localsizebytes[4*i+1]=240+(fersize>>16);
     localsizebytes[4*i+2]=(fersize>>8)&255; localsizebytes[4*i+3]=fersize&255;
    }
    // record length 8bits, start hexF ->12,16,20 bits
    for (s=0;s<mds;s++) moddegbytes[8*i+s]=(modulardeg>>(8*(mds-1-s)))&255;
    // MD if extant, then twdata (already done)
    numdone++;
   }
   if ((numdone==10) || ((i+1)==numcurves)) {writefile(filename); numdone=0;}
   avma=memspot;
  }
  free(curvedata);
  free(size16bytes); free(size256bytes); free(size4096bytes);
  free(leadbyte); free(moddegbytes); free(localsizebytes);
  for (i=0;i<numcurves;i++) free(twdata[i]); free(twdata);
 }
 printf("Nothing more to do!\n");
}

void mdforce()
{char filename[32]; int i,j; int64 modulardeg,ss; int mds,number;
 GEN C4,C6,p,l,X0CURVE; int which,s,numdone,r,domd;
 pari_sp memspot=avma; double L;

 PRINT=0; VERBOSE=0; numdone=0; TWPROD=gun; ansize=0;
 while (1)
 {scanf("%s",filename); if (!strcmp(filename,"END")) break; readfile(filename);
  printf("Opened %s ... %i curves here.\n",filename,numcurves);
  fflush(stdout);
  for (i=0;i<numcurves;i++)
  {memspot=avma; C4=gzero; C6=gzero;
   for (j=0;j<c4s;j++)
   {C4=gadd(gmul2n(C4,8),stoi((int) curvedata[i*recsize+j]));}
   for (j=0;j<c6s;j++)
   {C6=gadd(gmul2n(C6,8),stoi((int) curvedata[i*recsize+c4s+j]));}
   if (c4neg) C4=gneg(C4); if (c6neg) C6=gneg(C6);
   C4=gadd(stoi(c4r),gmul(C4,stoi(576)));
   C6=gadd(stoi(c6r),gmul(C6,stoi(1728)));
   CURVE=avecfromc4c6(C4,C6); TWCURVE=gcopy(CURVE);
   ISSQFREE=1; ISPRIME=0; ISSETZER=0;
   ROOTNO=-1; TAMA=gun; COND=gun; TWCOND=gun; ISOG=0; PLACE=0; CURVES=gzero;
   CF=factor(gabs((GEN) CURVE[12],-1)); number=itos((GEN) matsize(CF)[1]);
   for(j=1;j<=number;j++)
   {p=(GEN) ((GEN) CF[1])[j];
    ROOTNO*=ellrootno(CURVE,p); l=elllocalred(CURVE,p);
    COND=gmul(COND,gpow(p,(GEN) l[1],-1)); TAMA=gmul(TAMA,(GEN) l[4]);
    ((GEN) (CF)[2])[j]=lcopy((GEN) l[1]);
    if (mpcmp((GEN) l[1],gun)==1) ISSQFREE=0;
   }
   TWCF=gcopy(CF); TWCOND=gcopy(COND); symsq(); isogenies();
   antwsize=1; PREVD=0; ansize=1;
   anarray=malloc(u8((1+ansize)*sizeof(double))); anarray[1]=1.0;
   antwarray=malloc(u8((ansize+1)*sizeof(double))); antwarray[1]=1.0;
   which=(((int) leadbyte[i])>>6)+1;
   X0CURVE=ellinit0((GEN) CURVES[which],1,ELLACC);
   if (rootnotwist(1)==1) r=0; else r=1;
   qtunhandle(twdata[i][0],twdata[i][1],&r,&L); domd=0;
   ss=gint64(SSCOND); if (ss<=1000000) domd=1;
   if ((r>=4) && (ss<=10000000)) domd=1;
   if ((r>=5) && (ss<=100000000)) domd=1;
   mds=(leadbyte[i]>>3)&7; if (mds!=0) domd=0;
   if (domd)
   {printcurven(X0CURVE); outbrute(SSCOND); printf(" %i\n",r); fflush(stdout);
    degphi(X0CURVE,gun,1); modulardeg=gint64(MODDEG);
    mds=0; for (ss=1;ss<=modulardeg;ss*=256) mds++;
    leadbyte[i]=((which-1)<<6)+(mds<<3)+1;
    // WHICHX0 2 bits, MDsize 3bits, pointsearch? 1bit, done? 1bit
    for (s=0;s<mds;s++) moddegbytes[8*i+s]=(modulardeg>>(8*(mds-1-s)))&255;
    // MD if extant, then twdata (already done)
    writefile(filename);
   }
   free(anarray); free(antwarray); avma=memspot;
  }
  free(curvedata);
  free(size16bytes); free(size256bytes); free(size4096bytes);
  free(leadbyte); free(moddegbytes); free(localsizebytes);
  for (i=0;i<numcurves;i++) free(twdata[i]); free(twdata);
 }
 printf("Nothing more to do!\n");
}

