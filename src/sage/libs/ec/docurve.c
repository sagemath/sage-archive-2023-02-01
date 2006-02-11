#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void docurve
(char *INP,int dorank,int dox0isog,int domoddeg,int doanalrank,int doarderivs)
{unsigned char *TWFERRY; int fersize,b,j; GEN P,TA,T,P1,TA1,T1,Q,PV,E,p,l;
 GEN X0TWCURVE,X0TWEFF,X0CURVE; int i=1; int ar=-1; int terseoutput=1;
 double LRAT=0; pari_sp topmem=avma; GEN XMODDEG=gzero; int n;

 ELLACC=6; TRACE=0; VERBOSE=1; PRINT=1; DISK=0; READ=0; PCOUNT=0;
 INITANALRANK=0; plimit=1000000; DOTWIST=0; LOOKUP=0; HIRANK=0; MDFORCE=0;
 DOADD=0; SPIT=0; OFILE=0;

 if (!precdl)
   pari_init(4000000, plimit);

 PMAX=plimit; PSIZE=plimit;

 if (INP[0]=='[') CURVE=getcurve(INP); else CURVE=specialcurve(INP);
 if (gcmp0(CURVE)) {printf("Unknown entry\n"); return;}
 if (anarray) free(anarray); if (antwarray) free(antwarray);
 P1=gzero; TA1=gzero; T1=gzero; antwsize=1; PREVD=0; ansize=1;
 anarray=malloc(u8((1+ansize)*sizeof(double))); anarray[1]=1.0;
 antwarray=malloc(u8((ansize+1)*sizeof(double))); antwarray[1]=1.0;

 ROOTNO=-1; TAMA=gun; COND=gun; ISOG=0; PLACE=0; CURVES=gzero;
 ISPRIME=0; ISSQFREE=0; ISSETZER=0; MODDEG=gzero; CM=0; AISOG=0; NI=0;
 CURVE=minimalequation();
 if (PRINT)
 {printf("Minimal Equation is "); printcurven(CURVE);
  printf("of conductor "); output((GEN) (COND));
 }
 minimaltwist();
 symsq();
 isogenies();
 AISOG=absval(ISOG);
 NI=classsize();
 for (i=1;i<=NI;i++)
   CL[i] = (long) ellinit0((GEN) CL[i],1,ELLACC);
 sortcurves();

 if ((VERBOSE) && (CM!=0)) printf("Complex multiplication by %i\n",CM);
 if (PRINT)
 {printf("ISOG:%i PLACE:%i",ISOG,PLACE); n=itos((GEN) matsize(CF)[1]);
  if (VERBOSE>1)
  {printf(" ["); for (j=1;j<=n;j++)
   {p=(GEN) ((GEN) CF[1])[j]; outbrute(p); if (j!=n) printf(",");} printf("]");
  }
  printf("\n");
  for (i=1;i<=NI;i++)
  {E=gcopy((GEN) CL[i]); printcurven(E);
   if (VERBOSE>1)
   {printf(" "); if (gsigne((GEN) E[12])==1) printf("["); else printf("(");
    for (j=1;j<=n;j++)
    {p=(GEN) ((GEN) CF[1])[j]; printf("%li",ggval((GEN) E[12],p));
     l=elllocalred(E,p); if (j!=n) printf(",");
     else {if (gsigne((GEN) E[12])==1) printf("] "); else printf(") ");}
    }
    P=(GEN) periodvolvec(E)[2];
    Q=cgetg(20,t_VEC); for (j=1;j<=13;j++) Q[j]=E[j];
    b=4+((lgefint((GEN) Q[12])-2)>>1); PV=periodvolvec0(Q,b+2);
    Q[15]=PV[2]; Q[16]=PV[3]; Q[14]=(long) gzero; Q[17]=(long) gzero;
    Q[18]=(long) gzero; Q[19]=(long) gzero;
    if (gsigne((GEN) E[12])==1) P=gmul2n(P,1);
    TA=gettama(E); T=(GEN) elltors0(Q,0)[1];
    if (i==1) {printf("1 1 1 1\n"); P1=gcopy(P); TA1=gcopy(TA); T1=gcopy(T);}
    else
    {outbrute(volratio(P1,P)); printf(" "); outbrute(gdiv(TA1,TA));
     printf(" "); outbrute(gdiv(T1,T)); printf(" ");
     output(gdiv(gdiv(gsqr(gdiv(T,T1)),gdiv(TA,TA1)),volratio(P,P1)));
    }
   }
   else printf("\n");
  }
 }

printf("1\n");

 if (NI!=1) ISSETZER=KnownDiffOptimal((GEN) CL[1]);
 else X1VOL=myvol((GEN) CL[1]);
 if ((ISOG!=1) && (dox0isog))
   {
     X0NUM=x0isogeny();
     X0CURVE=ellinit0((GEN) CURVES[X0NUM],1,ELLACC);
     if ((X0NUM!=1) && (PRINT))
       {printf("Optimal X0-curve is "); printcurve(X0CURVE);
	if (!KnownDiffOptimal((GEN) CURVES[1])) printf("This curve is unknown.\n");
      }
 }
 else X0CURVE=ellinit0((GEN) CURVES[1],1,ELLACC);

 if (ISOG==1) X0NUM=1;

printf("2\n");

 if (doanalrank)
 {ar=analrank(qtwist(X0CURVE,TWPROD),&LRAT);
  if (PRINT)
  {printf("Analytic rank of "); printcurven(CURVE);
   printf("appears to be %i with L-ratio of %f\n",ar,LRAT);
  }
 }

printf("3\n");

 if (doarderivs) {firstderivs();}

printf("4\n");

 if (domoddeg)
 {
   X0TWCURVE = retminimaltwist(X0CURVE,&X0TWEFF);

   degphi(X0TWCURVE,X0TWEFF,1);

   printf("Modular Degree of ");
     printcurven(X0CURVE);
   printf("is ");
     output(MODDEG);

   printf("X0NUM is %d\n",X0NUM);
   i = isogdeg(X0NUM); XMODDEG=gmul(stoi(i),MODDEG);
   printf("Modular degree of "); printcurven(CURVE);
   printf("is "); output(XMODDEG);
 }

printf("5\n");

 if (DOTWIST)
 {TWFERRY=malloc(0); pile=malloc(0); fersize=whichtwists();
  free(TWFERRY); free(pile);
 }



/* if ((terseoutput) && (PRINT)) */
if (1)
 {printf("CURV "); printcurven(CURVE); printf(" "); outbrute(COND);
  printf(" "); outbrute(TWPROD); printf(" "); outbrute(EXEFF);
  printf(" %i",ISOG);
  if (dox0isog)
  {if (X0NUM==1) printf(" 1");
   else if (KnownDiffOptimal((GEN) CURVES[1])) printf(" %i",X0NUM);
   else printf(" %iX",X0NUM);
  }
  else
    printf(" 0");
  if (doanalrank) printf(" %i %f ",ar,LRAT); else printf(" X X ");
  if (domoddeg) output(XMODDEG); else printf("0\n");
 }

 avma=topmem;
}

#define append(key,val) ans = concat(ans, strtoGENstr(key)); ans = concat(ans, sep); \
                        ans = concat(ans, val); ans = concat(ans, sep);

GEN sage_docurve
(char *INP,int dorank,int dox0isog,int domoddeg,int doanalrank,int doarderivs)
{unsigned char *TWFERRY; int fersize,b,j; GEN P,TA,T,P1,TA1,T1,Q,PV,E,p,l;
 GEN X0TWCURVE,X0TWEFF,X0CURVE; int i=1; int ar=-1; int terseoutput=1;
 double LRAT=0; pari_sp topmem=avma; GEN XMODDEG=gzero; int n;
 GEN ans, sep;
PRINT=0;
VERBOSE=0;
 ans = strtoGENstr("");
 sep = strtoGENstr("|");

 ELLACC=6; TRACE=0; VERBOSE=1; PRINT=1; DISK=0; READ=0; PCOUNT=0;
 INITANALRANK=0; plimit=1000000; DOTWIST=0; LOOKUP=0; HIRANK=0; MDFORCE=0;
 DOADD=0; SPIT=0; OFILE=0;

 if (!precdl)
   pari_init(4000000, plimit);

 PMAX=plimit; PSIZE=plimit;

 if (INP[0]=='[') CURVE=getcurve(INP); else CURVE=specialcurve(INP);
 if (gcmp0(CURVE)) {printf("Unknown entry\n"); return;}
 if (anarray) free(anarray); if (antwarray) free(antwarray);
 P1=gzero; TA1=gzero; T1=gzero; antwsize=1; PREVD=0; ansize=1;
 anarray=malloc(u8((1+ansize)*sizeof(double))); anarray[1]=1.0;
 antwarray=malloc(u8((ansize+1)*sizeof(double))); antwarray[1]=1.0;

 ROOTNO=-1; TAMA=gun; COND=gun; ISOG=0; PLACE=0; CURVES=gzero;
 ISPRIME=0; ISSQFREE=0; ISSETZER=0; MODDEG=gzero; CM=0; AISOG=0; NI=0;
 CURVE=minimalequation();
 if (PRINT)
 {/* printf("Minimal Equation is "); printcurven(CURVE); */
   append("minimal equation", CURVE)
  /* printf("of conductor "); output((GEN) (COND)); */
   append("conductor",(GEN)COND)
 }
 minimaltwist();
 symsq();
 isogenies();
 AISOG=absval(ISOG);
 NI=classsize();
 for (i=1;i<=NI;i++)
   CL[i] = (long) ellinit0((GEN) CL[i],1,ELLACC);
 sortcurves();
 for (i=1;i<=NI;i++) {
   ans = concat(ans, strtoGENstr("Curve"));
   ans = concat(ans, stoi(i));
   ans = concat(ans, sep);
   ans = concat(ans, CL[i]);
   ans = concat(ans, sep);
 }

 if (NI!=1) ISSETZER=KnownDiffOptimal((GEN) CL[1]);
 else X1VOL=myvol((GEN) CL[1]);
 if ((ISOG!=1) && (dox0isog))
   {
     X0NUM=x0isogeny();
     X0CURVE=ellinit0((GEN) CURVES[X0NUM],1,ELLACC);
     if ((X0NUM!=1) && (PRINT))
       {/* printf("Optimal X0-curve is "); printcurve(X0CURVE); */
	append("Optimal X0-curve",X0CURVE)
	/* if (!KnownDiffOptimal((GEN) CURVES[1]))
	  printf("This curve is unknown.\n"); */
      }
 }
 else X0CURVE=ellinit0((GEN) CURVES[1],1,ELLACC);

 if (ISOG==1) X0NUM=1;

 if (doanalrank)
 {ar=analrank(qtwist(X0CURVE,TWPROD),&LRAT);
     append("analytic rank", stoi(ar));
    append("lratio", dbltor(LRAT));
 }

 if (doarderivs) {firstderivs();}

 if (domoddeg)
 {
   X0TWCURVE = retminimaltwist(X0CURVE,&X0TWEFF);

   degphi(X0TWCURVE,X0TWEFF,1);

   append("X0CURVE", X0CURVE);
   append("X0CURVE moddeg", MODDEG);

   append("X0NUM", stoi(X0NUM));

   i = isogdeg(X0NUM); XMODDEG=gmul(stoi(i),MODDEG);
   append("Modular degree",XMODDEG);
 }

 if (DOTWIST)
 {
   TWFERRY=malloc(0); pile=malloc(0); fersize=whichtwists();
   free(TWFERRY); free(pile);
 }

 avma=topmem;

 return ans;
}
