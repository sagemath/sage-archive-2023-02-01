#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

GEN avecfromc4c6(GEN c4,GEN c6)
{GEN a1,a2,a3,a4,a6,b2,AI,BI,DI,thevec,z; int T=0; int g;

 a1=gmod(c4,gdeux); a2=gmod(gneg(gadd(a1,c6)),stoi(3));
 if (gequal(a2,gdeux)) a2=gneg(gun);
 T=0; DI=gcd0(gdiv(gsub(gmul(c4,gsqr(c4)),gsqr(c6)),stoi(1728)),gsqr(c6),0);
 g=ggval(DI,gdeux)/12;
 AI=gmul2n(c4,-g*4); BI=gmul2n(c6,-g*6);
 if ((gequal(gmod(BI,gsqr(gdeux)),stoi(3))==1) &&
     (gequal(gmod(AI,gdeux),gun)==1))
   T=itos(gmod(gmul2n(gadd(BI,gun),-2),gsqr(gdeux)));
 z=gmod(BI,gmul2n(gun,5));
 if ((gequal(gmod(AI,gmul2n(gun,4)),gzero)==1) &&
     (gequal(z,gzero) || (gequal(z,gmul2n(gun,3))))==1)
   T=itos(gmod(gmul2n(BI,-3),gdeux));
 if (T&1) a3=gsub(gun,gmod(gmul(a1,a2),gdeux));
 else a3=gmod(gmul(a1,a2),gdeux);
 b2=gadd(a1,gmul2n(a2,2));
 a4=gdiv(gsub(gsqr(b2),gadd(c4,gmul(stoi(24),gmul(a1,a3)))),stoi(48));
 a6=gdiv(gsub(gmul(stoi(36),gmul(b2,gadd(gmul(a1,a3),gmul2n(a4,1)))),
	      gadd(gmul(b2,gsqr(b2)),gadd(c6,gmul(stoi(216),a3)))),
	 stoi(864));
 thevec=cgetg(6,t_VEC);
 thevec[1]=(long) a1; thevec[2]=(long) a2; thevec[3]=(long) a3;
 thevec[4]=(long) a4; thevec[5]=(long) a6;
 //output(a1); output(a2); output(a3); output(a4); output(a6);
 return(ellinit0(thevec,1,ELLACC));
}


GEN mineqfromc4c6(GEN c4,GEN c6)
{GEN AI,BI,F,disc,g,u,p,z; int i,j,number,d;

 u=gcd0(denom(c4),denom(c6),0); c4=gmul(c4,gsqr(gsqr(u)));
 c6=gmul(c6,gmul(gsqr(u),gsqr(gsqr(u))));
 disc=gdiv(gsub(gmul(c4,gsqr(c4)),gsqr(c6)),stoi(1728));
 g=gabs(gcd0(gsqr(c6),disc,0),-1); F=factor(gmul(stoi(6),g));
   ((GEN)F[2])[1] = (long) gsub((GEN) ((GEN) F[2])[1],gun);
   ((GEN)F[2])[2] = (long) gsub((GEN) ((GEN) F[2])[2],gun);

 number=itos((GEN) matsize(F)[1]); u=gun;
 for(i=1;i<=number;i++)
 {p=(GEN) ((GEN) F[1])[i];
  d=itos(gfloor(gdiv((GEN) ((GEN) F[2])[i],stoi(12))));
  if (gequal(p,gdeux))
  {AI=gmul2n(c4,-4*d); BI=gmul2n(c6,-6*d);
   z=gmod(BI,gmul2n(gun,5));
   if ((!gequal(gmod(AI,gdeux),gun) || (!gequal(gmod(BI,gsqr(gdeux)),stoi(3))))
       && ((!gequal(gmod(AI,gmul2n(gun,4)),gzero)) ||
	   (!gequal(z,gzero) && (!gequal(z,gmul2n(gun,3)))))==1) d--;
  }
  if (gequal(p,stoi(3))) {if (ggval(c6,stoi(3))==6*d+2) d--;}
  if (d==-1)
  {c4=gmul(c4,gsqr(gsqr(p))); c6=gmul(c6,gmul(gsqr(p),gsqr(gsqr(p))));}
  for (j=0;j<d;j++)
  {c4=gdiv(c4,gsqr(gsqr(p))); c6=gdiv(c6,gmul(gsqr(gsqr(p)),gsqr(p)));}
 }
 //printf("OUT\n"); output(c4); output(c6);
 return(avecfromc4c6(c4,c6));
}


GEN qtwist(GEN ellc,GEN p)
{GEN c4; GEN c6;

 if (gequal(p,gun)) return(ellc);
 c4=gmul((GEN) ellc[10],gsqr(p));
 c6=gmul(p,gmul((GEN) ellc[11],(GEN) gsqr(p)));
 return(mineqfromc4c6(c4,c6));
}

/* #define DBG(arg) printf("retminimaltwist: %s\n",arg); */
#define DBG(arg)

/* minimize twist at primes bigger than 3, then does 3, then does 2. */
GEN retminimaltwist(GEN EC,GEN *TWEF)
{
  GEN ELLAP,TEMP,currp,INVOL,OUTVOL;
  int number,i,kodsym,twocond,locv,toptwo;

 number = itos((GEN) matsize(CF)[1]); *TWEF=gun;

 DBG("1");

 for(i=1;i<=number;i++)
   {
     DBG("2");
     currp=(GEN) ((GEN) CF[1])[i];
     if (mpcmp(currp,gdeux)==0)
       {twocond=itos((GEN) ((GEN) CF[2])[i]);
	if ((twocond==4) || (twocond==6))
	  {
	    toptwo=twocond; INVOL=myvol(EC);
	    if (twocond==6) EC=qtwist(EC,gdeux);
	    twocond=itos(gel(elllocalred(EC,gdeux),1));
	    if (twocond==4) EC=qtwist(EC,gneg(gun));
	    toptwo-=itos(gel(elllocalred(EC,gdeux),1));
	    OUTVOL=myvol(EC); *TWEF=gmul2n(volratio(OUTVOL,INVOL),toptwo);
	    twocond=itos(gel(elllocalred(EC,gdeux),1));
	    if (twocond==0) {
	      ELLAP=ellap0(EC,gdeux,0); ELLAP=gsub(gmul(ELLAP,ELLAP),gdeux);
	      TEMP=gsub(gmul(gdeux,gsqr(gdeux)),gmul(gdeux,ELLAP));
	      TEMP=gsub(gadd(TEMP,ELLAP),gun); *TWEF=gmul2n(gmul(*TWEF,TEMP),-3);
	    }
	    if (twocond==1)
	      {TEMP=gmul2n(gadd(gdeux,gun),-2); *TWEF=gmul(*TWEF,TEMP);}
	  }
	if (twocond>=7)
	  {if ((ggval(gel(EC,10),gdeux)>=6) &&
	       (ggval((GEN) (EC)[11],gdeux)>=6) &&
	       (ggval((GEN) (EC)[12],gdeux)>=6))
	     {EC=qtwist(EC,gdeux); *TWEF=gmul(*TWEF,gdeux);}
	 }
	if ((twocond>=5) && (gsigne((GEN) EC[12])==-1)) EC=qtwist(EC,gneg(gun));
      }
     else {
       kodsym=itos(gel(elllocalred(EC,currp),2));
       if ((kodsym==-1) || (kodsym<=-5)) {
	 if (mod4(currp)==3) currp=gneg(currp);
	 EC=qtwist(EC,currp);
	 if (gsigne(currp)==-1) currp=gneg(currp);
	 locv=itos(gel(elllocalred(EC,currp),1));
	 if (locv==0) {
	   ELLAP=ellap0(EC,currp,0);
	   ELLAP=gsub(gmul(ELLAP,ELLAP),currp);
	   TEMP=gsub(gmul(currp,gsqr(currp)),gmul(currp,ELLAP));
	   TEMP=gsub(gadd(TEMP,ELLAP),gun);
	   *TWEF=gmul(*TWEF,TEMP); // (p+1+ap)(p+1-ap)(p-1) when good red
	 }
	 if (locv==1) {
	   TEMP=gsub(gsqr(currp),gun); *TWEF=gmul(*TWEF,TEMP);
	 }
       }
       else if (ggval((GEN) EC[10],currp)>=2) {
	 if (ggval((GEN) EC[11],currp)>=3)
	   {if (!gequal(currp,stoi(3)))
	      {if (mod4(currp)==3) currp=gneg(currp); EC=qtwist(EC,currp);
	       if (gsigne(currp)==-1) currp=gneg(currp);
	       *TWEF=gmul(*TWEF,currp);
	       }
	   else
	     {if ((ggval((GEN) EC[11],currp)!=5) &&
		  (ggval((GEN) EC[12],currp)>=6))
		{EC=qtwist(EC,gneg(currp)); *TWEF=gmul(*TWEF,currp);}
	    }
	  }
       }
     }
   }
 return(EC);
}

void minimaltwist()
{GEN ELLAP,TEMP,currp,INVOL,OUTVOL;
 int number,i,kodsym,twocond,locv,toptwo;

 number=itos((GEN) matsize(CF)[1]); TWCOND=gcopy(COND);
 TWCURVE=gcopy(CURVE); TWPROD=gun; TWEFF=gun; TWCF=gcopy(CF);

 for(i=1;i<=number;i++)
 {currp=(GEN) ((GEN) CF[1])[i];
  if (mpcmp(currp,gdeux)==0)
  {twocond=itos((GEN) ((GEN) CF[2])[i]);
   if ((twocond==4) || (twocond==6))
   {toptwo=twocond; INVOL=myvol(TWCURVE);
    if (twocond==6)
    {TWCURVE=qtwist(TWCURVE,gdeux); TWPROD=gdeux;}
    twocond=itos(gel(elllocalred(TWCURVE,currp),1));
    if (twocond==4)
    {TWCURVE=qtwist(TWCURVE,gneg(gun)); TWPROD=gneg(TWPROD);}
    twocond=itos(gel(elllocalred(TWCURVE,currp),1));
    TWCOND=gmul2n(TWCOND,twocond-toptwo); OUTVOL=myvol(TWCURVE);
    TWEFF=gmul2n(volratio(OUTVOL,INVOL),toptwo-twocond);
    TWPROD=gmul(TWPROD,gsqr(gdeux));
    if(twocond==0)
    {ELLAP=ellap0(TWCURVE,gdeux,0); ELLAP=gsub(gmul(ELLAP,ELLAP),gdeux);
     TEMP=gsub(gmul(gdeux,gsqr(gdeux)),gmul(gdeux,ELLAP));
     TEMP=gsub(gadd(TEMP,ELLAP),gun); TWEFF=gmul2n(gmul(TWEFF,TEMP),-3);
    }
    if (twocond==1)
    {TEMP=gmul2n(gadd(gdeux,gun),-2); TWEFF=gmul(TWEFF,TEMP);}
   }
   if (twocond>=7)
   {if ((ggval((GEN) TWCURVE[10],gdeux)>=6) &&
	(ggval((GEN) TWCURVE[11],gdeux)>=6) &&
	(ggval((GEN) TWCURVE[12],gdeux)>=6))
    {TWCURVE=qtwist(TWCURVE,gdeux); TWPROD=gmul(TWPROD,gdeux);
     TWEFF=gmul(TWEFF,gdeux);
    }
   }
   if ((twocond>=5) && (gsigne((GEN) (TWCURVE)[11])==-1))
   {TWCURVE=qtwist(TWCURVE,gneg(gun)); TWPROD=gneg(TWPROD);}
   ((GEN) TWCF[2])[i]=(long) stoi(twocond);
  }
  else
  {kodsym=itos(gel(elllocalred(TWCURVE,currp),2));
   if ((kodsym==-1) || (kodsym<=-5))
   {if (mod4(currp)==3) currp=gneg(currp);
    TWCURVE=qtwist(TWCURVE,currp); TWPROD=gmul(TWPROD,currp);
    if (gsigne(currp)==-1) currp=gneg(currp);
    locv=itos(gel(elllocalred(TWCURVE,currp),1));
    if (locv==0)
    {ELLAP=ellap0(TWCURVE,currp,0); ELLAP=gsub(gmul(ELLAP,ELLAP),currp);
     TEMP=gsub(gmul(currp,gsqr(currp)),gmul(currp,ELLAP));
     TEMP=gsub(gadd(TEMP,ELLAP),gun); TWEFF=gmul(TWEFF,TEMP);
     TWCOND=gdiv(TWCOND,gsqr(currp)); // (p+1+ap)(p+1-ap)(p-1) when good red
    }
    if (locv==1)
    {TEMP=gsub(gsqr(currp),gun); TWEFF=gmul(TWEFF,TEMP);
     TWCOND=gdiv(TWCOND,currp);
    }
    ((GEN) TWCF[2])[i]=(long) stoi(locv);
   }
   else if (ggval((GEN) TWCURVE[10],currp)>=2)
   {if (ggval((GEN) TWCURVE[11],currp)>=3)
    {if (!gequal(currp,stoi(3)))
     {if (mod4(currp)==3) currp=gneg(currp);
      TWCURVE=qtwist(TWCURVE,currp); TWPROD=gmul(TWPROD,currp);
      if (gsigne(currp)==-1) currp=gneg(currp); TWEFF=gmul(TWEFF,currp);
     }
     else
     {if ((ggval((GEN) TWCURVE[11],currp)!=5) &&
	  (ggval((GEN) TWCURVE[12],currp)>=6))
      {TWCURVE=qtwist(TWCURVE,gneg(currp));
       TWPROD=gmul(TWPROD,gneg(currp)); TWEFF=gmul(TWEFF,currp);
      }
     }
    }
   }
  }
 }

 if (gequal(gmod(TWPROD,gsqr(gdeux)),gdeux)) TWPROD=gmul2n(TWPROD,2);
 if (gequal(TWPROD,gneg(gun))) TWPROD=gneg(gsqr(gdeux));
 ISSQFREE=1; ISPRIME=0;
 for (i=1;i<=itos((GEN) matsize(TWCF)[1]);i++)
 {locv=itos((GEN) ((GEN) CF[2])[i]); if (locv>1) ISSQFREE=0; ISPRIME+=locv;}
 if (ISPRIME>1) ISPRIME=0;

 if (0 && VERBOSE)
 {if (gequal(TWPROD,gun)==1) return;
  printf("Minimal Rational Quadratic Twist is ");
  printcurven(TWCURVE); printf("of conductor "); output((GEN) (TWCOND));
  printf("Product of Twisted Primes is "); output(TWPROD);
  printf("Twists give a factor of "); output(TWEFF);
 }
}

GEN minimalequation()
{int i,number; GEN p,l;

 CURVE=mineqfromc4c6((GEN) CURVE[10],(GEN) CURVE[11]);
 CF=factor(gabs((GEN) CURVE[12],-1)); number=itos((GEN) matsize(CF)[1]);
 for(i=1;i<=number;i++)
 {p=(GEN) ((GEN) CF[1])[i];
  ROOTNO*=ellrootno(CURVE,p); l=elllocalred(CURVE,p);
  COND=gmul(COND,gpow(p,(GEN) l[1],-1)); TAMA=gmul(TAMA,(GEN) l[4]);
  ((GEN) CF[2])[i]=lcopy((GEN) l[1]);
 }
 return(CURVE);
}

GEN gettama(GEN C)
{int i,n; GEN p,l,T;

 T=gun; n=itos((GEN) matsize(CF)[1]);
 for(i=1;i<=n;i++)
 {p=(GEN) ((GEN) CF[1])[i]; l=elllocalred(C,p); T=gmul(T,(GEN) l[4]);}
 return(T);
}

void symsq()
{int n,i,myp,mypow; GEN p,pow;

 SSCOND=gcopy(TWCOND); EXPOS=gun; EXNEG=gun; EXEFF=gun;
 n=itos((GEN) matsize(TWCF)[1]);
 for(i=1;i<=n;i++)
 {p=(GEN) ((GEN) TWCF[1])[i];  pow=(GEN) ((GEN) TWCF[2])[i];
  if (itos(pow)>=2)
  {myp=itos(p);
   if (myp==2)
   {mypow=itos(pow);
    if (mypow==2)
    {SSCOND=gmul2n(SSCOND,-1); EXPOS=gmul2n(EXPOS,1);
     EXEFF=gmul2n(gmul(EXEFF,stoi(3)),-1);
    }
    if (mypow==3) {SSCOND=gmul2n(SSCOND,-1); /* Non-Exotic */}
    if (mypow==5) {SSCOND=gmul2n(SSCOND,-2); /* Non-Exotic */}
    if (mypow==7) {SSCOND=gmul2n(SSCOND,-3); /* Non-Exotic */}
    if (mypow==8)
    {if (ggval((GEN) TWCURVE[11],p)>=9)
     {SSCOND=gmul2n(SSCOND,-4); /* Non-Exotic */}
     else
     {SSCOND=gmul2n(SSCOND,-5);
      if (itos((GEN) gmod((GEN) TWCURVE[10],stoi(128)))==32)
      {EXPOS=gmul2n(EXPOS,1); EXEFF=gmul2n(gmul(EXEFF,stoi(3)),-1);}
      if (itos((GEN) gmod((GEN) TWCURVE[10],stoi(128)))==96)
      {EXNEG=gmul2n(EXNEG,1); EXEFF=gmul2n(EXEFF,-1);}
     }
    }
   }
   else if (myp==3)
   {mypow=itos(pow);
    if (mypow==2)
    {SSCOND=gdiv(SSCOND,stoi(3)); EXPOS=gmul(EXPOS,stoi(3));
     EXEFF=gdiv(gmul2n(EXEFF,2),stoi(3));
    }
    if (mypow==3) {SSCOND=gdiv(SSCOND,stoi(3)); /* Non-Exotic */}
    if (mypow==4)
    {SSCOND=gdiv(SSCOND,stoi(9));
     if (itos((GEN) gmod((GEN) TWCURVE[11],stoi(243)))==0)
     {if (itos((GEN) gmod((GEN) TWCURVE[10],stoi(81)))==27)
      {EXNEG=gmul(EXNEG,stoi(3)); EXEFF=gdiv(gmul2n(EXEFF,1),stoi(3));}
      if (itos((GEN) gmod((GEN) TWCURVE[10],stoi(81)))==54)
      {EXPOS=gmul(EXPOS,stoi(3)); EXEFF=gdiv(gmul2n(EXEFF,2),stoi(3));}
     }
     if (itos((GEN) gmod((GEN) TWCURVE[10],stoi(27)))==9)
     {if (itos((GEN) gmod((GEN) TWCURVE[11],stoi(243)))==54)
      {EXPOS=gmul(EXPOS,stoi(3)); EXEFF=gdiv(gmul2n(EXEFF,2),stoi(3));}
      if (itos((GEN) gmod((GEN) TWCURVE[11],stoi(243)))==108)
      {EXNEG=gmul(EXNEG,stoi(3)); EXEFF=gdiv(gmul2n(EXEFF,1),stoi(3));}
      if (itos((GEN) gmod((GEN) TWCURVE[11],stoi(243)))==135)
      {EXNEG=gmul(EXNEG,stoi(3)); EXEFF=gdiv(gmul2n(EXEFF,1),stoi(3));}
      if (itos((GEN) gmod((GEN) TWCURVE[11],stoi(243)))==189)
      {EXPOS=gmul(EXPOS,stoi(3)); EXEFF=gdiv(gmul2n(EXEFF,2),stoi(3));}
     }
    }
    if (mypow==5) {SSCOND=gdiv(SSCOND,stoi(9)); /* Non-Exotic */}
   }
   else
   {SSCOND=gdiv(SSCOND,p);
    if ((myp%12)==1)
    {EXNEG=gmul(EXNEG,p); EXEFF=gdiv(gmul(EXEFF,gsub(p,gun)),p);}
    if ((myp%12)==5)
    {if ((ggval((GEN) TWCURVE[10],p)<2) &&
	 (ggval((GEN) TWCURVE[11],p)>=2))
     {EXNEG=gmul(EXNEG,p); EXEFF=gdiv(gmul(EXEFF,gsub(p,gun)),p);}
     else
     {EXPOS=gmul(EXPOS,p); EXEFF=gdiv(gmul(EXEFF,gadd(p,gun)),p);}
    }
    if ((myp%12)==7)
    {if ((ggval((GEN) TWCURVE[11],p)<2) ||
	 (ggval((GEN) TWCURVE[10],p)>=2))
     {EXNEG=gmul(EXNEG,p); EXEFF=gdiv(gmul(EXEFF,gsub(p,gun)),p);}
     else
     {EXPOS=gmul(EXPOS,p); EXEFF=gdiv(gmul(EXEFF,gadd(p,gun)),p);}
    }
    if ((myp%12)==11)
    {EXPOS=gmul(EXPOS,p); EXEFF=gdiv(gmul(EXEFF,gadd(p,gun)),p);}
   }
  }
 }

 if (0 && (VERBOSE) && (PRINT))
 {if (gequal(SSCOND,TWCOND)==1) return;
  printf("Symmetric Square Conductor is "); output(SSCOND);
  if (gequal(EXPOS,gun)!=1)
  {printf("Plus product is "); output(EXPOS);}
  if (gequal(EXNEG,gun)!=1)
  {printf("Minus product is ");output(EXNEG);}
  printf("Even valuation gives a factor of "); output(EXEFF);
 }
}

