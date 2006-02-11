#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int mystop(double D,int r)
{if ((r==0) && (D<0.1)) return(0); if (D<0.004) return(0); return(1);}

GEN mytors(GEN C)
{byteptr pptr=diffptr; ulong pr=0; int i=0; GEN discr=(GEN) C[12]; GEN ap;
 GEN thegcd=gzero; GEN pp;

 return((GEN) elltors0(C,0)[1]);

/*
 pr+=*pptr; pptr++;
 for (i=0;i<24;i++)
 {pr+=*pptr; pptr++; pp=stoi(pr);
  if (ggval(discr,pp)==0)
  {ap=ellap0(C,pp,0); thegcd=gcd0(thegcd,gsub(gadd(pp,gun),ap),0);}
 }

 if (gequal(thegcd,gun)==1) return(gun);
 else return((GEN) elltors0(C,0)[1]);
*/
}

int analrank(GEN curv,double *lrat)
{int b,i,r,BND; double rpadj,SUM,alpha; GEN pv, curv2,u;

 if (!INITANALRANK) analinitw();
 Manal=2.0*rtodbl(gdiv(mppi(4),gsqrt(COND,4)));
 if (ROOTNO==1) r=0; else r=1;
 rpadj=rtodbl((GEN) periodvolvec0(CURVE,4)[2]);
 if (gsigne((GEN) CURVE[12])==1) rpadj=rpadj*2.0;
 rpadj*=(double) rtodbl(gmul(TAMA,dbltor(1.0)));
 curv2 = ellinit0(CURVE, 0, ELLACC);

/*todo cut out
 curv2=cgetg(20,t_VEC); for (i=1;i<=13;i++) curv2[i]=CURVE[i];
 b=4+((lgefint((GEN) CURVE[12])-2)>>1);
 pv=periodvolvec0(CURVE,b+2); curv2[15]=pv[2]; curv2[16]=pv[3];
 curv2[14]=(long) gzero; curv2[17]=(long) gzero;
 curv2[18]=(long) gzero; curv2[19]=(long) gzero;
*/

 rpadj/=(double) itos(gsqr(mytors(curv2)));
 u=gpow(gdeux,gsubgs(gdiv(glog(COND,4),glog(stoi(10),4)),6),4);
 alpha=1.74e-09/rpadj; if (mpcmp(gun,u)==-1) alpha*=rtodbl(u);
 // output(u); printf("%f\n",alpha);
 if (alpha>0.00001) BND=(int) ceil(3*rtodbl(gsqrt(COND,4)));
 else BND=(int) ceil(2*rtodbl(gsqrt(COND,4))); if (BND<1000) BND=1000;
 computeap(TWCURVE,BND);

 if (!gequal(TWPROD,gun))
 {computetwistarray(itos(TWPROD),BND);
  while (1)
  {SUM=0.0; for (i=1;i<=BND;i++) {SUM+=analw((int64) i,r)*antwarray[i];}
   SUM*=2.0; *lrat=SUM/rpadj; if (TRACE) printf("%f %f\n",SUM,SUM/rpadj);
   if (mystop(SUM/rpadj,r)) return(r); r+=2;
  }
 }
 else
 {while (1)
  {SUM=0.0; for (i=1;i<=BND;i++) {SUM+=analw((int64) i,r)*anarray[i];}
   SUM*=2.0; *lrat=SUM/rpadj; if (TRACE) printf("%f %f\n",SUM,SUM/rpadj);
   if (mystop(SUM/rpadj,r)) return(r); r+=2;
  }
 }
}

int qtwistar(int d,double *lrat,GEN thecond)
{int b,i,r,eps,BND,number; GEN l,p,pv,curv2,qcurv,twtama,u,FACT,TC;
 double SUM,rpadj,alpha;

 if (!INITANALRANK) analinitw(); TC=gcopy(thecond);
 Manal=2.0*rtodbl(gdiv(mppi(4),gsqrt(thecond,4)));
 qcurv=qtwist(TWCURVE,stoi(d)); eps=-1; twtama=gun;
 number=itos((GEN) matsize(TWCF)[1]);
 for(i=1;i<=number;i++)
 {p=(GEN) ((GEN) TWCF[1])[i];
  eps*=ellrootno(qcurv,p); l=elllocalred(qcurv,p);
  twtama=gmul(twtama,(GEN) l[4]);
  while (gequal(gmod(thecond,p),gzero)) thecond=gdiv(thecond,p);
 }
 FACT=factor(thecond); number=itos((GEN) matsize(FACT)[1]);
 for(i=1;i<=number;i++)
 {p=(GEN) ((GEN) FACT[1])[i]; eps*=ellrootno(qcurv,p);
  l=elllocalred(qcurv,p); twtama=gmul(twtama,(GEN) l[4]);
 }
 if (eps==1) r=0; else r=1;

 rpadj=rtodbl((GEN) periodvolvec0(qcurv,4)[2]);
 if (gsigne((GEN) qcurv[12])==1) rpadj*=2.0;
 rpadj*=gtodouble(twtama);
 curv2=cgetg(20,t_VEC); for (i=1;i<=13;i++) curv2[i]=qcurv[i];
 b=4+((lgefint((GEN) qcurv[12])-2)>>1);
 pv=periodvolvec0(qcurv,b+2); curv2[15]=pv[2]; curv2[16]=pv[3];
 curv2[14]=(long) gzero; curv2[17]=(long) gzero;
 curv2[18]=(long) gzero; curv2[19]=(long) gzero;
 rpadj/=(double) itos(gsqr(mytors(curv2)));

 u=gpow(gdeux,gsubgs(gdiv(glog(TC,4),glog(stoi(10),4)),6),4);
 alpha=1.74e-09/rpadj; if (mpcmp(gun,u)==-1) alpha*=rtodbl(u);
 if (alpha>0.00001) BND=(int) ceil(3*rtodbl(gsqrt(TC,4)));
 else BND=(int) ceil(2*rtodbl(gsqrt(TC,4))); if (BND<1000) BND=1000;
 computeap(TWCURVE,BND); computetwistarray(d,BND);
 while (1)
 {SUM=0.0; for (i=1;i<=BND;i++) {SUM+=analw((int64) i,r)*antwarray[i];}
  SUM*=2.0; *lrat=SUM; if (mystop(SUM/rpadj,r)) return(r); r+=2;
 }
}




