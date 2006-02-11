#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

double looktable[9][7]=
{{0.0000021,0.0000033,0.0000043,0.0000057,0.0000067,0.0000084,0.0000115},
 {0.0000044,0.0000064,0.0000078,0.0000098,0.0000113,0.0000136,0.0000178},
 {0.0000080,0.0000109,0.0000128,0.0000155,0.0000176,0.0000205,0.0000260},
 {0.0000131,0.0000169,0.0000195,0.0000231,0.0000257,0.0000294,0.0000363},
 {0.0000199,0.0000249,0.0000282,0.0000327,0.0000359,0.0000406,0.0000490},
 {0.0000286,0.0000349,0.0000390,0.0000446,0.0000485,0.0000542,0.0000643},
 {0.0000396,0.0000473,0.0000523,0.0000590,0.0000638,0.0000705,0.0000825},
 {0.0000530,0.0000623,0.0000682,0.0000762,0.0000818,0.0000898,0.0001037},
 {0.0000690,0.0000801,0.0000871,0.0000965,0.0001029,0.0001122,0.0001283}};

GEN linearcomb(GEN rp,GEN ip,GEN per)
{GEN TEMPi; GEN TEMPr; GEN vect;

 if (TRACE) {output(rp); output(ip); output(per);}
 vect=cgetg(3,t_VEC);
 TEMPi=bestappr(gdiv(gimag(per),gimag(ip)),gdeux);
 if (typ(TEMPi)!=t_INT) return(0);
 per=greal(gsub(per,gmul(TEMPi,ip)));
 TEMPr=bestappr(gdiv(per,rp),gdeux);
 if (typ(TEMPr)!=t_INT) return(0);
 vect[1]=(long) TEMPr; vect[2]=(long) TEMPi;
 return(vect);
}

int numterms(double Nsqrt,int d,double needacc)
{int L=(int) ceil(log(Nsqrt)/log(10.0)*2.0);
 int dshift=0; int Lshift=0; double haveacc,steps,isteps;

 //printf("NUMTERMS %f %i %f %i\n",Nsqrt,d,needacc,L);
 if (L<2) L=2; if (L>10) {Lshift=L-10; L=10;}
 if (d<=3) dshift=d-1; if (d==5) dshift=3; if (d==7) dshift=4;
 if (d==11) dshift=5; if (d>=13) dshift=6;
 haveacc=looktable[L-2][dshift];
 //printf("HAVEACC %f\n",haveacc);
 if (Lshift!=0) {haveacc+=Lshift*(looktable[8][dshift]-looktable[7][dshift]);}
 if (dshift==6) {haveacc+=((d/24)+1)*(looktable[L-2][6]-looktable[L-2][5]);}
 steps=log(haveacc/needacc)/log(1.9);
 //printf("STEPS %f\n",steps);
 if (steps>0) {isteps=ceil(steps); return(ceil(d*Nsqrt*pow(1.05,isteps)));}
 else {isteps=ceil(-steps); return(ceil(d*Nsqrt*pow(1.05,-isteps)));}
}

double minp(GEN rp,GEN ip)
{GEN rpp; GEN ipp;
 if (gequal(greal(ip),gzero))
 {if (mpcmp(gnorm(rp),gnorm(gimag(ip)))==1) return(rtodbl(gimag(ip)));
  else return(rtodbl(rp));
 }
 else
 {rpp=gmul2n(rp,-1); ipp=gimag(ip);
  if (mpcmp(gnorm(rpp),gnorm(ipp))==1) return(rtodbl(ipp));
  else return(rtodbl(rpp));
 }
}

GEN computeperiod(GEN curv,GEN d,double minp)
{GEN a,b,c; GEN barray; GEN carray; GEN currdexp=gun; GEN currnexp=gun;
 GEN TEMP; GEN total=gzero; GEN front; GEN ddecay; GEN ndecay;
 int dd; int i; int epsi; int n; int nbound;
 pari_sp memspot; pari_sp topmem=avma;

 nbound=numterms(rtodbl(gsqrt(COND,4)),itos(d),
		 minp/100.0/rtodbl(glog(COND,4)));
 //printf("%i %i %f\n",nbound,itos(d),minp);
 computeap(TWCURVE,nbound);
 if (!gequal(TWPROD,gun)) computetwistarray(itos(TWPROD),nbound);
 dd=itos(d); epsi=-ROOTNO;
 barray=cgetg(dd+1,t_VEC); carray=cgetg(dd+1,t_VEC);
 a=lift(ginv(gmodulcp(d,COND))); b=gdiv(gsub(gmul(a,d),gun),COND); c=gun;
 carray[dd]=(long) gexp(gdiv(gmul(gmul2n(mppi(3),1),gi),d),3);
 barray[dd]=(long) gpow((GEN) carray[dd],b,3);
 for (i=dd-1;i>=1;i--)
   carray[i]=(long) gmul((GEN) carray[i+1],(GEN) carray[dd]);
 for (i=dd-1;i>=1;i--)
   barray[i]=(long) gmul((GEN) barray[i+1],(GEN) barray[dd]);

 TEMP=gneg(gdiv(gmul2n(mppi(3),1),gsqrt(COND,3)));
 ddecay=gexp(gdiv(TEMP,d),3); ndecay=gexp(TEMP,3);

 memspot=avma;
 for (n=1;n<=nbound;n++)
 {if (gequal(TWPROD,gun)) front=gdiv(stoi((int) anarray[n]),stoi(n));
  else front=gdiv(stoi((int) antwarray[n]),stoi(n));
  currdexp=gmul(currdexp,ddecay);
  if (epsi==1)
  {TEMP=gsub((GEN) barray[1+(n%dd)],(GEN) carray[1+(n%dd)]);
   TEMP=gmul(currdexp,TEMP); total=gadd(total,gmul(TEMP,front));
  }
  else
  {TEMP=gadd((GEN) barray[1+(n%dd)],(GEN) carray[1+(n%dd)]);
   TEMP=gmul(currdexp,TEMP); currnexp=gmul(currnexp,ndecay);
   TEMP=gsub(TEMP,gmul2n(currnexp,1)); total=gadd(total,gmul(TEMP,front));
  }
  total=gcopy(total); currdexp=gcopy(currdexp); currnexp=gcopy(currnexp);
  avma=memspot; memspot=avma;
  total=gcopy(total); currdexp=gcopy(currdexp); currnexp=gcopy(currnexp);
 }
 avma=topmem; total=gcopy(total);
 return(total);
}
