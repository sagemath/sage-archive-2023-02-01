#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

GEN fundamentaltau(GEN *A,GEN *B)
{GEN TEMP; int OK=0;

 if (cmprs(gimag(gdiv(*A,*B)),0)==-1) {TEMP=*A; *A=*B; *B=TEMP;}
 while(!OK)
 {OK=1;
  TEMP=gdiv(*A,*B);
  if(mpcmp(gnorm(TEMP),dbltor(0.999))==-1)
  {OK=0; TEMP=*A; *A=gneg(*B); *B=TEMP; TEMP=gdiv(*A,*B);}
  TEMP=greal(TEMP);
  if(mpcmp(gnorm(TEMP),gnorm(dbltor(0.501)))==1)
  {if (gequal(ground(TEMP),gzero)==1) TEMP=stoi(gsigne(TEMP));
   else TEMP=ground(TEMP);
   *A=gsub(*A,gmul(*B,TEMP)); OK=0;
  }
 }
 return(gdiv(*A,*B));
}

void orderrealroots(GEN *theroots)
{int TEMP;
 if (cmprr((GEN) (*theroots)[1],(GEN) (*theroots)[2])==-1)
 {TEMP=(*theroots)[1]; (*theroots)[1]=(*theroots)[2]; (*theroots)[2]=TEMP;}
 if (cmprr((GEN) (*theroots)[2],(GEN) (*theroots)[3])==-1)
 {TEMP=(*theroots)[2]; (*theroots)[2]=(*theroots)[3]; (*theroots)[3]=TEMP;}
 if (cmprr((GEN) (*theroots)[1],(GEN) (*theroots)[2])==-1)
 {TEMP=(*theroots)[1]; (*theroots)[1]=(*theroots)[2]; (*theroots)[2]=TEMP;}
}

void allcubicroots(GEN poly,GEN *roots,int acc)
{GEN F=(GEN) poly[2]; GEN G=(GEN) poly[4];
 GEN R,Y,Z,T,U;  // leading coeff is 1, poly[3] is irrel

 R=cubicrealroot(poly,acc); (*roots)=cgetg(4,t_VEC); (*roots)[1]=(long) R;
 Y=gadd(F,R); Z=gneg(gdiv(G,R)); U=gneg(gmul2n(Y,-1));
 T=gmul2n(gsqrt(gsub(gsqr(Y),gmul2n(Z,2)),acc),-1);
 (*roots)[2]=(long) gadd(U,T); (*roots)[3]=(long) gsub(U,T);
 return;
}

GEN cubicrealroot(GEN poly,int acc)
{GEN A=(GEN) poly[2]; GEN B=(GEN) poly[3]; GEN C=(GEN) poly[4];
 GEN Q,R,S,T,U,V;

 Q=gdivgs(gsub(gmulsg(3,B),gsqr(A)),9);
 R=gmul2n(gmul(A,gsqr(A)),1); R=gadd(gmulsg(27,C),R);
 R=gsub(gmulsg(9,gmul(A,B)),R); R=gdivgs(R,54);
 U=gadd(gmul(Q,gsqr(Q)),gsqr(R));
 if (gequal(gsqr(R),U))
 {S=gmul2n(R,1);
  if (gsigne(S)>=0) S=gsqrtn(S,stoi(3),&V,acc);
  else S=gneg(gsqrtn(gneg(S),stoi(3),&V,acc));
  return(gsub(S,gdiv(A,stoi(3))));
 }
 if (gsigne(U)>=0)
 {U=gsqrt(U,acc+1);
  S=gadd(R,U); T=gsub(R,U);
  if (gsigne(S)>=0) S=gsqrtn(S,stoi(3),&V,acc);
  else S=gneg(gsqrtn(gneg(S),stoi(3),&V,acc));
  if (gsigne(T)>=0) T=gsqrtn(T,stoi(3),&V,acc);
  else T=gneg(gsqrtn(gneg(T),stoi(3),&V,acc));
  return(gsub(gadd(S,T),gdiv(A,stoi(3))));
 }
 U=gsqrt(U,acc+1); S=gsqrtn(gadd(R,U),stoi(3),&V,acc);
 V=gsub(gmul2n(greal(S),1),gdiv(A,stoi(3)));
 return(V);
}

GEN periodvolvec0(GEN ec,int z)
{GEN A,B,theroots,realroot,rp,ip,retvec,thevec;

 thevec=cgetg(5,t_VEC); thevec[1]=(long) gun;
 thevec[2]=(long) gmul2n((GEN) ec[6],-2);
 thevec[3]=(long) gmul2n((GEN) ec[7],-1);
 thevec[4]=(long) gmul2n((GEN) ec[8],-2); retvec=cgetg(4,t_VEC);
 retvec[1]=(long) gzero; retvec[2]=(long) gzero; retvec[3]=(long) gzero;

 if (gsigne((GEN) ec[12])==1)
 {allcubicroots(thevec,&theroots,z+2);  orderrealroots(&theroots);
  rp=gdiv(mppi(z),agm(gsqrt(gsub((GEN) theroots[1],(GEN) theroots[3]),z),
		      gsqrt(gsub((GEN) theroots[1],(GEN) theroots[2]),z),z));
  ip=gdiv(mppi(z),agm(gsqrt(gsub((GEN) theroots[1],(GEN) theroots[3]),z),
		      gsqrt(gsub((GEN) theroots[2],(GEN) theroots[3]),z),z));
  ip=gmul(gi,ip);
  retvec[1]=(long) gmul(rp,gimag(ip));
  retvec[2]=(long) rp; retvec[3]=(long) ip;
 }
 else
 {realroot=cubicrealroot(thevec,z+2);
  A=gadd(gmul(stoi(3),realroot),gmul2n((GEN) ec[6],-2));
  B=gsqrt(gadd(gadd(gmul(stoi(3),gsqr(realroot)),
		    gmul(realroot,gmul2n((GEN) ec[6],-1))),
	       gmul2n((GEN) ec[7],-1)),z);
  rp=gdiv(gmul2n(mppi(z),1),agm(gmul2n(gsqrt(B,3),1),
				gsqrt(gadd(gmul2n(B,1),A),z),z));
  ip=gdiv(mppi(z),agm(gmul2n(gsqrt(B,3),1),
		      gsqrt(gsub(gmul2n(B,1),A),z),z));
  ip=gadd(gmul2n(gneg(rp),-1),gmul(gi,ip));
  retvec[1]=(long) gmul(rp,gimag(ip));
  retvec[2]=(long) rp; retvec[3]=(long) ip;
 }
 return(retvec);
}

GEN volratio(GEN vola,GEN volb)
{if (mpcmp(vola,volb)==1) return(ground(gdiv(vola,volb)));
 else return(ginv(ground(gdiv(volb,vola))));
}

GEN periodvolvec(GEN ec) {return(periodvolvec0(ec,3));}
GEN myvol(GEN ec) {return((GEN) periodvolvec0(ec,3)[1]);}
