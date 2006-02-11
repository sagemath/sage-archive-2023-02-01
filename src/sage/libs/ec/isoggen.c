#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int doisogreal(GEN CURV1,GEN *CURV2,int p)
{return(doisoggeneric(CURV1,CURV2,p,0));}
int doisogimag(GEN CURV1,GEN *CURV2,int p)
{return(doisoggeneric(CURV1,CURV2,p,1));}

void shiftperiodi(int s,GEN *IP,GEN RP,int p)
{if (s==1) *IP=gdiv(*IP,stoi(p));
 else *IP=gadd(gmul2n(RP,-1),gmul(gi,gdiv(gimag(*IP),stoi(p))));
}

void shiftperiodr(int s,GEN *IP,GEN *RP,int p)
{*RP=gdiv(*RP,stoi(p));
 if (s==-1) *IP=gadd(gmul2n(*RP,-1),gmul(gi,gimag(*IP)));
}

int doisoggeneric(GEN CURV1,GEN *CURV2,int p,int w)
{GEN TEMP,realp,imagp,mytau,thevec,front; int n,j; int u=5;
 GEN mysum=gzero; GEN firstq,movingq,myc4,myc6,EARLYEXIT;

 EARLYEXIT=gsqr(dbltor(0.000000000001)); thevec=cgetg(6,t_VEC);
 TEMP=periodvolvec0(CURV1,u);
 realp=(GEN) TEMP[2]; imagp=(GEN) TEMP[3];
 if (w==1) shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
 else shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
 mytau=fundamentaltau(&realp,&imagp);

 front=gsqr(gdiv(gmul2n(mppi(u),1),imagp));
 u=itos(gceil(gdiv(gadd(dbltor(1.26),
			gmul(stoi(3),glogagm(gabs(front,3),3))),
		   dbltor(11.5125)))); u++;
 if (u>5)
 {TEMP=periodvolvec0(CURV1,u);
  realp=(GEN) TEMP[2]; imagp=(GEN) TEMP[3];
  if (w==1) shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
  else shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
  mytau=fundamentaltau(&realp,&imagp);
  front=gsqr(gdiv(gmul2n(mppi(u),1),imagp));
 }
 else u=5;

 firstq=gexp(gmul(gmul(gi,mytau),gmul2n(mppi(u),1)),u); movingq=gcopy(firstq);
 myc4=gsqr(front);
 for(n=1;n<=3*u;n++)
 {mysum=gadd(mysum,gdiv(gmul(movingq,gmul(stoi(n),gsqr(stoi(n)))),
			  gsub(gun,movingq)));
  movingq=gmul(movingq,firstq);
  if (mpcmp(gnorm(gmul(movingq,myc4)),EARLYEXIT)==-1) n=3*u;
 }
 mysum=gadd(gun,gmul(stoi(240),mysum));
 myc4=gmul(myc4,mysum);
 //output(myc4);
 if (mpcmp(gnorm(gimag(myc4)),dbltor(0.000001))==1) return(0);
 myc4=greal(myc4);
 if (mpcmp(gnorm(gsub(myc4,ground(myc4))),dbltor(0.000001))==1) return(0);
 myc4=ground(myc4);

 movingq=gcopy(firstq); mysum=gzero;
 myc6=gmul(gsqr(front),front);
 for(n=1;n<=3*u;n++)
 {mysum=gadd(mysum,gdiv(gmul(movingq,gmul(stoi(n),gsqr(gsqr(stoi(n))))),
			  gsub(gun,movingq)));
  movingq=gmul(movingq,firstq);
  if (mpcmp(gnorm(gmul(movingq,myc6)),EARLYEXIT)==-1) n=3*u;
 }
 mysum=gadd(gun,gmul(stoi(-504),mysum));
 myc6=gmul(myc6,mysum);
 //output(myc6);
 if (mpcmp(gnorm(gimag(myc6)),dbltor(0.000001))==1) return(0);
 myc6=greal(myc6);
 if (mpcmp(gnorm(gsub(myc6,ground(myc6))),dbltor(0.000001))==1) return(0);
 myc6=ground(myc6);

 if (mpcmp(gmod(gsub(gmul(myc4,gsqr(myc4)),gsqr(myc6)),stoi(1728)),gzero)!=0)
     return(0);
 if (ggval(myc6,stoi(3))==2) return(0);
 if (mod2(gabs(myc4,-1)))
 {if (mpcmp(gmod(myc6,gsqr(gdeux)),stoi(3))!=0) return(0);}
 else
 {if (mod16(gabs(myc4,-1))) return(0);
  j=itos(gmod(myc6,gmul2n(gun,5))); if ((j!=0) && (j!=8)) return(0);
 }
 *CURV2=mineqfromc4c6(myc4,myc6);
 if (gequal((GEN) (*CURV2)[10],(GEN) CURV1[10]) &&
     gequal((GEN) (*CURV2)[11],(GEN) CURV1[11])) return(0);
 return(checkconductor(*CURV2));
}

