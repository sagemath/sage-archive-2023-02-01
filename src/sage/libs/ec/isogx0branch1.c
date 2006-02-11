#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int x0isog1(GEN zn,int p)
{GEN listofd; GEN TEMP; GEN getperiods; GEN imagp; GEN realp; GEN linvec;
 int i; int j=1; GEN CURV1;
GEN dbg;
int dbg2;
 listofd=getd1(zn,p,stoi(p)); CURV1=ellinit0((GEN) CL[1],1,ELLACC);
 //27 unneeded as 18 divides phi(N)

 if ((ISOG>0) || (ISOG==-18))
 {if ((ISOG==18) && (p==2)) j=4;
  if ((ISOG==-18) && (p==2)) j=4;
  if ((ISOG==15) && (p==5)) j=3;
  if ((ISOG==21) && (p==7)) j=3;
  if (j==1) j=2;
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   if (gequal(gmod((GEN) linvec[1],stoi(p)),gzero)!=1) return(j);
   if (gequal(gmod((GEN) linvec[2],stoi(p)),gzero)!=1) return(j);
  }
  return(1);
 }
 if ((ISOG==-4) || (ISOG==-12))
 {getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,2);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,2);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod(((GEN) linvec[1]),gdeux),gzero)!=1)
   {if (gequal(gmod(((GEN) linvec[2]),gdeux),gzero)!=1)
     return(bentisog((GEN) CL[2],(GEN) CL[3],(GEN) CL[4]));
    else return(realisog2((GEN) CL[2],(GEN) CL[3],(GEN) CL[4],realp));
   }
   if (gequal(gmod(((GEN) linvec[2]),gdeux),gzero)!=1)
     return(imagisog2((GEN) CL[2],(GEN) CL[3],(GEN) CL[4],imagp));
  }
  return(1);
 }
 if ((ISOG==-9) || (ISOG==-25))
 {getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(p)),gzero)!=1)
     return(realisog((GEN) CL[2],(GEN) CL[3],realp));
   if (gequal(gmod((GEN) linvec[2],stoi(p)),gzero)!=1)
     return(imagisog((GEN) CL[2],(GEN) CL[3],imagp));
  }
  return(1);
 }
 printf("ERROR in x0isog1\n"); return(1);
}
