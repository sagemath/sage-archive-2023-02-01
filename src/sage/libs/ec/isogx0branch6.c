#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int x0isog6(GEN zn,int p,int q)
{GEN listofd,TEMP,getperiods,imagp,realp,linvec,CURV1; int i; int cman=1;

 listofd=getd1(zn,p,stoi(p));
 listofd=geval(gtoset(concat(listofd,getd1(zn,p*p,stoi(p)))));
 listofd=geval(gtoset(concat(listofd,getd1(zn,q,stoi(q)))));

 if (ISOG==12)
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,12);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,12);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(12)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[1],stoi(3)),gzero)!=1) cman+=4;
    if (gequal(gmod((GEN) linvec[1],gmul2n(gdeux,1)),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[1],gmul2n(gdeux,1)),gzero)==1) cman++;
     else{printf("Not implemented yet\n"); return(0);}
    }
   }
   if (gequal(gmod((GEN) linvec[2],stoi(12)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[2],stoi(3)),gzero)!=1) cman+=4;
    if (gequal(gmod((GEN) linvec[2],gmul2n(gdeux,1)),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[2],gmul2n(gdeux,1)),gzero)==1) cman++;
     else{printf("Not implemented yet\n"); return(0);}
    }
   }
  }
  return(cman);
 }
 if (ISOG==18)
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,18);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,18);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(18)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[1],gdeux),gzero)!=1) cman+=3;
    if (gequal(gmod((GEN) linvec[1],stoi(9)),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[1],stoi(3)),gzero)!=1) cman+=2; else cman++;}
   }
   if (gequal(gmod((GEN) linvec[2],stoi(6)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[2],gdeux),gzero)!=1) cman+=3;
    if (gequal(gmod((GEN) linvec[2],stoi(9)),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[2],stoi(3)),gzero)!=1) cman+=2; else cman++;}
   }
  }
  return(cman);
 }
 printf("ERROR in x0isog6\n"); return(0);
}
