#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int x0array4[4][4]={{0,1,0,1},{1,-1,5,-1},{0,5,6,5},{1,-1,5,-1}};
int x0array4B[4][4]={{0,5,0,5},{5,-1,1,-1},{0,1,6,1},{5,-1,1,-1}};

int x0isog4(GEN zn,int p)
{GEN listofd,TEMP,getperiods,imagp,realp,linvec,CURV1; int i,rm,im,tip;

 listofd=getd1(zn,2,gdeux);
 listofd=geval(gtoset(concat(listofd,getd1(zn,4,gdeux))));
 listofd=geval(gtoset(concat(listofd,getd1(zn,8,gdeux))));
 listofd=geval(gtoset(concat(listofd,getd1(zn,16,gdeux))));
 //nasty 16-isog structure
 return(1); // FIX THIS FIX THIS

 CURV1=ellinit0((GEN) CL[1],1,ELLACC);
 if (gsigne((GEN) CURV1[12])==-1)
 {CURV1=ellinit0((GEN) CL[6],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,2);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,2);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   output(linvec);
   rm=itos((GEN) linvec[1]); if (rm<0) rm=-rm;
   im=itos((GEN) linvec[2]); if (im<0) im=-im;
   tip=x0array4[rm&3][im&3];
   switch(tip)
   {case 1: return(6+whichbent((GEN) CL[8],(GEN) CL[7]));
    case -1: return(6+whichbent((GEN) CL[7],(GEN) CL[8]));
    case 3: {if (rm&3) return(5); else return(3);}
    case 6: return(6);
    case 0: continue;
   }
   TEMP=(GEN) periodvolvec(ellinit0((GEN) CL[2],1,ELLACC))[2];
   if (gequal(ground(gdiv(TEMP,gmul2n(realp,1))),gun)==1)
    {if (im&2) return(6); if (im&4) return(4);}
   else
   {if (rm&2) return(4); if (!(im&2)) return(2);}
  }
  return(1);
 }
 printf("Error in x0isog4\n"); return(0);
}
