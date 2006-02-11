#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int x0isog2(GEN zn,int p)
{GEN listofd,TEMP,getperiods,imagp,realp,linvec,CURV1; int i;

 listofd=getd1(zn,p,stoi(p));
 listofd=geval(gtoset(concat(listofd,getd1(zn,p*p,stoi(p)))));

 if ((ISOG==9) || (ISOG==25))
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(p*p)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[1],stoi(p)),gzero)!=1) return(3);
    else return(2);
   }
   if (gequal(gmod((GEN) linvec[2],stoi(p*p)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[2],stoi(p)),gzero)!=1) return(3);
    else return(2);
   }
  }
  return(1);
 }
 if (ISOG==27)
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(p*p)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[1],stoi(p)),gzero)!=1) return(4);
    else return(imagisog((GEN) CL[2],(GEN) CL[3],gmul(imagp,stoi(9))));
   }
   if (gequal(gmod((GEN) linvec[2],stoi(p*p)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[2],stoi(p)),gzero)!=1) return(4);
    else return(realisog((GEN) CL[2],(GEN) CL[3],gmul(realp,stoi(9))));
   }
  }
  return(1);
 }
 if ((ISOG==4) || (ISOG==8) || (ISOG==12) || (ISOG==16) || (ISOG==32))
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  if (gsigne((GEN) CURV1[12])==-1)
  {CURV1=ellinit0((GEN) CL[2],1,ELLACC);
   getperiods=periodvolvec(CURV1);
   realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
   shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,2);
   shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,2);
   for (i=1;i<=glength(listofd);i++)
   {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   //output(TEMP);
    linvec=linearcomb(realp,imagp,TEMP);
    //output(linvec);
    if (gequal(gmod(((GEN) linvec[1]),gdeux),gzero)!=1)
      return(1+realisog((GEN) CL[3],(GEN) CL[4],realp));
    if (gequal(gmod(((GEN) linvec[2]),gdeux),gzero)!=1)
      return(1+imagisog((GEN) CL[3],(GEN) CL[4],imagp));
    if (gequal(gmod(gadd((GEN) linvec[1],(GEN) linvec[2]),gmul2n(gdeux,1)),
	      gzero)!=1) return(2);
   }
   return(1);
  }
  else
  {CURV1=ellinit0((GEN) CL[2],1,ELLACC);
   getperiods=periodvolvec(CURV1);
   realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
   shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,2);
   shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,2);
   for (i=1;i<=glength(listofd);i++)
   {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
    linvec=linearcomb(realp,imagp,TEMP);
    //output(linvec);
    if ((gequal(gmod(((GEN) linvec[1]),gdeux),gzero)!=1))
    {if (gequal(gmod(((GEN) linvec[2]),gdeux),gzero)!=1)
      return(2+whichbent((GEN) CL[3],(GEN) CL[4]));
     else return(2+whichbent((GEN) CL[4],(GEN) CL[3]));
    }
    if ((gequal(gmod(((GEN) linvec[2]),gdeux),gzero)!=1))
      return(2+whichbent((GEN) CL[4],(GEN) CL[3]));
    if ((gequal(gmod(((GEN) linvec[1]),gmul2n(gdeux,1)),gdeux)==1) &&
	(gequal(gmod(((GEN) linvec[2]),gmul2n(gdeux,1)),gdeux)==1)) return(2);
    getperiods=periodvolvec(ellinit0((GEN) CL[1],1,ELLACC));
    if (gequal(ground(gdiv((GEN) getperiods[2],gmul2n(realp,1))),gun)==1)
    {if (gequal(gmod(((GEN) linvec[2]),gmul2n(gdeux,1)),gdeux)==1) return(2);}
    else
    {if (gequal(gmod(((GEN) linvec[1]),gmul2n(gdeux,1)),gdeux)==1) return(2);}
    //BUGGY
   }
   return(1);
  }
 }
 if ((ISOG==-8) || (ISOG==-16) || (ISOG==-32)) {return(1);}
 printf("Error in x0isog2\n"); return(0);
}
