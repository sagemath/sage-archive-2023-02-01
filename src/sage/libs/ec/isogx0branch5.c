#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int whichone5(int isog,int cman)
{if ((isog!=12) && (isog!=18))
 {if (cman==1) return(1); if (cman==2) return(2); if (cman==6) return(4);
  if (cman==3) {if (isog==6) return(3); else return(2);}
  if ((cman==5) || (cman==7)) return(3);
  if ((cman==10) || (cman==14)) return(4);
  if ((cman==15) || (cman==21)) return(4);
 }
 if (isog==12) {if (cman==3) return(5); return(cman);}
 if (isog==18)
 {if (cman==3) return(2); if (cman==2) return(3);
  if (cman==6) return(5); return(1);
 }
 return(-1);
}

int x0isog5(GEN zn,int p,int q)
{GEN listofd,TEMP,getperiods,imagp,realp,linvec,CURV1; int i; int cman=1;

 listofd=getd1(zn,p,stoi(p));
 listofd=geval(gtoset(concat(listofd,getd1(zn,q,stoi(q)))));

 if (ISOG>0)
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,p*q);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,p*q);
  //lame but it works --- DOES 2-isog work PROPERLY?!
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if ((gequal(gmod((GEN) linvec[1],stoi(p*q)),gzero)!=1))
   {if (gequal(gmod((GEN) linvec[1],stoi(p)),gzero)!=1) cman*=p;
    if (gequal(gmod((GEN) linvec[1],stoi(q)),gzero)!=1) cman*=q;
   }
   if ((gequal(gmod((GEN) linvec[2],stoi(p*q)),gzero)!=1))
   {if (gequal(gmod((GEN) linvec[2],stoi(p)),gzero)!=1) cman*=p;
    if (gequal(gmod((GEN) linvec[2],stoi(q)),gzero)!=1) cman*=q;
   }
  }
  return(whichone5(ISOG,cman));
 }
 if (ISOG==-12)
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,6);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,6);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(6)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[1],stoi(3)),gzero)!=1) cman+=4;
    if (gequal(gmod((GEN) linvec[1],gdeux),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[2],gdeux),gzero)!=1)
      cman+=bentisog((GEN) CL[2],(GEN) CL[3],(GEN) CL[4])-1;
     else cman+=realisog2((GEN) CL[2],(GEN) CL[3],
			  (GEN) CL[4],gmul(realp,stoi(6)))-1;
    }
   }
   if (gequal(gmod((GEN) linvec[2],stoi(6)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[2],stoi(3)),gzero)!=1) cman+=4;
    if (gequal(gmod((GEN) linvec[2],gdeux),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[1],gdeux),gzero)==1)
      cman+=imagisog2((GEN) CL[2],(GEN) CL[3],
		      (GEN) CL[4],gmul(imagp,stoi(6)))-1;
    }
   }
  }
  return(cman);
 }
 if (ISOG==-18)
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,6);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,6);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(6)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[1],gdeux),gzero)!=1) cman+=3;
    if (gequal(gmod((GEN) linvec[1],stoi(3)),gzero)!=1)
      cman+=realisog((GEN) CL[2],(GEN) CL[3],gmul(realp,stoi(6)))-1;
   }
   if (gequal(gmod((GEN) linvec[2],stoi(6)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[2],gdeux),gzero)!=1) cman+=3;
    if (gequal(gmod((GEN) linvec[2],stoi(3)),gzero)!=1)
      cman+=imagisog((GEN) CL[2],(GEN) CL[3],gmul(imagp,stoi(6)))-1;
   }
  }
  return(cman);
 }
 printf("ERROR in x0isog5\n"); return(0);
}
