#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int x0array3[4][4]={{0,1,0,1},{1,-1,3,-1},{0,3,4,3},{1,-1,3,-1}};
int x0array3B[4][4]={{0,3,0,3},{3,-1,1,-1},{0,1,4,1},{3,-1,1,-1}};
int x0array32[4][4]={{-1,4,3,4},{7,0,8,0},{6,5,0,5},{7,0,8,0}};

int x0isog3(GEN zn,int p)
{GEN listofd,TEMP,getperiods,imagp,realp,linvec,CURV1; int i,rm,im,tip;

 //8,+-16,32,-27
 listofd=getd1(zn,p,stoi(p));
 listofd=geval(gtoset(concat(listofd,getd1(zn,p*p,stoi(p)))));
 listofd=geval(gtoset(concat(listofd,getd1(zn,p*p*p,stoi(p)))));
 if ((ISOG==8) || (ISOG==16))
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  if (gsigne((GEN) CURV1[12])==-1)
  {CURV1=ellinit0((GEN) CL[4],1,ELLACC);
   getperiods=periodvolvec(CURV1);
   realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
   shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,2);
   shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,2);
   for (i=1;i<=glength(listofd);i++)
   {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
    linvec=linearcomb(realp,imagp,TEMP);
    //output(linvec);
    rm=itos((GEN) linvec[1]); if (rm<0) rm=-rm;
    im=itos((GEN) linvec[2]); if (im<0) im=-im;
    tip=x0array3B[rm&3][im&3];
    switch(tip)
    {case 1: return(4+whichbent((GEN) CL[6],(GEN) CL[5]));
     case -1: return(4+whichbent((GEN) CL[5],(GEN) CL[6]));
     case 3: return(3); case 4: return(4);
     case 0: continue;
    }
    TEMP=(GEN) periodvolvec(ellinit0((GEN) CL[2],1,ELLACC))[2];
    if (gequal(ground(gdiv(TEMP,gmul2n(realp,1))),gun)==1)
    {if (im&2) return(4); if (!(rm&2)) return(2);}
    else
    {if (rm&2) return(4); if (!(im&2)) return(2);}
   }
   return(1);
  }
  else
  {CURV1=ellinit0((GEN) CL[4],1,ELLACC);
   getperiods=periodvolvec(CURV1);
   realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
   shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,2);
   shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,2);
   for (i=1;i<=glength(listofd);i++)
   {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
    linvec=linearcomb(realp,imagp,TEMP);
    //output(linvec);
    rm=itos((GEN) linvec[1]); if (rm<0) rm=-rm;
    im=itos((GEN) linvec[2]); if (im<0) im=-im;
    tip=x0array3[rm&3][im&3];
    switch(tip)
    {case 1: return(4+whichbent((GEN) CL[6],(GEN) CL[5]));
     case -1: return(4+whichbent((GEN) CL[5],(GEN) CL[6]));
     case 3: return(3); case 4: return(4);
     case 0: continue;
    }
    TEMP=(GEN) periodvolvec(ellinit0((GEN) CL[2],1,ELLACC))[2];
    if (gequal(ground(gdiv(TEMP,gmul2n(realp,1))),gun)==1)
    {if (im&2) return(4); if (im&4) return(2);}
    else
    {if (rm&2) return(4); if (rm&4) return(2);}
   }
   return(1);
  }
 }
 if (ISOG==-16) {printf("Not implemented yet\n"); return(0);}
 if (ISOG==32)
 {CURV1=ellinit0((GEN) CL[2],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,4);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,4);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   rm=itos((GEN) linvec[1]); if (rm<0) rm=-rm;
   im=itos((GEN) linvec[2]); if (im<0) im=-im;
   tip=x0array3[rm&3][im&3];
   TEMP=(GEN) periodvolvec(ellinit0((GEN) CL[3],1,ELLACC))[2];
   if (gequal(ground(gdiv(gmul2n(realp,1),TEMP)),gun)!=1)
   {if (tip>=3) tip+=3; if (tip>=8) tip-=6;}
   switch(tip)
   {case 0: {printf("ERROR ISOG32\n"); exit(1);}
    case 3: return(3); case 6: return(6);
    case 4: return(3+whichbent((GEN) CL[5],(GEN) CL[4]));
    case 5: return(3+whichbent((GEN) CL[4],(GEN) CL[5]));
    case 7: return(6+whichbent((GEN) CL[8],(GEN) CL[7]));
    case 8: return(6+whichbent((GEN) CL[7],(GEN) CL[8]));
    case -1: if ((rm+im)&7) return(2);
   }
  }
  return(1);
 }
 if (ISOG==-27)
 {CURV1=ellinit0((GEN) CL[1],1,ELLACC);
  getperiods=periodvolvec(CURV1);
  realp=(GEN) getperiods[2]; imagp=(GEN) getperiods[3];
  shiftperiodi(gsigne((GEN) CURV1[12]),&imagp,realp,27);
  shiftperiodr(gsigne((GEN) CURV1[12]),&imagp,&realp,27);
  for (i=1;i<=glength(listofd);i++)
  {TEMP=computeperiod((GEN) CL[1],(GEN) listofd[i],minp(realp,imagp));
   linvec=linearcomb(realp,imagp,TEMP);
   //output(linvec);
   if (gequal(gmod((GEN) linvec[1],stoi(27)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[1],stoi(9)),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[1],stoi(3)),gzero)!=1)
      return(4); else return(3);
    }
    else return(2);
   }
   if (gequal(gmod((GEN) linvec[2],stoi(27)),gzero)!=1)
   {if (gequal(gmod((GEN) linvec[2],stoi(9)),gzero)!=1)
    {if (gequal(gmod((GEN) linvec[2],stoi(3)),gzero)!=1)
      return(4); else return(3);
    }
    else return(2);
   }
  }
  return(1);
 }
 printf("Error in x0isog3\n"); return(0);
}
