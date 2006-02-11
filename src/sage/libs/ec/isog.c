#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void isogenies()
{GEN CURV,CURV2,CURV3; GEN vect=gzero;
 int ptvec[30]; int OK=1; int i,j;
 pari_sp memspot=avma;
 int myprime[26]=
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};

 ISOG=1; PLACE=1; CURVES=cgetg(9,t_VEC);
 for(i=1;i<=8;i++)
 {CURVES[i]=lgetg(6,t_VEC);
  for(j=1;j<=5;j++) ((GEN) CURVES[i])[j]=(long) gzero;
 }
 cv(1,CURVE);
 if (ISPRIME)
 {if (gequal(TWCOND,stoi(11))) isogN11(); if (gequal(TWCOND,stoi(17))) isogN17();
  if (gequal(TWCOND,stoi(19))) isogN19(); if (gequal(TWCOND,stoi(37))) isogN37();
  if (ISOG!=1) return;
  if (carreparfait(gsub(TWCOND,gmul2n(gun,6)))) isogNsz(); return;
 }

 CURV=gcopy(CURVE); CURV2=gcopy(CURVE); CURV3=gcopy(CURVE);

 if (!ISSQFREE) {checkexotic(); if (ISOG!=1) return;}
 else
 {vect=cgetg(27,t_VEC);
  for(i=1;i<=25;i++) {vect[i]=(long) stoi(myprime[i]);}
  memspot=avma;
  for(i=1;i<=25;i++)
    ptvec[i]=myprime[i]+1-itos(ellap0(TWCURVE,(GEN) vect[i],0));
  for(i=3;i<=25;i++)
  {if (ggval(TWCOND,(GEN) vect[i])==0) if ((ptvec[i]%3)!=0) OK=0;}
 }

 if (OK) do3isog(&CURV,&CURV2,&CURV3);
 garbagecollect(&CURV,&CURV2,&CURV3,memspot);

 OK=1;
 if (ISSQFREE)
 {for(i=2;i<=25;i++)
  {if (ggval(TWCOND,(GEN) vect[i])==0) if ((ptvec[i]%2)!=0) OK=0;}
 }

 if (OK)
 {if ((ISOG==9) || (ISOG==-9)) do9twoisog(CURV,CURV2,CURV3);
  else if (ISOG==3) do3twoisog(CURV,CURV2); else dotwoisog(CURV,&CURV2);
  garbagecollect(&CURV,&CURV2,&CURV3,memspot);
 }
 if (ISOG==9) return; if (ISOG==-9) return; if (ISOG==3) return;
 if (ISOG==18) return; if (ISOG==-18) return;
 if (ISOG==12) return; if (ISOG==-12) return;
 if (ISOG==4) return; if (ISOG==-4) return;
 if (ISOG==8) return; if (ISOG==16) return;

 OK=1;
 if (ISSQFREE)
 {for(i=4;i<=25;i++)
  {if (ggval(TWCOND,(GEN) vect[i])==0) if ((ptvec[i]%5)!=0) OK=0;}
 }
 if (OK)
 {if (ISOG==2) do2fiveisog(CURV,CURV2); else do5isog(CURV);
  garbagecollect(&CURV,&CURV2,&CURV3,memspot);
 }

 if (ISOG==2)
 {if (mpcmp(myvol(CURV2),myvol(CURV))==1)
  {cv(2,CURV); cv(1,CURV2); PLACE=2; return;}
  else {cv(1,CURV); cv(2,CURV2); PLACE=1; return;}
 }
 if (ISOG!=1) return; OK=1;
 if (ISSQFREE)
 {for(i=5;i<=25;i++)
  {if (ggval(TWCOND,(GEN) vect[i])==0) if ((ptvec[i]%7)!=0) OK=0;}
 }
 if (OK) {do7isog(CURV); garbagecollect(&CURV,&CURV2,&CURV3,memspot);}
 if (ISOG==7) return;
 if (!ISSQFREE) {do13isog(CURV); garbagecollect(&CURV,&CURV2,&CURV3,memspot);}
 return;
}

void garbagecollect(GEN *A,GEN *B,GEN *C,pari_sp memspot)
     //     {return;}

{GEN TEMPA; GEN TEMPB; GEN TEMPC; GEN TEMPE;
 int i; int j;

 avma-=2048;
 TEMPA=cgetg(14,t_VEC); TEMPB=cgetg(14,t_VEC);
 TEMPC=cgetg(14,t_VEC); TEMPE=cgetg(9,t_VEC);
 for (i=1;i<=8;i++) TEMPE[i]=(long) cgetg(6,t_VEC);
 for (i=1;i<=13;i++) TEMPA[i]=lcopy((GEN) (*A)[i]);
 for (i=1;i<=13;i++) TEMPB[i]=lcopy((GEN) (*B)[i]);
 for (i=1;i<=13;i++) TEMPC[i]=lcopy((GEN) (*C)[i]);
 for (i=1;i<=8;i++)
 {for (j=1;j<=5;j++) {((GEN) TEMPE[i])[j]=lcopy((GEN) ((GEN) CURVES[i])[j]);}}
 avma=memspot;
 *A=cgetg(14,t_VEC); *B=cgetg(14,t_VEC);
 *C=cgetg(14,t_VEC); CURVES=cgetg(9,t_VEC);
 for (i=1;i<=8;i++) CURVES[i]=(long) cgetg(6,t_VEC);
 for (i=1;i<=13;i++) (*A)[i]=lcopy((GEN) (TEMPA)[i]);
 for (i=1;i<=13;i++) (*B)[i]=lcopy((GEN) (TEMPB)[i]);
 for (i=1;i<=13;i++) (*C)[i]=lcopy((GEN) (TEMPC)[i]);
 for (i=1;i<=8;i++)
 {for (j=1;j<=5;j++) {((GEN) CURVES[i])[j]=lcopy((GEN) ((GEN) TEMPE[i])[j]);}}
}

