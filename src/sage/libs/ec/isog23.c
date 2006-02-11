#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void do3twoisog(GEN CURV1,GEN CURV2)
{GEN CURV3,CURV4,CURV5,CURV6,CURV7,CURV8;
 int nroots; GEN root1,root2,root3; int i=0;

 nroots=findtwo(CURV1,&root1,&root2,&root3); if (nroots==0) return;
 if (nroots==1)
 {CURV3=twogetcurve(CURV1,root1); nroots=findtwo(CURV3,&root1,&root2,&root3);
  if (nroots==3)
  {CURV4=twogetcurve(CURV3,root1); CURV5=twogetcurve(CURV3,root2);
   if ((gequal((GEN) CURV4[11],(GEN) CURV1[11])==1) &&
       (gequal((GEN) CURV4[10],(GEN) CURV1[10])==1))
     CURV4=twogetcurve(CURV3,root3);
   if ((gequal((GEN) CURV5[11],(GEN) CURV1[11])==1) &&
       (gequal((GEN) CURV5[10],(GEN) CURV1[10])==1))
     CURV5=twogetcurve(CURV3,root3);
   nroots=findthree(CURV3,&root1,&root2); CURV6=threegetcurve(CURV3,root1);
   nroots=findthree(CURV4,&root1,&root2); CURV7=threegetcurve(CURV4,root1);
   nroots=findthree(CURV5,&root1,&root2); CURV8=threegetcurve(CURV5,root1);
   if (PLACE==2) i=1;
   sorttwelve(CURV3,CURV1,CURV4,CURV5,CURV6,CURV2,CURV7,CURV8,2);
   if (i==1) PLACE+=4;
   return;
  }
  ISOG=6; if (PLACE==2) PLACE=3;
  nroots=findtwo(CURV2,&root1,&root2,&root3);  CURV4=twogetcurve(CURV2,root1);
  if (mpcmp(myvol(CURV1),myvol(CURV3))==1)
  {cv(1,CURV1); cv(2,CURV3); cv(3,CURV2); cv(4,CURV4);}
  else {PLACE++; cv(2,CURV1); cv(1,CURV3); cv(4,CURV2); cv(3,CURV4);}
  return;
 }

 CURV3=twogetcurve(CURV1,root1); CURV4=twogetcurve(CURV1,root2);
 CURV5=twogetcurve(CURV1,root3); nroots=findthree(CURV3,&root1,&root2);
 CURV6=threegetcurve(CURV3,root1); nroots=findthree(CURV4,&root1,&root2);
 CURV7=threegetcurve(CURV4,root1); nroots=findthree(CURV5,&root1,&root2);
 CURV8=threegetcurve(CURV5,root1);

 if (PLACE==2) i=1;
 sorttwelve(CURV1,CURV3,CURV4,CURV5,CURV2,CURV6,CURV7,CURV8,1);
 if (i==1) PLACE+=4;
 return;
}

void do9twoisog(GEN CURV1,GEN CURV2,GEN CURV3)
{GEN CURV4,CURV5,CURV6; int nroots; GEN root1,root2,root3;

 nroots=findtwo(CURV1,&root1,&root2,&root3); if (nroots==0) return;
 ISOG*=2;
 CURV4=twogetcurve(CURV1,root1); nroots=findtwo(CURV2,&root1,&root2,&root3);
 CURV5=twogetcurve(CURV2,root1); nroots=findtwo(CURV3,&root1,&root2,&root3);
 CURV6=twogetcurve(CURV3,root1);
 if (mpcmp(myvol(CURV4),myvol(CURV1))==1)
 {cv(1,CURV4); cv(2,CURV5); cv(3,CURV6);
  cv(4,CURV1); cv(5,CURV2); cv(6,CURV3); PLACE+=3;
 }
 else {cv(4,CURV4); cv(5,CURV5); cv(6,CURV6);}
}
