#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void dotwoisog(GEN CURV1,GEN *RETCURV)
{GEN CURV2,CURV3,CURV4,CURV5,CURV6,CURV7,CURV8,root1,root2,root3,TEMP;
 int nroots;
 nroots=findtwo(CURV1,&root1,&root2,&root3); if (nroots==0) return;
 if (nroots==3) {do4twoisog(CURV1,root1,root2,root3); return;}

 CURV2=twogetcurve(CURV1,root1); nroots=findtwo(CURV2,&root1,&root2,&root3);
 if (nroots==0) {printf("ERROR: Factor bug in PARI\n"); nroots=1;}
 if (nroots==1)
 {*RETCURV=gcopy(CURV2); cv(1,CURV1); cv(2,CURV2); ISOG=2; return;}

 CURV3=twogetcurve(CURV2,root1); CURV4=twogetcurve(CURV2,root2);
 if ((gequal((GEN) CURV3[11],(GEN) CURV1[11])==1) &&
     (gequal((GEN) CURV3[10],(GEN) CURV1[10])==1))
   CURV3=twogetcurve(CURV2,root3);
 if ((gequal((GEN) CURV4[11],(GEN) CURV1[11])==1) &&
     (gequal((GEN) CURV4[10],(GEN) CURV1[10])==1))
   CURV4=twogetcurve(CURV2,root3);

 nroots=findtwo(CURV3,&root1,&root2,&root3);
 if (nroots==3) {TEMP=gcopy(CURV3); CURV3=gcopy(CURV4); CURV4=gcopy(TEMP);}
 else nroots=findtwo(CURV4,&root1,&root2,&root3);
 if (nroots==0) {printf("ERROR: Factor bug in PARI\n"); nroots=1;}
 if (nroots==1)
 {ISOG=4;
  if (mpcmp(myvol(CURV1),myvol(CURV2))==1)
  {PLACE=1; cv(1,CURV1); cv(2,CURV2); cv(3,CURV3); cv(4,CURV4); return;}
  if (mpcmp(myvol(CURV3),myvol(CURV2))==1)
  {PLACE=3; cv(3,CURV1); cv(2,CURV2); cv(1,CURV3); cv(4,CURV4); return;}
  if (mpcmp(myvol(CURV4),myvol(CURV2))==1)
  {PLACE=3; cv(3,CURV1); cv(2,CURV2); cv(4,CURV3); cv(1,CURV4); return;}
  else
  {PLACE=3; ISOG=-4; cv(3,CURV1); cv(1,CURV2);
   cv(2,CURV3); cv(4,CURV4); return;
  }
 }

 CURV5=twogetcurve(CURV4,root1); CURV6=twogetcurve(CURV4,root2);
 if ((gequal((GEN) CURV5[11],(GEN) CURV2[11])==1) &&
     (gequal((GEN) CURV5[10],(GEN) CURV2[10])==1))
   CURV5=twogetcurve(CURV4,root3);
 if ((gequal((GEN) CURV6[11],(GEN) CURV2[11])==1) &&
     (gequal((GEN) CURV6[10],(GEN) CURV2[10])==1))
   CURV6=twogetcurve(CURV4,root3);

 nroots=findtwo(CURV3,&root1,&root2,&root3);
 if (nroots==3)
 {CURV7=twogetcurve(CURV3,root1); CURV8=twogetcurve(CURV3,root2);
  if ((gequal((GEN) CURV7[11],(GEN) CURV2[11])==1) &&
      (gequal((GEN) CURV7[10],(GEN) CURV2[10])==1))
    CURV7=twogetcurve(CURV3,root3);
  if ((gequal((GEN) CURV8[11],(GEN) CURV2[11])==1) &&
      (gequal((GEN) CURV8[10],(GEN) CURV2[10])==1))
    CURV8=twogetcurve(CURV3,root3);
  sortsixteen(CURV2,CURV1,CURV3,CURV4,CURV7,CURV8,CURV5,CURV6,2); return;
 }

 nroots=findtwo(CURV5,&root1,&root2,&root3);
 if (nroots==3) {TEMP=gcopy(CURV6); CURV6=gcopy(CURV5); CURV5=gcopy(TEMP);}
 else nroots=findtwo(CURV6,&root1,&root2,&root3);
 if (nroots==3)
 {CURV7=twogetcurve(CURV6,root1);  CURV8=twogetcurve(CURV6,root2);
  if ((gequal((GEN) CURV7[11],(GEN) CURV4[11])==1) &&
      (gequal((GEN) CURV7[10],(GEN) CURV4[10])==1))
    CURV7=twogetcurve(CURV6,root3);
  if ((gequal((GEN) CURV8[11],(GEN) CURV4[11])==1) &&
      (gequal((GEN) CURV8[10],(GEN) CURV4[10])==1))
    CURV8=twogetcurve(CURV6,root3);
  sortsixteen(CURV4,CURV5,CURV2,CURV6,CURV1,CURV3,CURV7,CURV8,5); return;
 }
 sorteight(CURV2,CURV1,CURV3,CURV4,CURV5,CURV6,2); return;
}
