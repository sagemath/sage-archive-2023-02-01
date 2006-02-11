#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void do4twoisog(GEN CURV1,GEN root1,GEN root2,GEN root3)
{GEN CURV2,CURV3,CURV4,CURV5,CURV6,CURV7,CURV8,TEMP; int nroots;

 CURV2=twogetcurve(CURV1,root1); CURV3=twogetcurve(CURV1,root2);
 CURV4=twogetcurve(CURV1,root3); nroots=findtwo(CURV2,&root1,&root2,&root3);
 if (nroots==3) {TEMP=gcopy(CURV2); CURV2=gcopy(CURV4); CURV4=gcopy(TEMP);}
 else
 {nroots=findtwo(CURV3,&root1,&root2,&root3);
  if (nroots==3) {TEMP=gcopy(CURV3); CURV3=gcopy(CURV4); CURV4=gcopy(TEMP);}
  else nroots=findtwo(CURV4,&root1,&root2,&root3);
 }
 if (nroots==1)
 {ISOG=4; PLACE=2;
  if (mpcmp(myvol(CURV2),myvol(CURV1))==1)
  {cv(2,CURV1); cv(1,CURV2); cv(3,CURV3); cv(4,CURV4); return;}
  if (mpcmp(myvol(CURV3),myvol(CURV1))==1)
  {cv(2,CURV1); cv(3,CURV2); cv(1,CURV3); cv(4,CURV4); return;}
  if (mpcmp(myvol(CURV4),myvol(CURV1))==1)
  {cv(2,CURV1); cv(4,CURV2); cv(3,CURV3); cv(1,CURV4); return;}
  else
  {PLACE=1; ISOG=-4; cv(1,CURV1); cv(2,CURV2);
   cv(3,CURV3); cv(4,CURV4); return;
  }
 }

 CURV5=twogetcurve(CURV4,root1); CURV6=twogetcurve(CURV4,root2);
 if ((gequal((GEN) CURV5[11],(GEN) CURV1[11])==1) &&
     (gequal((GEN) CURV5[10],(GEN) CURV1[10])==1))
   CURV5=twogetcurve(CURV4,root3);
 if ((gequal((GEN) CURV6[11],(GEN) CURV1[11])==1) &&
     (gequal((GEN) CURV6[10],(GEN) CURV1[10])==1))
   CURV6=twogetcurve(CURV4,root3);
 nroots=findtwo(CURV2,&root1,&root2,&root3);
 if (nroots==3) {TEMP=gcopy(CURV2); CURV2=gcopy(CURV3); CURV3=gcopy(TEMP);}
 else nroots=findtwo(CURV3,&root1,&root2,&root3);

 if (nroots==3)
 {CURV7=twogetcurve(CURV3,root1); CURV8=twogetcurve(CURV3,root2);
  if ((gequal((GEN) CURV7[11],(GEN) CURV1[11])==1) &&
      (gequal((GEN) CURV7[10],(GEN) CURV1[10])==1))
    CURV7=twogetcurve(CURV3,root3);
  if ((gequal((GEN) CURV8[11],(GEN) CURV1[11])==1) &&
      (gequal((GEN) CURV8[10],(GEN) CURV1[10])==1))
    CURV8=twogetcurve(CURV3,root3);
  sortsixteen(CURV1,CURV2,CURV3,CURV4,CURV7,CURV8,CURV5,CURV6,1); return;
 }

 nroots=findtwo(CURV5,&root1,&root2,&root3);
 if (nroots==3) {TEMP=gcopy(CURV5); CURV5=gcopy(CURV6); CURV6=gcopy(TEMP);}
 else nroots=findtwo(CURV6,&root1,&root2,&root3);
 if (nroots==3)
 {CURV7=twogetcurve(CURV6,root1); CURV8=twogetcurve(CURV6,root2);
  if ((gequal((GEN) CURV7[11],(GEN) CURV4[11])==1) &&
      (gequal((GEN) CURV7[10],(GEN) CURV4[10])==1))
    CURV7=twogetcurve(CURV6,root3);
  if ((gequal((GEN) CURV8[11],(GEN) CURV4[11])==1) &&
      (gequal((GEN) CURV8[10],(GEN) CURV4[10])==1))
    CURV8=twogetcurve(CURV6,root3);
  sortsixteen(CURV4,CURV5,CURV1,CURV6,CURV2,CURV3,CURV7,CURV8,3); return;
 }

 ISOG=8; sorteight(CURV1,CURV2,CURV3,CURV4,CURV5,CURV6,1); return;
}


