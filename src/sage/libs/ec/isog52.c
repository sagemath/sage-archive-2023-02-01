#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void do2fiveisog(GEN CURV1,GEN CURV2)
{int a=0; int nroots; GEN A,B,root1,root2,root3,CURV3,CURV4;

 a=doisogimag(CURV1,&CURV3,5); if (a==0) a=doisogreal(CURV1,&CURV3,5);
 if (a==0) return;
 ISOG=10; nroots=findtwo(CURV3,&root1,&root2,&root3);
 CURV4=twogetcurve(CURV3,root1); A=myvol(CURV1); B=myvol(CURV2);
 if (mpcmp(A,B)==1)
 {if(mpcmp(A,myvol(CURV3))==1)
  {PLACE=1; cv(1,CURV1); cv(2,CURV2); cv(3,CURV3); cv(4,CURV4); return;}
  else
  {PLACE=3; cv(3,CURV1); cv(4,CURV2); cv(1,CURV3); cv(2,CURV4); return;}
 }
 else
 {if(mpcmp(B,myvol(CURV4))==1)
  {PLACE=2; cv(2,CURV1); cv(1,CURV2); cv(4,CURV3); cv(3,CURV4); return;}
  else
  {PLACE=4; cv(4,CURV1); cv(3,CURV2); cv(2,CURV3); cv(1,CURV4); return;}
 }
}
