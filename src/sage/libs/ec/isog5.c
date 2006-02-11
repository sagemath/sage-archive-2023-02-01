#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void do5isog(GEN CURV1)
{int a=0; int j=0; GEN CURV2; GEN CURV3; GEN A;

 a=doisogimag(CURV1,&CURV2,5); j=doisogreal(CURV1,&CURV3,5);
 if ((a==0) && (j==0)) return;
 if ((a==1) && (j==1))
 {ISOG=25; PLACE=2;  A=myvol(CURV1);
  if (mpcmp(myvol(CURV2),A)==1)
  {cv(2,CURV1); cv(1,CURV2); cv(3,CURV3); return;}
  if (mpcmp(myvol(CURV3),A)==1)
  {cv(2,CURV1); cv(3,CURV2); cv(1,CURV3); return;}
  else {ISOG=-25; PLACE=1; cv(1,CURV1); cv(2,CURV2); cv(3,CURV3); return;}
 }
 if (a==1)
 {j=doisogimag(CURV2,&CURV3,5);
  if (j==1)
  {ISOG=25; A=myvol(CURV2);
   if (mpcmp(myvol(CURV1),A)==1)
   {PLACE=1; cv(1,CURV1); cv(2,CURV2); cv(3,CURV3); return;}
   if (mpcmp(myvol(CURV3),A)==1)
   {PLACE=3; cv(3,CURV1); cv(2,CURV2); cv(1,CURV3); return;}
   else {ISOG=-25; PLACE=2; cv(2,CURV1); cv(1,CURV2); cv(3,CURV3); return;}
  }
  else
  {ISOG=5;
   if (mpcmp(myvol(CURV2),myvol(CURV1))==1)
   {cv(2,CURV1); cv(1,CURV2); PLACE=2; return;}
   else {cv(1,CURV1); cv(2,CURV2); PLACE=1; return;}
  }
 }
 else
 {a=doisogreal(CURV3,&CURV2,5);
  if (a==1)
  {ISOG=25; A=myvol(CURV3);
   if (mpcmp(myvol(CURV1),A)==1)
   {PLACE=1; cv(1,CURV1); cv(2,CURV3); cv(3,CURV2); return;}
   if (mpcmp(myvol(CURV2),A)==1)
   {PLACE=3; cv(3,CURV1); cv(2,CURV3); cv(1,CURV2); return;}
   else {ISOG=-25; PLACE=2; cv(2,CURV1); cv(1,CURV3); cv(3,CURV2); return;}
  }
  else
  {ISOG=5;
   if (mpcmp(myvol(CURV3),myvol(CURV1))==1)
   {cv(2,CURV1); cv(1,CURV3); PLACE=2; return;}
   else
   {cv(1,CURV1); cv(2,CURV3); PLACE=1; return;}
  }
 }
}
