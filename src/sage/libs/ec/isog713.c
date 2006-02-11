#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void do7isog(GEN CURV1)
{int a=0; GEN CURV2;
 a=doisogimag(CURV1,&CURV2,7); if (a==0) a=doisogreal(CURV1,&CURV2,7);
 if (a==0) return; ISOG=7;
 if (mpcmp(myvol(CURV2),myvol(CURV1))==1)
 {cv(2,CURV1); cv(1,CURV2); PLACE=2;}
 else
 {cv(1,CURV1); cv(2,CURV2); PLACE=1;}
}

void do13isog(GEN CURV1)
{int a=0; GEN CURV2;
 a=doisogimag(CURV1,&CURV2,13); if (a==0) a=doisogreal(CURV1,&CURV2,13);
 if (a==0) return; ISOG=13;
 if (mpcmp(myvol(CURV2),myvol(CURV1))==1)
 {cv(2,CURV1); cv(1,CURV2); PLACE=2;}
 else
 {cv(1,CURV1); cv(2,CURV2); PLACE=1;}
}
