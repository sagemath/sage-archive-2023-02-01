#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void do3isog(GEN *CURV1,GEN *CURV2,GEN *CURV3)
{GEN TEMP,root1,root2,A; int nroots;

 nroots=findthree(*CURV1,&root1,&root2);
 if (nroots==1)
 {ISOG=3;
  *CURV2=threegetcurve(*CURV1,root1);
  cv(2,*CURV2); nroots=findthree(*CURV2,&root1,&root2);
  if (nroots==2)
  {ISOG=9;
   *CURV3=threegetcurve(*CURV2,root1);
   if ((gequal((GEN) (*CURV3)[11],(GEN) (*CURV1)[11])==1) &&
       (gequal((GEN) (*CURV3)[10],(GEN) (*CURV1)[10])==1))
     *CURV3=threegetcurve(*CURV2,root2);
   cv(3,*CURV3); A=myvol(*CURV2);
   if (mpcmp(myvol(*CURV1),A)==-1)
   {if (mpcmp(A,myvol(*CURV3))==-1)
    {PLACE=3; cv(1,*CURV3); cv(3,*CURV1);
     TEMP=*CURV1; *CURV1=*CURV3; *CURV3=TEMP; return;
    }
    else
    {ISOG=-9; PLACE=2; cv(1,*CURV2); cv(2,*CURV1);
     TEMP=*CURV1; *CURV1=*CURV2; *CURV2=TEMP; return;
    }
   }
  }
  else
  {if (mpcmp(myvol(*CURV1),myvol(*CURV2))==-1)
   {PLACE=2; cv(2,*CURV1); cv(1,*CURV2);
    TEMP=*CURV1; *CURV1=*CURV2; *CURV2=TEMP;
   }
   return;
  }
 }
 else if (nroots==2)
 {ISOG=9; *CURV2=threegetcurve(*CURV1,root1); cv(2,*CURV2);
  *CURV3=threegetcurve(*CURV1,root2); cv(3,*CURV3);

  A=myvol(*CURV1);
  if (mpcmp(A,myvol(*CURV2))==-1)
  {PLACE=2; cv(1,*CURV2); cv(2,*CURV1);
   TEMP=*CURV1; *CURV1=*CURV2; (*CURV2)=TEMP;
  }
  else if (mpcmp(A,myvol(*CURV3))==-1)
  {PLACE=2; cv(1,*CURV3); cv(2,*CURV1); cv(3,*CURV2);
   TEMP=*CURV1; *CURV1=*CURV3; *CURV3=*CURV2; *CURV2=TEMP;
  }
  else {ISOG=-9;}
 }
}
