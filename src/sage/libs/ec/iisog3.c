#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int findthree(GEN CURV,GEN *root1,GEN *root2)
{GEN A,B,TEMP,TEMP2,poly,thepoly; int nroots=0; int i,number;

 *root1=gzero; *root2=gzero; poly=cgetg(6,t_VEC);
 A=gmul((GEN) CURV[10],stoi(-27)); B=gmul((GEN) CURV[11],stoi(-54));
 poly[1]=(long) gun; poly[2]=(long) gzero;
 poly[3]=(long) gmul2n(A,1); poly[4]=(long) gmul2n(B,2);
 poly[5]=(long) gdiv(gsqr(A),stoi(-3));
 thepoly=gtopoly(poly,fetch_var()); TEMP=factor(thepoly);
 number=itos((GEN) matsize(TEMP)[1]);
 for(i=1;i<=number;i++)
 {TEMP2=(GEN) ((GEN) TEMP[1])[i];
  if(degree(TEMP2)==1)
  {nroots++; if (nroots==1) *root1=gneg(polcoeff0(TEMP2,0,-1));
   else *root2=gneg(polcoeff0(TEMP2,0,-1));
  }
 }
 delete_var(); return(nroots);
}

GEN threegetcurve(GEN CURV, GEN root)
{GEN thevec=cgetg(6,t_VEC); GEN A,B,newcurve; long TH,TL; TH=avma;
 thevec[1]=(long) gzero; thevec[2]=(long) gzero; thevec[3]=(long) gzero;
 A=gmul(stoi(-3),gadd(gmul(stoi(-81),(GEN) CURV[10]),
		      gmul(stoi(10),gsqr(root))));
 B=gadd(gmul(stoi(-70),gmul(root,gsqr(root))),
	gadd(gmul(gmul((GEN) CURV[10],stoi(1134)),root),
	     gmul(stoi(1458),(GEN) CURV[11])));
 thevec[4]=(long) A; thevec[5]=(long) B;
 newcurve=ellinit0(thevec,1,ELLACC);
 newcurve=mineqfromc4c6((GEN) newcurve[10],(GEN) newcurve[11]);
 TL=avma; return(gerepile(TH,TL,gcopy(newcurve)));
}

