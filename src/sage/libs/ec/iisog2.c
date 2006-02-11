#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int findtwo(GEN CURV,GEN *root1,GEN *root2,GEN *root3)
{GEN poly,thepoly,TEMP,TEMP2; int nroots=0; int i,number;

 poly=cgetg(5,t_VEC); *root1=gzero; *root2=gzero; *root3=gzero;
 poly[1]=(long) gun; poly[2]=(long) CURV[6];
 poly[3]=(long) gmul2n((GEN) CURV[7],3);
 poly[4]=(long) gmul2n((GEN) CURV[8],4);
 thepoly=gtopoly(poly,fetch_var()); TEMP=factor(thepoly);
 number=itos((GEN) matsize(TEMP)[1]);
 for(i=1;i<=number;i++)
 {TEMP2=(GEN) ((GEN) TEMP[1])[i];
  if(degree(TEMP2)==1)
  {nroots++; if (nroots==1) *root1=gneg(gmul2n((polcoeff0(TEMP2,0,-1)),-2));
   if (nroots==2) *root2=gneg(gmul2n((polcoeff0(TEMP2,0,-1)),-2));
   if (nroots==3) *root3=gneg(gmul2n((polcoeff0(TEMP2,0,-1)),-2));
  }
 }
 delete_var(); return(nroots);
}

GEN twogetcurve(GEN CURV,GEN root)
{GEN thevec=cgetg(6,t_VEC); GEN T,W,newcurve; long TH,TL; TH=avma;
 thevec[1]=(long) CURV[1]; thevec[2]=(long) CURV[2]; thevec[3]=(long) CURV[3];
 T=gadd(gadd(gmul(stoi(3),gsqr(root)),gmul(gmul2n(root,1),(GEN) CURV[2])),
	gadd((GEN) CURV[4],gmul(gadd(gmul((GEN) CURV[1],root),(GEN) CURV[3]),
				gmul2n((GEN) CURV[1],-1))));
 W=gadd(gadd(gmul2n(gmul(root,gsqr(root)),2),
	     gmul((GEN) CURV[6],gsqr(root))),
	gadd(gmul2n(gmul(root,(GEN) CURV[7]),1),
	     gadd((GEN) CURV[8],gmul(root,T))));
 thevec[4]=(long) gadd((GEN) CURV[4],gmul(stoi(-5),T));
 thevec[5]=(long) gsub(gadd((GEN) CURV[5],gmul(stoi(-7),W)),
		       gmul((GEN) CURV[6],T));
 newcurve=ellinit0(thevec,1,ELLACC);
 newcurve=mineqfromc4c6((GEN) newcurve[10],(GEN) newcurve[11]);
 TL=avma; return(gerepile(TH,TL,gcopy(newcurve)));
}

