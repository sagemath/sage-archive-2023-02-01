#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int whichbent(GEN C1,GEN C2)
{GEN TEMP;
 TEMP=ellinit0(C1,1,ELLACC);
 if (gsigne((GEN) TEMP[12])==-1) {cgiv(TEMP); return(1);}
 TEMP=ellinit0(C2,1,ELLACC);
 if (gsigne((GEN) TEMP[12])==-1) {cgiv(TEMP); return(2);}
 printf("Error in whichbent\n"); return(0);
}
int bentisog(GEN C2,GEN C3,GEN C4)
{GEN TEMP;
 TEMP=ellinit0(C2,1,ELLACC);
 if (gsigne((GEN) TEMP[12])==-1) {cgiv(TEMP); return(2);}
 TEMP=ellinit0(C3,1,ELLACC);
 if (gsigne((GEN) TEMP[12])==-1) {cgiv(TEMP); return(3);}
 TEMP=ellinit0(C4,1,ELLACC);
 if (gsigne((GEN) TEMP[12])==-1) {cgiv(TEMP); return(4);}
 printf("Error in bentisog\n"); return(0);
}

int realisog(GEN C2,GEN C3,GEN rp)
{GEN TEMP,comp;
 TEMP=ellinit0(C2,1,ELLACC); comp=(GEN) periodvolvec(TEMP)[2];
 if (gequal(ground(gdiv(comp,rp)),gun)==1) return(2);
 TEMP=ellinit0(C3,1,ELLACC); comp=(GEN) periodvolvec(TEMP)[2];
 if (gequal(ground(gdiv(comp,rp)),gun)==1) return(3);
 printf("Error in realisog\n"); return(0);
}

int imagisog(GEN C2,GEN C3,GEN ip)
{GEN TEMP,comp;
 TEMP=ellinit0(C2,1,ELLACC); comp=gimag((GEN) periodvolvec(TEMP)[3]);
 if (gequal(ground(gdiv(comp,gimag(ip))),gun)==1) return(2);
 TEMP=ellinit0(C3,1,ELLACC); comp=gimag((GEN) periodvolvec(TEMP)[3]);
 if (gequal(ground(gdiv(comp,gimag(ip))),gun)==1) return(3);
 printf("Error in imagisog\n"); return(0);
}

int realisog2(GEN C2,GEN C3,GEN C4,GEN rp)
{GEN TEMP,comp;
 TEMP=ellinit0(C2,1,ELLACC); comp=(GEN) periodvolvec(TEMP)[2];
 if (gequal(ground(gdiv(comp,rp)),gun)==1) return(2);
 TEMP=ellinit0(C3,1,ELLACC); comp=(GEN) periodvolvec(TEMP)[2];
 if (gequal(ground(gdiv(comp,rp)),gun)==1) return(3);
 TEMP=ellinit0(C4,1,ELLACC); comp=(GEN) periodvolvec(TEMP)[2];
 if (gequal(ground(gdiv(comp,rp)),gun)==1) return(4);
 printf("Error in realisog2\n"); return(0);
}

int imagisog2(GEN C2,GEN C3,GEN C4,GEN ip)
{GEN TEMP,comp;
 TEMP=ellinit0(C2,1,ELLACC); comp=gimag((GEN) periodvolvec(TEMP)[3]);
 if (gsigne((GEN) TEMP[12])==1)
   if (gequal(ground(gdiv(comp,gimag(ip))),gun)==1) return(2);
 TEMP=ellinit0(C3,1,ELLACC); comp=gimag((GEN) periodvolvec(TEMP)[3]);
 if (gsigne((GEN) TEMP[12])==1)
   if (gequal(ground(gdiv(comp,gimag(ip))),gun)==1) return(3);
 TEMP=ellinit0(C4,1,ELLACC); comp=gimag((GEN) periodvolvec(TEMP)[3]);
 if (gsigne((GEN) TEMP[12])==1)
   if (gequal(ground(gdiv(comp,gimag(ip))),gun)==1) return(4);
 printf("Error in imagisog2\n"); return(0);
}

