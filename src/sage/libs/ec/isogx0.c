#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int x0isogeny()
{GEN TEMP=gun; GEN znstruc,zninner; int try,i;

 znstruc=myznstar(); zninner=(GEN) znstruc[2];
 for (i=1;i<=glength(zninner);i++) {TEMP=glcm(TEMP,(GEN) zninner[i]);}
 try=itos(gcd0(TEMP,stoi(ISOG),0));
 if ((ISOG==32) && (try>=8)) try=8;
 if ((ISOG==-4) && (try==4)) try=2; if ((ISOG==-8) && (try==8)) try=4;
 if ((ISOG==-16) && (try==16)) try=8; if ((ISOG==-32) && (try>=4)) try=4;
 if ((ISOG==-12) && (try==4)) try=2; if ((ISOG==-12) && (try==12)) try=6;
 if ((ISOG==-18) && (try==9)) try=3; if ((ISOG==-18) && (try==18)) try=6;
 if ((ISOG==-9) && (try==9)) try=3; if ((ISOG==-25) && (try==25)) try=5;
 if ((ISOG==27) && (try==27)) try=9;

 switch(try)
 {case 1: return(1);
  case 2: case 3: case 5: case 7: case 11: case 13:
  case 17: case 19: case 37: case 43: case 67: case 163:
    return(x0isog1(znstruc,try));
  case 4: return(x0isog2(znstruc,2)); case 6: return(x0isog5(znstruc,2,3));
  case 8: return(x0isog3(znstruc,2)); case 9: return(x0isog2(znstruc,3));
  case 10: return(x0isog5(znstruc,2,5)); case 12: return(x0isog6(znstruc,2,3));
  case 14: return(x0isog5(znstruc,2,7)); case 15: return(x0isog5(znstruc,3,5));
  case 16: return(x0isog4(znstruc,2)); case 18: return(x0isog6(znstruc,3,2));
  case 21: return(x0isog5(znstruc,3,7)); case 25: return(x0isog2(znstruc,5));
  case 27: return(x0isog3(znstruc,3));
 }
 return(-1);
}

GEN myznstar()
{GEN Z; GEN TEMP; int number; int i;

 number=itos((GEN) matsize(CF)[1]); Z=cgetg(5,t_VEC); Z[1]=(long) gun;
 if (mod8(COND)!=0)
 {Z[2]=lgetg(number+1,t_VEC); //phi of modulus
  Z[3]=lgetg(number+1,t_VEC); //modulus
  Z[4]=lgetg(number+1,t_VEC); //generators
  for(i=1;i<=number;i++)
  {TEMP=(GEN) ((GEN) CF[1])[i];
   ((GEN) Z[3])[i]=(long) gpow(TEMP,(GEN) ((GEN) CF[2])[i],-1);
   ((GEN) Z[2])[i]=
     (long) gmul((GEN) ((GEN) Z[3])[i],gdiv(gsub(TEMP,gun),TEMP));
   Z[1]=(long) glcm((GEN) Z[1],(GEN) ((GEN) Z[2])[i]);
   ((GEN) Z[4])[i]=(long) gener((GEN) ((GEN) Z[3])[i]);
  }
 }
 else
 {Z[2]=lgetg(number+2,t_VEC);
  Z[3]=lgetg(number+2,t_VEC);
  Z[4]=lgetg(number+2,t_VEC);
  i=ggval(COND,gdeux);
  Z[1]=(long) gmul2n(gun,i-2);
  ((GEN) Z[2])[1]=Z[1]; ((GEN) Z[2])[2]=(long) gdeux;
  ((GEN) Z[3])[1]=(long) gmul2n(gun,i); ((GEN) Z[3])[2]=((GEN) Z[3])[1];
  ((GEN) Z[4])[1]=(long) gmodulcp(stoi(5),(GEN) ((GEN) Z[3])[1]);
  ((GEN) Z[4])[2]=(long) gmodulcp(stoi(-1),(GEN) ((GEN) Z[3])[1]);
  for(i=2;i<=number;i++)
  {TEMP=(GEN) ((GEN) CF[1])[i];
   ((GEN) Z[3])[i+1]=(long) gpow(TEMP,(GEN) ((GEN) CF[2])[i],-1);
   ((GEN) Z[2])[i+1]=
     (long) gmul((GEN) ((GEN) Z[3])[i+1],gdiv(gsub(TEMP,gun),TEMP));
   Z[1]=(long) glcm((GEN) Z[1],(GEN) ((GEN) Z[2])[i+1]);
   ((GEN) Z[4])[i+1]=(long) gener((GEN) ((GEN) Z[3])[i+1]);
  }
 }
 if (TRACE) output(Z);
 return(Z);
}
