#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void sorteight(GEN A,GEN B,GEN C,GEN D,GEN E,GEN F,int w)
{GEN V1,V2; // B-A-C E-D-F A-D   w is 1 or 2

 V1=myvol(A); V2=myvol(D);
 if (mpcmp(V1,V2)==1)
 {if (mpcmp(myvol(B),V1)==1)
  {cv(2,A); cv(1,B); cv(3,C); cv(4,D); cv(5,E); cv(6,F);
   ISOG=8; PLACE=3-w; return;
  }
  else if (mpcmp(myvol(C),V1)==1)
  {cv(2,A); cv(3,B); cv(1,C); cv(4,D); cv(5,E); cv(6,F);
   ISOG=8; PLACE=w+1; return;
  }
  else
  {cv(1,A); cv(2,B); cv(3,C); cv(4,D); cv(5,E); cv(6,F);
   ISOG=-8; PLACE=w; return;
  }
 }
 else
 {if (mpcmp(myvol(E),V2)==1)
  {cv(4,A); cv(5,B); cv(6,C); cv(2,D); cv(1,E); cv(3,F);
   ISOG=8; PLACE=w+3; return;
  }
  else if (mpcmp(myvol(F),V2)==1)
  {cv(4,A); cv(5,B); cv(6,C); cv(2,D); cv(3,E); cv(1,F);
   ISOG=8; PLACE=w+3; return;
  }
  else
  {cv(4,A); cv(5,B); cv(6,C); cv(1,D); cv(2,E); cv(3,F); ISOG=-8;
   switch (w)
   {case 4: {PLACE=1; return;} case 5: {PLACE=2; return;}
    case 6: {PLACE=3; return;} default: {PLACE=w+3; return;}
   }
  }
 }
}

void sortsixteen(GEN A,GEN B,GEN C,GEN D,GEN E,GEN F,GEN G,GEN H,int w)
{GEN V1,V2; // C-A-D E-C-F G-D-H A-B   w is 1,2,3,5

 V1=myvol(A); V2=myvol(C);
 if (mpcmp(V1,V2)==1)
 {V2=myvol(D);
  if (mpcmp(V1,V2)==1)
  {if (mpcmp(V1,myvol(B))==1)
   {cv(1,A); cv(2,B); cv(3,C); cv(4,D); cv(5,E); cv(6,F); cv(7,G); cv(8,H);
    ISOG=-32; if (w==5) PLACE=4; else PLACE=w; return;
   }
   else
   {cv(2,A); cv(1,B); cv(3,C); cv(6,D); cv(4,E); cv(5,F); cv(7,G); cv(8,H);
    ISOG=32;
    switch (w)
    {case 1: {PLACE=2; return;} case 5: {PLACE=4; return;}
     case 2: {PLACE=1; return;} case 3: {PLACE=3; return;}
    }
   }
  }
  else
  {if (mpcmp(V2,myvol(G))==-1)
   {cv(4,A); cv(5,B); cv(6,C); cv(2,D); cv(8,E); cv(7,F); cv(1,G); cv(3,H);
    ISOG=16; PLACE=w+3; return;
   }
   else if (mpcmp(V2,myvol(H))==-1)
   {cv(4,A); cv(5,B); cv(6,C); cv(2,D); cv(8,E); cv(7,F); cv(3,G); cv(1,H);
    ISOG=16; PLACE=w+3; return;
   }
   else
   {cv(4,A); cv(5,B); cv(6,C); cv(1,D); cv(8,E); cv(7,F); cv(2,G); cv(3,H);
    ISOG=-16; PLACE=w+3; return;
   }
  }
 }
 else
 {if (mpcmp(V2,myvol(E))==-1)
  {cv(4,A); cv(5,B); cv(2,C); cv(6,D); cv(1,E); cv(3,F); cv(7,G); cv(8,H);
   ISOG=16;
   switch (w)
   {case 1: {PLACE=4; return;} case 5: {PLACE=1; return;}
    case 2: {PLACE=5; return;} case 3: {PLACE=2; return;}
   }
  }
  else if (mpcmp(V2,myvol(F))==-1)
  {cv(4,A); cv(5,B); cv(2,C); cv(6,D); cv(3,E); cv(1,F); cv(7,G); cv(8,H);
   ISOG=16;
   switch (w)
   {case 1: {PLACE=4; return;} case 5: {PLACE=3; return;}
    case 2: {PLACE=5; return;} case 3: {PLACE=2; return;}
   }
  }
  else
  {cv(4,A); cv(5,B); cv(1,C); cv(6,D); cv(2,E); cv(3,F); cv(7,G); cv(8,H);
   ISOG=-16;
   switch (w)
   {case 1: {PLACE=4; return;} case 5: {PLACE=2; return;}
    case 2: {PLACE=5; return;} case 3: {PLACE=1; return;}
   }
  }
 }
}

void sorttwelve(GEN A,GEN B,GEN C,GEN D,GEN E,GEN F,GEN G,GEN H,int w)
{GEN V1; // A-E 3, B-F,C-G,D-H   w is 1 or 2   A stronger than E

 V1=myvol(A);
 if (mpcmp(V1,myvol(B))==-1)
 {cv(2,A); cv(1,B); cv(3,C); cv(4,D); cv(6,E); cv(5,F); cv(7,G); cv(8,H);
  ISOG=12; PLACE=3-w; return;
 }
 if (mpcmp(V1,myvol(C))==-1)
 {cv(2,A); cv(3,B); cv(1,C); cv(4,D); cv(6,E); cv(7,F); cv(5,G); cv(8,H);
  ISOG=12; PLACE=w+1; return;
 }
 if (mpcmp(V1,myvol(D))==-1)
 {cv(2,A); cv(3,B); cv(4,C); cv(1,D); cv(6,E); cv(7,F); cv(8,G); cv(5,H);
  ISOG=12; PLACE=w+1; return;
 }
 else
 {cv(1,A); cv(2,B); cv(3,C); cv(4,D); cv(5,E); cv(6,F); cv(7,G); cv(8,H);
  ISOG=-12; PLACE=w; return;
 }
}

void swap(int i,int j)
{int A; A=CL[i]; CL[i]=CL[j]; CL[j]=A;}

int ordered(int i,int j)
{
  GEN Di,Dj; int k;
  Di = mpabs((GEN)(((GEN)CL[i])[12]));
  Dj = mpabs((GEN)(((GEN)CL[j])[12]));
  k=mpcmp(Di,Dj); if (k<0) return 1; if (k>0) return 0;
  k=mpcmp((GEN) (((GEN) CL[i])[10]), (GEN) (((GEN) CL[j])[10]));
  if (k<0) return 1; if (k>0) return 0;
  k=mpcmp((GEN) (((GEN) CL[i])[11]), (GEN) (((GEN) CL[j])[11]));
  if (k<0) return 1; if (k>0) return 0;
  return 0;
}

void sortcurves(void)
{
 if (classsize()<=2) return;
 if ((ISOG==9) || (ISOG==25) || (ISOG==14) ||
     (ISOG==15) || (ISOG==21) || (ISOG==27)) return;
 if ((ISOG==6) || (ISOG==10)) return;
 if (ISOG==-9) {if (!ordered(2,3)) swap(2,3); return;}
 if (ISOG==4) {if (!ordered(3,4)) swap(3,4); return;}
 if (ISOG==8) {if (!ordered(5,6)) swap(5,6); return;}
 if (ISOG==16) {if (!ordered(7,8)) swap(7,8); return;}
 if (ISOG==-4)
 {if (!ordered(2,3)) swap(2,3); if (!ordered(2,4)) swap(2,4);
  if (!ordered(3,4)) swap(3,4); return;
 }
 if (ISOG==18) {swap(2,4); swap(3,5); swap(3,4); return;}
 if (ISOG==12) {if (!ordered(3,4)) {swap(3,4); swap(7,8);}}
 return;
}
/*
[1,0,1,-498,4228]
[1,0,1,-1473,-16652]
[1,0,1,-578,2756]
[1,0,1,-21953,-1253644]
[1,0,1,-351233,-80149132]
[1,0,1,-4358,-109132]
[1,0,1,1922,20756]
[1,0,1,-20353,-1443724]
*/
