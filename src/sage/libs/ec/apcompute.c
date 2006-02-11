#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void computeap(GEN curv,int B)
{GEN T; GEN P; int64 i; double cur=0; int64 pk=0; int64 j=0;
 byteptr dp=diffptr; ulong p=0;

 if (ansize<B)
 {anarray=realloc(anarray,u8((B+1)*sizeof(double)));
  for (i=ansize+1;i<=B;i++) anarray[i]=1.0;
  p+=*dp++;
  while(p<=B)
  {if (p<=ansize) cur=anarray[p];
   else
   {P=stoi(p); T=ellap0(curv,P,0); cur=(double) itos(T); cgiv(T); cgiv(P);}
   for (pk=p;pk<=B;pk*=p)
   {j=ansize/pk+1; for (i=j*pk;i<=B;i+=pk) {if (j%p) anarray[i]*=cur; j++;}
    cur*=anarray[p]; P=stoi(p);
    if (ggval((GEN) curv[12],P)==0) cur-=((double) p)*anarray[pk/p];
    cgiv(P);
   }
   p+=*dp++;
  }
  ansize=B;
 }
}

void computetwistarray(int d,int B)
{GEN D=stoi(d); byteptr dp=diffptr; ulong p=0; int e; int64 i,pk; GEN P;

 if ((PREVD==d) && (B<=antwsize)) return;
 free(antwarray); antwarray=malloc(u8((B+1)*sizeof(double)));
 memcpy(antwarray,anarray,B*sizeof(double)+8);
 p+=*dp++;
 while(p<=B)
 {P=stoi(p); e=kronecker(D,P); cgiv(P);
  if (e==0) {for(i=p;i<=B;i+=p) antwarray[i]=0.0;}
  if (e==-1)
  {for (pk=p;pk<=B;pk*=p) {for(i=pk;i<=B;i+=pk) antwarray[i]=-antwarray[i];}}
  p+=*dp++;
 }
 PREVD=d; antwsize=B;
}
