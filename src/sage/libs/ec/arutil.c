#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void arround(double *L,int r)
{int T; double R=0.0; int i;
 if (r==0) {*L=floor(sqrt(*L)+0.5); return;}
 T=(int) ceil(log(*L)/log(2.0)); if (T<0) T=0; if (T>15) T=15;
 if (r==1) R=2048.0; if ((r==2) || (r==3)) R=1024.0;
 if ((r==4) || (r==5)) R=512.0; if ((r==6) || (r==7)) R=256.0;
 for (i=0;i<T;i++) *L=*L/2; *L*=R;
 *L=floor(*L); *L/=R;  for (i=0;i<T;i++) *L=*L*2;
}

int qtunhandle(unsigned char A,unsigned char B,int *ar, double *lrat)
{int i,mant,numb;
 int a=(int) A; int b=(int) B;

 if (((*ar)&1)==0)
 {if (a<112) {*ar=0; *lrat=(double) a; return(1);}
  if (a<128) {*ar=0; *lrat=256.0*(double) (a-112)+(double) b; return(2);}
 }
 if (a<128)
 {*ar=1; mant=(a>>3)&15; numb=((a&7)<<8)+b;
  *lrat=(double) numb/2048.0; for (i=0;i<mant;i++) *lrat*=2.0;
 }
 else if (a<192)
 {*ar=2+((*ar)&1); mant=(a>>2)&15; numb=((a&3)<<8)+b;
  *lrat=(double) numb/1024.0; for (i=0;i<mant;i++) *lrat*=2.0;
 }
 else if (a<224)
 {*ar=4+((*ar)&1); mant=(a>>1)&15; numb=((a&1)<<8)+b;
  *lrat=(double) numb/512.0; for (i=0;i<mant;i++) *lrat*=2.0;
 }
 if (a>=224)
 {*ar=6+((*ar)&1); mant=a&15; numb=b;
  *lrat=(double) numb/256.0; for (i=0;i<mant;i++) *lrat*=2.0;
 }
 return(2);
}

void qthandle(int ar,double LR,int *cursize)
{int j; int T; int A=*cursize;

 T=(int) floor(sqrt(LR)+0.5);
 if ((ar==0) && (T<112))
 {A++; TWFERRY=realloc(TWFERRY,u8(A)); TWFERRY[A-1]=(unsigned char) T;}
 else if (ar==0)
 {A+=2; TWFERRY=realloc(TWFERRY,u8(A));
  TWFERRY[A-2]=112+(T>>8); TWFERRY[A-1]=T&255;
 }
 else if (ar==1)
 {A+=2; TWFERRY=realloc(TWFERRY,u8(A)); TWFERRY[A-1]=0;
  T=(int) ceil(log(LR)/log(2.0)); if (T<0) T=0; if (T>15) T=15;
  TWFERRY[A-2]=(T<<3); for (j=0;j<T;j++) LR/=2.0;
  if (LR<1.0)
  {for (j=0;j<3;j++) {if (LR>0.5) {TWFERRY[A-2]+=1<<(2-j); LR-=0.5;} LR*=2.0;}
   for (j=0;j<8;j++) {if (LR>0.5) {TWFERRY[A-1]+=1<<(7-j); LR-=0.5;} LR*=2.0;}
  }
 }
 else if ((ar==2) || (ar==3))
 {A+=2; TWFERRY=realloc(TWFERRY,u8(A)); TWFERRY[A-1]=0;
  T=(int) ceil(log(LR)/log(2.0)); if (T<0) T=0; if (T>15) T=15;
  TWFERRY[A-2]=128+(T<<2); for (j=0;j<T;j++) LR/=2.0;
  if (LR<1.0)
  {for (j=0;j<2;j++) {if (LR>0.5) {TWFERRY[A-2]+=1<<(1-j); LR-=0.5;} LR*=2.0;}
   for (j=0;j<8;j++) {if (LR>0.5) {TWFERRY[A-1]+=1<<(7-j); LR-=0.5;} LR*=2.0;}
  }
 }
 else if ((ar==4) || (ar==5))
 {A+=2; TWFERRY=realloc(TWFERRY,u8(A)); TWFERRY[A-1]=0;
  T=(int) ceil(log(LR)/log(2.0)); if (T<0) T=0; if (T>15) T=15;
  TWFERRY[A-2]=192+(T<<1); for (j=0;j<T;j++) LR/=2.0;
  if (LR<1.0)
  {for (j=0;j<1;j++) {if (LR>0.5) {TWFERRY[A-2]+=1; LR-=0.5;} LR*=2.0;}
   for (j=0;j<8;j++) {if (LR>0.5) {TWFERRY[A-1]+=1<<(7-j); LR-=0.5;} LR*=2.0;}
  }
 }
 else if ((ar==6) || (ar==7))
 {A+=2; TWFERRY=realloc(TWFERRY,u8(A)); TWFERRY[A-1]=0;
  T=(int) ceil(log(LR)/log(2.0)); if (T<0) T=0; if (T>15) T=15;
  TWFERRY[A-2]=224+T; for (j=0;j<T;j++) LR/=2.0;
  if (LR<1.0)
  {for (j=0;j<8;j++) {if (LR>0.5) {TWFERRY[A-1]+=1<<(7-j); LR-=0.5;} LR*=2.0;}}
 }
 *cursize=A;
}

/* int Ordering(struct MYST *A, struct MYST *B) */
int Ordering(__const void *a, __const void *b)
{
  struct MYST *A=(struct MYST*)a, *B=(struct MYST*)b;

  if ((*A).cost<(*B).cost) return(-1); if ((*A).cost>(*B).cost) return(1);
  if ((*A).which<(*B).which) return(-1); return(1);
}

int istwminat(GEN C,int p)
{if (ggval((GEN) C[12],stoi(p))<6) return(1); if (p>3) return(0);
 if (ggval((GEN) C[11],stoi(3))==5) return(1); return(0);
}

int blocktwists()
{int M=1; int i;
 if ((CM==-1) && (ggval(TWCOND,gdeux)==8)) return(2);
 if (CM) return(-CM); if (ISOG==37) return(1); if (ISOG==17) return(17);
 if (ISOG==21) return(3); if (ISOG==15) return(5);
 if (!(AISOG&1)) {if (ggval(TWCOND,gdeux)>=7) M=2;}
 if (((AISOG%3)==0) && (ggval(TWCOND,stoi(3))>=2))
 {for (i=2;i<=NI;i++) {if (!istwminat((GEN) CL[i],3)) {M*=3; i=NI;}}}
 if (((AISOG%5)==0) && (ggval(TWCOND,stoi(5))>=2))
 {for (i=2;i<=NI;i++) {if (!istwminat((GEN) CL[i],5)) {M*=5; i=NI;}}}
 if (((AISOG%7)==0) && (ggval(TWCOND,stoi(7))>=2))
 {for (i=2;i<=NI;i++) {if (!istwminat((GEN) CL[i],7)) {M*=7; i=NI;}}}
 if (((AISOG%13)==0) && (ggval(TWCOND,stoi(13))>=2))
 {for (i=2;i<=NI;i++) {if (!istwminat((GEN) CL[i],13)) {M*=13; i=NI;}}}
 return(M);
}

void createpile(int block)
{int num,numb; int64 i,j,k,l,c; int temp1,temp2; GEN P,c2,pv; int vecsize;
 int64 vec2[16]; int64 vec1[16]; int v1c=0; int v2c=0; int twocond;
 double condsqrt; double COST,S1,S2,T; int skip;

 twocond=ggval(TWCOND,gdeux); condsqrt=rtodbl(gsqrt(TWCOND,4));
 vecsize=itos((GEN) matsize(TWCF)[1]);
 for (i=1;i<=vecsize;i++)
 {P=(GEN) ((GEN) TWCF[1])[i];
  if (!gequal(P,gdeux))
  {l=mpcmp((GEN) ((GEN) TWCF[2]) [i],gun);
   if (l==0) {vec1[v1c]=gint64(P); v1c++;}
   if (l==1) {vec2[v2c]=gint64(P); v2c++;}
  }
 }

 COST=rtodbl(gsqrt(gdivsg(100000000,TWCOND),4));
 numb=(int) ceil(COST); numb=1+(numb>>1);
 pile=realloc(pile,u8(sizeof(struct MYST)*numb));
 for (i=1;i<numb;i++)
 {pile[i].cost=2.0*i-1.0; pile[i].which=2*i-1; pile[i].condcorr=1;}

 num=numb;
 for (j=3;j*j<2*num;j+=2)
 {k=j*j; for (i=1;i<numb;i++)
  {if (((pile[i].which)%k)==0)
   {pile[i].which=pile[numb-1].which; pile[i].condcorr=pile[numb-1].condcorr;
    pile[i].cost=pile[numb-1].cost; i--; numb--;
   }
  }
 }
 for(i=0;i<v2c;i++)
 {for (j=1;j<numb;j++)
  {if ((pile[j].which)%(vec2[i])==0)
   {pile[j].which=pile[numb-1].which; pile[j].condcorr=pile[numb-1].condcorr;
    pile[j].cost=pile[numb-1].cost; j--; numb--;
   }
  }
 }
 for(i=0;i<v1c;i++)
 {for (j=1;j<numb;j++)
  {if ((pile[j].which)%(vec1[i])==0)
   {pile[j].which=pile[numb-1].which; pile[j].condcorr=pile[numb-1].condcorr;
    pile[j].cost=pile[numb-1].cost; j--; numb--;
   }
  }
 }
 pile[0].cost=100000000000000.0; pile[0].which=0;

 pile=realloc(pile,u8(sizeof(struct MYST)*numb*3));
 if (twocond<=3) temp1=1<<twocond; else temp1=16;
 if (twocond<=5) temp2=1<<twocond; else temp2=64;
 S1=4.0/sqrt(temp1); S2=8.0/sqrt(temp2); c=numb;
 for (i=1;i<numb;i++)
 {T=pile[i].cost*S1; if (T<COST)
  {pile[c].cost=T; pile[c].which=pile[i].which*4;
   pile[c].condcorr=pile[i].condcorr*temp1; c++;
  }
  T=pile[i].cost*S2; if (T<COST)
  {pile[c].cost=T; pile[c].which=pile[i].which*8;
   pile[c].condcorr=pile[i].condcorr*temp2; c++;
  }
 }
 numb=c;

 if (CM==-1)
 {for (i=1;i<numb;i++)
  {if ((pile[i].which&7)==4)
   {pile[i].which=pile[numb-1].which; pile[i].condcorr=pile[numb-1].condcorr;
    pile[i].cost=pile[numb-1].cost; i--; numb--;
   }
  }
 }
 if (!(block&1))
 {for (i=1;i<numb;i++)
  {if ((pile[i].which&7)==0)
   {pile[i].which=pile[numb-1].which; pile[i].condcorr=pile[numb-1].condcorr;
    pile[i].cost=pile[numb-1].cost; i--; numb--;
   }
  }
 }

 for (j=0;j<v2c;j++)
 {if ((block%vec2[j])!=0)
  {pile=realloc(pile,u8(sizeof(struct MYST)*numb*2)); c=numb;
   for (i=1;i<numb;i++)
   {pile[c].cost=pile[i].cost; pile[c].which=pile[i].which*vec2[j];
    pile[c].condcorr=pile[i].condcorr*vec2[j]*vec2[j]; c++;
   }
   numb=c;
  }
 }

 for (j=0;j<v1c;j++)
 {pile=realloc(pile,u8(sizeof(struct MYST)*numb*2)); c=numb;
  for (i=1;i<numb;i++)
  {T=pile[i].cost*sqrt(vec1[j]); if (T<COST)
   {pile[c].cost=T; pile[c].which=pile[i].which*vec1[j];
    pile[c].condcorr=pile[i].condcorr*vec1[j]; c++;
   }
  }
  numb=c;
 }
 qsort((struct MYST *) pile,numb,sizeof(struct MYST),Ordering);
 //for (i=0;i<numb;i++)
 //printf("%i %i %f\n",pile[i].which,pile[i].condcorr,pile[i].cost);
}
