#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ellcurv.h"

char INPUT[1024]; char *a; int64 BOUND=9765625; int NB,nb;

int my_getline() {a=fgets(INPUT,128,stdin); if (NULL==a) return(0); return(1);}

void topcheck()
{char *GO; char GO2[128]; int64 N; int OK=1; int i,cnt;
 double f; int l=0; int r=0;
 GO=malloc(128);
 sscanf(INPUT,"%s",GO); N=atoll(GO); if ((N<11) || (N>BOUND)) OK=0;
 sscanf(INPUT,"%s %s",GO,GO); if ((GO[0]!='[')||(GO[strlen(GO)-1]!=']')) OK=0;
 for (i=0;i<strlen(GO);i++) {if (GO[i]=='[') l++; if (GO[i]==']') r++;}
 if ((l!=1) || (r!=1)) OK=0;
 sscanf(INPUT,"%s %s %s",GO,GO,GO); i=atoi(GO); if ((i<0) || (i>6)) OK=0;
 sscanf(INPUT,"%s %s %s %s",GO,GO,GO,GO); f=atof(GO);
 if ((f<0.0) || (f>100000.0)) OK=0;
 sscanf(INPUT,"%s %s %s %s %s",GO,GO,GO,GO,GO); i=atoi(GO);
 if ((i<1) || (i>163)) OK=0;
 sscanf(INPUT,"%s %s %s %s %s %s",GO,GO,GO,GO,GO,GO);
 if ((!strcmp(GO,"X")) && (!strcmp(GO,"+X")) && (!strcmp(GO,"+*X")))
 {if (GO[0]=='+') GO++; if (GO[0]=='*') GO++; N=atoll(GO); if (N<1) OK=0;}
 cnt=0; for (i=0;i<strlen(INPUT);i++) if (INPUT[i]==' ') cnt++;
 if (cnt!=5) OK=0; if (!OK) {printf("FAILTOP\n%s\n",INPUT); exit(-2);}
 free(GO);
}

int ISOK(int i)
{int j; if (i<1) return(0);
 if (i==1) return(1); if (i==4) return(1); if (i==9) return(1);
 if (i==16) return(1); j=(int) sqrt(i); if (i!=j*j) return(0); return(1);
}

void bigcheck(int CC)
{char GO[128]; char GO2[128]; int l=0; int r=0; int i=0; int OK=1; int cnt;
 GEN Q,PV; int b,j,t;

 sscanf(INPUT,"%s",GO); if ((GO[0]!='[') || (GO[strlen(GO)-1]!=']')) OK=0;
 for (i=0;i<strlen(GO);i++) {if (GO[i]=='[') l++; if (GO[i]==']') r++;}
 if ((l!=1) || (r!=1)) {printf("FAILBIG1\n%s\n",INPUT); exit(-2);}
 sscanf(INPUT,"%s %s",GO,GO); l=0; r=0;
 if (GO[0]=='[')
 {for (i=0;i<strlen(GO);i++) {if (GO[i]=='[') l++; if (GO[i]==']') r++;}
  if ((l!=1) || (r!=1)) OK=0; if (GO[strlen(GO)-1]!=']') OK=0;
 }
 else if (GO[0]=='(')
 {for (i=0;i<strlen(GO);i++) {if (GO[i]=='(') l++; if (GO[i]==')') r++;}
  if ((l!=1) || (r!=1)) OK=0; if (GO[strlen(GO)-1]!=')') OK=0;
 }
 else OK=0;
 sscanf(INPUT,"%s %s %s",GO,GO,GO);
 if (strcmp(GO,"X")) if (!ISOK(atoi(GO))) OK=0;
 sscanf(INPUT,"%s %s %s %s",GO,GO,GO,GO);
 if (GO[strlen(GO)-1]=='x') GO[strlen(GO)-1]=0;
 i=atoi(GO); if ((i<1) || (i>12)) OK=0;
 cnt=0; for (i=0;i<strlen(INPUT);i++) {if (INPUT[i]==' ') cnt++; nb++;}
 if (cnt!=3) OK=0; if (!OK) {printf("FAILBIG2\n%s\n",INPUT); exit(-2);}
 if (nb>CC)
 {sscanf(INPUT,"%s",GO); CURVE=getcurve(GO);
  Q=cgetg(20,t_VEC); for (j=1;j<=13;j++) Q[j]=CURVE[j];
  b=4+((lgefint((GEN) Q[12])-2)>>1); PV=periodvolvec0(Q,b+2);
  Q[15]=PV[2]; Q[16]=PV[3]; Q[14]=(long) gzero; Q[17]=(long) gzero;
  Q[18]=(long) gzero; Q[19]=(long) gzero;
  t=itos((GEN) elltors0(Q,0)[1]); sscanf(INPUT,"%s %s %s %s",GO,GO,GO,GO);
  if (GO[strlen(GO)-1]=='x') {GO[strlen(GO)-1]=0; i=2*atoi(GO);}
  else i=atoi(GO); if (t!=i) {printf("FAILTORS\n%s\n",INPUT);}
 }
}

void processhead(int CC)
{int64 N; char GO[128]; int I,J,b,i;
 pari_sp memspot;
 sscanf(INPUT,"%s",GO); topcheck(); N=atoll(GO); nb+=strlen(INPUT);
 if (GO[0]=='[') {printf("FAILPROC\n%s\n",INPUT); exit(-2);}
 sscanf(INPUT,"%s %s %s %s %s",GO,GO,GO,GO,GO);
 ISOG=atoi(GO); I=classsize(); J=ISOG;
 memspot=avma;
 for (i=0;i<I;i++)
 {b=my_getline(); bigcheck(CC);
  if ((!i) &&  (nb>CC) && 0)
  {sscanf(INPUT,"%s",GO); CURVE=getcurve(GO);
   ROOTNO=-1; TAMA=gun; COND=gun; ISOG=0; PLACE=0; CURVES=gzero;
   ISPRIME=0; ISSQFREE=0; ISSETZER=0; MODDEG=gzero; CM=0; AISOG=0; NI=0;
   CURVE=minimalequation(); minimaltwist(); symsq(); isogenies();
   AISOG=absval(ISOG);
   if (J!=AISOG) {printf("FAILISOG\n%s\n",INPUT); exit(-2);}
  }
 }
 avma=memspot; //if (nb>NB) {printf("%i\n",nb); NB+=1000000;}
}

void checkit(int CC)
{int b;
 if (CC==-1) CC++; CC*=1000000;
 BOUND<<=10; a=malloc(128); VERBOSE=0; NB=CC+1000000;
 while(1) {b=my_getline(); if (!b) break; processhead(CC);}
}
