#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

uint MM[512];
uint ALL=(1<<30)-1;

int dimF2(int m,int n)
{int i; int j; int r=0; int k; int d; int c[m+1]; uint f=1;

 for (i=0;i<m;i++) c[i]=-1;
 for (k=0;k<n;k++)
 {for (j=0;j<m;j++)
  {if (((MM[j]&f)!=0) && (c[j]==-1))
   {for (i=0;i<m;i++)
    {if (i!=j)
     {d=MM[i]&f; MM[i]&=ALL-f; if (d!=0) MM[i]^=(MM[j]&((ALL+1-f)<<1));}}
    c[j]=k; j=m+5;
   }
  }
  if (j==m) r++;
  f<<=1;
 }
 return(r);
}


void do2divismd(int r)
{int c,i,h,L,s,jj,n,B,S,G,j,k; GEN l,q,p,ppow,GOODPOW,BADPOW,GOOD,BAD;
 int OK=0; char *A; pari_sp memspot=avma;

 GOOD=cgetg(1,t_VEC); BAD=cgetg(1,t_VEC);
 GOODPOW=cgetg(1,t_VEC); BADPOW=cgetg(1,t_VEC); c=0;
 n=itos((GEN) matsize(CF)[1]);
 for (j=1;j<=n;j++)
 {p=(GEN) ((GEN) CF[1])[j]; ppow=gpow(p,(GEN) ((GEN) CF[2])[j],-1);
  if (ellrootno(CURVE,p)==1)
  {GOOD=concat(GOOD,p); GOODPOW=concat(GOODPOW,ppow);}
  else {BAD=concat(BAD,p); BADPOW=concat(BADPOW,ppow);}
 }
 G=glength(GOOD); B=glength(BAD);

 jj=1<<(B+G);
 for (j=jj-1;j>0;j--)
 {k=1; s=0; q=gun; l=gun; for (i=0;i<B+G;i++)
  {L=1<<i; if (L&j)
   {if (i<B) {q=gmul(q,(GEN) BADPOW[i+1]); s++;}
    else q=gmul(q,(GEN) GOODPOW[i-B+1]);
   }
  }
  if (!(s&1))
  {if (gequalgs(q,4)) OK=1;
   for (i=0;i<B+G;i++)
   {L=1<<i; if (!(L&j))
    {if (i<B)
     {k*=(1+kronecker(gneg(q),(GEN) BAD[i+1])); l=gmul(l,(GEN) BADPOW[i+1]);}
     else
     {k*=(1+kronecker(gneg(q),(GEN) GOOD[i-B+1]));
      l=gmul(l,(GEN) GOODPOW[i-B+1]);
     }
    }
   }
  }
  else k=0;
  if (k)
  {h=ggval(l,gdeux);
   if (h<=1) OK=1; if (h==2) {if (mod4(q)!=1) OK=1;}
   if (h>=3) {if (mod8(q)==7) OK=1;}
  }
  if (OK) {MM[c]=(uint) j; c++;} OK=0;
 }
 if (c>0) S=B+G-dimF2(c,B+G); else S=0;

 A=GENtostr((GEN) CURVE[1]); fprintf(OUTFILE,"[%s",A); free(A);
 A=GENtostr((GEN) CURVE[2]); fprintf(OUTFILE,",%s",A); free(A);
 A=GENtostr((GEN) CURVE[3]); fprintf(OUTFILE,",%s",A); free(A);
 A=GENtostr((GEN) CURVE[4]); fprintf(OUTFILE,",%s",A); free(A);
 A=GENtostr((GEN) CURVE[5]); fprintf(OUTFILE,",%s] ",A); free(A);
 A=GENtostr(COND); fprintf(OUTFILE,"%s ",A); free(A); A=GENtostr(MODDEG);
 fprintf(OUTFILE,"%s %li %i %i %i\n",A,ggval(MODDEG,gdeux),r,S,B+G);
 fflush(OUTFILE); free(A); avma=memspot;
}

void sortthem(GEN WT,int *I)
{int i; GEN D,V,W;
 D=cgetg(NI+1,t_VEC); V=cgetg(NI+1,t_VEC);
 for (i=1;i<=NI;i++)
 {D[i]=(long) qtwist((GEN) CL[I[i]],WT);
  V[i]=(long) gdiv(gun,myvol((GEN) D[i]));
 }
 W=indexsort(V); for (i=1;i<=NI;i++) I[i]=itos((GEN) W[i]);
}

void twistdump(int wt,int r,double L,GEN C)
{GEN E,T,p,c2,pv,CCF,WT,l,tam,TWE,ap; int64 pp,cond;
 int i,j,n,pow,cpow; double ra; int I[9]={0,1,2,3,4,5,6,7,8};
 pari_sp memspot=avma;

 if (wt==0) return;
 CCF=factor(C); WT=stoi(wt); n=itos((GEN) matsize(CCF)[1]); TWE=gun;
 cond=gint64(C); printf("%lli [",cond);
 for (j=1;j<=n;j++)
 {p=(GEN) ((GEN) CCF[1])[j]; pp=gint64(p); printf("%lli",pp);
  if (j!=n) printf(","); else printf("] ");
  if ((wt%pp)==0)
  {cpow=ggval(TWCOND,p);
   if (pp!=2)
   {if (cpow>=2) TWE=gmul(TWE,p);
    if (cpow==1) TWE=gmul(TWE,gsub(gsqr(p),gun));
    if (cpow==0)
    {ap=ellap0(TWCURVE,p,0); TWE=gmul(TWE,gsub(p,gun));
     TWE=gmul(TWE,gadd(gadd(p,gun),ap)); TWE=gmul(TWE,gsub(gadd(p,gun),ap));
    }
   }
   else
   {if ((wt%8)==0) TWE=gmul(TWE,volratio(X1VOL,myvol(qtwist(TWCURVE,gdeux))));
    else TWE=gmul(TWE,volratio(X1VOL,myvol(qtwist(TWCURVE,gneg(gun)))));
    pow=itos((GEN) ((GEN) CCF[2])[j]);
    if (cpow>=2) TWE=gmul2n(TWE,pow-cpow);
    if (cpow==1) {TWE=gmul2n(TWE,pow-cpow-2); TWE=gmulgs(TWE,3);}
    if (cpow==0)
    {TWE=gmul2n(TWE,pow-cpow-3); ap=ellap0(TWCURVE,gdeux,0);
     TWE=gmul(TWE,gaddsg(3,ap)); TWE=gmul(TWE,gsubsg(3,ap));
    }
   }
  }
 }
 printf("%i %f %i ",r,L,AISOG);
 if (wt==1) printf("+"); if ((wt==1) && (ISSETZER)) printf("*");
 else if (ISSETZER) TWE=gdivgs(TWE,ISSETZER);
 if ((wt==-4) && (ISSETZER==2) && (ggval(COND,gdeux)>0))
 {printf("*"); TWE=gmul2n(TWE,1);}
 if ((wt==8) && (gequal(COND,gmul2n(gun,5)))) {printf("*"); TWE=gmul2n(TWE,1);}

 if (((wt%8)==0) && (ggval(COND,gdeux)>=5))
 {sortthem(WT,I); if (I[1]!=1) TWE=gmul2n(TWE,-1);}
 else if ((ISOG!=8) && (ISOG!=16)) sortthem(WT,I);
 if (!gcmp0(MODDEG))
 {if ((wt==1) && OFILE) do2divismd(r); output(gmul(TWE,MODDEG));}
 else printf("X\n");
 for (i=1;i<=NI;i++)
 {E=qtwist((GEN) CL[I[i]],WT); printcurven(E); tam=gun;
  if (gsigne((GEN) E[12])==1) printf("["); else printf("(");
  for (j=1;j<=n;j++)
  {p=(GEN) ((GEN) CCF[1])[j]; printf("%li",ggval((GEN) E[12],p));
   l=elllocalred(E,p); tam=gmul(tam,(GEN) l[4]);
   if (j!=n) printf(",");
   else {if (gsigne((GEN) E[12])==1) printf("] "); else printf(") ");}
  }
  c2=cgetg(20,t_VEC); for (j=1;j<=13;j++) c2[j]=E[j];
  c2[17]=(long) gzero; c2[18]=(long) gzero; c2[19]=(long) gzero;
  c2[14]=(long) gzero; c2[15]=(long) gzero; c2[16]=(long) gzero;
  j=4+((lgefint((GEN) c2[12])-2)>>1);
  pv=periodvolvec0(c2,j+2); c2[15]=pv[2]; c2[16]=pv[3];

  T=elltors0(c2,0);
  if (r==0)
  {ra=rtodbl((GEN) pv[2]); if (gsigne((GEN) c2[12])==1) ra*=2.0;
   ra*=gtodouble(tam); ra/=(double) itos(gsqr((GEN) T[1]));
   printf("%i ",(int) (0.5+(L/ra)));
  }
  else printf("X ");
  if (glength((GEN) T[2])>=1) outbrute((GEN) ((GEN) T[2])[1]);
  else printf("1");
  if (glength((GEN) T[2])==2) printf("x"); printf("\n");
 }
 avma=memspot;
}

int whichtwists()
{/* ar = analytic rank? */
 int i,ar,cursize,M; GEN CC; int DONE=0;
 pari_sp memspot;
 double LR; GEN ib=stoi(100000000);

 if (mpcmp(ib,TWCOND)==-1)
 {ar=qtwistar(1,&LR,COND); twistdump(1,ar,LR,COND); return(1);}
 createpile(blocktwists()); cursize=0; i=0;
 while(!DONE)
 {memspot=avma; if (pile[i].which==0) break;
  if ((pile[i].which&3)==3) pile[i].which=-pile[i].which;
  if ((pile[i].which&15)==4) pile[i].which=-pile[i].which;
  CC=gdiv(gmul(TWCOND,gsqr(stoi(pile[i].which))),stoi(pile[i].condcorr));
  if (i==0) {ar=qtwistar(1,&LR,COND); twistdump(1,ar,LR,COND);}
  else if ((mpcmp(ib,CC)==1) && (pile[i].cost<100000000000.0))
  {ar=qtwistar(pile[i].which,&LR,CC); twistdump(pile[i].which,ar,LR,CC);
   if (((pile[i].which&7)==0) && (CM!=-1))
   {ar=qtwistar(-pile[i].which,&LR,CC); twistdump(-pile[i].which,ar,LR,CC);}
  }
  else DONE=1;
  i++; avma=memspot;
 }
 return(cursize);
}

int rootnotwist(int d)
{int eps,i,j,number; GEN condmul,qcurv,p,FACT;

 condmul=gabs(gmul(TWCOND,stoi(d)),-1); qcurv=qtwist(TWCURVE,stoi(d)); eps=-1;
 number=itos((GEN) matsize(TWCF)[1]);
 for(i=1;i<=number;i++)
 {p=(GEN) ((GEN) TWCF[1])[i];
  j=itos((GEN) ((GEN) TWCF[2])[i]); if (j!=0) eps*=ellrootno(qcurv,p);
  while (gequal(gmod(condmul,p),gzero)) condmul=gdiv(condmul,p);
 }
 FACT=factor(condmul); number=itos((GEN) matsize(FACT)[1]);
 for(i=1;i<=number;i++)
 {p=(GEN) ((GEN) FACT[1])[i]; eps*=ellrootno(qcurv,p);}
 return(eps);
}
