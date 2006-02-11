#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int WRITE=0; int R0=0; double SUM[16];
int64 *apsm; int64 *ansm; unsigned char *auxsieve;
int SQBND; int64 BOUND; byteptr DP;

void aradd(int64 n,int i,int64 a,int64 la,int64 pj)
{int j=0; int j0; int64 nexta; int r; int g; GEN T; double A=(double) a;
//printf("A %lli %i %lli %lli %lli %f\n",n,i,a,la,pj,analw(n,R0));
//printf("A %lli %lli %f\n",n,a,analw(n,R0));
 if (a==0) j0=i;
 else {for (r=R0;r<10;r+=2) SUM[r]+=A*analw(n,r); pj=0; j0=1;}

 if (n<SQBND) ansm[n]=a;
 for (j=j0;j<=i;j++)
 {pj+=(int64) *(DP+j-1);
  if (pj>BOUND/n) j=i+1;
  else
  {T=stoi(pj); g=ggval(COND,T); cgiv(T);
   if ((g==0) || (apsm[j]!=0))
   {nexta=a*apsm[j];
    if ((j==i) && (g==0))
    {nexta+=-pj*la; aradd(pj*n,j,nexta,a,pj-(int64) *(DP+j-1));}
    else aradd(pj*n,j,nexta,a,pj-(int64) *(DP+j-1));
   }
  }
 }
}


void firstderivs()
{int b,i,r; byteptr dp=diffptr; pari_sp memspot;
 int NUMPRIMES,STEPSIZE,ap,PIX; int64 pi,ppi,mult,multpi;
 GEN pv,curv2,TEMP,BUMP,u; byteptr myptr; double apd,alpha,rpadj;
 FILE *STUFF; char filename[16]="STUFF";

 if (WRITE) STUFF=fopen(filename,"w");
 memspot=avma; auxp=malloc(0); if (!INITANALRANK) analinitw();
 PIX=0; pi=0; while(pi<PMAX) {pi+=*dp++; PIX++;} dp=diffptr; DP=dp;
 Manal=2.0*rtodbl(gdiv(mppi(4),gsqrt(COND,4)));
 if (ROOTNO==1) R0=0; else R0=1;
 rpadj=rtodbl((GEN) periodvolvec0(CURVE,4)[2]);
 if (gsigne((GEN) CURVE[12])==1) rpadj=rpadj*2.0;
 rpadj*=(double) rtodbl(gmul(TAMA,dbltor(1.0)));
 curv2=cgetg(20,t_VEC); for (i=1;i<=13;i++) curv2[i]=CURVE[i];
 b=4+((lgefint((GEN) CURVE[12])-2)>>1);
 pv=periodvolvec0(CURVE,b+2); curv2[15]=pv[2]; curv2[16]=pv[3];
 curv2[14]=(long) gzero; curv2[17]=(long) gzero;
 curv2[18]=(long) gzero; curv2[19]=(long) gzero;
 rpadj/=(double) itos(gsqr(mytors(curv2)));

 u=gpow(gdeux,gsubgs(gdiv(glog(COND,4),glog(stoi(10),4)),6),4);
 alpha=1.74e-09/rpadj; if (mpcmp(gun,u)==-1) alpha*=rtodbl(u);
output(u); printf("%f\n",alpha);
 if (alpha>0.00001) BOUND=(int64) ceil(3*rtodbl(gsqrt(COND,4)));
 else BOUND=(int64) ceil(2*rtodbl(gsqrt(COND,4)));

 // BOUND=(int64) ceil(15.0/Manal);
 if (BOUND<1000) BOUND=1000;
 printf("Need to compute to %lli\n",BOUND); fflush(NULL);
 SQBND=(int) ceil(sqrt(1.01*(double) BOUND));
 if (SQBND<990) SQBND=990;
 NUMPRIMES=(int) floor((double) SQBND/log((double) SQBND)*
                       (1.0+1.5/log((double) SQBND)));
 if (NUMPRIMES<25) NUMPRIMES=25;
 STEPSIZE=(int) ((double) BOUND/(double) log((double) BOUND)/20.0);
 apsm=malloc(u8(NUMPRIMES*8+32)); apsm[0]=0;
 ansm=malloc(u8(SQBND*8+32)); ansm[0]=0; ansm[1]=1;

 for (i=2;i<SQBND+1;i++) ansm[i]=0;
 pi=0; for(i=1;i<NUMPRIMES;i++)
 {pi+=*(dp+i-1); BUMP=stoi(pi); TEMP=ellap0(CURVE,BUMP,0);
  ap=itos(TEMP); cgiv(TEMP); cgiv(BUMP); apsm[i]=ap;
 }

 for (r=R0;r<10;r+=2) SUM[r]=analw((int64) 1,r); pi=0;
 for (i=1;i<NUMPRIMES;i++)
 {ppi=pi; pi+=*(dp+i-1); aradd(pi,i,(int64) apsm[i],1,(int64) ppi);}
 printf("Done with small products  %lli\n",pi); printf("%lli ",pi);
 for (r=R0;r<10;r+=2) printf("%f ",2.0*SUM[r]/rpadj);
 printf("\n"); fflush(NULL);
 dp=dp+NUMPRIMES-1; myptr=diffptr;
 auxp=realloc(auxp,u8(PIX+64));
 getprimes(pi,PSIZE); PMAX=pi+PSIZE; dp=auxp+1; myptr=dp; memspot=avma;

 while (pi<BOUND)
 {pi+=((*dp)<<1); dp++;
  if (pi>PMAX)
  {auxp=realloc(auxp,u8(PIX+64)); getprimes(pi,PSIZE); PMAX+=PSIZE;
   dp=auxp+1; myptr=dp; printf("%lli  ",pi);
   for (r=R0;r<10;r+=2) printf("%f ",2.0*SUM[r]/rpadj);
   printf("\n"); fflush(NULL);
  }
  BUMP=lltoi(pi); TEMP=ellap0(CURVE,BUMP,0);
  ap=itos(TEMP); cgiv(TEMP); cgiv(BUMP);
  //  if (WRITE) fprintf(STUFF,"%i\n",ap);
  if (ap!=0)
  {mult=1; multpi=pi; apd=(double) ap;
   while (multpi<BOUND)
   {for (r=R0;r<10;r+=2) SUM[r]+=analw(multpi,r)*apd*(double) ansm[mult];
    //printf("A %lli %lli %f\n",multpi,ap*ansm[mult],analw(multpi,R0));
    mult++; multpi+=pi;
   }
  }
  avma=memspot;
 }
 printf("FINAL "); for (r=R0;r<10;r+=2) printf("%f ",2.0*SUM[r]/rpadj);
 printf("\n"); fflush(NULL);
 free(auxp); free(apsm); free(ansm); avma=memspot; PMAX=PSIZE;
 //if (WRITE) fclose(STUFF);
}
