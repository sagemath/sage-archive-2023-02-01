#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int *bpsm; int *bnsm; unsigned char *auxsieve;
int NUMPRIMES; int64 BOUND; int SQBND; int PIX;
double M; double MSQPI; double S=0.0;
int INITIALISED=0; double invfactorial[21];
double Ea[40][21]; double Fa[40][21]; double Eb[320][21]; double Fb[320][21];

void goprimes(int64 start,int size,int *j,int *n)
{byteptr mp=diffptr; int i,p,s,l;
//printf("%lli %i %i %i\n",start,size,*j,*n);
 p=0; p+=*mp++; for (i=0;i<size;i++) auxsieve[i]=1;
 while(p<=(int) sqrt(start+size+100))
 {s=(int) (start%p); if (s!=0) l=p-s; else l=0;
  for (i=l;i<size;i+=p) auxsieve[i]=0; p+=*mp++;
 }
 for (i=0;i<size;i++) {(*j)++; if (auxsieve[i]) {auxp[*n]=*j/2; (*n)++; *j=0;}}
}

void getprimes(int64 start,int size)
{int i; int j=0; int pnum=0;
//printf("getprimes %lli %i\n",start,size);
 size+=1000; auxsieve=malloc(u8((size/10)+72));
 for (i=0;i<10;i++) {goprimes(start,size/10,&j,&pnum); start+=size/10;}
 free(auxsieve);
}

void recomp1(double *E,double t)
{int i;
 for (i=4;i<=20;i++) E[i]=(-(i-2)*E[i-1]-2*E[i-3])/t;
 for (i=2;i<=20;i++) E[i]*=invfactorial[i];
}

void recomp2(double *F,double t)
{double tsq=t*t; double n=-tsq*t; int i;
F[4]=(-tsq*F[3]-2*F[1]+2*t*F[2]+2*tsq*F[1])/n;
F[5]=(2*tsq*F[4]+2*tsq*F[2]+4*t*F[1])/n;
F[6]=(5*tsq*F[5]+4*t*F[4]+2*tsq*F[3]+8*t*F[2]+4*F[1])/n;
F[7]=(8*tsq*F[6]+14*t*F[5]+4*F[4]+2*tsq*F[4]+12*t*F[3]+12*F[2])/n;
F[8]=(11*tsq*F[7]+30*t*F[6]+18*F[5]+2*tsq*F[5]+16*t*F[4]+24*F[3])/n;
F[9]=(14*tsq*F[8]+52*t*F[7]+48*F[6]+2*tsq*F[6]+20*t*F[5]+40*F[4])/n;
F[10]=(17*tsq*F[9]+80*t*F[8]+100*F[7]+2*tsq*F[7]+24*t*F[6]+60*F[5])/n;
F[11]=(20*tsq*F[10]+114*t*F[9]+180*F[8]+2*tsq*F[8]+28*t*F[7]+84*F[6])/n;
F[12]=(23*tsq*F[11]+154*t*F[10]+294*F[9]+2*tsq*F[9]+32*t*F[8]+112*F[7])/n;
F[13]=(26*tsq*F[12]+200*t*F[11]+448*F[10]+2*tsq*F[10]+36*t*F[9]+144*F[8])/n;
F[14]=(29*tsq*F[13]+252*t*F[12]+648*F[11]+2*tsq*F[11]+40*t*F[10]+180*F[9])/n;
F[15]=(32*tsq*F[14]+310*t*F[13]+900*F[12]+2*tsq*F[12]+44*t*F[11]+220*F[10])/n;
F[16]=(35*tsq*F[15]+374*t*F[14]+1210*F[13]+2*tsq*F[13]+48*t*F[12]+264*F[11])/n;
F[17]=(38*tsq*F[16]+444*t*F[15]+1584*F[14]+2*tsq*F[14]+52*t*F[13]+312*F[12])/n;
F[18]=(41*tsq*F[17]+520*t*F[16]+2028*F[15]+2*tsq*F[15]+56*t*F[14]+364*F[13])/n;
F[19]=(44*tsq*F[18]+602*t*F[17]+2548*F[16]+2*tsq*F[16]+60*t*F[15]+420*F[14])/n;
F[20]=(47*tsq*F[19]+690*t*F[18]+3150*F[17]+2*tsq*F[17]+64*t*F[16]+480*F[15])/n;
for (i=2;i<=20;i++) F[i]*=invfactorial[i];
}

void InitialiseWeights()
{int i,j,p; byteptr dp=diffptr;

 PIX=0; p=0; while(p<PMAX) {p+=*dp++; PIX++;}

 invfactorial[1]=1.0;
 for (i=2;i<=20;i++) invfactorial[i]=invfactorial[i-1]/(double) i;

 Ea[0][0]=0.571220809781621776397779175329;
 Ea[0][1]=-1.17552038116431292642263284582;
 Ea[0][2]=4.17235578464720035925283455107;
 Ea[0][3]=-26.5739870210522195602419645087;

 Fa[0][0]=0.942023312587544175347448939461;
 Fa[0][1]=-.41671112528252129492049970875;
 Fa[0][2]=-.60449401084343923398903234101;
 Fa[0][3]=5.37038302981623065266017263208;

 Eb[0][0]=0.16930874045437966393905180;
 Eb[0][1]=-.21092648721548305615222730;
 Eb[0][2]=0.33075813243002310703316760;
 Eb[0][3]=-.66937561333878243491127120;

 Fb[0][0]=0.59321869178228331361583278;
 Fb[0][1]=-.37385746452305603724210338;
 Fb[0][2]=0.21239606102085506448239242;
 Fb[0][3]=-.01393033247674442378267187;

 recomp1(Ea[0],0.2); recomp2(Fa[0],0.2);
 for (j=0;j<39;j++)
 {Ea[j+1][0]=0.0; Ea[j+1][1]=0.0; Ea[j+1][2]=0.0; Ea[j+1][3]=0.0;
  for (i=20;i>=0;i--) Ea[j+1][0]=0.025*Ea[j+1][0]+Ea[j][i];
  for (i=20;i>=1;i--) Ea[j+1][1]=0.025*Ea[j+1][1]+i*Ea[j][i];
  for (i=20;i>=2;i--) Ea[j+1][2]=0.025*Ea[j+1][2]+i*(i-1)*Ea[j][i];
  for (i=20;i>=3;i--) Ea[j+1][3]=0.025*Ea[j+1][3]+i*(i-1)*(i-2)*Ea[j][i];
  recomp1(Ea[j+1],0.2+0.025*(double) (j+1));
  Fa[j+1][0]=0.0; Fa[j+1][1]=0.0; Fa[j+1][2]=0.0; Fa[j+1][3]=0.0;
  for (i=20;i>=0;i--) Fa[j+1][0]=0.025*Fa[j+1][0]+Fa[j][i];
  for (i=20;i>=1;i--) Fa[j+1][1]=0.025*Fa[j+1][1]+i*Fa[j][i];
  for (i=20;i>=2;i--) Fa[j+1][2]=0.025*Fa[j+1][2]+i*(i-1)*Fa[j][i];
  for (i=20;i>=3;i--) Fa[j+1][3]=0.025*Fa[j+1][3]+i*(i-1)*(i-2)*Fa[j][i];
  recomp2(Fa[j+1],0.2+0.025*(double) (j+1));
 }

 recomp1(Eb[0],1.0); recomp2(Fb[0],1.0);
 for (j=0;j<319;j++)
 {Eb[j+1][0]=0.0; Eb[j+1][1]=0.0; Eb[j+1][2]=0.0; Eb[j+1][3]=0.0;
  for (i=20;i>=0;i--) Eb[j+1][0]=0.125*Eb[j+1][0]+Eb[j][i];
  for (i=20;i>=1;i--) Eb[j+1][1]=0.125*Eb[j+1][1]+i*Eb[j][i];
  for (i=20;i>=2;i--) Eb[j+1][2]=0.125*Eb[j+1][2]+i*(i-1)*Eb[j][i];
  for (i=20;i>=3;i--) Eb[j+1][3]=0.125*Eb[j+1][3]+i*(i-1)*(i-2)*Eb[j][i];
  recomp1(Eb[j+1],1.0+0.125*(double) (j+1));
  Fb[j+1][0]=0.0; Fb[j+1][1]=0.0; Fb[j+1][2]=0.0; Fb[j+1][3]=0.0;
  for (i=20;i>=0;i--) Fb[j+1][0]=0.125*Fb[j+1][0]+Fb[j][i];
  for (i=20;i>=1;i--) Fb[j+1][1]=0.125*Fb[j+1][1]+i*Fb[j][i];
  for (i=20;i>=2;i--) Fb[j+1][2]=0.125*Fb[j+1][2]+i*(i-1)*Fb[j][i];
  for (i=20;i>=3;i--) Fb[j+1][3]=0.125*Fb[j+1][3]+i*(i-1)*(i-2)*Fb[j][i];
  recomp2(Fb[j+1],1.0+0.125*(double) (j+1));
 }
 // for (i=32;i<40;i++) for (j=0;j<20;j++)
 // printf("%i %i %.30f %.30f\n",i,j,Ea[i][j],Fa[i][j]);
 INITIALISED=1;
}

double csm1(double h)
{double ans=1.0;
 double LOG=log(h);
 ans=ans+h*
   (-.151401970301401363712811834249+h*
    (-1.00000000000000000000000000000+h*
     (0.275984587738125243706837284175+h*
      (0.055555555555555555555555555555+h*
       (-.010190720597481709433201395571+h*
        (-.000740740740740740740740740741+h*
         (0.000098647366279703548252545702+h*
          (0.000003779289493575207860922147+h*
           (-.000000389312193549144569237960)))))))));
 ans=ans+LOG*h*
   (1.12837916709551257389615890312+h*h*
    (-.18806319451591876231602648385+h*h*
     (0.00470157986289796905790066210+h*h*
      (-.00003731412589601562744365605+h*h*
       (0.00000012956293713894315084602)))));
 return(ans);
}

double csm2(double h)
{double ans=1.0;
 double LOG=log(h);
 ans=ans+h*h*
   (0.365823497352299290909768135130+h*
    (-1.18163590060367735153211165556+h*
     (0.346044125661925177272557966217+h*
      (0.078775726706911823435474110371+h*
       (-.014820670157275699368682165728+h*
        (-.001125367524384454620506773005+h*
         (0.000151958516742729593504450147+h*
          (0.000005954325525843675240776577+h*
           (-.0000006187376637259239072626810)))))))));
 ans=ans+LOG*h*h*
   (1.0000000000000000000000000000000+h*h*
    (-.2500000000000000000000000000000+h*h*
     (0.0069444444444444444444444444444+h*h*
      (-.0000578703703703703703703703704+h*h*
       (0.000000206679894179894179894180)))));
 return(ans);
}

double cmd1(double h, int n)
{double ans=0.0; int i;
 for (i=8;i>=0;i--) ans=h*ans+Ea[n][i];
 return(ans);
}

double cmd2(double h, int n)
{double ans=0.0; int i;
 for (i=8;i>=0;i--) ans=h*ans+Fa[n][i];
 return(ans);
}

double clg1(double h, int n)
{double ans=0.0; int i;
 for (i=8;i>=0;i--) ans=h*ans+Eb[n][i];
 return(ans);
}

double clg2(double h, int n)
{double ans=0.0; int i;
 for (i=8;i>=0;i--) ans=h*ans+Fb[n][i];
 return(ans);
}

double WEIGHT(int64 n)
{double nd=(double) n;
 double xval1=nd*M; double xval2=xval1;
 double sum1=0.0; double sum2=0.0;
 int A;

 //printf("W %f\n",xval1);
 if (xval1<=0.19) sum1=csm1(xval1)*MSQPI/nd;
 else if (xval1<=1.00)
 {A=(int) floor(xval1*40-7.5);
 sum1=cmd1(xval1-0.20-(double) A/40.0,A)*MSQPI/nd;
 }
 else if (xval1<=40.00)
 {A=(int) floor(xval1*8-7.5);
 sum1=clg1(xval1-1.00-(double) A/8.0,A)*MSQPI/nd;
 }
 else sum1=0.0;

 if (xval2<=0.19) sum2=csm2(xval2)/nd/nd;
 else if (xval2<=1.00)
 {A=(int) floor(xval2*40-7.5);
  sum2=cmd2(xval2-0.20-(double) A/40.0,A)/nd/nd;
 }
 else if (xval2<=40.00)
 {A=(int) floor(xval2*8-7.5);
  sum2=clg2(xval2-1.00-(double) A/8.0,A)/nd/nd;
 }
 else sum2=0.0;
 // printf("S %lli %f %f %f %f\n",n,xval1,sum1,sum2,sum1+sum2);

 return(sum1+sum2);
}

void ADDLL(int64 n,int i,int64 a,int64 la,int64 lla,int64 pj,int64 sscond)
{int j=0; int j0; int64 nexta; byteptr dp=diffptr;

//printf("A %lli %i %lli %lli %lli %lli %lli %.12f\n",
//     n,i,a,la,lla,pj,sscond,S);
 if (a==0) j0=i; else {S+=(double) a*WEIGHT(n); pj=0; j0=1;}
 //printf("R %f\n",S);
 if (n<SQBND) bnsm[n]=a;
 for (j=j0;j<=i;j++)
 {pj+=(int64) *(dp+j-1);
  if (pj>BOUND/n) j=i+1;
  else if (((sscond%pj)==0) && (bpsm[j]==0)) ;
  else
  {nexta=a*bpsm[j];
   if ((j==i) && ((sscond%pj)!=0))
   {nexta+=pj*pj*pj*lla-pj*bpsm[j]*la;
    ADDLL(pj*n,j,nexta,a,la,pj-(int64) *(dp+j-1),sscond);
   }
   else ADDLL(pj*n,j,nexta,a,0,pj-(int64) *(dp+j-1),sscond);
  }
 }
}

/* #define DBG(arg) printf("degphi: %s\n",arg); */
#define DBG(arg)

void degphi(GEN CURV,GEN tweffect,int goon)
{int64 pi=0; int i=1; int bnd=0; int64 ppi; int ap=0;
 GEN TEMP,BUMP; GEN LAST=gzero; GEN GUARD,fundvol;
 int64 bp=0; int STEPSIZE=0; int64 mult; int64 multpi; double bpd;
 byteptr myptr; byteptr dp=diffptr; int64 sscond,explus,exminus;
 pari_sp myspot,memspot;

 DBG("in degphi");
 memspot=avma; myspot=avma; MODDEG=gzero; auxp=malloc(0);
 sscond=gint64(SSCOND); explus=gint64(EXPOS); exminus=gint64(EXNEG);
 //printf("%lli %lli %lli\n",sscond,explus,exminus);

 GUARD=stoi(100); fundvol=(GEN) periodvolvec0(CURV,4)[1];
 M=11.1366559936634156905696359642/(double) sscond;
 MSQPI=1.77245385090551602729816748334*M;
 bnd=2*(15-((int) 2.0*log(rtodbl(fundvol))));
 BOUND=(int64) ceil((double) bnd/M); bnd/=2;
 if (BOUND<1000) BOUND=1000;
 SQBND=(int) ceil(sqrt(1.01*(double) BOUND));
 if (SQBND<990) SQBND=990;
 NUMPRIMES=(int) floor((double) SQBND/log((double) SQBND)*
		       (1.0+1.5/log((double) SQBND)));
 if (NUMPRIMES<25) NUMPRIMES=25;
 STEPSIZE=(int) ((double) BOUND/(double) log((double) BOUND)/20.0);

 DBG("1");

 if (!INITIALISED) InitialiseWeights(); timer();
 bpsm=malloc(u8(NUMPRIMES*4+32)); bpsm[0]=0;
 bnsm=malloc(u8(SQBND*4+32)); bnsm[0]=0; bnsm[1]=1;

 for (i=2;i<SQBND+1;i++) bnsm[i]=0;
 pi=0; for(i=1;i<NUMPRIMES;i++)
 {pi+=*(dp+i-1); BUMP=stoi(pi);
  if (pi>ansize) {TEMP=ellap0(CURV,BUMP,0); ap=itos(TEMP); cgiv(TEMP);}
  else ap=(int) anarray[pi];
  cgiv(BUMP);
  if ((sscond%pi)==0)
  {if ((explus%pi)==0) bpsm[i]=-pi; //REVERSED
   else if ((exminus%pi)==0) bpsm[i]=pi;
   else bpsm[i]=ap*ap;
  }
  else bpsm[i]=ap*ap-pi;
 }

 S=WEIGHT(1); pi=0;
 for (i=1;i<NUMPRIMES;i++)
 {ppi=pi; pi+=*(dp+i-1);
  ADDLL(pi,i,(int64) bpsm[i],1,0,(int64) ppi,sscond);
 }

 if (TRACE) printf("%li\n",timer());
 dp=dp+NUMPRIMES-1; myptr=diffptr;
 auxp=realloc(auxp,u8(PIX+64));
 getprimes(pi,PSIZE); PMAX=pi+PSIZE; dp=auxp+1; myptr=dp;
 if (TRACE) printf("%lli %lli\n",pi,pi+((*dp)<<1));
 memspot=avma; myspot=avma;
 while (1){
   for (i=0;i<STEPSIZE;i++)
  {pi+=(int64) ((*dp)<<1); dp++;
   if (pi>PMAX)
   {if (!goon) {MODDEG=gzero; break;}
    auxp=realloc(auxp,u8(PIX+64)); getprimes(pi,PSIZE); PMAX+=PSIZE;
    dp=auxp+1; myptr=dp;
    if (TRACE) printf("%lli %lli %lli\n",pi,pi+((*dp)<<1),PMAX);
   }
   if (pi>ansize) ap=itos(ellap0(CURV,lltoi(pi),0)); else ap=(int) anarray[pi];
   if ((sscond%pi)==0)
   {if ((exminus%pi)==0) bp=pi; //REVERSED
    else if ((explus%pi)==0) bp=-pi;
    else bp=(int64) ap*(int64) ap;
   }
   else bp=(int64) ap*(int64) ap-pi;
   if (bp!=0)
   {mult=1; multpi=pi; bpd=(double) bp;
    while (multpi<BOUND)
    {S+=WEIGHT(multpi)*bpd*bnsm[mult];
      mult++; multpi+=pi;
    }
   }
   avma=memspot;
  }
/*  if (pi>PMAX) break; */
  if ((pi>2*BOUND) && (goon))
  {MODDEG=gdiv(dbltor(S),gmul2n(mppi(4),1));
   MODDEG=gmul(MODDEG,gdiv(gmul(EXEFF,TWCOND),fundvol));
   MODDEG=gmul(bestappr(MODDEG,gun),tweffect);
   break;
  }
  if (pi>BOUND/6)
  {MODDEG=gdiv(dbltor(S),gmul2n(mppi(4),1));
   MODDEG=gmul(MODDEG,gdiv(gmul(EXEFF,TWCOND),fundvol));
   MODDEG=gmul(bestappr(MODDEG,GUARD),tweffect);
   if(typ(MODDEG)==t_INT)
   {if (gequal(MODDEG,gneg(LAST))) LAST=gcopy(MODDEG);
    else if (gequal(MODDEG,LAST)) break;
    else LAST=gneg(MODDEG);
   }
   else LAST=gzero;
   avma=myspot; LAST=gcopy(LAST); memspot=avma;
   if (TRACE) {printf("%lli %li %.11f  ",pi,timer(),S); output(MODDEG);}
  }
 }

 free(auxp);
 free(bpsm);
 free(bnsm);

 PMAX=PSIZE;
}
