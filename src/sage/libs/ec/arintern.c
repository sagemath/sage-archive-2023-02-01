#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

double GfuncM[40][10][21]; double GfuncH[320][10][21];
double invfact[21]; double psarray[10];

void recompM(double t,int w)
{int i; int j;

 for (i=1;i<=20;i++) GfuncM[w][0][i]=-GfuncM[w][0][i-1];
 for (j=1;j<=9;j++)
 {GfuncM[w][j][1]=-GfuncM[w][j-1][0]/t;
  for (i=2;i<=20;i++)
    GfuncM[w][j][i]=((-i+1)*GfuncM[w][j][i-1]-GfuncM[w][j-1][i-1])/t;
 }
 for (j=0;j<=9;j++) for (i=2;i<=20;i++) GfuncM[w][j][i]*=invfact[i];
}

void recompH(double t,int w)
{int i; int j;

 for (i=1;i<=20;i++) GfuncH[w][0][i]=-GfuncH[w][0][i-1];
 for (j=1;j<=9;j++)
 {GfuncH[w][j][1]=-GfuncH[w][j-1][0]/t;
  for (i=2;i<=20;i++)
    GfuncH[w][j][i]=((-i+1)*GfuncH[w][j][i-1]-GfuncH[w][j-1][i-1])/t;
 }
 for (j=0;j<=9;j++) for (i=2;i<=20;i++) GfuncH[w][j][i]*=invfact[i];
}

void analinitw()
{int i; int j; int k;

 psarray[0]=1.0;
 psarray[1]=-0.5772156649015328606065120900;
 psarray[2]=0.9890559953279725553953956514;
 psarray[3]=-0.9074790760808862890165601672;
 psarray[4]=0.9817280868344001873363802939;
 psarray[5]=-0.9819950689031452021047014136;
 psarray[6]=0.9931491146212761931538672531;
 psarray[7]=-0.9960017604424315339700784194;
 psarray[8]=0.9981056937831289219785754027;
 psarray[9]=-0.9990252676219548677946780593;

 invfact[1]=1.0;
 for (i=2;i<=20;i++) invfact[i]=invfact[i-1]/(double) i;

 GfuncM[0][0][0]=0.8187307530779818586699602572;
 GfuncM[0][1][0]=1.222650544183893088334770926;
 GfuncM[0][2][0]=1.160064331844838079811398545;
 GfuncM[0][3][0]=0.8291340003216042069564831197;
 GfuncM[0][4][0]=0.4819084443966072531913594061;
 GfuncM[0][5][0]=0.2379318188879429847848837427;
 GfuncM[0][6][0]=0.1026495824390527947845122350;
 GfuncM[0][7][0]=0.03947110292115808449906886803;
 GfuncM[0][8][0]=0.01372738860259250798508162512;
 GfuncM[0][9][0]=0.004367095636610454091471201768;

 recompM(0.2,0);
 for (j=0;j<39;j++)
 {for (k=0;k<=9;k++)
  {GfuncM[j+1][k][0]=0.0;
   for (i=20;i>=0;i--)
     GfuncM[j+1][k][0]=0.025*GfuncM[j+1][k][0]+GfuncM[j][k][i];
  }
  recompM(0.2+0.025*(double) (j+1),j+1);
 }

 GfuncH[0][0][0]=0.3678794411714423215955609022;
 GfuncH[0][1][0]=0.2193839343955202736772025982;
 GfuncH[0][2][0]=0.09784319721667017932558040432;
 GfuncH[0][3][0]=0.03560349192847501782584410099;
 GfuncH[0][4][0]=0.01107089544600878118841673686;
 GfuncH[0][5][0]=0.003027611195878790266071323620;
 GfuncH[0][6][0]=0.0007426583004868970395963771984;
 GfuncH[0][7][0]=0.0001657562560603850063446475627;
 GfuncH[0][8][0]=0.00003403139486808645492679345187;
 GfuncH[0][9][0]=0.000006482609817428667767728767327;
 recompH(1.0,0);
 for (j=0;j<319;j++)
 {for (k=0;k<=9;k++)
  {GfuncH[j+1][k][0]=0.0;
   for (i=20;i>=0;i--)
     GfuncH[j+1][k][0]=0.125*GfuncH[j+1][k][0]+GfuncH[j][k][i];
  }
  recompH(1.0+0.125*(double) (j+1),j+1);
 }
 INITANALRANK=1;
}

double csm(double h,int k)
{double ans; int i; int n; double hinit;
 double LOG=log(1.0/h); double LOGINIT=LOG;

 i=0; ans=psarray[k];
 for (i=1;i<=k;i++) {ans+=LOG*psarray[k-i]*invfact[i]; LOG*=LOGINIT;}
 hinit=-h; if ((k%2)==0) h=-h;
 for (n=1;n<=10;n++)
 {ans+=h*invfact[n]/pow((double) n,(double) k); h*=hinit;}
 return(ans);
}

double cmd(double h,int n,int k)
{double ans=0.0;
 int i;
 for (i=8;i>=0;i--) ans=h*ans+GfuncM[n][k][i];
 return(ans);
}

double clg(double h,int n,int k)
{double ans=0.0;
 int i;
 for (i=8;i>=0;i--) ans=h*ans+GfuncH[n][k][i];
 return(ans);
}

double analw(int64 n,int k)
{double nd=(double) n;
 double xval=nd*Manal; double sum=0.0; //M is 2*Pi/sqrt(N)
 int A;

 if (xval<=0.19) sum=csm(xval,k)/nd;
 else if (xval<=1.00)
 {A=(int) floor(xval*40-7.5);
  sum=cmd(xval-0.20-(double) A/40.0,A,k)/nd;
 }
 else if (xval<=40.00)
 {A=(int) floor(xval*8-7.5);
  sum=clg(xval-1.00-(double) A/8.0,A,k)/nd;
 }
 else sum=0.0;
 return(sum);
}
