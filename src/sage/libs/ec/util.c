#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

int plwh8[6][6]={{1,2,4,4,8,8},{2,1,2,2,4,4},{4,2,1,4,8,8},
		 {4,2,4,1,2,2},{8,4,8,2,1,4},{8,4,8,2,4,1}};
int plwh18[6][6]={{1,3,9,2,6,18},{3,1,3,6,2,6},{9,3,1,18,6,2},
		  {2,6,18,1,3,9},{6,2,6,3,1,3},{18,6,2,9,3,1}};
int plwhn8[6][6]={{1,2,2,2,4,4},{2,1,4,4,8,8},{2,4,1,4,8,8},
		  {2,4,4,1,2,2},{4,8,8,2,1,4},{4,8,8,2,4,1}};
int plwhn18[6][6]={{1,3,3,2,6,6},{3,1,9,6,2,18},{3,9,1,6,18,2},
		   {2,6,6,1,3,3},{6,2,18,3,1,9},{6,18,2,3,9,1}};
int plwh12[8][8]={{1,2,4,4,3,6,12,12},{2,1,2,2,6,3,6,6},
		  {4,2,1,4,12,6,3,6},{4,2,4,1,12,6,12,3},
		  {3,6,12,12,1,2,4,4},{6,3,6,6,2,1,2,2},
		  {12,6,3,6,4,2,1,4},{12,6,12,3,4,2,4,1}};
int plwh16[8][8]={{1,2,4,4,8,8,16,16},{2,1,2,2,4,4,8,8},
		  {4,2,1,4,8,8,16,16},{4,2,4,1,2,2,4,4},
		  {8,4,8,2,1,4,8,8},{8,4,8,2,4,1,2,2},
		  {16,8,16,4,8,2,1,4},{16,8,16,4,8,2,4,1}};
int plwh32[8][8]={{1,2,4,8,8,4,8,8},{2,1,2,4,4,2,4,4},
		  {4,2,1,2,2,4,8,8},{8,4,2,1,4,8,16,16},
		  {8,4,2,4,1,8,16,16},{4,2,4,8,8,1,2,2},
		  {8,4,8,16,16,2,1,4},{8,4,8,16,16,2,4,1}};
int plwhn12[8][8]={{1,2,2,2,3,6,6,6},{2,1,4,4,6,3,12,12},
		   {2,4,1,4,6,12,3,12},{2,4,4,1,6,12,12,3},
		   {3,6,6,6,1,2,2,2},{6,3,12,12,2,1,4,4},
		   {6,12,3,12,2,4,1,4},{6,12,12,3,2,4,4,1}};
int plwhn16[8][8]={{1,2,2,2,4,4,8,8},{2,1,4,4,8,8,16,16},
		   {2,4,1,4,8,8,16,16},{2,4,4,1,2,2,4,4},
		   {4,8,8,2,1,4,8,8},{4,8,8,2,4,1,2,2},
		   {8,16,16,4,8,2,1,4},{8,16,16,4,8,2,4,1}};

int absval(int INP) {if (INP<0) return(-INP); return(INP);}
GEN lltoi(int64 lli)
{int n=0; GEN RET;
 if (lli<0) {n=1; lli=-lli;}
 RET=stoi(lli&((1<<22)-1));
 RET=gadd(RET,gmul2n(stoi((lli>>22)&((1<<22)-1)),22));
 RET=gadd(RET,gmul2n(stoi(lli>>44),44));
 if (n) RET=gneg(RET); return(RET);
}

int64 gint64(GEN inp)
{int64 temp; GEN TEMP; int n=0;
 if (gsigne(inp)==-1) {n=1; inp=gneg(inp);}
 // if(mpcmp(gmul2n(gun,31),inp)==1) return((int64) itos(inp));
 TEMP=gmod(inp,gmul2n(gun,31)); temp=(int64) itos(TEMP);
 inp=gmul2n(gsub(inp,TEMP),-31);
 temp+=(((int64) itos(inp))<<31); if (n) temp=-temp;
 return(temp);
}

int classsize()
{
 switch (ISOG)
 {case 1: return(1);
  case 2: case 3: case 5: case 7: case 11: case 13: case 17:
   case 19: case 37: case 43: case 67: case 163: return(2);
  case 9: case 25: case -9: case -25: return(3);
  case 4: case 6: case 10: case 14: case 15:
   case 21: case 27: case -4: return(4);
  case 8: case 18: case -8: case -18: return(6);
  case 12: case 16: case 32: case -12: case -16: case -32: return(8);
 }
  return(0);
}

int isogdeg(int which)
{int i;

 if (ISOG==1) return(1); if (which==PLACE) return(1);
 i=classsize(); if (i==2) return(ISOG);
 if (i==3)
 {if (ISOG==-9) return(3); if (ISOG==-25) return(5);
  if (ISOG==9)
  {if ((PLACE-which)==1) return(3); if ((PLACE-which)==2) return(9);
   if ((PLACE-which)==-2) return(9); if ((PLACE-which)==-1) return(3);
  }
  if (ISOG==25)
  {if ((PLACE-which)==1) return(5); if ((PLACE-which)==2) return(25);
   if ((PLACE-which)==-2) return(25); if ((PLACE-which)==-1) return(5);
  }
 }
 if (i==4)
 {if (ISOG==4)
  {if ((PLACE==2) || (which==2)) return(2); else return(4);}
  if (ISOG==-4)
  {if ((PLACE==1) || (which==1)) return(2); else return(4);}
  if (ISOG==6)
  {if ((PLACE+which)==5) return(6);
   if ((PLACE-which)==2) return(3); if ((PLACE-which)==-2) return(3);
   if ((PLACE-which)==1) return(2); if ((PLACE-which)==-1) return(2);
  }
  if (ISOG==10)
  {if ((PLACE+which)==5) return(10);
   if ((PLACE-which)==2) return(5); if ((PLACE-which)==-2) return(5);
   if ((PLACE-which)==1) return(2); if ((PLACE-which)==-1) return(2);
  }
  if (ISOG==14)
  {if ((PLACE+which)==5) return(14);
   if ((PLACE-which)==2) return(7); if ((PLACE-which)==-2) return(7);
   if ((PLACE-which)==1) return(2); if ((PLACE-which)==-1) return(2);
  }
  if (ISOG==15)
  {if ((PLACE+which)==5) return(15);
   if ((PLACE-which)==2) return(5); if ((PLACE-which)==-2) return(5);
   if ((PLACE-which)==1) return(3); if ((PLACE-which)==-1) return(3);
  }
  if (ISOG==21)
  {if ((PLACE+which)==5) return(21);
   if ((PLACE-which)==2) return(7); if ((PLACE-which)==-2) return(7);
   if ((PLACE-which)==1) return(3); if ((PLACE-which)==-1) return(3);
  }
  if (ISOG==27)
  {if ((PLACE+which)==5) return(9);
   if ((PLACE+which)==6) return(27);
   if ((PLACE+which)==3) return(3);
   if ((PLACE+which)==4) return(3);
   if ((PLACE+which)==7) return(3);
  }
 }
 if (i==6)
 {if (ISOG==8) return(plwh8[PLACE-1][which-1]);
  if (ISOG==18) return(plwh18[PLACE-1][which-1]);
  if (ISOG==-8) return(plwhn8[PLACE-1][which-1]);
  if (ISOG==-18) return(plwhn18[PLACE-1][which-1]);
 }
 if (i==8)
 {if (ISOG==12) return(plwh12[PLACE-1][which-1]);
  if (ISOG==16) return(plwh16[PLACE-1][which-1]);
  if (ISOG==32) return(plwh32[PLACE-1][which-1]);
  if (ISOG==-12) return(plwhn12[PLACE-1][which-1]);
  if (ISOG==-16) return(plwhn16[PLACE-1][which-1]);
  if (ISOG==-32) return(plwh32[PLACE-1][which-1]);
 }
 return(-1);
}

void printcurve(GEN curv)
{
 printf("[");
 outbrute((GEN) (curv)[1]); printf(",");
 outbrute((GEN) (curv)[2]); printf(",");
 outbrute((GEN) (curv)[3]); printf(",");
 outbrute((GEN) (curv)[4]); printf(",");
 outbrute((GEN) (curv)[5]); printf("]\n");
}

void printcurven(GEN curv)
{
 printf("[");
 outbrute((GEN) (curv)[1]); printf(",");
 outbrute((GEN) (curv)[2]); printf(",");
 outbrute((GEN) (curv)[3]); printf(",");
 outbrute((GEN) (curv)[4]); printf(",");
 outbrute((GEN) (curv)[5]); printf("] ");
}

int u8(int x) {return(((x+7)>>3)<<3);}

int checkconductor(GEN C)
{GEN D,p,l; int number,g,i;

 D=gabs((GEN) C[12],-1); number=itos((GEN) matsize(CF)[1]);
 for(i=1;i<=number;i++)
 {p=(GEN) ((GEN) CF[1])[i]; g=ggval(D,p); if (g<1) return(0);
  l= gel(elllocalred(C,p),1);
  if (!gequal(l,(GEN) ((GEN) CF[2])[i])) return(0);
  D=gdiv(D,gpow(p,stoi(g),-1));
 }
 if (gequal(D,gun)) return(1); return(0);
}

int gettheint(unsigned char *IN,int i)
{int A=0;
 A+=((int) IN[i])<<24; A+=((int) IN[i+1])<<16;
 A+=((int) IN[i+2])<<8; A+=(int) IN[i+3]; return(A);
}

void get4bytes(int w,unsigned char *OUT)
{OUT[3]=(unsigned char) (w&255); w=w>>8;
 OUT[2]=(unsigned char) (w&255); w=w>>8;
 OUT[1]=(unsigned char) (w&255); w=w>>8;
 OUT[0]=(unsigned char) (w&255);
}

int uwlsb(unsigned char *IN,int i)
{int A=(int) IN[4*i];
 //printf("%3i %3i %3i %3i  ",IN[4*i],IN[4*i+1],IN[4*i+2],IN[4*i+3]);
 if (A>=240) {A-=240; A<<=8; A+=(int) IN[4*i+1];}
 if (IN[4*i]==255)
 {A-=3840; A<<=8; A+=(int) IN[4*i+2];
  if (IN[4*i+1]>=240) {A-=61440; A<<=8; A+=(int) IN[4*i+3];}
 }
 //printf("%i\n",A);
 return(A);
}

void cv(int w,GEN C)
{int j; for (j=1;j<=5;j++) {((GEN) CURVES[w])[j]=C[j];}}

void docount()
{char filename[80]; FILE *infile; unsigned char headbytes[32];
 int nc; int running=0; pari_sp memspot; int i,j; GEN C4,C6; int total=0;

 while(1)
 {scanf("%s",filename); if (!strcmp(filename,"END")) break; nc=0;
  if (PCOUNT)
  {readfile(filename); total+=numcurves;
   for (i=0;i<numcurves;i++)
   {memspot=avma; C4=gzero; C6=gzero;
    for (j=0;j<c4s;j++)
    {C4=gadd(gmul2n(C4,8),stoi((int) curvedata[i*recsize+j]));}
    for (j=0;j<c6s;j++)
    {C6=gadd(gmul2n(C6,8),stoi((int) curvedata[i*recsize+c4s+j]));}
    if (c4neg) C4=gneg(C4); if (c6neg) C6=gneg(C6);
    C4=gadd(stoi(c4r),gmul(C4,stoi(576)));
    C6=gadd(stoi(c6r),gmul(C6,stoi(1728))); CURVE=avecfromc4c6(C4,C6);
    if (isprime(gabs((GEN) CURVE[12],-1))) nc++; avma=memspot;
   }
   free(curvedata);
   free(size16bytes); free(size256bytes); free(size4096bytes);
   free(leadbyte); free(moddegbytes); free(localsizebytes);
   for (i=0;i<numcurves;i++) free(twdata[i]); free(twdata); running+=nc;
   printf("%s %i %i %i %i\n",filename,numcurves,total,nc,running);
  }
  else
  {infile=fopen(filename,"r"); fread(headbytes,1,8,infile); fclose(infile);
   nc=(((int) headbytes[5])<<16)+
     (((int) headbytes[6])<<8)+((int) headbytes[7]); running+=nc;
   printf("%s %i %i\n",filename,nc,running);
  }
 }
}

GEN isprimepower(GEN X)
{int i,u,g,sz; GEN G,p,L;

 if (gequal(X,gzero)) return(gzero); if (gequal(X,gun)) return(gun);
 if (!mod2(X))
 {g=ggval(X,gdeux);
  if (gequal(gpow(gdeux,stoi(g),-1),X)) return(gdeux); else return(gzero);
 }
 p=stoi(3); g=ggval(X,p);
 if (g!=0) {if (gequal(gpow(p,stoi(g),-1),X)) return(p); else return(gzero);}
 p=stoi(5); g=ggval(X,p);
 if (g!=0) {if (gequal(gpow(p,stoi(g),-1),X)) return(p); else return(gzero);}
 p=stoi(7); g=ggval(X,p);
 if (g!=0) {if (gequal(gpow(p,stoi(g),-1),X)) return(p); else return(gzero);}

 if (isprime(X)) return(X); sz=taille(X);
 L=glog(X,3);  u=itos(gceil(gdiv(L,glog(stoi(11),3))));
 for (i=2;i<=u;i++)
 {G=ground(gpow(X,gdiv(gun,stoi(i)),(sz+3)/2));
  if (gequal(gpow(G,stoi(i),-1),X)) {if (isprime(G)) return(G);}
 }
 return(gzero);
}

int KnownDiffOptimal(GEN C)
{GEN p,D,E,c4,F,P,Q,CF,S; int i,n,OK;

 C=ellinit0(C,1,ELLACC); X1VOL=myvol(C); D=(GEN) C[12]; c4=(GEN) C[10];
 if (gequal(gun,c4) && gequal(D,stoi(-15))) return(4);
 if (gequal(c4,stoi(96)) && gequal(D,stoi(80))) return(2);
 if (gequal(c4,gmul2n(gun,4)) && gequal(D,stoi(-11))) return(5);
 if (gequal(c4,stoi(48)) && gequal(D,stoi(64))) return(2);
 if (gequal(c4,stoi(-48)) && gequal(D,stoi(-64))) return(2);
 if (gequal(c4,stoi(-32)) && gequal(D,stoi(-48))) return(2);
 if (gequal(c4,stoi(112)) && gequal(D,stoi(128))) return(2);
 if (gequal(c4,stoi(33)) && gequal(D,stoi(17))) return(4);

 E=gsub(D,gmul2n(gun,6));
 if (gsigne(D)==1)
 {if (Z_issquarerem(E, NULL))
  {if (isprime(D)) {if (gequal(c4,gsub(D,gmul2n(gun,4)))) return(2);}
   if (mod16(D)==0)
   {if (isprime(gmul2n(D,-4)))
    {if (gequal(c4,gsub(D,gmul2n(gun,4)))) return(2);}
   }
  }
 }
 E=gadd(D,gmul2n(gun,6));
 if (Z_issquarerem(E,&S))
 {p=gsub(S,gmul2n(gun,3));
  P=isprimepower(gabs(p,-1)); Q=isprimepower(gadd(p,gmul2n(gun,4)));
  if (!gequal(P,gzero) && !gequal(Q,gzero))
  {if (gequal(P,gun)) return(1);
   if ((gequal(gmodulcp(P,gsqr(gdeux)),stoi(3))) ||
       (gequal(gmodulcp(Q,gsqr(gdeux)),stoi(3)))) return(2);
  }
 }
 E=gadd(D,stoi(27)); if (gequal(E,gzero)) return(1);
 F=icbrt(E);
 if (!gequal(F,gzero))
 {E=gmul(F,gsub(gmul(gsqr(F),F),stoi(24)));
  if (gequal(E,c4))
  {F=gsub(F,stoi(3));
   if (!gequal(isprimepower(gadd(gsqr(F),
                                gadd(gmul(F,stoi(9)),stoi(27)))),gzero))
   {CF=factor(F); n=itos((GEN) matsize(CF)[1]); OK=1;
    for(i=1;i<=n;i++)
    {E=(GEN) ((GEN) CF[1])[i]; if (gequal(gmodulcp(E,stoi(3)),gun)) OK=0;}
    if (OK) return(3);
   }
  }
 }
 return(0);
}

GEN icbrt(GEN X)
{int n=0; int sz; GEN G;

 if (gequal(X,gzero)) return(gzero);
 if (gsigne(X)==-1) {n=1; X=gneg(X);}
 sz=taille(X); //output(X); printf("ICBRT %i\n",sz);
 G=ground(gpow(X,gdiv(gun,stoi(3)),(sz+3)/2));
 if (!gequal(gpow(G,stoi(3),-1),X)) return(gzero);
 if (n) G=gneg(G);
 return(G);
}

GEN getcurve(char *INP)
{GEN A;
 A=gp_read_str(INP); if (glength(A)==5) return(ellinit0(A,1,ELLACC));
 if (glength(A)==2)
   return(mineqfromc4c6(gmul(stoi(20736),(GEN) A[1]),
			gmul(stoi(2985984),(GEN) A[2])));
 return(gzero);
}

GEN jcurve(char *INP)
{GEN j,k;
 j=gp_read_str(INP); k=gsub(j,stoi(1728));
 if (gcmp0(j)) return(mineqfromc4c6(gzero,stoi(-864)));
 if (gcmp0(k)) return(mineqfromc4c6(stoi(-48),gzero));
 return(mineqfromc4c6(gdiv(gmul(stoi(144),j),k),gdiv(gmul(stoi(1728),j),k)));
}
