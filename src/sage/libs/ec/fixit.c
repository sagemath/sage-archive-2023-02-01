#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

char INP[128]; int DEAD=0; char *A;

//int my_getline() {A=fgets(INP,128,OUTFILE); if (NULL==A) return(0); return(1);}

int getcurrent()
{int i,k; char G[128]; int r,s; GEN H,T;
 k=my_getline(); if (!k) {DEAD=1; return(-1);}
 //printf("IINP1 %s",INP);
 sscanf(INP,"%s %s %s %s %s %s",G,G,G,G,G,G); if (G[0]=='+') return(-1);
 k=my_getline(); sscanf(INP,"%s",G); H=getcurve(G);
 //printf("IINP2 %s",INP);
 //printf("HI\n"); output(CURVE); output(H);
 //outbrute((GEN) H[11]); printf(" "); outbrute((GEN) CURVE[11]);
 if (gcmp0((GEN) H[10])) T=icbrt(gdiv((GEN) H[11],(GEN) CURVE[11]));
 else if (gcmp0((GEN) H[11])) T=racine(gdiv((GEN) H[10],(GEN) CURVE[10]));
 else
 T=gmul(gdiv((GEN) CURVE[10],(GEN) CURVE[11]),gdiv((GEN) H[11],(GEN) H[10]));
 if (typ(T)!=t_INT) T=gmul(T,gsqr(denom(T))); r=itos(T);
 //printf(" "); output(T);
 s=r%4; if (s<0) s+=4; if (s>1) r=r*4;
 if (CM==-1) {s=r%8; if (s<0) s+=8; if ((s&7)==4) r=-(r>>2);}
 for (i=1;i<NI;i++) k=my_getline(); if ((r%8)==0) for (i=0;i<=NI;i++) k=my_getline();
 return(r);
}

void fixtwists()
{int i; GEN CC; int ar; int DONE=0; int CURR=1;
 double LR; pari_sp memspot; GEN ib=stoi(100000000);

 if (mpcmp(ib,TWCOND)==-1) {i=my_getline(); return;}
 createpile(blocktwists()); i=0;
 while(!DONE)
 {memspot=avma;
  if ((pile[i].which&3)==3) pile[i].which=-pile[i].which;
  if ((pile[i].which&15)==4) pile[i].which=-pile[i].which;
  //printf("CURR %i %i\n",CURR,pile[i].which);
  if (CURR==pile[i].which) CURR=getcurrent();
  else if (pile[i].which)
  {CC=gdiv(gmul(TWCOND,gsqr(stoi(pile[i].which))),stoi(pile[i].condcorr));
   if ((mpcmp(ib,CC)==1) && (pile[i].cost<100000000000.0))
   {ar=qtwistar(pile[i].which,&LR,CC); twistdump(pile[i].which,ar,LR,CC);
    if (((pile[i].which&7)==0) && (CM!=-1))
    {ar=qtwistar(-pile[i].which,&LR,CC); twistdump(-pile[i].which,ar,LR,CC);}
   }
  }
  else DONE=1; i++; avma=memspot;
 }
}

void fixit()
{pari_sp memspot; int i,j,k; char *G; char *GS; int64 N;

 VERBOSE=0; G=malloc(128); GS=G; A=malloc(128); j=my_getline();
 while (!DEAD)
 {memspot=avma; G=GS; sscanf(INP,"%s %s %s %s %s %s",G,G,G,G,G,G);
 //printf("INP %s",INP);
  while (((G[0]<48) || (G[0]>57)) && (G[0]!='X')) G++;
  if (G[0]=='X') MODDEG=gzero; else MODDEG=gp_read_str(G);
  sscanf(INP,"%s %s %s %s %s",G,G,G,G,G); AISOG=atoi(G); ISOG=AISOG;
  sscanf(INP,"%s",G); N=atoll(G); NI=classsize();
  j=my_getline(); sscanf(INP,"%s",G); CURVE=getcurve(G);
  CL=cgetg(9,t_VEC); for(i=1;i<=8;i++)
  {CL[i]=lgetg(6,t_VEC); for(j=1;j<=5;j++) ((GEN) CL[i])[j]=(long) gzero;}

  if (((N%27)==0) && ((N%81)!=0) && (ISOG==9))
  {if (!gcmp0(icbrt((GEN) CURVE[12]))) ISOG=-9;}
  if (((N&31)==0) && ((N&63)!=0) && (ISOG==4)) ISOG=-4;

  if (anarray) free(anarray); if (antwarray) free(antwarray);
  antwsize=1; PREVD=0; ansize=1;
  anarray=malloc(u8((1+ansize)*sizeof(double))); anarray[1]=1.0;
  antwarray=malloc(u8((ansize+1)*sizeof(double))); antwarray[1]=1.0;

  ROOTNO=-1; TAMA=gun; COND=gun; PLACE=1;
  ISPRIME=0; ISSQFREE=0; ISSETZER=0; CM=0; CURVE=minimalequation();
  minimaltwist(); symsq(); if (!ISSQFREE) checkexotic();
  CL[1]=(long) ellinit0((GEN) CURVE,1,ELLACC);

  for (j=2;j<=NI;j++)
  {k=my_getline(); sscanf(INP,"%s",G); CL[j]=(long) gp_read_str(G);
   CL[j]=(long) ellinit0((GEN) CL[j],1,ELLACC);
  }
  ISSETZER=KnownDiffOptimal((GEN) CL[1]); fixtwists(); avma=memspot;
 }
 fclose(OUTFILE);
}
