#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void readfile(char *filename)
{FILE *WORKFILE; int i,j,query,wp,loc,mds; unsigned char *wd;
 int twdatsize;

 WORKFILE=fopen(filename,"r"); fread(headbytes,1,8,WORKFILE);
 residues=
   (((int) headbytes[0])<<16)+(((int) headbytes[1])<<8)+((int) headbytes[2]);
 c6r=residues%1728; c4r=(residues-c6r)/1728; c4neg=0; c6neg=0;
 c4s=(int) headbytes[3]; if (c4s>128) {c4s=256-c4s; c4neg=1;}
 c6s=(int) headbytes[4]; if (c6s>128) {c6s=256-c6s; c6neg=1;}
 numcurves=
   (((int) headbytes[5])<<16)+(((int) headbytes[6])<<8)+((int) headbytes[7]);
 recsize=c4s+c6s; cdatsize=numcurves*recsize;
 curvedata=malloc(u8(cdatsize+4*recsize));
 fread(curvedata,1,cdatsize,WORKFILE); query=fread(backbytes,1,4,WORKFILE);
 if (query==0) {twdatsize=0; wd=malloc(0);}
 else twdatsize=(((int) backbytes[0])<<24)+(((int) backbytes[1])<<16)+
	(((int) backbytes[2])<<8)+((int) backbytes[3]);
 wd=malloc(u8(twdatsize)); fread(wd,1,twdatsize,WORKFILE); fclose(WORKFILE);
 size16bytes=malloc(u8((numcurves+10)>>2));
 size256bytes=malloc(u8((numcurves+10)>>6));
 size4096bytes=malloc(u8((numcurves+10)>>10));
 leadbyte=malloc(u8(numcurves+16)); moddegbytes=malloc(u8(8*numcurves+80));
 localsizebytes=malloc(u8(4*numcurves+80)); twdata=malloc(u8(4*numcurves+80));
 if (!twdatsize)
 {for (i=0;i<numcurves;i++)
  {leadbyte[i]=0; moddegbytes[8*i]=0; moddegbytes[8*i+1]=0;
   moddegbytes[8*i+2]=0; moddegbytes[8*i+3]=0;
   moddegbytes[8*i+4]=0; moddegbytes[8*i+3]=0;
   moddegbytes[8*i+6]=0; moddegbytes[8*i+5]=0;
   localsizebytes[4*i]=0; localsizebytes[4*i+1]=0;
   localsizebytes[4*i+2]=0; localsizebytes[4*i+3]=0; twdata[i]=malloc(0);
  }
 }
 else
 {wp=0; for (i=0;i<numcurves;i++)
  {if (((i&4095)==0) && (i+4096<numcurves))
   {size4096bytes[(i>>10)+2]=wd[wp+2]; size4096bytes[(i>>10)+3]=wd[wp+3];
    size4096bytes[i>>10]=wd[wp]; size4096bytes[(i>>10)+1]=wd[wp+1]; wp+=4;
   }
   if (((i&255)==0) && (i+256<numcurves))
   {size256bytes[i>>6]=wd[wp]; size256bytes[(i>>6)+1]=wd[wp+1];
    size256bytes[(i>>6)+2]=wd[wp+2]; size256bytes[(i>>6)+3]=wd[wp+3]; wp+=4;
   }
   if (((i&15)==0) && (i+16<numcurves))
   {size16bytes[i>>2]=wd[wp]; size16bytes[(i>>2)+1]=wd[wp+1];
    size16bytes[(i>>2)+2]=wd[wp+2]; size16bytes[(i>>2)+3]=wd[wp+3]; wp+=4;
   }
   leadbyte[i]=wd[wp]; wp++; localsizebytes[4*i]=wd[wp]; loc=wd[wp]; wp++;
   if (localsizebytes[4*i]>=240)
   {localsizebytes[4*i+1]=wd[wp]; loc-=240; loc<<=8; loc+=wd[wp]; wp++;}
   if (localsizebytes[4*i]==255)
   {localsizebytes[4*i+2]=wd[wp]; loc-=3840; loc<<=8; loc+=wd[wp]; wp++;
    if (localsizebytes[4*i+1]>=240)
    {localsizebytes[4*i+3]=wd[wp]; loc-=61440; loc<<=8; loc+=wd[wp]; wp++;}
   }
   mds=(leadbyte[i]&63)>>3;
   for (j=0;j<mds;j++) moddegbytes[8*i+j]=wd[wp+j]; wp+=mds;
   twdata[i]=malloc(u8(loc));
   for (j=0;j<loc;j++) twdata[i][j]=wd[wp+j]; wp+=loc;
  }
 }
 free(wd);
}

void writefile(char *filename)
{int wp,s,k,j,cliff,loc,mds,u; unsigned char fb[4]; FILE *WORKFILE;
 unsigned char *wd;

 wp=8; s=12+numcurves*(recsize+10);
 for(k=0;k<numcurves;k++) s+=uwlsb(localsizebytes,k); wd=malloc(u8(s));
 for (k=0;k<8;k++) wd[k]=headbytes[k];
 for (k=0;k<cdatsize;k+=recsize)
 {for (j=0;j<recsize;j++) {wd[wp+j]=curvedata[k+j];} wp+=recsize;}
 cliff=wp; wp+=4;
 for(k=0;k<numcurves;k++)
 {if ((k&4095)==0)
  {j=k>>10; if (k!=0)
   {u=gettheint(size4096bytes,j-4); get4bytes(wp-u,fb);
    wd[u]=fb[0]; wd[u+1]=fb[1]; wd[u+2]=fb[2]; wd[u+3]=fb[3];
   }
   if (k+4096<numcurves)
   {get4bytes(wp,fb); size4096bytes[j]=fb[0]; size4096bytes[j+1]=fb[1];
    size4096bytes[j+2]=fb[2]; size4096bytes[j+3]=fb[3]; wp+=4;
   }
  }
  if ((k&255)==0)
  {j=k>>6; if (k!=0)
   {u=gettheint(size256bytes,j-4); get4bytes(wp-u,fb);
    wd[u]=fb[0]; wd[u+1]=fb[1]; wd[u+2]=fb[2]; wd[u+3]=fb[3];
   }
   if (k+256<numcurves)
   {get4bytes(wp,fb); size256bytes[j]=fb[0]; size256bytes[j+1]=fb[1];
    size256bytes[j+2]=fb[2]; size256bytes[j+3]=fb[3]; wp+=4;
   }
  }
  if ((k&15)==0)
  {j=k>>2; if (k!=0)
   {u=gettheint(size16bytes,j-4); get4bytes(wp-u,fb);
    wd[u]=fb[0]; wd[u+1]=fb[1]; wd[u+2]=fb[2]; wd[u+3]=fb[3];
   }
   if (k+16<numcurves)
   {get4bytes(wp,fb); size16bytes[j]=fb[0]; size16bytes[j+1]=fb[1];
    size16bytes[j+2]=fb[2]; size16bytes[j+3]=fb[3]; wp+=4;
   }
  }
  wd[wp]=leadbyte[k]; wp++; wd[wp]=localsizebytes[4*k]; loc=wd[wp]; wp++;
  if (localsizebytes[4*k]>=240)
  {wd[wp]=localsizebytes[4*k+1]; loc-=240; loc<<=8; loc+=wd[wp]; wp++;}
  if (localsizebytes[4*k]==255)
  {wd[wp]=localsizebytes[4*k+2]; loc-=3840; loc<<=8; loc+=wd[wp]; wp++;
   if (localsizebytes[4*k+1]>=240)
   {wd[wp]=localsizebytes[4*k+3]; loc-=61440; loc<<=8; loc+=wd[wp]; wp++;}
  }
  mds=(leadbyte[k]&63)>>3;
  for (j=0;j<mds;j++) wd[wp+j]=moddegbytes[8*k+j]; wp+=mds;
  for (j=0;j<loc;j++) wd[wp+j]=twdata[k][j]; wp+=loc;
 }
 s=wp-12-cdatsize; wd[cliff]=s>>24; wd[cliff+1]=(s>>16)&255;
 wd[cliff+2]=(s>>8)&255; wd[cliff+3]=s&255;
 WORKFILE=fopen(filename,"w"); fwrite(wd,1,wp,WORKFILE); fclose(WORKFILE);
 free(wd);
}
