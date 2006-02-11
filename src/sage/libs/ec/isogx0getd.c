#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

GEN createcycgenlist(GEN firstpow,int num)
{GEN rvec; int i;

 rvec=cgetg(num,t_VEC); rvec[1]=lcopy(firstpow);
 for (i=2;i<=num-1;i++) rvec[i]=(long) gmul((GEN) rvec[i-1],(GEN) rvec[1]);
 return(rvec);
}

GEN findincycgenlist(GEN find,GEN list,GEN p)
{int i;

 if (gequal(lift(find),gun)==1) return(gmodulcp(gzero,p));
 for (i=1;i<=glength(list);i++)
 {if(gequal(find,(GEN) list[i])==1) return(gmodulcp(stoi(i),p));}
 if(gequal(find,gmodulcp(stoi(3),stoi(8)))==1) return(gmodulcp(gzero,p));
 if(gequal(find,gmodulcp(stoi(5),stoi(8)))==1) return(gmodulcp(gun,p));
 if(gequal(find,gmodulcp(stoi(7),stoi(8)))==1) return(gmodulcp(gun,p));
 printf("ERROR in findincycgenlist\n");
 output(find); output(list);
 return(gzero);
}

GEN getd1(GEN zn,int p,GEN pm)
{GEN cycparts,pp,TEMP,cycgen,cycgenlist,powraise;
 GEN makevec,currpow,currmat,listofd,newmat,whichmodulus;
 byteptr pptr=diffptr; ulong pr=0; int i,k,newrank; int j=0; int currrank=0;

 currmat=gzero; cycparts=gcopy((GEN) zn[2]); pp=stoi(p);
 whichmodulus=cgetg(glength(cycparts)+1,t_VEC);
 powraise=cgetg(glength(cycparts)+1,t_VEC);
 cycgenlist=cgetg(glength(cycparts)+1,t_VEC);
 for (i=1;i<=glength(cycparts);i++)
 {TEMP=gdiventres((GEN) cycparts[i],pp);
  if (gequal((GEN) TEMP[2],gzero)==1)
  {j++;
   cycgen=gpow((GEN) ((GEN) zn[4])[i],(GEN) TEMP[1],-1);
   whichmodulus[j]=((GEN) zn[3])[i];
   cycgenlist[j]=(long) createcycgenlist(cycgen,p);
   powraise[j]=lcopy((GEN) TEMP[1]);
  }
 }
 listofd=cgetg(j+1,t_VEC); i=1;
 while(1)
 {pr+=*pptr++; pp=stoi(pr);
  if (gequal(gmod(COND,pp),gzero)!=1)
  {makevec=cgetg(j+1,t_VEC);
   for (k=1;k<=j;k++)
   {currpow=gpow(gmodulcp(pp,(GEN) whichmodulus[k]),(GEN) powraise[k],-1);
    if (mod8((GEN) whichmodulus[k])==0)
    {if ((k==1) && (p==2))
     {if (mod4(pp)==3) makevec[2]=(long) gmodulcp(gun,gdeux);
      else makevec[2]=(long) gmodulcp(gzero,gdeux);
      makevec[k]=(long) findincycgenlist(currpow,(GEN) cycgenlist[k],pm);
     }
     else if (k!=2)
     {if (mod4(lift(currpow))==3) currpow=gneg(currpow);
      makevec[k]=(long) findincycgenlist(currpow,(GEN) cycgenlist[k],pm);
     }
    }
    else makevec[k]=(long) findincycgenlist(currpow,(GEN) cycgenlist[k],pm);
   }
   if (currrank==0)
   {newmat=cgetg(j+1,t_MAT);
    for (k=1;k<=j;k++)
    {newmat[k]=lgetg(2,t_VEC); ((GEN) newmat[k])[1]=makevec[k];}
   }
   else newmat=concat(currmat,makevec);
   newrank=rank(newmat);
   if (newrank!=currrank)
   {currrank++; listofd[currrank]=(long) pp; currmat=newmat;}
  }
  i++; if (currrank==j) break;
 }
 return(listofd);
}

