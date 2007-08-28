#include  "lie.h"
#define local  static
/* #include  <string.h> */
#define two_lengths(type) (strchr("BCFG",type)!=NULL)
#define grp_less(x,y) \
 ((x)->lietype<(y)->lietype \
      || (x)->lietype==(y)->lietype && (x)->lierank<(y)->lierank)
#define opposite(a,b) adj[b][adj[b][0]==(a)]


local void fundam(matrix* roots, index first, index* last);


local index s; /* the current semisimple rank */


local simpgrp* simp_type(entry** m, entry n)
{ matrix* adjs=mkmatrix(n,3);
  entry** adj=adjs->elm /* |adj[i]| lists up to 3 neighbours of node |i| */
  ,* norm=mkintarray(3*n) /* norms of roots */
  ,* valency=&norm[n] /* valencies in Dynkin diagram */
  ,* p=&valency[n]; /* permutation of |n| */
  simpgrp* result;
  index i,j,k, a_val[4]={-1,-1,-1,-1};
    /* |a_val[i]| is index of a node of valency |i|, if any */

  if (n==0) error("empty input in simp_type\n");

  { for (i=0;i<n;i++) valency[i]=0;
      /* |valency[i]| is also index of next slot in |adj[i]| */
    for (i=n; --i>=0;)
    { norm[i]= Norm(m[i]); /* where |Norm(x)==Inprod(x,x)/2| */
      for (j=i; --j>=0;)
        if (Inprod(m[i],m[j])!=0) /* then valencies increase */
        { if (valency[i]>=3 || valency[j]>=3) error ("valency >3 found\n");
  	adj[i][valency[i]++]=j; adj[j][valency[j]++]=i;
  	/* update valencies and adjacencies */
        }
      a_val[valency[i]]=i; /* valency of node |i| is now known */
    }
  }
  if (a_val[3]<0)

  { index e; /* index of end node (|valency[e]<=1|) */
    if (a_val[0]>=0) p[0]=e=a_val[0]; /* must be type $A_1$ */
    else
    { if (a_val[1]>=0) p[0]=e=a_val[1]; /* other linear types */
      else error("no end node found\n");


      { k=p[1]=adj[e][0]; /* the unique neighbour of node |e| */
        for(i=2;i<n;i++)  p[i]=k=opposite(p[i-2],k); /* here |k==p[i-1]| */
      }

      if ( n==2 && norm[p[0]]+2*norm[p[1]]==5
        || n>=3 && norm[p[0]]!=norm[p[1]]
        || n==4 && norm[p[1]]<norm[p[2]]
         )
      { for (i=0; i<n-1-i; i++) swap(&p[i],&p[n-1-i]); e=p[0]; }
    }

    { entry norm0=norm[p[0]], norm1=norm[p[n-1]];
      if (norm0==norm1) result = mksimpgrp('A',n);
      else if (norm1==3) result=mksimpgrp('G',2);
      else if (norm1==2) result=mksimpgrp('C',n);
      else if (norm0!=2) error("I don't recognize this Cartan Type\n");
      else if (n==4 && norm[p[2]]==1) result=mksimpgrp('F',4);
      else result=mksimpgrp('B',n);
    }
  }
 /* no nodes of valency 3 */
  else
       { entry* branch=adj[a_val[3]], end[3], end_count=0;
         for (j=2; j>=0; j--)
           if (valency[branch[j]]==1) end[end_count++]=branch[j];
         if (end_count>1)
                         { p[n-1]=end[1]; p[n-2]=end[0]; p[n-3]=a_val[3];
                           k=p[n-4]=branch[0]+branch[1]+branch[2]-p[n-1]-p[n-2];
                             /* the remaining branch */
                           for(i=n-5; i>=0; i--)
                           { if (valency[k]!=2) error("unlinear Dn tail.\n");
                             p[i]=k=opposite(p[i+2],k);
                           }
                           result=mksimpgrp('D',n);
                         }
         else if (end_count==1)
                              { p[3]=a_val[3]; p[1]=end[0];
                                for (j=2; j>=0; j--)
                                  if (valency[branch[j]]==2)
                                    if (valency[opposite(a_val[3],branch[j])]==1) break;
                                if (j<0) error("type E not recognised\n");
                                p[2]=branch[j]; p[0]=opposite(p[3],p[2]);
                                p[4]=k=branch[0]+branch[1]+branch[2]-p[1]-p[2]; /* remaining branch */
                                for(i=5;i<n;i++)
                                { if (valency[k]!=2) error("wrong type E system.\n");
                                  p[i]=k=opposite(p[i-2],k);
                                }
                                result=mksimpgrp('E',n);
                              }
         else error("no end node adjacent to valency 3 node\n");
       }

  for (i=0; i<n; i++) if (p[i]>=0) /* then |p[i]| starts an untreated cycle */
  { entry* mi=m[j=i]; /* record beginning of cycle */
    while (p[j]!=i)
      { k=j; j=p[j]; m[k]=m[j]; p[k]= -1; }
        /* assign |m[j]=m[p[j]]| and advance */
    m[j]=mi; p[j]= -1; /* close the cycle */
  }
  freemem(adjs); freearr(norm);
  return result;
}

static void cycle_block(matrix* m, index first, index last, index amount)
{ entry** row=&m->elm[first];  index modulus=last-first,i,min_done=amount;
  if (amount>0 && modulus>amount) /* otherwise there is nothing to do */
    for (i=0; i<min_done; ++i) /* |min_done| is fixed after first time round */
    { entry* row_i=row[i]; index j=i+amount,old_j=i;
      do /* perform cycle containing |i| */
      { row[old_j]=row[j]; old_j=j;
	if ((j+=amount)>=modulus) j-=modulus; /* wrap downwards */
      } while (j>=min_done || j>i && (min_done=j,true));
      row[old_j]=row_i; /* close the cycle */
      if (j!=i) error("System error cycling.\n");
    }
}

local void long_close(matrix* m, index first, index last)
{ index i,j;  entry* root_i,* root_j,* t=mkintarray(s);
  for (i=first; i<last; ++i)
  { root_i=m->elm[i]; if (Norm(root_i)>1) continue;
    for (j=i+1; j<last; ++j)
    { root_j=m->elm[j]; if (Norm(root_j)>1) continue;
      subrow(root_i,root_j,t,s);
      if (isroot(t))
	if (isposroot(t))
	{ copyrow(t,root_i,s); break;  } /* need not consider more |j|'s */
	else add_xrow_to(root_j,-1,root_i,s);
    }
  }
  freearr(t);
}


matrix* Closure(matrix* m, boolean close, group* lie_type)
{ matrix* result;  index i,j;
  group* tp=(s=Ssrank(grp), lie_type==NULL ? mkgroup(s) : lie_type);

  tp->toraldim=Lierank(grp); tp->ncomp=0; /* start with maximal torus */
  m=copymatrix(m);

  if (close)
    if (type_of(grp)==SIMPGRP) close = two_lengths(grp->s.lietype);
    else
    { for (i=0; i<grp->g.ncomp; i++)
        if (two_lengths(Liecomp(grp,i)->lietype)) break;
      close= i<grp->g.ncomp;
    }

  { entry* t;
    for (i=0; i<m->nrows; i++)
    if (!isroot(t=m->elm[i]))
      error("Set of root vectors contains a non-root\n");
    else if (!isposroot(t=m->elm[i]))
      for (j=0; j<m->ncols; j++) t[j]= -t[j]; /* make positive root */
    Unique(m,cmpfn);
  }
  { index next;
    for (i=0; i<m->nrows; i=next)

    { index d,n=0;  simpgrp* c;
      next=isolcomp(m,i);
      fundam(m,i,&next);
      if (close) long_close(m,i,next),fundam(m,i,&next);
      c=simp_type(&m->elm[i],d=next-i);

      { j=tp->ncomp++;
        while(--j>=0 && grp_less(tp->liecomp[j],c))
          n += (tp->liecomp[j+1]=tp->liecomp[j])->lierank;
        tp->liecomp[++j]=c; tp->toraldim -= d;
          /* insert component and remove rank from torus */
        cycle_block(m,i-n,next,n);
          /* move the |d| rows down across |n| previous rows */
      }
    }
  }
  if (lie_type==NULL)
    return result=copymatrix(m),freemem(m),freemem(tp),result;
  else return freemem(m),(matrix*)NULL; /* |Cartan_type| doesn't need |m| */
}

group* Carttype(matrix* m)
{ group* type=mkgroup(s=Ssrank(grp)); /* rank bounds number of components */
  Closure(m,false,type); return type;
}

local void fundam(matrix* roots, index first, index* last)
{ index i,j,d;  boolean changed;
  entry* t=mkintarray(s);  matrix mm,* m=&mm;
  mm.elm=&roots->elm[first]; mm.nrows=*last-first; mm.ncols=roots->ncols;
  for (i=m->nrows-1; i>0; changed ? Unique(m,cmpfn),i=m->nrows-1 : --i)
  { entry* root_i=m->elm[i]; changed=false;

    for (j=i-1; j>=0; j--)
    { entry* root_j=m->elm[j]; entry c=Cart_inprod(root_j,root_i);
      if (c==2 && eqrow(root_j,root_i,s))

        { cycle_block(m,j,m->nrows--,1); root_i=m->elm[--i];
        }
      else if (c>0)
      { changed=true;

        { copyrow(root_j,t,s); add_xrow_to(t,-c,root_i,s);
          if (isposroot(t)) copyrow(t,root_j,s);
          else
          { j=i; c=Cart_inprod(root_i,root_j);
            copyrow(root_i,t,s); add_xrow_to(t,-c,root_j,s);
            if (isposroot(t)) copyrow(t,root_i,s);
            else
                 { index k;  entry* ln,* sh; /* the longer and the shorter root */
                   if (Norm(root_i)>Norm(root_j))
                     ln=root_i, sh=root_j;  else ln=root_j, sh=root_i;
                   switch (Norm(ln))
                   { case 2: subrow(ln,sh,sh,s); /* |sh=ln-sh| */
                     add_xrow_to(ln,-2,sh,s); /* |ln=ln-2*sh| */
                   break; case 3: /* |grp=@t$G_2$@>| now */
                     for (k=0; sh[k]==0; ++k) {} /* find the place of this $G_2$ component */
                     sh[k]=1; sh[k+1]=0; ln[k]=0; ln[k+1]=1;
                       /* return standard basis of full system */
                   break; default: error("problem with norm 1 roots\n");
                   }
                 }
          }
        }
      }
    }
  }
  cycle_block(roots,first+mm.nrows,roots->nrows,d=*last-first-mm.nrows);
*last-=d; roots->nrows-=d;
  freearr(t);
}

matrix* Resmat(matrix* roots)
{ index i,j,k,r=Lierank(grp),s=Ssrank(grp), n=roots->nrows;
  vector* root_norms=Simproot_norms(grp);
  entry* norms=root_norms->compon;
    /* needed to compute $\<\lambda,\alpha[i]>$ */
  matrix* root_images=Matmult(roots,Cartan()),* result=mkmatrix(r,r);
  entry** alpha=roots->elm,** img=root_images->elm,** res=result->elm;

  for (i=0; i<r; i++) for (j=0; j<r; j++) res[i][j]= i==j;
    /* initialise |res| to identity */
  for (j=0; j<n; j++) /* traverse the given roots */

  { entry* v=img[j], norm=(checkroot(alpha[j]),Norm(alpha[j]));
    for (k=s-1; v[k]==0; k--) {}
    if (k<j)
      error("Given set of roots is not independent; apply closure first.\n");

    if (v[k]<0)
    { for (i=j; i<n; i++) img[i][k]= -img[i][k];
      for (i=k-j; i<s; i++) res[i][k]= -res[i][k];
    }
    while(--k>=j)
      /* clear |v[k+1]| by unimodular column operations with column~|j| */
    {
        entry u[3][2];  index l=0;
        u[0][1]=u[1][0]=1; u[0][0]=u[1][1]=0;
        u[2][1]=v[k]; u[2][0]=v[k+1];
        if (v[k]<0) u[2][1]= -v[k], u[0][1]= -1; /* make |u[2][1]| non-negative */
        do /* subtract column |l| some times into column |1-l| */
        { entry q=u[2][1-l]/u[2][l];  for (i=0; i<3; i++) u[i][1-l]-=q*u[i][l];
        } while (u[2][l=1-l]!=0);
        if (l==0)  for (i=0; i<2; i++) swap(&u[i][0],&u[i][1]);

      { for (i=j; i<n; i++) /* combine columns |k| and |k+1| */
        { entry img_i_k=img[i][k];
          img[i][k]  =img_i_k*u[0][0]+img[i][k+1]*u[1][0];
          img[i][k+1]=img_i_k*u[0][1]+img[i][k+1]*u[1][1];
        }
        for (i=k-j; i<s; i++)
        { entry res_i_k=res[i][k];
          res[i][k]=res_i_k*u[0][0]+res[i][k+1]*u[1][0];
          res[i][k+1]=res_i_k*u[0][1]+res[i][k+1]*u[1][1];
        }
      }
     }
    for (i=0; i<s; i++)
                    { index inpr= norms[i]*alpha[j][i]; /* this is $(\omega_i,\alpha[j])$ */
                      if (inpr%norm!=0) error("Supposed root has non-integer Cartan product.\n");
                      res[i][j]=inpr/norm; /* this is $\<\omega_i,\alpha[j]>$ */
                    }
  }
  freemem(root_norms); freemem(root_images); return result;
}

