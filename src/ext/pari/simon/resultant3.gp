\\ Ce fichier gp contient des fonctions pour calculer
\\ le resultant de trois polynomes p1,p2,p3 homogenes
\\ en trois variables (toujours x,y,z), ainsi que
\\ le discriminant d'un polynome homogene en ces trois variables.
\\ L'algorithme utilise est celui du sous-resultant.
\\
\\ ***********************************************
\\ *       Auteur: Denis SIMON                   *
\\ *       mail: desimon@math.unicaen.fr.fr      *
\\ *       date: 12 Dec 2003                     *
\\ ***********************************************
\\
\\ exemple d'utilisation :
\\
\\ ? p1=x^2-3*z^2+y*z;p2=x-y+15*z;p3=y^2*x+z^3-x*y*z+x^2*z;
\\ ? resultant3([p1,p2,p3])
\\ %2 = 521784
\\ ? discriminant3(p3)
\\ %3 = -63
\\
\\ la fonction hom sert a rendre homogene un polynome en x et y:
\\ ? ell=y^2-y-x^3+x^2;
\\ ? discriminant3(hom(ell))
\\ %5 = -11
\\
\\

{hom(p)=
p=subst(subst(z^100*p,x,x/z),y,y/z);
p/=z^myvaluation(p,z);
return(p);
}
\\
\\ degre d'un polynome homogene.
\\
{degreetot(p)=
local(auxdeg,auxp);
auxdeg=poldegree(p,x);auxp=pollead(p,x);
auxdeg+=poldegree(auxp,y);auxp=pollead(auxp,y);
auxdeg+=poldegree(auxp,z);auxp=pollead(auxp,z);
return(auxdeg);
}
\\
{ex(vec,i,j)=
local(aux);
aux=vec[i];vec[i]=vec[j];vec[j]=aux;
return(vec);
}
{myvaluation(p,var)=
local(valp);
valp=0;
while(polcoeff(p,valp,var)==0,valp++);
return(valp);
}
{myvaluation2(p,var)=
local(pp,valp);
pp=p;
valp=0;
while(subst(pp,var,0)==0,valp++;pp/=var);
return(valp);
}
{mycontent(p)=
local(vz,co,dco);
vz=myvaluation(p,z);
p=subst(p,z,1);
co=content(p);
dco=poldegree(co,y);
if(dco,co=subst(co,y,y/z)*z^dco);
return(co*z^vz);
}
{resultant2(p,q)=
local(dp,dq,valp,valq,auxr,auxp,auxq,res);
if(p==0 || q==0,return(0));
dp=degreetot(p);dq=degreetot(q);
valp=myvaluation(p,y);
valq=myvaluation(q,y);
if(valp && valq,return(0));
auxr=1;
if(valp,
  if(dq%2 && valp%2,auxr=-1);
  auxr*=pollead(q,x)^valp);
if(valq,auxr=pollead(p,x)^valq);
auxp=subst(p,y,1);
auxq=subst(q,y,1);
res=auxr*polresultant(auxp,auxq,x);
return(res);
}
{resultantcomp(vp=[x,y,z])=
\\ les vp[i] sont des pol homogenes en x,y et z.
\\ vp[3] ne depend que de y et z.
local(p,q,rt,dp,dq,drt,vy,vz,lrt,lp,res,aux,res2,dres2);
p=vp[1];q=vp[2];rt=vp[3];
if(p==0 || q==0 || rt==0,return(0));
dp=degreetot(p);dq=degreetot(q);drt=degreetot(rt);
if(drt==0,return(rt^(dp*dq)));
vy=myvaluation(rt,y);vz=myvaluation(rt,z);
rt=subst(rt,y,1);if(vz,rt/=z^vz);
lrt=polcoeff(rt,drt,z);lp=polcoeff(p,dp,x);
res=1;
if(lp==0,
  aux=p;p=q;q=aux;
  aux=dp;dp=dq;dq=aux;
  if(dp%2 && dq%2 && drt%2,res=-res);
  lp=polcoeff(p,dp,x);
  if(lp==0,return(0)));
if(vy,
  if(dp%2 && dq%2, res2=-1,res2=1);
  res2*=resultant2(subst(subst(p,y,0),z,y),subst(subst(q,y,0),z,y));
  if(res2==0,return(0));
  res*=res2^vy);
if(vz,
  res2=resultant2(subst(p,z,0),subst(q,z,0));
  if(res2==0,return(0));
  res*=res2^vz);
drt-=(vy+vz);
if(drt==0,return(res*rt^(dp*dq)));
res2=polresultant(subst(p,y,1),subst(q,y,1),x);
dres2=poldegree(res2,z);
res2=polresultant(rt,res2,z);
res*=res2;
if(dq!=poldegree(subst(q,y,1),x),res*=lp^(drt*(dq-poldegree(subst(q,y,1),x))));
if(dres2!=dp*dq,res*=pollead(rt,z)^(dp*dq-dres2));
return(res);
}
{resultant3(vp=[x,y,z],fl)=
\\ resultant de 3 polynomes homogenes en x,y,z.
\\ fl=1 pour afficher des resultats intermediaires.
\\ on travaille sur la variable x.
local(vdt,vdx,prodd,s,denom,nume,lp,p0,q0,r0,dd,cq,rm,res);

for(i=1,3,if(vp[i]==0,return(0)));
vdt=vector(3,i,degreetot(vp[i]));
vdx=vector(3,i,poldegree(vp[i],x));
if(vecmin(vdt-vdx),return(0));
\\ en effet, dans ce cas, les polynomes sont dans l'ideal <y,z>

prodd=prod(i=1,3,vdt[i]);
s=1;
denom=1;nume=1;
\\ on echange pour que vdx[1] >= vdx[2] >= vdx[3]
if(vdx[1]<vdx[2] || (vdx[1]==vdx[2] && vdt[1]<vdt[2]),
  vp=ex(vp,1,2);vdx=ex(vdx,1,2);vdt=ex(vdt,1,2);if(prodd%2,s=-s));
if(vdx[1]<vdx[3] || (vdx[1]==vdx[3] && vdt[1]<vdt[3]),
  vp=ex(vp,1,3);vdx=ex(vdx,1,3);vdt=ex(vdt,1,3);if(prodd%2,s=-s));
if(vdx[2]<vdx[3] || (vdx[2]==vdx[3] && vdt[2]<vdt[3]),
  vp=ex(vp,2,3);vdx=ex(vdx,2,3);vdt=ex(vdt,2,3);if(prodd%2,s=-s));

\\ cas particulier ou vp[3] est constant
if(vdt[3]==0,return(vp[3]^(vdt[1]*vdt[2])*s));

\\ on fait un echange pour que vp[1] contienne le monome x^vdt[1]
lp=polcoeff(vp[1],vdt[1],x);
if(lp==0,
  lp=polcoeff(vp[2],vdt[2],x);
  if(lp==0,
    lp=polcoeff(vp[3],vdt[3],x);
    if(lp==0,return(0));
    vp=ex(vp,1,3);vdx=ex(vdx,1,3);vdt=ex(vdt,1,3);if(prodd%2,s=-s)
  , vp=ex(vp,1,2);vdx=ex(vdx,1,2);vdt=ex(vdt,1,2);if(prodd%2,s=-s)));

\\ on supprime de vp[1] les puissances de x
p0=polcoeff(vp[1],0,x);q0=polcoeff(vp[2],0,x);r0=polcoeff(vp[3],0,x);
while(p0==0,
  vp[1]/=x;vdx[1]-=1;vdt[1]-=1;p0=polcoeff(vp[1],0,x);
  nume*=resultant2(subst(subst(q0,y,x),z,y),subst(subst(r0,y,x),z,y)));

if(nume==0,return(0));
dd=vecmax(vdx)+1;
while(vdx[3],

  if(fl,print("vp="vp));
  if(fl,print("vdx="vdx));
  if(fl,print("vdt="vdt));
  if(fl,print([nume,denom]));

  if(denom==0,print(" ERREUR ");return(x));
  if(nume==0,return(0));
  for(i=2,3,if(vp[i]==0,return(0)));

  q0=polcoeff(vp[2],0,x);
  if(q0==0,  \\ si vp[2] est divisible par x
    vp[2]/=x;vdx[2]-=1;vdt[2]-=1;
    nume*=resultant2(subst(subst(r0,y,x),z,y),subst(subst(p0,y,x),z,y));next);
  r0=polcoeff(vp[3],0,x);
  if(r0==0,  \\ si vp[3] est divisible par x
    vp[3]/=x;vdx[3]-=1;vdt[3]-=1;
    nume*=resultant2(subst(subst(p0,y,x),z,y),subst(subst(q0,y,x),z,y));next);

  \\ on enleve le contenu
  cq=mycontent(vp[2]);
  if(cq!=1,
    nume*=resultantcomp([vp[3],vp[1],cq]);
    vp[2]/=cq;vdt[2]-=poldegree(cq);next);
  cq=mycontent(vp[3]);
  if(cq!=1,
    nume*=resultantcomp([vp[1],vp[2],cq]);
    vp[3]/=cq;vdt[3]-=poldegree(cq);next);

  rm=pollead(vp[3],x);\\ le coefficient dominant de vp[3].
  res=resultantcomp([vp[1],vp[3]-rm*x^vdx[3],rm]);
  if(res==0,print1("cas non traite");return(x));\\ ce cas est-il possible ?
  if(vdt[1]%2 && vdt[3]%2 && degreetot(rm)%2, res=-res);
  \\ on baisse le degre de vp[2] en retranchant des multiples de vp[3].
  while(vdx[3]<=vdx[2],
    delta=vdx[2]-vdx[3];
    lp=pollead(vp[2],x);
    vp[2]=simplify(rm*vp[2]-lp*vp[3]*x^delta);
    vdx[2]=poldegree(vp[2],x);
    denom*=res
  );
  vdt[2]=degreetot(vp[2]);

  prodd=prod(i=1,3,vdt[i]);
  vp=ex(vp,2,3);vdx=ex(vdx,2,3);
  vdt=ex(vdt,2,3);
  if(prodd%2,s=-s);
);
vp[3]=polcoeff(vp[3],0,x);
nume*=resultantcomp(vp);
if(fl,print([nume,denom]));
return(simplify(s*nume/denom));
}
\\
\\ discriminant d'un pol homogene en 3 variables x, y et z.
\\
{
discriminant3(p,fl)=
local(dp,normal,re);
dp=degreetot(p);
normal=dp^(dp^2-3*dp+3)*(-1)^(dp*(dp-1)/2);
re=resultant3([deriv(p,x),deriv(p,y),deriv(p,z)],fl);
if(re==x,return(0));
return(re/normal);
}






