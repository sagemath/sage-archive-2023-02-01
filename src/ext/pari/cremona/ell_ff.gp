/* Functions for working in the function field K=k(u,v) of an elliptic curve
 (over an arbitrary ground field k).

Author: John Cremona (john.cremona@nottingham.ac.uk)
Last edited: 13/10/04

 Here
          v^2+a1*u*v+a3*v-(u^3+a2*u^2+a4*u+a6)=0.

 So [u,v] is a generic point on E.

 The generic point is obtained by calling ellgenpt(e) when needed;
 we could store this and pass it to the functions which need it, but
 currently it is recomputed as needed, which should be essentially
 cost-free.

requires ell_utils.gp

*/

global(ellffV,ellffU);

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ The polynomial defining the curve (in universal global variables)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffeqn(e)=
local(U,V);  \\ for readability
U=ellffU; V=ellffV;
V^2+e.a1*U*V+e.a3*V-(U^3+e.a2*U^2+e.a4*U+e.a6)
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ The generic point on the curve (in terms of universal global variables)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffgenpt(e)=
local(eqn);
eqn=ellffeqn(e);
[Mod(ellffU,eqn),Mod(ellffV,eqn)]
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ evaluation at a point (user's responsibility to avoid poles!)
\\ and at a divisor (p1)-(p2) (ditto)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffeval(h,p)=
if(p==[0],print("Cannot evaluate functions at [0]");return(0));
if(h==0,0,subst(subst(lift(h),ellffV,p[2]),ellffU,p[1]))
}

{
ellffeval2(h,p1,p2)=ellffeval(h,p1)/ellffeval(h,p2);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ elffvert(e,t) returns a function h whose divisor is (t)+(-t)-2(0),
\\           i.e. the vertical line through t
\\ if t==[0] this degenerates to 1
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffvert(e,t)=
local(eqn);
eqn=ellffeqn(e);
if(t==[0], Mod(1,eqn), Mod(ellffU-t[1],eqn));
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellfftang(e,t) returns a function h whose divisor is 2(t)+(-2t)-3(0)
\\             i.e. the tangent to e at t
\\ if t has order 2 this degenerates to ellffvert(e,t)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellfftang(e,t)=
local(gpt,dyf,dxf,u,v);
gpt=ellffgenpt(e);
if(t==[0],return(Mod(1,gpt[1].mod)));
u=gpt[1]; v=gpt[2];
dyf=2*t[2]+e[1]*t[1]+e[3];
if(dyf==0,return(u-t[1]));
dxf=e[1]*t[2]-(3*t[1]^2+2*e[2]*t[1]+e[4]);
v-t[2]+(dxf/dyf)*(u-t[1]);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellffchord(e,t,s) returns a function h whose divisor is
\\     (t)+(s)+(-s-t)-3(0)
\\ i.e. the chord through t and s
\\ if s=t  this degenerates to ellfftang(e,t)
\\ if s=-t this degenerates to ellffvert(e,t)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffchord(e,t,s)=
local(gpt);
if(s==[0],return(ellffvert(e,t)));
if(t==[0],return(ellffvert(e,s)));
if(s==t,return(ellfftang(e,t)));
if(s[1]==t[1],return(ellffvert(e,t))); \\s==-t!=t here
\\ Now s[1]!=t[1]:
gpt=ellffgenpt(e);
gpt[2]-t[2]-((t[2]-s[2])/(t[1]-s[1]))*(gpt[1]-t[1]);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellffaddpol(e,t,s) returns [s+t,h] where h's divisor is
\\     (t)+(s)-(s+t)-(0)
\\ i.e. the function "proving" s+t is correct
\\ if s=-t this degenerates to ellffvert(e,t)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffaddpol(e,t,s)=
local(st);st=elladd(e,s,t);
return([st,ellffchord(e,s,t)/ellffvert(e,st)]);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellffmultpol(e,t,m) : recursive function for m>0
\\ returns a function h whose divisor is
\\     m(t)-(m*t)-(m-1)(0)
\\ i.e. the function "proving" m*t is correct
\\ for recursive efficiency also returns m*t
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffmultpol(e,t,m)=
local(m2,m2th,m2t,mt,h);
if(m==1, return([t,1]));
if(m==2, return(ellffaddpol(e,t,t)));
m2=m>>1; \\ = m/2 rounded down
mth=ellffmultpol(e,t,m2);
m2t=mth[1];
h=mth[2];
mth=ellffaddpol(e,m2t,m2t);
mt=mth[1];
h=h*h*mth[2];
if(m%2==0, return([mt,h]));  \\ even m
mth=ellffaddpol(e,mt,t);    \\ odd m
h=h*mth[2];
mt=mth[1];
[mt,h]
}

\\ Redundant and inefficient old version:
\\ ellffmultpol(e,t,m)[1] is a much faster way to get the function
\\ (uses repeated doubling instead of loop!)

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellffweilpol(e,t,m) : provided m*t=0,
\\ returns a function h whose divisor is
\\     m(t)-m(0)
\\ No check is made that m*t=0 or that m>0
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellffweilpol(e,t,m)=
local(kt,h);
if(m==1,return(1));
if(m==2,return(ellffvert(e,t)));
if(m==3,return(ellfftang(e,t)));
kt=elladd(e,t,t);
h=ellfftang(e,t);
for(k=2,m-2,
	  h*=ellffchord(e,t,kt);
	  h/=ellffvert(e,kt);
          kt=elladd(e,kt,t)
   );
h
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ elldivpoint(e,q,m) : returns a list of points p such that m*p=q
\\
\\ q=[0] is allowed (hence one can find m-torsion).
\\
\\ all integers m are allowed.
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
elldivpoint(e,q,m)=
local(pt,mpt,xmpt,fx,fxj,xp,yps,both,ans);
pt=ellffgenpt(e);
\\print("genpt=",pt);
if(m==0,if(q==[0],return([[0]]),return([])));
if(m==1,return([q]));
if(m==-1,return([ellsub(e,[0],q)]));
mpt=ellpow(e,pt,m);
\\print("mpt=",mpt);
xmpt=polcoeff(lift(mpt[1]),0);
\\print("xmpt=",xmpt);
if(q==[0],fx=factor(denominator(xmpt)),
fx=factor(numerator(xmpt)-q[1]*denominator(xmpt)));
ans=[];if(q==[0],ans=[[0]]);
both=ellpow(e,q,2)==[0];
for(j=1,#(fx~),if(poldegree(fxj=fx[j,1])==1,
xp=-polcoeff(fxj,0)/polcoeff(fxj,1);yps=ellordinate(e,xp);
if(#yps,ans=concat(ans,[[xp,yps[1]]]);
if(both&&(#yps>1),ans=concat(ans,[[xp,yps[2]]])))));
ans;
}
