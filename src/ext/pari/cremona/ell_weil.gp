/* Functions for computing Weil Pairing and related things
   -- mostly over an arbitrary ground field k, but at present we need
   to be able to generate random points on the curves which is only
   implemented for k=Z/pZ

Author: John Cremona (john.cremona@nottingham.ac.uk)
Last edited: 13/10/04

*/

\\ Two "dummy" functions, not needed as they shadow built-in functions
\\ elladd() and ellpow(), but included as templates for later versions to
\\ follow

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellsum(e,t,s) returns s+t
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellsum(e,t,s)=elladd(e,s,t);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellmult(e,t,m) : recursive function for m>0
\\ returns m*t
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellmult(e,t,m)=
local(m2,m2t,mt);
if(m==1, return(t));
if(m==2, return(ellpow(e,t,2)));
m2=m>>1; \\ = m/2 rounded down
m2t=ellmult(e,t,m2); \\ the recursive call
mt=elladd(e,m2t,m2t); \\ double
if(m%2==0, return(mt));  \\ even m: return
mt=elladd(e,mt,t);    \\ odd m: add
mt
}

\\ in ell_ff.gp there are functions ellffsumpol(e,s,t) and
\\ ellffmultpol(e,s,t) with return [s+t,h] resp. [m*t,h] where
\\ h is in the function field and satisfies
\\
\\ div(h) = (s)+(t)-(s+t)-(0), resp.
\\ div(h) = m*(t)-(m*t)-(m-1)*(0).
\\
\\ In practice it is more efficient to do the following:  instead of
\\ constructing h and then evaluating it later, we evaluate as we go
\\ along.  For this to work we must only evaluate at a point (or divisor)
\\ whose support is disjoint from all points in the support of the
\\ functions arising in the contruction of h, which will (for the mult
\\ function) be some subset of the multiples of t of size of the order of
\\ log(m).
\\
\\ To achive this we first give version of the above which return as
\\ second component a Set of the x-coordinates of the points to be
\\ avoided, lift-ed in case the are of Mod() type.
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
{
ellsum_supp(e,t,s)=
local(st);st=elladd(e,s,t);
return([st,Set(lift([t[1],s[1],st[1]]))]);
}
\\
{
ellmult_supp(e,t,m)=
local(m2,mth,m2t,mt,supp);
\\print("In ellmult_supp with m=",m);
if(m==1, return([t,Set([lift(t[1])])]));
if(m==2, return(ellsum_supp(e,t,t)));
m2=m>>1; \\ = m/2 rounded down
mth=ellmult_supp(e,t,m2);
m2t=mth[1];
supp=mth[2];
mth=ellsum_supp(e,m2t,m2t);
mt=mth[1];
supp=setunion(supp,mth[2]);
if(m%2==0, return([mt,supp]));  \\ even m
mth=ellsum_supp(e,mt,t);    \\ odd m
supp=setunion(supp,mth[2]);
mt=mth[1];
[mt,supp]
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Now the evaluation versions.  We evaluate at a divisor of the form
\\ (s1)-(s2) : the functions in question are only defined up to a
\\ constant but this is well-defined
\\
\\ NB It is the caller's responsibility to ensure that none of the
\\ evaluation points are in the support of the partial functions
\\ evaluated!  See the function ellweilpairing() for an example
\\
\\ We use the function ellffaddpol() from ell_ff.gp
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellsumpoleval2(e,t1,t2,s1,s2)=
local(mth);
mth=ellffaddpol(e,t1,t2);
[mth[1],ellffeval2(mth[2],s1,s2)]
}

{
ellmultpoleval2(e,t,m,s1,s2)=
local(m2,mth,mt,h);
if(m==1, return([t,1]));
if(m==2, return(ellsumpoleval2(e,t,t,s1,s2)));
m2=m>>1; \\ = m/2 rounded down
mth=ellmultpoleval2(e,t,m2,s1,s2);
mt=mth[1];
h=mth[2]^2;
mth=ellsumpoleval2(e,mt,mt,s1,s2);
mt=mth[1];
h*=mth[2];
if(m%2==0, return([mt,h]));          \\ even m
mth=ellsumpoleval2(e,mt,t,s1,s2);    \\ odd m
mt=mth[1];
h*=mth[2];
[mt,h]
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ ellffweilpairing(e,s,t,m) : provided m*t=m*s=0,
\\ returns e_m(s,t) = Weil Pairing of s,t, an m'th root of unity
\\
\\ NB This uses ellzprandompoint() so only works over Z/pZ.  It could
\\ work elsewhere if we had a way of finding  suitable auxiliary point r
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

debug_weil=0;

{
ellweilpairing2(e,s,t)=local(p=e.disc.mod);
if(debug_weil,print("Weil pairing for m = 2"));
if(s==t,return(Mod(1,p)),return(Mod(-1,p)));
}

{
ellweilpairing3(e,s,t)=local(p=e.disc.mod,t2,s2,ts,tt,w);
t2=elladd(e,t,t);
if(debug_weil,print("Weil pairing for m = 3"));
if((s==t)||(s==t2),return(Mod(1,p)));
ts=ellfftang(e,s);
tt=ellfftang(e,t);
s2=elladd(e,s,s);
w=(ellffeval(ts,t)*ellffeval(tt,s2))/(ellffeval(ts,t2)*ellffeval(tt,s));
w=w^2;
if(debug_weil,print("returning = ",w));
return(w);
}

{
ellweilpairing4(e,s,t)=local(p=e.disc.mod,t2,s2,t3,ts,tt,w);
if(debug_weil,print("Weil pairing for m = 4"));
s2=elladd(e,s,s);
t2=elladd(e,t,t);
if(s2==[0],return(ellweilpairing2(e,s,t2)));
if(t2==[0],return(ellweilpairing2(e,s2,t)));
if(s2==t2,return(ellweilpairing2(e,s2,ellsub(e,s,t))));
t3=elladd(e,t2,t);
if((s==t)||(s==t3),return(Mod(1,p)));
ws1=ellffweilpol(e,s,4);
ws2=ellffweilpol(e,s2,4);
wt1=ellffweilpol(e,t,4);
wt2=ellffweilpol(e,t2,4);
w=(ellffeval(ws2,t2)*ellffeval(ws1,t)*ellffeval(wt2,s)*ellffeval(wt1,s2))/
  (ellffeval(wt2,s2)*ellffeval(wt1,s)*ellffeval(ws2,t)*ellffeval(ws1,t2));
if(debug_weil,print("returning = ",w));
return(w);
}

{
ellweilpairing(e,s,t,m)=
local(badxt,badxs,t2,s2,ts,tt,r,rs,tr,w);
if(debug_weil,print("In ellweilpairing with m=",m,", s=",s,", t=",t));
if(ellpow(e,s,m)!=[0],print("Error in Weil Pairing: s must be in E[m]");
                      return(0));
if(ellpow(e,t,m)!=[0],print("Error in Weil Pairing: t must be in E[m]");
                      return(0));
if(s==[0],return(Mod(1,e.disc.mod)));
if(m==2,return(ellweilpairing2(e,s,t)));
if(m==3,return(ellweilpairing3(e,s,t)));
if(m==4,return(ellweilpairing4(e,s,t)));
if(debug_weil,print("m > 3"));
badxt=ellmult_supp(e,t,m)[2];
badxs=ellmult_supp(e,s,m)[2];
if(debug_weil,print("badxt=",badxt));
if(debug_weil,print("badxs=",badxs));
\\
r=ellzprandompoint(e);
rs=elladd(e,r,s); \\ rs=r+s
tr=ellsub(e,t,r); \\ tr=t-r
if(debug_weil,print("r=",r));
if(debug_weil,print("r+s=",rs));
if(debug_weil,print("t-r=",tr));
while( (rs==[0]) || (rs==t) || (tr==[0]) || (r==[0]) ||
      setsearch(badxt,Str(lift(r[1]))) ||
      setsearch(badxt,Str(lift(rs[1]))) ||
      setsearch(badxs,Str(lift(tr[1]))) ||
      setsearch(badxs,Str(lift(r[1]))) ,
if(debug_weil,print("No good, trying another random point r"));
r=ellzprandompoint(e);
rs=elladd(e,r,s);
tr=ellsub(e,t,r);
if(debug_weil,print("r=",r));
if(debug_weil,print("r+s=",rs));
if(debug_weil,print("t-r=",tr));
);
w = ellmultpoleval2(e,t,m,rs,r)[2] /
    ellmultpoleval2(e,s,m,tr,ellsub(e,[0],r))[2];
if(m%znorder(w)!=0,print("Error in Weil Pairing (order ",m,"):
returned value has order ",znorder(w)));
w
}

