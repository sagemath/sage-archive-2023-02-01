/* Some simple utility functions for elliptic curves in GP

Author: John Cremona (john.cremona@nottingham.ac.uk)
Last edited: 13/10/04

*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Negation of a point:
\\ e is an ellinit, p a point on e;  returns -p
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

ellneg(e,p)=ellsub(e,[0],p)

/* Example:
           e=ellinit([0,0,1,-7,6]); p=[-2,3];
           ellneg(e,p)
   returns [-2, -4]
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Member function to extract [a1,a2,a3,a4,a6] from a longer vector:
\\ e is an ellinit;  returns [a1,a2,a3,a4,a6]
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

e.ai=vecextract(e,[1,2,3,4,5])

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Inverse of a standard transformation [u,r,s,t]:
\\ urst=[u,r,s,t] encoding a standard transformation of Weierstrass
\\ equations (so u!=0);  returns the inverse trabsformation
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellinvert(urst)=
local(u=urst[1],r=urst[2],s=urst[3],t=urst[4]);
[1/u,-r/u^2,-s/u,(r*s-t)/u^3]
}

/* Example:
 e1=ellinit([1,2,3,4,6]);
 v=[u,r,s,t] \\generic
 e2=ellchangecurve(e1,v);
 e1==ellchangecurve(e2,ellinvert(v)) \\ =1 (true)
 P1=[-1,1]; ellisoncurve(e1,P1)       \\ =1 (true)
 P2=ellchangepoint(P1,v)
 P1==ellchangepoint(P2,ellinvert(v))   \\ =1 (true)
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Composition of two standard transformations:
\\ urst1, urts2 are two standard transformations; returns their composite
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellcompose(urst1,urst2)=
local(u1=urst1[1],r1=urst1[2],s1=urst1[3],t1=urst1[4],
      u2=urst2[1],r2=urst2[2],s2=urst2[3],t2=urst2[4]);
[u1*u2,u1^2*r2+r1,u1*s2+s1,u1^3*t2+s1*u1^2*r2+t1]
}

/* Example:
  ellcompose([u,r,s,t],ellinvert([u,r,s,t])) == [1,0,0,0]

  e1=ellinit([1,2,3,4,6]);
  v1=[u1,r1,s1,t1] \\generic 1
  v2=[u2,r2,s2,t2] \\generic 2
  e2=ellchangecurve(e1,v1);
  e3=ellchangecurve(e2,v2);
  e3==ellchangecurve(e1,ellcompose(v1,v2))
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\ Roots of a polynomial!  (not elliptic curve specific)
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
{
polratroots(pol)=
fx=factor(pol); ans=[];
for(j=1,#(fx~),if(poldegree(fxj=fx[j,1])==1,
ans=concat(ans,[-polcoeff(fxj,0)/polcoeff(fxj,1)])));
ans;
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Test of isomorphism between to curves e1,e2:  either returns urst
\\ such that ellchangecurve(e1,urst)==e2, or 0
\\ Requires field of definition such that issquare(), sqrt() and
\\ (for the case j=0 only) factorization of x^3-a to get a cube root
\\
\\ Do not expect this to work in characteristics 2,3!!
\\
\\ We don't use the form issquare(u2,&u) since it does not work (for
\\ example) over Z/pZ...
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
elliso(e1,e2)=
local(u6,u3,uu,u,u4,u2,r,s,t);

\\ Check we have at least 13-vectors:
if(#e1<13,e1=ellinit(e1,1));
if(#e2<13,e2=ellinit(e2,1));

\\ Check j-invariants equal:
if(!(e1.j==e2.j),return(0));

\\ Special case where j==c4==0:
if(e1.j==0,u6=e1.c6/e2.c6;
           if(!issquare(u6),return(0));
           u3=polratroots(x^2-u6)[1];
           uu=polratroots(x^3-u3);
           if(#uu==0,return(0),u=uu[1]),

\\ Special case where j-1728==c6==0:
if(e1.j==1728,u4=e1.c4/e2.c4;
              if(!issquare(u4),return(0));
              u2=polratroots(x^2-u4)[1];
              if(!issquare(u2),return(0));
              u=polratroots(x^2-u2)[1],

\\ Generic case, we have u^2:
                u2=(e1.c6*e2.c4)/(e2.c6*e1.c4);
                if(!issquare(u2),return(0));
                u=polratroots(x^2-u2)[1]));

\\Given u, solve for r,s,t:
 s=(u*e2.a1-e1.a1)/2;
 r=(u^2*e2.a2-e1.a2+s*e1.a1+s^2)/3;
 t=(u^3*e2.a3-e1.a3-r*e1.a1)/2;
 [u,r,s,t];
}

compsq(e) = if((e.a1==0)&&(e.a3==0),e,ellchangecurve(e,[1/2,0,-e.a1/2,-e.a3/2]));

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Linear combinations of points
\\ e is an ellinit, pts a vector of points, v a vector of integers
\\ with #v=#pts;  returns sum(v[j]*pts[j]) on E
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ (May become obsolete when factorback() is adapted)
{
ellcomb(E,pts,v)=
local(P); P=[0];
if(#pts==#v,
for(i=1,#pts, if(v[i],  P=elladd(E,P,ellpow(E,pts[i],v[i])))),
print("Error in ellcomb:  arguments 2 and 3 must have the same length!"));
P
}

/* Example:
  e=ellinit([0,0,1,-7,6]);
  pts=[[-2,3],[-1,3],[0,2]];
  ellcomb(e,pts,[-1,2,3])
  returns [13961/100, -1649791/1000]
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ The primes of bad reduction for an integral elliptic curve over Q
\\ e is an ellinit; returns a vector of primes dividing the
\\ discriminant (so these are bad reduction for this model but may not be
\\ bad for the underlying elliptic curve if this is not a minimal model)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

ellbadprimes(e)=factor(abs(e.disc))[,1]~

/* Example:
  ellbadprimes(ellinit([1,0,1,1,2]))
returns [2,3,5]
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\ Construct a vector of all points on e from a given vector of
\\ (possible) x-coordinates);  at most one from each x unless flag!=0
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
{
ellpointsfromx(e,xi,flag=0)=
local(yij,ans);
yij=ellordinate(e,xi);maxny=if(flag,2,1);
ans=vector(#xi,i,vector(min(#yij[i],maxny),j,[xi[i],yij[i][j]]));
if(#ans,concat(ans),[])
}

/* Example:
  e=ellinit([0,0,1,-7,6]);
  ellpointsfromx(e,vector(1004,j,j-4))
returns 18 integral points on e with x-coordinates between -3 and +1000
(these are all the integral points on e, up to sign, but that's another story!
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\ 2- and 3-torsion points (see ell_ff.gp for general case)
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ell2torsion(e)=
ellpointsfromx(e,polratroots(4*(x^3+e.a2*x^2+e.a4*x+e.a6)+(e.a1*x+e.a3)^2));
}

{
ell3torsion(e)=
ellpointsfromx(e,polratroots(3*x^4+(e.a1^2+4*e.a2)*x^3+(3*e.a1*e.a3+6*e.a4)*x^2+(3*e.a3^2+12*e.a6)*x+e.a1^2*e.a6-e.a1*e.a3*e.a4+e.a2*e.a3^2+4*e.a2*e.a6-e.a4^2));
}

/* Examples:
 e=ellinit([0,1,0,4,4]);
 ell3torsion(e)
 ell2torsion(e)
returns
 [[0, 2], [0, -2]]
 [[-1, 0]]

BUT
 e=ellinit(1.0*[0,0,0,-2,0]);
 ell2torsion(e)
returns 4 points, since ellordinate(e,sqrt(2)) returns two values,
both approximately zero!

 e=ellinit(Mod(1,47)*[0,0,0,-2,0]);
 lift(ell2torsion(ellinit(Mod(1,47)*[0,0,0,-2,0])))
 lift(ell3torsion(ellinit(Mod(1,47)*[0,0,0,-2,0])))
returns
[[0, 0], [40, 0], [7, 0]]
[[21, 30], [21, 17]]
(we "lift"ed the results for easier readability)
*/
