/* Division polynomial functions for elliptic curves over arbitary fields

Author: John Cremona (john.cremona@nottingham.ac.uk)
Last edited: 13/10/04

*/

\\ elldivpol0(e,n) returns a polynomial in x whose zeros are the
\\ (x-coordinates of the) non-zero points P on e satisfying n*P=0
\\ and 2*P!=0

\\ elldivpol(e,n)  returns a polynomial in x,y whose zeros are the
\\ non-zero points P on e satisfying n*P=0.  Same as elldivpol0(e,n) when
\\ n is odd, but when n is even is equals elldivpol0(e,n)*(2*y+a1*x+a3).

\\ Users are only expected to use elldivpol normally

global(x,y);

{
elldivpol(e,n)=
local(ans=elldivpol0(e,n));
if(n%2==1,ans,ans*(2*y+e[1]*x+e[3]))
}

{
elldivpol0(e,n)=
local(m,a1,a2,a3,a4,a6,t1,t2,f1,f2,psi24);
a1=e[1];a2=e[2];a3=e[3];a4=e[4];a6=e[5];
f1=x^3+a2*x^2+a4*x+a6; f2=a1*x+a3;
n=abs(n);
if(n==0,return(0));
if(n==1,return(1));
if(n==2,return(1));
if(n==3,return(3*x^4+(a1^2+4*a2)*x^3+(3*a1*a3+6*a4)*x^2+(3*a3^2+12*a6)*x+a1^2*a6-a1*a3*a4+a2*a3^2+4*a2*a6-a4^2));
if(n==4,return(2*x^6+(a1^2+4*a2)*x^5+(5*a1*a3+10*a4)*x^4+(10*a3^2+40*a6)*x^3+(10*a1^2*a6-10*a1*a3*a4+10*a2*a3^2+40*a2*a6-10*a4^2)*x^2+(a1^4*a6-a1^3*a3*a4+a1^2*a2*a3^2+8*a1^2*a2*a6-a1^2*a4^2-4*a1*a2*a3*a4-a1*a3^3-4*a1*a3*a6+4*a2^2*a3^2+16*a2^2*a6-4*a2*a4^2-2*a3^2*a4-8*a4*a6)*x+a1^3*a3*a6-a1^2*a3^2*a4+2*a1^2*a4*a6+a1*a2*a3^3+4*a1*a2*a3*a6-3*a1*a3*a4^2+2*a2*a3^2*a4+8*a2*a4*a6-a3^4-8*a3^2*a6-2*a4^3-16*a6^2));
\\ general case, use recursion
\\ If n is odd, n=2m+1:
if(n%2==1,m=(n-1)/2;
t1=elldivpol0(e,m+2)*elldivpol0(e,m)^3;
t2=elldivpol0(e,m-1)*elldivpol0(e,m+1)^3;
psi24=(4*f1+f2^2)^2;
if(m%2==1,return(t1-psi24*t2),return(psi24*t1-t2)));
\\ Now n is even, n=2m:
m=n/2;
t1=elldivpol0(e,m+2)*elldivpol0(e,m-1)^2;
t2=elldivpol0(e,m-2)*elldivpol0(e,m+1)^2;
elldivpol0(e,m)*(t1-t2);
}

/* Examples:
  e=ellinit([0, -1, 1, 0, 0]);
  for(j=1,5,print(j,": ",elldivpol(e,j)))
displays
 1: 1
 2: (2*y + 1)
 3: 3*x^4 - 4*x^3 + 3*x - 1
 4: (4*y + 2)*x^6 + (-8*y - 4)*x^5 + (20*y + 10)*x^3 + (-20*y - 10)*x^2
     + (8*y + 4)*x + (-2*y - 1)
 5: 5*x^12 - 20*x^11 + 16*x^10 + 95*x^9 - 285*x^8 + 360*x^7 - 255*x^6 +
    94*x^5 + 15*x^4 - 45*x^3 + 25*x^2 - 5*x

 factor(elldivpol(e,5))
shows that e has a rational five-isogeny


NB Readers might expect elldivpol(e,2) to be 4*x^3 - 4*x^2 + 1
since the roots of this are indeed the x-coordinates of the 2-torsion
points.  But as a function on the curve, the latter has degree 6, with
double zeros at the 2-torsion points (and a pole of order 6 at
infinity), while 2*y+1 has degree 3 with a simple zero at each
2-torsion point.  Of course, (2*y+1)^2=4*x^3 - 4*x^2 + 1.

*/


\\ For convenience we therefore also provide the function elldivpol2():

\\ the 2-division x-polynomial of [a1,a2,a3,a4,a6]

elldivpol2(e)=4*(x^3+e.a2*x^2+e.a4*x+e.a6)+(e.a1*x+e.a3)^2;


{
elltorsionpoints(e,m)=
ellpointsfromx(e,polratroots(elldivpol(e,m)))
}

