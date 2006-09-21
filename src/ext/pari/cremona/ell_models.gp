\\ Some functions for working with models of elliptic curves in GP \\

\\ For curves defined over Q these can test integrality locally or
\\ globally, and transform to an integral model by scaling (as in the
\\ static function ellintegralmodel() in elliptic.c); and test
\\ minimality locally or globally;  transforming to a minimal model is
\\ already provided in minimalmodel().

\\ In all the followign functions it is assumed that the elliptic
\\ curve e has rational coefficients (not necessarily integral) but we
\\ do not check this!

\\ In the functions which return a new model we would also like to
\\ return the transformation [u,0,0,0] used, as in the built-in function
\\ ellminimalmodel; but user-defined functions cannot have reference
\\ parameters :(  -- but for the future, we do compute the factor u such
\\ that output = ellchangecurve(e,[u,0,0,0])

\\
\\ INTEGRALITY FUNCTIONS
\\

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Checking integrality:  cheats by using type-checking
\\ e is an ellinit structure of any length (5,13 or 19)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellisintegral(e)=
prod(j=1,5,(type(e[j])=="t_INT"));
}

\\ Examples:
\\ ellisintegral([0,0,1,-7,6])
\\ returns 1 (true)
\\ ellisintegral(ellchangecurve(ellinit([0,0,1,-7,6]),[2,0,0,0]))
\\ returns 0 (false)
\\ ellisintegral(ellchangecurve(ellinit([0,0,1,-7,6]),[1/2,0,0,0]))
\\ returns 1 (true)

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Checking local integrality at a prime:
\\ e is an ellinit structure of any length (5,13 or 19),
\\ p should be a prime
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellislocallyintegral(e,p)=
prod(j=1,5,(valuation(e[j],p)>=0));
}

\\ Example:
\\ e = ellchangecurve(ellinit([0,0,1,-7,6]),[2,0,0,0]);
\\ (e.ai is [0, 0, 1/8, -7/16, 3/32])
\\ ellislocallyintegral(e,2)
\\ returns 0 (false)
\\ ellislocallyintegral(e,3)
\\ returns 1 (true)

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Making locally integral at a prime:
\\ e is an ellinit structure of any length (5,13 or 19),
\\ returns same length structure;
\\ p should be a prime
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
elllocalintegralmodel(e,p)=
local(ai,u);
u=p^vecmax(vector(5,k,-valuation(e[k],p)/if(k==5,6,k)));
ai=vector(5,k,e[k]*u^(if(k==5,6,k)));
if(#e==5,ai,ellinit(ai,#e==13));
}

\\ Example:
\\ e = ellchangecurve(ellinit([0,0,1,-7,6]),[6,0,0,0]).ai;
\\ returns [0, 0, 1/216, -7/1296, 1/7776]
\\ elllocalintegralmodel(e,2)
\\ returns [0, 0, 1/27, -7/81, 2/243]
\\ elllocalintegralmodel(e,3)
\\ returns [0, 0, 1/8, -7/16, 3/32]

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Making globally integral:
\\ e is an ellinit structure of any length (5,13 or 19)
\\ returns same length structure
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Method:  for each ai in turn and each p dividing denom(ai) we call
\\ the local version to scale up the short ai vector
\\
\\ To get the global u take the product of the local ones.

{
ellintegralmodel(e)=
local(ai,plist);
ai=vecextract(e,[1,2,3,4,5]);
for(i=1,5,plist=factor(denominator(ai[i]))[,1]~;
for(j=1,#plist,ai=elllocalintegralmodel(ai,plist[j])));
if(#e==5,ai,ellinit(ai,#e==13));
}

\\ Example:
\\ e = ellchangecurve(ellinit([0,0,1,-7,6]),[6,0,0,0]).ai;
\\ returns [0, 0, 1/216, -7/1296, 1/7776]
\\ ellintegralmodel(e)
\\ returns [0, 0, 1, -7, 6]

\\
\\ MINIMALITY FUNCTIONS
\\

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  test whether (integer) c4,c6 is a valid pair of invariants
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellvalidc4c6(c4,c6)=local(d=c4^3-c6^2);
if(d==0,return(0));
if(d%1728!=0,return(0));
if(valuation(c4,3)==2,return(0));
if(c6%4==3, return(1));
if(c4%16!=0, return(0));
d=c6%32;
(d==0)||(d==8)
}

\\Examples:
\\ ellvalidc4c6(336,-5400) returns 1
\\ ellvalidc4c6(random(),random()) probably returns 0

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  test local minimality at 2
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellislocallyminimalat2(e)=
ellislocallyintegral(e,2) &&
(valuation(e.c4,2)<4) || (valuation(e.c6,2)<6) ||
!ellvalidc4c6(e.c4>>4,e.c6>>6)
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  test local minimality at 3
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellislocallyminimalat3(e)=
ellislocallyintegral(e,3) &&
(valuation(e.c4,3)<4) || (valuation(e.c6,3)<6) ||
!ellvalidc4c6(e.c4/3^4,e.c6/3^6)
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  test local minimality at arbitrary p
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellislocallyminimal(e,p)=
ellislocallyintegral(e,p) &&
if(p==2,ellislocallyminimalat2(e),
 if(p==3,ellislocallyminimalat3(e),
  (valuation(e.c4,p)<4)||(valuation(e.c6,p)<6)))
}

\\ Examples:
\\
\\ e=ellinit([0,0,1,-7,6]);
\\ e1=ellchangecurve(e,[1/120,0,0,0,0]); (ellisintegral(e1)==1)
\\ ellislocallyminimal(e1,2) returns 0
\\ ellislocallyminimal(e1,3) returns 0
\\ ellislocallyminimal(e1,5) returns 0
\\ ellislocallyminimal(e1,7) returns 1

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  test global minimality
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellisminimal(e)=
local(badp);
if(!ellisintegral(e),return(0));
badp=factor(abs(e.disc))[,1]~;
for(i=1,#badp,if(!ellislocallyminimal(e,badp[i]),return(0)));
1
}

\\ Examples:
\\
\\ e=ellinit([0,0,1,-7,6]);
\\ e1=ellchangecurve(e,[1/120,0,0,0,0]); (ellisintegral(e1)==1)
\\ ellisminimal(e) returns 1
\\ ellisminimal(e1) returns 0

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  Construct e from c4, c6
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellinitfromc4c6(c4,c6)=
local(a1,a2,a3,a4,a6,b2,b4,b6);
if(!ellvalidc4c6(c4,c6),
print("Error in ellinitfromc4c6(): ",c4,",",c6," are not valid!");
return(0),
b2 = centerlift(Mod(-c6,12));
b4 = (b2^2-c4)/24;
b6 = (-b2^3+36*b2*b4-c6)/216;
a1 = b2%2;
a3 = b6%2;
a2 = (b2-a1)/4;
a4 = (b4-a1*a3)/2;
a6 = (b6-a3)/4;
ellinit([a1,a2,a3,a4,a6]));
}

\\ Examples:
\\
\\ e=ellinitfromc4c6(336, -5400);
\\ e.ai returns [0, 0, 1, -7, 6]

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  Given valid integral c4, c6 return those of a minimal model
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellminimisec4c6(c4,c6)=
local(u,d,g,plist,p,a,b);
u=1;
d = (c4^3-c6^2)/1728;
g=gcd(c4,c6); if(g==1,return([c4,c6]));
g=gcd(c6^2,d); if(g==1,return([c4,c6]));
plist = factor(g)[,1]~;
for(i=1,#plist,
    p = plist[i];
    d = floor(valuation(g,p)/12);
    if(p==2,
	a = (c4 >> (4*d)) % 16;
	b = (c6 >> (6*d)) % 32;
	if ((( (b%4)!=3) && !( (a==0) && (( b==0) || (b==8) ))),d-=1));
    if(p==3, if (valuation(c6,3)==(6*d + 2), d-=1));
    if(d>0, u *= p^d));
 [c4/u^4,c6/u^6];
}

\\ Examples:
\\
\\ e=ellinit([0,0,1,-7,6]);
\\ e1=ellchangecurve(e,[1/120,0,0,0]);
\\ ellminimisec4c6(e1.c4,e1.c6) returns [336, -5400]

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  Minimal model without Tate's algorithm:
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellminimalmodel2(e)=
local(c4c6);
e=ellintegralmodel(e);
c4c6=ellminimisec4c6(e.c4,e.c6);
ellinitfromc4c6(c4c6[1],c4c6[2]);
}



