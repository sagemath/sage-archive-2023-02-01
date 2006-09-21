\\ Buhler-Gross algorithm for computing L-series derivatives

\\ Implementation by tom@womack.net
\\ Revised by john.cremona@nottingham.ac.uk

\\ This is practical for any rank and conductors <= 10^10 in 10
\\ minutes per L' calculation (on P3/500E).

\\ Section 0: prime-number utility

pith(p)=primepi(p)

\\ Section 1: G_r

\\ We need a load of Magic Constants, which are the coefficients
\\ of the Taylor expansion of GAMMA(1+s). Now, [BGZ] we have
\\ log(Gamma(1+s)) = -gamma s + sum(n=2,infinity,(-1)^n*zeta(n)/n)
\\ and Pari can do polynomial exponential

\\ These need to be replaced with better functions from Cohen's new
\\ (2004) book...

global(ellld_ll,ellld_mm,ellld_pn_coeffs,ellld_s);
global(ellld_pnum,ellld_ap,ellld_p,ellld_E);
global(ellld_rootbnd,ellld_an,ellld_sum,ellld_r,ellld_X,ellld_bnd,ellld_N);

ellld_ll = -Euler*ellld_s + sum(n=2,30,(-1)^n * zeta(n) / n * ellld_s^n);
ellld_mm = exp(ellld_ll);
ellld_pn_coeffs=vector(17,i,polcoeff(ellld_mm,i-1));

ellld_P(r,t)=sum(m=0,r,ellld_pn_coeffs[r-m+1]*t^m/factorial(m));

\\ Warning: G2 does not work properly for large x
{G2(x)=
local(eps,ans,p,n,term);
  eps = 10.0^(1-precision(x));
  ans = -log(x) - Euler;
  ans=ans*ans/2 + Pi*Pi/12;
  n = p = 1;
  while(n<500,
      p /=n;  p*= -x;
      term = (p/n)/n;
      ans += term;
      if(abs(term/ans)<eps,break);
      n += 1);
   ans;
}


\\ Warning: G3 does not work properly for large x

{G3a(x)=
local(eps,ans,p,n,term);
\\  if (x>40, return(0));
\\  print("G3a(",x,")....");
  eps = 10.0^(-1-precision(x));
  ans = -log(x) - Euler();
  ans = (Pi()^2 + 2*ans^2) * ans/12 - zeta(3)/3;
  p = -1;
  n = 1;
  while(n<500,
      p *= -x/n;
      term = p/(n^3);
      ans += term;
      if(abs(term/ans)<eps,break);
      n += 1);
if(n==500,print("Reached end of loop in G3()!"));
  ans;
}

{G3(x)=if(x<30,G3a(x),0);
}

{ellld_G(r,x)=
   local(ss,oldss,n,sn,eps);
   eps = 10.0^(1-precision(x));
   if(r==0,return(2*exp(-x)));
   if(r==1,return(2*eint1(x)));
   if(r==2,return(2*G2(x)));
   if(r==3,return(2*G3(x)));
   if (x>30,return(0));
   ss = ellld_P(r,log(1/x));
   oldss = ss-1;n=1;sn=-((-1)^r);
   while (abs(ss-oldss)>eps,
     oldss = ss;
     ss = ss + sn*x^n/(n^r*factorial(n));
     n=n+1;sn=-sn;
  );
  2*ss;
}

\\ Section 2: the BG recursion
\\ Some global variables are necessary, but I've named them
\\ to try to avoid contaminating the namespace *too* much.

ellld_getap(r)=if(r<=ellld_pnum,ellld_ap[r],ellap(ellld_E,prime(r)));
ellld_getapp(p)=ellap(ellld_E,p);
ellld_getp(r)=if(r<=ellld_pnum,ellld_p[r],prime(r));

{
BGadd(n,i,a,last_a) =
	local(j,j0,next_a,p);
        if(n<=ellld_rootbnd,ellld_an[n]=a);
	if(a==0,j0=i,ellld_sum+=a*ellld_G(ellld_r,ellld_X*n)/n;j0=1;);
	j=j0;
        p=ellld_getp(j);
	if((a==0)&&(p>ellld_rootbnd),return);
	while(j<=i && p*n <= ellld_bnd,
		next_a=a*ellld_getap(j);
		if(j==i && ellld_N%p!=0,
			next_a = next_a - p*last_a);
		BGadd(p*n,j,next_a,a);
		j=j+1;
                p=nextprime(p+1);
	)
}

{
elllderiv(e,r,bnd,v) =
  local(cached_pmax,N,p,ap,i,n,aovern,mmax);
  if(v>1,print("in elllderiv with r = ",r,", bnd = ",bnd));
  ellld_E      = e;
  N	       = ellglobalred(e)[1];
if(v>1,print("N = ",N));
  ellld_N      = N;
  ellld_X      = 2*Pi/sqrt(N);
  ellld_r      = r;
  ellld_bnd    = bnd;
  ellld_rootbnd= sqrt(bnd);
if(v>2,print("ellld_rootbnd = ",ellld_rootbnd));
  ellld_sum    = ellld_G(r,ellld_X);
if(v>2,print("ellld_sum = ",ellld_sum));
  ellld_an     = vector(floor(ellld_rootbnd)+1,j,0);
  ellld_an[1]  = 1;
if(v>2,print("ellld_an initially = ",ellld_an));

  cached_pmax  = 10+floor(sqrt(ellld_bnd));
if(v>2,print("cached_pmax = ",cached_pmax,", calling pith..."));
  ellld_pnum   = pith(cached_pmax);
if(v>2,print("ellld_pnum = ",ellld_pnum));
  ellld_p      = primes(ellld_pnum);
if(v>2,print("ellld_p = ",ellld_p));
  ellld_ap     = vector(ellld_pnum,i,ellap(e,ellld_p[i]));
if(v>2,print("ellld_ap = ",ellld_ap));
\\ Version 2.0.19 has acknowledged bug in ellap for p=2
  ellld_ap[1]  = ellak(e,2);
if(v>2,print("ellld_ap = ",ellld_ap));

\\ Recursion for n up to rootbnd:
if(v>1,print("First stage, using recursion for p up to ",cached_pmax));
\\  for(i=1,pith(bnd),
  for(i=1,ellld_pnum,
    p=ellld_getp(i); ap=ellld_getap(i);
    if(v>1,if((i<100)||(i%100==0),print1("p=",p,"\tap=",ap)));
    BGadd(ellld_getp(i),i,ellld_getap(i),1);
    if(v>1,if((i<100)||(i%100==0),print("\tsum=",ellld_sum)))
     );
if(v>1,print("Second stage, looping for p up to ",bnd));
\\ Looping for large n, using cached values
  p=prime(ellld_pnum); i=ellld_pnum;
\\p=bnd+1;
  while(p<ellld_bnd,p=nextprime(p+1);i=i+1;
    ap=ellld_getapp(p);
    if(ap!=0,
    mmax=floor(ellld_bnd/p);
    mmax=min(mmax,ellld_rootbnd);
    if(v>1,if((i<100)||(i%100==0),print1("p=",p,"\tap=",ap)));
    for(m=1,mmax,  n=p*m;  aovern=ap*ellld_an[m]*1.0/n;
    if(v>2,print("\tUsing n = ",n, "\tan = ",ap*ellld_an[m]));
    ellld_sum+=aovern*ellld_G(ellld_r,ellld_X*n));
    if(v>1,if((i<100)||(i%100==0),print("\tsum=",ellld_sum))));
  );

  ellld_sum;
}


{
  ellanalyticrank(e,v)=
  local(rksgn,egrd,cond,localbnd,Lr1,om,ntors,prodcp,c,prelocalbnd,rk,regsha);
  egrd=ellglobalred(e);
\\ We had to compute this anyway to get the conductor, and
\\ ellap gives wrong answers on non-minimal curves ...
  if(egrd[2][1]!=1,
    print("INPUT CURVE NOT MINIMAL");
    e=ellchangecurve(e,egrd[2]);
    print("MINIMISED TO ",[e[1],e[2],e[3],e[4],e[5]]));
  cond=egrd[1]; prodcp=egrd[3];
  om=e[15]; if(e[12]>0,om=2*om);
  ntors=elltors(e)[1];
  if(v,print("Conductor = ",cond));
  rksgn=ellrootno(e);
  if(v,if(rksgn==-1,print("Rank is odd"),print("Rank is even")));
\\ Initial computation of L^(r)(1) to low precision...
  c=2*Pi/sqrt(cond);
  prelocalbnd=round(5*log(10)/c);
  if(v,print("Step 1 (low precision): Summing ",prelocalbnd," a_n terms"));
\\  print("prelocalbnd = ",prelocalbnd);
  rk=(1-rksgn)/2;
  Lr1=0;
\\ And I'm sure it's conjectural that this is close-enough to 0
  while(abs(Lr1)<0.0001,
    Lr1=elllderiv(e,rk,prelocalbnd,v);
    if(v,print("L^(",rk,")=",Lr1," (approx)"));
    rk+=2);
  rk-=2;
  if(v,print("Analytic rank = ",rk));
  return(rk);  \\ Lr1 and regsha are totally bogus for rank 2 curves, so don't bother

  regsha=(Lr1*ntors^2)/(om*prodcp);
  if(v,print("Approx R*S = ",regsha));
  if(v,print("\nNow recomputing L^(",rk,")(E,1) to higher precision"));
  localbnd=round((15+precision(0.0))*log(10)/c);
  if(localbnd>10^7,print("bound in ellanalyticrank()=",localbnd,", reducing to 10^7");
                   localbnd=10^7);
  if(v,print("Step 2 (high precision): Summing ",localbnd," a_n terms"));
  Lr1=elllderiv(e,rk,localbnd,v);
  if(v,print("L^(",rk,")=",Lr1));

\\ Also compute the quantity reg*sha from the BSD conjecture
  regsha=(Lr1*ntors^2)/(om*prodcp);
  [rk,Lr1,regsha];
}

{ellRegSha(e,r)=
  local(egrd,om,Lfn);
  om=e[15];if(e[12]>0,om=2*om);
  egrd=ellglobalred(e);
\\  Lfn=elllderiv(e,r,2*sqrt(egrd[1]));
  Lfn=elllderiv(e,r,precision(1.0)*log(10)*sqrt(egrd[1])/(2*Pi));
  (Lfn*elltors(e)[1]^2)/(om*egrd[3]);
}

ellSha(e,pts)=ellRegSha(e,length(pts))/matdet(ellheightmatrix(e,pts));
