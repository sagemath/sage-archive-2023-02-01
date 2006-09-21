\r bgc.gp
\r lambda.gp
\r forms.gp

\\ The functions below require some global variables and functions:

global(H_E, H_ratp, H_qi, H_classno, H_classno_0, H_N, H_D, H_p, H_ap, H_an, H_z, H_bnd, H_rootbnd, H_indexbnd, H_sum, H_indexsum, H_indexq, H_pnum, H_ht, H_lambdas, H_ind);

{
printproj(p)=
local(z);
z=denominator(p[2]);
print("[",z*p[1],":",z*p[2],":",z,"]");
}

\\test if (z1+m*x)/k mod L gives a good point on e;  e[15] is the real period.

{
testz(e,z,v)=
local(p,ok,rp,d2,nd,d,hprec,ratp);
p=ellztopoint(e,z);
ok=0;
rp=[real(p[1]),real(p[2])];
nd=1;hprec=precision(1.0);
while((nd<hprec)&&(ok==0),d2=denominator(bestappr(rp[1],10^nd));
if(v>1,print("nd=",nd,", d2=",d2));
  if(issquare(d2),ok=1;d=sqrtint(d2);
if(v>1,print("d2=",d2,", d=",d));
  ratp=[round(d2*rp[1])/d2,round(d*d2*rp[2])/(d*d2)];
if(v>1,print("ratp=",ratp));
  ok=ellisoncurve(e,ratp));
  nd=nd+1);
H_ratp=[0];
if(ok,H_ratp=ratp;print();print("***RATIONAL P = ",ratp);
                  print1(" = ");printproj(ratp);
                  print("h(P) = ",ellheight(e,ratp));
                  if(v,print("x(P) = ",rp)));
H_ratp;
}

{test(e,z1,k,m,v)=local(z2,P);
  if(v,print1([k,m]," : "));
  z2=(z1+m*e[15])/k; P=testz(e,z2);
  if(P!=[0],return(P));
  if(e[12]<0,return(P));
  if(v,print1([k,m],"b: "));
  testz(e,z2+e[16]/2);
}

{testk(e,z1,k)=local(P);
for(m=0,k-1,P=test(e,z1,k,m);if(P!=[0],return(P)));
P;
}

{testrun(e,z1,kk)=local(P);
for(k=1,kk,P=testk(e,z1,k);if(P!=[0],return(P)));
P;
}

\\test if (z1+m*x)/k mod L gives a good point on e; given the expected
\\ height (ht) and a list lambdas=lambdalist(e), using Silverman's trick

{testzht(e,ht,lambdas,z,v)=local(p,xrp,ratp);
H_ratp=[0];
p=ellztopoint(e,z);
xrp=real(p[1]);
if(v>1,print("xrp=",xrp));
ratp=makept(e,lambdas,ht,xrp,v);
if(v>1,print("ratp=",ratp));
if(((length(ratp)>1) && ellisoncurve(e,ratp)),
     H_ratp=ratp;
     print();
     print("***RATIONAL P = ",ratp);
     print1(" = ");printproj(ratp);
     print("h(P) = ",ellheight(e,ratp));
     if(v>0,print("real x(P)     = ",xrp))
);
H_ratp;
}

{testht(e,ht,lambdas,z1,k,m,v)=local(z2,P);
z2=(z1+m*e[15])/k;
if(v,print1([k,m]," : "));
P=testzht(e,ht,lambdas,z2,v);
if(P!=[0],return(P));
if(e[12]<0,return(P));
if(v,print1([k,m],"b: "));
testzht(e,ht,lambdas,z2+e[16]/2);
}

{testkht(e,ht,lambdas,z1,k)=local(P);
for(m=0,k-1,P=testht(e,ht,lambdas,z1,k,m,1);if(P!=[0],return(P)));
return(P);
}

{testrunht(e,ht,lambdas,z1,kk)=local(P);
for(k=1,kk,P=testkht(e,ht,lambdas,z1,k);if(P!=[0],return(P)));
return(P);
}

\\The following assumes odd rank; superseded by elllderiv(e,1) in bg.gp
{ldash1(e,an,n)=local(c,l,ll);
l=length(an);
c=2*Pi/sqrt(n);
ll=round(precision(0.0)*log(10)/c);
if(l>ll,l=ll);
if(l<10,l=ll;an=ellan(e,ll));
print("ldash1 using ",l," terms");
2*sum(j=1,l,an[j]*eint1(c*j)/j,0.0);
}

\\ Superseded by ellanalyticrank(e)[3] in bg.gp;
{getht(e,an)=local(gr,l1,t,om,ht);
gr=ellglobalred(e);
l1=ldash1(e,an,gr[1]);
print("L'(E,1)=",l1);
t=elltors(e)[1];
om=e[15];if(e[12]>0,om=2*om,);
ht=l1*t*t/om/gr[3];
ht}


{ltwist1_old(e,an,n,d)=local(l,t);
l=length(an);t=exp(2*Pi/(d*sqrt(n)));
2*sqrt(abs(d))*sum(j=1,l,(an[j]*kronecker(d,j)/j)*t^j)
}

{twistcurve(e,d)=local(ed,gr);
ed=ellinit([0,0,0,-27*d*d*(e.c4),-54*d*d*d*(e.c6)]);
gr=ellglobalred(ed);
ellchangecurve(ed,gr[2])
}

{ltwist1(e,d)=local(ed,ll,localbnd,c,maxbnd); maxbnd=10^7;
ed=twistcurve(e,d);
c=2*Pi/sqrt(ellglobalred(ed)[1]);
localbnd=round((precision(0.0))*log(10)/c);
if(localbnd>maxbnd,print("bound in ltwist1()=",localbnd,", reducing to ",maxbnd);
                 localbnd=maxbnd);
ll=elllderiv(ed,0,localbnd);
if(abs(ll)<0.001,0,ll)
}

{heegindex_old(e,an,n,d,v)=local(gr,t,om,ind2);
gr=ellglobalred(e);
if(!isfundamental(d),print("Warning: ",d," is not a fundmental discriminant!");
                     return(0));
if(gcd(n,d)>1,print("heegindex requires gcd(n,d)=1!");
              return(0));
if(!isprime(abs(d)),print("Warning: non-prime discriminant, may not be correct"));
t=elltors(e)[1];
om=2*abs(imag(e[16]));
if(e[12]<0,om=2*om);
ind2=ltwist1_old(e,an,n,d)*gr[3]/(om*t*t);
if(d==-3,ind2=ind2*9,if(d==-4,ind2=ind2*4));
if(v,print("Square of index before rounding = ",ind2));
sqrtint(round(ind2))
}

{heegindex(e,d,v)=local(gr,n,t,om,ind2);
gr=ellglobalred(e);n=gr[1];
if(!isfundamental(d),print("Warning: ",d," is not a fundmental discriminant!");
                     return(0));
if(gcd(n,d)>1,print("heegindex requires gcd(n,d)=1!");
              return(0));
if(v&&!isprime(abs(d)),print("Warning: non-prime discriminant, may not be correct"));
t=elltors(e)[1];
om=2*abs(imag(e[16]));
if(e[12]<0,om=2*om);
ind2=sqrt(abs(d))*ltwist1(e,d)*gr[3]/(om*t*t);
if(d==-3,ind2=ind2*9,if(d==-4,ind2=ind2*4));
if(v,print("Square of index before rounding = ",ind2));
sqrtint(round(ind2))
}

\\ test how close z is to the lattice
{ex1(e,z)=local(w1,w2,a,b);
w1=e[15];w2=e[16];
\\We assume that w1 is real here
b=imag(z)/imag(w2);a=real(z-b*w2)/w1;
\\Now z=a*w1+b*w2 with a,b real
\\We return a,b and a measure of how close to L z is
[a,b,abs(round(a)*w1+round(b)*w2-z)]
}

\\ reduce z mod lattice
{ex(e,z)=local(w1,w2,a,z1);
w1=e[15];w2=e[16];
\\We assume that w1 is real here
a=imag(z)/imag(w2);z1=z-round(a)*w2;
a=real(z1)/w1;z1=z1-round(a)*w1;
z1;
}


\\ NB We are adding the real parts only here, so we may be out by
\\ half a period at the end;
\\ hence we must double up at the end (and double the index)

{
H_add(n,i,a,last_a) =
	local(j,j0,next_a,p,k,ap,aovern);
\\print("\tUsing n = ",n, "\tan = ",a);
        if(n<=H_rootbnd,H_an[n]=a
\\;print("\tStoring a[",n,"]=",a)
           );
	if(a==0,j0=i,
                aovern=(1.0*a)/n;
                for(k=1,H_classno_0,H_sum[k]+=(aovern*real(H_qi[k]^n)));
		k=kronecker(H_D,n);
                if(k!=0,H_indexsum+=aovern*k*H_indexq^n);
	        j0=1;);
	j=j0;
        if(j<=H_pnum,p=H_p[j],
                     print("Calling prime(",j,")");prime(j));
        if((a==0)&&(p>H_rootbnd),return);
	while(j<=i && ((p*n) <= H_bnd),
		ap=if(j<=H_pnum,H_ap[j],
                                if((H_E.c4==0)&&(p%3==2),0,
                                if((H_E.c6==0)&&(p%4==3),0,
                                ellap(H_E,p))));
		next_a=a*ap;
		if((j==i) && ((H_N%p)!=0),next_a = next_a - p*last_a);
		H_add(p*n,j,next_a,a);
		j=j+1;
                p=nextprime(p+1);
	)
}

{
H_addup(e,D,b,alist,v) =
  local(cached_pmax,alist_0,amax,i,ap,z,theprec,xprec,gr,ind,cp,t,om,indexfactor,H_XD,tau,p,mmax,qip,qin,H_indexqp,H_indexqn,n,aovern,k);
  H_E	   = ellinit(e);
  gr	   = ellglobalred(H_E);
  H_N	   = gr[1];
  H_XD     = 2*Pi/(sqrt(H_N)*D);
  H_D      = D;
  H_classno=length(alist);
 if(v,print("Class number = ",H_classno));
  cp=conjpairs(H_N,D,b,alist);
 if(v,print("Conjugate pairs = ",cp));
  H_classno_0=length(cp);
  alist_0=vector(H_classno_0,k,alist[cp[k][1]]);
  amax=vecmax(alist_0);

  t=elltors(e)[1];
  om=abs(imag(e[16]));  if(H_E.disc<0,om=2*om);
  indexfactor=if(D==-3,9,if(D==-4,4,1))*gr[3]*sqrt(abs(D))/(om*t*t);

  theprec=precision(1.0);
  H_bnd    = round(amax*theprec*log(10)*H_N/(Pi*sqrt(abs(D))));
  if(H_bnd>10^8,print("H_bnd=",H_bnd,", reducing to 10^8");
                H_bnd=10^8;
                xprec=H_bnd*Pi*sqrt(abs(D))/(amax*log(10)*H_N);
                print("Effective precision will be at most ",xprec," d.p."));
  H_rootbnd=floor(sqrt(H_bnd))+1;
  cached_pmax  = H_rootbnd;
  if(cached_pmax>default(primelimit,,1),default(primelimit,cached_pmax+100));
  H_pnum       = pith(cached_pmax);
  if(v,print1("Caching ",H_pnum," [p,ap] pairs, p< ",cached_pmax," ... "));
  H_p      = primes(H_pnum);
  H_ap     = vector(H_pnum,i,ellap(H_E,H_p[i]));
\\ Version 2.0.19 has acknowledged bug in ellap for p=2
  H_ap[1]  = ellak(H_E,2);
  if(v,print("done"));

  H_indexbnd = floor(-10*log(10)/H_XD);
  print("Summing up to ",H_bnd,", first stage up to ",H_p[#H_p],", index sum up to ",H_indexbnd);

\\ This will store cached an as they are computed
  H_an=vector(H_rootbnd,j,0);
  H_an[1]=1;

  tau      =(-b+sqrt(D))/H_N;
  H_qi=vector(H_classno_0,k,exp(Pi*I*tau/alist_0[k]));
  H_sum    = vector(H_classno_0,k,real(H_qi[k]));
  H_indexq = precision(exp(H_XD),9);  \\ low precision sufficient
  H_indexsum = H_indexq;
\\ print("Starting recursion...");
\\ Phase 1, up to the sqrt bound:
  for(i=1,H_pnum,
    p=H_p[i];
    ap=H_ap[i];
    if(v,if((i<100)||(i%100==0),print1("p=",p,"\tap=",ap)));
    H_add(p,i,ap,1);
    if(v,if((i<100)||(i%100==0),  z=0;
          for(k=1,H_classno_0,if(cp[k][1]==cp[k][2],z+=H_sum[k],z+=2*H_sum[k]));
          print("\tz=",z,"\tindex^2 ~ ",indexfactor*H_indexsum)));
  );
\\print("H_an=",H_an);
  if(v,print("------------"));
\\ Phase 2, the rest
  i=H_pnum;
  while(p<H_bnd,
    i+=1;
    p=nextprime(p+1);
    ap= if((H_E.c4==0)&&(p%3==2),0,
        if((H_E.c6==0)&&(p%4==3),0,
        ellap(H_E,p)));
    if(ap==0,next); \\ N.B. p^2 > H_bnd!
    mmax=floor(H_bnd/p);
    if(v,if((i<100)||(i%100==0),print1("p=",p,"\tap=",ap)));
    qip=vector(H_classno_0,k,H_qi[k]^p);
    qin=vector(H_classno_0,k,1.0);
    H_indexqp=H_indexq^p;
    H_indexqn=precision(1.0,9);
    for(m=1,mmax,  n=p*m;  aovern=ap*H_an[m]*1.0/n;
\\print("\tUsing n = ",n, "\tan = ",ap*H_an[m]);
        for(k=1,H_classno_0,qin[k]*=qip[k]);
        H_indexqn*=H_indexqp;
        if(H_an[m]==0,next);
        for(k=1,H_classno_0,
                            H_sum[k]+=(aovern*real(qin[k])));
        if(n<=H_indexbnd,
        	k=kronecker(H_D,n);
                if(k!=0,H_indexsum+=(aovern*k*H_indexqn)))
       );
    if(v,if((i<100)||(i%100==0),  z=0;
          for(k=1,H_classno_0,if(cp[k][1]==cp[k][2],z=z+H_sum[k],z=z+2*H_sum[k]));
          print1("\tz=",z);
          if(p<=H_indexbnd,print("\tindex^2 ~ ",indexfactor*H_indexsum),print())
  ))
  );
  z=0;
  for(k=1,H_classno_0,if(cp[k][1]==cp[k][2],z=z+H_sum[k],z=z+2*H_sum[k]));

\\ Now compute the index:
  if(v,print("H_indexsum = ",H_indexsum));
  ind=indexfactor*H_indexsum;
  if(v,print("Square of index = ",ind));
  ind=sqrtint(round(ind));
  if(v,print("Index = ",ind));

  if(ind==0, print("Twist has positive rank!"));

  [z,ind];
}

{H_getpoint(e,D,b,alist,ind,ht,lambdas,v) =
local(zind,z,P);
\\print("Calling H_addup with [e,D,b,alist,v] = ",[e,D,b,alist,v]);
zind=H_addup(e,D,b,alist,v);
\\print("H_addup returns [z,ind] = ",zind);
z=zind[1]; if(ind==0,ind=zind[2]);
H_ind=ind;
if(ind==0,print("Twist has positive rank, this discriminant is no good!"),
if(v,print("H_z=",z,"\nreducing modulo real period..."));
z=z-e[15]*round(z/e[15]);
H_z=z;
if(v,print("H_z=",z));
nt=elltors(e)[1];
if(v,print("Dividing 2*|T|*H_z by ",2*nt*ind,"..."));
P=testkht(e,ht,lambdas,2*nt*z,2*nt*ind));
print("---------------------------------------------------------");
P;
}

{H_point(e,dmin,dmax,v)=local(e0,gr,N,ht,lambdas,fn,np,plist,d,h,bai,P);
e0=vector(5,j,e[j]);print("E=",e0);
e=ellinit(e0);
gr=ellglobalred(e);
\\if(gr[2]!=[1,0,0,0],print("Non-minimal model input!");
if(gr[2][1]!=1,print("Non-minimal model input!");
              e=ellchangecurve(e,gr[2]);
              print("replacing with minimal model ",vector(5,j,e[j]));
              gr=ellglobalred(e));
N=gr[1];print("N=",N);
ht=ellanalyticrank(e)[3];print("Expected height=",ht);
lambdas=lambdalist(e);
H_ht=ht; H_lambdas=lambdas;
fn=mattranspose(factor(N));np=length(fn);
plist=vector(np,j,[fn[1,j],fn[2,j]]);
for(p=dmin,dmax,d=-p;
 if(isfundamental(d)&&(gcd(N,d)==1),
  if(good_disc(d,plist),h=qfbclassno(d);
    print("[D,h]=[",d,",",h,"]");
    bai=findb(N,d);
    if(v,print("[b,alist]=",bai));
    P=H_getpoint(e,d,bai[1],bai[2],0,ht,lambdas,v)
    )
   )
)
}

\\ Thanks to the global variables, one should be able
\\ to recompute the point using:
\\
\\ testkht(e,H_ht,H_lambdas,H_z,H_ind);

\\ In the next function we already know the index for the given
\\ discriminant and the height:

{H_point_ind(e,ht,d,ind,v)=local(e0,gr,N,lambdas,h,bai);
e0=vector(5,j,e[j]);print("E=",e0);
e=ellinit(e0);
gr=ellglobalred(e);N=gr[1];print("N=",N);
lambdas=lambdalist(e);
H_ht=ht; H_lambdas=lambdas;
h=qfbclassno(d);
H_ind=ind;
print("\n\n[D,h,index]=[",d,",",h,",",ind,"]");
bai=findb(N,d);
if(v,print("[b,alist]=",bai));
H_getpoint(e,d,bai[1],bai[2],ind,ht,lambdas,v)
}

\\ function to look for good discriminants

\\ We try all fundamental discriminants in the range dmin<=|D|<=dmax
\\ satisfying the conditions, such that the rank of the twist of e by
\\ D is 0, returning the best of the first H_MAXND found.  We do not store
\\ the quadratic form data for the returned D but it is trivial to
\\ recompute.

H_MAXND=5;

{H_disc(e,N,dmin,dmax,v)=local(fn,np,plist,d,h,bai,b,alist,cp,h0,alist_0,amax,best_d,ht,best_ht,nd,oldprec);
oldprec=default(realprecision,,1);
if(v,print("Saving old precision ",oldprec));
default(realprecision,19);  \\ for computing L(ed,1) approximately
if(v,print("N=",N));
fn=mattranspose(factor(N));np=length(fn);
plist=vector(np,j,[fn[1,j],fn[2,j]]);
best_d=0; best_ht=0;  nd=0;
for(p=dmin,dmax,d=-p;
 if(isfundamental(d)&&(gcd(N,d)==1),
  if(good_disc(d,plist),
    if(v,print1("\nTesting rank condition for D=",d," ... "));
    if(ltwist1(e,d)>0,
    if(v,print("OK"));
    nd+=1;
    h=qfbclassno(d);
    bai=findb(N,d); b=bai[1]; alist=bai[2];
    cp=conjpairs(N,d,b,alist);
    h0=length(cp);
    if(v,print("[D,h,h0]=[",d,",",h,",",h0,"]"));
    alist_0=vector(h0,k,alist[cp[k][1]]);
    if(v,print("[b,alist0]=",[b,alist_0]));
    amax=vecmax(alist_0);
    ht=sqrt(-d)/amax;
    if(ht>best_ht,best_ht=ht;best_d=d);
    if(v,print("max a = ",amax, ":\tsqrt(|d|)/maxa = ",ht));
    if(nd==H_MAXND,default(realprecision,oldprec);return(best_d)),
    if(v,print("no, twist has rank >=2"))
        )
      )
    )
   );
default(realprecision,oldprec);
best_d;
}

DISCLIM=1000;

{H_auto(e0,v) = local(e,N,d);
	   e=ellinit(e0);
	   N=ellglobalred(e)[1];
	   d=-H_disc(e,N,1,DISCLIM,v);
	   if(d==0,print("No suitable discriminant found with |d|<",DISCLIM),
	   print("Using discriminant -",d);
	   H_point(e,d,d,v));
}
