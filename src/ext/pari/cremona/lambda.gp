
\\Returns lambda_bad list for one prime p
\\e must be a long elliptic curve structure, p a prime
{lambda1(e,p) = local(kod,m,n,nl);
   kod=elllocalred(e,p)[2];
   if(kod>4,
      m=kod-4;n=valuation(e[12],p);nl=1+floor(m/2);
      vector(nl,j,(j-1)^2/n-(j-1))*log(p),
   if(kod<-4,
      m=-kod-4;[0,-1,-1-m/4]*log(p),
   if((kod==2)|(kod==-2)|(kod==-1)|(p>3),
      [0],if(kod==3,[0,-1/2]*log(p),
   if(kod==4,
      [0,-2/3]*log(p),
   if(kod==-3,
      [0,-3/2]*log(p),
   if(kod==-4,
      [0,-4/3]*log(p),
   )))))))
}

{lambdalist(e)=local(n,fn,dis,np,plist,ans,nans,p,ll,nll,nnewans,newans,t);
   n=ellglobalred(e)[1];
   fn=mattranspose(factor(n)); np=length(fn);
   plist=vector(np,j,fn[1,j]);
\\print(plist);
   dis=e[12];ans=[0];nans=1;
   for(j=1,np,
       p=plist[j];
       if(dis%(p*p)!=0,,
          ll=lambda1(e,p);
          nll=length(ll);
\\print("p=",p,"; no. lambdas = ",nll);
          nnewans=nans*nll;newans=vector(nnewans,k,0);
          for(k=1,nans,
              t=ans[k];
              for(l=1,nll,
                  newans[k+nans*(l-1)]=t+ll[l]
                 )
              );
          nans=nans*nll;ans=newans)
       );
   ans}


{makept(e,lambdas,ht,x,v)= local(rh,nl,p,il,logd,approxd,d,d2,n,nn,ydmax,yd,ydiff,ylist,ddiff,ndiff,approxn);
   rh=realheight(e,x); if(v>1,print("rh=",rh));
   nl=length(lambdas);p=[0];il=1;
   while((il<=nl)&&(p==[0]),
         logd=(ht-rh-lambdas[il])/2;
         if(v>1,print("Trying lambda number ",il,", logd=",logd));
         il=il+1;
         approxd=exp(logd); if(v>1,print("approxd=",approxd));
         d=round(approxd); if(v>1,print("rounded d=",d));
         if((d>0)&&(abs(d-approxd)<0.01),
            d2=d*d;approxn=d2*x;n=round(approxn); if(v>1,print("rounded n=",n));
	    ddiff=abs(d-approxd); if(v>1,print("ddiff=",ddiff));
	    ndiff=abs(n-approxn); if(v>1,print("ndiff=",ndiff));
if(v>1,print("[approxd,d,n,x] = ",[approxd,d,n,x]));
            if(abs(d-approxd)<0.00001,
if(v>1,print("[il,d,ddiff,ndiff,x] = ",[il,d,ddiff,ndiff,x,abs(x-n/d2)]));
		   yd=0; ydmax=100000;
		   while((yd<ydmax)&&(p==[0]),
                    yd+=1;
                    ydiff=if(yd%2==0,yd/2,(1-yd)/2);
                    nn=n+ydiff;
                    ylist=ellordinate(e,nn/d2);
                    if(length(ylist)>0,p=[nn/d2,ylist[1]]);
                    )
                  )
              )
        );
   p
}

{realheight(e,x) = local(b2,b4,b6,b8,b2dash,b4dash,b6dash,b8dash,h,t,nlim,b,ans,f,w,z,zw,dmu,hprec);
b2=1.0*e[6];b4=1.0*e[7];b6=1.0*e[8];b8=1.0*e[9];
b2dash=b2-12;
b4dash=b4-b2+6;
b6dash=b6-2*b4+b2-4;
b8dash=b8-3*b6+3*b4-b2+3;
h=4;
t=abs(b2);if(t>h,h=t,);
t=2*abs(b4);if(t>h,h=t,);
t=2*abs(b6);if(t>h,h=t,);
t=abs(b8);if(t>h,h=t,);
hprec=precision(1.0);
nlim=ceil(5*hprec/3 + 1/2 + (3/4)*log(7+(4/3)*log(h)));
if(abs(x)<1/2,t=1/(x+1);b=0,t=1/x;b=1);
ans=-log(abs(t)); f=1;
for(n=0,nlim,
\\print([n,z,dmu,ans]);
f=f/4;
if(b==1,w = (((b6*t + 2*b4)*t + b2)*t + 4)*t;
        z = 1 - t*t*(b4 + t*(2*b6 + t*b8));
        zw = z + w,
        w = (((b6dash*t + 2*b4dash)*t + b2dash)*t + 4)*t;
        z = 1 - t*t*(b4dash + t*(2*b6dash + t*b8dash));
        zw = z - w);
if(abs(w)<=2*abs(z),dmu = f*log(abs(z)); ans = ans+dmu; t = w/z,
                    dmu = f*log(abs(zw));ans = ans+dmu; t = w/zw;b=1-b)
);
ans}

