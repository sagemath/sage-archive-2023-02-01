\\Use this to find a suitable discriminant for conductor N:

\\vecin(V,e) returns 1 if e is in the vector, and 0 otherwise
{vecin(V,e)=
	for (i=1,length(V),if(V[i]==e,return(1)));
	return(0);
}

\\vecadduniq(V,e) adds e to the vector iff e wasn't there already
{vecadduniq(V,e)=
	local(i,lv,absent);
	i=1;lv=length(V);absent=1;
	while(i<=lv&&absent,if(V[i]==e,absent=0);i=i+1;);
	if(absent,concat(V,[e]),V);
}

\\good_disc(D,pl) checks whether discriminant D is good wrt a list of
\\prime factors in pl (listed as [prime,exp] pairs)

{good_disc(D,pl) =
	local(dd,np,OK,ip,p);
	dd = quaddisc(D);
	np = length(pl);
	OK=1;ip=1;
	while (OK&(ip<=np),
		p=pl[ip][1];
		if(p==2,OK=(dd%8)==1,
		if(pl[ip][2]>1,OK=kronecker(dd,p)==1,OK=kronecker(dd,p)!=-1));
		ip=ip+1;
	);
	OK
}

\\ This tells you which pairs of Heegner points will be (complex)
\\ conjugate, so you only need to compute (approximately) half of them.
\\ Within each pair you use the one with smallest a, for better accuracy.

{conjpairs(N,D,b,alist)=
local(ans,h,flist,f0);
ans=[];
h=length(alist);
flist=vector(h,j,Qfb(alist[j]*N,b,(b*b-D)/(4*N*alist[j])));
\\assumes alist[1]=1
f0=qfbred(flist[1]);
for(j=1,h,for(k=j,h,if(qfbred(qfbcompraw(flist[j],flist[k]))==f0,ans=concat(ans,[[j,k]]))));
ans
}

\\listforms(n,d,b) returns a list of a_i such that [n*a_i,b,c] are a
\\set of quadratic forms with discriminant d containing at most one
\\element from each class in Cl[Q(sqrt d)]

{listforms(n,d,b)=
	local(FL,ac,c,f,al);
	FL=[];al=[];
	ac=(b^2-d)/n/4;
	fordiv(ac,a,
		c=ac/a;
		f=qfbred(Qfb(a*n,b,c));
		if(vecin(FL,f)==0,
			FL=concat(FL,f);
			al=concat(al,a);
		)
	);
	al;
}

\\ findb1(n,d) solves b^2=d (mod 4n) for b (if possible)

{findb1(n,d,v)=
	local(b,bp,plist,fn,np,e,p);
	fn=mattranspose(factor(4*n));
	np=length(fn);
	plist=vector(np,j,[fn[1,j],fn[2,j]]);
	if(v,print("plist = ",plist));
	e=plist[1][2];
	if(e==2,b=Mod(d,2);  if(v,print("(e=2) b2 = ",b)),
          if(e==3,b=Mod(1,4); if(v,print("(e=3) b2 = ",b)),
        	b=sqrt(d+O(2^e));
	        if(v,print("(e>3) b2 = ",b));
	        b=Mod(b%(2^(e-1)),2^(e-1))));
        	if(v,print("b2 = ",b));
	        for(j=2,np,p=plist[j][1];e=plist[j][2];
	           bp=Mod(sqrt(d+O(p^e))%(p^e),p^e);
                   if(v,print("p=",p,"; bp = ",bp));
		   b=chinese(b,bp);
                   if(v,print("so far, b = ",b))
            );
	b=lift(b);
	if((b*b-d)%(4*n)==0,
          if(v,print("Success: b=",b)),
          print("Failure in findb1()"));
        b
}

\\ findb(n,d) prints out (if(v)) a list of entries of the form [b, [a_1, ..., a_h]]
\\ where [a_1, ..., a_h] are such that Qfb(a_i*n,b,c), where c is chosen to
\\ make the discriminant =d, cover the classes in Q(\sqrt d);
\\
\\ The returned set is the  one which minimises max(a_i).
\\ This could be improved, since the measure of "best" is just that
\\ the maximum a is least, not taking into account the conjugation pairing
\\ which means that around half the ai are not used.

{findb(n,d,v)=
	local(b,bb,L,h,maxa,maxal,bestb,bestL);
	h=qfbclassno(quaddisc(d));maxa=-1;
	b=findb1(n,d);if(v,print(b));
        for(t=-50,50,bb=b+t*2*n;
		     L=listforms(n,d,bb);
                     if (length(L)==h,maxal=vecmax(L);
                     if((maxa==-1)||(maxa>maxal),bestb=bb;bestL=L;maxa=maxal);
                     if(v,print([bb,L])))
        	);
[bestb,bestL];
}

