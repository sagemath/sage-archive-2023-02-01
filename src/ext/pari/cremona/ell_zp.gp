/* Functions for elliptic curves over the finite field Z/pZ;  could be
   extended to more general finite fields one day.

Author: John Cremona (john.cremona@nottingham.ac.uk)
Last edited: 13/10/04

*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ A random point on an elliptic curve defined over Z/pZ
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellzprandompoint(ep)=
local(p,xp,yp);
p=ep.a1.mod; yp=[];
while(#yp==0,yp=ellordinate(ep,xp=Mod(random(p),p)));
[xp,yp[1]];
}

/* Example:
 e=ellinit(Mod(1,101)*[0,0,1,-7,6]);
 lift(vector(5,j,ellzprandompoint(e)))
returns something like
 [[53, 6], [44, 6], [63, 54], [21, 5], [99, 3]]
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Reduction of rational points mod p
\\ [0] -> [0]
\\ [x,y] -> [0] if val(x),val(y)<0 else [Mod(x,p),Mod(y,p)]
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellreducepoint(P,p)=if(P==[0],[0],if(valuation(P[1],p)<0,[0],Mod(1,p)*P));
}

/* Examples:
ellreducepoint([0],5)
\\ returns [0]
ellreducepoint([1/25,1/125],5)
\\ returns  [0]
ellreducepoint([1/25,1/125],3)
\\ returns  [Mod(1, 3), Mod(2, 3)]
*/



\\ Extension to built-in data type: the return vector from ellinit()
\\ on input of type t_INTMOD has length 19 but entries 14-19 are unused
\\ and set to 0.  We use these as follows (provisional):
\\
\\ E[14] : the group order n = #E(Fp)
\\ E[15] : factorization of the group order n
\\ E[16] : the group structure: [n] if cyclic, or
\\                              [n1,n2] with n=n1*n2 and n2|n1
\\ E[17] : the group generators: [P] is cyclic, or
\\                               [P1,P2] with order(Pi)=ni
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Input can be either a vector of INT_MODs (in which case the
\\ parameter p is ignored so can be omitted) or a vector of INTs, in
\\ which case we need to also give the prime p.  No check is made that p
\\ is prime.  Singular data will cause a run-time error
\\
\\ Output is an extended ellinit structure as explained above
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellzpinit(ai,p)=
local(ep,gp);
ai=vecextract(ai,[1,2,3,4,5]);
if(type(ai[1])=="t_INTMOD",
  ep=ellinit(ai);
  p=ai[1].mod;
  ap=ellap(ellinit(lift(ep)),p),
if(p==0,
  print("Error in ellzpinit: p=0");
  return(0),
ep=ellinit(Mod(1,p)*ai);
ap=ellap(ellinit(ai,1),p)));
ep[15]=factor(ep[14]=1+p-ap);
gp=ellzpgroup(ep);
if(#gp<2,print("Problem finding group structure"),
ep[16]=gp[1];
ep[17]=gp[2]);
ep;
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Extracting the group order / factored group order:
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

ellzp.grouporder=ellzp[14];
ellzp.factoredgrouporder=ellzp[15];

/* Example:
 e=ellzpinit([0,0,1,-7,6],101);
 e.grouporder
 e.factoredgrouporder
returns
 98
[2 1]
[7 2]
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Extracting the group structure:
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

ellzp.groupstr=[ellzp[16],ellzp[17]];
ellzp.iscyclic=(#ellzp[16])==1;
ellzp.isotype=ellzp[16];
ellzp.generators=ellzp[17];

/* Example 1 :
 e=ellzpinit([0,0,1,-7,6],101);
 e.groupstr
 e.iscyclic
 e.isotype
 e.generators
returns (possibly with a different generator)
[[98], [[Mod(38, 101), Mod(8, 101)]]]
1
[98]
[[Mod(38, 101), Mod(8, 101)]]

/* Example 2 :
 e=ellzpinit([0,0,1,-7,6],103);
 e.groupstr
 e.iscyclic
 e.isotype
 e.generators
returns (possibly with a different first generator)
 [[52, 2], [[Mod(88, 103), Mod(25, 103)], [Mod(83, 103), Mod(51, 103)]]]
 0
 [52, 2]
 [[Mod(88, 103), Mod(25, 103)], [Mod(83, 103), Mod(51, 103)]]
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Order of a point (using the factored group order):
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
ellzppointorder(ellzp,t)=
local(ans,fm,p,m,s);
\\print("Finding order of ",t);
\\print("Group order =  ",ellzp.grouporder);
ans=1; fm=(ellzp.factoredgrouporder)~;
for(j=1,#fm,p=fm[1,j];
\\print("p=",p);
m=factorback(vecextract(fm,concat("^",Str(j)))~);
\\print("m=",m);
s=ellpow(ellzp,t,m);
\\print(s);
while(s!=[0],s=ellpow(ellzp,s,p);ans*=p));
ans;
}

/* Example:
 ellzp=ellzpinit([0,0,1,-7,6],101);
 vector(10,j,ellzppointorder(ellzp,ellzprandompoint(ellzp)))
returns something like
 [7, 49, 98, 49, 98, 98, 14, 98, 98, 49]
which shows (assuming that 98 does appear in the list, as it does with
probability more than 99.6%) that the group is cyclic.  We can then
get a generator:
 P=[0];while(ellzppointorder(ellzp,P)<98,P=ellzprandompoint(ellzp));P
returns (perhaps)
 Mod(3, 101), Mod(3, 101)]
*/

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Group structure:
\\ returns either
\\     [[n],[P]] (if cyclic of order n) where P is a generator
\\ or
\\     [[n1,n2],[P1,P2]] is (Z/n1)x(Z/n2) with n2|n1
\\                          and P1, P2 generators of order n1,n2
\\
\\ Internal function used in ellzpinit, need/should not be called by user
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

debug_group=0;

\\ Utility function: given positive integers m,n returns [l,m',n']
\\ where gcd(m',n')=1, lcm(m',n')=m'n'=lcm(m,n), m'|m, n'|n:
\\
\\ (Nothing to do with elliptic curves, may be useful elsewhere)

{
tidy_lcm(m,n)=
  local(g,l,m0,n0);m0=m;n0=n;
  g=gcd(m,n);
  l=m*n/g;  \\ = lcm(m,n)
  g=gcd(m,n/g);  \\ divisible by primes dividing n to a higher power than m
  while(g!=1, m/=g; g=gcd(m,g));
  n=l/m;
  if((l==m*n)&&
     (gcd(m,n)==1)&&
     (m0%m==0)&&
     (n0%n==0),,
  print("error in tidy_lcm(",m0,",",n0,")"));
  [l,m,n];
}

\\ Utility function: given a positive integer n returns the
\\ largest d such that d^2|n.
\\
\\ (Nothing to do with elliptic curves, may be useful elsewhere)

{
maxsquarefreediv(n)=
  local(fn);
  fn=factor(n);
  fn[,2]\=2;
  factorback(fn)
}

\\ Utility function: given P1 of order n1 and P2 of order n2,
\\ construct a point Q of order n3=lcm(n1,n2); returns [Q,n3]
\\
\\ (algorithm may be useful elsewhere)

{
ellzpmerge1(ellzp,P1,n1,P2)=
 local(n2,abc,m,m1,m2);
 if(ellpow(ellzp,P2,n1)==[0],return([P1,n1]));
 n2=ellzppointorder(ellzp,P2);
 if(n2%n1==0,return([P2,n2]));
 abc=tidy_lcm(n1,n2);
 m=abc[1];
 m1=n1/abc[2];
 m2=n2/abc[3];
 [elladd(ellzp,ellpow(ellzp,P1,m1),ellpow(ellzp,P2,m2)), m];
}

\\ Utility function: given independent points P1,P2 of order n1,n2,
\\ and a further point Q, replace P2 by a point still independent of
\\ P1 whose order is the lcm of n2 and the order of Q mod <P1>.
\\ n2target is n/n1 (n=group order) so n2|n2target and we are aiming
\\ to achieve n2=n2target
\\
\\ All of P1,P2,n1,n2 may be changed;  as we have no writeable
\\ variable parameters in gp we return [P1,n1,P2,n2]
\\
\\ (algorithm may be useful elsewhere)

{
ellzpmerge2(ellzp,P1,n1,P2,n2,n2target,Q)=
local(Q1,Q2,P1n1,a,w,m,l,ordQ,lmn);
if(debug_group>2,print("merge2 called with ",[P1,n1,P2,n2,Q,n2target]));
if(n2==0,print("Error: n2=0 at place A");1/0);
Q1 = ellpow(ellzp,Q,n2);
if(Q1==[0], if(debug_group>1,print("Order(Q) divides n2"));
            return([P1,n1,P2,n2])
  );
Q2 = ellpow(ellzp,Q1,n1/n2);
if(Q2!=[0],
      \\ P1 needs updating and we discard P2
      if(debug_group>1, print("Order(Q) does not divide n1, updating P1"));
      P1n1=ellzpmerge1(ellzp,P1,n1,Q);
      P1=P1n1[1]; n1=P1n1[2];
      if(debug_group>1, print("New P1 has order ",n1));
      if(n2>1, P2 = [0]; n2=1);
      return([P1,n1,P2,n2])
  );

\\ Now we find a multiple a*P1 such that Q-a*P1 is killed by
\\ n2target so we can apply the Weil Pairing of order n2target

Q1 = ellpow(ellzp,Q,n2target);
Q2 = ellpow(ellzp,P1,n2target);  \\ has exact order n1/n2target
a = bg_algorithm(ellzp,Q2,Q1,0,n1/n2target);
  if(debug_group,
      print("Dlog of ",Q1," w.r.t. ",Q2," (order ",n1/n2target,") is ",a);
      print("Check: a*Q2-Q1 =  ",ellsub(ellzp,ellpow(ellzp,Q2,a),Q1))
                                  \\a*Q2-Q1
    );
  Q = elladd(ellzp,Q,ellpow(ellzp,P1,-a)); \\ = Q-a*P1
  if(Q==[0],
      if(debug_group, print("Q-a*P1 = ",Q,", no use"));
      return([P1,n1,P2,n2])
    );
  if(debug_group,
      print("Replacing Q by Q-a*P1 = ",Q," where a = ", a);
      if(ellpow(ellzp,Q,n2target)==[0],
	print1("order divides n2target, OK; computing Weil pairing of ");
	print("order " ,n2target),
        print("order does not divide n2target -- WRONG!")
        ));
  \\ compute Weil pairing of P1 and Q, which is an n2target'th root of unity:
  Q1 = ellpow(ellzp,P1,(n1/n2target));
  if(debug_group,
      print("order((n1/n2target)*P1) = ",Q1," is ",ellzppointorder(ellzp,Q1));
      print("order(Q) =                ",Q," is ",ellzppointorder(ellzp,Q))
    );
  if(debug_group,print("Computing Weil pairing of those two points..."));
  w = ellweilpairing(ellzp,Q1,Q,n2target);
  m = znorder(w);
  if(debug_group, print("w = ", w," of order ", m ));
  \\ Compare this with n2 to see if we have gained:
  l = lcm(n2,m);
  if(l==n2, return([P1,n1,P2,n2])    ); \\ no gain
  ordQ = ellzppointorder(ellzp,Q);
  if(debug_group, print("ordQ = ", ordQ, "; ordQ/m = ",ordQ/m));
  Q1 = ellpow(ellzp,Q,ordQ/m);   \\ of order m
  if(l==m, \\ replace P2, n2 by Q1, m
      P2 = Q1;
      n2 = m;
      return([P1,n1,P2,n2])
    );
  \\ Now P2,Q1 have orders n2,m & both are independent of P1,
  \\ so we combine them to get a point of order l=lcm(n2,m)
  \\ still independent of P1
  lmn = tidy_lcm(n2,m);
  n2=lmn[2]; m=lmn[3]; if(n2==0,print("Error: n2=0 at place B");1/0);
  P2 = elladd(ellzp,ellpow(ellzp,P2,l/n2),ellpow(ellzp,Q,l/m));
  n2 = l; if(n2==0,print("Error: n2=0 at place C");1/0);
  if(debug_group,
      print("Changed P2 = ",P2,":\t order(P2) = ",n2);
      if(ellzppointorder(ellzp,P2)!=n2, print("that's wrong!"))
    );
  return([P1,n1,P2,n2]);
}

{
ellzpgroup(ellzp)=
 local(N,p,P1,n1,n2,n,Q,Pm,n2target);
 p=ellzp.a1.mod;
if(debug_group,
   print("Curve is  ", lift(vecextract(ellzp,[1,2,3,4,5])), " mod ",p);
   print("Factorization of p-1 is ",factor(p-1)));
 fN=ellzp.factoredgrouporder;      \\ factorization of the order
 N=ellzp.grouporder; \\ the order
if(debug_group,print("N = ", N, " is the group order"));
if(N==1,return([[],[]]));
 fN[,2]\=2;
 n2max = gcd(factorback(fN),p-1);
 iscyclic=(n2max==1);  \\ 0 means we don't yet know
if(debug_group,
 print("A priori maximal possible n2 = ", n2max);
 if(iscyclic,print("Group must be cyclic"))
);

\\ Start finding random points...

 P1 = ellzprandompoint(ellzp);
 n1 = ellzppointorder(ellzp,P1);
if(debug_group,print("P1 = ", P1, ", order ",n1));
 if(n1==N, return([[N],[P1]])); \\ we hit the jackpot, quick exit
 n2max = gcd(n2max,N/n1);
if(debug_group, print("maximal possible n2 = ", n2max));
 iscyclic=(n2max==1);  \\ 0 means we don't yet know
if(debug_group, if(iscyclic,print("Group must be cyclic")));

\\ use more random points to increase the order; after a few we are
\\ highly likely to be holding a point of maximal order (hence a
\\ generator when E is cyclic, which happens often)

 n=0;
 while((n<10)||((n1*gcd(n1,p-1))%N!=0),n+=1;
 Q = ellzprandompoint(ellzp);
 Pm=ellzpmerge1(ellzp,P1,n1,Q);
 P1=Pm[1]; n1=Pm[2];
 if(debug_group,print("P1 = ", P1, ", order ",n1));
 if(n1==N, return([[N],[P1]])); \\ cyclic, generator P1, exit
 n2max = gcd(n2max,N/n1);
 if(debug_group, print("maximal possible n2 = ", n2max));
 iscyclic=(n2max==1);  \\ 0 means we don't yet know
 if(debug_group, if(iscyclic,print("Group must be cyclic")));
);

 if(n1==N, return([[N],[P1]])); \\ cyclic, generator P1, exit

 \\ See whether [n1,N/n1] is possible group structure...

 if(N%n1,print("problem 1");return([])); \\contradicts Lagrange...
 n2=N/n1;
 if(n1%n2,print("problem 2");return([])); \\ should not happen
                              \\ since once we get here, (n1*gcd(n1,p-1))%N==0
 if((p-1)%n2,print("problem 3");return([])); \\ should not happen
 if(debug_group,
    print1("ellzpgroup found probable group structure, now finding second ");
    print("generator of order ",n2)
  );

  P2=[0];
  n2target=n2;
  n2=1; \\ holds order(P2) always
  while(n2<n2target,
    Q = ellzprandompoint(ellzp);
    if(debug_group, print("Using Q = ",Q));
    P1n1P2n2=ellzpmerge2(ellzp,P1,n1,P2,n2,n2target,Q);
    P1=P1n1P2n2[1]; newn1=P1n1P2n2[2];
    if(newn1>n1,
    n2target=n2target*n1/newn1;
    if(debug_group,
    print("Increased n1 to ",newn1,"; decreased n2target to ",n2target));
    n1=newn1);
    P2=P1n1P2n2[3]; n2=P1n1P2n2[4];
   if(n2==0,print("Error: n2=0 at place D");1/0);
    if(debug_group, print("Group order now ",(n1*n2)))
  );
  return( [[n1,n2],[P1,P2]]);
}

/* OLD CODE for second phase
 Q1=ellpow(ellzp,P1,n1/n2);
\\ At the moment this is done in a slightly stupid way....
 n2pts=elldivpoint(ellzp,[0],n2);
if(debug_group,
    print(#n2pts," points of order dividing ",n2," (should be ",n2^2,"): ");
    print(n2pts));
\\ Look for one of those whose Weil pairing with Q1 has order n2:
for(j=1,#n2pts,P2=n2pts[j];if(ellzppointorder(ellzp,P2)!=n2,next);print(P2);
               if((znorder(ellweilpairing(ellzp,Q1,P2,n2))==n2),
                  return( [[n1,n2],[P1,P2]])));
*/
