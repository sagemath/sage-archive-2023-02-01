\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\       Copyright (C) 2008 Denis Simon
\\
\\ Distributed under the terms of the GNU General Public License (GPL)
\\
\\    This code is distributed in the hope that it will be useful,
\\    but WITHOUT ANY WARRANTY; without even the implied warranty of
\\    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
\\    General Public License for more details.
\\
\\ The full text of the GPL is available at:
\\
\\                 http://www.gnu.org/licenses/
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\
\\ Auteur:
\\ Denis SIMON -> simon@math.unicaen.fr
\\ adresse du fichier:
\\ www.math.unicaen.fr/~simon/ellQ.gp
\\
\\  *********************************************
\\  *          VERSION 29/04/2008               *
\\  *********************************************
\\
\\ Programme de calcul du rang des courbes elliptiques sur Q.
\\ langage: GP
\\ pour l'utiliser, lancer gp, puis taper
\\ \r ellQ.gp
\\
\\ Ce programme utilise le module de resolution des formes quadratiques
\\ situe a l'adresse
\\ www.math.unicaen.fr/~simon/qfsolve.gp
\\ Il faut donc aussi taper:
\\ \r qfsolve.gp
\\
\\ Explications succintes :
\\ La fonction ellrank() accepte toutes les courbes sous la forme
\\ [a1,a2,a3,a4,a6]
\\ Les coefficients peuvent etre entiers ou non.
\\ L'algorithme utilise est celui de la 2-descente.
\\ La 2-torsion peut etre quelconque.
\\ Il suffit de taper :
\\
\\ gp > ell = [a1,a2,a3,a4,a6];
\\ gp > ellrank(ell)
\\
\\ Retourne un vecteur [r,s,vec]
\\ ou r est le rang probable (c'est toujours une minoration du rang),
\\ s est le 2-rang du groupe de Selmer,
\\ vec est une liste de points independants dans E(Q)/2E(Q).
\\
\\ Courbes de la forme: k*y^2 = x^3+A*x^2+B*x+C
\\ sans 2-torsion, A,B,C entiers.
\\ gp > bnf = bnfinit(x^3+A*x^2+B*x+C);
\\ gp > ell = ellinit([0,A,0,B,C],1);
\\ gp > rank = ell2descent_gen(ell,bnf,k);
\\
\\ Courbes avec #E[2](Q) >= 2:
\\ ell doit etre sous la forme
\\ y^2 = x^3 + A*^2 + B*x
\\ avec A et B entiers.
\\ gp > ell = [0,A,0,B,0]
\\ gp > ell2descent_viaisog(ell)
\\ = algorithme de la 2-descente par isogenies
\\ Attention A et B doivent etre entiers
\\
\\ Courbes avec #E[2](Q) = 4: y^2 = (x-e1)*(x-e2)*(x-e3)
\\ gp > ell2descent_complete(e1,e2,e3)
\\ = algorithme de la 2-descente complete
\\ Attention: les ei doivent etre entiers.
\\
\\
\\ On peut avoir plus ou moins de details de calculs avec
\\ DEBUGLEVEL_ell = 0;
\\ DEBUGLEVEL_ell = 1; 2; 3;...
\\
\\
\\

{
\\
\\ Variables globales usuelles
\\

  DEBUGLEVEL_ell = 1; \\ pour avoir plus ou moins de details
  LIM1 = 5;       \\ limite des points triviaux sur les quartiques
  LIM3 = 50;      \\ limite des points sur les quartiques ELS
  LIMTRIV = 50;   \\ limite des points triviaux sur la courbe elliptique

\\
\\  Variables globales techniques
\\

  MAXPROB = 20;
  LIMBIGPRIME = 30; \\ pour distinguer un petit nombre premier d'un grand
                    \\ utilise un test probabiliste pour les grands
                    \\ si LIMBIGPRIME = 0, n'utilise aucun test probabiliste
  ELLREDGENFLAG = 1;\\ pour reduire les genereteurs a la fin de l'algorithme
}

\\
\\  Programmes
\\

\\
\\ Fonctions communes ell.gp et ellQ.gp
\\
{
ellinverturst(urst) =
local(u = urst[1], r = urst[2], s = urst[3], t = urst[4]);
  [1/u,-r/u^2,-s/u,(r*s-t)/u^3];
}
{
ellchangecurveinverse(ell,v) = ellchangecurve(ell,ellinverturst(v));
}
{
ellchangepointinverse(pt,v) = ellchangepoint(pt,ellinverturst(v));
}
{
ellcomposeurst(urst1,urst2) =
local(u1 = urst1[1], r1 = urst1[2], s1 = urst1[3], t1 = urst1[4],
      u2 = urst2[1], r2 = urst2[2], s2 = urst2[3], t2 = urst2[4]);
  [u1*u2,u1^2*r2+r1,u1*s2+s1,u1^3*t2+s1*u1^2*r2+t1];
}
{
ellinverturst(urst) =
local(u = urst[1], r = urst[2], s = urst[3], t = urst[4]);
  [1/u,-r/u^2,-s/u,(r*s-t)/u^3];
}
{
ellchangecurveinverse(ell,v) = ellchangecurve(ell,ellinverturst(v));
}
{
ellchangepointinverse(pt,v) = ellchangepoint(pt,ellinverturst(v));
}
{
ellcomposeurst(urst1,urst2) =
local(u1 = urst1[1], r1 = urst1[2], s1 = urst1[3], t1 = urst1[4],
      u2 = urst2[1], r2 = urst2[2], s2 = urst2[3], t2 = urst2[4]);
  [u1*u2,u1^2*r2+r1,u1*s2+s1,u1^3*t2+s1*u1^2*r2+t1];
}
{
polratroots(pol) =
local(f,ans);
  f = factor(pol)[,1];
  ans = [];
  for( j = 1, #f,
    if( poldegree(f[j]) == 1,
      ans = concat(ans,[-polcoeff(f[j],0)/polcoeff(f[j],1)])));
  return(ans);
}
if( DEBUGLEVEL_ell >= 4, print("mysubst"));
{
mysubst(polsu,subsx) =
  if( type(lift(polsu)) == "t_POL",
    return(simplify(subst(lift(polsu),variable(lift(polsu)),subsx)))
  , return(simplify(lift(polsu))));
}
if( DEBUGLEVEL_ell >= 4, print("degre"));
{
degre(idegre) =
local(ideg,jdeg);

  ideg = idegre; jdeg = 0;
  while( ideg >>= 1, jdeg++);
  return(jdeg);
}
if( DEBUGLEVEL_ell >= 4, print("nfissquare"));
{
nfissquare(nf, a) = #nfsqrt(nf,a) > 0;
}
if( DEBUGLEVEL_ell >= 4, print("nfsqrt"));
{
nfsqrt( nf, a) =
\\ si a est un carre, renvoie [sqrt(a)], sinon [].
local(alift,ta,res,pfact,r1,py);

if( DEBUGLEVEL_ell >= 5, print("entree dans nfsqrt ",a));
  if( a==0 || a==1,
if( DEBUGLEVEL_ell >= 5, print("fin de nfsqrt"));
    return([a]));

  alift = lift(a);
  ta = type(a);
  if( !poldegree(alift), alift = polcoeff(alift,0));

  if( type(alift) != "t_POL",
    if( issquare(alift),
if( DEBUGLEVEL_ell >= 5, print("fin de nfsqrt"));
      return([sqrtrat(alift)])));

  if( poldegree(nf.pol) <= 1,
if( DEBUGLEVEL_ell >= 5, print("fin de nfsqrt"));
    return([]));
  if( ta == "t_POL", a = Mod(a,nf.pol));

\\ tous les plgements reels doivent etre >0
\\
  r1 = nf.sign[1];
  for( i = 1, r1,
    py = mysubst(alift,nf.roots[i]);
    if( sign(py) < 0,
if( DEBUGLEVEL_ell >= 5, print("fin de nfsqrt"));
      return([])));
\\ factorisation sur K du polynome X^2-a :
  if( variable(nf.pol) == x,
    py = subst(nf.pol,x,y);
    pfact = lift(factornf(x^2-mysubst(alift,Mod(y,py)),py)[1,1])
  ,
    pfact = lift(factornf(x^2-a,nf.pol)[1,1]));
  if( poldegree(pfact) == 2,
if( DEBUGLEVEL_ell >= 5, print("fin de nfsqrt"));
    return([]));
if( DEBUGLEVEL_ell >= 5, print("fin de nfsqrt"));
  return([subst(polcoeff(pfact,0),y,Mod(variable(nf.pol),nf.pol))]);
}
if( DEBUGLEVEL_ell >= 4, print("sqrtrat"));
{
sqrtrat(a) =
  sqrtint(numerator(a))/sqrtint(denominator(a));
}
\\
\\ Fonctions propres a ellQ.gp
\\

if( DEBUGLEVEL_ell >= 4, print("ellhalf"));
{
ellhalf(ell,P)=
\\ renvoie tous les points Q sur ell tels que 2Q = P.
local(pol2,ratroots,half,x2,y2,P2);

  if(#ell < 13, ell=ellinit(ell,1));

  pol2 = Pol([4,ell.b2,2*ell.b4,ell.b6]); \\ polynome de 2-division

  if( P == [0],
    ratroots = polratroots(pol2);
    half = vector(#ratroots,i,[ratroots[i],ellordinate(ell,ratroots[i])[1]]);
    half = concat( [[0]], half);
    return(half)
  );

  x2=Pol([1,0,-ell.b4,-2*ell.b6,-ell.b8]); \\ x(2P) = x2/pol2

  half = [];
  ratroots = polratroots(x2-P[1]*pol2);
  if( #ratroots == 0, return(half));
  for( i = 1, #ratroots,
    y2 = ellordinate(ell,ratroots[i]);
    for( j = 1, #y2,
      P2 = [ratroots[i],y2[j]];
      if( ellpow(ell,P2,2) == P, half = concat(half,[P2]))
    )
  );
  return(half);
}
if( DEBUGLEVEL_ell >= 4, print("elltors2"));
{
elltors2(ell)=
\\ Calcule le sous-groupe de 2-torsion de la courbe elliptique ell.
local(pol2,ratroots,tors2);

if( DEBUGLEVEL_ell >= 4, print("calcul de la 2-torsion"));
  if(#ell < 13, ell=ellinit(ell,1));
  tors2 = ellhalf(ell,[0]);
  if( #tors2 == 1,
    tors2 = [1, [], []],
  if( #tors2 == 2,
    tors2 = [2, [2], [tors2[2]]]
  , tors2 = [4, [2,2], [tors2[2],tors2[3]]]
  ));
if( DEBUGLEVEL_ell >= 4, print("E[2] = ",tors2));
  return(tors2);
}
if( DEBUGLEVEL_ell >= 4, print("elltorseven"));
{
elltorseven(ell)=
\\ Calcule le 2-Sylow sous-groupe de torsion de la courbe elliptique ell.
local(torseven,P2);

if( DEBUGLEVEL_ell >= 4, print("calcul de la 2^n-torsion"));
  if(#ell < 13, ell=ellinit(ell,1));
  torseven = elltors2(ell);

  while( torseven[1] != 1,
    P2 = ellhalf(ell,torseven[3][1]);
    if( #P2 > 0,
       torseven[1] *= 2;
       torseven[2][1] *= 2;
       torseven[3][1] = P2[1];
       next
    );
    if( #torseven[3] == 1, break());

    P2 = ellhalf(ell,torseven[3][2]);
    if( #P2 > 0,
       torseven[1] *= 2;
       torseven[2][2] *= 2;
       torseven[3][2] = P2[1];
       next
    );
    P2 = ellhalf(ell,elladd(ell,torseven[3][1],torseven[3][2]));
    if( #P2 > 0,
       torseven[1] *= 2;
       torseven[2][1] *= 2;
       torseven[3][1] = P2[1];
       next
    );
    break()
  );

if( DEBUGLEVEL_ell >= 4, print("E[2^n] = ",torseven));
  return(torseven);
}
if( DEBUGLEVEL_ell >= 4, print("polratroots"));
{
polratroots(pol) =
local(f,ans);
  f = factor(pol)[,1];
  ans=[];
  for( j = 1, #f,
    if( poldegree(f[j]) == 1,
      ans = concat(ans,[-polcoeff(f[j],0)/polcoeff(f[j],1)])));
  return(ans);
}
if( DEBUGLEVEL_ell >= 4, print("ratpoint"));
{
ratpoint(pol,lim=1,singlepoint=1,tryhard=0) =
\\ Recherche de points sur y^2=pol(x).
\\ Les coeff de pol sont entiers.
\\ Si singlepoint >= 1, cherche un seul point, sinon plusieurs.
\\ Si tryhard == 1, on essaye une autre strategie quand pol est imprimitif.

local(listpoints,point1,odd,deg4,pol16,tab16,pol9,tab9,pol5,tab5,pol0,vecz,vecx,lead,zz,xx,x16,x9,x5,evpol,ix,iz,K,factK,e,cont,ind,p,sol,factpol,r,M,U);

if( DEBUGLEVEL_ell >= 4, print("entree dans ratpoint avec pol = ",pol); print("lim = ",lim););
  if( !singlepoint, listpoints = []);
  point1 = [];

\\          cas triviaux
  if( issquare(pollead(pol)),
    point1 = [ 1, sqrtrat(pollead(pol)), 0];
if( DEBUGLEVEL_ell >= 3, print("solution triviale: a est un carre"));
    if( singlepoint,
if( DEBUGLEVEL_ell >= 4, print("fin de ratpoint"));
      return(point1));
    listpoints = concat(listpoints,[point1]));
  if( issquare(polcoeff(pol,0)),
    point1 = [ 0, sqrtrat(polcoeff(pol,0)) ];
if( DEBUGLEVEL_ell >= 3, print("solution triviale: e est un carre"));
    if( singlepoint,
if( DEBUGLEVEL_ell >= 4, print("fin de ratpoint"));
      return(point1));
    listpoints = concat(listpoints,[point1]));
  odd = poldegree(pol)%2;
  deg4 = poldegree(pol) == 4;

\\ initialisation du crible modulo 16, 9 et 5
  if( deg4,
    pol16 = (Vec(pol)*Mod(1,16))~;
    tab16 = matrix(16,16);
    for(xx = 0, 16-1,
      for(zz = 0, 16-1,
        tab16[xx+1,zz+1] = !issquare([xx^4,xx^3*zz,xx^2*zz^2,xx*zz^3,zz^4]*pol16)));
    pol9 = (Vec(pol)~)*Mod(1,9);
    tab9 = matrix(9,9);
    for(xx = 0, 9-1,
      for(zz = 0, 9-1,
        tab9[xx+1,zz+1] = !issquare([xx^4,xx^3*zz,xx^2*zz^2,xx*zz^3,zz^4]*pol9)));
    pol5 = (Vec(pol)~)*Mod(1,5);
    tab5 = matrix(5,5);
    for(xx = 0, 5-1,
      for(zz = 0, 5-1,
        tab5[xx+1,zz+1] = !issquare([xx^4,xx^3*zz,xx^2*zz^2,xx*zz^3,zz^4]*pol5)))
  );

  lead = pollead(pol);
  pol0 = polcoeff(pol,0);

  if( odd,
    vecz = vector(lim,i,i^2);
  ,
\\ si le degre de pol est pair, il faut que le coeff dominant soit
\\ un carre mod zz.
    vecz = vector(lim);
    zz = 0;
    for( i = 1, lim,
      zz++; while( !issquare(Mod(lead,zz)),zz++); vecz[i] = zz
  ));
\\ le coeff constant doit etre un carre mod xx.
  vecx = vector(lim);
  xx = 0;
  for( i = 1, lim,
    xx++; while( !issquare(Mod(pol0,xx)),xx++); vecx[i] = xx);

if( DEBUGLEVEL_ell >= 4, print("xmax = ",vecx[lim]));
if( DEBUGLEVEL_ell >= 4, print("zmax = ",vecz[lim]));

if( DEBUGLEVEL_ell >= 5, print("vecx = ",vecx));
if( DEBUGLEVEL_ell >= 5, print("vecz = ",vecz));

\\ boucle sur x = xx/zz
  for( somme = 2, 2*lim,
    for( ix = max(1,somme-lim), min(lim,somme-1),
      xx = vecx[ix]; iz = somme-ix; zz = vecz[iz];
      if( gcd(zz,xx) > 1, next);
      if( odd && !issquare(lead*Mod(xx,zz)), next);
      for( eps = 1, 2, if( eps == 2, zz = -zz);
      if( deg4 &&
        (tab16[xx%16+1,zz%16+1] || tab9[xx%9+1,zz%9+1] || tab5[xx%5+1,zz%5+1])
      , next);
      evpol = subst(pol,variable(pol),xx/zz);
      if( issquare(evpol),
        point1 = [xx/zz,sqrtrat(evpol)];
        if( singlepoint, break(3));
        listpoints = concat(listpoints,[point1])
  ))));

  if( point1 != [],
if( DEBUGLEVEL_ell >= 3, print("point trouve par ratpoint = ",point1));
if( DEBUGLEVEL_ell >= 4, print("sortie de ratpoint "));
    if( singlepoint, return(point1), return(listpoints))
  );

\\
\\ Essaye une autre strategie quand pol a un content non trivial
\\

  if( !odd && tryhard,
if( DEBUGLEVEL_ell >= 4, print(" Autre strategie dans ratpoint **********"));
    K = content(pol);
    if( K != 1,
      pol /= K;
      factK = factor(K);
      e = factK[,2]\2;
      cont = factorback(factK[,1],e);
      K /= cont^2;
      if(K != 1,
        e = factK[,2]%2;
        ind = #e; while( !e[ind], ind--);
        p = factK[ind,1];
        if( valuation( pollead(pol), p) == 1 ||
            ( valuation( pollead(pol), p) >= 2 && valuation( polcoeff(pol,poldegree(pol)-1), p) == 0),
if( DEBUGLEVEL_ell >= 4, print(" utilise une racine de pol mod p = ",p));
          sol = ratpoint(K/p^2*subst(polrecip(pol),variable(pol),p*variable(pol)),lim,singlepoint,1);
          if( #sol > 0,
            point1 = [ 1/(sol[1]*p),
              sol[2]*cont*p/(p*sol[1])^(poldegree(pol)/2) ];
if( DEBUGLEVEL_ell >= 4, print("sortie de ratpoint ",point1));
            return(point1)
          )
        );
        factpol = factormod(pol,p)[,1];
        for( i = 1, #factpol,
          if( poldegree(factpol[i]) !=1, next);
if( DEBUGLEVEL_ell >= 4, print(" utilise une racine de pol mod p = ",p));
          r = -centerlift(polcoeff(factpol[i],0));
          if( valuation(subst(pol,variable(pol),r),p) > 2, next);
          M = [p,r;0,1];
          U = redquartique(subst(K*pol,variable(pol),p*variable(pol)+r));
          if( content(U[1]) != p, next);
          sol = ratpoint(K/p^2*U[1],lim,singlepoint,1);
          if( #sol > 0,
            M = (M*U[2])*[sol[1], #sol == 2]~;
            point1 = [ M[1]/M[2], sol[2]*cont*p/M[2]^(poldegree(pol)/2) ];
if( DEBUGLEVEL_ell >= 4, print("sortie de ratpoint ",point1));
            return(point1)
          )
        )
      )
    )
  );

  return([]);
}
if( DEBUGLEVEL_ell >= 5, print("psquare"));
{
psquare( a, p) =
local(ap,v);

if( DEBUGLEVEL_ell >= 5, print("entree dans psquare ",[a,p]));
\\ a est un entier
\\ renvoie 1 si a est un carre dans Zp 0 sinon
  if( a == 0,
if( DEBUGLEVEL_ell >= 5, print("fin de psquare 1"));
    return(1));
\\
  v = valuation(a,p);
  if( v%2,
if( DEBUGLEVEL_ell >= 5, print("fin de psquare 0"));
    return(0));
  if( p == 2,
    ap = (a>>v)%8-1,
    ap = kronecker(a/p^v,p)-1
  );
if( DEBUGLEVEL_ell >= 5, print("fin de psquare ", !ap));
  return(!ap);
}
if( DEBUGLEVEL_ell >= 4, print("lemma6"));
{
lemma6(pol, p, nu, xx) =
local(gx,gpx,lambda,mu);

\\ pour les p <> 2
  gx = subst( pol, variable(pol), xx);
  if( psquare(gx,p), return(1));
  gpx = subst( pol', variable(pol), xx);
  lambda = valuation(gx,p); mu = valuation(gpx,p);

  if( lambda > 2*mu, return(1));
\\  if( (lambda >= mu+nu) && (nu > mu), return(1));
  if( (lambda >= 2*nu) && (mu >= nu), return(0));
  return(-1);
}
if( DEBUGLEVEL_ell >= 4, print("lemma7"));
{
lemma7( pol, nu, xx) =
local(gx,gpx,lambda,mu,q);

\\ pour p = 2
  gx = subst( pol, variable(pol), xx);
  if( psquare(gx,2), return(1));
  gpx = subst( pol', variable(pol), xx);
  lambda = valuation(gx,2); mu = valuation(gpx,2);
  if( lambda > 2*mu, return(1));
  if( nu > mu,
    if( lambda%2, return(-1));
    q = mu+nu-lambda;
    if( q == 1, return(1));
    if( q == 2 && (gx>>lambda)%4 == 1, return(1));
    return(-1));
  q = lambda-2*nu;
  if( q >= 0, return(0));
  if( q == -2 && (gx>>lambda)%4 == 1, return(0));
  return(-1);
}
if( DEBUGLEVEL_ell >= 4, print("zpsoluble"));
{
zpsoluble(pol, p, nu, pnu, x0, pnup) =
local(result,pol2,fact,x1);

if( DEBUGLEVEL_ell >= 5, print("entree dans zpsoluble ",[pol,p,x0,nu]));
  if( p == 2,
    result = lemma7(pol,nu,x0),
    result = lemma6(pol,p,nu,x0));
  if( result == +1,
if( DEBUGLEVEL_ell >= 5, print("fin de zpsoluble 1 lemma"));
    return(1));
  if( result == -1,
if( DEBUGLEVEL_ell >= 5, print("fin de zpsoluble 0 lemma"));
    return(0));
  pnup = pnu*p;
  nu++;
  if( p< LIMBIGPRIME || !LIMBIGPRIME,
    for( i = 0, p-1,
      if( zpsoluble(pol,p,nu,pnup,x0+pnu*i),
if( DEBUGLEVEL_ell >= 5, print("fin de zpsoluble"));
        return(1)))
  ,
    pol2 = subst(pol,variable(pol),x0+pnu*variable(pol));
    pol2 /= content(pol2);
    pol2 = pol2*Mod(1,p);
    if( !poldegree(pol2), return(0));
    fact = factormod(pol2,p)[,1];
    for( i = 1, #fact,
      x1 = -centerlift(polcoeff(fact[i],0));
      if( zpsoluble(pol,p,nu,pnup,x0+pnu*x1),
if( DEBUGLEVEL_ell >= 5, print("fin de zpsoluble"));
        return(1)));
    for( i = 1, MAXPROB,
      x1 = random(p);
      if( zpsoluble(pol,p,nu,pnup,x0+pnu*x1),
if( DEBUGLEVEL_ell >= 5, print("fin de zpsoluble"));
        return(1)))
  );
if( DEBUGLEVEL_ell >= 2,
  if( p >= LIMBIGPRIME,
    print("******* test probabiliste en p = ",p,"*******")));
if( DEBUGLEVEL_ell >= 5, print("fin de zpsoluble"));
  return(0);
}
if( DEBUGLEVEL_ell >= 4, print("qpsoluble"));
{
qpsoluble(pol, p) =
if( DEBUGLEVEL_ell >= 5, print("entree dans qpsoluble ",p); print("pol = ",pol));
  if( psquare(pollead(pol),p),
if( DEBUGLEVEL_ell >= 5, print("fin de qpsoluble 1"));
    return(1));
  if( psquare(polcoeff(pol,0),p),
if( DEBUGLEVEL_ell >= 5, print("fin de qpsoluble 1"));
    return(1));
  if( zpsoluble(pol,p,0,1,0),
if( DEBUGLEVEL_ell >= 5, print("fin de qpsoluble 1"));
    return(1));
  if( zpsoluble(polrecip(pol),p,1,p,0),
if( DEBUGLEVEL_ell >= 5, print("fin de qpsoluble 1"));
    return(1));
if( DEBUGLEVEL_ell >= 5, print("fin de qpsoluble 0"));
  return(0);
}
if( DEBUGLEVEL_ell >= 4, print("locallysoluble"));
{
locallysoluble(pol) =
\\ teste l'existence locale de solutions de y^2 = pol(x,z)
local(plist,disc0,p,c,vc);

if( DEBUGLEVEL_ell >= 4, print("entree dans locallysoluble :",pol));

\\ place reelle
  if( !(poldegree(pol)%2) && sign(pollead(pol)) < 0
         && sign(polcoeff(pol,0)) < 0 && polsturm(pol) == 0,
if( DEBUGLEVEL_ell >= 3, print(" non ELS a l'infini"));
if( DEBUGLEVEL_ell >= 4, print("fin de locallysoluble"));
    return(0));

\\
\\ places finies de plist */
\\
  pol *= denominator(content(pol))^2;
  c = content(pol);

  disc0 = poldisc(pol);
  plist = factor (abs(2*disc0));
if( DEBUGLEVEL_ell >= 4, print("liste de premiers = ",plist));
  for( i = 1, #plist[,1],
    p = plist[i,1];
if( DEBUGLEVEL_ell >= 4, print("p = ",p));
    vc = valuation(c,p);
    if( vc >= 2,
      pol /= p^(2*(vc\2));
      plist[i,2] -= 2*(vc\2)*(2*poldegree(pol)-2)
    );
    if( poldegree(pol) == 4 && p != 2 && plist[i,2] < 2, next);
    if( !qpsoluble(pol,p),
if( DEBUGLEVEL_ell >= 3, print(" non ELS en ",p));
if( DEBUGLEVEL_ell >= 4, print("fin de locallysoluble"));
      return(0)));

if( DEBUGLEVEL_ell >= 2, print(" quartique ELS"));
if( DEBUGLEVEL_ell >= 4, print("fin de locallysoluble"));
  return(1);
}
if( DEBUGLEVEL_ell >= 4, print("redquartique"));
{
redquartique(pol) =
\\ reduction d'une quartique.
\\ ou plus generalement d'un polynome de degre deg.
local(prec,prec0,d,disc2,test,normderiv,disc2v,r,q,M);

if( DEBUGLEVEL_ell >= 4, print("entree dans redquartique"));
if( DEBUGLEVEL_ell >= 3, print(" reduction de la quartique ",pol));

\\ choix de la precision des calculs.
  prec = prec0 = default(realprecision);
  d = poldegree(pol);
  disc2 = poldisc(pol)^2;
  test = 0;
  while( test == 0,
if( DEBUGLEVEL_ell >= 4, print(" precision = ",prec));
    r = polroots(pol);
    normderiv = vector( d, i, norm(subst(pol',variable(pol),r[i])));
    disc2v = prod( i = 1, d, normderiv[i]) * pollead(pol)^(2*d-4);
    test = abs(disc2v-disc2) < 10^(-prec\2);
    if( !test, default(realprecision, prec *= 2))
  );

\\ On n'utilise plus
\\  q = Vec(sum( i = 1, d, norm(x-r[i])));
\\ mais la normalisation de Cremona-Stoll
  q = Vec(sum( i = 1, d, norm(x-r[i]) / normderiv[i]^(1/(d-2))));
  M = QfbReduce([q[1],q[2]/2;q[2]/2,q[3]]);
  pol = subst(pol,variable(pol),Pol(M[1,])/Pol(M[2,]))*Pol(M[2,])^poldegree(pol);

  if( prec != prec0, default(realprecision,prec0));

if( DEBUGLEVEL_ell >= 3, print(" quartique reduite = ",pol));
if( DEBUGLEVEL_ell >= 4, print("sortie de redquartique"));

  return([pol,M]);
}
if( DEBUGLEVEL_ell >= 4, print("reducemodsquares"));
{
reducemodsquares(delta,d) =
\\ reduction du coefficient de x^d dans ( delta modulo les carres )
local(deg,xx,z,qd,Qd,reduc);

  deg = poldegree(delta.mod);
  xx = Mod(x,delta.mod);
  z = subst(Pol(vector(deg,i,eval(Str("a"i)))),x,xx);
  qd = polcoeff(lift(delta*z^2),d,x);
  Qd = simplify(matrix(deg,deg,i,j,deriv(deriv(qd,eval(Str("a"i))),eval(Str("a"j)))/2));

  reduc = IndefiniteLLL(Qd);
  if( #reduc == 2, reduc = reduc[2][,1]);

  return(delta*subst(Pol(reduc),x,xx)^2);
}
if( DEBUGLEVEL_ell >= 4, print("ellsort"));
{
ellsort(listpts) =
\\ tri des points listpts sur une courbe elliptique
\\ suivant la hauteur naive.
local(n,v,aux,ord);

  v = vector(n = #listpts);
  for( i = 1, n,
    if( listpts[i] == [0], v[i] = [0,0,0]; next);
    aux = denominator(listpts[i][2])/denominator(listpts[i][1]);
    v[i] = vecsort(abs([listpts[i][1]*aux^2, listpts[i][2]*aux^3,aux]),,4);
  );
  ord = vecsort(v,,3);
  return(vector(n,i,listpts[ord[i]]));
}
if( DEBUGLEVEL_ell >= 4, print("ellredgen"));
{
ellredgen(ell,listgen,K=1) =
\\ reduction des generateurs de listgen
\\ sur la courbe ell = [a1,a2,a3,a4,a6]
\\ ou K*y^2 = x^3 + a2*x^2 + a4*x + a6 (lorsque a1 = a3 = 0);
local(d,sqrtK,urst,M,U,limgoodrelations,listgen2);

if( DEBUGLEVEL_ell >= 3, print("Reduction des generateurs ",listgen));
if( DEBUGLEVEL_ell >= 5, print("ell=",ell));
  d = #listgen;
  if( d == 0, return([]));

  if( #ell < 19, ell = ellinit(ell));
  if( K != 1,
    if( ell.a1 != 0 || ell.a3 != 0, error(" ellredgen : a1*a3 != 0"));
    ell[2] *= K; ell[4] *= K^2; ell[5] *= K^3;
    ell[6] *= K; ell[7] *= K^2; ell[8] *= K^3; ell[9] *= K^4;
    ell[10] *= K^2; ell[11] *= K^3; ell[12] *= K^6;
    sqrtK = sqrt(K);
    ell[14] *= K;
    ell[15] /= sqrtK; ell[16] /= sqrtK;
    ell[17] *= sqrtK; ell[18] *= sqrtK;
    ell[19] /= K;

    for( i = 1, d,
      for( j = 1, #listgen[i],
        listgen[i][j] *= K^j))
  );

  ell = ellminimalmodel(ell,&urst);
  listgen = ellchangepoint(listgen,urst);

\\ Recherche des relations entre les points de listgen
\\ par recherche du noyau de la matrice des hauteurs

  M = ellheightmatrix(ell,listgen);
if( DEBUGLEVEL_ell >= 4, print("matrice des hauteurs = ",M));
  M = round( M*10^(default(realprecision)-10) );
  U = qflll(M,4);
  U = concat(U[1],U[2]);

  /* BEGIN patch to work with PARI 2.4.4 */
  /* AUTHORS: John Cremona, Jeroen Demeyer (Sage Trac #11130) */
if( DEBUGLEVEL_ell >= 4, print("    change of basis proposed by LLL = ",U));
  \\ The columns of U that have very small coefficients (coeff < 20)
  \\ are either exact relations or reductions.  These are the ones we
  \\ want to keep, the other ones are irrelevant.
  keep = 0;
  for( i = 1, d,
    if( vecmax(abs(U[,i])) < 20, keep += 1<<(i-1))
  );
  U = vecextract(U, keep);
  /* END patch from Sage Ticket #11130 to work with PARI 2.4.4 */

  U = completebasis(U);
if( DEBUGLEVEL_ell >= 4, print("changement de base = ",U));

  listgen2 = vector(d);
  for( i = 1, d,
    listgen2[i] = [0];
    for( j = 1, d,
      listgen2[i] = elladd(ell,listgen2[i],ellpow(ell,listgen[j],U[j,i]))));
  listgen = listgen2;

\\ Tri des points d'ordre infini

  listgen2 = [];
  for( i = 1, d,
    if( !ellorder(ell,listgen[i]),
      listgen2 = concat(listgen2,[listgen[i]])));
  listgen = listgen2;
if( DEBUGLEVEL_ell >= 3, print("points d'ordre infini = ",listgen));
  d = #listgen;
  if( d == 0, return([]));

\\ Reduction des points d'ordre infini

  if( d > 1,
    M = ellheightmatrix(ell,listgen);
if( DEBUGLEVEL_ell >= 4, print("matrice des hauteurs = ",M));
    U = qflllgram(M);
if( DEBUGLEVEL_ell >= 4, print("changement de base = ",U));

    listgen2 = vector(d);
    for( i = 1, d,
      listgen2[i] = [0];
      for( j = 1, d,
        listgen2[i] = elladd(ell,listgen2[i],ellpow(ell,listgen[j],U[j,i]))));
    listgen = listgen2
  );

  listgen = ellsort(listgen);
if( DEBUGLEVEL_ell >= 4, print("generateurs tries = ",listgen));

  listgen = ellchangepointinverse(listgen,urst);
  if( K != 1,
    for( i = 1, d,
      for( j = 1, 2,
        listgen[i][j] /= K^j)));

\\ on ne garde que les points (x,y) avec y >= 0

  if( ell.a1 == 0 && ell.a3 == 0,
    for( i = 1, d,
      if( #listgen[i] == 2,
        listgen[i][2] = abs(listgen[i][2]))));

if( DEBUGLEVEL_ell >= 2, print("generateurs reduits = ",listgen));
  return(listgen);
}
if( DEBUGLEVEL_ell >= 4, print("ell2descent_gen"));
{
ell2descent_gen(ell,ext,K=1,help=[],redflag=0) =
\\ si ell= K*y^2=P(x), alors ext est le buchinitfu de l'extension.
\\ theta est une racine de P.
\\ dans la suite ext est note L = Q(theta).
\\ help est une liste de points deja connus sur ell.
\\ ell est de la forme K*y^2=x^3+A*x^2+B*x+C */
\\ ie ell=[0,A,0,B,C], avec A,B et C entiers */
\\
\\ si redflag != 0, on utilise la reduction modulo les carres.
\\

local(A,B,C,S,SLprod,SLlist,aux,oddclass,LS2gen,fact,trouve,i,polrel,ttheta,polprime,KS2gen,LS2genunit,normLS2gen,normcoord,LS2coordtilda,LS2,listgen,listpoints,listpointstriv,listpointsmwr,list,m1,m2,lastloc,maskwhile,iwhile,zc,j,iaux,liftzc,den,ispointtriv,point,found,idealfactorzc,idealzc,baseidealzc,q2,sol,q1,param,pol,redq,q0,pointxx,point2,v,rang);

if( DEBUGLEVEL_ell >= 4, print("entree dans ell2descent_gen"));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      construction de L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  if( #ell < 13, ell = ellinit(ell,1));

  if( ell.a1 != 0 || ell.a3 != 0,
    error(" ell2descent_gen : la courbe n'est pas sous la forme [0,a,0,b,c]"));
  if( denominator(ell.a2) > 1 || denominator(ell.a4) > 1 || denominator(ell.a6) >1,
    error(" ell2descent_gen : coefficients non entiers"));

  A = ell.a2; if( DEBUGLEVEL_ell >= 2, print("A = ",A));
  B = ell.a4; if( DEBUGLEVEL_ell >= 2, print("B = ",B));
  C = ell.a6; if( DEBUGLEVEL_ell >= 2, print("C = ",C));

  polrel = Pol([1,A,B,C]);
  if( !ext,
if( DEBUGLEVEL_ell >= 2, print("bnfinit(",polrel,")"));
    ext = bnfinit(polrel,1));

  S = -abs(K*ext.index);
  SLprod = idealmul(ext,K,idealadd(ext,ext.pol',ext.index));
  SLlist = idealfactor(ext,SLprod)[,1]~;
  aux = []; SLprod = 1;
   for( i = 1, #SLlist,
     if( !(K%SLlist[i][1]) || valuation(ext.index,SLlist[i][1]),
       SLprod = idealmul(ext,SLprod,SLlist[i]);
       aux = concat(aux,[SLlist[i]])));
  SLlist = aux;
  oddclass = 0;
  while( !oddclass,
\\ Constructoin de S:
if( DEBUGLEVEL_ell >= 4, print("SLlist = ",SLlist));
\\ Construction des S-unites
    LS2gen = bnfsunit(ext,SLlist);
if( DEBUGLEVEL_ell >= 4, print("LS2gen = ",LS2gen));
\\ on ajoute la partie paire du groupe de classes.
    oddclass = LS2gen[5][1]%2;
    if( !oddclass,
if( DEBUGLEVEL_ell >= 3, print("Groupe de classes pair"));
if( DEBUGLEVEL_ell >= 4, print(LS2gen[5]));
      S *= LS2gen[5][3][1][1,1];
      SLprod = idealmul(ext,SLprod,LS2gen[5][3][1]);
      fact = idealfactor(ext,LS2gen[5][3][1])[,1];
      trouve = 0; i = 0;
      while( !trouve,
        i++; trouve = 1;
        for( j = 1, #SLlist,
          if( SLlist[j] == fact[i], trouve = 0; break)));
      SLlist = concat(SLlist,[fact[i]]))
  );

if( DEBUGLEVEL_ell >= 4, print("S = ",S));

  ttheta = Mod(x,polrel);            \\ ttheta est la racine de P(x)
  polprime = Mod(polrel',polrel);

  KS2gen = factor(S)[,1]~;

if( DEBUGLEVEL_ell >= 3, print("#KS2gen = ",#KS2gen));
if( DEBUGLEVEL_ell >= 3, print("KS2gen = ",KS2gen));

  LS2genunit = ext.tufu;
  LS2genunit = concat(LS2gen[1],LS2genunit);

  LS2genunit = subst(LS2genunit,x,ttheta);
  LS2genunit = LS2genunit*Mod(1,polrel);
if( DEBUGLEVEL_ell >= 3, print("#LS2genunit = ",#LS2genunit));
if( DEBUGLEVEL_ell >= 3, print("LS2genunit = ",LS2genunit));

\\ dans LS2gen, on ne garde que ceux dont la norme
\\ est un carre.

  normLS2gen = norm(LS2genunit);
if( DEBUGLEVEL_ell >= 4, print("normLS2gen = ",normLS2gen));

\\ matrice de l'application norme

  normcoord = matrix(#KS2gen,#normLS2gen);
  for( j = 1, #normLS2gen,
    normcoord[1,j] = (sign(normLS2gen[j]) < 0);
    for( i = 2, #KS2gen,
      normcoord[i,j] = valuation(normLS2gen[j],KS2gen[i])));
if( DEBUGLEVEL_ell >= 4, print("normcoord = ",normcoord));

\\ construction du noyau de la norme

  LS2coordtilda = lift(matker(normcoord*Mod(1,2)));
if( DEBUGLEVEL_ell >= 4, print("LS2coordtilda = ",LS2coordtilda));
  LS2 = vector(#LS2coordtilda[1,],i,0);
  for( i = 1, #LS2coordtilda[1,],
    aux = 1;
    for( j = 1, #LS2coordtilda[,i],
      if( sign(LS2coordtilda[j,i]),
        aux *= LS2genunit[j]));
    LS2[i] = aux
  );
if( DEBUGLEVEL_ell >= 3, print("LS2 = ",LS2));
if( DEBUGLEVEL_ell >= 3, print("norm(LS2) = ",norm(LS2)));

\\ Reduction des generateurs de LS2

  if( redflag,
    for( i = 1, #LS2,
      LS2[i] = reducemodsquares(LS2[i],2)));

\\ Fin de la construction de LS2

  listgen = LS2;
if( DEBUGLEVEL_ell >= 2, print("LS2gen = ",listgen));
if( DEBUGLEVEL_ell >= 2, print("#LS2gen = ",#listgen));
  listpoints = [];

if( DEBUGLEVEL_ell >= 3, print("(A,B,C) = ",[A,B,C]));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Recherche de points triviaux.  \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if( DEBUGLEVEL_ell >= 2, print("Recherche de points triviaux sur la courbe"));
  listpointstriv = ratpoint(K^3*subst(polrel,x,x/K),LIMTRIV,0,0);
  for( i = 1, #listpointstriv,
    if( #listpointstriv[i] == 3,
      listpointstriv[i] = [0]
    , for( j = 1, 2, listpointstriv[i][j] /= K^j))
   );
  listpointstriv = concat(help,listpointstriv);
if( DEBUGLEVEL_ell >= 1, print("Points triviaux sur la courbe = ",listpointstriv));


\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          parcours de L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  listpointsmwr = [];
  list = [ 6, ell.disc, 0 ];
  m1 = 0; m2 = 0; lastloc = -1;
  maskwhile = 1<<#listgen;
  iwhile = 1;
  while( iwhile < maskwhile,
if( DEBUGLEVEL_ell >= 4, print("iwhile = ",iwhile); print("listgen = ",listgen));
    zc = Mod(1,polrel); j = 1; iaux = iwhile;
    while( iaux,
      if( iaux%2, zc *= listgen[j]);
      iaux >>= 1; j++);
if( DEBUGLEVEL_ell >= 2, print(); print("zc = ",zc));
    liftzc = lift(zc);
    if( redflag,
      zc = reducemodsquares(zc,2);
      liftzc = lift(zc);
      den = denominator(content(liftzc))^2;
      zc *= den; liftzc *= den;
if( DEBUGLEVEL_ell >= 2, print("zcred = ",zc))
    );

\\ Est-ce un point trivial ?
    ispointtriv = 0;
    for( i = 1, #listpointstriv,
      point = listpointstriv[i];
      if( #point == 2,
        if( nfissquare(ext.nf,K*(point[1]-x)*liftzc),
if( DEBUGLEVEL_ell >= 2, print(" vient du point trivial ",point));
          listpointsmwr = concat(listpointsmwr,[point]);
          m1++;
          if( degre(iwhile) > lastloc, m2++);
          found = (ispointtriv = 1);
          break
    )));

\\ Ce n'est pas un point trivial
    if( !ispointtriv,
\\ Il faut resoudre une forme quadratique
\\      q2 = matrix(3,3,i,j,trace(zc*ttheta^(i+j-2)/polprime));
\\if( DEBUGLEVEL_ell >= 4, print("q2 = ",q2));
      idealfactorzc = idealfactor(ext,zc);
      idealfactorzc[,2] *= -1;
      idealfactorzc[,2] \= 2;
      idealzc = idealfactorback(ext,idealfactorzc);
      if( idealzc == 1, idealzc = matid(3));
      baseidealzc = vector(3,i,nfbasistoalg(ext,idealzc[,i]));
      q2 = matrix(3,3,i,j,trace(zc*baseidealzc[i]*baseidealzc[j]/polprime));
if( DEBUGLEVEL_ell >= 4, print("q2 = ",q2));
\\      q2 *= ext.index;
\\ if( DEBUGLEVEL_ell >= 4, print("q2 = ",q2));
if( DEBUGLEVEL_ell >= 4, print("q2/content(q2) = ",q2/content(q2)));
      sol = Qfsolve(q2/content(q2));
if( DEBUGLEVEL_ell >= 4,print("sol = ",sol));
      if( type(sol) == "t_INT",
if( DEBUGLEVEL_ell >= 3, print("non ELS en ",sol));
        iwhile++; next
      );

\\ \\\\\\\\\\\
\\ Construction de la quartique
\\ \\\\\\\\\\\
      q1 = -matrix(3,3,i,j,trace(zc*baseidealzc[i]*baseidealzc[j]*(ttheta+A)/polprime));
      param = Qfparam(q2,sol)*[x^2,x,1]~;
      param /= content(param);
      pol = param~*q1*param;
if( DEBUGLEVEL_ell >= 2, print(" quartique: ",K,"*Y^2 = ",pol));
      redq = redquartique(pol);
if( DEBUGLEVEL_ell >= 2, print(" reduite: ",K,"*Y^2 = ",redq[1]));
      pol = redq[1];
      den = denominator(content(K*pol));
      pol *= den^2;

\\ \\\\\\\\\\\
\\ Recherche de points sur la quartique
\\ \\\\\\\\\\\

      point = ratpoint(K*pol,LIM1,1,0);
      if( point != [],
        m1++;
        if( #point == 2, point = concat(point,[1]));
        point = concat(redq[2]*[point[1],point[3]]~,[point[2]/den]~);
        param = subst(param,x,x/y)*y^2;
        param = subst(subst(param,x,point[1]),y,point[2]);
if( DEBUGLEVEL_ell >= 2, print(" point sur la quartique = ",point));
        param *= K/point[3];
if( DEBUGLEVEL_ell >= 3, print("reconstruction du point sur la courbe"));
        q0 = matrix(3,3,i,j,trace(zc*baseidealzc[i]*baseidealzc[j]*(ttheta^2+A*ttheta+B)/polprime));
        pointxx = param~*q0*param/K;
        point2 = [ pointxx, sqrtrat(subst(x^3+A*x^2+B*x+C,x,pointxx)/K)];
if( DEBUGLEVEL_ell >= 1, print("point sur la courbe = ",point2));
        listpointsmwr = concat(listpointsmwr,[point2]);
        if( degre(iwhile) > lastloc, m2++);
        found = 1
      ,
        if( locallysoluble(K*pol),
          if( degre(iwhile) > lastloc, m2++; lastloc = degre(iwhile));
          point = ratpoint(K*pol,LIM3,1,1);
          if( #point > 0,
            m1++;
            if( #point == 2, point = concat(point,[1]));
            point = concat(redq[2]*[point[1],point[3]]~,[point[2]/den]~);
            param = subst(param,x,x/y)*y^2;
            param = subst(subst(param,x,point[1]),y,point[2]);
if( DEBUGLEVEL_ell >= 3, print(" point sur la quartique = ",point));
            param *= K/point[3];
if( DEBUGLEVEL_ell >= 3, print("reconstruction du point sur la courbe"));
            q0 = matrix(3,3,i,j,trace(zc*baseidealzc[i]*baseidealzc[j]*(ttheta^2+A*ttheta+B)/polprime));
            pointxx = param~*q0*param/K;
            point2 = [ pointxx, sqrtrat(subst(x^3+A*x^2+B*x+C,x,pointxx)/K)];
if( DEBUGLEVEL_ell >= 1, print("point sur la courbe = ",point2));
            listpointsmwr = concat(listpointsmwr,[point2]);
            found = 1
          )
        )
      )
    );
    if( found || ispointtriv,
      found = 0; lastloc = -1;
      v = 0; iaux = (iwhile>>1);
      while( iaux, iaux >>= 1; v++);
      maskwhile >>= 1;
      listgen = vecextract(listgen,1<<#listgen-1<<v-1);
      iwhile = 1<<v
    , iwhile ++
    )
  );

if( DEBUGLEVEL_ell >= 2,
  print();
  print("rang des points trouves = ",m1);
  print("rang du groupe de Selmer = ",m2));
if( DEBUGLEVEL_ell >= 1,
  print("#S(E/Q)[2]    = ",1<<m2));
  if( m1 == m2,
if( DEBUGLEVEL_ell >= 1,
    print("#E(Q)/2E(Q)   = ",1<<m1);
    print("#III(E/Q)[2]  = 1");
    print("rang(E/Q)     = ",m1));
    rang = m1
  ,
if( DEBUGLEVEL_ell >= 1,
    print("#E(Q)/2E(Q)  >= ",1<<m1);
    print("#III(E/Q)[2] <= ",1<<(m2-m1));
    print("rang(E/Q)    >= ",m1));
    rang = m1;
    if( (m2-m1)%2,
if( DEBUGLEVEL_ell >= 1,
      print(" III devrait etre un carre, donc ");
      if( m2-m1 > 1,
        print("#E(Q)/2E(Q)  >= ",1<<(m1+1));
        print("#III(E/Q)[2] <= ",1<<(m2-m1-1));
        print("rang(E/Q)    >= ",m1+1)
      ,
        print("#E(Q)/2E(Q)  = ",1<<(m1+1));
        print("#III(E/Q)[2] = 1");
        print("rang(E/Q)    = ",m1+1)));
        rang = m1+1
    )
  );
if( DEBUGLEVEL_ell >= 1, print("listpointsmwr = ",listpointsmwr));
  for( i = 1, #listpointsmwr,
    if( subst(polrel,x,listpointsmwr[i][1])-K*listpointsmwr[i][2]^2,
      error(" ell2descent_gen : MAUVAIS POINT = ",listpointsmwr[i])));
  if( #listpointsmwr >= 2,
    listpointsmwr = ellredgen(ell,listpointsmwr,K));
if( DEBUGLEVEL_ell >= 4, print("fin de ell2descent_gen"));
  return([rang,m2,listpointsmwr]);
}
if( DEBUGLEVEL_ell >= 4, print("ellrank"));
{
ellrank(ell,help=[]) =
\\ Algorithme de la 2-descente sur la courbe elliptique ell.
\\ help est une liste de points connus sur ell.
local(urst,urst1,den,eqell,tors2,bnf,rang,time1);

if( DEBUGLEVEL_ell >= 3, print("entree dans ellrank"));
  if( #ell < 13, ell = ellinit(ell,1));

\\ supprime les coefficients a1 et a3
  urst = [1,0,0,0];
  if( ell.a1 != 0 || ell.a3 != 0,
    urst1 = [1,0,-ell.a1/2,-ell.a3/2];
    ell = ellchangecurve(ell,urst1);
    urst = ellcomposeurst(urst,urst1)
  );

\\ supprime les denominateurs
  while( (den = denominator([ell.a2,ell.a4,ell.a6])) > 1,
    den = factor(den); den[,2] = vectorv(#den[,2],i,1);
    den = factorback(den);
    urst1 = [1/den,0,0,0];
    ell = ellchangecurve(ell,urst1);
    urst = ellcomposeurst(urst,urst1)
  );

  help = ellchangepoint(help,urst);
  eqell = Pol([1,ell.a2,ell.a4,ell.a6]);
if( DEBUGLEVEL_ell >= 1, print("courbe elliptique : Y^2 = ",eqell));

\\ choix de l'algorithme suivant la 2-torsion

  tors2 = ellhalf(ell,[0]);
if( DEBUGLEVEL_ell >= 1, print("E[2] = ",tors2));

  if( #tors2 == 1,                              \\ cas 1: 2-torsion triviale
if( DEBUGLEVEL_ell >= 3, print1("bnfinit "));
if( DEBUGLEVEL_ell >= 4, gettime());
    bnf = bnfinit(eqell,1);
if( DEBUGLEVEL_ell >= 4, time1 = gettime());
if( DEBUGLEVEL_ell >= 3, print("ok"));
    rang = ell2descent_gen(ell,bnf,1,help);
if( DEBUGLEVEL_ell >= 4, print("temps dans bnfinit = ",time1));
if( DEBUGLEVEL_ell >= 4, print("temps pour le reste = ",gettime()));
  ,
  if( #tors2 >= 2,                              \\ cas 2: 2-torsion >= Z/2Z
    if( ell.a6 != 0,
      urst1 = [1,tors2[2][1],0,0];
      ell = ellchangecurve(ell,urst1);
      urst = ellcomposeurst(urst,urst1)
    );
    eqell = Pol([1,ell.a2,ell.a4,ell.a6]);
if( DEBUGLEVEL_ell >= 1, print("courbe elliptique : Y^2 = ",eqell));

    rang = ell2descent_viaisog(ell,help)
  ,                                             \\ cas 3: 2-torsion = Z/2Z*Z/2Z
\\    rang = ell2descent_complete(tors2[2][1],tors2[3][2],tors2[4][3])
  ));

  rang[3] = ellchangepointinverse(rang[3],urst);
if( DEBUGLEVEL_ell >= 3, print("fin de ellrank"));

return(rang);
}
if( DEBUGLEVEL_ell >= 4, print("ell2descent_complete"));
{
ell2descent_complete(e1,e2,e3) =
\\ calcul du rang d'une courbe elliptique
\\ par la methode de 2-descente complete.
\\ Y^2 = (x-e1)*(x-e2)*(x-e3);
\\ en suivant la methode decrite par J.Silverman
\\ renvoie [r,s,v] avec
\\   r est une borne inferieure du rang de E(Q)
\\   s est le rang du 2-groupe de Selmer
\\   v est un systeme de points independants dans E(Q)/2E(Q)

\\ e1, e2 et e3 sont des entiers.

local(d32,d31,d21,G1,G2,vect1,vect2,selmer,rang,listepoints,b1,b2,q1,q2,sol1,param1,param1x,sol2,quart,point,z1,solx,soly,strange);

if( DEBUGLEVEL_ell >= 2, print("Algorithme de la 2-descente complete"));

\\ calcul des groupes G1 et G2

  d32 = e3-e2; d31 = e3-e1; d21 = e2-e1;
  G1 = factor(-abs(d31*d21))[,1];
  G2 = factor(-abs(d32*d21))[,1];

if( DEBUGLEVEL_ell >= 3, print("G1 = ",G1));
if( DEBUGLEVEL_ell >= 3, print("G2 = ",G2));

\\ parcours de G1*G2

  vect1 = vector(#G1,i,[0,1]);
  vect2 = vector(#G2,i,[0,1]);
  selmer = 0;
  rang = 0;
  listepoints = [];

  forvec( X = vect1,
    b1 = prod( i = 1, #G1, G1[i]^X[i]);
    forvec( Y = vect2,
      b2 = prod( i = 1, #G2, G2[i]^Y[i]);

if( DEBUGLEVEL_ell >= 3, print("[b1,b2] = ",lift([b1,b2])));

\\ points triviaux provenant de la 2-torsion

      if( b1==1 && b2==1,
if( DEBUGLEVEL_ell >= 4, print(" point trivial [0]"));
        selmer++; rang++; next);
      if( issquare(-d21*b2) && issquare(d31*d21*b1),
if( DEBUGLEVEL_ell >= 3, print(" point trivial [e1,0]"));
        selmer++; rang++; listepoints = concat(listepoints,[[e1,0]]); next);
      if( issquare(d21*b1) && issquare(-d32*d21*b2),
if( DEBUGLEVEL_ell >= 3, print(" point trivial [e2,0]"));
        selmer++; rang++; listepoints = concat(listepoints,[[e2,0]]); next);
      if( issquare(d31*b1) && issquare(d32*b2),
if( DEBUGLEVEL_ell >= 3, print(" point trivial [e3,0]"));
        selmer++; rang++; listepoints = concat(listepoints,[[e3,0]]); next);

\\ il faut resoudre 2 equations quadratiques:
\\ (1) b1*z1^2-b2*z2^2 = e2-e1
\\ (2) b1*z1^2-b1*b2*z3^2 = e3-e1
\\ on aura alors (x,y) = (b1*z1^2+e1,b1*b2*z1*z2*z3)

      q1 = matdiagonal([b1,-b2,-d21]);
if( DEBUGLEVEL_ell >= 3, print(" q1 = ",q1));
      q2 = matdiagonal([b1,-b1*b2,-d31]);
if( DEBUGLEVEL_ell >= 3, print(" q2 = ",q2));

\\ solution de la premiere forme quadratique

      sol1 = Qfsolve(q1);
      if( type(sol1) == "t_INT",
if( DEBUGLEVEL_ell >= 3, print(" q1 non ELS en ",sol1));
        next);
if( DEBUGLEVEL_ell >= 3, print(" sol part de q1 = ",sol1));
      param1 = Qfparam(q1,sol1,1);
if( DEBUGLEVEL_ell >= 3, print(" param de q1 = ",param1));
      param1x = param1*[x^2,x,1]~;

      sol2 = Qfsolve(q2);
      if( type(sol2) == "t_INT",
if( DEBUGLEVEL_ell >= 3, print(" q2 non ELS en ",sol2));
        next);

      quart = b1*b2*(b1*param1x[1]^2-d31*param1x[3]^2);
if( DEBUGLEVEL_ell >= 3, print(" quart = ",quart));

\\ la quartique est-elle localement soluble ?

      if( !locallysoluble(quart),
if( DEBUGLEVEL_ell >= 3, print(" quartique non ELS "));
        next);
if( DEBUGLEVEL_ell >= 2, print(" y^2 = ",quart));
      selmer++;

\\ recherche de points sur la quartique.

      point = ratpoint(quart,LIM3,1);
      if( point != [],
if( DEBUGLEVEL_ell >= 2, print("point trouve sur la quartique !!"));
if( DEBUGLEVEL_ell >= 3, print(point));
        if( #point == 2,
          z1 = subst(param1x[1],x,point[1])/subst(param1x[3],x,point[1])
        , z1 = param1[1,1]/param1[3,1]);
        solx = b1*z1^2+e1;
        soly = sqrtrat((solx-e1)*(solx-e2)*(solx-e3));
        listepoints = concat(listepoints,[[solx,soly]]);
if( DEBUGLEVEL_ell >= 1, print("point sur la courbe elliptique = ",[solx,soly]));
        rang++
      ,
if( DEBUGLEVEL_ell >= 2, print("aucun point trouve sur la quartique"))
      )
    )
  );

\\ fin

if( DEBUGLEVEL_ell >= 1,
  print("#S^(2)      = ",selmer));
  if( rang > selmer/2, rang = selmer);
if( DEBUGLEVEL_ell >= 1,
  strange = rang != selmer;
  if( strange,
  print("#E[K]/2E[K]>= ",rang)
, print("#E[K]/2E[K] = ",rang));
  print("#E[2]       = 4");
);
  rang = ceil(log(rang)/log(2))-2;
  selmer = valuation(selmer,2);
if( DEBUGLEVEL_ell >= 1,
  if( strange,
  print(selmer-2," >= rang  >= ",rang)
, print("rang        = ",rang));
  if( rang, print("points = ",listepoints));
);
  ell = ellinit([0,-(e1+e2+e3),0,e1*e2+e2*e3+e3*e1,-e1*e2*e3]);
  if( ELLREDGENFLAG, listepoints = ellredgen(ell,listepoints));
  listepoints = concat(ellsort(elltorseven(ell)[3]),listepoints);

return([rang,selmer,listepoints]);
}
if( DEBUGLEVEL_ell >= 4, print("ellcount"));
{
ellcount( c, d, KS2gen, listpointstriv=[]) =
local(found,listgen,listpointscount,m1,m2,lastloc,mask,i,d1,iaux,j,triv,pol,point,deuxpoints,v);

if( DEBUGLEVEL_ell >= 4, print("entree dans ellcount ",[c,d]));
if( DEBUGLEVEL_ell >= 4, print("KS2gen = ",KS2gen));
if( DEBUGLEVEL_ell >= 4, print("listpointstriv = ",listpointstriv));

  found = 0;
  listgen = KS2gen;
  listpointscount = [];

  m1 = m2 = 0; lastloc = -1;

  mask = 1 << #KS2gen;
  i = 1;
  while( i < mask,
    d1 = 1; iaux = i; j = 1;
    while( iaux,
      if( iaux%2, d1 *= listgen[j]);
      iaux >>= 1; j++);
if( DEBUGLEVEL_ell >= 2, print("d1 = ",d1));
    triv = 0;
    for( j = 1, #listpointstriv,
      if( listpointstriv[j][1] && issquare(d1*listpointstriv[j][1]),
        listpointscount = concat(listpointscount,[listpointstriv[j]]);
if( DEBUGLEVEL_ell >= 2, print("point trivial"));
        triv = 1; m1++;
        if( degre(i) > lastloc, m2++);
        found = 1; lastloc = -1; break));
    if( !triv,
    pol = Pol([d1,0,c,0,d/d1]);
if( DEBUGLEVEL_ell >= 3, print("quartique = y^2 = ",pol));
    point = ratpoint(pol,LIM1,1);
    if( point != [],
if( DEBUGLEVEL_ell >= 2, print("point sur la quartique"));
if( DEBUGLEVEL_ell >= 3, print(point));
      m1++;
      listpointscount = concat(listpointscount,[d1*point[1]*point]);
      if( degre(i) > lastloc, m2++);
      found = 1; lastloc = -1
    ,
      if( locallysoluble(pol),
        if( degre(i) > lastloc, m2++; lastloc = degre(i));
        point = ratpoint(pol,LIM3,1);
        if( point != [],
if( DEBUGLEVEL_ell >= 2, print("point sur la quartique"));
if( DEBUGLEVEL_ell >= 3, print(point));
          m1++;
          listpointscount = concat(listpointscount,[d1*point[1]*point]);
          if( degre(i) > lastloc, m2++);
          found = 1; lastloc = -1
        ,
if( DEBUGLEVEL_ell >= 2, print(" pas de point trouve sur la quartique"))
          ))));
    if( found,
      found = 0;
      v = 0; iaux = (i>>1);
      while( iaux, iaux >>= 1; v++);
      mask >>= 1;
      listgen = vecextract(listgen,(1<<#listgen)-(1<<v)-1);
      i = (1<<v)
    , i++)
  );

  for( i = 1, #listpointscount,
   if( #listpointscount[i] > 1,
    if( subst(x^3+c*x^2+d*x,x,listpointscount[i][1])-listpointscount[i][2]^2 != 0,
      error(" ellcount : MAUVAIS POINT "))));
if( DEBUGLEVEL_ell >= 4, print("fin de ellcount"));
  return([listpointscount,[m1,m2]]);
}
if( DEBUGLEVEL_ell >= 4, print("ell2descent_viaisog"));
{
ell2descent_viaisog(ell,help=[]) =
\\ Calcul du rang des courbes elliptiques avec 2-torsion
\\ par la methode des 2-isogenies.
\\
\\ ell = [a1,a2,a3,a4,a6]
\\ y^2+a1xy+a3y=x^3+a2x^2+a4x+a6
\\
\\ ell doit etre sous la forme
\\ y^2=x^3+ax^2+bx -> ell = [0,a,0,b,0]
\\ avec a et b entiers.
local(i,P,Pfact,tors,listpointstriv,apinit,bpinit,plist,KS2prod,oddclass,KS2gen,listpoints,pointgen,n1,n2,certain,np1,np2,listpoints2,aux1,aux2,certainp,rang,strange);

if( DEBUGLEVEL_ell >= 2, print("Algorithme de la 2-descente par isogenies"));
  if( #ell < 13, ell = ellinit(ell,1));

  if( ell.disc == 0,
    error(" ell2descent_viaisog : courbe singuliere !!"));
  if( ell.a1 != 0 || ell.a3 != 0 || ell.a6 != 0,
    error(" ell2descent_viaisog : la courbe n'est pas sous la forme [0,a,0,b,0]"));
  if( denominator(ell.a2) > 1 || denominator(ell.a4) > 1,
    error(" ell2descent_viaisog : coefficients non entiers"));

  P = Pol([1,ell.a2,ell.a4]);
  Pfact = factor(P)[,1];
  tors = #Pfact;
  listpointstriv = concat(help,elltorseven(ell)[3]);

  apinit = -2*ell.a2; bpinit = ell.a2^2-4*ell.a4;

  plist = factor(6*ell.disc)[,1];

if( DEBUGLEVEL_ell >= 3, print(" Recherche de points triviaux sur la courbe"));
  P *= x;
if( DEBUGLEVEL_ell >= 3, print("Y^2 = ",P));
  listpointstriv = concat( listpointstriv, ratpoint(P,LIMTRIV,0));
if( DEBUGLEVEL_ell >= 1, print("points triviaux sur E(Q) = ",listpointstriv); print());

  KS2prod = -abs(ell.a4);
  if( bpinit < 0, KS2prod = -KS2prod);
  KS2gen = factor(KS2prod)[,1];

if( DEBUGLEVEL_ell >= 2,
  print("#K(b,2)gen          = ",#KS2gen);
  print("K(b,2)gen = ",KS2gen));

  listpoints = ellcount(ell.a2,ell.a4,KS2gen,listpointstriv);
  pointgen = listpoints[1];
if( DEBUGLEVEL_ell >= 1, print("points sur E(Q) = ",pointgen); print());
  n1 = listpoints[2][1]; n2 = listpoints[2][2];

  certain = (n1 == n2);
if( DEBUGLEVEL_ell >= 1,
  if( certain,
    print("[E(Q):phi'(E'(Q))]  = ",1<<n1);
    print("#S^(phi')(E'/Q)     = ",1<<n2);
    print("#III(E'/Q)[phi']    = 1"); print()
  ,
    print("[E(Q):phi'(E'(Q))] >= ",1<<n1);
    print("#S^(phi')(E'/Q)     = ",1<<n2);
    print("#III(E'/Q)[phi']   <= ",1<<(n2-n1)); print())
);

  KS2prod = -abs(bpinit);
  if( ell.a4 < 0, KS2prod = -KS2prod);
  KS2gen = factor(KS2prod)[,1];

if( DEBUGLEVEL_ell >= 2,
  print("#K(a^2-4b,2)gen     = ",#KS2gen);
  print("K(a^2-4b,2)gen     = ",KS2gen));

  P = Pol([1,apinit,bpinit]);
  listpointstriv = elltorseven([0,apinit,0,bpinit,0])[3];

if( DEBUGLEVEL_ell >= 3, print(" Recherche de points triviaux sur la courbe"));
  P *= x;
if( DEBUGLEVEL_ell >= 3, print("Y^2 = ",P));
  listpointstriv = concat( listpointstriv, ratpoint(P,LIMTRIV,0));
if( DEBUGLEVEL_ell >= 1, print("points triviaux sur E'(Q) = ",listpointstriv); print());

  listpoints = ellcount(apinit,bpinit,KS2gen,listpointstriv);
if( DEBUGLEVEL_ell >= 1, print("points sur E'(Q) = ",listpoints[1]));
  np1 = listpoints[2][1]; np2 = listpoints[2][2];
  listpoints2 = vector(#listpoints[1],i,0);
  for( i = 1, #listpoints[1],
    listpoints2[i] = [0,0];
    aux1 = listpoints[1][i][1]^2;
    if( aux1 != 0,
      aux2 = listpoints[1][i][2];
      listpoints2[i][1] = aux2^2/aux1/4;
      listpoints2[i][2] = aux2*(bpinit-aux1)/aux1/8
    , listpoints2[i] = listpoints[1][i]));
if( DEBUGLEVEL_ell >= 1, print("points sur E(Q) = ",listpoints2); print());
  pointgen = concat(pointgen,listpoints2);

  certainp = (np1 == np2);
if( DEBUGLEVEL_ell >= 1,
  if( certainp,
    print("[E'(Q):phi(E(Q))]   = ",1<<np1);
    print("#S^(phi)(E/Q)       = ",1<<np2);
    print("#III(E/Q)[phi]      = 1"); print()
  ,
    print("[E'(Q):phi(E(Q))]  >= ",1<<np1);
    print("#S^(phi)(E/Q)       = ",1<<np2);
    print("#III(E/Q)[phi]     <= ",1<<(np2-np1)); print());

  if( !certain && (np2 > np1), print1(1<<(np2-np1)," <= "));
  print1("#III(E/Q)[2]       ");
  if( certain && certainp, print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));

  print("#E(Q)[2]            = ",1<<tors);
);
  rang = n1+np1-2;
if( DEBUGLEVEL_ell >= 1,
  if( certain && certainp,
    print("#E(Q)/2E(Q)         = ",(1<<(rang+tors)));
    print("rang                = ",rang); print()
  ,
    print("#E(Q)/2E(Q)        >= ",(1<<(rang+tors))); print();
    print(rang," <= rang          <= ",n2+np2-2); print()
  ));

  strange = (n2+np2-n1-np1)%2;
  if( strange,
if( DEBUGLEVEL_ell >= 1,
      print(" !!! III doit etre un carre !!!"); print("donc"));
    if( certain,
      np1++;
      certainp = (np1 == np2);
if( DEBUGLEVEL_ell >= 1,
        if( certainp,
          print("[E'(Q):phi(E(Q))]   = ",1<<np1);
          print("#S^(phi)(E/Q)       = ",1<<np2);
          print("#III(E/Q)[phi]      = 1"); print()
        ,
          print("[E'(Q):phi(E(Q))]  >= ",1<<np1);
          print("#S^(phi)(E/Q)       = ",1<<np2);
          print("#III(E/Q)[phi]     <= ",1<<(np2-np1)); print())
      )
    ,
    if( certainp,
      n1++;
      certain = (n1 == n2);
if( DEBUGLEVEL_ell >= 1,
        if( certain,
          print("[E(Q):phi'(E'(Q))]   = ",1<<n1);
          print("#S^(phi')(E'/Q)       = ",1<<n2);
          print("#III(E'/Q)[phi']      = 1"); print()
        ,
          print("[E(Q):phi'(E'(Q))]  >= ",1<<n1);
          print("#S^(phi')(E'/Q)      = ",1<<n2);
          print("#III(E'/Q)[phi']    <= ",1<<(n2-n1)); print())
      )
    , n1++)
  );

if( DEBUGLEVEL_ell >= 1,
  print("#S^(2)(E/Q)           = ",1<<(n2+np2-1));
  if( !certain && (np2 > np1), print1(1<<(np2-np1)," <= "));
  print1("#III(E/Q)[2]       ");
  if( certain && certainp, print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));
  print("#E(Q)[2]            = ",1<<tors);
);
  rang = n1+np1-2;
if( DEBUGLEVEL_ell >= 1,
  if( certain && certainp,
    print("#E(Q)/2E(Q)         = ",(1<<(rang+tors))); print();
    print("rang                = ",rang); print()
  ,
    print("#E(Q)/2E(Q)        >= ",(1<<(rang+tors))); print();
    print(rang," <= rang          <= ",n2+np2-2); print())
  ));

\\ fin de strange

  if( ELLREDGENFLAG, pointgen = ellredgen(ell,pointgen));
  pointgen = concat(ellsort(elltorseven(ell)[3]),pointgen);
if( DEBUGLEVEL_ell >= 1, print("points = ",pointgen));
if( DEBUGLEVEL_ell >= 3, print("fin de ell2descent_viaisog"));
  return([rang,n2+np2-2+tors,pointgen]);
}


