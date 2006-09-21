\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\       Copyright (C) 2005 Denis Simon
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
\\  *          VERSION 25/10/2005               *
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
\\ definition du corps :
\\ Courbes sans 2-torsion: y^2 = x^3+A*x^2+B*x+C
\\ A,B,C entiers.
\\ -> main(A,B,C);
\\
\\ retourne un vecteur [r,s,vec]
\\ ou r est le rang probable,
\\ s est le 2-rang du groupe de Selmer
\\ vec est une liste de points independants
\\
\\ Courbes de la forme: k*y^2 = x^3+A*x^2+B*x+C
\\ sans 2-torsion, A,B,C entiers.
\\ bnf = bnfinit(x^3+A*x^2+B*x+C);
\\ ell = ellinit([0,A,0,B,C]);
\\ rank = ellrank(ell,bnf,k);
\\
\\ On peut avoir plus ou moins de details de calculs avec
\\ DEBUGLEVEL = 0;
\\ DEBUGLEVEL = 1; 2; 3;...
\\
\\
\\

{
\\
\\ Variables globales usuelles
\\

global(DEBUGLEVEL,LIM1,LIM3,LIMTRIV,MAXPROB,LIMBIGPRIME);

  DEBUGLEVEL = 1; \\ pour avoir plus ou moins de details
  LIM1 = 5;       \\ limite des points triviaux sur les quartiques
  LIM3 = 50;      \\ limite des points sur les quartiques ELS
  LIMTRIV = 10;   \\ limite des points triviaux sur la courbe elliptique

\\
\\  Variables globales techniques
\\

  MAXPROB = 20;
  LIMBIGPRIME = 30; \\ pour distinguer un petit 1nombre premier d'un grand
                    \\ utilise un test probabiliste pour les grands
                    \\ si LIMBIGPRIME = 0, n'utilise aucun test probabiliste
}

\\
\\  Programmes
\\

if(DEBUGLEVEL >= 4, print("mysubst "));
{
mysubst(polsu,subsx) =
  if ( type(lift(polsu)) == "t_POL" ,
    return ( simplify(subst(lift(polsu),variable(lift(polsu)),subsx)) )
  , return( simplify(lift(polsu))));
}
if(DEBUGLEVEL >= 4, print("degre "));
{
degre(idegre)=
local(ideg,jdeg);
  ideg = idegre; jdeg = 0;
  while (ideg>>=1, jdeg++);
  return (jdeg);
}
if(DEBUGLEVEL >= 4, print("ratpoint "));
{
ratpoint(pol,lim,flag=0)=
local(listpoints,point1,odd,deg4,pol16,tab16,pol9,tab9,pol0,vecz,vecx,lead,zz,xx,x16,x9,step,evpol);
\\ Recherche de points sur y^2=pol.
\\ Les coeff de pol sont entiers.
\\ Si flag ==0, cherche un seul point, sinon plusieurs.
if (DEBUGLEVEL >= 4, print("entree dans ratpoint avec pol=",pol);print("lim = ",lim););
\\
  if (flag, listpoints = []);
\\
  point1 = [];
\\          cas triviaux
  if ( issquare(polcoeff(pol,0)),
    point1 = [ 0 , sqrtrat(polcoeff(pol,0)) ];
if (DEBUGLEVEL >= 3, print("solution triviale: e est un carre"));
    if ( !flag ,
if (DEBUGLEVEL >= 4, print("fin de ratpoint"));
      return (point1));
    listpoints = concat(listpoints,[point1]));
  if ( issquare(pollead(pol)),
    point1 = [ 1 , sqrtrat(pollead(pol)) , 0];
if (DEBUGLEVEL >= 3, print("solution triviale: a est un carre"));
    if ( !flag ,
if (DEBUGLEVEL >= 4, print("fin de ratpoint"));
      return (point1));
    listpoints = concat(listpoints,[point1]));
  odd = poldegree(pol)%2;
  deg4 = poldegree(pol) == 4;

\\ initialisation du crible modulo 16  et 9
  if ( deg4 ,
    pol16 = (Vec(pol)*Mod(1,16))~;
    tab16 = matrix(16,16);
    for(xx = 0, 16-1,
      for(zz = 0, 16-1,
        tab16[xx+1,zz+1] = issquare([xx^4,xx^3*zz,xx^2*zz^2,xx*zz^3,zz^4]*pol16)));
    pol9 = (Vec(pol)~)*Mod(1,9);
    tab9 = matrix(9,9);
    for(xx = 0, 9-1,
      for(zz = 0, 9-1,
        tab9[xx+1,zz+1] = issquare([xx^4,xx^3*zz,xx^2*zz^2,xx*zz^3,zz^4]*pol9)));
  );

pol0=polcoeff(pol,0);

  if ( odd,
    vecz=vector(lim,i,i);
    vecx=vector(lim,i,i);
  ,
\\ si le degre de pol est pair, il faut que le coeff constant soit
\\ un carre mod xx. Idem pour le dominant mod zz.
    vecz = vector(lim);
    vecx = vector(lim);
    lead = pollead(pol);
    zz = 0;xx = 0;
    for( i = 1, lim,
      xx++; while( !issquare(Mod(pol0,xx)),xx++);vecx[i] = xx;
      zz++; while( !issquare(Mod(lead,zz)),zz++);vecz[i] = zz;
  ));

if( DEBUGLEVEL >= 5, print("vecx = ",vecx));
if( DEBUGLEVEL >= 5, print("vecz = ",vecz));

\\ boucle sur x=xx/zz
  for(eps = 0,1,
    if(eps == 1,vecx = -vecx);
    for(ix = 1,lim, xx = vecx[ix];
      if (deg4 , x16 = xx%16+1; x9 = xx%9+1;);
      if(odd || xx%2 , step = 0, step = 1 );
      for( iz = 1, lim , zz = vecz[iz];
        if(gcd(zz,xx)>1,next);
        if(  !deg4 ||
             (tab16[x16,zz%16+1] && tab9[x9,zz%9+1])
        ,
          evpol = subst(pol,x,xx/zz);
          if ( issquare(evpol) ,
            point1=[xx/zz,sqrtrat(evpol)];
            if ( !flag , break(2));
            if ( !flag ,if (DEBUGLEVEL >= 4, print("fin de ratpoint"));return (point1));
            listpoints = concat(listpoints,[point1])));
  )));
if (DEBUGLEVEL >= 3, print("point sur la quartique=",point1));
if (DEBUGLEVEL >= 4, print("sortie de ratpoint"));
  if (!flag ,return (point1),return (listpoints));
}
if(DEBUGLEVEL >= 4, print("nfissquare"));
{
nfissquare( nf, a)=
local(ta,pfact,alift,r1,py);
if (DEBUGLEVEL >= 5, print("entree dans nfissquare ",a));
\\ si a est un carre, renvoit 1, sinon 0.
  if ( a==0 || a==1,
if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
    return (1));

  alift=lift(a);
  if(!poldegree(alift),alift=polcoeff(alift,0));
  ta = type(alift);

  if ( ta != "t_POL" ,
if (DEBUGLEVEL >= 5, print("fin de nfissquare "));
    return (issquare(alift)));

\\
\\ tous les plongements reels doivent etre >0
\\
  r1=nf.sign[1];
  for ( i = 1 ,r1,
    py=mysubst(alift,nf.roots[i]);
    if ( sign(py) < 0 ,
if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
      return ([])));

\\ factorisation sur nf du polynome y^2-a :

   pfact=factor(polresultant(y^2-alift,nf.pol))[1,1];
if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
  return (poldegree(pfact) == 3);
}
if(DEBUGLEVEL >= 5, print("psquare "));
{
psquare( a, p)=
local(ap,v);
if (DEBUGLEVEL >= 5, print("entree dans psquare",[a,p]));
\\ a est un entier
\\ renvoie 1 si a est un carre dans Zp 0 sinon
  if ( a == 0, if (DEBUGLEVEL >= 5, print("fin de psquare 1"));return (1));
\\
  v = valuation(a,p);
  if ( v%2 ,
if (DEBUGLEVEL >= 5, print("fin de psquare 0"));
    return (0));
  if ( p == 2 ,
    ap = (a>>v)%8-1,
    ap = kronecker(a/p^v,p)-1
  );
if (DEBUGLEVEL >= 5, print("fin de psquare ", !ap));
  return (!ap);
}
if(DEBUGLEVEL >= 4, print("lemma6 "));
{
lemma6(pol, p, nu, xx)=
local(gx,gpx,lambda,mu);
\\ pour les p <>2
  gx = subst( pol , x , xx);
  if ( psquare(gx,p) , return (1));
  gpx = subst( pol' , x , xx);
  lambda = valuation(gx,p);  mu = valuation(gpx,p);

  if ( lambda> 2*mu , return (1));
\\  if ( (lambda>=mu+nu) && (nu>mu) , return(1));
  if ( (lambda >= 2*nu)  && (mu >= nu) , return (0));
  return (-1);
}
if(DEBUGLEVEL >= 4, print("lemma7 "));
{
lemma7( pol, nu, xx)=
local(gx,gpx,lambda,mu,q);
\\ pour p=2
  gx = subst( pol , x , xx);
  if ( psquare(gx,2) , return (1));
  gpx = subst( pol' , x , xx);
  lambda = valuation(gx,2);mu = valuation(gpx,2);
  if ( lambda>2*mu , return (1));
  if (nu > mu ,
    if ( lambda%2 , return (-1));
    q = mu+nu-lambda;
    if ( q == 1 , return (1));
    if ( q == 2 && (gx>>lambda)%4 == 1 , return(1));
    return (-1));
  q = lambda-2*nu;
  if ( q >= 0 , return (0));
  if ( q == -2 && (gx>>lambda)%4 == 1 , return (0));
  return (-1);
}
if(DEBUGLEVEL >= 4, print("zpsoluble "));
{
zpsoluble(pol, p, nu, pnu, x0, pnup)=
local(result,pol2,fact,x1);
if (DEBUGLEVEL >= 5, print("entree dans zpsoluble",[pol,p,x0,nu]));
  if ( p == 2 ,
    result = lemma7(pol,nu,x0),
    result = lemma6(pol,p,nu,x0));
  if ( result == +1 ,
if (DEBUGLEVEL >= 5, print("fin de zpsoluble 1 lemma"));
    return (1));
  if ( result == -1 ,
if (DEBUGLEVEL >= 5, print("fin de zpsoluble 0 lemma"));
    return (0));
  pnup=pnu*p;
  nu++;
  if ( p< LIMBIGPRIME || !LIMBIGPRIME,
    for ( i = 0 ,p-1,
      if ( zpsoluble(pol,p,nu,pnup,x0+pnu*i) ,
if (DEBUGLEVEL >= 5, print("fin de zpsoluble"));
        return (1)));
  ,
    pol2 = subst(pol,x,x0+pnu*x);
    pol2 /= content(pol2);
    pol2 = pol2*Mod(1,p);
    if ( !poldegree(pol2), return(0));
    fact = factormod(pol2,p)[,1];
    for ( i = 1 , length(fact),
      x1 = -centerlift(polcoeff(fact[i],0));
      if ( zpsoluble(pol,p,nu,pnup,x0+pnu*x1) ,
if (DEBUGLEVEL >= 5, print("fin de zpsoluble"));
        return (1)));
    for ( i = 1 ,MAXPROB,
      x1 = random(p);
      if ( zpsoluble(pol,p,nu,pnup,x0+pnu*x1) ,
if (DEBUGLEVEL >= 5, print("fin de zpsoluble"));
        return (1)));
  );
if (DEBUGLEVEL >= 2, print("******* test probabiliste en p = ",p,"*******"));
if (DEBUGLEVEL >= 5, print("fin de zpsoluble"));
  return (0);
}
if(DEBUGLEVEL >= 4, print("qpsoluble "));
{
qpsoluble(pol, p) =
if (DEBUGLEVEL >= 5, print("entree dans qpsoluble ",p);print("pol = ",pol));
  if ( psquare(pollead(pol),p) ,
    if (DEBUGLEVEL >= 5, print("fin de qpsoluble 1"));
    return (1));
  if ( psquare(polcoeff(pol,0),p) ,
    if (DEBUGLEVEL >= 5, print("fin de qpsoluble 1"));
    return (1));
  if ( zpsoluble(pol,p,0,1,0) ,
    if (DEBUGLEVEL >= 5, print("fin de qpsoluble 1"));
    return (1));
  if ( zpsoluble(polrecip(pol),p,1,p,0) ,
    if (DEBUGLEVEL >= 5, print("fin de qpsoluble 1"));
    return (1));
if (DEBUGLEVEL >= 5, print("fin de qpsoluble 0"));
  return (0);
}
if(DEBUGLEVEL >= 4, print("locallysoluble "));
{
locallysoluble(pol)=
local(plist,pol0,pro,p);
\\ teste l'existence locale de solutions de y^2 = pol(x,z)
if (DEBUGLEVEL >= 4, print("entree dans locallysoluble :",pol));
pol0=pol;

\\ place reelle
  if ( !(poldegree(pol)%2) &&  sign(pollead(pol)) < 0
         && sign(polcoeff(pol,0)) < 0 && polsturm(pol)==0 ,
    if (DEBUGLEVEL >= 2, print(" not ELS at infty "));
    if (DEBUGLEVEL >= 4, print("fin de locallysoluble"));
    return (0));

\\
\\ places finies de plist */
\\
  pol *= denominator(content(pol))^2;

  pro = poldisc(pol0);
  plist = factor (abs(2*pro))[,1] ;
if(DEBUGLEVEL >= 4, print("liste de premiers = ",plist));
  for ( i = 1 ,length(plist),
    p =  plist[i];
if (DEBUGLEVEL >= 4, print("p=",p));
    if ( !qpsoluble(pol,p) ,
if (DEBUGLEVEL >= 2, print(" not ELS at ",p));
if (DEBUGLEVEL >= 4, print("fin de locallysoluble"));
      return (0)));

  if (DEBUGLEVEL >= 2, print(" quartic ELS "));
  if (DEBUGLEVEL >= 4, print("fin de locallysoluble"));
  return (1);
}
if(DEBUGLEVEL >= 4, print("sqrtrat "));
{
sqrtrat(a) =
if (DEBUGLEVEL >= 5, print("entree dans sqrtrat ",a));
  return(sqrtint(numerator(a))/sqrtint(denominator(a)));
}
if(DEBUGLEVEL >= 4, print("redquartique "));
{
redquartique(pol,deg=4)=
local(M,t);
\\ reduction d'une quartique.
\\ ou plus generalement d'un polynome de degre deg.
if (DEBUGLEVEL >= 4, print("entree dans redquartique"));
if (DEBUGLEVEL >= 3, print(" reduction de la quartique ",pol));

  M=[1,0;0,1];t=1;
  if(abs(polcoeff(pol,0))<abs(polcoeff(pol,deg)),
    pol=polrecip(pol);
    M=[0,1;1,0]);
  while( abs(polcoeff(pol,0))>=abs(polcoeff(pol,deg)) && t,
    t=round(polcoeff(pol,deg-1)/deg/polcoeff(pol,deg));
    pol=polrecip(subst(pol,x,x-t));
if (DEBUGLEVEL >= 4, print("maxpol="vecmax(abs(Vec(pol)))));
    M=M*[1,-t;0,1]*[0,1;1,0];
  );

if (DEBUGLEVEL >= 3, print(" quartique reduite = ",pol));
if (DEBUGLEVEL >= 4, print("sortie de redquartique"));

  return([pol,M]);
}
if(DEBUGLEVEL >= 4, print("minquartique "));
{
minquartique( pol, r,a,b)=
local(gcc,rp,pol2);
\\ minimisation d'une quartique issue de la 2-descente
if (DEBUGLEVEL >= 4, print("entree dans minquartique"));
if (DEBUGLEVEL >= 4, print([r,a,b]));
if (DEBUGLEVEL >= 3, print(" minimisation de la quartique ",pol));
  if (a==0,
    rp=0
  ,
    gcc=gcd(b,a);
    while(gcc!=1 ,
      b/=gcc;
if(DEBUGLEVEL >= 5, print("********** gcc = ",gcc));
      gcc=gcd(b,gcc));
    rp=(r/a)%b;
    if( rp > (b>>1) , rp-=b);
  );
 pol2=subst(pol,x,rp+b*x)/b^4;
 if (DEBUGLEVEL >= 3, print(" quartique minime = ",pol2));
 if (DEBUGLEVEL >= 4, print("fin de minquartique"));
 return([pol2,[b,rp;0,1]])
}
if(DEBUGLEVEL >= 4, print("ellrank "));
{
ellrank(ell,ext,K=1,help=[],redflag=1) =
local(A,B,C,discP,S,SLprod,fact,trouve,SLlist,oddclass,LS2gen,polrel,alpha,ttheta,KS2gen,LS2genunit,normLS2gen,normcoord,LS2coordtilda,LS2,i,aux,j,listgen,listpoints,listpointstriv,listpointsmwr,list,m1,m2,lastloc,maskwhile,iwhile,zc,iaux,liftzc,den,ispointtriv,point,c,b,a,sol,found,r2,normzc,r,q2,mattr,vec,z1,ff,cont,d,e,polorig,pol,redq,transl,UVW,pointxx,point2,v,rang);
if (DEBUGLEVEL >= 4, print("entree dans ellrank"));
\\ ell est donne par ellinit.
\\ si ell= K*y^2=P(x), alors ext est le buchinitfu de l'extension.
\\ theta est une racine de P.
\\ dans la suite ext est note L = Q(theta).
\\ help est une liste de points deja connus sur ell.
\\ ell est de la forme K*y^2=x^3+A*x^2+B*x+C */
\\ ie ell=[0,A,0,B,C], avec A,B et C entiers */
\\
\\ si redflag != 0, on utilise la reduction modulo les carres.
\\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      construction de L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
  A = ell.a2;if(DEBUGLEVEL >= 2, print("A= ",A));
  B = ell.a4;if(DEBUGLEVEL >= 2, print("B= ",B));
  C = ell.a6;if(DEBUGLEVEL >= 2, print("C= ",C));
  discP = ell.disc/16;
  S = -abs(K*discP);
  SLprod = K*(3*x^2+2*A*x+B);
  SLlist = idealfactor(ext,SLprod)[,1]~;
  aux = []; SLprod = 1;
   for (i = 1 , length(SLlist),
     if ( !(K%SLlist[i][1]) || valuation(discP,SLlist[i][1]) > 1,
       SLprod = idealmul(ext,SLprod,SLlist[i]);
       aux = concat(aux,[SLlist[i]])));
  SLlist = aux;
  oddclass = 0;
  while (!oddclass,
\\ Constructoin de S:
if (DEBUGLEVEL >= 4, print("SLlist=",SLlist));
\\ Construction des S-unites
    LS2gen = bnfsunit(ext,SLlist);
if (DEBUGLEVEL >= 4, print("LS2gen=",LS2gen));
\\ on ajoute la partie paire du groupe de classes.
    oddclass = LS2gen[5][1]%2;
    if (!oddclass,
if (DEBUGLEVEL >= 4,print(LS2gen[5]));
      S *= LS2gen[5][3][1][1,1];
      SLprod = idealmul(ext,SLprod,LS2gen[5][3][1]);
      fact = idealfactor(ext,LS2gen[5][3][1])[,1];
      trouve = 0; i = 0;
      while (!trouve ,
        i++; trouve = 1;
        for (j = 1, length(SLlist),
          if (SLlist[j] ==  fact[i] , trouve = 0;break)));
      SLlist = concat(SLlist,[fact[i]]))
  );

if (DEBUGLEVEL >= 4, print("S = ",S));

  polrel = x^3+A*x^2+B*x+C;
  ttheta = Mod(x,polrel);            \\ ttheta est la racine de P(x)

  KS2gen = factor(S)[,1]~;

if (DEBUGLEVEL >= 3, print("#KS2gen = ",length(KS2gen)));
if (DEBUGLEVEL >= 3, print("KS2gen = ",KS2gen));

  LS2genunit = ext.tufu;
  LS2genunit = concat(lift(nfbasistoalg(ext,LS2gen[1])),LS2genunit);

  LS2genunit = subst(LS2genunit,x,ttheta);
  LS2genunit = LS2genunit*Mod(1,polrel);
if (DEBUGLEVEL >= 3, print("#LS2genunit = ",length(LS2genunit)));
if (DEBUGLEVEL >= 3, print("LS2genunit = ",LS2genunit));

\\ dans LS2gen, on ne garde que ceux dont la norme
\\ est un carre.

  normLS2gen = norm(LS2genunit);
if (DEBUGLEVEL >= 4, print("normLS2gen=",normLS2gen));

\\ matrice de l'application norme

  normcoord = matrix(length(KS2gen),length(normLS2gen));
  for( j = 1 ,length(normLS2gen),
    normcoord[1,j] = (sign(normLS2gen[j]) < 0);
    for( i = 2 ,length(KS2gen),
      normcoord[i,j] = valuation(normLS2gen[j],KS2gen[i])));
if (DEBUGLEVEL >= 4, print("normcoord=",normcoord));

\\ construction du noyau de la norme

  LS2coordtilda = lift(matker(normcoord*Mod(1,2)));
if(DEBUGLEVEL >= 4, print("LS2coordtilda=",LS2coordtilda));
  LS2 = vector(length(LS2coordtilda[1,]),i,0);
  for( i = 1 ,length(LS2coordtilda[1,]),
    aux = 1;
    for ( j = 1 ,length(LS2coordtilda[,i]),
      if ( sign(LS2coordtilda[j,i]),
        aux *= LS2genunit[j]));
    LS2[i] = aux;
  );
if (DEBUGLEVEL >= 3, print("LS2 = ",LS2));
if (DEBUGLEVEL >= 3, print("norm(LS2)= ",norm(LS2)));

\\ Reduction des generateurs de LS2

  if(redflag,
    for( i = 1,length(LS2),
      LS2[i]=reducemodsquares(LS2[i],2)));

\\ Fin de la construction de LS2

  listgen = LS2;
  if (DEBUGLEVEL >= 2, print("LS2gen = ",listgen));
  if (DEBUGLEVEL >= 2, print("#LS2gen = ",length(listgen)));
  listpoints = [];
\\
  if(DEBUGLEVEL >= 3, print("(A,B,C) = ",[A,B,C]));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Recherche de points triviaux.  \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if(DEBUGLEVEL >= 2, print(" recherche de points triviaux sur la courbe "));
  listpointstriv = ratpoint(K*polrel,LIMTRIV,1);
  for(i=1,length(listpointstriv),listpointstriv[i][2]/=K);
  listpointstriv = concat(help,listpointstriv);
if(DEBUGLEVEL >= 1, print("points triviaux sur la courbe = ",listpointstriv));
\\

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          parcours de L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  listpointsmwr = [];
  list = [ 6 , ell.disc, 0 ];
  m1 = 0;m2 = 0;lastloc = -1;
  maskwhile = 1<<length(listgen);
  iwhile = 1;
  while ( iwhile < maskwhile ,
    if ( DEBUGLEVEL >= 4, print("iwhile = ",iwhile);print("listgen = ",listgen));
    zc = Mod(1,polrel); j = 1; iaux = iwhile;
    while ( iaux ,
      if( iaux%2, zc *= listgen[j] );
      iaux >>= 1; j++ );
if(DEBUGLEVEL >= 2, print();print("zc=",zc));
    if(redflag,
      zc=reducemodsquares(zc,2);
    liftzc = lift(zc);
    den = denominator(content(liftzc))^2;
    zc *= den; liftzc *= den;
if(DEBUGLEVEL >= 2, print("zcred=",zc)));


\\ Est-ce un point trivial ?
    ispointtriv = 0;
    for ( i = 1,length(listpointstriv),
      point = listpointstriv[i];
      if ( length(point)==2 || point[3] != 0,
        if ( nfissquare(ext.nf,K*(point[1]-x)*liftzc),
          if(DEBUGLEVEL >= 2, print(" vient du point trivial ",point));
          listpointsmwr = concat(listpointsmwr,[point]);
          m1++;
          if ( degre(iwhile) > lastloc ,m2++);
          found = (ispointtriv = 1); break)));

\\ Ce n'est pas un point trivial
    if ( !ispointtriv ,
      c = polcoeff(liftzc,2);
      b = -polcoeff(liftzc,1);
      a = polcoeff(liftzc,0);
      r2 = sqrtint(norm(zc));
      sol = 0;
\\ \\\\\\\\\\\\\
\\ On cherche a ecrire zc sous la forme a-b*theta
\\ \\\\\\\\\\\\\
      if ( c == 0 || ispointtriv, sol = 1);
\\ Il faut resoudre une forme quadratique
\\ Existence (locale = globale) d'une solution :
      alpha = A*b*c+B*c^2-a*c+b^2;
      if (!alpha,
        zc = norm(zc)/zc;
if(DEBUGLEVEL >= 2, print("zc = ",zc));
        sol = 1);
      if (!sol ,
        q2 = [c,-(b+A*c),(A^2*c+A*b-B*c+a);
          -(b+A*c),A^2*c+A*b-B*c+a,-A^3*c-A^2*b+2*A*B*c-A*a+B*b-C*c;
          (A^2*c+A*b-B*c+a),-A^3*c-A^2*b+2*A*B*c-A*a+B*b-C*c,
          A^4*c+A^3*b-3*A^2*B*c+A^2*a-2*A*B*b+2*A*C*c+B^2*c-B*a+C*b];
        sol=Qfsolve(q2);
        if(type(sol) == "t_INT",
if(DEBUGLEVEL >= 2, print("hilbert(",alpha,",",c,") =  -1"));
          iwhile++;next
        );
        z1 = (sol[3]*ttheta+sol[2])*ttheta+sol[1];
        zc *= z1^2; sol = 1;
if(DEBUGLEVEL >= 2, print("zc = ",zc));
      );
\\ \\\\\\\\\\
\\ Maintenant zc est de la forme a-b*theta
\\ \\\\\\\\\\
      if (sol && !ispointtriv,
        liftzc = lift(zc);
        normzc = norm(zc);
if(DEBUGLEVEL >= 3, print(" zc = ",liftzc));
if(DEBUGLEVEL >= 4, print(" N(zc) = ",normzc));
        if ( poldegree(liftzc) < 2 ,, print("****** ERROR  c <> 0 *******");1/0);
        b = -polcoeff(liftzc,1); \\ rem: b n'est jamais nul
        a = polcoeff(liftzc,0);
        ff = factor(gcd(a,b));
if(DEBUGLEVEL >= 4, print(" ff=",ff));
        for(i=1,length(ff[,1]),
          cont=ff[i,1]^(2*(ff[i,2]\2));
          a/=cont;b/=cont;zc/=cont;normzc/=cont^3);
if(DEBUGLEVEL >= 4, print(" [a,b]=",[a,b]));
        if(issquare(K*b),
if (DEBUGLEVEL >= 3, print("b est un carre"));
          point=[a/b,sqrtrat(subst(polrel,x,a/b)/K)];
if (DEBUGLEVEL >= 2, print("point trouve = ",point));
          listpointsmwr=concat(listpointsmwr,[point]);
          m1++;
          if ( degre(iwhile) > lastloc ,m2++);
          ispointtriv=1
        );
      );
\\ \\\\\\\\\\\
\\ Construction de la quartique
\\ \\\\\\\\\\\
      if (sol && !ispointtriv,
        r = sqrtint(normzc); \\ r est entier car a et b le sont.
if(DEBUGLEVEL >= 4, print(" r=",r));
        c = -2*(A*b+3*a);
if(DEBUGLEVEL >= 4, print(" c=",c));
        d = 8*r;
if(DEBUGLEVEL >= 4, print(" d=",d));
        e = (A^2*b^2 - 2*A*a*b-4*B*b^2-3*a^2);
if(DEBUGLEVEL >= 4, print(" e=",e));
        polorig = b*(x^4+c*x^2+d*x+e);
if(DEBUGLEVEL >= 2, print(" quartique : ",b*K,"*Y^2 = ",lift(polorig/b)));
        list[3] = b;
        pol = polorig;
\\ minimisation du discriminant de la quartique
        redq = minquartique(pol,r,a,b);
if(DEBUGLEVEL >= 2, print(" minime: ",K,"*Y^2 = ",lift(redq[1])));
        pol = redq[1]; transl = redq[2];
\\ reduction des coefficients de la quartique
        redq = redquartique(pol,4);
if(DEBUGLEVEL >= 2, print(" reduite: ",K,"*Y^2 = ",lift(redq[1])));
        pol = redq[1];
        transl = transl*redq[2];
        point = ratpoint(K*pol,LIM1);
        if ( point != 0,
          m1++;
          if( length(point) == 2, point = concat(point,[1]));
if(DEBUGLEVEL >= 2, print("point=",point));
          point[1] = (point[1]*transl[1,1]+transl[1,2]*point[3])/(point[1]*transl[2,1]+transl[2,2]*point[3]);
          mattr = matid(3);
          mattr[1,1] = -2*b^2; mattr[1,2] = (A*b+a)*b;
          mattr[1,3] = a^2+(2*B-A^2)*b^2; mattr[2,2] = -b;
          mattr[2,3] = a+A*b; mattr[3,3] =r;
          UVW = [point[1]^2,1,point[1]]~;
          vec = (mattr^(-1))*UVW;
          z1 = (vec[3]*ttheta+vec[2])*ttheta+vec[1];
          zc *= z1^2;if(DEBUGLEVEL >= 3, print("zc*z1^2 = ",zc));
          zc /= -polcoeff(lift(zc),1);
if(DEBUGLEVEL >= 3, print("zc*z1^2 = ",zc));
          pointxx = polcoeff(lift(zc),0);
          point2 = [ pointxx , sqrtrat(subst(x^3+A*x^2+B*x+C,x,pointxx)/K)];
if(DEBUGLEVEL >= 1, print(" point trouve = ",point2));
          listpointsmwr = concat(listpointsmwr,[point2]);
          if ( degre(iwhile) > lastloc ,m2++);
          found = 1;
        ,
          if ( locallysoluble(K*pol),
            if ( degre(iwhile) > lastloc ,m2++;lastloc=degre(iwhile));
            point = ratpoint(K*pol,LIM3);
            if ( point!=0,
              m1++;
              if( length(point) == 2, point = concat(point,[1]));
              point[1] =    (point[1]*transl[1,1]+transl[1,2]*point[3])/(point[1]*transl[2,1]+transl[2,2]*point[3]);
              point[2] = sqrtrat(subst(polorig,x,point[1]/point[3])/K);
              mattr = matid(3);
              mattr[1,1] = -2*b^2; mattr[1,2] = (A*b+a)*b;
              mattr[1,3] = a^2+(2*B-A^2)*b^2; mattr[2,2] = -b;
              mattr[2,3] = a+A*b; mattr[3,3] =r;
              UVW = [point[1]^2,point[3]^2,point[1]*point[3]]~;
              vec = (mattr^(-1))*UVW;
              z1 = (vec[3]*ttheta+vec[2])*ttheta+vec[1];
              zc *= z1^2;
              zc /= -polcoeff(lift(zc),1);
if(DEBUGLEVEL >= 3, print(" zc*z1^2 = ",zc));
              pointxx = polcoeff(lift(zc),0);
              point2 = [ pointxx , sqrtrat(subst(x^3+A*x^2+B*x+C,x,pointxx)/K)];
if(DEBUGLEVEL >= 1, print(" point trouve = ",point2));
              listpointsmwr=concat(listpointsmwr,[point2]);
              found = 1;
            )
          )
        )
      )
    );
    if ( found || ispointtriv,
      found = 0; lastloc = -1;
      v = 0; iaux = (iwhile>>1);
      while( iaux , iaux>>=1;v++);
      maskwhile >>= 1;
      listgen = vecextract(listgen,1<<length(listgen)-1<<v-1);
      iwhile = (1<<v)
    ,iwhile ++
    )
  );

if(DEBUGLEVEL >= 2,
  print("");
  print("rang des points trouves = ",m1);
  print("rang du groupe de Selmer = ",m2));
if(DEBUGLEVEL >= 1,
  print("#S(E/Q)[2]    = ",1<<m2));
  if (m1==m2,
if(DEBUGLEVEL >= 1,
    print("#E(Q)/2E(Q)   = ",1<<m1);
    print("#III(E/Q)[2]  = 1");
    print("rang(E/Q)     = ",m1));
    rang = m1;
  ,
if(DEBUGLEVEL >= 1,
    print("#E(Q)/2E(Q)  >= ",1<<m1);
    print("#III(E/Q)[2] <= ",1<<(m2-m1));
    print("rang(E/Q)    >= ",m1));
    rang = m1;
    if ((m2-m1)%2,
if(DEBUGLEVEL >= 1,
      print(" III devrait etre un carre, donc ");
      if(m2-m1>1,
        print("#E(Q)/2E(Q)  >= ",1<<(m1+1));
        print("#III(E/Q)[2] <= ",1<<(m2-m1-1));
        print("rang(E/Q)    >= ",m1+1);
      ,
        print("#E(Q)/2E(Q)  = ",1<<(m1+1));
        print("#III(E/Q)[2] = 1");
        print("rang(E/Q)    = ",m1+1)));
      rang = m1+1);
  );
if(DEBUGLEVEL >= 1, print("listpointsmwr = ",listpointsmwr));
  for ( i = 1 ,length(listpointsmwr) ,
    if ( subst(polrel,x,listpointsmwr[i][1])-K*listpointsmwr[i][2]^2,
      print("****** MAUVAIS POINT = ",listpointsmwr[i])));
if (DEBUGLEVEL >= 4, print("fin de ellrank"));
  return ([rang,m2,listpointsmwr]);
}
if(DEBUGLEVEL >= 4, print("main"));
{
main(A,B,C,help=[])=
local(eqtheta,bnf,ell,rang,time1);
if (DEBUGLEVEL >= 3, print("entree dans main"));
  eqtheta = x^3+A*x^2+B*x+C;
if (DEBUGLEVEL >= 1, print("courbe elliptique : Y^2 = ",eqtheta));
if (DEBUGLEVEL >= 3, print1("bnfinit "));
  gettime();
  bnf = bnfinit(eqtheta,1);
  time1 = gettime();
if (DEBUGLEVEL >= 3, print("done"));
  ell = ellinit([0,A,0,B,C]);
  rang = ellrank(ell,bnf,1,help);
if (DEBUGLEVEL >= 3, print("fin de main"));
if (DEBUGLEVEL >= 2, print("temps dans bnfinit = ",time1));
if (DEBUGLEVEL >= 2, print("temps pour le reste = ",gettime()));
  return (rang);
}
if( DEBUGLEVEL >= 4, print("reducemodsquares "));
{reducemodsquares(delta,d)=
\\ reduction du coefficient de x^d dans ( delta modulo les carres )
local(deg,xx,z,qd,Qd,reduc);

  deg=poldegree(delta.mod);
  xx=Mod(x,delta.mod);
  z=subst(Pol(vector(deg,i,eval(Str("a"i)))),x,xx);
  qd=polcoeff(lift(delta*z^2),d,x);
  Qd=simplify(matrix(deg,deg,i,j,deriv(deriv(qd,eval(Str("a"i))),eval(Str("a"j)))/2));

  reduc=IndefiniteLLL(Qd);
  if(length(reduc)==2, reduc=reduc[2][,1] );

  return(delta*subst(Pol(reduc),x,xx)^2);
}
