\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\       Copyright (C) 2007 Denis Simon
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
\\
\\  *********************************************
\\  *          VERSION 16/03/2007               *
\\  *********************************************
\\
\\ Programme de calcul du rang des courbes elliptiques
\\ dans les corps de nombres.
\\ langage: GP
\\ pour l'utiliser, lancer gp, puis taper
\\ \r ell.gp
\\
\\
\\ Explications succintes :
\\ definition du corps :
\\ bnf=bnfinit(y^2+1);
\\ (il est indispensable que la variable soit y).
\\ on peut ensuite poser :
\\ X = Mod(y,bnf.pol);
\\
\\ Courbes avec 2-torsion: y^2 = x^3+a*x^2+b*x
\\ ( E doit toujours etre sous cette forme )
\\ a, b entiers.
\\ -> main2(bnf,[0,a,0,b,0]);
\\
\\ Courbes sans 2-torsion: y^2 = x^3+A*x^2+B*x+C
\\ A,B,C entiers.
\\ -> main(bnf,A,B,C);
\\ Attention A,B,C ne doivent pas etre dans un sous-corps de bnf.
\\ au besoin, on peut utiliser une translation du type x -> x+X
\\
\\ Courbes avec E[2]=4: y^2 = (x-e1)*(x-e2)*(x-e3)
\\ ->complete(bnf,e1,e2,e3);
\\ attention: les ei doivent etre entiers algebriques.
\\
\\
\\
\\ On peut avoir plus ou moins de details de calculs avec
\\ DEBUGLEVEL = 0;
\\ DEBUGLEVEL = 1; 2; 3;...
\\

{
\\
\\ Variables globales usuelles
\\

  DEBUGLEVEL = 2;   \\ pour avoir plus ou moins de details
  LIM1 = 2;    \\ limite des points triviaux sur les quartiques
  LIM3 = 4;    \\ limite des points sur les quartiques ELS
  LIMTRIV = 2; \\ limite des points triviaux sur la courbe elliptique

\\
\\  Variables globales techniques
\\

  BIGINT = 32000;   \\ l'infini
  MAXPROB = 20;
  LIMBIGPRIME = 100; \\ pour distinguer un petit 1er d'un grand
  NBIDEAUX = 10;

}

\\
\\  Programmes
\\

if (DEBUGLEVEL >= 4, print("ellchangecurveinverse "));
{
ellchangecurveinverse(ell,v)=
ellchangecurve(ell,[1/v[1],-v[1]/v[2]^2,-v[3]/v[1],(v[2]*v[3]-v[4])/v[1]^3]);
}
if (DEBUGLEVEL >= 4, print("ellchangepointinverse "));
{
ellchangepointinverse(pt,v)=
ellchangepoint(pt,[1/v[1],-v[2]/v[1]^2,-v[3]/v[1],(v[2]*v[3]-v[4])/v[1]^3]);
}
if (DEBUGLEVEL >= 4, print("mysubst "));
{
mysubst(polsu,subsx) =
  if ( type(lift(polsu)) == "t_POL" ,
    return ( simplify(subst(lift(polsu),variable(lift(polsu)),subsx)) )
  , return( simplify(lift(polsu))));
}
if (DEBUGLEVEL >= 4, print("nfmodid2 "));
{
nfmodid2(nf,a,ideal) =
if (DEBUGLEVEL >= 5, print("entree dans nfmodid2"));
\\ ideal doit etre sous la forme primedec
  if ( length(nf.zk) == 1,
if (DEBUGLEVEL >= 5, print("fin de nfmodid2"));
    return (a*Mod(1,ideal[1])));
  a = mynfeltmod(nf,a,nfbasistoalg(nf,ideal[2]));
  if (gcd(denominator(content(lift(a))),ideal[1]) == 1,
if (DEBUGLEVEL >= 5, print(" fin de nfmodid2"));
    return (a*Mod(1,ideal[1])));
if (DEBUGLEVEL >= 5, print("fin de nfmodid2"));
  return (a);
}
if (DEBUGLEVEL >= 4, print("hilb2nf "));
{
hilb2nf(nf,a,b,p)=
local(res);
if (DEBUGLEVEL >= 5, print("entree dans hilb2nf"));
  if (qpsolublenf(nf,a*x^2+b,initp(nf,p)),res=1,res=-1);
if (DEBUGLEVEL >= 5, print("fin de hilb2nf"));
  return (res);
}
if (DEBUGLEVEL >= 4, print("mynfhilbertp "));
{
mynfhilbertp(nf,a,b,p)=
\\ calcule le symbole de Hilbert quadratique local (a,b)_p
\\ * en l'ideal premier p du corps nf,
\\ * a et b sont des elements non nuls de nf, sous la forme
\\ * de polymods ou de polynomes, et p renvoye par primedec.
local(alpha,beta,sig,aux,aux2,rep);

if (DEBUGLEVEL >= 5, print("entree dans mynfhilbertp",p));
  if ( a==0 || b==0 , print("0 argument in mynfhilbertp"));
  if ( p[1]==2,
if (DEBUGLEVEL >= 5, print("fin de mynfhilbertp"));
    return (hilb2nf(nf,a,b,p)));
  if ( type(a)!="t_POLMOD",a=Mod(a,nf.pol));
  if ( type(b)!="t_POLMOD",b=Mod(b,nf.pol));

  alpha = idealval(nf,a,p); beta = idealval(nf,b,p);
if (DEBUGLEVEL >= 5, print("[alpha,beta]=",[alpha,beta]));
  if ( (alpha%2 == 0) && (beta%2 == 0) ,
if (DEBUGLEVEL >= 5, print("fin de mynfhilbertp"));
    return (1));
  aux2 = idealnorm(nf,p)\2;
  if ( alpha%2 && beta%2 && aux2%2 ,sig = 1 ,sig = -1);
  if ( beta , aux = nfmodid2(nf,a^beta/b^alpha,p) , aux = nfmodid2(nf,b^alpha,p));
  aux = aux^aux2 + sig;
  aux = lift(lift(aux));
  if ( aux == 0 ,rep = 1 ,rep = (idealval(nf,aux,p)>= 1) );
if (DEBUGLEVEL >= 5, print("fin de mynfhilbertp"));
  if ( rep , return (1) , return (-1));
}
if (DEBUGLEVEL >= 4, print("ideallistfactor "));
{
ideallistfactor(nf,listfact)=
local(Slist,S1,test,i,j,k);
if (DEBUGLEVEL >= 5, print("entree dans ideallistfactor"));
  Slist = [];test=1;
  for ( i = 1 ,length(listfact),
    if ( listfact[i] == 0 ,next);
    S1 = idealfactor(nf,listfact[i])[,1];
    for( j = 1 ,length(S1),k=length(Slist);
      for( k = 1 ,length(Slist),
        if ( Slist[k] == S1[j] ,test=0;break));
      if ( test, Slist=concat(Slist,[S1[j]]) , test=1);
     ));
if (DEBUGLEVEL >= 5, print("fin de ideallistfactor"));
  return(Slist);
}
if (DEBUGLEVEL >= 4, print("mynfhilbert "));
{
mynfhilbert(nf,a,b)=
local(al,bl,r1,i,S);
if (DEBUGLEVEL >= 4, print("entree dans mynfhilbert",[a,b]));
\\ calcule le symbole de Hilbert quadratique global (a,b):
\\ =1 si l'equation X^2-aY^2-bZ^2=0 a une solution non triviale,
\\ =-1 sinon,
\\ a et b doivent etre non nuls.

  if ( a == 0 || b == 0 , print("0 argument in mynfhilbert"));
  al = lift(a); bl = lift(b);

\\ solutions locales aux places infinies reelles

  r1 = nf.sign[1];
  for ( i = 1, r1,
    if ( sign(mysubst(al,nf.roots[i])) < 0,
      if ( sign(mysubst(bl,nf.roots[i]))< 0 ,
if (DEBUGLEVEL >= 3, print("mynfhilbert not soluble at a infty"));
if (DEBUGLEVEL >= 4, print("fin de mynfhilbert"));
        return (-1))));

  if ( type(a) == "t_POLMOD" ,, a=Mod(a,nf.pol));
  if ( type(b) == "t_POLMOD" ,, b=Mod(b,nf.pol));

\\  solutions locales aux places finies (celles qui divisent 2ab)

  S = ideallistfactor (nf,[2,a,b]);
  forstep ( i = length(S), 2 ,-1 ,
\\ d'apres la formule du produit on peut eviter un premier
    if ( mynfhilbertp(nf,a,b, S[i]) == -1,
if (DEBUGLEVEL >= 3, print("mynfhilbert not soluble at finite place: ",S[i]));
if (DEBUGLEVEL >= 4, print("fin de mynfhilbert"));
      return (-1)));
if (DEBUGLEVEL >= 4, print("fin de mynfhilbert"));
  return(1);
}
if (DEBUGLEVEL >= 4, print("initp "));
{
initp( nf, p)=
\\   pp[1] est l'ideal sous forme reduite
\\   pp[2] est un entier de Zk avec une valuation 1 en p
\\   pp[3] est la valuation de 2 en p
\\   pp[4] sert a detecter les carres dans Qp
\\     si p|2 il faut la structure de Zk/p^(1+2v) d'apres Hensel
\\     sinon il suffit de calculer x^(N(p)-1)/2
\\   pp[5] est un systeme de representants de Zk/p
\\     c'est donc un ensemble de cardinal p^f .
local(idval,pp);

if (DEBUGLEVEL >= 5, print("entree dans initp"));
  idval=idealval(nf,2,p);
  pp=[ p , nfbasistoalg(nf,p[2]) , idval , 0 , repres(nf,p) ];
  if ( idval ,
    pp[4] = idealstar(nf,idealpow(nf,p,1+2*idval)),
    pp[4] = p[1]^p[4]\2 );
if (DEBUGLEVEL >= 5, print("fin de initp"));
  return (pp);
}
if (DEBUGLEVEL >= 4, print("deno "));
{
deno(num) =
\\ calcul un denominateur du polynome num

if (DEBUGLEVEL >= 5, print("entree dans deno"));
  if ( num == 0 ,
if (DEBUGLEVEL >= 5, print("fin de deno"));
    return (1));
  if (type(num) == "t_POL" ,
if (DEBUGLEVEL >= 5, print("fin de deno"));
    return (denominator(content(num))));
if (DEBUGLEVEL >= 5, print("fin de deno"));
  return (denominator(num));
}
if (DEBUGLEVEL >= 4, print("degre "));
{
degre(idegre)=
local(ideg,jdeg);
if (DEBUGLEVEL >= 5, print("entree dans degre ",idegre));
  ideg = idegre; jdeg = 0;
  while (ideg>>= 1, jdeg++);
if (DEBUGLEVEL >= 5, print("fin de degre ",jdeg));
  return (jdeg);
}
if (DEBUGLEVEL >= 4, print("square "));
{
square( nf, a, myplist)=
local(r1,i,ar,degsq);
if (DEBUGLEVEL >= 5, print("entree dans square",a));
  if ( a == 0 ,
if (DEBUGLEVEL >= 5, print("fin de square"));
    return (1));
  if ( a == 1 ,
if (DEBUGLEVEL >= 5, print("fin de square"));
    return (1));
  if (type(lift(a)) == "t_INT" && issquare(lift(a)) ,
if (DEBUGLEVEL >= 5, print("fin de square"));
    return (1));
\\ tous les plgments reels doivent etre >0
  r1 = nf.sign[1];
  for ( i = 1 ,r1,
    ar = mysubst(a,nf.roots[i]);
    if ( sign(simplify(ar)) < 0 ,
if (DEBUGLEVEL >= 5, print("fin de square"));
      return (0)));
\\ a doit etre un carre modulo les ideaux p
  for ( i = 1 ,length(myplist),
  if ( !psquare(nf,a,myplist[i][1],myplist[i][4]) ,
if (DEBUGLEVEL >= 5, print("fin de square"));
    return (0)));
\\ la norme de a doit etre un carre */
  if (!issquare(norm(a)) ,
if (DEBUGLEVEL >= 5, print("fin de square"));
    return (0));
\\ factorisation sur K du polynome X^2-a :
\\ renvoit 1 si a est un carre, 0 sinon
  degsq = poldegree( factornf(x^2-a,nf.pol)[1,1]);
if (DEBUGLEVEL >= 5, print("fin de square"));
  return ( degsq == 1);
}
if (DEBUGLEVEL >= 4, print("ratpoint "));
{
ratpoint(nf,pol,lim,myplist,flag)=
\\ si flag ==0, cherche un seul point, sinon plusieurs.
local(compt1,compt2,Np,deg,n,AA,point,listpoints,denoz,vectx,xx,evpol);

if (DEBUGLEVEL >= 4, print("entree dans ratpoint avec pol=",pol);print("lim = ",lim););
  compt1 = 0; compt2 = 0;
  Np = length(myplist); deg = poldegree(pol); n = poldegree(nf.pol);
  AA = lim<<1;
  if (flag, listpoints = []);
  point = [];

\\ cas triviaux
  if (square(nf,polcoeff(pol,0),myplist),
    point = [ 0 , mysqrtnf(nf,polcoeff(pol,0)) ,1 ];
if (DEBUGLEVEL >= 3, print("solution triviale: e est un carre"));
    if ( !flag ,
if (DEBUGLEVEL >= 4, print("fin de ratpoint"));
      return (point));
    listpoints = concat(listpoints,[point])
  );
  if ( square(nf,pollead(pol),myplist),
    point = [ 1 , mysqrtnf(nf,Mod(lift(pollead(pol)),nf.pol)) , 0];
if (DEBUGLEVEL >= 3, print("solution triviale: a est un carre"));
    if ( !flag ,
if (DEBUGLEVEL >= 4, print("fin de ratpoint"));
      return (point));
    listpoints = concat(listpoints,[point])
  );

\\ boucle generale
  vectx = vector(n,i,[-lim,lim]);
  for( denoz = 1 ,lim,
    forvec( xx = vectx ,
      if(denoz==1 || gcd(content(xx),denoz)==1 ,
        xpol = nfbasistoalg(nf,xx~);
        evpol = subst(pol,x,xpol/denoz);
        if ( square(nf,evpol,myplist) ,
          point=[xpol/denoz,mysqrtnf(nf,evpol),1];
          if ( !flag , break(2));
          if ( !flag ,
if (DEBUGLEVEL >= 4, print("fin de ratpoint"));
            return (point));
          listpoints = concat(listpoints,[point])));
  ));

if (DEBUGLEVEL >= 4, print("sortie de ratpoint"));
if (DEBUGLEVEL >= 3, print("points sur la quartique=",point));
  if (!flag ,return (point),return (listpoints));
}
if (DEBUGLEVEL >= 4, print("repres "));
{
repres(nf,p)=
\\ calcule un systeme de representants Zk/p
local(fond,mat,i,j,k,f,rep,pp,ppi,pp2,jppi,gjf);

if (DEBUGLEVEL >= 5, print("entree dans repres"));
  fond = [];
  mat = idealhnf(nf,p);
  for( i = 1 ,length(mat),
    if ( mat[i,i] != 1 ,fond = concat(fond,nf.zk[i])));
  f = length(fond);
  pp = p[1];
  rep = vector(pp^f,i,0);
  rep[1] = 0;
  ppi = 1;
  pp2 = pp\2;
  for ( i = 1 ,f,
    for ( j = 1 ,pp-1,
      if ( j <= pp2 ,gjf = j*fond[i], gjf = (j-pp)*fond[i]);
      jppi = j*ppi;
      for ( k = 0 ,ppi-1,rep[jppi+k+1] = rep[k+1]+gjf ));
    ppi *= pp);
if (DEBUGLEVEL >= 5, print("fin de repres"));
  return (Mod(rep,nf.pol));
}
if (DEBUGLEVEL >= 4, print("val "));
{
val(nf,num,p)=
local(res);
if (DEBUGLEVEL >= 5, print("entree dans val"));
  if ( num == 0 ,
if (DEBUGLEVEL >= 5, print("fin de val"));
    return (BIGINT));
  res = idealval(nf,lift(num),p);
if (DEBUGLEVEL >= 5, print("fin de val"));
  return (res);
}
if (DEBUGLEVEL >= 4, print("nfissquarep "));
{
nfissquarep(nf,a,p,q)=
\\ suppose que a est un carre modulo p^q
\\ et renvoit sqrt(a) mod p^q (ou plutot p^(q/2))
local(pherm,f,aaa,n,pp,qq,e,z,xx,yy,r,aux,b,m,vp,inv2x,zinit,zlog,expo);

if (DEBUGLEVEL >= 5, print("entree dans nfissquarep",a,p,q));
  if ( a == 0 || a == 1 ,
if (DEBUGLEVEL >= 4, print("fin de nfissquarep"));
    return (a));
  pherm = idealhnf(nf,p);
if (DEBUGLEVEL >= 5, print("pherm=",pherm));
  f = idealval(nf,a,p);
  if ( f >= q ,
    if ( f > q ,aaa = nfbasistoalg(nf,p[2])^((q+1)>>1) , aaa = 0);
if (DEBUGLEVEL >= 4, print("fin de nfissquarep"));
    return (aaa));
  if ( f , aaa = a*nfbasistoalg(nf,p[5]/p[1])^f , aaa = a);
  if ( pherm[1,1] != 2 ,
\\ cas ou p ne divise pas 2
\\ algorithme de Shanks
    n = nfrandintmodid(nf,pherm);
    while ( psquarenf(nf,n,p) , n=nfrandintmodid(nf,pherm));
    pp = Mod(1,p[1]);
    n *= pp;
    qq=idealnorm(nf,pherm)\2;
    e = 1; while ( !(qq%2) ,e++; qq \= 2);
    z = mynfeltreduce(nf,lift(lift(n^qq)),pherm);
    yy = z;r = e;
    xx = mynfeltreduce(nf,lift(lift((aaa*pp)^(qq\2))),pherm);
    aux = mynfeltreduce(nf,aaa*xx,pherm);
    b = mynfeltreduce(nf,aux*xx,pherm);
    xx = aux;
    aux = b;m = 0;
    while ( !val(nf,aux-1,p) ,m++; aux = mynfeltreduce(nf,aux^2,pherm));
    while( m,
      if( m==r , print("error in nfissquarep");1/0;break);
      yy *= pp;
      aux = mynfeltreduce(nf,lift(lift(yy^(1<<(r-m-1)))),pherm);
      yy = mynfeltreduce(nf,aux^2,pherm);
      r = m;
      xx = mynfeltreduce(nf,xx*aux,pherm);
      b = mynfeltreduce(nf,b*yy,pherm);
      aux = b;m = 0;
      while ( !val(nf,aux-1,p) , m++; aux = mynfeltreduce(nf,aux^2,pherm));
    );
\\ lift de Hensel
\\
    if( q > 1,
      vp = idealval(nf,xx^2-aaa,p);
      if( vp < q-f,
        yy = 2*xx;
        inv2x = nfbasistoalg(nf,idealaddtoone(nf,yy,p)[1])/yy;
        while (vp<q , vp++; xx-=(xx^2-aaa)*inv2x);
      );
      if (f ,xx*=nfbasistoalg(nf,p[2])^(f>>1));
    );
    xx = mynfeltreduce(nf,xx,idealpow(nf,p,q));
  ,
\\ cas ou p divise 2 */
    if( q-f>1 , id = idealpow(nf,p,q-f), id = pherm);
    zinit = idealstar(nf,id,2);
    zlog = ideallog(nf,aaa,zinit);
    xx = 1;
    for ( i = 1 , length(zlog),
      expo = zlog[i];
      if ( expo ,
        if ( !expo%2 , expo = expo>>1, aux = zinit[2][i]; expo = expo*((aux+1)>>1)%aux);
        if ( expo==1, xx *= nfbasistoalg(nf,(zinit[2][3][i])),
                     xx *= nfbasistoalg(nf,(zinit[2][3][i]))^expo);
      )
    );
    if ( f, xx*=nfbasistoalg(nf,p[2])^(f>>1);  id = idealpow(nf,p,q));
    xx = mynfeltreduce(nf,xx,id);
  );
if (DEBUGLEVEL >= 4, print("fin de nfissquarep",xx));
  return (xx);
}
if (DEBUGLEVEL >= 4, print("nfissquare"));
{
nfissquare( nf, a)=
\\ si a n'est pas un carre, alors renvoit [].
\\   sinon, renvoit [sqrt(a)]
local(ta,res,pfact,i,alift,r1,py);

if (DEBUGLEVEL >= 5, print("entree dans nfissquare ",a));
if ( a==0 || a==1,
if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
    return ([a]));
  res=[];

  ta = type(a);
  alift=lift(a);
  if(!poldegree(alift),alift=polcoeff(alift,0));

  if ( type(alift) != "t_POL" ,
    if( issquare(alift) ,
      res = [round(sqrt(alift))] ;
if (DEBUGLEVEL >= 5, print("fin de nfissquare 1"));
      return (res)));
  if( poldegree(nf.pol)<=1,
    res = [];
if (DEBUGLEVEL >= 5, print("fin de nfissquare 0"));
    return (res));
  if ( ta == "t_POL" ,a = Mod(a,nf.pol));

\\ a doit etre un carre modulo les ideaux p

\\  pfact = idealfactor(nf,a);
\\  for ( i = 1 ,length(pfact[,1]),
\\    if (pfact[i,2]%2 ,
\\if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
\\      return (res)));
\\  for ( i = 1 ,length(pfact[,1]),
\\    if ( pfact[i,1][1]%2 && !psquarenf(nf,a,pfact[i,1]) ,
\\if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
\\      return (res)));
\\

\\ tous les plgements reels doivent etre >0
\\
  r1=nf.sign[1];
  for ( i = 1 ,r1,
    py=mysubst(alift,nf.roots[i]);
    if ( sign(py) < 0 ,
if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
      return (res)));
\\ factorisation sur K du polynome X^2-a :
\\ renvoit 1 si a est un carre, 0 sinon
\\
  if (variable(nf.pol) == x,
    py = subst(nf.pol,x,y);
    pfact = lift(factornf( x^2-mysubst(alift,Mod(y,py)),py)[1,1]);
  ,
    pfact = lift(factornf(x^2-a,nf.pol)[1,1]));
  if ( poldegree(pfact) == 2 ,
if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
    return ([]));
if (DEBUGLEVEL >= 5, print("fin de nfissquare"));
  return ([subst(polcoeff(pfact,0),y,Mod(variable(nf.pol),nf.pol))]);
}
if (DEBUGLEVEL >= 4, print("psquarenf "));
{
psquarenf( nf, a, p)=
local(v,ap,norme,den);
\\ renvoie 1 si a est un carre dans ZK_p 0 sinon
\\ seulement pour p premier avec 2

if (DEBUGLEVEL >= 5, print("entree dans psquarenf"));
if ( a == 0 ,
if (DEBUGLEVEL >= 5, print("fin de psquarenf"));
    return (1));
  v = idealval(nf,lift(a),p);
  if ( v%2 ,
if (DEBUGLEVEL >= 5, print("fin de psquarenf"));
    return (0));
  ap = a*(1/nfbasistoalg(nf,p[2])^v);
\\
  norme = idealnorm(nf,p)\2;
  den = denominator(content(lift(ap)))%p[1];
  if (sign(den), ap*=Mod(1,p[1]));
  ap = ap^norme-1;
  if ( ap == 0 ,
if (DEBUGLEVEL >= 5, print("fin de psquarenf"));
    return (1));
  ap = lift(lift(ap));
  if ( idealval(nf,ap,p) > 0 ,
if (DEBUGLEVEL >= 5, print("fin de psquarenf"));
    return (1));
if (DEBUGLEVEL >= 5, print("fin de psquarenf"));
  return (0);
}
if (DEBUGLEVEL >= 4, print("psquare "));
{
psquare( nf, a, p, zinit)=
\\ a est un entier de K
\\ renvoie 1 si a est un carre dans ZKp 0 sinon
local(valap,zlog,i);

if (DEBUGLEVEL >= 5, print("entree dans psquare",[a,p,zinit]));
  if ( a == 0,
if (DEBUGLEVEL >= 5, print("fin de psquare"));
    return (1));
\\
  if ( (p[1]%2) ,
if (DEBUGLEVEL >= 5, print("fin de psquare"));
    return (psquarenf(nf,a,p)));
\\
  valap=idealval(nf,a,p);
  if ( valap%2 ,
if (DEBUGLEVEL >= 5, print("fin de psquare"));
    return (0));
  if(valap,
    zlog = ideallog(nf,a*(nfbasistoalg(nf,p[5])/p[1])^valap,zinit)
  ,
    zlog = ideallog(nf,a,zinit));
  for ( i = 1 ,length(zinit[2][2]),
    if ( !(zinit[2][2][i]%2) && (zlog[i]%2) ,
if (DEBUGLEVEL >= 5, print("fin de psquare"));
      return (0)));
if (DEBUGLEVEL >= 5, print("fin de psquare"));
  return (1);
}
if (DEBUGLEVEL >= 4, print("psquareq "));
{
psquareq( nf, a, p, q)=
\\ cette fonction renvoie 1 si a est un carre
\\ ?inversible? modulo P^q et 0 sinon.
\\ P divise 2, et ?(a,p)=1?.
local(vala,zinit,zlog,i);

if (DEBUGLEVEL >= 5, print("entree dans psquareq",[a,p,q]));
  if(a==0,
if (DEBUGLEVEL >= 5, print("fin de psquareq"));
    return (1));
  vala=idealval(nf,a,p);
  if(vala>=q,
if (DEBUGLEVEL >= 5, print("fin de psquareq"));
    return (1));
  if(vala%2,
if (DEBUGLEVEL >= 5, print("fin de psquareq"));
    return (0));
  zinit = idealstar(nf,idealpow(nf,p,q-vala),2);
  zlog = ideallog(nf,a*nfbasistoalg(nf,p[5]/2)^vala,zinit);
  for ( i = 1 ,length(zinit[2][2]),
    if ( !(zinit[2][2][i]%2) && (zlog[i]%2) ,
if (DEBUGLEVEL >= 5, print("fin de psquareq"));
      return (0)));
if (DEBUGLEVEL >= 5, print("fin de psquareq"));
  return (1);
}
if (DEBUGLEVEL >= 4, print("lemma6nf "));
{
lemma6nf( nf, pol, p, nu, xx)=
local(gx,gpx,lambda,mu);

if (DEBUGLEVEL >= 5, print("entree dans lemma6nf"));
  gx = subst( pol , x , xx);
  if ( psquarenf(nf,gx,p) ,
if (DEBUGLEVEL >= 5, print("fin de lemma6nf"));
    return (1));
  gpx = subst( pol' , x , xx);
  lambda = val(nf,gx,p);mu = val(nf,gpx,p);

  if ( lambda>(2*mu) ,
if (DEBUGLEVEL >= 5, print("fin de lemma6nf"));
    return (1));
  if ( (lambda >= 2*nu)  && (mu >= nu) ,
if (DEBUGLEVEL >= 5, print("fin de lemma6nf"));
    return (0));
if (DEBUGLEVEL >= 5, print("fin de lemma6nf"));
  return (-1);
}
if (DEBUGLEVEL >= 4, print("lemma7nf "));
{
lemma7nf( nf, pol, p, nu, xx, zinit)=
local(gx,gpx,v,lambda,mu,q);

if (DEBUGLEVEL >= 5, print("entree dans lemma7nf",[xx,nu]));
  gx = subst( pol , x , xx);
  if ( psquare(nf,gx,p,zinit) ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
    return (1));
  gpx = subst( pol' , x , xx);
  v = p[3];
  lambda = val(nf,gx,p);mu = val(nf,gpx,p);
  if ( lambda>2*mu ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
    return (1));
  if (nu > mu ,
    if ( lambda%2 ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
      return (-1));
    q = mu+nu-lambda;
    if ( q>2*v ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
      return (-1));
    if ( psquareq(nf,gx*nfbasistoalg(nf,p[5]/2)^lambda,p,q) ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
      return (1));
  ,
    if ( lambda>= 2*nu ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
      return (0));
    if ( lambda%2 ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
      return (-1));
    q = 2*nu-lambda;
    if ( q > 2*v ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
      return (-1));
    if ( psquareq(nf,gx*nfbasistoalg(nf,p[5]/2)^lambda,p,q) ,
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
      return (0));
  );
if (DEBUGLEVEL >= 5, print("fin de lemma7nf"));
  return (-1);
}
if (DEBUGLEVEL >= 4, print("zpsolublenf "));
{
zpsolublenf( nf, pol, p, nu, pnu, x0)=
local(result,pnup,lrep);

if (DEBUGLEVEL >= 5, print("entree dans zpsolublenf",[lift(x0),nu]));
  if ( p[3] == 0 ,
    result = lemma6nf(nf,pol,p[1],nu,x0),
    result = lemma7nf(nf,pol,p[1],nu,x0,p[4]));
  if ( result == +1 ,
if (DEBUGLEVEL >= 5, print("fin de zpsolublenf"));
    return (1));
  if ( result == -1 ,
if (DEBUGLEVEL >= 5, print("fin de zpsolublenf"));
    return (0));
  pnup=pnu*p[2];
  lrep=length(p[5]);
  nu++;
  for ( i = 1 ,lrep,
    if ( zpsolublenf(nf,pol,p,nu,pnup,x0+pnu*p[5][i]) ,
if (DEBUGLEVEL >= 5, print("fin de zpsolublenf"));
      return (1)));
if (DEBUGLEVEL >= 5, print("fin de zpsolublenf"));
  return (0);
}

if (DEBUGLEVEL >= 4, print("mynfeltmod "));
{
mynfeltmod(nf,a,b)=
local(qred);
qred = round(nfalgtobasis(nf,a/b));
qred = a-b*nfbasistoalg(nf,qred);
return (qred);
}
if (DEBUGLEVEL >= 4, print("mynfeltreduce "));
{
mynfeltreduce(nf,a,id) = nfbasistoalg(nf,nfeltreduce(nf,nfalgtobasis(nf,a),id))
}
if (DEBUGLEVEL >= 4, print("nfrandintmodid "));
{
nfrandintmodid( nf, id)=
local(res);

if (DEBUGLEVEL >= 5, print("entree dans nfrandintmodid"));
  res = 0;
  while( !res,
    res=nfrandint(nf,0);
    res=mynfeltreduce(nf,res,id));
if (DEBUGLEVEL >= 5, print("fin de nfrandintmodid"));
  return (res);
}
if (DEBUGLEVEL >= 4, print("nfrandint "));
{
nfrandint( nf, borne)=
local(l,res,i);

if (DEBUGLEVEL >= 5, print("entree dans nfrandint"));
  l = length(nf.zk);
  res = vectorv(l,i,0);
  for ( i = 1 ,l,
      if ( borne , res[i] = random(borne<<1)-borne ,res[i] = random() ));
  res = nfbasistoalg(nf,res);
if (DEBUGLEVEL >= 5, print("fin de nfrandint"));
  return (res);
}
if (DEBUGLEVEL >= 4, print("qpsolublenfbig "));
{
qpsolublenfbig( nf, pol, p,ap=0,b=1) =
local(deg,i,xx,z,Px,j,cont,pi,pol2,Roots);
if (DEBUGLEVEL >= 4, print("entree dans qpsolublenfbig avec ",p[1]));
  deg = poldegree(pol);
\\
  if ( psquarenf(nf,polcoeff(pol,0),p) ,
if (DEBUGLEVEL >= 4, print("fin de qpsolublenfbig"));
    return (1));
  if ( psquarenf(nf,pollead(pol),p) ,
if (DEBUGLEVEL >= 4, print("fin de qpsolublenfbig"));
    return (1));

\\ on tient compte du contenu de pol
  cont=idealval(nf,polcoeff(pol,0),p);
  for(i=1,deg,
    if(cont, cont=min(cont,idealval(nf,polcoeff(pol,i),p))));
  if(cont, pi = nfbasistoalg(nf,p[5]/p[1]));
  if(cont>1, pol *= pi^(2*(cont\2)));

\\ On essaye des valeurs de x au hasard
  if( cont%2,
    pol2 = pol*pi
  , pol2 = pol;
    for(i=1 ,MAXPROB ,
      xx = nfrandint(nf,0);
      z = 0; while (!z, z=random());
\\      xx /= z;
      xx = -ap*z+b*xx;
      Px=polcoeff(pol,deg);
      forstep (j=deg-1,0,-1,Px=Px*xx+polcoeff(pol,j));
      Px *= z^(deg);
      if ( psquarenf(nf,Px,p) ,
if (DEBUGLEVEL >= 4, print("fin de qpsolublenfbig"));
        return (1));
    )
  );

\\ On essaye les racines de pol
  Roots = nfpolrootsmod(nf,pol2,p);
  pi = nfbasistoalg(nf,p[2]);
  for(i=1,length(Roots),
    if(qpsolublenfbig(nf,subst(pol,x,pi*x+Roots[i]),p),return(1)));

if (DEBUGLEVEL >= 4, print("fin de qpsolublenfbig"));
  return (0);
}
if (DEBUGLEVEL >= 4, print("nfpolrootsmod "));
{nfpolrootsmod(nf,pol,p)=
\\ calcule les racines modulo l'ideal p du polynome pol.
\\ p est un ideal premier de nf, sous la forme idealprimedec
local(factlist,sol);
  factlist = nffactormod(nf,pol,p);
\\ CETTE LIGNE NE DOIT PAS RESTER TRES LONGTEMPS
  if(type(factlist)=="t_VEC",factlist=factlist[1],factlist=factlist[,1]);
  sol=[];
  for(i=1,length(factlist),
    if(poldegree(factlist[i])==1,
      sol = concat(sol, [-polcoeff(factlist[i],0)/polcoeff(factlist[i],1)])));
  return(sol);
}
if (DEBUGLEVEL >= 4, print("qpsolublenf "));
{
qpsolublenf( nf, pol, p) =

if (DEBUGLEVEL >= 4, print("entree dans qpsolublenf ",p));
if (DEBUGLEVEL >= 5, print("pol = ",pol));
  if ( psquare(nf,pollead(pol),p[1],p[4]) ,
if (DEBUGLEVEL >= 5, print("fin de qpsolublenf"));
    return (1));
  if ( psquare(nf,polcoeff(pol,0),p[1],p[4]) ,
if (DEBUGLEVEL >= 5, print("fin de qpsolublenf"));
    return (1));
  if ( zpsolublenf(nf,pol,p,0,1,0) ,
if (DEBUGLEVEL >= 5, print("fin de qpsolublenf"));
    return (1));
  if ( zpsolublenf(nf,polrecip(pol),p,1, p[2],0) ,
if (DEBUGLEVEL >= 5, print("fin de qpsolublenf"));
    return (1));
if (DEBUGLEVEL >= 5, print("fin de qpsolublenf"));
  return (0);
}
if (DEBUGLEVEL >= 4, print("locallysoluble "));
{
locallysoluble( nf, pol, r=0,a=1,b=1)=
local(pol0,plist,add,ff,p,r1,Delta,vecpol,ar,er,Deltar,vecpolr,Sturmr);

if (DEBUGLEVEL >= 4, print("entree dans locallysoluble",[pol,r,a,b]));
  pol0=pol;
\\
\\ places finies de plist */
\\
  pol *= deno(content(lift(pol)))^2;
  for( ii = 1 , 3,
    if ( ii == 1 , plist = idealprimedec(nf,2));
    if ( ii == 2 && r , plist = idealfactor(nf,poldisc(pol0/pollead(pol0))/pollead(pol0)^6/2^12)[,1]);
    if ( ii == 2 && !r , plist = idealfactor(nf,poldisc(pol0))[,1]);
    if ( ii == 3 ,
      add=idealadd(nf,a,b);
      ff=factor(idealnorm(nf,add))[,1];
      addprimes(ff);
if (DEBUGLEVEL >= 4, print("liste de premiers = ",ff));
      plist = idealfactor(nf,add)[,1]);
      for ( i =1 ,length(plist),
        p =  plist[i];
if (DEBUGLEVEL >= 3, print("p=",p));
        if ( p[1] < LIMBIGPRIME ,
          if ( !qpsolublenf(nf,pol,initp(nf,p)) ,
if (DEBUGLEVEL >= 2, print(" not ELS at ",p));
if (DEBUGLEVEL >= 4, print("fin de locallysoluble"));
            return (0)),
          if ( !qpsolublenfbig(nf,pol,p,r/a,b) ,
if (DEBUGLEVEL >= 2, print(" not ELS at big ",p[1]));
if (DEBUGLEVEL >= 4, print("fin de locallysoluble"));
            return (0))));
);
\\ places reelles
  r1 = nf.sign[1];
  if ( r1,
    Delta = poldisc(pol); vecpol = Vec(pol);
    for ( i = 1 ,r1,
      ar = mysubst(pollead(pol),nf.roots[i]);
      if ( ar > 0 ,next);
      er = mysubst(polcoeff(pol,0),nf.roots[i]);
      if ( er > 0 ,next);
      Deltar = mysubst(Delta,nf.roots[i]);
      if ( Deltar < 0 ,next);
      vecpolr = vector(poldegree(pol)+1,j,mysubst(vecpol[j],nf.roots[i]));
      Sturmr = polsturm(Pol(vecpolr));
      if ( Sturmr == 0 ,
if (DEBUGLEVEL >= 2, print(" not ELS at infty "));
if (DEBUGLEVEL >= 4, print("fin de locallysoluble"));
        return (0));
  ));
if (DEBUGLEVEL >= 2, print(" quartic ELS "));
if (DEBUGLEVEL >= 4, print("fin de locallysoluble"));
  return (1);
}
if (DEBUGLEVEL >= 4, print("mysqrtnf "));
{
mysqrtnf( nf, a)=
local( aux);
if (DEBUGLEVEL >= 5, print("entree dans mysqrtnf"));
  if ( a == 0 ,
if (DEBUGLEVEL >= 5, print("fin de mysqrtnf"));
    return (0));
  if ( a == 1 ,
if (DEBUGLEVEL >= 5, print("fin de mysqrtnf"));
    return (1));
  aux = factornf(x^2-a,nf.pol)[1,1];
  if (poldegree(aux)==2, print("********** BUG : NON CARRE DANS mysqrtnf");1/0;);
if (DEBUGLEVEL >= 5, print("fin de mysqrtnf"));
  return (polcoeff(aux,0));
}
if (DEBUGLEVEL >= 4, print("count "));
{
count( nf, c, d, KS2gen, myplist, pointstriv)=
local(found,listgen,listpointscount,m1,m2,lastloc,mask,i,d1,iaux,j,triv,pol,point,deuxpoints);

if (DEBUGLEVEL >= 4, print("entree dans count ",[c,d]));
  found = 0;
  listgen = KS2gen;
  listpointscount = [];

  m1 = m2 = 0;lastloc = -1;

  mask = 1 << length(KS2gen) ;
  i=1;
  while (i < mask ,
    d1 = 1; iaux = i; j = 1;
    while ( iaux ,
      if ( iaux%2 ,d1 *= listgen[j]);
      iaux >>= 1; j++);
if (DEBUGLEVEL >= 2, print("d1=",d1));
    triv = 0;
    for(j=1,length(pointstriv),
      if(pointstriv[j][3]*pointstriv[j][1],
        if(square(nf,d1*pointstriv[j][1]*pointstriv[j][3],myplist),
          listpointscount=concat(listpointscount,[pointstriv[j]]);
if (DEBUGLEVEL >= 2, print("point trivial"));
          triv = 1; m1++;
          if ( degre(i) > lastloc ,m2++);
          found = 1; lastloc = -1;break)));
    if(!triv ,
    pol = (d1*x^2+c)*x^2+d/d1;
if (DEBUGLEVEL >= 3, print("quartic = y^2 = ",pol));
    point = ratpoint(nf,pol,LIM1,myplist,0);
    if ( point != [],
if (DEBUGLEVEL >= 2, print ("point sur la quartique"));
if (DEBUGLEVEL >= 3, print (point));
      m1++;
      if ( point[3] != 0 ,
        aux = d1*point[1]/point[3]^2;
        deuxpoints = [ aux*point[1] , aux*point[2]/point[3] ];
      ,
        deuxpoints = [0]);
      listpointscount = concat(listpointscount,[deuxpoints]);
      if ( degre(i) > lastloc ,m2++);
      found = 1;lastloc = -1;
    ,
      if ( locallysoluble(nf,pol),
        if ( degre(i) > lastloc ,m2++;lastloc=degre(i));
        point = ratpoint(nf,pol,LIM3,myplist,0);
        if ( point != [],
if (DEBUGLEVEL >= 2, print ("point sur la quartique"));
if (DEBUGLEVEL >= 3, print (point));
          m1++;
          aux = d1*point[1]/point[3]^2;
          deuxpoints = [ aux*point[1] , aux*point[2]/point[3] ];
          listpointscount = concat(listpointscount,[deuxpoints]);
          if ( degre(i) > lastloc ,m2++);
          found = 1;lastloc = -1);
        ,if (DEBUGLEVEL >= 2, print("pas de point trouve sur la quartique"));
          )));
    if ( found,
      found = 0;
      v = 0;iaux = (i>>1);
      while ( iaux ,iaux >>= 1;v++);
      mask >>= 1;
      listgen = vecextract(listgen,(1<<(length(listgen)))-(1<<v)-1);
      i=(1<<v);
    ,i++);
  );
  for(i=1,length(listpointscount),
   if(length(listpointscount[i])>1,
    if(subst(x^3+c*x^2+d*x,x,listpointscount[i][1])-listpointscount[i][2]^2 != 0,
      print("********* MAUVAIS POINT DANS COUNT ****** =",listpointscount[i]))));
if (DEBUGLEVEL >= 4, print("fin de count"));
  return ([listpointscount,[m1,m2]]);
}
if (DEBUGLEVEL >= 4, print("makemyplist "));
{
makemyplist( nf, nbid)=
local(i,p,myplist,v,pdec,j);
if (DEBUGLEVEL >= 4, print("entree dans makemyplist"));
  i=0;p=1;
  myplist=[];v=[0];

  while ( i<=nbid ,
    p+=2;
    if (isprime(p),
      pdec = idealprimedec(nf,p);
      j=1;
      while (j<=length(pdec) && j<=2 ,
        if ( pdec[j][4] == 1,
          i++;
          v[1] = initp(nf, pdec[j]);
          myplist = concat(myplist,v));
        j++);
  ));
if (DEBUGLEVEL >= 4, print("fin de makemyplist"));
  return (myplist);
}
if (DEBUGLEVEL >= 4, print("main2 "));
{
main2( bnf, Ell )=
\\ Calcul du rang des courbes elliptiques avec 2-torsion
\\ dans le corps de nombres bnf
\\ par la methode des 2-isogenies.
\\
\\ Ell = [a1,a2,a3,a4,a6]
\\ y^2+a1xy+a3y=x^3+a2x^2+a4x+a6
\\
\\ Ell doit etre sous la forme
\\ y^2=x^3+ax^2+bx -> Ell = [0,a,0,b,0]
\\ avec a et b entiers.
local(i,P,Pfact,tors,pointstriv,apinit,bpinit,plist,myplist,KS2prod,oddclass,KS2gen,listpoints,pointgen,n1,n2,certain,np1,np2,listpoints2,aux1,aux2,certainp,rang,strange);

if (DEBUGLEVEL >= 3, print("entree dans main2"));
  if ( variable(bnf.pol)!=y,
    print(" ********** la variable du corps de nombres doit etre y ");return(0));
  Ell=ellinit(Mod(lift(Ell),bnf.pol),1);

  if ( Ell.disc == 0 ,
    print(" ************* discriminant=0 !!");return(0));
  if ( Ell.a6 != 0 ,
    print(" ************* la courbe n'est pas sous la forme y^2 = y^2=x^3+a*x^2+b*x dans main2 !!");return(0));
\\  a=Ell.a2;
  if( denominator(nfalgtobasis(bnf,Ell.a2)) > 1,
    print(" *** WARNING ****");
    print(" ************ coefficients non entiers dans main2"));
\\  b=Ell.a4;
  if( denominator(nfalgtobasis(bnf,Ell.a4)) > 1,
    print(" *** WARNING ****");
    print(" ************ coefficients non entiers dans main2"));

  P = (x^2+Ell.a2*x+Ell.a4)*Mod(1,bnf.pol);
  Pfact= factornf(P,bnf.pol)[,1];
  tors=length(Pfact);
  if(length(Pfact)>1,
    pointstriv=[[0,0,1],[-polcoeff(Pfact[1],0),0,1],[-polcoeff(Pfact[2],0),0,1]]
    , pointstriv = [[0,0,1]]);
\\
  apinit = -2*Ell.a2;bpinit = Ell.a2^2-4*Ell.a4;

\\ calcul des ideaux premiers de plist
\\ et de quelques renseignements associes
  plist = idealfactor(bnf,6*Ell.disc)[,1];

  myplist = makemyplist(bnf.nf,NBIDEAUX);

if (DEBUGLEVEL >= 3, print(" recherche de points triviaux sur la courbe "));
  P *= x;
if (DEBUGLEVEL >= 3, print("Y^2 = ",P));
  pointstriv = concat(pointstriv , ratpoint(bnf.nf,P,LIMTRIV,myplist,1));
if (DEBUGLEVEL >= 1, print("points triviaux sur E(K)= ");
  print(lift(pointstriv));print());

  KS2prod = Ell.a4;
  oddclass = 0;
  while ( !oddclass ,
    KS2gen = bnfsunit(bnf,idealfactor(bnf,KS2prod)[,1]~);
    oddclass = (KS2gen[5][1]%2);
    if (!oddclass,
      KS2prod = idealmul(bnf,KS2prod,(KS2gen[5][3][1])));
 );
  KS2gen = KS2gen[1];
  for( i = 1 ,length(KS2gen),
    KS2gen[i] = nfbasistoalg(bnf, KS2gen[i]));
  KS2gen = concat(Mod(lift(bnf.tufu),bnf.pol),KS2gen);
if (DEBUGLEVEL >= 2,
  print("#K(b,2)gen          = ",length(KS2gen));
  print("K(b,2)gen=",KS2gen));

  listpoints = count(bnf.nf,Ell.a2,Ell.a4,KS2gen,myplist,pointstriv);
  pointgen = listpoints[1];
if (DEBUGLEVEL >= 1, print("points sur E(K) = ",lift(pointgen));print());
  n1=listpoints[2][1];n2=listpoints[2][2];
\\
  certain = (n1==n2);
if (DEBUGLEVEL >= 1,
  if (certain,
    print("[E(K):phi'(E'(K))]  = ",1<<n1);
    print("#S^(phi')(E'/K)     = ",1<<n2);
    print("#III(E'/K)[phi']    = 1");print();
  ,
    print("[E(K):phi'(E'(K))] >= ",1<<n1);
    print("#S^(phi')(E'/K)     = ",1<<n2);
    print("#III(E'/K)[phi']   <= ",1<<(n2-n1));print());
);
\\
  KS2prod = bpinit;
  oddclass = 0;
  while (!oddclass ,
    KS2gen = bnfsunit(bnf,idealfactor(bnf,KS2prod)[,1]~);
    oddclass = (KS2gen[5][1]%2);
    if ( !oddclass ,
      KS2prod = idealmul(bnf,KS2prod,(KS2gen[5][3][1]))));
  KS2gen = KS2gen[1];
  for ( i = 1 ,length(KS2gen),
    KS2gen[i]= nfbasistoalg(bnf, KS2gen[i]));
  KS2gen = concat(Mod(lift(bnf.tufu),bnf.pol),KS2gen);
if (DEBUGLEVEL >= 2,
  print("#K(a^2-4b,2)gen     = ",length(KS2gen));
  print("K(a^2-4b,2)gen     = ",KS2gen));

  P = (x^2+apinit*x+bpinit)*Mod(1,bnf.pol);
  Pfact= factornf(P,bnf.pol)[,1];
  if(length(Pfact)>1,
    pointstriv=[[0,0,1],[-polcoeff(Pfact[1],0),0,1],[-polcoeff(Pfact[2],0),0,1]]
    , pointstriv = [[0,0,1]]);

if (DEBUGLEVEL >= 3, print(" recherche de points triviaux sur la courbe "));
  P *= x;
if (DEBUGLEVEL >= 3, print("Y^2 = ",P));
  pointstriv = concat(pointstriv , ratpoint(bnf.nf,P,LIMTRIV,myplist,1));
if (DEBUGLEVEL >= 1, print("points triviaux sur E'(K)= ");
  print(lift(pointstriv));print());

  listpoints = count(bnf.nf,apinit,bpinit,KS2gen,myplist,pointstriv);
if (DEBUGLEVEL >= 1, print("points sur E'(K) = ",lift(listpoints[1])));
  np1=listpoints[2][1];np2=listpoints[2][2];
  listpoints2 = vector(length(listpoints[1]),i,0);
  for ( i = 1 ,length(listpoints[1]),
    listpoints2[i]=[0,0];
    aux1 = listpoints[1][i][1]^2;
    if ( aux1 != 0,
      aux2 = listpoints[1][i][2];
      listpoints2[i][1] = aux2^2/aux1/4;
      listpoints2[i][2] = aux2*(bpinit-aux1)/aux1/8;
    , listpoints2[i] = listpoints[1][i]));
if (DEBUGLEVEL >= 1, print("points sur E(K) = ",lift(listpoints2));print());
  pointgen = concat(pointgen,listpoints2);
\\
  certainp = (np1==np2);
if (DEBUGLEVEL >= 1,
  if (certainp,
   print("[E'(K):phi(E(K))]   = ",1<<np1);
   print("#S^(phi)(E/K)       = ",1<<np2);
   print("#III(E/K)[phi]      = 1");print();
  ,
   print("[E'(K):phi(E(K))]  >= ",1<<np1);
   print("#S^(phi)(E/K)       = ",1<<np2);
   print("#III(E/K)[phi]     <= ",1<<(np2-np1));print());
\\
  if (!certain && (np2>np1) , print1(1<<(np2-np1)," <= "));
  print1("#III(E/K)[2]       ");
  if (certain && certainp , print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));
\\
  print("#E(K)[2]            = ",1<<tors);
);
  rang=n1+np1-2;
if (DEBUGLEVEL >= 1,
  if (certain && certainp,
    print("#E(K)/2E(K)         = ",(1<<(rang+tors)));
    print("rang                = ",rang);print();
  ,
    print("#E(K)/2E(K)        >= ",(1<<(rang+tors)));print();
    print(rang," <= rang          <= ",n2+np2-2);print()
  ));
\\
  strange = (n2+np2-n1-np1)%2;
  if (strange,
if (DEBUGLEVEL >= 1,
      print(" !!! III doit etre un carre !!!");print("donc"));
    if (certain,
      np1++;
      certainp = (np1==np2);
if (DEBUGLEVEL >= 1,
        if (certainp,
          print("[E'(K):phi(E(K))]   = ",1<<np1);
          print("#S^(phi)(E/K)       = ",1<<np2);
          print("#III(E/K)[phi]      = 1");print();
        ,
          print("[E'(K):phi(E(K))]  >= ",1<<np1);
          print("#S^(phi)(E/K)       = ",1<<np2);
          print("#III(E/K)[phi]     <= ",1<<(np2-np1));print());
      );
    ,
    if (certainp,
      n1++;
      certain = (n1==n2);
if (DEBUGLEVEL >= 1,
        if (certain,
          print("[E(K):phi'(E'(K))]   = ",1<<n1);
          print("#S^(phi')(E'/K)       = ",1<<n2);
          print("#III(E'/K)[phi']      = 1");print();
        ,
          print("[E(K):phi'(E'(K))]  >= ",1<<n1);
          print("#S^(phi')(E'/K)      = ",1<<n2);
          print("#III(E'/K)[phi']    <= ",1<<(n2-n1));print())
      );
    ,n1++)
  );
\\
if (DEBUGLEVEL >= 1,
  if (!certain && (np2>np1) , print1(1<<(np2-np1)," <= "));
  print1("#III(E/K)[2]       ");
  if (certain && certainp , print1(" "), print1("<"));
  print("= ",1<<(n2+np2-n1-np1));
  print("#E(K)[2]            = ",1<<tors);
);
  rang=n1+np1-2;
if (DEBUGLEVEL >= 1,
  if (certain && certainp,
    print("#E(K)/2E(K)         = ",(1<<(rang+tors)));print();
    print("rang                = ",rang);print();
  ,
    print("#E(K)/2E(K)        >= ",(1<<(rang+tors)));print();
    print(rang," <= rang          <= ",n2+np2-2);print())
  ));

\\ fin de strange

if (DEBUGLEVEL >= 1, print ("points = ",pointgen));
if (DEBUGLEVEL >= 3, print("fin de main2"));
  return (rang);
}
if (DEBUGLEVEL >= 4, print("chinremain "));
{
chinremain( nf, b, fact)=
\\ Chinese Ramainder Theorem
local(l,fact2,i);

if (DEBUGLEVEL >= 4, print("entree dans chinremain"));
  l = length(fact[,1]) ;
  fact2=vector(l,i,0);
  for ( i = 1 ,l,
    fact2[i] = idealdiv(nf,b,idealpow(nf,fact[i,1],fact[i,2])));
  fact2 = idealaddtoone(nf,fact2);
  for ( i = 1 ,l,
    fact2[i] = nfbasistoalg(nf,fact2[i]));
if (DEBUGLEVEL >= 4, print("fin de chinremain"));
  return (fact2);
}
if (DEBUGLEVEL >= 4, print("bnfqfsolve2 "));
{
bnfqfsolve2(bnf, aleg, bleg, auto=[y])=
\\ Solves Legendre Equation x^2-aleg*Y^2=bleg*Z^2
\\ Using quadratic norm equations
\\ auto contient les automorphismes de bnf sous forme de polynomes
\\ en y, avec auto[1]=y .
local(aux,solvepolrel,auxsolve,solvepolabs,exprxy,rrrnf,bbbnf,SL0,i,SL1,SL,sunL,fondsunL,normfondsunL,SK,sunK,fondsunK,vecbleg,matnorm,matnormmod,expsolution,solution,reste,carre,verif);

if (DEBUGLEVEL >= 3, print("entree dans bnfqfsolve2"));
  solvepolrel = x^2-aleg;
if (DEBUGLEVEL >= 4, print("aleg=",aleg));
if (DEBUGLEVEL >= 4, print("bleg=",bleg));

  if (length(auto)>1,
if (DEBUGLEVEL >= 4, print("factorisation du discriminant avec les automorhpismes de bnf"));
    for ( i = 2, length(auto),
      aux = abs(polresultant(lift(aleg)-subst(lift(aleg),y,auto[i]),bnf.pol));
      if(aux, addprimes(factor(aux)[,1]))));

  auxsolve = rnfequation(bnf,solvepolrel,1);
  solvepolabs = auxsolve[1];
  exprxy = auxsolve[2];
if ( auxsolve[3],
if( DEBUGLEVEL >=5, print(" CECI EST LE NOUVEAU CAS auxsolve[3]!=0 ")));
if (DEBUGLEVEL >= 4, print(" bbbnfinit ",solvepolabs));
  rrrnf = rnfinit(bnf,solvepolrel);
  bbbnf = bnfinit(solvepolabs,1);
if (DEBUGLEVEL >= 4, print(" done"));
  SL0 = 1;
if (DEBUGLEVEL >= 4, print("bbbnf.clgp=",bbbnf.clgp));
  for ( i = 1 ,length(bbbnf.clgp[2]),
    if ( bbbnf.clgp[2][i]%2 == 0,
      SL0 = idealmul(bbbnf,SL0,bbbnf.clgp[3][i])));
  SL1 = idealmul(bbbnf,SL0,rnfeltup(rrrnf,bleg));
  SL = idealfactor(bbbnf,SL1)[,1]~;
  sunL = bnfsunit(bbbnf,SL);
  fondsunL = concat(bbbnf.futu,nfbasistoalg(bbbnf,sunL[1]));
  normfondsunL = norm(rnfeltabstorel( rrrnf,fondsunL));
  SK = idealfactor(bnf,idealnorm(bbbnf,SL1))[,1]~;
  sunK = bnfsunit(bnf,SK);
  fondsunK = concat(bnf.futu,nfbasistoalg(bnf,sunK[1]));
  vecbleg = bnfissunit(bnf,sunK,bleg);
  matnorm = matrix(length(fondsunK),length(normfondsunL),i,j,0);
  for(i=1,length(normfondsunL),
    matnorm[,i] = lift(bnfissunit( bnf,sunK,normfondsunL[i] )));
  matnormmod = matnorm*Mod(1,2);
  expsolution = lift(matinverseimage(matnormmod , vecbleg*Mod(1,2)));
if (expsolution == []~, print("IL N'Y A PAS DE SOLUTION DANS bnfqfsolve2 ");1/0);
  solution = prod (i=1,length(expsolution) , fondsunL[i]^expsolution[i]);
  solution = rnfeltabstorel(rrrnf,solution);
  reste = (lift(vecbleg) - matnorm*expsolution)/2;
  carre = prod (i=1,length(vecbleg) , fondsunK[i]^reste[i]);
  solution *= carre;
  x1=polcoeff(lift(solution),1,x);x0=polcoeff(lift(solution),0,x);
  verif = x0^2 - aleg*x1^2-bleg;
if (verif , print("ERREUR DANS bnfqfsolve2");1/0);
if (DEBUGLEVEL >= 3, print("fin de bnfqfsolve2"));
  return ([x0,x1,1]);
}
if (DEBUGLEVEL >= 4, print("bnfqfsolve "));
{
bnfqfsolve(bnf, aleg, bleg, myplist,flag3,auto=[y])=
\\ cette fonction resout l'equation X^2-aleg*Y^2=bleg*Z^2
\\ dans le corps de nombres nf.
\\ la solution est [X,Y,Z],
\\ [0,0,0] sinon.
local(nf,aa,bb,na,nb,maxna,maxnb,mat,resl,t,sq,pol,vecrat,alpha,xx,yy,borne,test,sun,fact,suni,k,f,l,aux,alpha2,maxnbiter,idbb,rem,nbiter,mask,oldnb,newnb,bor,testici,de,xxp,yyp,rap,verif);

if (DEBUGLEVEL >=5 , print("entree dans bnfqfsolve"));
if (DEBUGLEVEL >= 3, print("(a,b)=(",aleg,",",bleg,")"));
  nf=bnf.nf;
  aleg = Mod(lift(aleg),nf.pol); aa = aleg;
  bleg = Mod(lift(bleg),nf.pol); bb = bleg;
\\
  if ( aa == 0,
if (DEBUGLEVEL >= 5, print("fin de bnfqfsolve"));
    return ([0,1,0]~));
  if ( bb == 0,
if (DEBUGLEVEL >= 5, print("fin de bnfqfsolve"));
    return ([0,0,1]~));
\\
  na = abs(norm(aa)); nb = abs(norm(bb));
  if( na > nb ,maxnb = na ,maxnb = nb);
  maxnb <<= 20;
  mat = Mod(matid(3),nf.pol); borne = 1;
  test = 0; nbiter = 0;
\\
  while (1,
    if(flag3 && bnf.clgp[1]>1, resl = bnfqfsolve2(bnf,aa,bb,auto)~;break);
if (DEBUGLEVEL >= 4, print("(na,nb,a,b)=",lift([na,nb,aa,bb,norm(aa),norm(bb)])));
if (DEBUGLEVEL >=5 , print("***",nb,"*** "));
    if ( nb >= maxnb ,
      mat=Mod(matid(3),nf.pol);
      aa = aleg; bb = bleg; na = abs(norm(aleg)); nb = abs(norm(bleg)));
    if ( aa == 1 ,resl = [1,1,0]~; break);
    if ( bb == 1 ,resl = [1,0,1]~; break);
    if ( aa+bb == 1 ,resl = [1,1,1]~; break);
    if ( aa+bb == 0 ,resl = [0,1,1]~; break);
    if ( aa == bb  && aa != 1 ,
      t = aa*mat[,1];
      mat[,1] = mat[,3]; mat[,3] = t;
      aa = -1; na = 1);
    if ( issquare(na) ,
      sq = nfissquare(nf,aa);
      if ( sq != [] ,resl = [sq[1],1,0]~; break));
    if ( issquare(nb) ,
      sq = nfissquare(nf,bb);
      if ( sq != [] ,resl = [sq[1],0,1]~; break));
    if ( na > nb ,
      t = aa; aa = bb; bb = t;
      t = na; na = nb; nb = t;
      t = mat[,3]; mat[,3] = mat[,2]; mat[,2] = t);
    if( nb == 1 ,
      if (DEBUGLEVEL >= 4, print("(a,b)=",lift([aa,bb])));
      if (DEBUGLEVEL >= 4, print("(na,nb)=",lift([na,nb])));
      if ( aleg == aa && bleg == bb ,mat = Mod(matid(3),nf.pol));
      if(flag3, resl = bnfqfsolve2(bnf,aa,bb,auto)~;break);
      pol = aa*x^2+bb;
      vecrat = ratpoint(nf,pol,borne++,myplist,0);
      if ( vecrat != 0 ,resl=[vecrat[2],vecrat[1],vecrat[3]]~; break);
  \\
      alpha = 0;
if (DEBUGLEVEL >= 4,     print("borne = ",borne));
      while ( alpha==0,
        xx = nfrandint(nf,borne); yy = nfrandint(nf,borne);
        borne++;
        alpha = xx^2-aa*yy^2 );
      bb *= alpha; nb *= abs(norm(alpha));
      t = xx*mat[,1]+yy*mat[,2];
      mat[,2] = xx*mat[,2]+aa*yy*mat[,1];
      mat[,1] = t;
      mat[,3] *= alpha;
    ,
      test = 1;
if (DEBUGLEVEL >= 4, print("on factorise bb=",bb));
      sun=bnfsunit(bnf,idealfactor(bnf,bb)[,1]~);
      fact=lift(bnfissunit(bnf,sun,bb));
if (DEBUGLEVEL >= 4, print("fact=",fact));
      suni=concat(bnf.futu,vector(length(sun[1]),i,nfbasistoalg(bnf,sun[1][i])));
      for(i=1,length(suni),
        if ( (f=fact[i]>>1) ,
          test =0;
          for(k=1,3,mat[k,3] /= suni[i]^f);
          nb /= abs(norm(suni[i]))^(2*f);
          bb/=suni[i]^(2*f)));
if (DEBUGLEVEL >= 4, print("on factorise bb=",bb));
      fact = idealfactor(nf,bb);
if (DEBUGLEVEL >= 4, print("fact=",fact));
      l = length(fact[,1]);
\\\\
      if ( test ,
        aux = 1;
        for ( i = 1 ,l,
	  if ( (f=fact[i,2]>>1) &&
	       !(fact[i,1][1]%2) && !psquarenf(nf,aa,fact[i,1]) ,
	    aux=idealmul(nf,aux,idealpow(nf,fact[i,1],f))));
        if ( aux != 1 ,
	  test = 0;
	  alpha = nfbasistoalg(nf,idealappr(nf,idealinv(nf,aux)));
          alpha2 = alpha^2;
          bb *= alpha2; nb *= abs(norm(alpha2));
          mat[,3] *= alpha));
      if ( test ,
	maxnbiter = 1<<l;
	sq = vector(l,i,nfissquarep(nf,aa,fact[i,1],fact[i,2]));
        l = length(sq);
if (DEBUGLEVEL >= 4, print("sq=",sq);print("fact=",fact);print("l=",l));
        if ( l > 1,
	  idbb = idealhnf(nf,bb);
	  rem = chinremain(nf,idbb,fact));
        test = 1; nbiter = 1;
        while (test && nbiter<=maxnbiter ,
	  if ( l > 1 ,
	    mask = nbiter; xx = 0;
	    for ( i = 1 ,l,
	      if ( mask%2 ,xx += rem[i]*sq[i] , xx -= rem[i]*sq[i] ); mask >>= 1)
          ,
            test = 0; xx = sq[1]);
          xx = mynfeltmod(nf,xx,bb);
          alpha = xx^2-aa;
          if ( alpha == 0, resl=[xx,1,0]~; break(2));
          t = alpha/bb;
if (DEBUGLEVEL >= 4, print("[alpha,bb]=",[alpha,bb]));
          oldnb = nb;
          newnb=abs(norm(t));
if (DEBUGLEVEL >= 4, print("[oldnb,newnb,oldnb/newnb]=",[oldnb,newnb,oldnb/newnb+0.]));
          while ( nb > newnb ,
            mat[,3] *= t;
            bb = t; nb = newnb;
            t = xx*mat[,1]+mat[,2];
            mat[,2] = aa*mat[,1] + xx*mat[,2];
            mat[,1] = t;
            xx = mynfeltmod(nf,-xx,bb);
            alpha = xx^2-aa;
            t = alpha/bb;
            newnb = abs(norm(t));
          );
        if ( nb == oldnb ,nbiter++ , test = 0);
        );
        if (nb == oldnb ,
          if(flag3,resl=bnfqfsolve2(bnf,aa,bb,auto)~;break);
          pol = aa*x^2+bb;
          vecrat =ratpoint(nf,pol,borne++<<1,myplist,0);
          if ( vecrat != 0 ,resl=[vecrat[2],vecrat[1],vecrat[3]]~; break);
\\
          bor=1000;yy=1;testici=1;
          for(i=1,10000,de=nfbasistoalg(nf,vectorv(poldegree(nf.pol),j,random(bor)));
            if(idealadd(bnf,de,bb)!=matid(poldegree(bnf.pol)),next);
            xxp=mynfeltmod(bnf,de*xx,bb);yyp=mynfeltmod(bnf,de*yy,bb);
            rap=(norm(xxp^2-aa*yyp^2)/nb^2+0.);
            if(abs(rap)<1,
if (DEBUGLEVEL >= 4, print("********** \n \n MIRACLE ",rap," \n \n ***"));
              t=(xxp^2-aa*yyp^2)/bb;
              mat[,3] *= t;
              bb = t; nb = abs(norm(bb));
if (DEBUGLEVEL >= 4, print("newnb = ",nb));
              t = xxp*mat[,1]+yyp*mat[,2];
              mat[,2] = aa*yyp*mat[,1] + xxp*mat[,2];
              mat[,1] = t;
              xx=xxp;yy=-yyp;testici=0;
              ));
\\
          if(testici,
            alpha = 0;
            while ( alpha == 0,
              xx = nfrandint(nf,4*borne); yy = nfrandint(nf,4*borne);
              borne++;
              alpha = xx^2-aa*yy^2);
            bb *= alpha; nb *= abs(norm(alpha));
            t = xx*mat[,1] + yy*mat[,2];
            mat[,2] = xx*mat[,2]+aa*yy*mat[,1];
            mat[,1] = t;
            mat[,3] *= alpha;)))));
  resl = lift(mat*resl);
if (DEBUGLEVEL >= 5, print("resl1=",resl));
if (DEBUGLEVEL >= 5, print("content=",content(resl)));
  resl /= content(resl);
  resl = Mod(lift(resl),nf.pol);
if (DEBUGLEVEL >=5, print("resl3=",resl));
  fact = idealadd(nf,idealadd(nf,resl[1],resl[2]),resl[3]);
  fact = bnfisprincipal(bnf,fact,3);
  resl *=1/nfbasistoalg(nf,fact[2]);
if (DEBUGLEVEL >= 5, print("resl4=",resl));
\\  if(resl[1] && resl[2] && resl[3] ,
\\if (DEBUGLEVEL >= 4, print("on essaye mordellreduce"));
\\    for(i=1,10,resl=mordellreduce(bnf,aleg,bleg,-1,resl));
\\if (DEBUGLEVEL >= 4, print("resl=",resl);print("on enleve les unites"));
\\    listfa = ideallistfactor(nf,resl);
\\    sun=bnfsunit(bnf,listfa);
\\    moy=round(sum(i=1,3,lift(bnfissunit(bnf,sun,resl[i])))/3);\\print("moy = ",moy);
\\    for(i=1,length(bnf.fu),resl/=bnf.fu[i]^moy[i]);
\\  );
if (DEBUGLEVEL >= 3, print("resl=",resl));
  verif = (resl[1]^2-aleg*resl[2]^2-bleg*resl[3]^2==0);
  if ( !verif && DEBUGLEVEL >= 0 , print("***** erreur dans bnfqfsolve : mauvaise solution *******",[a,b]);1/0);
if (DEBUGLEVEL >= 3, print("fin de bnfqfsolve"));return (resl);
}
if (DEBUGLEVEL >= 4, print("redquartique2 "));
{
redquartique2( bnf, pol, r,a,b)=
local(gcc,princ,den,ap,rp,pol2);
\\ reduction d'une quartique issue de la 2-descente
\\ en connaissant les valeurs de r, a et b.
if (DEBUGLEVEL >= 4, print("entree dans redquartique2"));
if (DEBUGLEVEL >= 4, print([r,a,b]));
if (DEBUGLEVEL >= 3, print(" reduction de la quartique ",pol));

  if (a==0,
    rp=0
  ,
    gcc=idealadd(bnf,b,a);
    if(gcc==1,
      rp=nfbasistoalg(bnf,idealaddtoone(bnf.nf,a,b)[1])/a;
      rp=mynfeltmod(bnf,r*rp,b)
    ,
      princ=bnfisprincipal(bnf,gcc,3);
      if (princ[1]==0,gcc=nfbasistoalg(bnf,princ[2])
      ,
if (DEBUGLEVEL >= 3, print(" quartique non reduite"));
if (DEBUGLEVEL >= 4, print("fin de redquartique2"));
        return([pol,0,1]));
      rp=nfbasistoalg(bnf,idealaddtoone(bnf.nf,a/gcc,b/gcc)[1])/(a/gcc);
      rp=mynfeltmod(bnf,r*rp,b)/gcc;
      b/=gcc;
    )
  );
 pol2=subst(pol/b,x,rp+b*x)/b^3;
if (DEBUGLEVEL >= 3, print(" quartique reduite = ",pol2));
if (DEBUGLEVEL >= 4, print("fin de redquartique2"));
 return([pol2,rp,b])
}
if (DEBUGLEVEL >= 4, print("redquartique22 "));
{
redquartique22( bnf, pol, r,a,b)=
\\ reduction d'une quartique issue de la 2-descente
\\ en connaissant les valeurs de r, a et b.
local(gcc,be,al,bl,princ,g,l,a1,rp,pol2);

if (DEBUGLEVEL >= 4, print("entree dans redquartique22"));
if (DEBUGLEVEL >= 4, print([r,a,b]));
if (DEBUGLEVEL >= 3, print(" reduction de la quartique ",pol));

  gcc=idealadd(bnf,b,a);
  be=idealaddtoone(bnf,idealdiv(bnf,a,gcc),idealdiv(bnf,b,gcc));
  al=nfbasistoalg(bnf,be[1]);bl=1-al;
  princ=bnfisprincipal(bnf,gcc,3);
  g=nfbasistoalg(bnf,princ[2]); print("g=",g);
  g*=idealdiv(bnf,gcc,g)[1,1];print("g=",g);
  l=al*g/a;a1=bl*g/b;print("l=",l);print("a1=",a1);
  rp=mynfeltmod(bnf,l*r,b);print("rp=",rp);print("b=",b);
  pol2=subst(subst(pol,x,x/g)*g^4,x,b*x+rp)/b^4;

if (DEBUGLEVEL >= 3, print(" quartique reduite = ",pol2));
if (DEBUGLEVEL >= 4, print("fin de redquartique22"));
  return([pol2,rp,b,g])
}
if (DEBUGLEVEL >= 4, print("redquartique3 "));
{
redquartique3( nf, pol, r,a,b)=
\\ reduction d'une quartique issue de la 2-descente
local(l,gcc,a0,r1,id1,lp,pol2,mat,a0p,u,v);

if (DEBUGLEVEL >= 4, print("entree dans redquartique3"));
if (DEBUGLEVEL >= 3, print(" reduction de la quartique ",pol));

  if(abs(norm(b))==1,return([pol/b^3,matid(2)]));
  print("r=",r);  print("a=",a);  print("b=",b);

  gcc=idealadd(nf,r,a);  print("gcc=",gcc);
  if(gcc==1,    l=1
  , print("1/gcc=",idealinv(nf,gcc));
    l=nfbasistoalg(nf,idealappr(nf,idealinv(nf,gcc)));print("l=",l);
    if(idealadd(nf,idealadd(nf,l*r,l*a),b)!=1,
      print("CAS No 1 NON ENCORE TRAITE DANS redquartique3");1/0));
  print("l=",l);
  a0=-nfbasistoalg(nf,round(nfalgtobasis(nf,a*l/b)));a=l*a+b*a0;print("new a = ",a);
  r1=-nfbasistoalg(nf,round(nfalgtobasis(nf,r*l/b)));r=l*r+b*r1;print("new r = ",r);


  gcc=idealadd(nf,a,b);
  if(gcc==1,
    print("CAS PARTICULIER DANS redquartique3");
    id1=idealaddtoone(nf,a,b);print("id1=",id1);
    lp=nfbasistoalg(nf,id1[1])/a;print("lp=",lp);
    r1=-nfbasistoalg(nf,round(nfalgtobasis(nf,r*lp/b)));r=lp*r+b*r1;print("new r = ",r);
    pol2=subst(pol,x,b*x+r)/b^3;print("pol2=",pol2);
    mat=[b,r;0,1];
    return([pol2,mat]);
  );

  gcc=idealadd(nf,r,a);
  print("gcc=",gcc);
  if(gcc==1,a0p=0
  , print("CAS No 2 NON ENCORE TRAITE DANS redquartique3");1/0);
  a=a+b*a0p;


  id1=idealaddtoone(nf,r,a);
  print("id1=",id1);
  v=nfbasistoalg(nf,id1[1])/r;
  u=-nfbasistoalg(nf,id1[2])/a;
  mat=[b*u,r;b*v,a];
  pol2=subst(pol,x,(b*u*x+r)/(b*v*x+a))*(b*v*x+a)^4/b^3;
if (DEBUGLEVEL >= 3, print(" quartique reduite = ",pol2));
if (DEBUGLEVEL >= 4, print("fin de redquartique3"));
  return([pol2,mat]);
}
if (DEBUGLEVEL >= 4, print("redquartique "));
{
redquartique( bnf, pol, list)=
\\ reduction d'une quartique de la forme y^2=ax^4+cx^2+dx+e
local(nf,pol2,multiplicateur,den,fact,i,j,cont,bnfisprinc,aux,factprinc,factprinc2,ordre,test,ordres);

if (DEBUGLEVEL >= 4, print("entree dans redquartique"));
if (DEBUGLEVEL >= 3, print("reduction de la quartique ",pol));

  nf = bnf.nf;
  pol2 = pol;
  multiplicateur = 1;
  den = deno(content(lift(pol2)));
  if ( den != 1 ,
      if (type(den) == "t_INT",fact = factor(den,0),fact = Mat([den,1]));
      den = 1;
      for (j=1,length(fact[,1]), den *= fact[j,1]^((fact[j,2]+1)>>1));
      pol2 *= den^2);
  cont=pollead(pol2);
  for ( i = 0, poldegree(pol2)-1 ,
    if ( polcoeff(pol2,i) ,
      cont = idealadd(nf,cont,polcoeff(pol2,i))));
if (DEBUGLEVEL >= 4, print("cont = ",cont));
  fact = idealfactor(nf,cont);
  for ( i = 1 ,length(fact[,1]),
    if ( fact[i,2]>1,
      bnfisprinc = bnfisprincipal(bnf,fact[i,1],3);
      if ( bnfisprinc[1] == 0,
        aux=nfbasistoalg(nf,bnfisprinc[2])^((fact[i,2]\2));
        pol2 /= aux^2;
        if (DEBUGLEVEL >= 4, print("on simplifie par un contenu ->",pol2)))));
  fact = ideallistfactor(nf,list);
  factprinc = vector(length(fact),i,0);
  factprinc2 = vector(length(fact),i,0);
  for ( j = 1 ,length(fact),\\print("j=",j,fact[j]);
    if ( fact[j][1] == nfbasistoalg(nf,fact[j][2]) ,
      factprinc[j] = fact[j][1]; factprinc2[j] = 1; next);
    aux = bnfisprincipal(bnf,fact[j],3);
    if ( aux[1] != 0 ,
      ordre=1;
      for( i = 1 ,length(aux[1]),
        ordre = lcm(ordre,bnf[8][1][2][i]/gcd(bnf[8][1][2][i],aux[1][i])));
      factprinc2[j] = ordre;
      factprinc[j] = nfbasistoalg(nf,bnfisprincipal(bnf,idealpow(nf,fact[j],ordre),3)[2]);
      next
    );
    factprinc[j] = nfbasistoalg(nf,aux[2]);
    factprinc2[j] = 1;
  );

  for( j = 1 ,length(fact),
	test = 0;
	ordres = factprinc2[j];
	while ( val(nf,polcoeff(pol2,0),fact[j])>= 2*ordres
               && val(nf,polcoeff(pol2,1),fact[j])>= 2*ordres
               && val(nf,polcoeff(pol2,2),fact[j])>= 2*ordres
               && val(nf,polcoeff(pol2,4),fact[j])>= 2*ordres ,
          test=1;pol2/=factprinc[j]^2);
        while ( val(nf,polcoeff(pol2,2),fact[j])>= 2*ordres
               && val(nf,polcoeff(pol2,4),fact[j])>=6*ordres ,
          test=1;pol2=subst(pol2,x,x/factprinc[j]^2)*factprinc[j]^2;
          multiplicateur/=factprinc[j]^2);
        while ( val(nf,polcoeff(pol2,1),fact[j])>=ordres
               && val(nf,polcoeff(pol2,2),fact[j])>= 2*ordres
               && val(nf,polcoeff(pol2,4),fact[j])>= 4*ordres ,
          test=1;pol2=subst(pol2,x,x/factprinc[j]);multiplicateur/=factprinc[j]);
        while ( val(nf,polcoeff(pol2,2),fact[j])>= 2*ordres
               && val(nf,polcoeff(pol2,1),fact[j])>= 3*ordres
               && val(nf,polcoeff(pol2,0),fact[j])>= 4*ordres ,
          test=1;pol2=subst(pol2,x,x*factprinc[j])/factprinc[j]^4;
          multiplicateur*=factprinc[j]);
        while ( val(nf,polcoeff(pol2,1),fact[j])>=ordres
               && val(nf,polcoeff(pol2,0),fact[j])>= 2*ordres ,
          test=1;pol2=subst(pol2,x,x*factprinc[j])/factprinc[j]^2;
          multiplicateur*=factprinc[j]);
        if (test, j--);
  );
if (DEBUGLEVEL >= 3, print("quartique reduite = ",pol2));
if (DEBUGLEVEL >= 4, print("multiplicateur = ",multiplicateur));
if (DEBUGLEVEL >= 4, print("fin de redquartique"));
  return ([pol2,multiplicateur]);
}
if (DEBUGLEVEL >= 4, print("ellranknf "));
{
ellranknf( bnf, ell, ext , help=[],bigflag=1,flag3=1,auto=[y])=
\\ bnf a un polynome en y.
\\ ell est donne par smallinitell.
\\ si ell= y^2=P(x), alors ext est
\\ ext[1] est une equation relative du corps (=P(x)),
\\ ext[2] est le resultat rnfequation(bnf,P,2);
\\ ext[3] est le buchinitfu (sur Q) de l'extension.
\\ dans la suite ext est note L = K(theta).
\\ help est une liste de points deja connus sur ell.
\\ si bigflag !=0 alors on applique redquartique.
\\ si flag3 ==1 alors on utilise bnfqfsolve2 (equation aux normes) pour resoudre Legendre
\\ auto est une liste d'automorphismes connus de bnf
\\ ca peut aider a factoriser certains discriminiants
\\ ell est de la forme y^2=x^3+A*x^2+B*x+C
\\ ie ell=[0,A,0,B,C], avec A,B et C entiers
\\
local(nf,unnf,ellnf,AAA,BBB,CCC,S,plist,Lrnf,SLprod,oddclass,SLlist,LS2gen,polrel,alpha,ttheta,KS2gen,LS2genunit,normLS2,normcoord,LS2coordtilda,LS2tilda,i,aux,j,myplist,listgen,listpoints,listpointstriv,listpointsmwr,list,m1,m2,loc,lastloc,maskwhile,iwhile,zc,iaux,liftzc,ispointtriv,point,c,b,a,sol,found,alphac,r,denc,dena,cp,alphacp,beta,mattr,vec,z1,ff,cont,d,e,polorig,pol,redq,transl,multip,UVW,pointxx,point2,v,rang);

if (DEBUGLEVEL >= 4, print("entree dans ellranknf"));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\      construction de L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  nf = bnf.nf;
  unnf = Mod(1,nf.pol);
  ellnf = ell*unnf;
  AAA = ellnf.a2;
if (DEBUGLEVEL >= 2, print("A= ",AAA));
  BBB = ellnf.a4;
if (DEBUGLEVEL >= 2, print("B= ",BBB));
  CCC = ellnf.a6;
if (DEBUGLEVEL >= 2, print("C= ",CCC));
  S = 6*lift(ellnf.disc);
  plist = idealfactor(nf,S)[,1];
  Lrnf = ext[3];
  SLprod = subst(lift(ext[1]'),y,lift(ext[2][2]));
if (DEBUGLEVEL >= 3, print(ext[2]));
  oddclass = 0;
  while (!oddclass,
\\ Constructoin de S:
    SLlist = idealfactor(Lrnf,SLprod)[,1]~;
\\ Construction des S-unites
    LS2gen = bnfsunit(Lrnf,SLlist);
if (DEBUGLEVEL >= 4, print("LS2gen=",LS2gen));
\\ on ajoute la partie paire du groupe de classes.
    oddclass = LS2gen[5][1]%2;
    if (!oddclass,
 if (DEBUGLEVEL >= 3, print("2-class group",LS2gen[5][3][1][1,1]));
      S*=LS2gen[5][3][1][1,1];
      SLprod=idealmul(Lrnf,SLprod,(LS2gen[5][3][1])));
  );

  polrel = ext[1];
  alpha = Mod(Mod(y,nf.pol),polrel); \\ alpha est l'element primitif de K
  ttheta = Mod(x,polrel);            \\ ttheta est la racine de P(x)

  KS2gen = bnfsunit(bnf,idealfactor(nf,S)[,1]~);

if (DEBUGLEVEL >= 3, print("#KS2gen = ",length(KS2gen[1])));
if (DEBUGLEVEL >= 3, print("KS2gen = ",KS2gen[1]));

  LS2genunit = lift(Lrnf.futu);
  LS2genunit = concat(LS2genunit,lift(nfbasistoalg(Lrnf,LS2gen[1])));

  LS2genunit = subst(LS2genunit,x,ttheta);
  LS2genunit = LS2genunit*Mod(1,polrel);
if (DEBUGLEVEL >= 3, print("#LS2genunit = ",length(LS2genunit)));
if (DEBUGLEVEL >= 3, print("LS2genunit = ",LS2genunit));

\\ dans LS2gen, on ne garde que ceux dont la norme est un carre.

  normLS2gen = norm(LS2genunit);
if (DEBUGLEVEL >= 4, print("normLS2gen=",normLS2gen));

\\ matrice de l'application norme

  normcoord = matrix(length(KS2gen[1])+length(bnf[8][5])+1,length(normLS2gen),i,j,0);
  for( i = 1 ,length(normLS2gen),
    normcoord[,i] = bnfissunit(bnf,KS2gen,normLS2gen[i]));
if (DEBUGLEVEL >= 4, print("normcoord=",normcoord));

 \\ construction du noyau de la norme

  LS2coordtilda = lift(matker(normcoord*Mod(1,2)));
if (DEBUGLEVEL >= 4, print("LS2coordtilda=",LS2coordtilda));
  LS2tilda = vector(length(LS2coordtilda[1,]),i,0);
  for( i = 1 ,length(LS2coordtilda[1,]),
    aux = 1;
    for ( j = 1 ,length(LS2coordtilda[,i]),
      if ( sign(LS2coordtilda[j,i]),
        aux *= LS2genunit[j]));
    LS2tilda[i] = aux;
  );

if (DEBUGLEVEL >= 3, print("LS2tilda = ",LS2tilda));
if (DEBUGLEVEL >= 3, print("norm(LS2tilda)= ",norm(LS2tilda)));

  myplist = makemyplist(nf,NBIDEAUX);
if (DEBUGLEVEL >= 4, print("myplist = ",myplist));

\\ Fin de la construction de LS2

  listgen = LS2tilda;
if (DEBUGLEVEL >= 2, print("LS2gen = ",listgen));
if (DEBUGLEVEL >= 2, print("#LS2gen = ",length(listgen)));
  listpoints = [];

if (DEBUGLEVEL >= 3, print("(A,B,C) = ",[AAA,BBB,CCC]));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\   Recherche de points triviaux   \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if (DEBUGLEVEL >= 2, print(" recherche de points triviaux sur la courbe "));
  listpointstriv = ratpoint(nf,x^3+AAA*x^2+BBB*x+CCC,LIMTRIV,myplist,1);
  listpointstriv = concat(help,listpointstriv);
if (DEBUGLEVEL >= 1, print("points triviaux sur la courbe = ",listpointstriv));

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\          parcours de L(S,2)         \\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  listpointsmwr = [];
  list = [ 6 , ellnf.disc,0 ];
  m1 = 0; m2 = 0; lastloc = -1;
  maskwhile = 1<<length(listgen);
  listELS=[0]; listnotELS=[];
  iwhile = 1;
  while ( iwhile < maskwhile,
if (DEBUGLEVEL >= 4, print("iwhile = ",iwhile);print("listgen = ",listgen));

\\ utilise la structure de groupe pour detecter une eventuelle solubilite locale.
if (DEBUGLEVEL >= 4, print("listELS = ",listELS);
                   print("listnotELS = ",listnotELS));
    sol=1;loc=0;
    for(i=1,length(listELS),
      for(j=1,length(listnotELS),
        if(bitxor(listELS[i],listnotELS[j])==iwhile,
if (DEBUGLEVEL >= 3, print(" Not ELS from group structure"));
          listnotELS=concat(listnotELS,[iwhile]);
          sol=0;break(2))));
    if(!sol, iwhile++; next);

    for(i=1,length(listELS),
      for(j=i+1,length(listELS),
        if(bitxor(listELS[i],listELS[j])==iwhile,
if (DEBUGLEVEL >= 3, print(" ELS from group structure"));
          listELS=concat(listELS,[iwhile]);
          loc=2;
          break(2))));

    iaux=vector(length(listgen),i,bittest(iwhile,i-1))~;
    iaux=(LS2coordtilda*iaux)%2;

    zc=unnf*prod(i=1,length(LS2genunit),LS2genunit[i]^iaux[i]);

if (DEBUGLEVEL >= 2, print("zc=",zc));
    liftzc = lift(zc);

\\ Est-ce un point trivial ?
    ispointtriv=0;
    for ( i = 1, length(listpointstriv),
      point = listpointstriv[i];
      if ( length(point)==2 || point[3] != 0,
        if ( nfissquare(Lrnf.nf,subst((lift(point[1])-x)*lift(liftzc),y,lift(ext[2][2]))),
if (DEBUGLEVEL >= 2, print(" vient du point trivial ",point));
          listpointsmwr=concat(listpointsmwr,[point]);
          m1++;
          listELS=concat(listELS,[iwhile]);
          if ( degre(iwhile) > lastloc ,m2++);
          sol = found = ispointtriv = 1; break
          )));

\\ Ce n'est pas un point trivial
    if ( !ispointtriv,
      c = polcoeff(liftzc,2);
      b = -polcoeff(liftzc,1);
      a = polcoeff(liftzc,0);
      sol = 0; found = 0;
\\ \\\\\\\\\\\\\
\\ On cherche a ecrire zc sous la forme a-b*theta
\\ \\\\\\\\\\\\\
      if ( c == 0, sol = 1,
        alphac = (AAA*b+BBB*c-a)*c+b^2;
if (DEBUGLEVEL >= 3, print("alphac = ",alphac));
        r = mysqrtnf(nf,norm(zc));
        if ( alphac == 0,
\\ cas particulier
if (DEBUGLEVEL >= 3, print(" on continue avec 1/zc "));
          sol = 1; zc = norm(zc)*(1/zc);
if (DEBUGLEVEL >= 2, print(" zc = ",zc))
        ,
\\ Il faut resoudre une forme quadratique
\\ Existence (locale = globale) d'une solution :
          denc = deno(lift(c));
          if ( denc != 1 ,cp=c*denc^2, cp=c);
          dena = deno(lift(alphac));
          if ( dena != 1 ,alphacp = alphac*dena^2, alphacp=alphac);
if (DEBUGLEVEL >= 2, print1(" symbole de Hilbert (",alphacp,",",cp,") = "));
          sol = loc || (mynfhilbert(nf, alphacp,cp)+1);
if (DEBUGLEVEL >= 2, print(sol-1));
          if (sol,
            beta = AAA*(AAA*b*c+BBB*c^2+b^2)-CCC*c^2+a*b;
            mattr = matid(3);
            mattr[1,1] = c ;mattr[2,2] = alphac ;
            mattr[3,3] = r ;mattr[2,3] = -beta;
            mattr[1,2] = -(b +AAA*c) ;mattr[1,3] = a-BBB*c+AAA*(AAA*c+b);
if (DEBUGLEVEL >= 2, print1(" sol de Legendre = "));
            vec = bnfqfsolve(bnf,alphacp,cp,myplist,flag3,auto);
if (DEBUGLEVEL >= 2, print(lift(vec)));
            aux = vec[2]*dena;
            vec[2] = vec[1];vec[1] = aux;
            vec[3] = vec[3]*denc;
            vec = (mattr^(-1))*vec;
            vec /= content(lift(vec));
            z1 = (vec[3]*ttheta+vec[2])*ttheta+vec[1];
if (DEBUGLEVEL >= 3, print(" z1 = ",z1));
            zc *= z1^2;
if (DEBUGLEVEL >= 2, print(" zc*z1^2 = ",zc));
          )
        )
      )
    );

    if (!sol,
      listnotELS=concat(listnotELS,[iwhile]);
      iwhile++;
      next
    );

\\ \\\\\\\\\\
\\ Maintenant zc est de la forme a-b*theta
\\ \\\\\\\\\\
    if (!ispointtriv,
      liftzc = lift(zc);
if (DEBUGLEVEL >= 3, print(" zc = ",liftzc));
if (DEBUGLEVEL >= 4, print(" N(zc) = ",norm(zc)));
      if ( poldegree(liftzc) < 2 ,, print("****** ERROR  c <> 0 *******");1/0);
      b = -polcoeff(liftzc,1);
      a = polcoeff(liftzc,0);
if (DEBUGLEVEL >= 4, print(" on factorise (a,b) = ",idealadd(nf,a,b)));
      ff = idealfactor(nf,idealadd(nf,a,b));
if (DEBUGLEVEL >= 4, print(" ff=",ff));
      cont = 1;
      for( i = 1, length(ff[,1]),
        cont = idealmul(nf,cont,idealpow(nf,ff[i,1],ff[i,2]\2)));
      cont = idealinv(nf,cont);
      cont = nfbasistoalg(nf,bnfisprincipal(bnf,cont,3)[2])^2;
      a *= cont; b *= cont; zc *= cont;
if (DEBUGLEVEL >= 4, print(" [a,b]=",[a,b]));
      if(nfissquare(nf,b),
if (DEBUGLEVEL >= 3, print("b est un carre"));
        point=[a/b,mysqrtnf(nf,subst(x^3+AAA*x^2+BBB*x+CCC,x,a/b))];
if (DEBUGLEVEL >= 2, print("point trouve = ",point));
        listpointsmwr=concat(listpointsmwr,[point]);
        m1++;
        if ( degre(iwhile) > lastloc ,m2++);
        found = (ispointtriv=1)
      )
    );

\\ \\\\\\\\\\\
\\ Construction de la quartique
\\ \\\\\\\\\\\
    if (!ispointtriv,
      r = mysqrtnf(nf,norm(zc));
if (DEBUGLEVEL >= 4, print(" r=",r));
      c = -2*(AAA*b+3*a);
if (DEBUGLEVEL >= 4, print(" c=",c));
      d = 8*r;
if (DEBUGLEVEL >= 4, print(" d=",d));
      e = (AAA^2*b^2 - 2*AAA*a*b-4*BBB*b^2-3*a^2);
if (DEBUGLEVEL >= 4, print(" e=",e));
      polorig = b*(x^4+c*x^2+d*x+e)*unnf;
if (DEBUGLEVEL >= 2, print(" quartique : (",lift(b),")*Y^2 = ",lift(polorig/b)));
      list[3] = b;
      pol = polorig;
      if ( bigflag ,
        redq = redquartique2(bnf,pol,r,a,b);
if (DEBUGLEVEL >= 2, print(" reduite: Y^2 = ",lift(redq[1])));
        pol = redq[1];transl=redq[2];multip=redq[3]
      );
      point = ratpoint(nf,pol,LIM1,myplist);
      if ( point != 0,
if (DEBUGLEVEL >= 2, print("point=",point));
        m1++;
        if ( bigflag ,
          point[1] = point[1]*multip+transl;
          point[2] = mysqrtnf(nf,subst(polorig,x,point[1]/point[3])));
        mattr = matid(3);
        mattr[1,1] = -2*b^2; mattr[1,2] = (AAA*b+a)*b;
        mattr[1,3] = a^2+(2*BBB-AAA^2)*b^2; mattr[2,2] = -b;
        mattr[2,3] = a+AAA*b; mattr[3,3] =r;
        UVW = [point[1]^2,point[3]^2,point[1]*point[3]]~;
        vec = (mattr^(-1))*UVW;
        z1 = (vec[3]*ttheta+vec[2])*ttheta+vec[1];
        zc *= z1^2;
        zc /= -polcoeff(lift(zc),1);
if (DEBUGLEVEL >= 3, print("zc*z1^2 = ",zc));
        pointxx = polcoeff(lift(zc),0);
        point2 = [ pointxx , mysqrtnf(nf,subst(x^3+AAA*x^2+BBB*x+CCC,x,pointxx))];
if (DEBUGLEVEL >= 1, print(" point trouve = ",point2));
        listpointsmwr = concat(listpointsmwr,[point2]);
        if ( degre(iwhile) > lastloc ,m2++);
        found = 1; lastloc = -1;
      ,
        if ( loc
             || (!bigflag && locallysoluble(nf,pol,r,a,b))
             || (bigflag && locallysoluble(nf,pol,0,1,1)),

          if (!loc, listlistELS=concat(listELS,[iwhile]));
          if ( degre(iwhile) > lastloc ,m2++; lastloc=degre(iwhile));
          point = ratpoint(nf,pol,LIM3,myplist);
          if ( point!=0,
            m1++;
            if ( bigflag ,
              point[1] = point[1]*multip+transl;
              point[2] = mysqrtnf(nf,subst(polorig,x,point[1]/point[3]));
            );
            mattr = matid(3);
            mattr[1,1] = -2*b^2; mattr[1,2] = (AAA*b+a)*b;
            mattr[1,3] = a^2+(2*BBB-AAA^2)*b^2; mattr[2,2] = -b;
            mattr[2,3] = a+AAA*b; mattr[3,3] =r;
            UVW = [point[1]^2,point[3]^2,point[1]*point[3]]~;
            vec = (mattr^(-1))*UVW;
            z1 = (vec[3]*ttheta+vec[2])*ttheta+vec[1];
            zc *= z1^2;
            zc /= -polcoeff(lift(zc),1);
if (DEBUGLEVEL >= 3, print(" zc*z1^2 = ",zc));
            pointxx = polcoeff(lift(zc),0);
            point2 = [ pointxx , mysqrtnf(nf,subst(x^3+AAA*x^2+BBB*x+CCC,x,pointxx))];
if (DEBUGLEVEL >= 1, print(" point trouve = ",point2));
            listpointsmwr=concat(listpointsmwr,[point2]);
            found = 1;lastloc = -1;
          )
        , listnotELS=concat(listnotELS,[iwhile])
        )
      )
    );
    if ( found,
      found = 0; lastloc = -1;
      v = degre(iwhile);
      iwhile = 1<<v;
      maskwhile >>= 1;
      LS2coordtilda = vecextract(LS2coordtilda,1<<length(listgen)-iwhile-1);
      listgen = vecextract(listgen,1<<length(listgen)-iwhile-1);
      while( listELS[length(listELS)] >= iwhile,
        listELS = vecextract(listELS,1<<(length(listELS)-1)-1));
      while( length(listnotELS) && listnotELS[length(listnotELS)] >= iwhile,
        listnotELS = vecextract(listnotELS,1<<(length(listnotELS)-1)-1));
    , iwhile ++
    )
  );

if (DEBUGLEVEL >= 2,
  print("m1 = ",m1);
  print("m2 = ",m2));
if (DEBUGLEVEL >= 1,
  print("#S(E/K)[2]    = ",1<<m2));
  if (m1==m2,
if (DEBUGLEVEL >= 1,
    print("#E(K)/2E(K)   = ",1<<m1);
    print("#III(E/K)[2]  = 1");
    print("rang(E/K)     = ",m1));
    rang = m1;
  ,
if (DEBUGLEVEL >= 1,
    print("#E(K)/2E(K)  >= ",1<<m1);
    print("#III(E/K)[2] <= ",1<<(m2-m1));
    print("rang(E/K)    >= ",m1));
    rang = m1;
    if ((m2-m1)%2,
if (DEBUGLEVEL >= 1,
      print(" III devrait etre un carre, donc ");
      if(m2-m1>1,
        print("#E(K)/2E(K)  >= ",1<<(m1+1));
        print("#III(E/K)[2] <= ",1<<(m2-m1-1));
        print("rang(E/K)    >= ",m1+1);
      ,
        print("#E(K)/2E(K)  = ",1<<(m1+1));
        print("#III(E/K)[2] = 1");
        print("rang(E/K)    = ",m1+1)));
      rang = m1+1)
  );
if (DEBUGLEVEL >= 1, print("listpointsmwr = ",listpointsmwr));
  for ( i = 1 ,length(listpointsmwr) ,
    if ( !ellisoncurve(ellnf,listpointsmwr[i]) ,
      print("****** MAUVAIS POINT = ",listpointsmwr[i])));
if (DEBUGLEVEL >= 4, print("fin de ellranknf"));
  return ([rang,m2,listpointsmwr]);
}
if (DEBUGLEVEL >= 4, print("main"));
{
main(bnf,Amain,Bmain,Cmain,help=[],bigflag=1,flag3=1)=
\\ attention bnf a un polynome en y.
\\ si bigflag !=0, on reduit les quartiques
\\ si flag3 != 0, on utilise bnfqfsolve2
local(eqtheta,rnfeq,bbnf,ellmain,extmain,rang);

if (DEBUGLEVEL >= 3, print("entree dans main"));
  eqtheta = x^3+Amain*x^2+Bmain*x+Cmain;
if (DEBUGLEVEL >= 1, print("courbe elliptique : Y^2 = ",eqtheta));
  rnfeq = rnfequation(bnf,eqtheta,1);
if (DEBUGLEVEL >= 3, print1("bbnfinit "));
  bbnf = bnfinit(rnfeq[1],1);
if (DEBUGLEVEL >= 3, print("done"));
  ellmain = ellinit([0,Amain,0,Bmain,Cmain])*Mod(1,bnf.pol);
  extmain = [eqtheta, rnfeq , bbnf];
  rang = ellranknf(bnf,ellmain, extmain,help,bigflag,flag3);
if (DEBUGLEVEL >= 3, print("fin de main"));
  return (rang);
}
if (DEBUGLEVEL >= 4, print("mordellreduce "));
{
mordellreduce(bnf,amor,bmor,cmor=-1,vecx)=
local(k,nf,x0,y0,z0,pgcd,vectone,v,u,w,aux1,aux2,solution,fact);

  k=2*cmor;nf=bnf.nf;
  x0=vecx[2];y0=vecx[3];z0=vecx[1];pgcd=idealadd(nf,x0,y0);
  vectone=idealaddtoone(nf,idealdiv(nf,x0,pgcd),idealdiv(nf,y0,pgcd));
  v=-cmor*nfbasistoalg(nf,vectone[1])/x0;u=cmor*nfbasistoalg(nf,vectone[2])/y0;
  w=-(amor*u*x0+bmor*v*y0)/(cmor*z0);
  w=nfbasistoalg(nf,round(nfalgtobasis(nf,w)));
  aux1=amor*u^2+bmor*v^2+cmor*w^2;
  aux2=2*(amor*u*x0+bmor*v*y0+cmor*w*z0);
  solution=[z0*aux1-w*aux2,x0*aux1-u*aux2,y0*aux1-v*aux2]~/(k);
  solution/=content(lift(solution));
  fact = idealadd(nf,idealadd(nf,solution[1],solution[2]),solution[3]);
  fact = bnfisprincipal(bnf,fact,3);
  if ( fact[1] == 0,
    solution = solution*(1/nfbasistoalg(nf,fact[2])));
if (DEBUGLEVEL >= 4, print(amor*solution[2]^2+bmor*solution[3]^2+cmor*solution[1]^2));
  return (solution);
}
if (DEBUGLEVEL >= 4, print("complete"));
{
complete(bnf,e1,e2,e3,flag3=1,auto=[y])=
\\ calcul du rang d'une courbe elliptique
\\ par la methode de 2-descente complete.
\\ Y^2 = (x-e1)*(x-e2)*(x-e3);
\\ en suivant la methode decrite par J.Silverman
\\ si flag3 ==1 alors on utilise bnfqfsolve2 (equation aux normes)
\\ pour resoudre Legendre
\\ on pourra alors utiliser auto=nfgaloisconj(bnf.pol)

\\ e1, e2 et e3 sont des entiers algebriques de bnf.

local(KS2prod,oddclass,KS2gen,i,myplist,vect,selmer,rang,X,Y,b1,b2,vec,z1,z2,d31,quart0,quart,cont,fa,point,solx,soly,listepoints);

\\ calcul de K(S,2)

  KS2prod = (e1-e2)*(e2-e3)*(e3-e1)*2;
  oddclass = 0;
  while ( !oddclass ,
    KS2gen = bnfsunit(bnf,idealfactor(bnf,KS2prod)[,1]~);
    oddclass = (KS2gen[5][1]%2);
    if (!oddclass,
      KS2prod = idealmul(bnf,KS2prod,(KS2gen[5][3][1])));
  );
  KS2gen = KS2gen[1];
  for( i = 1 ,length(KS2gen),
    KS2gen[i] = nfbasistoalg(bnf, KS2gen[i]));
  KS2gen = concat(Mod(lift(bnf.tufu),bnf.pol),KS2gen);
if (DEBUGLEVEL >= 2,
  print("#K(S,2)gen          = ",length(KS2gen));
  print("K(S,2)gen=",KS2gen));

  myplist = makemyplist(bnf.nf,NBIDEAUX);
if (DEBUGLEVEL >= 4, print("myplist = ",myplist));

\\ parcours de K(S,2)*K(S,2)

  vect=vector(length(KS2gen),i,[0,1]);
  selmer = 0;
  rang = 0;
  listepoints = [];

  forvec(X = vect ,
    b1 = prod(i=1,length(KS2gen),KS2gen[i]^X[i]);
    forvec(Y = vect ,
      b2 = prod(i=1,length(KS2gen),KS2gen[i]^Y[i]);

if (DEBUGLEVEL >= 3, print("[b1,b2]=",lift([b1,b2])));

\\ points triviaux provenant de la 2-torsion

      if ( b1==1 && b2==1 ,
if (DEBUGLEVEL >= 2, print("point trivial infty"));
        selmer++;rang++;next);
      if ( nfissquare(bnf.nf,(e2-e1)*b1)
        && nfissquare(bnf.nf,(e2-e3)*(e2-e1)*b2) ,
if (DEBUGLEVEL >= 2, print(" point trivial [e2,0]"));
        selmer++;rang++;next);
      if ( nfissquare(bnf.nf,(e1-e2)*b2)
        && nfissquare(bnf.nf,(e1-e3)*(e1-e2)*b1) ,
if (DEBUGLEVEL >= 2, print(" point trivial [e1,0]"));
        selmer++;rang++;next);
      if ( nfissquare(bnf.nf,(e3-e1)*b1)
        && nfissquare(bnf.nf,(e3-e2)*b2) ,
if (DEBUGLEVEL >= 2, print(" point trivial [e3,0]"));
        selmer++;rang++;next);

\\ premier critere local : sur les formes quadratiques

      if(mynfhilbert(bnf.nf,b1*b2,b1*(e2-e1)) < 0
        || mynfhilbert(bnf.nf,b2,b1*(e3-e1)) < 0
        || mynfhilbert(bnf.nf,b1,b2*(e3-e2)) < 0
        ,
if (DEBUGLEVEL >= 3, print("not ELS"));
        next);

if (DEBUGLEVEL >= 2, print("[b1,b2]=",lift([b1,b2])));
if (DEBUGLEVEL >= 2, print("qf loc soluble"));

\\ solution de la premiere forme quadratique

      if(b1!=b2,
        vec = bnfqfsolve(bnf,b1*b2,b1*(e2-e1),myplist,flag3);
if (DEBUGLEVEL >= 3, print("sol part = ",vec));
        if(vec[3]==0, print(" BUG !!! : vec[3]=0 "));
        z1=vec[1]/vec[3]/b1;z2=vec[2]/vec[3]
      ,
        z1=(1+(e2-e1)/b1)/2;z2=z1-1;
      );
      d31=e3-e1;
      quart0=b2^2*(z1^2*b1 - d31)*x^4 - 4*z1*b2^2*z2*b1*x^3
           + 2*b1*b2*(z1^2*b1 + 2*b2*z2^2 + d31)*x^2 - 4*z1*b2*z2*b1^2*x
           + b1^2*(z1^2*b1 - d31);
      quart=quart0*b1*b2;
if (DEBUGLEVEL >= 4, print("quart=",quart));
      quart*=denominator(simplify(content(quart)))^2;
      cont=simplify(content(lift(quart)));
      fa=factor(cont);
      for(i=1,length(fa[,1]),quart/=fa[i,1]^(2*(fa[i,2]\2)));
if (DEBUGLEVEL >= 3, print("quart red=",quart));

\\ la quartique est-elle localement soluble ?

      if ( !locallysoluble(bnf.nf,quart) ,
if (DEBUGLEVEL >= 2, print (" quartic not ELS "));
        next);
      selmer++;

\\ recherche de points sur la quartique.

      point = ratpoint(bnf.nf,quart,LIM3,myplist);
      if ( point ,
if (DEBUGLEVEL >= 2, print("point trouve sur la quartique !!"));
if (DEBUGLEVEL >= 3, print(point));
        if( point[3] ,
          point/=point[3];
          z1=(2*b2*point[1]*z2-z1*(b1+b2*point[1]^2))/(b1-b2*point[1]^2);
          solx=b1*z1^2+e1;
          soly=mysqrtnf(bnf.nf,(solx-e1)*(solx-e2)*(solx-e3));
          listepoints=concat(listepoints,[[solx,soly]]);
if (DEBUGLEVEL >= 1, print("point sur la courbe elliptique =",[solx,soly]))
        );
        rang++
      ,
if (DEBUGLEVEL >= 2, print("aucun point trouve sur la quartique"))
      )
    )
  );

\\ fin

if (DEBUGLEVEL >= 1,
  print("# Selmer    = ",selmer);
  print("#E[K]/2E[K] = ",rang);
  print("#E[2]       = 4");
);
  rang = ceil(log(rang)/log(2))-2;
if (DEBUGLEVEL >= 1,
  print("rang        = ",rang);
  if (rang , print("points = ",listepoints));
);
return(rang);
}


