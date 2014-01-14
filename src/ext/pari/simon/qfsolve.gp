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
\\ adresse du fichier:
\\ www.math.unicaen.fr/~simon/qfsolve.gp
\\
\\  *********************************************
\\  *          VERSION 02/10/2009               *
\\  *********************************************
\\
\\ Programme de resolution des equations quadratiques
\\ langage: GP
\\ pour l'utiliser, lancer gp, puis taper
\\ \r qfsolve.gp
\\
\\ Ce fichier contient 4 principales fonctions:
\\
\\ - Qfsolve(G,factD): pour resoudre l'equation quadratique X^t*G*X = 0
\\ G doit etre une matrice symetrique n*n, a coefficients dans Z.
\\ S'il n'existe pas de solution, la reponse est un entier
\\ indiquant un corps local dans lequel aucune solution n'existe
\\ (-1 pour les reels, p pour Q_p).
\\ Si on connait la factorisation de -abs(2*matdet(G)),
\\ on peut la passer par le parametre factD pour gagner du temps.
\\
\\ - Qfparam(G,sol,fl): pour parametrer les solutions de la forme
\\ quadratique ternaire G, en utilisant la solution particuliere sol.
\\ si fl>0, la 'fl'eme forme quadratique est reduite.
\\
\\ - IndefiniteLLL(G,c): pour resoudre ou reduire la forme quadratique
\\ G a coefficients entiers. Il s'agit d'un algorithme type LLL, avec la
\\ constante 1/4<c<1.
\\
\\ - class2(d,factd): determine le 2-Sylow du (narrow) groupe de classes de
\\ discriminant d, ou d est un discriminant fondamental.
\\ Si on connait la factorisation de abs(2*d),
\\ on peut la donner dans factd, et dans ce cas le reste
\\ de l'algorithme est polynomial.
\\

\\
DEBUGLEVEL_qfsolve = 0;

\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\  DEBUT DES SOURCES                \\
\\
\\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
{QfbReduce(M) =
\\ M = [a,b;b;c] est a coefficients entiers.
\\ Reduction de la forme quadratique binaire
\\   qf = (a,b,c)=a*X^2+2*b*X*Y+c*Y^2
\\ Renvoit la matrice de reduction de det = +1.

local(a,b,c,H,test,di,q,r,nexta,nextb,nextc,aux);

if( DEBUGLEVEL_qfsolve >= 5, print("entree dans QfbReduce avec ",M));

  a = M[1,1]; b = M[1,2]; c = M[2,2];

  H = matid(2); test = 1;
  while( test && a,
    di = divrem(b,a); q = di[1]; r = di[2];
    if( 2*r > abs(a), r -= abs(a); q += sign(a));
    H[,2] -= q*H[,1];
    nextc = a; nextb = -r; nexta= (nextb-b)*q+c;

    if( test = abs(nexta) < abs(a),
      c = nextc; b = nextb; a = nexta;
      aux = H[,1]; H[,1] = -H[,2]; H[,2] = aux
    )
  );

if( DEBUGLEVEL_qfsolve >= 5, print("sortie de QfbReduce avec ",H));
return(H);
}
{IndefiniteLLL(G,c=1,base=0) =
\\ Performs first a LLL reduction on a positive definite
\\ quadratic form QD bounding the indefinite G.
\\ Then finishes the reduction with IndefiniteLLL2.

local(n,M,QD,M1,S,red);

  n = length(G);
  M = matid(n);
  QD = G;
  for( i = 1, n-1,
    if( !QD[i,i],
return(IndefiniteLLL2(G,c,base))
    );
    M1 = matid(n);
    M1[i,] = -QD[i,]/QD[i,i];
    M1[i,i] = 1;
    M = M*M1;
    QD = M1~*QD*M1
  );
  M = M^(-1);
  QD = M~*abs(QD)*M;
  S = qflllgram(QD/content(QD));
  red = IndefiniteLLL2(S~*G*S,c,base);
  if( type(red) == "t_COL",
return(S*red));
  if( length(red) == 3,
return([red[1],S*red[2],S*red[3]]));
return([red[1],S*red[2]]);
}
{IndefiniteLLL2(G,c=1,base=0) =
\\ following Cohen's book p. 86
\\ but without b and bstar: works on G
\\ returns [H~*G*H,H] where det(H) = 1 and H~*G*H is reduced.
\\ Exit with a norm 0 vector if one such is found.
\\ If base == 1 and norm 0 is obtained, returns [H~*G*H,H,sol] where
\\   sol is a norm 0 vector and is the 1st column of H.

local(n,H,M,A,aux,sol,k,nextk,swap,q,di,HM,aux1,aux2,Mkk1,bk1new,Mkk1new,newG);

  n = length(G);
if( DEBUGLEVEL_qfsolve >= 3, print("LLL dim ",n," avec G=",log(vecmax(abs(G)))/log(10)));
if( DEBUGLEVEL_qfsolve >= 4, print("LLL avec ");printp(G));

 if( n <= 1, return([G,matid(n)]));

  H = M = matid(n); A = matrix(n,n);

\\ compute Gram-Schmidt

  for( i = 1, n,
    if( !(A[i,i] = G[i,i]),
      if( base,
        aux = H[,1]; H[,1] = H[,i]; H[,i] = -aux;
        return([H~*G*H,H,H[,1]])
      , return(M[,i])));
    for( j = 1, i-1,
      A[i,j] = G[i,j] - sum( k = 1, j-1, M[j,k]*A[i,k]);
      M[i,j] = A[i,j]/A[j,j];
      A[i,i] -= M[i,j]*A[i,j];
      if( !A[i,i],
        sol = (M^(-1))~[,i]; sol /= content(sol);
        if( base,
          H = completebasis(sol);
          aux = H[,1]; H[,1] = H[,n]; H[,n]= -aux;
          return([H~*G*H,H,H[,1]])
        , return(sol)))
    )
  );

\\ LLL loop

  k = 2; nextk = 1;
  while( k <= n,

    swap = 1;
    while( swap,
      swap = 0;

\\ red(k,k-1);
      if( q = round(M[k,k-1]),
        for( i = 1, k-2, M[k,i] -= q*M[k-1,i]);
        M[k,k-1] -= q;
        for( i = 1, n,
          A[k,i] -= q*A[k-1,i];
          H[i,k] -= q*H[i,k-1]
        )
      );

\\ preparation du swap(k,k-1)

      if( issquare( di = -A[k-1,k-1]*A[k,k]),
\\ di est le determinant de matr
\\ On trouve une solution
        HM = (M^(-1))~;
        aux1 = sqrtint(numerator(di));
        aux2 = sqrtint(denominator(di));
        sol = aux1*HM[,k-1]+aux2*A[k-1,k-1]*HM[,k];
        sol /= content(sol);
        if( base,
          H = H*completebasis(sol,1);
          aux = H[,1]; H[,1] = H[,n]; H[,n] = -aux;
          return([H~*G*H,H,H[,1]])
        , return(H*sol)
        )
      );

\\ On reduit [k,k-1].
      Mkk1 = M[k,k-1];
      bk1new = Mkk1^2*A[k-1,k-1] + A[k,k];
      if( swap = abs(bk1new) < c*abs(A[k-1,k-1]),
        Mkk1new = -Mkk1*A[k-1,k-1]/bk1new
      );

\\ Calcul des nouvelles matrices apres le swap.
      if( swap,
        for( j = 1, n,
          aux = H[j,k-1]; H[j,k-1] = H[j,k]; H[j,k] = -aux);
        for( j = 1, k-2,
          aux = M[k-1,j]; M[k-1,j] = M[k,j]; M[k,j] = -aux);
        for( j = k+1, n,
          aux = M[j,k]; M[j,k] = -M[j,k-1]+Mkk1*aux; M[j,k-1] = aux+Mkk1new*M[j,k]);
        for( j = 1, n, if( j != k && j != k-1,
          aux = A[k-1,j]; A[k-1,j] = A[k,j]; A[k,j] =- aux;
          aux = A[j,k-1]; A[j,k-1] = Mkk1*aux+A[j,k]; A[j,k] = -aux-Mkk1new*A[j,k-1]));

        aux1 = A[k-1,k-1];
        aux2 = A[k,k-1];
        A[k,k-1]  =-A[k-1,k] - Mkk1*aux1;
        A[k-1,k-1]= A[k,k]   + Mkk1*aux2;
        A[k,k]    = aux1     - Mkk1new*A[k,k-1];
        A[k-1,k]  =-aux2     - Mkk1new*A[k-1,k-1];

        M[k,k-1] = Mkk1new;

if( DEBUGLEVEL_qfsolve >=4, newG=H~*G*H;print(vector(n,i,matdet(vecextract(newG,1<<i-1,1<<i-1)))));

        if( k != 2, k--)
      )
    );

    forstep( l = k-2, 1, -1,
\\ red(k,l)
      if( q = round(M[k,l]),
        for( i = 1, l-1, M[k,i] -= q*M[l,i]);
        M[k,l] -= q;
        for( i = 1, n,
          A[k,i] -= q*A[l,i];
          H[i,k] -= q*H[i,l]
        )
      )
    );
    k++
  );
return([H~*G*H,H]);
}
{kermodp(M,p) =
\\ Compute the kernel of M mod p.
\\ returns [d,U], where
\\ d = dim (ker M mod p)
\\ U in SLn(Z), and its first d columns span the kernel.

local(n,U,d);

  n = length(M);
  U = centerlift(matker(M*Mod(1,p)));
  d = length(U);
  U = completebasis(U);
  U = matrix(n,n,i,j,U[i,n+1-j]);
return([d,U]);
}
{Qfparam(G,sol,fl=3) =
\\ G est une matrice symetrique 3*3, et sol une solution de sol~*G*sol=0.
\\ Renvoit une parametrisation des solutions avec de bons invariants,
\\ sous la forme d'une matrice 3*3, dont chaque ligne contient
\\ les coefficients de chacune des 3 formes quadratiques.

\\ si fl!=0, la 'fl'eme forme quadratique est reduite.

local(U,G1,G2);

if( DEBUGLEVEL_qfsolve >= 5, print("entree dans Qfparam"));
  sol /= content(sol);
\\ construction de U telle que U[,3] est sol, et det(U) = +-1
  U = completebasis(sol,1);
  G1 = U~*G*U; \\ G1 a un 0 en bas a droite.
  G2 = [-2*G1[1,3],-2*G1[2,3],0;
      0,-2*G1[1,3],-2*G1[2,3];
      G1[1,1],2*G1[1,2],G1[2,2]];
  sol = U*G2;
  if(fl && !issquare(matdet( U = [sol[fl,1],sol[fl,2]/2;
                                  sol[fl,2]/2,sol[fl,3]])),
    U = QfbReduce(U);
    U = [U[1,1]^2,2*U[1,2]*U[1,1],U[1,2]^2;
         U[2,1]*U[1,1],U[2,2]*U[1,1]+U[2,1]*U[1,2],U[1,2]*U[2,2];
         U[2,1]^2,2*U[2,1]*U[2,2],U[2,2]^2];
    sol = sol*U
  );
if( DEBUGLEVEL_qfsolve >= 5, print("sortie de Qfparam"));
return(sol);
}
{LLLgoon3(G,c=1) =
\\ reduction LLL de la forme quadratique G (matrice de Gram)
\\ en dim 3 seulement avec detG = -1 et sign(G) = [1,2];

local(red,U1,G2,bez,U2,G3,cc,U3);

  red = IndefiniteLLL(G,c,1);
\\ On trouve toujours un vecteur isotrope.
  U1 = [0,0,1;0,1,0;1,0,0];
  G2 = U1~*red[1]*U1;
\\ la matrice G2 a un 0 dans le coin en bas a droite.
  bez = bezout(G2[3,1],G2[3,2]);
  U2 = [bez[1],G2[3,2]/bez[3],0;bez[2],-G2[3,1]/bez[3],0;0,0,-1];
  G3 = U2~*G2*U2;
\\ la matrice G3 a des 0 sous l'anti-diagonale.
  cc = G3[1,1]%2;
  U3 = [1,0,0;  cc,1,0;
        round(-(G3[1,1]+cc*(2*G3[1,2]+G3[2,2]*cc))/2/G3[1,3]),
        round(-(G3[1,2]+cc*G3[2,2])/G3[1,3]),1];
return([U3~*G3*U3,red[2]*U1*U2*U3]);
}
{completebasis(v,redflag=0) =
\\ Donne une matrice unimodulaire dont la derniere colonne est v.
\\ Si redflag <> 0, alors en plus on reduit le reste.

local(U,n,re);

  U = (mathnf(Mat(v~),1)[2]~)^-1;
  n = length(v~);
  if( n==1 || !redflag, return(U));
  re = qflll(vecextract(U,1<<n-1,1<<(n-1)-1));
return( U*matdiagonalblock([re,Mat(1)]));
}
{LLLgoon(G,c=1) =
\\ reduction LLL de la forme quadratique G (matrice de Gram)
\\ ou l'on continue meme si on trouve un vecteur isotrope

local(red,U1,G2,U2,G3,n,U3,G4,U,V,B,U4,G5,U5,G6);

  red = IndefiniteLLL(G,c,1);
\\ si on ne trouve pas de vecteur isotrope, il n'y a rien a faire.
  if( length(red) == 2, return(red));
\\ sinon :
  U1 = red[2];
  G2 = red[1]; \\ On a G2[1,1] = 0
  U2 = mathnf(Mat(G2[1,]),4)[2];
  G3 = U2~*G2*U2;
\\ la matrice de G3 a des 0 sur toute la 1ere ligne,
\\ sauf un 'g' au bout a droite, avec g^2| det G.
  n = length(G);
  U3 = matid(n); U3[1,n] = round(-G3[n,n]/G3[1,n]/2);
  G4 = U3~*G3*U3;
\\ le coeff G4[n,n] est reduit modulo 2g
  U = vecextract(G4,[1,n],[1,n]);
  if( n == 2,
    V = matrix(2,0)
  , V = vecextract(G4,[1,n],1<<(n-1)-2));
  B = round(-U^-1*V);
  U4 = matid(n);
  for( j = 2, n-1,
    U4[1,j] = B[1,j-1];
    U4[n,j] = B[2,j-1]
  );
  G5 = U4~*G4*U4;
\\ la derniere colonne de G5 est reduite
  if( n < 4, return([G5,U1*U2*U3*U4]));

  red = LLLgoon(matrix(n-2,n-2,i,j,G5[i+1,j+1]),c);
  U5 = matdiagonalblock([Mat(1),red[2],Mat(1)]);
  G6 = U5~*G5*U5;
return([G6,U1*U2*U3*U4*U5]);
}
{QfWittinvariant(G,p) =
\\ calcule l'invariant c (=invariant de Witt) d'une forme quadratique,
\\ p-adique (reel si p = -1)

local(n,det,diag,c);

  n = length(G);
\\ On diagonalise d'abord G
  det  = vector( n+1, i, matdet(matrix(i-1,i-1,j,k,G[j,k])));
  diag = vector( n, i, det[i+1]/det[i]);

\\ puis on calcule les symboles de Hilbert
  c = prod( i = 1, n,
        prod( j = i+1, n,
          hilbert( diag[i], diag[j], p)));
return(c);
}
{Qflisteinvariants(G,fa=[]) =
\\ G est une forme quadratique, ou une matrice symetrique,
\\ ou un vecteur contenant des formes quadratiques de meme discriminant.
\\ fa = factor(-abs(2*matdet(G)))[,1] est facultatif.

local(l,sol,n,det);

if( DEBUGLEVEL_qfsolve >= 4, print("entree dans Qflisteinvariants",G));
  if( type(G) != "t_VEC", G = [G]);
  l = length(G);
  for( j = 1, l,
    if( type(G[j]) == "t_QFI" || type(G[j]) == "t_QFR",
      G[j] = mymat(G[j])));

  if( !length(fa),
    fa = factor(-abs(2*matdet(G[1])))[,1]);

  if( length(G[1]) == 2,
\\ En dimension 2, chaque invariant est un unique symbole de Hilbert.
    det = -matdet(G[1]);
    sol = matrix(length(fa),l,i,j,hilbert(G[j][1,1],det,fa[i])<0);
if( DEBUGLEVEL_qfsolve >= 4, print("sortie de Qflisteinvariants"));
    return([fa,sol])
  );

  sol = matrix(length(fa),l);
  for( j = 1, l,
    n = length(G[j]);
\\ En dimension n, on calcule un produit de n symboles de Hilbert.
    det = vector(n+1, i, matdet(matrix(i-1,i-1,k,m,G[j][k,m])));
    for( i = 1, length(fa),
      sol[i,j] = prod( k = 1, n-1, hilbert(-det[k],det[k+1],fa[i]))*hilbert(det[n],det[n+1],fa[i]) < 0;
    )
  );
if( DEBUGLEVEL_qfsolve >= 4, print("sortie de Qflisteinvariants"));
return([fa,sol]);
}
{Qfsolvemodp(G,p) =
\\ p a prime number.
\\ finds a solution mod p for the quadatic form G
\\ such that det(G) !=0 mod p and dim G = n>=3;
local(n,vdet,G2,sol,x1,x2,x3,N1,N2,N3,s,r);

  n = length(G);
  vdet = [0,0,0];
  for( i = 1, 3,
    G2 = vecextract(G,1<<i-1,1<<i-1)*Mod(1,p);
    if( !(vdet[i] = matdet(G2)),
      sol = kermodp(lift(G2),p)[2][,1];
      sol = vectorv(n, j, if( j <= i, sol[j]));
      return(sol)
    )
  );
\\ maintenant, on resout en dimension 3...
\\ d'abord, on se ramene au cas diagonal:
  x1 = [1,0,0]~;
  x2 = [-G2[1,2],G2[1,1],0]~;
  x3 = [G2[2,2]*G2[1,3]-G2[2,3]*G2[1,2],G2[1,1]*G2[2,3]-G2[1,3]*G2[1,2],G2[1,2]^2-G2[1,1]*G2[2,2]]~;
  while(1,
    if( issquare( N1 = -vdet[2]),                 s = sqrt(N1); sol = s*x1+x2; break);
    if( issquare( N2 = -vdet[3]/vdet[1]),         s = sqrt(N2); sol = s*x2+x3; break);
    if( issquare( N3 = -vdet[2]*vdet[3]/vdet[1]), s = sqrt(N3); sol = s*x1+x3; break);
    r = 1;
    while( !issquare( s = (1-N1*r^2)/N3), r = random(p));
    s = sqrt(s); sol = x1+r*x2+s*x3; break
  );
  sol = vectorv(n, j, if( j <= 3, sol[j]));
return(sol);
}
{Qfminim(G,factdetG=0) =
\\ Minimisation de la forme quadratique entiere non degeneree G,
\\   de dimension n >=2. On suppose que G est symetrique a coefficients entiers.
\\ Renvoit [G',U,factd] avec U \in GLn(Q) telle que G'=U~*G*U*constante est entiere
\\   de determinant minimal.
\\ Renvoit p si la reduction est impossible a cause de la non-solubilite locale
\\   en p (seulement en dimension 3-4).
\\ factd est la factorisation de abs(det(G'))
\\ Si on connait la factorisation de matdet(G), on peut la donner en parametre
\\   pour gagner du temps.
local(n,U,factd,detG,i,vp,Ker,dimKer,Ker2,dimKer2,sol,aux,p,di,m);

  n = length(G);
  U = matid(n);
  factd = matrix(0,2);
  detG = matdet(G);
  if( !factdetG, factdetG = factor(detG));
  i = 1;
  while(i <= length(factdetG[,1]),
    p = factdetG[i,1];
    if( p == -1, i++; next);
if( DEBUGLEVEL_qfsolve >= 4, print("p=",p,"^",factdetG[i,2]));
    vp = factdetG[i,2];
    if( vp == 0, i++; next);
\\ Le cas vp=1 n'est minimisable que si n est impair.
    if( vp == 1 && n%2 == 0,
      factd = concat(factd~, Mat([p,1])~)~;
      i++; next
    );
    Ker = kermodp(G,p); dimKer = Ker[1]; Ker = Ker[2];
\\ Rem: on a toujours dimKer <= vp
if( DEBUGLEVEL_qfsolve >= 4, print("dimKer = ",dimKer));
\\ cas trivial: G est divisible par p.
    if( dimKer == n,
      G /= p;
      factdetG[i,2] -= n;
      next
    );
    G = Ker~*G*Ker;
    U = U*Ker;
\\ 1er cas: la dimension du noyau est plus petite que la valuation
\\ alors le noyau mod p contient un noyau mod p^2
    if( dimKer < vp,
if( DEBUGLEVEL_qfsolve >= 4, print("cas 1"));
      Ker2 = kermodp(matrix(dimKer,dimKer,j,k,G[j,k]/p),p);
      dimKer2 = Ker2[1]; Ker2 = Ker2[2];
      for( j = 1, dimKer2, Ker2[,j] /= p);
      Ker2 = matdiagonalblock([Ker2,matid(n-dimKer)]);
      G = Ker2~*G*Ker2;
      U = U*Ker2;
      factdetG[i,2] -= 2*dimKer2;
if( DEBUGLEVEL_qfsolve >= 4, print("fin cas 1"));
      next
    );

\\ Maintenant, vp = dimKer
\\ 2eme cas: la dimension du noyau est >=2 et contient un element de norme 0 mod p^2
    if( dimKer > 2 ||
       (dimKer == 2 && issquare(di=Mod((G[1,2]^2-G[1,1]*G[2,2])/p^2,p))),
\\ on cherche dans le noyau un elt de norme p^2...
      if( dimKer > 2,
if( DEBUGLEVEL_qfsolve >= 4, print("cas 2.1"));
        dimKer = 3;
        sol = Qfsolvemodp(matrix(3,3,j,k,G[j,k]/p),p)
      ,
if( DEBUGLEVEL_qfsolve >= 4, print("cas 2.2"));
        if( G[1,1]%p^2 == 0,
          sol = [1,0]~
        , sol = [-G[1,2]/p+sqrt(di),Mod(G[1,1]/p,p)]~
        )
      );
      sol = centerlift(sol); sol /= content(sol);
if( DEBUGLEVEL_qfsolve >= 4, print("sol=",sol));
      Ker = vectorv(n, j, if( j<= dimKer, sol[j], 0)); \\ on complete avec des 0
      Ker = completebasis(Ker,1);
      Ker[,n] /= p;
      G = Ker~*G*Ker;
      U = U*Ker;
      factdetG[i,2] -= 2;
if( DEBUGLEVEL_qfsolve >= 4, print("fin cas 2"));
      next
    );

\\ Maintenant, vp = dimKer <= 2
\\   et le noyau ne contient pas de vecteur de norme p^2...

\\ Dans certains cas, on peut echanger le noyau et l'image
\\   pour continuer a reduire.

    m = (n-1)\2-1;
    if( ( vp == 1 && issquare(Mod(-(-1)^m*matdet(G)/G[1,1],p)))
     || ( vp == 2 && n%2 == 1 && n >= 5)
     || ( vp == 2 && n%2 == 0 && !issquare(Mod((-1)^m*matdet(G)/p^2,p)))
    ,
if( DEBUGLEVEL_qfsolve >= 4, print("cas 3"));
      Ker = matid(n);
      for( j = dimKer+1, n, Ker[j,j] = p);
      G = Ker~*G*Ker/p;
      U = U*Ker;
      factdetG[i,2] -= 2*dimKer-n;
if( DEBUGLEVEL_qfsolve >= 4, print("fin cas 3"));
      next
    );

\\ On n'a pas pu minimiser.
\\ Si n == 3 ou n == 4 cela demontre la non-solubilite locale en p.
    if( n == 3 || n == 4,
if( DEBUGLEVEL_qfsolve >= 1, print("pas de solution locale en ",p));
      return(p));

if( DEBUGLEVEL_qfsolve >= 4, print("plus de minimisation possible"));
    factd = concat(factd~,Mat([p,vp])~)~;
    i++
  );
\\print("Un=",log(vecmax(abs(U))));
  aux = qflll(U);
\\print("Ur=",log(vecmax(abs(U*aux))));
return([aux~*G*aux,U*aux,factd]);
}
{mymat(qfb) = qfb = Vec(qfb);[qfb[1],qfb[2]/2;qfb[2]/2,qfb[3]];
}
{Qfbsqrtgauss(G,factdetG) =
\\ pour l'instant, ca ne marche que pour detG sans facteur carre
\\ sauf en 2, ou la valuation est 2 ou 3.
\\ factdetG contient la factorisation de 2*abs(disc G)
local(a,b,c,d,m,n,p,aux,Q1,M);
if( DEBUGLEVEL_qfsolve >=3, print("entree dans Qfbsqrtgauss avec",G,factdetG));
  G = Vec(G);
  a = G[1]; b = G[2]/2; c = G[3]; d = a*c-b^2;

\\ 1ere etape: on resout m^2 = a (d), m*n = -b (d), n^2 = c (d)
  m = n = Mod(1,1);
  factdetG[1,2] -= 3;
  for( i = 1, length(factdetG[,1]),
    if( !factdetG[i,2], next);
    p = factdetG[i,1];
    if( gcd(a,p) == 1,
      aux = sqrt(Mod(a,p));
      m = chinese(m,aux);
      n = chinese(n,-b/aux)
    ,
      aux = sqrt(Mod(c,p));
      n = chinese(n,aux);
      m = chinese(m,-b/aux)
    )
  );
  m = centerlift(m);  n = centerlift(n);
if( DEBUGLEVEL_qfsolve >=4, print("m=",m); print("n=",n));

\\ 2eme etape: on construit Q1 de det -1 tq Q1(x,y,0) = G(x,y)
  Q1 = [(n^2-c)/d, (m*n+b)/d, n ;
        (m*n+b)/d, (m^2-a)/d, m ;
        n,         m,         d ];
  Q1 = -matadjoint(Q1);

\\ 3eme etape: on reduit Q1 jusqu'a [0,0,-1;0,1,0;-1,0,0]
  M = LLLgoon3(Q1)[2][3,];
  if( M[1] < 0, M = -M);
if( DEBUGLEVEL_qfsolve >=3, print("fin de Qfbsqrtgauss "));
  if( M[1]%2,
    return(Qfb(M[1],2*M[2],2*M[3]))
  , return(Qfb(M[3],-2*M[2],2*M[1])));
}
{class2(D,factdetG,Winvariants,U2) =
\\ On suit l'algorithme de Shanks/Bosma-Stevenhagen
\\ pour construire tout le 2-groupe de classes
\\ seulement pour D discriminant fondamental.
\\ lorsque D = 1(4), on travaille avec 4D.
\\ Si on connait la factorisation de abs(2*D),
\\ on peut la donner dans factdetG, et dans ce cas le reste
\\ de l'algorithme est polynomial.
\\ Si Winvariants est donne, l'algorithme s'arrete des qu'un element d'invariants=Winvariants
\\ est trouve.

local(factD,n,rang,m,listgen,vD,p,vp,aux,invgen,im,Ker,Kerim,listgen2,G2,struct,E,compteur,red);

if( DEBUGLEVEL_qfsolve >= 1, print("Construction du 2-groupe de classe de discriminant ",D));
  if( D%4 == 2 || D%4 == 3, print("Discriminant not congruent to 0,1 mod 4");return(0));

  if( D==-4, return([[1],[Qfb(1,0,1)]]));

  if( !factdetG, factdetG = factor(2*abs(D)));
  factD = concat([-1],factdetG[,1]~);
  if( D%4 == 1, D *= 4; factdetG[1,2] += 2);

  n = length(factD); rang = n-3;
  if(D>0, m = rang+1, m = rang);
if( DEBUGLEVEL_qfsolve >= 3, print("factD = ",factD));
  listgen = vector(max(0,m));

  if( vD = valuation(D,2),
    E = Qfb(1,0,-D/4)
  , E = Qfb(1,1,(1-D)/4)
  );
if( DEBUGLEVEL_qfsolve >= 3, print("E = ",E));

  if( type(Winvariants) == "t_COL" && (Winvariants == 0 || length(matinverseimage(U2*Mod(1,2),Winvariants))>0), return([[1],[E]]));

  for( i = 1, m, \\ on  ne regarde pas factD[1]=-1, ni factD[2]=2
    p = factD[i+2];
    vp = valuation(D,p);
    aux = p^vp;
    if( vD,
      listgen[i] = Qfb(aux,0,-D/4/aux)
    , listgen[i] = Qfb(aux,aux,(aux-D/aux)/4))
  );
  if( vD == 2 && D%16 != 4,
    m++; rang++; listgen = concat(listgen,[Qfb(2,2,(4-D)/8)]));
  if( vD == 3,
    m++; rang++; listgen = concat(listgen,[Qfb(2^(vD-2),0,-D/2^vD)]));

if( DEBUGLEVEL_qfsolve >= 3, print("listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 2, print("rang = ",rang));

  if( !rang, return([[1],[E]]));

 invgen = Qflisteinvariants(listgen,factD)[2]*Mod(1,2);
if( DEBUGLEVEL_qfsolve >= 3, printp("invgen = ",lift(invgen)));

  struct = vector(m,i,2);
  im = lift(matinverseimage(invgen,matimage(invgen)));
  while( (length(im) < rang)
  || (type(Winvariants) == "t_COL" && length(matinverseimage(concat(invgen,U2),Winvariants) == 0)),
    Ker = lift(matker(invgen));
    Kerim = concat(Ker,im);
    listgen2 = vector(m);
    for( i = 1, m,
      listgen2[i] = E;
      for( j = 1, m,
        if( Kerim[j,i],
          listgen2[i] = qfbcompraw(listgen2[i],listgen[j])));
      if( norml2(Kerim[,i]) > 1,
        red = QfbReduce(aux=mymat(listgen2[i]));
        aux = red~*aux*red;
        listgen2[i] = Qfb(aux[1,1],2*aux[1,2],aux[2,2]))
    );
    listgen = listgen2;
    invgen = invgen*Kerim;

if( DEBUGLEVEL_qfsolve >= 4, print("listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 4, printp("invgen = ",lift(invgen)));

    for( i = 1, length(Ker),
      G2 = Qfbsqrtgauss(listgen[i],factdetG);
      struct[i] <<= 1;
      listgen[i] = G2;
      invgen[,i] = Qflisteinvariants(G2,factD)[2][,1]*Mod(1,2)
    );

if( DEBUGLEVEL_qfsolve >= 3, print("listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 3, printp("invgen = ",lift(invgen)));
if( DEBUGLEVEL_qfsolve >= 3, print("struct = ",struct));

    im = lift(matinverseimage(invgen,matimage(invgen)))
  );

  listgen2 = vector(rang);
  for( i = 1, rang,
    listgen2[i] = E;
    for( j = 1, m,
      if( im[j,i],
        listgen2[i] = qfbcompraw(listgen2[i],listgen[j])));
    if( norml2(im[,i]) > 1,
      red = QfbReduce(aux=mymat(listgen2[i]));
      aux = red~*aux*red;
      listgen2[i] = Qfb(aux[1,1],2*aux[1,2],aux[2,2]))
  );
  listgen = listgen2;
\\ listgen = vector(rang,i,listgen[m-rang+i]);
  struct = vector(rang,i,struct[m-rang+i]);

if( DEBUGLEVEL_qfsolve >= 2, print("listgen = ",listgen));
if( DEBUGLEVEL_qfsolve >= 2, print("struct = ",struct));

return([struct,listgen]);
}
{Qfsolve(G,factD) =
\\ Resolution de la forme quadratique X^tGX=0 de dimension n >= 3.
\\ On suppose que G est entiere et primitive.
\\ La solution peut etre un vectorv ou une matrice.
\\ S'il n'y a pas de solution, alors renvoit un p tel
\\ qu'il n'existe pas de solution locale en p.
\\
\\ Si on connait la factorisation de -abs(2*matdet(G)),
\\ on peut la passer par le parametre factD pour gagner du temps.

local(n,M,signG,d,Min,U,codim,aux,G1,detG1,M1,subspace1,G2,subspace2,M2,solG2,Winvariants,dQ,factd,U2,clgp2,V,detG2,dimseti,solG1,sol,Q);

if( DEBUGLEVEL_qfsolve >= 1, print("entree dans Qfsolve"));
\\
\\ 1ere reduction des coefficients de G
\\

  n = length(G);
  M = IndefiniteLLL(G);
  if( type(M) == "t_COL",
if( DEBUGLEVEL_qfsolve >= 1, print("solution ",M));
    return(M));
  G = M[1]; M = M[2];

\\ Solubilite reelle
  signG = qfsign(G);
  if( signG[1] == 0 || signG[2] == 0,
if( DEBUGLEVEL_qfsolve >= 1, print("pas de solution reelle"));
    return(-1));
  if( signG[1] < signG[2], G = -G; signG = signG*[0,1;1,0]);

\\ Factorisation du determinant
  d = matdet(G);
  if( !factD,
if( DEBUGLEVEL_qfsolve >= 1, print("factorisation du determinant"));
    factD = factor(-abs(2*d));
if( DEBUGLEVEL_qfsolve >= 1, print(factD))
  );
  factD[1,2] = 0;
  factD[2,2] --;

\\
\\ Minimisation et solubilite locale.
\\

if( DEBUGLEVEL_qfsolve >= 1, print("minimisation du determinant"));
  Min = Qfminim(G,factD);
  if( type(Min) == "t_INT",
if( DEBUGLEVEL_qfsolve >= 1, print("pas de solution locale en ",Min));
    return(Min));

  M = M*Min[2];
  G = Min[1];
\\  Min[3] contient la factorisation de abs(matdet(G));

if( DEBUGLEVEL_qfsolve >= 4, print("G minime = ",G));
if( DEBUGLEVEL_qfsolve >= 4, print("d=",d));

\\ Maintenant, on sait qu'il y a des solutions locales
\\ (sauf peut-etre en 2 si n==4),
\\ si n==3, det(G) = +-1
\\ si n==4, ou si n est impair, det(G) est sans facteur carre.
\\ si n>=6, det(G) toutes les valuations sont <=2.

\\ Reduction de G et recherche de solution triviales
\\ par exemple quand det G=+-1, il y en a toujours.

if( DEBUGLEVEL_qfsolve >= 1, print("reduction"));
  U = IndefiniteLLL(G);
  if(type(U) == "t_COL",
if( DEBUGLEVEL_qfsolve >= 1, print("solution ",M*U));
    return(M*U));
  G = U[1]; M = M*U[2];

\\
\\ Quand n >= 6 est pair, il faut passer en dimension n+1
\\ pour supprimer tous les carres de det(G).
\\

  if( n >= 6 && n%2 == 0 && matsize(Min[3])[1] != 0,
if( DEBUGLEVEL_qfsolve >= 1, print("On passe en dimension ",n+1));
    codim = 1; n++;
\\ On calcule le plus grand diviseur carre de d.
    aux = prod( i = 1, matsize(Min[3])[1], if( Min[3][i,1] == 2, Min[3][i,1], 1));
\\ On choisit le signe de aux pour que la signature de G1
\\ soit la plus equilibree possible.
    if( signG[1] > signG[2],
      signG[2] ++; aux = -aux
    , signG[1] ++
    );
    G1 = matdiagonalblock([G,Mat(aux)]);
    detG1 = 2*matdet(G1);
    for( i = 2, length(factD[,1]),
      factD[i,2] = valuation(detG1,factD[i,1]));
    factD[2,2]--;
    Min = Qfminim(G1,factD);
    G1 = Min[1];
    M1 = Min[2];
    subspace1 = matrix(n,n-1,i,j, i == j)
  , codim = 0;
    G1 = G;
    subspace1 = M1 = matid(n)
  );
\\ Maintenant d est sans facteur carre.

\\
\\ Si d n'est pas +-1, il faut passer en dimension n+2
\\

  if( matsize(Min[3])[1] == 0, \\ if( abs(d) == 1,
if( DEBUGLEVEL_qfsolve >= 2, print(" detG2 = 1"));
     G2 = G1;
     subspace2 = M2 = matid(n);
     solG2 = LLLgoon(G2,1)
  ,
if( DEBUGLEVEL_qfsolve >= 1, print("On passe en dimension ",n+2));
    codim += 2;
    subspace2 = matrix( n+2, n, i, j, i == j);
    d = prod( i = 1, matsize(Min[3])[1],Min[3][i,1]);    \\ d = abs(matdet(G1));
    if( signG[2]%2 == 1, d = -d);                        \\ d = matdet(G1);
    if( Min[3][1,1] == 2, factD = [-1], factD = [-1,2]); \\ si d est pair ...
    factD = concat(factD, Min[3][,1]~);
if( DEBUGLEVEL_qfsolve >= 4, print("factD=",factD));

\\ Solubilite locale en 2 (c'est le seul cas qui restait a traiter !!)
    if( n == 4 && d%8 == 1,
      if( QfWittinvariant(G,2) == 1,
if( DEBUGLEVEL_qfsolve >= 1, print("pas de solution locale en 2"));
        return(2)));

\\
\\ Construction d'une forme Q de dim 2, a partir de ces invariants.
\\
    Winvariants = vectorv(length(factD));

\\ choix de la signature de Q.
\\ invariant reel et signe du discriminant.
    dQ = abs(d);
    if( signG[1] ==signG[2], dQ = dQ; Winvariants[1] = 0); \\ signQ = [1,1];
    if( signG[1] > signG[2], dQ =-dQ; Winvariants[1] = 0); \\ signQ = [2,0];
    if( n == 4 && dQ%4 != 1, dQ *= 4);
    if( n >= 5, dQ *= 8);

\\ invariants p-adiques
\\ pour p = 2, on ne choisit pas.
    if( n == 4,
if( DEBUGLEVEL_qfsolve >= 1, print("calcul des invariants locaux de G1"));
      aux = Qflisteinvariants(-G1,factD)[2][,1];
      for( i = 3, length(factD), Winvariants[i] = aux[i])
    ,
      aux = (-1)^((n-3)/2)*dQ/d; \\ ici aux = 8 ou -8
      for( i = 3, length(factD), Winvariants[i] = hilbert(aux,factD[i],factD[i]) > 0)
    );
    Winvariants[2] = sum( i = 1, length(factD), Winvariants[i])%2;

if( DEBUGLEVEL_qfsolve >= 1,
  print("Recherche d'un forme binaire de discriminant = ",dQ);
  print("et d'invariants de Witt = ",Winvariants));

\\ On construit le 2-groupe de classes de disc dQ,
\\ jusqu'a obtenir une combinaison des generateurs qui donne les bons invariants.
\\ En dim 4, il faut chercher parmi les formes du type q, 2*q
\\ car Q peut etre imprimitive.

    factd = matrix(length(factD)-1,2);
    for( i = 1, length(factD)-1,
      factd[i,1] = factD[i+1];
      factd[i,2] = valuation(dQ,factd[i,1]));
    factd[1,2]++;
    U2 = matrix(length(factD), n == 4, i,j, hilbert(2,dQ,factD[i])<0);
    clgp2 = class2(dQ,factd,Winvariants,U2);
if( DEBUGLEVEL_qfsolve >= 4, print("clgp2=",clgp2));

    clgp2 = clgp2[2];
    U = Qflisteinvariants(clgp2,factD)[2];
    if( n == 4, U = concat(U,U2));
if( DEBUGLEVEL_qfsolve >= 4, printp("U=",U));

    V = lift(matinverseimage(U*Mod(1,2),Winvariants*Mod(1,2)));
    if( !length(V), next);
if( DEBUGLEVEL_qfsolve >= 4, print("V=",V));

    if( dQ%2 == 1, Q = qfbprimeform(4*dQ,1), Q = qfbprimeform(dQ,1));
    for( i = 1, length(clgp2),
      if( V[i], Q = qfbcompraw(Q,clgp2[i])));
    Q = mymat(Q);
    if( norml2(V) > 1, aux = QfbReduce(Q); Q = aux~*Q*aux);
    if( n == 4 && V[length(V)], Q*=  2);
if( DEBUGLEVEL_qfsolve >= 2, print("Q=",Q));
if( DEBUGLEVEL_qfsolve >= 3, print("invariants de Witt de Q=",Qflisteinvariants(Q,factD)));

\\
\\ Construction d'une forme de dim=n+2 potentiellement unimodulaire
\\

    G2 = matdiagonalblock([G1,-Q]);
if( DEBUGLEVEL_qfsolve >= 4, print("G2=",G2));

if( DEBUGLEVEL_qfsolve >= 2, print("minimisation de la forme quadratique de dimension ",length(G2)));
\\ Minimisation de G2
    detG2 = matdet(G2);
    factd = matrix(length(factD)-1,2);
    for( i = 1, length(factD)-1,
      factd[i,2] = valuation(detG2, factd[i,1] = factD[i+1]));
if( DEBUGLEVEL_qfsolve >= 3, print("det(G2) = ",factd));
    Min = Qfminim(G2,factd);
    M2 = Min[2]; G2 = Min[1];
if( abs(matdet(G2)) > 2, print("******* ERREUR dans Qfsolve: det(G2) <> +-1 *******",matdet(G2));return(0));
if( DEBUGLEVEL_qfsolve >= 4, print("G2=",G2));

\\ Maintenant det(G2) = +-1

\\ On construit un seti pour G2 (Sous-Espace Totalement Isotrope)
if( DEBUGLEVEL_qfsolve >= 2, print("recherche d'un espace de solutions pour G2"));
    solG2 = LLLgoon(G2,1);
    if( matrix(codim+1,codim+1,i,j,solG2[1][i,j]) != 0,
if( DEBUGLEVEL_qfsolve >= 2, print(" pas assez de solutions dans G2"));return(0))
  );

\\ G2 doit avoir un espace de solutions de dimension > codim
  dimseti = 0;
  while( matrix(dimseti+1,dimseti+1,i,j,solG2[1][i,j]) == 0, dimseti ++);
  if( dimseti <= codim,
print("dimseti = ",dimseti," <= codim = ",codim);
print("************* ERREUR : pas assez de solutions pour G2"); return(0));
  solG2 = matrix(length(G2),dimseti,i,j,solG2[2][i,j]);
if( DEBUGLEVEL_qfsolve >= 3, print("solG2=",solG2));

\\ La solution de G1 se trouve a la fois dans solG2 et dans subspace2
if( DEBUGLEVEL_qfsolve >= 1, print("Reconstruction de la solution de G1"));
  solG1 = matintersect(subspace2,M2*solG2);
  solG1 = subspace2~*solG1;
if( DEBUGLEVEL_qfsolve >= 3, print("solG1 = ",solG1));

\\ La solution de G se trouve a la fois dans solG et dans subspace1
if( DEBUGLEVEL_qfsolve >= 1, print("Reconstruction de la solution de G"));
  sol = matintersect(subspace1,M1*solG1);
  sol = subspace1~*sol;
  sol = M*sol;
  sol /= content(sol);
  if( length(sol) == 1, sol = sol[,1]);
if( DEBUGLEVEL_qfsolve >= 3, print("sol = ",sol));
if( DEBUGLEVEL_qfsolve >= 1, print("fin de Qfsolve"));
  return(sol);
}
{matdiagonalblock(v) =
local(lv,lt,M);
  lv = length(v);
  lt = sum( i = 1, lv, length(v[i]));
  M = matrix(lt,lt);
  lt = 0;
  for( i = 1, lv,
    for( j = 1, length(v[i]),
      for( k = 1, length(v[i]),
        M[lt+j,lt+k] = v[i][j,k]));
    lt += length(v[i])
  );
return(M);
}






