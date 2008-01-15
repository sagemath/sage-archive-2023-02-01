/****************************************************************************/
/*                                                                          */
/*   Description: computes the cardinality of an elliptic curve / Fp        */
/*   Time-stamp: <2004-02-26 11:10:56 kb>, version 1.3                      */
/*   Original authors: Christophe Doche   <cdoche@math.u-bordeaux.fr>       */
/*                     Sylvain Duquesne <duquesne@math.u-bordeaux.fr>       */
/*   Universite Bordeaux I, Laboratoire A2X                                 */
/*   For the AREHCC project, see http://www.arehcc.com                      */
/*   Created: Apr 30 2003                                                   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                  TOOLS                                   */
/****************************************************************************/

\r sea_init.gp
global(VERBOSE, EARLY_ABORT, BOUND_ONE_ROOT, MAXL = 199);

\\Implement Vecrev for GP 2.3
vecrev(w) = w=Vec(w);vecextract(w, Str(#w,"..",1));

\\Builds the modular equation corresponding to the vector list
list_to_pol(list) = Pol( vector(#list, i, Pol(list[i], 'J)), 'x );

\\Performs the inverse operation
pol_to_list(T) =
{local(list, tmp);
  list = [];
  forstep (i = poldegree(T), 0, -1,
    tmp = Vec(polcoeff(T, i, 'x));
    list = concat(list, if (#tmp == 1, tmp[1], [tmp]))
  ); list;}

\\Given power series s1 and s2, finds a polynomial P
\\such that s2 = P(s1)
find_transformation(s2, s1) =
{ local(vs1, vs2, degP, s1pl, s1i, P, invs1coeff, Pcoeff, d);
  vs1 = valuation(s1, 'z);
  vs2 = valuation(s2, 'z); degP = vs2 \ vs1;
  if (vs2%vs1 || degP <= 0, error("transformation cannot exist"));
  s1i = 1;
  s1pl = vector(degP, i, s1i *= s1);
  P = 0;
  invs1coeff = 1/polcoeff(s1, vs1);
  for (i = 0, degP - 1,

    Pcoeff = polcoeff(s2, vs2 - i*vs1) * invs1coeff;
    d = degP - i; /* > 0 */
    P += Pcoeff*'x^d;
    s2 -= Pcoeff*s1pl[d]; /* * s1^d */
    if (!truncate(s2), error("not enough terms to determine the polynomial"))
  );
  Pcoeff = polcoeff(s2, vs2 - degP*vs1) * invs1coeff;
  s2 -= Pcoeff;
  if (truncate(s2), error("polynomial expression does not match"));
  P + Pcoeff;}

/****************************************************************************/
/*                             ELLIPTIC CURVES                              */
/****************************************************************************/

\\ Finds a random point on E
find_pt_aff(E, p) =
{ local(x, y);
  until (#y, x = random(p); y = ellordinate(E, x));
  [Mod(x,p), y[1]];}

\\ E has order o[1], ..., or o[#o], draw random points until all solutions
\\ but one are eliminated
choose_card(o, E, p)=
{ local(P, t, lasto, lastgood, nbo = #o);
  if (nbo == 1, return (o[1]));
  o = vecsort(o); \\ minimize max( o[i+1] - o[i] )
  lastgood = o[#o];
  while (1,
    lasto = 0; t = [0];
    P = find_pt_aff(E, p);
    for (i = 1, #o, \\ t = ellpow(E, P, lasto)
      if (!o[i], next);
      t = elladd(E, t, ellpow(E, P, o[i] - lasto)); /* P^o[i] */
      lasto = o[i];
      if (t != [0],
        nbo--; if (nbo == 1, return (lastgood));
        o[i] = 0;
      ,
        lastgood = o[i];
      );
    );
  );
}

\\ Returns the double of the point [P[1], P[2]*Y]
\\ Compute modulo Y^2 - (X^3 + a4X + a6)
my_elldbl(E, P) =
{ local(lambda, C);

  if (P == [0], return([0]));
  if (P == [X,1], \\ save 1 square (lambda*DRHS is free)
    lambda = DRHS/(4*RHS); \\ 1/2 x the "standard" value
    C = lambda*DRHS - 2 * X;
    [C, 2*lambda*(X - C) - 1]
  ,
    lambda = (3*P[1]^2 + E.a4)/(2*P[2]*RHS);
    C = sqr(lambda)*RHS - 2*P[1];
    [C, lambda*(P[1] - C) - P[2]]
  );
}

\\Returns the addition of [P[1], P[2]*Y] and of [Q[1], Q[2]*Y]
\\Computations are done modulo Y^2 - (X^3 + a4X + a6)
\\An inversion is equivalent to 4M, so that this function requires about 7M
\\which is the same as with the method using ell-division polynomials
\\Working in mixed projective coordinates would require 11M
my_elladd(E, P, Q) =
{ local(lambda, C);
  if (P == [0], return(Q));
  if (Q == [0], return(P));
  if (P[1] == Q[1],
    return ( if (P[2] == Q[2], my_elldbl(E, P), [0]) )
  );
  lambda = (P[2] - Q[2])/(P[1] - Q[1]);
  C = sqr(lambda) * RHS - P[1] - Q[1];
  [C, lambda*(P[1] - C) - P[2]];}

my_ellpow(E, P, m) =
{ local(v, Ptmp);

  v = binary(m); Ptmp = [0];
  for (nb = 1, #v,
    Ptmp = my_elldbl(E, Ptmp);
    if (v[nb], Ptmp = my_elladd(E, Ptmp, P))
  ); Ptmp; }

\\Gives the first precS terms of the Weierstrass series related to
\\E: y^2 = x^3 + a4x + a6
find_coeff(a4, a6, precS) =
{ local(res);

  res = vector(precS); if (!precS, return(res));
  res[1] = -a4/5;     if (precS == 1, return(res));
  res[2] = -a6/7;
  for (k = 3, precS,
    res[k] = 3*sum(h=1, k-2, res[h]*res[k-1-h]) / ((k-2)*(2*k+3))
  ); res;}

\\Writes P(x, y) as A(x)y + B(x) replacing y^2 by x^3 + a4x + a6 in P(x, y)
reduce(P, E, h, p) = substpol(lift(P), 'y^2, 'x^3 + E.a4*'x + E.a6) * Mod(1,p)*Mod(1,h);

\\Computes the n-division polynomial modulo the polynomial h \in Fp[x]
elldivpol(E, n, h, p) =
{local(f, g, a4, a6, N, inv2y2);
  a4 = E.a4;
  a6 = E.a6; N = max(5, n+1);
  f = vector(N);
  g = vector(N);
  inv2y2 = 1/Mod('x^3+a4*'x+a6,h) * Mod(2,p)^-1;
  f[1] = Mod(1,h);
  g[1] = Mod(1,h);
  f[2] = 2*'y*Mod(1,h);
  g[2] = reduce(f[2]^2, E, h, p);
  f[3] = (3*'x^4 + 6*a4*'x^2 + 12*a6*'x - a4^2)*Mod(1,h);
  g[3] = reduce(f[3]^2, E, h, p);
  f[4] = 4*'y*('x^6 + 5*a4*'x^4 + 20*a6*'x^3 - 5*a4^2*'x^2 - 4*a4*a6*'x - 8*a6^2 - a4^3)*Mod(1,h);
  g[4] = reduce(f[4]^2, E, h, p);
  f[5] = reduce(f[4]*f[2]*g[2]-f[3]*g[3]*f[1], E, h, p);
  g[5] = reduce(f[5]^2, E, h, p);
  for (m = 3, n\2,
    f[2*m] = reduce(f[m]*('y*inv2y2)*(f[m+2]*g[m-1] - f[m-2]*g[m+1]), E, h, p);
    g[2*m] = reduce(f[2*m]^2, E, h, p);
    f[2*m+1] = reduce(f[m+2]*f[m]*g[m]-f[m+1]*g[m+1]*f[m-1], E, h, p);
    g[2*m+1] = reduce(f[2*m+1]^2, E, h, p);
  );
  lift(f[n]);}

\\ Finds E' q-isogenous to E and the trace term p1 from canonical modular
\\   equation meqn
\\ E: elliptic curve, q: a prime, meqn: canonical modular equation
\\ g: root of meqn defining isogenous curve Eb.
find_isogenous_from_canonical(E, q, meqn, g, p) =
{local(j, E4, E6, delta, is, Dx, DJ, Dxg, DJg, Dxx, px, pJ, ExJ, dx, dJ, deltal, E4l,
  a4tilde, jl, a6tilde, p1, E2s, gd, jd, E0b, Dgd, Djd, E0bd, f, fd,
  Dgs, Djs, jld, E6l);

  E4 = -E.a4/3;
  E6 = -E.a6/2; delta = (E4^3 - E6^2)/1728; j = E.j;
  Dx = deriv(meqn, 'x);
  DJ = deriv(meqn, 'J);
  Dxg = subst(Dx, 'x, g); px = subst(Dxg, 'J, j); dx = g * px;
  DJg = subst(DJ, 'x, g); pJ = subst(DJg, 'J, j); dJ = j * pJ;
  Dxx = deriv(Dx, 'x);
  ExJ = subst(deriv(Dxg,'J), 'J, j);
  is = 12 / gcd(12,q-1);
  deltal = delta * g^(12/is)/q^12;
  if (!dJ,
    if (VERBOSE, print("Division by zero ", lift(E)," "p));
    E4l = E4/q^2;
    jl = E4l^3/deltal;
    E6l = -sqrt((jl - 1728)*deltal);
    p1 = 0
  ,
    E2s = (-12*E6*dJ)/(is*E4*dx);
    gd = -(is/12)*E2s*g;
    jd = -E4^2*E6/delta;
    E0b = E6/(E4*E2s);
    Dgd = gd*px + g * (gd*subst(subst(Dxx, 'x, g), 'J, j) + jd*ExJ);
    Djd = jd*pJ + j * (jd*subst(deriv(DJg,'J), 'J, j) + gd*ExJ);
    E0bd = ((-is*Dgd)/12 - E0b*Djd)/dJ;

    E4l = (E4 - E2s*(12*E0bd/E0b + 6*E4^2/E6 - 4*E6/E4 - E2s))/q^2;
    jl = E4l^3/deltal;
    f = q^is/g;
    fd = is*E2s*f/12;
    Dgs = subst(subst(Dx, 'x, f), 'J, jl);
    Djs = subst(subst(DJ, 'x, f), 'J, jl);
    jld = -fd*Dgs/(q*Djs);
    E6l = -E4l*jld/jl;
    p1 = -q*E2s/2;
  );
  a4tilde = -3*q^4*E4l;
  a6tilde = -2*q^6*E6l;
  [ellinit([0, 0, 0, a4tilde, a6tilde]*Mod(1,p)), p1];}

\\Finds the isogenous EC, and the sum of the x-coordinates of the points in
\\ the kernel of the isogeny E -> Eb
\\E: elliptic curve, q: a prime, meqn: Atkin modular equation
\\g: root of meqn defining isogenous curve Eb.
find_isogenous_from_Atkin(E, q, meqn, g, p) =
{ local(Roots, jt, j, dx, dJ, E4, E6, dxstar,
    Dx, DJ, Dxg, DJg, Dxxg, DJJg, DxJg, gprime, E4t, E6t, a4t, a6t, px, pJ,
    pxstar, pJstar, dJstar, u1, u2, p1, Eb, check, u,v);

  E4 = -E.a4/3;
  E6 = -E.a6/2; j = E.j;
  Dx = deriv(meqn, 'x);
  DJ = deriv(meqn, 'J);
  Dxg = subst(Dx, 'x, g); px = subst(Dxg, 'J, j); dx = g * px;
  DJg = subst(DJ, 'x, g); pJ = subst(DJg, 'J, j); dJ = j * pJ;
  Dxxg = subst(deriv(Dx, 'x), 'x, g);
  DJJg = deriv(DJg, 'J);
  DxJg = deriv(Dxg, 'J);
  if (!dx || !E4,
    if (VERBOSE,
      print("find_isogenous_from_Atkin: division by zero at prime ", q));
    return([]);
  );
  gprime = (g*E6*dJ)/(E4 * dx);
  u1 = (-gprime*subst(Dxxg, 'J,j)
        + 2*j*subst(DxJg, 'J, j)*E6/E4
        - (E6^2/(gprime*E4^2))*j*(pJ + j*subst(DJJg, 'J, j))) / px
      + E6/(3*E4) - E4^2/(2 * E6);

  Roots = polrootsmod(subst(meqn, 'x, g), p);
  forstep (k = #Roots, 1, -1,
    jt = Roots[k];
    pxstar = subst(Dxg, 'J, jt); dxstar = g * pxstar;
    pJstar = subst(DJg, 'J, jt); dJstar = q * jt * pJstar;
    u = dxstar*dJ*E6;
    v = dJstar*dx*E4;
    E4t = (sqr(u)*jt)/(sqr(v)*(jt - 1728));
    E6t = (u * E4t)/v;
    u2 = (-gprime*subst(Dxxg, 'J, jt)
          + 2*q*jt*subst(DxJg, 'J, jt)*E6t/E4t
          - q^2*(E6t^2/(gprime*E4t^2))*jt*(pJstar + jt*subst(DJJg, 'J, jt))
         ) / pxstar
       + q*(E6t/(3*E4t) - E4t^2/(2*E6t));
    p1 = (u1 - u2)*6*q/2;
    a4t = -3*q^4*E4t;
    a6t = -2*q^6*E6t;
    Eb = ellinit([0, 0, 0, a4t, a6t]*Mod(1,p));
    check = find_kernel(E, q, Eb, p1, p);
    if (check[2], return([Eb, check[1]]))
  );
  error(" kernel not found");}

\\Finds the kernel polynomial h, dividing the ell-division polynomial from the
\\isogenous curve Eb and trace term pp1.
\\Uses CCR algorithm and returns [h, 1] if the result is correct or [h, 0]
\\otherwise. This last case meant that E and Eb were not isogenous
find_kernel(E, ell, Eb, pp1, p) =
{ local(ext, deg, a4, a6, Coeff, Coefftilde, psi2, Dpsi2, list, tsil, o, r,
        M, dim, N, V, K,K1,K2, v, tlist, tmp, h);
  ext = 2;
  deg = (ell - 1)\2;
  dim = deg+ext;
  if ((2*dim + 3) >= p, error("division by zero"));
  a4 = E.a4;
  a6 = E.a6;
  Coeff      = find_coeff(a4,       a6, dim + 1);
  Coefftilde = find_coeff(Eb.a4, Eb.a6, dim + 1);
  psi2 = 4*'x^3 + 4*a4*'x + 4*a6;
  Dpsi2 = 6*'x^2 + 2*a4;
  list = vector(dim); list[1] = r = Dpsi2;
  for (k = 2, dim,
    tsil = vecrev(r);
    r = tsil[2]*Dpsi2;
    for (j = 3, #tsil,
      o = j - 1;
      r += o*tsil[o+1]*(Dpsi2*'x + (o-1)*psi2)*'x^(o-2)
    );
    list[k] = r;
  );
  M = matrix(dim, dim + 2);
  for (k = 1, dim,
    tsil = vecrev(list[k]) * 2/(2*k)!;
    for (i = 1, #tsil, M[k,i] = lift(tsil[i]));
  );
  N = vecextract(M, Str("..",dim));
  V = vectorv(dim, i, lift(Coefftilde[i] - Coeff[i]));
  v = concat(FpM_gauss(N, V, p), [0,0]~);
  K = FpM_ker(M, p) * Mod(1,p);
  if (#K != 2, error("trace not determined in a unique way"));
  K1 = K[,1] / K[1,1];
  K2 = K[, 2] - K[1, 2]*K1;
  K2 /= K2[2];
  K1 -= K1[2]*K2;
  v += K1*(deg - v[1]);
  v += K2*(pp1 - v[2]);
  tlist = vector(dim+1); tlist[dim+1] = 1;
  for (k = 1, dim,
    tmp = -v[k+1] - sum(i = 1, k - 1, tlist[dim-i+1]*v[k-i+1]);
    tlist[dim-k+1] = tmp/k
  );
  h = Pol(vector(deg+1 , i, tlist[dim+2-i]));
  [h, !vector(ext, i, tlist[i])];}

find_kernel_power(Eb, Ec, ell, meqn, meqntype, mpoly, kpoly, Ib, p) =
{ local(num_iso, mroots, tmp, Etmp, p1c, gtmp, check, Ic, kpoly_new);

  num_iso = find_numerator_isogeny(Eb, Ec, kpoly, ell+1);
  mroots = polrootsmod(mpoly, p);
  for (i = 1, #mroots,
    if (meqntype == "C",
      tmp = find_isogenous_from_canonical(Ec, ell, meqn, mroots[i], p);
      Etmp = tmp[1];
      p1c  = tmp[2];
      tmp = find_kernel(Ec, ell, Etmp, p1c, p);
      gtmp  = tmp[1];
      check = tmp[2] ,
    if (meqntype == "A",
      tmp = find_isogenous_from_Atkin(Ec, ell, meqn, mroots[i], p);
      if (!#tmp, return ([]));
      Etmp = tmp[1];
      gtmp = tmp[2];
      check = 1));
    if (check && \\check that the kernel kpoly is the good one
        elldivpol(Eb, ell, numerator(subst(gtmp, 'x, num_iso/kpoly^2)), p),
      Ic = subst(num_iso, 'x, Ib) / subst(kpoly, 'x, Ib)^2;
      kpoly_new = numerator(subst(gtmp, 'x, Ic));
      return([Etmp, kpoly_new, gtmp, Ic])
    )
  );
  error("failed to find kernel polynomial");}

compute_W(E, precS)=  1/'z^2 + subst(Polrev(find_coeff(E.a4, E.a6, precS),'z),'z,'z^2)*'z^2 + O('z^(2*precS));

\\Finds numerator phi of the isogeny between Eb and Ec whose denominator is h.
find_numerator_isogeny(Eb, Ec, h, precS) =
{ local(WEb, WEc, den);
  WEb = compute_W(Eb, precS);
  WEc = compute_W(Ec, precS);
  den = subst(h, 'x, WEb);
  find_transformation(den^2*WEc, WEb);}

\\Assume E defined over Fp
\\Checks that #E(F_p) = ord by computing ord*Pt for N random Pt on E.
\\Note that this is only probabilistic
check_order(E, ord, N, p) =
{
  for (i = 1, N,
    if (ellpow(E, find_pt_aff(E, p), ord) != [0], return (0));
  ); 1;}

/****************************************************************************/
/*                              EIGENVALUE                                  */
/****************************************************************************/

init_eigen(E,h,p)=
{
  X = Mod('x*Mod(1,p), h);
  RHS = Mod('x^3 + E.a4*'x + E.a6, h);
  DRHS= Mod(3*'x^2 + E.a4, h);
  Gr = FpXQ_pow('x^3 + lift(E.a4)*'x + lift(E.a6), (p-1)/2, lift(h), p)
       * Mod(Mod(1,p), h);
}

\\Finds the eigenvalue of the Frobenius given E, ell, h a factor of the
\\ell-division polynomial, p and tr the possible values for the trace
\\(useful for primes with one root less than BOUND_ONE_ROOT)
find_eigen_value(E, ell, h, p, tr = []) =
{ local(DRHS, RHS, X, Gr, BP, Dr, bont);
  init_eigen(E,h,p);
  Dr = BP = [X, 1];
  \\[0,Gr], BP, Dr are not points on the curve.
  \\To obtain the corresponding points, multiply the y-coordinates by Y
  if (!#tr,
    for (t = 1, (ell-1)/2,
      if (Dr[2] == Gr,  bont = t;     break);
      if (Dr[2] == -Gr, bont = ell-t; break);
      Dr = my_elladd(E, Dr, BP)
    )
  ,
    bont = centerlift(Mod(tr[1]/2, ell));
    Dr = my_ellpow(E, BP, bont);
    if (Dr[2] != Gr, bont = ell - bont);
  ); bont;}

\\Finds the eigenvalue of the Frobenius modulo ell^k given E, ell, k, h a factor
\\of the ell-division polynomial, lambda the previous eigen value and p
find_eigen_value_power(E, ell, k, h, lambda, p) =
{ local(DRHS, RHS, X, Gr, BP, Dr, bont, Ell);
  init_eigen(E,h,p);
  \\[0,Gr], BP, Dr are not points on the curve.
  \\To obtain the corresponding points, multiply the y-coordinates by Y
  Ell = ell^(k-1);
  BP = my_ellpow(E, [X, 1], Ell);
  Dr = my_ellpow(E, [X, 1], lambda);
  for (t = 0, ell - 1,
    if (Dr[2] ==  Gr, bont = t*Ell+ lambda; break);
    if (Dr[2] == -Gr, bont = Ell*ell - (t*Ell+ lambda); break);
    Dr = my_elladd(E, Dr, BP);
  ); bont;}

/****************************************************************************/
/*                                  TRACE                                   */
/****************************************************************************/

\\ this is the all-important routine: > 3/4 of total time is spent here
\\ Cantor-Zassenhaus variant, with blocks (size 2) to reduce the number of gcd
\\ computations
\\Determines the type of the prime q
/*
study_modular_eqn(q, mpoly, p) =
{local(r, R, g, XP, G, lg, L, L2, XP2, T);

  r = 0;
  if (poldegree(FpX_gcd(mpoly, deriv(mpoly), p)) > 0,
    T = "P"; g = [0];
  ,
    XP = FpXQ_pow(x, p, mpoly, p);
    G = FpX_gcd(XP - x, mpoly, p);
    lg = poldegree(G); \\ # of roots of G = #g below
    if (!lg, \\ compute r = degree of smallest divisor of mpoly
      T = "A"; g = [0];
      L = bkinit(XP, q, mpoly, p); \\ deg(mpoly) = q+1
      XP2 = brent_kung(XP, L, mpoly, p); \\ Frob^2 x
      if (poldegree(FpX_gcd(XP2 - x, mpoly, p)), return (["A",0,2]));
      r = 2; R = (q+1) >> 1;
      while(1,
        if (r > R, return (["A",0,q+1]));
        r += 2;
        XP = brent_kung(XP2, L, mpoly, p); \\ Frob^(r-1) x
        XP2 = brent_kung(XP, L, mpoly, p); \\ Frob^r x
        \\ 2 gcd at at time, to save half of them
        G = FpX_gcd((XP2 - x) * (XP - x), mpoly, p);
        lg = poldegree(G);
        if (lg,
\\ G has irreducible factors of the same degree, r-1 or r. Decide which is the
\\ case, without an extra gcd unless r(r-1) | lg.
          if (lg % r,
            r--  \\ it's r-1
          , \\ dubious case is rare
            if (lg % (r-1) == 0 && poldegree(FpX_gcd(XP - x, G, p)), r--)
          );
          break
        )
      );
    ,
      g = polrootsmod(G, p);
      if (lg == 1, T = "1",
      if (lg == 2, T = "E",
      if (lg == q+1, T = "2",
                     T = "P"))))
  ); [T, g[1], r];}
*/

\\ Berlekamp variant [override the above]. Faster, simpler, but uses more space.
study_modular_eqn(q, mpoly, p) =
{local(r, g, XP, G, lg, L, T, s);

  r = 0;
  if (poldegree(FpX_gcd(mpoly, deriv(mpoly), p)) > 0,
    T = "P"; g = 0;
  ,
    XP = FpXQ_pow('x, p, mpoly, p);
    G = FpX_gcd(XP - 'x, mpoly, p);
    lg = poldegree(G); \\ # of roots of G = #g below
    if (!lg, \\ compute r = degree of smallest divisor of mpoly ~ CanZass
      T = "A"; g = 0;
      L = FpXQ_matrix_pow(XP,q+1,q+1,mpoly,p); \\ deg(mpoly) = q+1
      s = q+1 - FpM_rank(L - 1, p); \\ # of irreducible factors
      r = (q+1) / s \\ equal degree
    ,
      g = polrootsmod(G, p)[1];
      if (lg == 1, T = "1",
      if (lg == 2, T = "E",
      if (lg == q+1, T = "2",
                     T = "P"))))
  ); [T, g, r];}

\\Returns the possible values of the trace when ell is an Atkin prime,
\\given r the splitting degree of the modular equation at J = E.j
find_trace_Atkin(ell, r, p) =
{local(val_pos, tmp, a, T, P, invp, pell);
  val_pos = []; pell = p % ell;
  P = factor(r)[,1]; invp = 1/pell % ell;
  for (teta = 0, ell - 1,
    if (kronecker(teta^2 - 4*pell, ell) < 0,
      T = 't^2 - teta*'t + pell;
      tmp = invp*teta*'t - 1;
      a = FpXQ_pow(tmp, r/P[1], T, ell);
      if (a != 1 && FpXQ_pow(a,P[1], T,ell) == 1,
        for (i = 2, #P,
          if (FpXQ_pow(tmp, r/P[i], T,ell) == 1, next(2))
        );
        val_pos = concat(val_pos, teta);
      );
    );
  );
  val_pos;}

\\Returns the possible traces when there is only one root
find_trace_one_root(ell, p) =
  local(a); a = 2*Fp_sqrt(p%ell,ell); [a, ell - a];

\\Returns the possible traces when there are l + 1 roots
find_trace_lp1_roots(ell, p) =
  local(a); a = Mod(Fp_sqrt(p%ell,ell),ell^2); a=lift(p/a+a); [a, ell^2 - a];

\\Returns the trace modulo ell^k when ell is an Elkies prime
find_trace_Elkies_power(E, ell, k, meqn, meqntype, g, tr, p) =
{ local(Eb, Ib, tmp, Ec, p1, kpoly, lambda, mpoly);

  Eb = E;
  Ib = 'x;
  if (VERBOSE, print1("Compute trace mod ", ell));
  if (meqntype == "C",
    tmp = find_isogenous_from_canonical(E, ell, meqn, g, p);
    Ec = tmp[1];
    p1 = tmp[2];
    tmp = find_kernel(E, ell, Ec, p1, p);
    kpoly = tmp[1],
  if (meqntype == "A",
    tmp = find_isogenous_from_Atkin(E, ell, meqn, g, p);
    if (!#tmp, return ([]));
    Ec = tmp[1];
    kpoly = tmp[2]));
  lambda = find_eigen_value(E, ell, kpoly, p, tr);
  if (EARLY_ABORT && (p + 1 - lambda - p/lambda)%ell == 0,
    return( [(lambda + p/lambda)%ell] )
  );
  for (cnt = 2, k,
    if (VERBOSE, print1(", ", ell^cnt));
    mpoly = subst(meqn, 'J, Ec.j);
    tmp = find_kernel_power(Eb, Ec, ell, meqn, meqntype, mpoly, kpoly, Ib, p);
    if (!#tmp, return ([]));
    lambda = find_eigen_value_power(E, ell, cnt, tmp[2], lambda, p);
    Eb = Ec;
    Ec    = tmp[1];
    kpoly = tmp[3];
    Ib    = tmp[4]
  ); [(lambda + p/lambda) % (ell^k)];}

\\Returns [ell^k, v] where v is a vector containing the possible values of the
\\trace modulo ell^k: [], [t] or [t1,...,td]
find_trace(E, ell, k = 1, nb, p) =
{ local(kt, tmp, g, tr, meqn, meqntype);
  if (VERBOSE, print1("\nProcess prime ", ell, ".\tType: "));
  meqn = list_to_pol(modular_eqn[nb]);
  meqntype = modular_eqn_type[nb];
  kt = k;
  tmp = study_modular_eqn(ell, lift(subst(meqn, 'J, E.j)), p);
  g = tmp[2];
/* If l is an Elkies prime, search for a factor of the l-division polynomial.
 * Then deduce the trace by looking for eigenvalues of the Frobenius by
 * computing modulo this factor */
  if (tmp[1] == "E",
    if (VERBOSE, print1("Elkies.\t "));
    tr = find_trace_Elkies_power(E, ell, k, meqn, meqntype, g, [], p);
    if (!#tr, return ([]));
  );
  if (tmp[1] == "1",
    if (VERBOSE, print1("One root.\t "));
    tr = find_trace_one_root(ell, p);
    kt = 1;
    if (ell < BOUND_ONE_ROOT,
      tr = find_trace_Elkies_power(E, ell, 1, meqn, meqntype, g, tr, p);
      if (!#tr, return ([]));
    ,
      if (VERBOSE, print1("Compute possible values for the trace"));
    );
  );
  if (tmp[1] == "2",
    if (VERBOSE, print1("l+1 roots.\t Compute possible values for the trace"));
    tr = find_trace_lp1_roots(ell, p);
    kt = 2
  );
  if (tmp[1] == "A",
    if (VERBOSE, print1("Atkin.\t Compute possible values for the trace"));
    tr = find_trace_Atkin(ell, tmp[3], p);
    kt = 1;
  );
  if (tmp[1] == "P",
    if (VERBOSE, print1("Pathological."));
    tr = []
  ); [ell^kt, tr];}

/****************************************************************************/
/*                              MATCH AND SORT                              */
/****************************************************************************/

cost_without_precomp(ell, compile_atkin) =
{ local(fact, res, j, P, E);

  fact = factor(ell); P = fact[,1]; E = fact[,2];
  res = j = 1;
  for (i = 1, #P,
    while ( P[i]^E[i] != compile_atkin[j][1],
      j++; if (j > #compile_atkin, return("inf"))
    );
    res *= #compile_atkin[j][2];
  ); res;}

\\ assume cost_vec precomputed
cost(ell)=
{ local(res, fact, P);
  fact = factor(ell); P = fact[,1];
  res = prod(i = 1, #P, cost_vec[ P[i] ]);
  if (!res, "inf", res);}

champion(compile_atkin) =
{ local(k, i, n, B, Bp, i1, i2, b, costb, cost_vec, res);

  k = #compile_atkin;
  cost_vec = vector(compile_atkin[k][1]);
  for (i = 1, k,
    n = compile_atkin[i][1];
    cost_vec[ if(isprime(n), n, round(sqrt(n))) ] = #compile_atkin[i][2]
  );
  n = 2;
  B  = vector(2^k); B[1] = 1;
  Bp = B;
  Bp[2] = compile_atkin[1][1];
  for (j = 2, k,
    i = 1; i1 = 2; i2 = 1;
    while(i1 <= n,
      if(Bp[i1] < compile_atkin[j][1] * Bp[i2],
        b = Bp[i1];  i1++
      ,
        b = compile_atkin[j][1] * Bp[i2]; i2++
      );
      costb = cost(b);
      while(costb < cost(B[i]), i--);
      B[i++] = b;
    );
    while(i2 <= n,
      b = compile_atkin[j][1] * Bp[i2];  i2++;
      costb = cost(b);
      while(costb < cost(B[i]), i--);
      B[i++] = b;
    );
    n = i;
    for (i = 1, n, Bp[i] = B[i])
  );
  res = [];
  for (i = 1, 2^k,
    if (B[i] > 1, res = concat(res, [[B[i], cost(B[i])]]))
  ); res;}

\\ A partition of compile_atkin in baby and giant is represented as the binary
\\ developpement of an integer; if the i-th bit is 1, the i-th prime in
\\ compile-atkin is a baby. The optimum is obtained when the ratio between
\\ the number of possibilities for traces modulo giants (p_g) and babies (p_b)
\\ is near 3/4.
separation(cnt) =
{ local(k, best_i, best_r, P, p_b, r, v);
  k = #cnt; P = prod(j=1, k, cnt[j]); \\ p_b * p_g = P is constant
  best_i = 0;
  best_r = 3/4;
  for (i = 1, (1<<k)-2, \\ scan all possibilities
    v = vecextract(cnt, i);
    p_b = prod(j=1,#v, v[j]);
    r = abs(p_b^2/P - 3/4); \\ |p_b/p_g - 3/4|
    if (!r, return(i));
    if (r < best_r, best_i = i; best_r = r);
  ); best_i;}

\\ chinese(Mod(a,A), Mod(b,B)), A,B coprime
crt(A, a, B, b) =
{ local(u, v, M);
  u = A * (1/A % B);
  v = B * (1/B % A); M = A * B;
  [M, (u * b + v * a) % M];
}

global(global_P);
\\ x vector defined modulo P (= global_P), y vector modulo q, (q,P) = 1
\\ return the vector mod q P congruent to x (resp. y) mod P (resp. q).
\\ update global_P ( <-- qP )
multiple_crt(x, y, q) =
{ local(t, k, a1, a2);

  t = vector(#x * #y);
  a1 = global_P * (1/global_P % q);
  a2 = q * (1/q % global_P);
  k = 0; global_P *= q;
  for (i = 1, #x,
    for (j = 1, #y, t[k++] = (a1*y[j] + a2*x[i]) % global_P)
  ); t; }

\\ update global_P
possible_traces(C) =
{ local(v);
  global_P = C[1][1]; v = C[1][2];
  for (i=2, #C, v = multiple_crt(v, C[i][2], C[i][1]));
  v; }

\\ u = trace_elkies, Mu = prod_elkies
match_and_sort(compile_atkin, Mu, u, E, p) =
{ local(baby, giant, Mb, Mg, den, Sg, dec_inf, div, P, Pb,
        Pg, point, diff, pre, d, table, table_ind, r, s, card, best_i, k);
  k = #compile_atkin;
  if (!k, \\no Atkin prime: Mu >= 4*sqrt(p).
    card = p + 1 - u;
    return( choose_card([card, card + Mu], E, p))
  );
  if (k == 1, \\only one Atkin prime, check the cardinality with random points
    s = Mod(u, Mu);
    r = compile_atkin[1];
    card = vector(#r[2], i, p + 1 - centerlift(chinese(s, Mod(r[2][i],r[1]))));
    return ( choose_card(card, E, p) );
  );
  best_i = separation( vector(k, j, #compile_atkin[j][2]) );
  giant= vecextract(compile_atkin, (1<<k)-1-best_i);
  baby = vecextract(compile_atkin, best_i);
  giant = possible_traces(giant); Mg = global_P;
  baby  = possible_traces(baby);  Mb = global_P;
  \\ the variant from Lercier's thesis, section 11.2.3
  den = 1/(Mu * Mg) % Mb;
  for (i = 1, #baby, baby[i] = (baby[i] - u) * den % Mb);
  den = 1/(Mu * Mb) % Mg; Sg = (- u*den) % Mg;
  for (i = 1, #giant, giant[i] = giant[i] * den % Mg);
  dec_inf = ceil(-Mb/2 - Sg*Mb/Mg);
  div = (dec_inf \ Mb) * Mb;
  for (i = 1, #baby,
    baby[i] += div;
    if(baby[i] < dec_inf, baby[i] += Mb)
  );
  P = find_pt_aff(E, p);
  point = ellpow(E, P, Mu);
  Pb = ellpow(E, point, Mg);
  Pg = ellpow(E, point, Mb);

  \\Some precomputations
  baby = vecsort(baby);
  diff = ZV_sort_uniq( vector(#baby-1, i, baby[i+1]-baby[i]) );
  pre = vector(#diff); pre[1] = ellpow(E, Pb, diff[1]);
  \\ what we'd _really_ want here is a hashtable diff[i] -> pre[i]
  for (i = 2, #diff,
    d = diff[i] - diff[i-1];
    pre[i] = elladd(E, pre[i-1], ellpow(E,Pb,d))
  );

  \\Now we compute the table of babies, this table contains only the
  \\lifted x-coordinate of the points in order to use less memory
  table = vector(#baby);
  \\ (p+1 - u - Mu*Mb*Sg) P - (baby[1]) Pb
  point = ellpow(E, P, p + 1 - u - Mu * (Sg*Mb + Mg*baby[1]));
  table[1] = lift(point[1]);
  for (i = 2, #baby,
    d = baby[i] - baby[i-1];
    point = ellsub(E, point, pre[ ZV_search(diff, d) ]);
    table[i] = lift(point[1]);
  );

  \\ Same precomputation for giants
  giant = vecsort(giant);
  diff = ZV_sort_uniq( vector(#giant-1, i, giant[i+1]-giant[i]) );
  pre = vector(#diff); pre[1] = ellpow(E, Pg, diff[1]);
  for (i = 2, #diff,
    d = diff[i] - diff[i-1];
    pre[i] = elladd(E, pre[i-1], ellpow(E,Pg,d))
  );

  \\ Look for a collision among the x-coordinates
  table_ind = vecsort(table,,1);
  table = vecextract(table, table_ind);
  point = ellpow(E, Pg, giant[1]);
  for (i = 1, #giant,
    s = ZV_search(table, lift(point[1]));
    if (s, s = table_ind[s]; r = i; break);
    d = giant[i+1] - giant[i]; \\ error if i = #giant ==> no match
    point = elladd(E, point, pre[ ZV_search(diff, d) ]);
  );
  card = p + 1 - u - Sg * Mu * Mb - Mu * (giant[r] * Mb +  Mg * baby[s]);
  choose_card([card, card + 2 * Mu * Mb * giant[r]], E, p);
}

/* E is an elliptic curve defined over Z or over Fp in ellinit format, defined
 * by the equation E: y^2 + a1*x*y + a2*y = x^3 + a2*x^2 + a4*x + a6
 * p is a prime number
 * set VERBOSE to have information on the computation process
 * set EARLY_ABORT to stop whenever a small factor of the order is detected.
 *   Useful when searching for a good curve for cryptographic applications */
ellsea(E, p, verbose = 0, early_abort = 0) =
{ local(change, Etmp, power_max, i, xp, tr, bound, product,
  l, compile_atkin, nb, trace_mod, M, bound_bsgs,
  growth_factor, best_champ, fact, bound_champ, champ, compile_atkin_tmp, j, T, lp);
  local(x, J, y, z);
  VERBOSE     = verbose;
  EARLY_ABORT = early_abort;
  if (E.a1 || E.a2 || E.a3,
    Etmp = [0,0,0,-27*lift(E.c4), -54*lift(E.c6)];
    change = 1
  ,
    Etmp = [0,0,0,lift(E.a4), lift(E.a6)];
    change = 0;
  );
  Etmp = ellinit(Etmp * Mod(1,p));
  BOUND_ONE_ROOT = 23;

  power_max = vector(47, i, 1);
  \\ sets which power of ell to investigate. Depends on the size of p
  power_max[primepi(3)] = 3;
  power_max[primepi(5)] = 2;
  power_max[primepi(7)] = 2;
  lp = #binary(p);
  \\Uses ellap if p is small
  if (lp < 63, return(p + 1 - ellap(Etmp, p)));
  if (lp > 160, power_max[primepi(3)] = 4);
  if (lp > 260,
    power_max[primepi(5)]  = 3;
    power_max[primepi(11)] = 2;
    power_max[primepi(13)] = 2
  );
  if (lp > 350, power_max[primepi(3)] = 5);
  if (lp > 390, power_max[primepi(7)] = 3);

  if (VERBOSE,
    print1("Computing the order of the elliptic curve\nE:  ");
    if (change,
      print("y^2 + ", lift(E.a1), "*x*y + ", lift(E.a3), "*y = x^3 + ",
            lift(E.a2), "*x^2 + ",lift(E.a4), "*x + ", lift(E.a6));
      print("which has same order as ")
    );
    print("y^2 = x^3 + ", lift(Etmp.a4), "*x + ", lift(Etmp.a6));
    print("defined over the prime field of order ", p, " (", lp, "bits)")
  );
  i = CM_CardEFp(Etmp, p); if (i, return(i));

  x = 'x; J = 'J; y = 'y; z = 'z;

  T = x^3 + lift(Etmp.a4) * x + lift(Etmp.a6);
  xp = FpXQ_pow(x, p, T, p) - x;
  \\Fisrt compute the trace modulo 2
  if (poldegree(FpX_gcd(T, xp, p)) > 0,
    tr = [2, 0];
    if (EARLY_ABORT,
      if (VERBOSE, print("\nAborting: the number of points is divisible by 2"));
      return(0)
    ),
    tr = [2, 1];
  );
  \\Product is the product of the primes we use, the computation
  \\stops if product > 4*sqrt(p).
  \\l is the current prime,
  \\compile_atkin is a vector containing informations about Atkin primes,
  \\informations about Elkies primes lies in tr.
  bound = 4*sqrt(p);
  product = 2;
  l = 3;
  compile_atkin = [];
  while (product <= bound,
    nb = primepi(l);
    trace_mod = find_trace(Etmp, l, power_max[nb], nb, p);
    if (!#trace_mod, next);
    if (#trace_mod[2] == 1,
      if (EARLY_ABORT && (p+1 - trace_mod[2][1]) % l == 0,
        if (VERBOSE, print("\nAborting: #E(Fp) divisible by ", l));
        return(0)
      );
      tr = crt(trace_mod[1], trace_mod[2][1], tr[1], tr[2])
    , \\ else compute possible values for the trace using Atkin method.
      compile_atkin = concat(compile_atkin, [trace_mod]);
    );
\\We can now increase the product with this prime and go to the next prime
    if (#trace_mod[2], product *= trace_mod[1]);
    l = nextprime(l+1);
    if (l > MAXL,
      if (VERBOSE,
        print();
        print("Warning: no more modular polynomials available!");
        print("Match and sort may be very long : it remains ",
              prod(i = 1, #compile_atkin, #compile_atkin[i][2]),
              "possibilities for the trace")
      );
      return( match_and_sort(compile_atkin, tr[1], tr[2], Etmp, p) );
    );
  );
  M = 1000000;
  if(lp <= 160, bound_bsgs = 1.048^lp/9 * M,
  if(lp <= 192, bound_bsgs = 1.052^lp/16.65 * M,
  if(lp <= 306, bound_bsgs = 1.035^lp * 1.35 * M,
                bound_bsgs = 50000*M)));
  growth_factor = 1.26;
  best_champ = [prod(i = 1, #compile_atkin, compile_atkin[i][1]),
                prod(i = 1, #compile_atkin, #compile_atkin[i][2])];
  \\If the number of possible traces is too large, we treat a new prime
  if (VERBOSE && best_champ[2] >= bound_bsgs,
    print("\n\nToo many possibilities for the trace: ", best_champ[2],
          ". Look for new primes")
  );
  while (best_champ[2] >= bound_bsgs,
    trace_mod = find_trace(Etmp, l, 1, primepi(l), p);
    if (!#trace_mod, next);
    if (#trace_mod[2] == 1,
      tr = crt(trace_mod[1], trace_mod[2][1], tr[1], tr[2])
    ,
      compile_atkin = concat(compile_atkin, [trace_mod]);
    );
    if ((best_champ[2]/bound_bsgs < l^2/25 && #trace_mod[2]<=2)
      || best_champ[2]/bound_bsgs < 5 || l == MAXL,
      l = nextprime(l + 1);
      bound_bsgs *= growth_factor,
      \\Let us now treat an other prime if we are too far from the bound_bsgs
      l = nextprime(l + 1);
      bound_bsgs *= growth_factor;
      trace_mod = find_trace(Etmp, l, 1, primepi(l), p);
      if (!#trace_mod, next);
      if (#trace_mod[2] == 1,
        tr = crt(trace_mod[1], trace_mod[2][1], tr[1], tr[2])
      ,
        compile_atkin = concat(compile_atkin, [trace_mod]);
      );
      l = nextprime(l + 1);
      bound_bsgs *= growth_factor;
    );
    bound_champ = 4*sqrt(p) / tr[2];
    champ = champion(compile_atkin);
    best_champ = champ[#champ];
    for (i = 1, #champ,
      if (champ[i][1] > bound_champ && champ[i][2] < best_champ[2],
        best_champ = champ[i]
      );
    );
    \\best_champ is the champion of lowest cost among champions less than the
    \\required bound. We now recover the corresponding sets of traces.
    compile_atkin_tmp = [];
    fact = factor(best_champ[1]);
    i = 1;
    j = 1;
    while (j <= #compile_atkin && i <= matsize(fact)[1],
      if (fact[i, 1]^fact[i, 2] == compile_atkin[j][1],
        compile_atkin_tmp = concat(compile_atkin_tmp, [compile_atkin[j]]);
        i++;
      ); j++
    );
    compile_atkin  = compile_atkin_tmp;

    if (l > MAXL && best_champ[2] > bound_bsgs,
      if(VERBOSE,
        print("\nWarning: no more modular polynomials available, ",
          "match and sort may be very long : it remains ",
          best_champ[2], " possibilities for the trace"));
      break;
    );
  );
  if (VERBOSE,
    print("\n\nComputation of traces done. Entering match and sort algorithm.")
  );
  match_and_sort(compile_atkin, tr[1], tr[2], Etmp, p); }

\\Ensures that E is not supersingular, anomalous and fullfills the MOV condition
bad_curve(p, res) =
{local(d);
  if (p == res, return(1)); \\ anomalous
  if (p == res + 1, return(2)); \\ supersingular
  d = 1;
  for (i = 0, 30,
    d = (d * p) % res; if (d == 1, return(3)) \\ MOV condition not satisfied
  ); 0; }

is_singular(A, B, p) = (4*A^3 + 27*B^2) % p == 0;

\\Finds a curve whose order is a prime of size approximately lg bits
ellcrypto(lg) =
{local(p, nbessai, E, res, a4, a6, t = 1<<(lg-1));
  p = nextprime(t + random(t));
  nbessai = 0;
  while (1, print1("*"); nbessai++;
    a4 = random(p);
    a6 = random(p);
    if (is_singular(a4,a6,p), next);
    E = ellinit([0, 0, 0, a4, a6] * Mod(1,p));
    res = ellsea(E,p,0,1);
    if (isprime(res) && !bad_curve(p, res), break)
  );
  [lift(E.a4), lift(E.a6), res]; }

addhelp(ellsea,"ellsea(E,p,{flag1=0},{flag2=0}): returns the cardinality of the group of rational points E(Fp) where E is an elliptic curve in ellinit format defined over Z or Fp by the equation E: y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6. If flag1 is equal to 1 the program displays information about the computation process. If flag2 is set to 1 the early abort technique is used and the computation is interrupted as soon as a small divisor of the order is detected.");

addhelp(ellcrypto,"ellcrypto(lg): returns a curve whose order is a prime of size approximately lg bits. Checks also that the curve is not supersingular, anomalous and that it fulfills the MOV condition.");
