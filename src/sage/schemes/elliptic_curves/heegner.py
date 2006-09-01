"""nodoctest
Heegner points
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import ellcurve

"""



Heegner points package.
> E :=

> E := EllipticCurve([0,0,1,-1,0]);   // 37A
> A := [-D : D in [1..100] | HeegnerHypothesis(E,-D)]; A;
[ -3, -4, -7, -11, -40, -47, -67, -71, -83, -84, -95 ]
P := [* *]; for D in A do time Y := HeegnerPoint(E,D); Append(~P,Y); end for;
Time: 0.040
Time: 0.030
Time: 0.040
Time: 0.050
Time: 0.090
Time: 0.180
Time: 0.050
Time: 0.250
Time: 0.120
Time: 0.140
WARNING: Precision overflow or point is 0.  Use L-series to decide or try
AssertAttribute(FldPr,"Precision",n) or tweaking parameters.
Time: 0.270

> for i in [1..10] do Recognize(P[i]); end for;
(-1 : 0 : 1)
(1 : -1 : 1)
(0 : 0 : 1)
(0 : -1 : 1)
(1 : -1 : 1)
(0 : 0 : 1)
(6 : -15 : 1)
(0 : -1 : 1)
(0 : 0 : 1)
(0 : 0 : 1)


* BUG!?

> E := EllipticCurve([ 0, 1, 1, -26, -61 ]);
> time P := HeegnerPoint(E,-107); P;
> (3.665647323987715053919669397 + 6.116738926164336860764880743 E-25*i : -0.49999999999999999999999867820523179988 + 9.6625393082244036119397240285848446462*i : 1)
> IdentifyPoint(P);
> false 6938055/1892723
// This time it doesn't identify it as a rational point.  In fact,
// the coordinates of P generate the Hilbert class field of Q(sqrt(-107)),
// which is really weird because we just took a trace to get down
// from the Hilbert class field.  What is going on!
// Note: Pete Green's Heegner point package gives the same wrong answer.

"""

"""
intrinsic tau_to_z(E::CrvEll, tau::FldQuadElt, qexp_prec::RngIntElt) -> FldPrElt
{}
   C<i> := ComplexField();
   tau := CoerceIntoC(tau);
   if Imaginary(tau) lt 0 then
      tau := ComplexConjugate(tau);
   end if;
   Pi := Pi(C);
   f,an_over_n := qExpansion(E,qexp_prec);
   e := Exp(2*Pi*i*tau);
   V := VectorSpace(C,qexp_prec-1);
   z := InnerProduct(V!an_over_n, V![e^n : n in [1..qexp_prec-1]]);
   return z;
end intrinsic;

intrinsic ideal_to_tau(a::RngQuadIdl, NN::RngQuadIdl, K::FldQuad) -> FldQuadElt
{Given a nonzero integral fractional ideal a of O_K and an integral
ideal NN such that O_K/NN = Z/NZ, compute tau in K such that
(E(tau),C(tau)) = (C/a^(-1), N^(-1)*a^(-1)/a^(-1)).}
   // reduce ideal: gives equivalent point, since just
   // multiplies a be a scalar.
    a := Ideal(Reduction(QuadraticForm(a)));

   // compute I1 and I2 in point (I1, I2/I1)
   I1 := a^(-1);
   I2 := I1*NN^(-1);

   // Find a basis for I2 of the form omega_1, omega_2/N for
   // some basis om1 and om2 of I1.   Such a basis
   // exists because I2/I1 = Z/NZ.  We compute it using
   // Smith Normal Form.
   V := VectorSpace(RationalField(),2);
   W := VectorSpaceWithBasis([V | Eltseq(b) : b in Basis(I2)]);
   B1 := Basis(I1);
   A := MatrixAlgebra(Integers(),2)!&cat[Coordinates(W,V!Eltseq(b)) : b in B1];

   // The smith form S of A, together with unimodular
   // matrices P and Q such that P * A * Q = S
   S, P, Q := SmithForm(A);
   // Apply P to the basis B1 to get om1, om2.
   om1 := P[1,1]*B1[1] + P[1,2]*B1[2];
   om2 := P[2,1]*B1[1] + P[2,2]*B1[2];

   // Finally, the point z is 1/om2;
   z := K!(1/om2);
   return z;
end intrinsic;

intrinsic ideal_to_tau(f::QuadBinElt) -> FldQuadElt
{}
    D := Discriminant(f);
    K := QuadraticField(D);
    sqD := K.1;
    if D mod 4 eq 0 then
       sqD := 2*sqD;
    end if;
    z := (-f[2] + sqD)/(2*f[1]);
    return z;
end intrinsic;

intrinsic levelNforms(N::RngIntElt, D::RngIntElt) -> SeqEnum
{As in Step 11 of Stephens's paper.}
    beta := Integers()!Sqrt(Integers(4*N)!D);
    Q := QuadraticForms(D);
    G,m := ClassGroup(Q);
    F := Sort([Reduction(m(g)) : g in G]);
    levelforms := [Q!1 : i in [1..#F]];
    known     := [false : i in [1..#F]];
    found := 0;
    A := 1;
    N4 := 4*N;
    alpha := 0;
    while found lt #F do
       B := beta + N4*alpha;
       C := (B^2 - D)/(N4*A);
       if Denominator(C) eq 1 then
          test := Q![Integers()|A*N,B,C];
          i := Index(F, Reduction(test));
          if not known[i] then
             found := found + 1;
             known[i] := true;
             levelforms[i] := test;
          end if;
       end if;
       alpha := alpha + 1;
       if alpha gt A+2 then
          alpha := 0;
          A := A + 1;
       end if;
    end while;

    return levelforms;
end intrinsic;

intrinsic HeegnerPointUseIdeals(E::CrvEll, D::RngIntElt,
                  q_prec::RngIntElt, wp_prec::RngIntElt) -> SeqEnum, CrvEll
{}
   require IsFundamentalDiscriminant(D) : "Argument 2 must be a fundamental discriminant.";
   require D lt 0 : "Argument 2 must be negative.";
   N := Conductor(E);
   require GCD(D, 4*N) eq 1 : "D and 4*N must be coprime";
   require Type(BaseRing(E)) eq FldRat : "Argument 1 must be defined over the rational field.";

   K := QuadraticField(D);
   OK := MaximalOrder(K);
   // Find an ideal NN with O/NN = Z/NZ.
   NN := 1*OK;
   for p in Factorization(N) do
      Q := Factorization(p[1]*OK)[1][1];
      if Degree(Q) ne 1 then
         require false : "Each prime dividing argument 3 must split completely in Q(sqrt(D)).";
      end if;
      NN := NN * Q^p[2];
   end for;
   G, m := ClassGroup(OK);
   vprint Heegner :  "Class number =", #G;
   taus := [ideal_to_tau(m(g), NN, K) : g in G];
   vprint Heegner :  "taus =", taus;
   // taus2 := [ideal_to_tau_2(m(g), NN, K) : g in G];
   // vprint Heegner :  "taus2 =", taus2;
   z := &+[tau_to_z(E, tau, q_prec) : tau in taus];
   vprint Heegner :  "z    =", z;
   // z2 := &+[tau_to_z(E, tau, q_prec) : tau in taus2];
   // vprint Heegner :  "z2    =", z2;
   p := z_to_point(E, z, wp_prec) ;
   vprint Heegner :  "p = ", p;
   // p2 := z_to_point(E, z2, wp_prec) ;
   // vprint Heegner :  "p2 = ", p2;
   return p;
end intrinsic;

intrinsic HeegnerHypothesis(E::CrvEll, D::RngIntElt) -> BoolElt
{True if D satisfies the Heegner Hypotheses, so there is a Heegner
point attached to discriminant D.}
   if not IsFundamentalDiscriminant(D) or D ge 0 then
      return false;
   end if;
   N := Conductor(E);
   return IsSquare(Integers(4*N)!D) and GCD(N,D) eq 1;
end intrinsic;

intrinsic RootNumber(E::CrvEll) -> RngIntElt
{The sign in the functional equation for E, so the analytic rank of E is even
if and only if the root number is +1.}
   require Type(BaseRing(E)) eq FldRat : "Argument 1 must be defined over the rational field.";
   eps := -1;
/*
   for p in Factorization(Conductor(E)) do
   end for;
*/
end intrinsic;

intrinsic MyHeegnerPoint(E::CrvEll, D::RngIntElt :
                        q_prec:= -1, wp_prec := 40, t := 6) -> PtEll, PtEll
{The Heegner point on E corresponding to D (over the complex numbers), and the a corresponding Q-rational point on
 either E or on the twist E^D (this point might be a trace from Q(sqrt(D)) to Q).  Returns 0 if the Heegner point is torsion.}
   require E eq MinimalModel(E) : "Argument 1 must be minimal.";
   require HeegnerHypothesis(E,D) : "The Heegner Hypotheses must hold.";
   EE := EllipticCurve([ComplexField()|a : a in aInvariants(E)]);
   tiny := 10^(-6);

   /*
   if Abs(LOne(E)) lt tiny then  // E has rank > 0
      if Abs(LPrimeOne(E)) lt tiny or Abs(LOneChi(E,D)) lt tiny then
         return EE!0, E!0;
      end if;
   else      // E has rank 0, so need twist to have rank at most 1
      if LPrimeOneChi(E,D) lt tiny then
         return EE!0, E!0;
      end if;
   end if;
   */

   N := Conductor(E);
   require Type(BaseRing(E)) eq FldRat : "Argument 1 must be defined over the rational field.";

   if q_prec le 0 then
      q_prec := 2*N*Ceiling(Log(-D)/Log(5)) + 100;
      vprint Heegner : "Using q-expansion to precision", q_prec;
   end if;

   forms := levelNforms(N, D);
   taus := [ideal_to_tau(f) : f in forms];
   vprint Heegner :  "h_K   =", #taus;
   vprint Heegner :  "tau's =", taus;
   z := &+[tau_to_z(E, tau, q_prec) : tau in taus];
   vprint Heegner :  "z     =", z;
 p := z_to_point(E, z, wp_prec : t := t) ;
   vprint Heegner :  "p     =", p;
   if not IsPoint(EE,p) then
      error "*WARNING*: Precision Overflow!";
   end if;
   p := EE!p;

   if Abs(LOne(E)) lt tiny then  // E has rank >=1
      if Max(Abs(Imaginary(p[1])), Abs(Imaginary(p[2]))) gt 10^(-2) then
         q := TracePoint(p);
      else
         q := p;
      end if;
      t, Q := Recognize(q);
   else   // E has rank 0
      phi := TwistMap(E, D);
      q0 := phi(p);
print "q0 = ", q0;
      if Max(Abs(Imaginary(q0[1])), Abs(Imaginary(q0[2]))) gt 10^(-2) then
         q := TracePoint(q0);
      else
         q := q0;
      end if;
print "q = ", q;
      t, Q := Recognize(q);
   end if;
   if not t then
      print "WARNING: Couldn't correctly recognize point for D=", D;
   end if;

   return EE!p, Q;
end intrinsic;

intrinsic ApproximateRational(x::FldPrElt) -> FldRatElt
{Try to approximate x using continued fractions.}
   vprint Heegner: "Approximating x = ", x;
   cf := ContinuedFraction(Real(x));
   vprint Heegner: "Continued fraction = ", cf;
   for i in [2..#cf] do
      if cf[i] gt 1000 then
         v := Convergents([cf[j] : j in [1..i-1]]);
         return v[1,1]/v[2,1];
      end if;
   end for;
   v := Convergents(cf);
   p := v[1,1]/v[2,1];
   return p;
end intrinsic

intrinsic Recognize(P::PtEll) ->  BoolElt, PtEll
{}
   EC := Scheme(Parent(P));
   E := EllipticCurve([Round(Real(a)) : a in aInvariants(EC)]);
   if P eq EC!0 then
      return true, E!0;
   end if;
   x := ApproximateRational(P[1]);
   R<z> := PolynomialRing(Rationals());
   roots := Roots(Evaluate(DefiningPolynomial(E),[x,z,1]));
   if #roots eq 0 then
      vprint Heegner : "Unable to correctly identify rational point!";
      print "P = ", P;
      return false, x;
   end if;
   if #roots eq 1 then
      return true, E![x,roots[1][1]];
   end if;
   if Abs(roots[1][1]-P[2]) lt Abs(roots[2][1]-P[2]) then
      return true, E![x,roots[1][1]];
   else
      return true, E![x,roots[2][1]];
   end if;
end intrinsic;

intrinsic TwistMap(E::CrvEll, D::RngIntElt) -> Map
{The map phi: E_K --> F_K, where F is the quadratic twist of E
by K = Q(sqrt(D)), and the model for F is Z-minimal.  The second
argument is the base extension of phi to the complex numbers.}
   F := MinimalModel(QuadraticTwist(E, D));
   K := QuadraticField(D);
   C := ComplexField();
   E_K := BaseExtend(E,K);
   F_K := BaseExtend(F,K);
   phi_K := Isomorphism(E_K, F_K);
   E_C := EllipticCurve([CoerceIntoC(a) : a in aInvariants(E)]);
   F_C := EllipticCurve([CoerceIntoC(a) : a in aInvariants(F)]);
   r,s,t,u := Explode([CoerceIntoC(a) : a in IsomorphismData(phi_K)]);
   function doit(P, FF)
      Q := [u^2*P[1]+r, u^3*P[2]+s*u^2*P[1]+t];
      if not IsPoint(FF,Q) then
         print "*WARNING*: Possible precision problem!";
         return FF!0;
      end if;
      return FF!Q;
   end function;
   phi_C := hom< E_C -> F_C | P :-> doit(P, F_C)> ;
   return phi_C, [r,s,t,u];
end intrinsic;

intrinsic TracePoint(P::PtEll) -> PtEll
{P + Pbar}
   if P[3] eq 0 then
      return P;
   end if;
   E := Scheme(Parent(P));
   Pbar := E![ComplexConjugate(P[1]), ComplexConjugate(P[2])];
   return P + Pbar;
end intrinsic;


intrinsic Rank0ShaBound(E::CrvEll)
{Compute an upper bound on Sha for the rank 0 elliptic curve E.}

   A := [-D : D in [1..200] | HeegnerHypothesis(E,-D)];
   print "First few Heegner primes: ", A;
   for D in [A[1]] do
      F := MinimalModel(QuadraticTwist(E,D));
      time P := HeegnerPoint(E,D : wp_prec := 100);
      phi := TwistMap(E,D);
      T := TracePoint(phi(P));
print "T = ", T;
      t, TQ := Recognize(T);
print "TQ = ", TQ;
      if not t then
         print "Could not recognize T!";
      else
         h := Height(TQ);
         R := Regulator(F);
         s := Sqrt(h/R);
         print "s = ", s;
         a := ApproximateRational(s);
         print "a = ", a;
      end if;
   end for;

end intrinsic;
"""
