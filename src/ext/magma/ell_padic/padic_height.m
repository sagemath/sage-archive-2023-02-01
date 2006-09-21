freeze;

/*****************************************************************

 padic_height.m - Computation of the canonical global p-adic
                  cyclotomic height pairing.

 (c) 2004, William Stein, was@math.harvard.edu

 ****************************************************************/

/*****************************************************************
1. You *must* put a copy of kedlaya.m in the same directory as this file.
   The following line imports five functions from the kedlaya.m
   that is distributed with MAGMA (as package/Geometry/CrvHyp/kedlaya.m)
   The functions FrobYInv, Convert, ReduceA, and ReduceB implement
   arithmetic in Monsky-Washnitzer cohomology, following Kedlaya's
   algorithms.  These functions are used below in the intrinsic E2
   but nowhere else.
 ****************************************************************/

import "kedlaya.m": myXGCD, FrobYInv, Convert, ReduceA, ReduceB;

/*****************************************************************
   2. How do you use this package?

First try the following calculations and see if they work:

> Attach("kedlaya.m");
> Attach("padic_height.m");
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> P := E![0,0];
> h := height_function(E, 5, 20);
> h(P);
-201147061754703 + O(5^21)

Also, compute some regulators:

> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> regulator(E,5,20);
138360994885922 + O(5^21)
> regulator(E,7,20);
4709403600911866 + O(7^20)
> E := EllipticCurve(CremonaDatabase(), "389A");
> regulator(E,7,20);
22975764581280320 + O(7^20)

Now you should know enough to be able to make use of the package.
However, you should probably read through the source code below, To
make this easier, I've put a description of *every* function with a
complete example of usage.  I've also grouped the functions into
logical sections.

****************************************************************/


/*****************************************************************

  Formal groups of elliptic curves.

 ****************************************************************/

R<t> := LaurentSeriesRing(RationalField());

intrinsic formal_w(E::CrvEll, prec::RngIntElt) -> RngSerLaurElt
{
Compute and return the formal power series w(t) defined in Chapter IV, Section 1
of Silverman's AEC book.  This function is useful because Silverman gives
expressions in w(t) for other formal power series related to formal groups.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_w(E,19);
t^3 + t^6 - t^7 + 2*t^9 - 4*t^10 + 2*t^11 + 5*t^12 - 15*t^13 +
    15*t^14 + 9*t^15 - 56*t^16 + 84*t^17 - 14*t^18 + O(t^19)
}
   w := t^3 + O(t^prec);
   a1, a2, a3, a4, a6 := Explode(aInvariants(E));
   for n in [4..prec-1] do
      w := w +  t^n*(a1*Coefficient(w,n-1) +
                a2*Coefficient(w,n-2) +
                a3*&+[Integers()|Coefficient(w,i)*Coefficient(w,n-i) : i in [3..n-3]] +
                a4*&+[Integers()|Coefficient(w,i)*Coefficient(w,n-1-i) : i in [3..n-4]] +
                a6*&+[Integers()|Coefficient(w,i)*Coefficient(w,j)*Coefficient(w,n-i-j) :
                         i in [3..n-3], j in [3..n-3] | n-i-j ge 3]);
   end for;
   return w;
end intrinsic;

intrinsic formal_x(E::CrvEll, prec::RngIntElt) -> RngSerLaurElt
{
Compute and return the formal group power series defined by the
function x on E.

EXAMPLES:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_x(E,15);
t^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 -
    3*t^10 + 20*t^11 - 30*t^12 + 6*t^13 + 65*t^14 + O(t^15)
}
   return t/formal_w(E,prec+5);
end intrinsic;

intrinsic formal_y(E::CrvEll, prec::RngIntElt) -> RngSerLaurElt
{
Compute and return the formal group power series defined by the
function y on E.

EXAMPLES:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_y(E,15);
-t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9
    - 20*t^10 + 30*t^11 - 6*t^12 - 65*t^13 + 140*t^14 + O(t^15)
}
   return -1/formal_w(E,prec+6);
end intrinsic;


intrinsic formal_omega(E::CrvEll, prec::RngIntElt) -> RngSerLaurElt
{
Compute and return the formal differential in terms of t, i.e., return
the function f(t) such that f(t)dt is the formal differential omega.

ALGORITHM:
    Use that omega = dx/(2*y + a1*x + a3) = f(t) dt.

EXAMPLES:
> formal_omega(E,10);
1 + 2*t^3 - 2*t^4 + 6*t^6 - 12*t^7 + 6*t^8 + 20*t^9 - 60*t^10 +
    60*t^11 + O(t^12)
}

   a1, a2, a3, a4, a6 := Explode(aInvariants(E));
   x := formal_x(E,prec);
   y := formal_y(E,prec);

   return Derivative(x)/(2*y + a1*x + a3);

end intrinsic;

intrinsic formal_log(E::CrvEll, prec::RngIntElt) -> RngSerLaurElt
{
The formal log in terms of t.  This is sometimes written z=z(s).  It is the
integral with respect to t of the formal differential.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_log(E,10);
t + 1/2*t^4 - 2/5*t^5 + 6/7*t^7 - 3/2*t^8 + 2/3*t^9 + 2*t^10 -
    60/11*t^11 + 5*t^12 + O(t^13)
}
   y := formal_y(E,prec);
   x := formal_x(E,prec);
   a1, a2, a3, a4, a6 := Explode(aInvariants(E));
print 1/(2*y + a1*x + a3) * Derivative(x);
   return Integral(1/(2*y + a1*x + a3) * Derivative(x));
end intrinsic;

intrinsic formal_inverse(E::CrvEll, prec::RngIntElt) -> RngSerLaurElt
{
Compute the power series in Z[[t]] corresponding to the inverse
map in the formal group on E.

EXAMPLES:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_inverse(E,10);
-t - t^4 - 2*t^7 + t^8 - 5*t^10 + 6*t^11 - 2*t^12 + O(t^13)
}
   x := formal_x(E,prec);
   y := formal_y(E,prec);
   a1, a2, a3, a4, a6 := Explode(aInvariants(E));
   return x/(y + a1*x + a3);
end intrinsic;

intrinsic formal_group_law(E::CrvEll, prec::RngIntElt) -> RngSerLaurElt
{
Compute and return the formal group law on E as an element of Z[[t,t2]].
The precision must be at least 4.

ALGORITHM:
  Use the formula at the end of Section 1 of Chapter 4 of Silverman AEC,
on page 114.  This is the formula for "z_3" there.  Note though that
the FORMULA IS *WRONG* in Silverman's book!!  The a_1 and a_3 in that
formula should have minus signs in front, but they don't.


EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_group_law(E,4);
t + O(t^4) + (1 - 2*t^3 + 2*t^4 + O(t^6))*t2 + (-3*t^2 + 4*t^3 +
    O(t^5))*t2^2 + (-2*t + 4*t^2 + O(t^4))*t2^3 + O(t2^4)
}
   require prec ge 4: "The precision must be at least 4.";
   t1 := t;
   S<t2> := LaurentSeriesRing(R);
   w := formal_w(E,prec);
   function tsum(n)
      return &+[t2^m * t1^(n-m-1) : m in [0..n-1]];
   end function;
   lam := &+[tsum(n)*Coefficient(w,n) : n in [3..prec-1]];

   w1 := Evaluate(w,t1+O(t^prec));
   nu := w1 - lam*(t1+O(t^prec));
   a1, a2, a3, a4, a6 := Explode(aInvariants(E));
   t3 := -t1 - t2 + \
         (-a1*lam - a3*lam^2 - a2*nu - 2*a4*lam*nu - 3*a6*lam^2*nu)/
         (1 + a2*lam + a4*lam^2 + a6*lam^3) + O(t1^prec) + O(t2^prec);
   inv := formal_inverse(E, prec);
   return Evaluate(inv, t3);
end intrinsic;

intrinsic formal_mult(E::CrvEll, n::RngIntElt, prec::RngIntElt) -> RngSerLaurElt
{
Compute and return the formal power series corresponding to multiplication by n on E,
to precision prec.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
>  formal_mult(E,2,5);
2*t - 7*t^4 + 12*t^5 + 4*t^7 + O(t^8)
>  formal_mult(E,3,5);
3*t - 39*t^4 + 96*t^5 + 234*t^7 + O(t^8)
>  formal_mult(E,7,5);
7*t - 1197*t^4 + 6720*t^5 + 115254*t^7 + O(t^8)
}

   if n eq -1 then
      return formal_inverse(E,prec);
   end if;

   // The following is less direct, but it works more quickly...
   L := FunctionField(E);
   EL := BaseExtend(E,L);
   Q := EL![L.1,L.2];
   x,y := Explode(Eltseq(n*Q));
   a := -x/y;
   xx := formal_x(E,prec);
   yy := formal_y(E,prec);
   phi := hom<Parent(x) -> Parent(xx) | xx,yy>;
   return phi(a);
end intrinsic;

intrinsic formal_divpoly(E::CrvEll, m::RngIntElt, prec::RngIntElt) -> RngSerElt
{
Compute and return the formal division "polynomial" f_m for E to
precision prec.

EXAMPLES:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_divpoly(E,2,5);
2*t^-3 - 3 + 2*t - 2*t^3 + O(t^4)
> formal_divpoly(E,3,5);
3*t^-8 - 12*t^-5 + 6*t^-4 + 9*t^-2 + O(t^-1)
> formal_divpoly(E,7,5);
7*t^-48 - 168*t^-45 - 140*t^-44 + 2750*t^-42 + O(t^-41)
}
   if m eq 1 then
      return R!1;
   end if;
   if m eq -1 then
      return R!(-1);
   end if;
   if m eq 2 then
      return Sqrt(Evaluate(DivisionPolynomial(E,m),formal_x(E,prec)));
   end if;
   assert IsOdd(m);
   return Evaluate(DivisionPolynomial(E,m),formal_x(E,prec));
end intrinsic;

intrinsic weierstrass_p(E::CrvEll, prec::RngIntElt) -> FldPrElt
{
   Compute and return the Weierstrass p-series for E to precision prec.

EXAMPLES:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> weierstrass_p(E,9);
t^-2 + 1/5*t^2 - 1/28*t^4 + 1/75*t^6 - 3/1540*t^8 + O(t^10)
}
   L<w> := PolynomialRing(Rationals());
   R<t> := LaurentSeriesRing(L);
   p := t^(-2);
   c4, c6 := Explode(cInvariants(E));
   a := 4; b := -c4/12; c := -c6/216;
   for n in [1..prec] do
      ptest := p + w*t^n + O(t^(n+1));
      ptest_prime := Derivative(ptest);
      f := ptest_prime^2 - a*ptest^3 - b*ptest - c;
      x := Coefficient(f,-4+n);
      coeff := -Coefficient(x,0)/Coefficient(x,1);
      p := p + coeff*t^n;
   end for;
   return LaurentSeriesRing(Rationals())!(p + O(t^(prec+1)));
end intrinsic;


intrinsic t_val(Q::PtEll)  -> FldRatElt
{
Given a point (x,y) on an elliptic curve, return t = -x/y.  (This function doesn't
do much!)

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> P := E![0,0]; P;
(0 : 0 : 1)
> Q := 5*P;
(1 : 0 : 1)
> Q := 5*P;
> Q;
(1/4 : -5/8 : 1)
> t_val(Q);
2/5
}
   require Q[2] ne 0 : "The denominator y must be nonzero.";
   return -Q[1]/Q[2];
end intrinsic;


/*****************************************************************

    Computing sigma: Integrality algorithm

 ****************************************************************/

intrinsic formal_sigma_in_s2(E::CrvEll, prec::RngIntElt) -> RngSerElt
{
Compute a formal power series for sigma with coefficient polynomials
in the coefficient s_2 of t^3.

NOTE: This is adapted from Christian Wuthrich's PARI code, which is also the
algorithm Mazur-Tate and I found on 2004-06-09.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> formal_sigma_in_s2(E, 10);
t + s2*t^3 + 1/2*t^4 + (1/2*s2^2 - 5/12)*t^5 + 3/2*s2*t^6 +
    (1/6*s2^3 - 73/60*s2 + 103/120)*t^7 + (5/4*s2^2 - 37/24)*t^8
    + (1/24*s2^4 - 121/120*s2^2 + 2791/840*s2 + 1411/2016)*t^9 +
    (7/12*s2^3 - 691/120*s2 + 481/240)*t^10 + (-7/15*s2^3 +
    95/28*s2^2 + 379/150*s2 - 85793/15400)*t^11 + (3/16*s2^4 -
    463/80*s2^2 + 4873/560*s2 + 6977/1344)*t^12 + O(t^13)
}

   R<s2> := FieldOfFractions(PolynomialRing(RationalField()));
   S<z> := LaurentSeriesRing(R);
   w := S!weierstrass_p(E,prec);
   w := &+[Coefficient(w,n)*z^n : n in [0..prec-1]] + O(z^prec);
   sigma_of_z := z*Exp( s2*z^2 - Integral(Integral(w)) + O(z^prec) );
   z_of_s := S!formal_log(E, prec);
   sigma_of_s := Evaluate(sigma_of_z, z_of_s);
   R<s2> := PolynomialRing(RationalField());
   S<t> := PowerSeriesRing(R);
   return S!sigma_of_s;
end intrinsic;

intrinsic solve_for_s2(sigma::RngSerElt, p::RngIntElt) -> .
{
Given a prime p and the formal power series for sigma in terms of s_2,
as output by formal_sigma_in_s2, find the element s_2 of Z_p to as high
a precision as possible by using that the coefficients of sigma must be
p-adic integers.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> sig := formal_sigma_in_s2(E, 40);
> solve_for_s2(sig, 5);
53 + O(5^3)
> print_padic($1);
3*5^0 + 2*5^2 + O(5^3)
}
   a := 0;   // a is what we know about s2 so far.
   n := 0;
   S<t> := Parent(sigma);
   R<s2> := BaseRing(S);
   prec := AbsolutePrecision(sigma);
   for i in [3..prec-1] do
      c := Coefficient(sigma,i);
      d := Coefficients(c);
      m := Max([0] cat [Valuation(Denominator(x),p) : x in d]);
      if m gt 0 then
         f := p^m*c;
         // Now viewing f as a poly in s2, it must be 0 modulo p^m, so
         // as to cancel the denominator.
         R<x> := PolynomialRing(GF(p));
         g := R!f;
         X := Roots(g);
         assert #X le 1;
         if #X eq 1 then
            b := Integers()!X[1][1];
            a +:= b*p^n;
            n +:= 1;
            // Now s2 = a + O(p^(n+1)).
            z := b+ p*s2;
            sigma := &+[Evaluate(Coefficient(sigma,r),z)*t^r : r in [0..prec-1]] + O(t^prec);
         end if;
      end if;
   end for;
   return pAdicField(p,n)!a;
end intrinsic;


intrinsic sigma_using_integrality(E::CrvEll, p::RngIntElt,
				  prec::RngIntElt) -> .
{
Compute the function sigma for E at p using the classical integrality algorithm.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> sigma_using_integrality(E, 5, 40);
O(5^3) + t + O(5^3)*t^2 + (53 + O(5^3))*t^3 - (62 + O(5^3))*t^4 -
    (23 + O(5^3))*t^5 + (17 + O(5^3))*t^6 - (6 + O(5^2))*t^7 -
    (58 + O(5^3))*t^8 - (5 + O(5^2))*t^9 + (11 + O(5^2))*t^10 -
    (2 + O(5))*t^11 + (11 + O(5^2))*t^12 - (1 + O(5))*t^13 +
    O(5)*t^14 - (1 + O(5))*t^15 + (1 + O(5))*t^16 + O(1)*t^17 +
    (2 + O(5))*t^18 + O(1)*t^19 + O(1)*t^20 + O(5^-1)*t^21 +
    O(1)*t^22 + O(5^-1)*t^23 + O(5^-1)*t^24 + O(5^-1)*t^25 +
    O(5^-1)*t^26 + O(5^-2)*t^27 + O(5^-1)*t^28 + O(5^-2)*t^29 +
    O(5^-2)*t^30 + O(5^-3)*t^31 + O(5^-2)*t^32 + O(5^-3)*t^33 +
    O(5^-3)*t^34 + O(5^-3)*t^35 + O(5^-3)*t^36 + O(5^-4)*t^37 +
    O(5^-3)*t^38 + O(5^-4)*t^39 + O(t^40)
}

   sigma := formal_sigma_in_s2(E, prec);
   s2 := solve_for_s2(sigma, p);
   R<t> := PowerSeriesRing(Parent(s2));
   return &+[Evaluate(Coefficient(sigma,n),s2) * t^n : n in [0..prec-1]] + O(t^prec);
end intrinsic;

/**************************************************************************

 Computation of E2(E,omega) and the matrix of absolute Frobenius.

 **************************************************************************/

intrinsic E2(E::CrvEll, p::RngIntElt, prec::RngIntElt) -> .
{
Returns value of E2 on the elliptic curve and the matrix
of frob on omega=dx/y, eta = x*dx/y, where the elliptic
curve is first put in WeiestrassModel form y^2=x^3+ax+b.

INPUT:
   E -- elliptic curve
   p -- prime
   prec -- p-adic precision for computations.
}
   require p ge 5 : "Argument 2 must be at least 5.";
   require IsPrime(p) : "Argument 2 must be prime.";
   require Conductor(E) mod p ne 0 : "Argument 1 must have good reduction at argument 2.";

   c4, c6 := Explode(cInvariants(MinimalModel(E)));
   a4 := - c4/(2^4 * 3);
   a6 := - c6/(2^5 * 3^3);
   assert IsIsomorphic(E, EllipticCurve([a4,a6]));

   n := 1;
   d := 2;
   g := 1;

   /* Set cube = false, so that the following basis is used
      in the computation of Frob on cohomology:

        omega = dx/y  and  eta = x(dx/y)

      This is the basis that Katz uses in his computation of E2.
   */
   cube := false;

   N1 := Ceiling((n*g/2)+Log(p,2*Binomial(2*g,g))) + prec;
   N := N1 + Floor(Log(p,2*N1))+1;
   K := pAdicField(p,N);
   R1 := quo<Integers(K)|p^N>;  // R1 = Z / p^N Z

   L<T> := LaurentSeriesRing(R1);
   P1 := PolynomialRing(L);
   X := P1.1;
   S := quo<P1 | X^3+a4*X+a6 - T^-1>;
   x := S.1;

   W<X> := PolynomialRing(K);
   precs := [X^3+a4*X+a6];
   Append(~precs,Derivative(precs[1]));

   A,B := myXGCD(precs[1],precs[2]);
   // A,B satisfy A*Q+B*Q'=1 where Q is the lift of poly to char 0.
   Append(~precs,A);
   Append(~precs,B);

   /*
     Compute
         p*x^(p-1)*(Frob(y))^-1  (or 3 if cube, which isn't true here!)
   */

   // print "Computing (y^Frobenius)^(-1)";
   tyme :=Cputime();
   x_pow := x^(p-1);
   poly := PolynomialRing(R1)![a6,a4,0,1];
   difl := FrobYInv(poly, p, N, x_pow, S,cube) * x_pow;
   x_pow := x_pow*x;
   // printf "Expansion time: %o\n",Cputime(tyme);

   /* Calculate the rows of the Frobenius matrix */
   // print "Reducing differentials modulo cohomology relations.";
   R1 := pAdicField(p,N1);
   M := MatrixAlgebra(R1,d)!0;
   i := 1;
   boundu := p*N + (p div 2) - 1;
   S1 := PolynomialRing(BaseRing(S));
   while true do
       boundl := (p div 2) - Floor((i*p-1)/(d+1));
       polys,bot := Convert(S1!difl, boundu, boundl, K);
       diffr := ReduceA([polys[k] : k in [1..Min(1-bot,#polys)]], precs, -bot)+
                ReduceB([polys[k] : k in [Max(2-bot,1)..#polys]], precs, Max(bot,1), cube);
       M[i] := RSpace(R1,d)![R1!Coefficient(diffr,k) : k in [0..(d-1)]];
       if i eq d then
          break;
       end if;
       i +:= 1;
       difl *:= x_pow;
   end while;

   /*
      We know Frob=M to precision N1.  Compute the N1-th power of
      Frob.  The second *row* (since MAGMA matrices act from left)
      contains linear combination of omega and eta that equals
      Frob^N1(eta).
   */
   A := M^N1;
   return -12 * (A[2,1]/A[2,2]), M;

end intrinsic;


/*****************************************************************

    Computing sigma: Cohomology algorithm

 ****************************************************************/

intrinsic sigma_using_e2(E::CrvEll, p::RngIntElt, prec::RngIntElt :
                    e2prec := 50) -> .
{
   Compute and return the formal series for sigma(t) using
   the Katz/Kedlaya algorithm for computing E2.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> sigma_using_e2(E, 5, 10);
t + O(5^52)*t^2 + (331883312126563673413235306321062928 +
    O(5^52))*t^3 - (1110223024625156540423631668090820312 +
    O(5^52))*t^4 - (114496958818289158695517187452183148 +
    O(5^51))*t^5 + (445053157775380010065972176906892 +
    O(5^49))*t^6 - (246425046817409275170929836993681 +
    O(5^48))*t^7 - (248945902282571925665450930262558 +
    O(5^47))*t^8 - (9730291437202192499870687559126*5 +
    O(5^47))*t^9 + O(t^10)
}
   prec +:= 2;

   // 1. Compute value of p-adic E_2, using Kedlaya's MW algorithm
   e2 := E2(E,p,e2prec);

   // 2. Define some notation.
   K := pAdicField(p,Precision(Parent(e2)));
   a1, a2, a3, a4, a6 := Explode(aInvariants(E));
   S<t> := LaurentSeriesRing(K);

   // 3. Compute formal Weierstrass p function in terms of t.
   x  := S!formal_x(E, prec);
   wp := x + (a1^2 + 4*a2)/12;

   // 4. Compute wp(z), where z = FormalLog(t)
   z_of_t := S!formal_log(E, prec);
   t_of_z := Reverse(z_of_t);
   wp_of_z := Evaluate(wp, t_of_z);
   R<z> := S;

   // 5. Compute
   pr := AbsolutePrecision(wp_of_z);

   // Let g be    1/z^2 - \wp + E2/12.  Notice that we compute this
   // by just ignoring the coefficients of wp_of_z for n=-2.
   g := e2/12 +  &+[z^n*Coefficient(-wp_of_z,n) : n in [-1..pr-1]] + O(z^pr);
   /* It's *really* weird that it isn't -E2/12, since the canonical
      eks function is \wp + E2/12.

      There's a factor of -1 in my definition of E2 in the other file,
      which I think Tate and I got out of Katz's paper.   Maybe that
      sign is wrong?   There's some -1 normalization that is different
      in two papers. */

   // 6. Integrate twice, exponentiate, etc.
   sigma_of_z := z*my_exp(my_integral(my_integral(g)));
   sigma_of_t := Evaluate(sigma_of_z, z_of_t);
   R<t> := PowerSeriesRing(K);
   sigma := R!sigma_of_t + O(t^(prec-2));
   return sigma;
end intrinsic;


/*****************************************************************

    Computing the p-adic Height Function

 ****************************************************************/

intrinsic prep_multiple(P::PtEll, p::RngIntElt) -> PtEll, RngIntElt
{
Returns an integer n such that n*P is in the subgroup of E that
reduces to the identity mod p and lands in the connected component of
the Neron model modulo all primes.

EXAMPLE:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> P := E![0,0];
> prep_multiple(P,5);
8
> prep_multiple(P,11);
17
}
   require IsPrime(p) : "Argument 2 must be prime.";

   function prime_divisors(n)
      return [F[1] : F in Factorization(n)];
   end function;
   E := Curve(P);
   c := LCM([TamagawaNumber(E,q) : q in prime_divisors(Conductor(E))]);
   N := #ChangeRing(E,GF(p));
   n := LCM(c,N);
   return n;
end intrinsic;

intrinsic log(x::FldPadElt) -> FldPadElt
{
Compute and return the value of the unique extension of log to a
group homomorphism on all Q_p.

EXAMPLES:
> K := pAdicField(5);
> log(K!(1+5+2*5^2));
771685135946*5 + O(5^20)
> log(K!(5+2*5^2));
-159535800608*5 + O(5^20)
> log(K!(1/5+2*5^2));
-69652953373*5^3 + O(5^20)
}
   p := Prime(Parent(x));
   u := x * p^(-Valuation(x));
   return Log(u^(p-1))/(p-1);
end intrinsic;

intrinsic x_coord_d(Q::PtEll) -> RngIntElt
{
Computes and returns the positive squareroot of the denominator
of the x coordinate of Q.  Bombs if the denominator is not a perfect square.
}
   t, d := IsSquare(Denominator(Q[1]));
   assert t;
   return d;
end intrinsic;

intrinsic my_exp(f::.) -> RngSerLaurElt
{
Computes exp(f) to the precision of f.  This is here because the Exp function
built into MAGMA does not work as it should.
}
   return &+[f^n/Factorial(n) : n in [0..AbsolutePrecision(f)]];
end intrinsic;

intrinsic my_is_zero(a::.) -> .
{
Determines whether or not the p-adic number a is 0.
This is here because testing for equality with 0
for p-adic numbers in MAGMA does not work as it should.
}
   return Valuation(a) ge Precision(Parent(a));
end intrinsic;

intrinsic my_integral(f::RngSerLaurElt) -> RngSerLaurElt
{
 Integrate Laurent series, without stupid check on coefficient
 of t^(-1) being 0, which doesn't work in MAGMA, since nonzero
 for p-adic elements is messed up.
}
   require my_is_zero(Coefficient(f,-1)) : "Coefficient of 1/t must be zero.";
   pr := AbsolutePrecision(f);
   t := Parent(f).1;
   return &+[(Coefficient(f,n)/(n+1))*t^(n+1) : n in [0..pr-1] |
                not my_is_zero(Coefficient(f,n)) ] + O(t^pr);
end intrinsic;

intrinsic height_function(E::CrvEll, p::RngIntElt, prec::RngIntElt :
            use_integrality := false) -> .
{
Returns the CANONICAL p-adic global height function.
Note that we normalize this height so it takes values in Z_p, unless
p is anomalous in which case it takes values in (1/p)*Z_p.

INPUT:
  E -- an elliptic curve over Q,
  p -- a prime such that E has good ordinary reduction at p,
  prec -- precision parameter (see below).

OPTIONS:
  use_integrality -- this defaults to false.

             * false: Use the E2 algorithm for computing sigma.
                      The heights are computed to precision very
                      close to O(p^prec), but maybe slightly less.
                      I haven't figured out why the precision is
                      slightly less then O(p^prec), since I think
                      it should be exactly O(p^prec), but see that
                      is is by computing the height with larger
                      and larger precision and consider the valuations
                      of the differences.

             * true: Use the integrality of sigma algorithm to
                     compute sigma.  The prec argument then determines
                     the number of terms of sigma used in computing s_2.
                     For example, maybe around 40 or 50 terms gives s_2 to
                     precision O(p^3), at least when p=5 and E is 37A.

EXAMPLES:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> P := E![0,0];
> h := height_function(E, 5, 10);
> h(P);
10120297 + O(5^11)
> h := height_function(E, 5, 20);
> h(P);
-201147061754703 + O(5^21)
> Valuation(-201147061754703 - 10120297, 5);
8
> h := height_function(E, 5, 30);
> h(P);
-1305176752965909410953 + O(5^31)
> Valuation(-1305176752965909410953  + 201147061754703, 5);
18
> h := height_function(E, 5, 60 : use_integrality := true);  // slow
> h(P);
> -3 + O(5^2)                                                // pitiful precision
}
   require Conductor(E) mod p ne 0 and ap(E,p) mod p ne 0 : "Curve must have good ordinary reduction at p.";

   if use_integrality then
      sig := sigma_using_integrality(E, p, prec) ;
   else
      sig := sigma_using_e2(E, p, prec : e2prec := prec) ;
   end if;
   function h(P)
      assert Curve(P) eq E;
      n := prep_multiple(P, p);
      Q := n*P;
      d := x_coord_d(Q);  // positive sqrt of denominator
      t := t_val(Q);
      v := Evaluate(sig,t);
      HQ := log(v/d);
      HP := HQ/n^2;
      return HP/p;           // note the normalization
   end function;
   return h;
end intrinsic;

intrinsic regulator(E::CrvEll, p::RngIntElt, prec::RngIntElt) -> .
{
   Computes and returns the p-adic regulator of E, which is the
   discriminant of the height pairing on the Mordell-Weil group E(Q).

INPUT:
   E -- an elliptic curve over Q.
   p -- a good ordinary prime for E.
   prec -- the precision that is input to the height_function command, which
           is computed using the E2 algorithm.

EXAMPLES:
> E := EllipticCurve([ 0, 0, 1, -1, 0 ]);
> regulator(E,5,20);
138360994885922 + O(5^21)
> regulator(E,7,20);
4709403600911866 + O(7^20)
> E := EllipticCurve(CremonaDatabase(), "389A");
> regulator(E,7,20);
22975764581280320 + O(7^20)
}
   require Conductor(E) mod p ne 0 and ap(E,p) mod p ne 0 :
                    "Curve must have good ordinary reduction at p.";

   h := height_function(E, p, prec);
   G, f := MordellWeilGroup(E);
   B := [f(G.i) : i in [1..Ngens(G)] | Order(G.i) eq 0];
   function pairing(P,Q)
      return (h(P+Q)-h(P)-h(Q))/2;
   end function;
   K := Parent(h(B[1]));
   M := MatrixAlgebra(K,#B)!0;
   for i in [1..#B] do
      for j in [i..#B] do
         aij := pairing(B[i],B[j]);
         M[i,j] := aij;
         M[j,i] := aij;
      end for;
   end for;
   return Determinant(M);
end intrinsic;


/***********************************************************************

   Miscellaneous functions that are useful to have around.

 ***********************************************************************/
intrinsic ap(E::CrvEll, p::RngIntElt) ->RngIntElt
{
  Returns #E(F_p) - (p + 1).
}
   return TraceOfFrobenius(ChangeRing(E,GF(p)));
end intrinsic;

intrinsic good_ordinary_primes(E::CrvEll, pmax::RngIntElt) -> SeqEnum
{
   Compute and return the primes p>=5 that are good ordinary for E and are
   less than pmax.
}
   N := Conductor(E);
   return [p : p in [5..pmax] | IsPrime(p) and
                                N mod p ne 0 and
                                ap(E,p) mod p ne 0];
end intrinsic;

intrinsic print_padic(x::.)  -> .
{
I can't figure out how to print p-adics correctly in MAGMA (as power
series in p) anymore.   This function isn't very refined though, since
e.g., it outputs 3*5^0 + 2*5^2 instead of 3+2*5^2.

EXAMPLE:
> K := pAdicField(5);
> a := K!9484945; a;
1896989*5 + O(5^21)
> print_padic(a);
4*5^1 + 2*5^2 + 4*5^3 + 2*5^6 + 5^7 + 4*5^8 + 4*5^9 + O(5^10)
}
   z := Integers()!x;
   p := Prime(Parent(x));
   i := 0;
   s := "";
   while z ne 0 and i lt Precision(Parent(x)) do
      c := z mod p;
      if c ne 0 then
         if c ne 1 then
            s *:= Sprintf("%o*", c);
         end if;
         s *:= Sprintf("%o^%o + ", p, i);
      end if;
      z := z div p;
      i := i + 1;
   end while;
   s *:= Sprintf("O(%o^%o)", p, i);
   return s;
end intrinsic;


