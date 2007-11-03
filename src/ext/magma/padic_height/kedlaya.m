freeze;
/*******************************************************************
 *******************************************************************

      KEDLAYA'S ALGORITHM (in ODD CHARACTERISTIC)

              Mike Harrison.

*******************************************************************
*******************************************************************/

LINEAR_LIMIT := 7;
ALWAYS_CUBE := false;

function GetRing(R1,Q,prec)

    L<T> := LaurentSeriesRing(ChangePrecision(R1,prec));
    P1 := PolynomialRing(L);
    return quo<P1|P1!Q-T^-1>;

end function;

function myXGCD(p1,p2)

    // p1 and p2 are relatively prime integral polynomials over an
    // unramified pAdic field of bounded precision n.
    // It is assumed that the leading coeff of p1 is prime to p.
    // The function calculates and returns integral polynomials of
    // the same precision s.t. A*p1+B*p2 = 1 mod p^n.
    // Begins by lifting the Xgcd result over the Residue field and
    // iterates up mod p^i for 1 <= i <= n.
    R := Parent(p1);
    F,mp := ResidueClassField(IntegerRing(BaseRing(R)));
    S := PolynomialRing(F);
    p1r := S![mp(j) : j in Coefficients(p1)];
    p2r := S![mp(j) : j in Coefficients(p2)];
    _,Ar,Br := Xgcd(p1r,p2r);
    u := R!Ar; v := R!Br;
    A := u; B := v;
    delta := R!-1;
    p := Characteristic(F);
    for i in [1..BaseRing(R)`DefaultPrecision-1] do
        delta := (delta+u*p1+v*p2)/p;
        delr := S![mp(j) : j in Coefficients(delta)];
        v1 := (-Br*delr) mod p1r;
        v := R!v1;
        u := R!((-(delr+v1*p2r)) div p1r);
        /* NB. the simpler
          u := R!((-Ar*delr) mod p2r);
          v := R!((-Br*delr) mod p1r);
           doesn't necessarily work if p2 has leading
           coeff div by p, when deg p2r < deg p2.
            In this case, u*p1+v*p2 != -delta mod p if
           deg p1r+deg p2r <= deg delr (< deg p1+deg p2)
	*/
        A +:= u*(p^i);
        B +:= v*(p^i);
    end for;
    return A,B;

end function;

function SigmaPowMat(M,m,n)

    // returns s^(n-1)(M)*s^(n-2)(M)*..*s(M)*M where M is
    // an m*m matrix over an unramified pAdic field and
    // s is the absolute Frobenius of that field. n >= 1.
    // Uses a basic left-right binary powering algorithm.

    if n eq 1 then return M; end if; //special case
    bits := Intseq(n,2);
    N := M;
    r := 1;
    for i := (#bits-1) to 1 by -1 do
        N := Matrix(m,[GaloisImage(x,r) : x in Eltseq(N)])*N;
        r *:= 2;
	if bits[i] eq 1 then
	    N := Matrix(m,[GaloisImage(x,1) : x in Eltseq(N)])*M;
            r +:= 1;
	end if;
    end for;
    return N;

end function;

function FrobYInv(Q,p,N,x_pow,S,cube)

    // Computes p*(Frob y)^-1 (cube=false) or p*(Frob y)^(-3) (cube=true)
    // mod p^N.
    // This means calculating
    //     (1+((FrobQ)(x^p)-(Q(x))^p)*X^p)^(-1/2){resp(-3/2)}
    // to the required precision in S (S defined in Kedlaya fn).
    //
    // Starts by computing terms in the binomial expansion of the above
    // then uses Newton iteration. The first part computes the most
    // appropriate number <= LINEAR_LIMIT of terms for the Newton phase.
    // In the Newton phase, the power series (in T) coefficients of powers
    // of x are truncated noting that result mod p^t contains no non-zero
    // coefficients of T beyond T^(p*(t-1)).
    //
    // Q(x) is the defining polynomial for the hyperelliptic curve (y^2=Q(x))
    // and x_pow = x^(p-1) in S.

    R1 := BaseRing(BaseRing(S));
    T := BaseRing(S).1;
    E := 1 + (T^p)*(Evaluate(Parent(Q1)![GaloisImage(j,1):
            j in Coefficients(Q1)], x_pow*S.1) - Evaluate(Q1,S.1)^p)
                where Q1 is PolynomialRing(R1)!Q;
    // now compute E^(-1/2) (E^(-3/2) if cube) mod p^N
    // ( this is weaker than mod T^(p(N-1)+1) )
    // by "linear lift" followed by Newton iteration.

    //first backtrack to find starting point for Newton phase
    prec := N;
    adjs := [Booleans()|];
    while prec gt LINEAR_LIMIT do
        Append (~adjs,IsOdd(prec));
        if p eq 3 then
	    prec := prec div 2;
	else
            prec := (prec+1) div 2;
        end if;
    end while;
    // now do the linear part
    Sc := GetRing(R1,Q,prec);
    Rc := BaseRing(BaseRing(Sc));
    u := Sc!1;
    Epow := u;
    E1 := Sc![BaseRing(Sc)!a : a in Coefficients(E)];
    half := cube select -(1/2) else (1/2);
    bincoeff := 1/1;
    for i in [1..prec] do
        bincoeff *:= (half-i)/i;
        Epow := (E1-1)*Epow;
	u := u + (Rc!bincoeff)*Epow;
    end for;
    delete Epow;
    // u is now the answer mod p^prec. Finish by Newton iteration
    // x -> (3*x - E*x^3)/2.
    if cube then E := E^3; end if;
    half := BaseRing(Parent(Q))!(1/2);
    if prec eq N then PolR := PolynomialRing(BaseRing(S)); end if;
    while prec lt N do
        tyme := Cputime();
	prec *:= 2;
        if adjs[#adjs] then
            prec +:= ((p eq 3) select 1 else -1);
        end if;
        Prune(~adjs);
        if prec lt N then
            Sc := GetRing(R1,Q,prec);
            E1 := Sc![BaseRing(Sc)!a : a in Coefficients(E)];
        else
            Sc := S; E1 := E;
        end if;
        T := BaseRing(Sc).1;
        PolR := PolynomialRing(BaseRing(Sc));
        u := Sc![BaseRing(Sc)!a : a in Coefficients(u)];
        v := Sc![j+O(T^(p*prec-(p-1))) : j in Coefficients(PolR!(u^2))];
        u := (3*u - E1*(u*v))*(BaseRing(BaseRing(Sc))!(1/2));
        // remove O(T^r) terms
        u := Sc![ elt<Parent(j)|v,coeffs,-1> where coeffs,v is
                  Coefficients(j) : j in Coefficients(PolR!u)];
        vprintf JacHypCnt, 3:
           "Reached precision %o   time = %o\n",prec,Cputime(tyme);
    end while;
    // now "clean" the result mod T^(pN-(p-1))
    u := S![ elt<Parent(j)|v,coeffs,-1> where coeffs,v is
       Coefficients(j+O(T^((p*(N-1))+1))) : j in Coefficients(PolR!(p*u))];
    return (u * T^(((cube select 3 else 1)*p) div 2));

end function;

function Convert(powTx,bdu,bdl,K)
    // Takes a differential powTx*(dx/y) where powTx is of the form
    //   p0(T) + p1(T)*x +... pr(T)*x^r
    //    (T := 1/y^2, pi(T) are finite-tailed Laurent series)
    // and returns [Ar,A(r+1),...],r where Ai = Coefficients(ai)
    //  powTx = ar(x)*T^r + a(r+1)(x)*T^(r+1) + ... (ai in K[x]).
    // K is the unramified pAdic coefficient field.
    // bdu,bdl are upper and lower bounds for non-zero powers of
    // T in {p0,p1,...}.

    vec := [PowerSequence(K)|];
    found_start := false;
    start_adj := 0;
    for i in [bdl..bdu] do
        v1 := [K!Coefficient(ser,i) : ser in Coefficients(powTx)];
	if not found_start then
	    if &and[RelativePrecision(k) eq 0 : k in v1] then
	        start_adj +:= 1;
		continue;
	    else
	        found_start := true;
	    end if;
	end if;
	Append(~vec,v1);
    end for;
    return vec,bdl+start_adj;

end function;

function PrecDiv(pol,d)

    // Used by ReduceB to avoid losing padic precision when
    // dividing polynomial pol by integer d (p may divide d).
    // If d isn't a padic unit, arbitrary "noise" is added to
    // restore full precision after the division.

    K := BaseRing(Parent(pol));
    pold := pol/d;
    if Valuation(d,Prime(K)) ne 0 then
        pold := Parent(pol)!
           [(RelativePrecision(x) eq 0) select
           O(UniformizingElement(K)^K`DefaultPrecision) else
           ChangePrecision(x,Max(0,K`DefaultPrecision-Valuation(x))) :
            x in Coefficients(pold)];
    end if;
    return pold;

end function;

function ReduceA(polys,precs,n0)

    // Performs the reduction of
    // {a_(n0-1)(x)*y^(2*(n0-1)) + .. + a1(x)*y^2 + a0(x)}*(dx/y)
    // to canonical return_val *(dx/y) form.

    PolR := Parent(precs[1]);
    if IsEmpty(polys) then
        return PolR!0;
    end if;
    d := Degree(precs[1]);
    pol := PolR!polys[1];
    N := BaseRing(PolR)`DefaultPrecision;
    for k := (n0-1) to 0 by -1 do
        pol := ((k le (n0-#polys)) select PolR!0 else PolR!(polys[n0+1-k])) +
	              pol*precs[1];
	for j := (2*d-1) to d by -1 do
	    lc := Coefficient(pol,j);
	    if RelativePrecision(lc) ne 0 then
	        pol := pol - ((PolR.1)^(j-d))*
                (ChangePrecision(lc,N)/((2*k-1)*d+2*(j+1))) *
                ((2*k+1)*precs[2]*PolR.1+2*(j+1-d)*precs[1]);
	    end if;
	end for;
        pol := PolR![Coefficient(pol,i):i in [0..(d-1)]];
    end for;
    lc := Coefficient(pol,d-1);
    if RelativePrecision(lc) ne 0 then
        pol := pol - (ChangePrecision(lc,N)/d)*precs[2];
    end if;
    return pol;

end function;

function ReduceB(polys,precs,n0,cube)

    // Performs the reduction of
    // {a1(x)*(1/y^2) + .. a_n0(x)*(1/y^2n0)}*(dx/y^r)
    // to canonical return_val *(dx/y^r) form.
    // r = 1 if cube = false, else r = 3.

    PolR := Parent(precs[1]);
    if IsEmpty(polys) then
        return PolR!0;
    end if;
    divcase := (Valuation(LeadingCoefficient(precs[2])) gt 0);
    pol := PolR!polys[#polys];
    for k := (n0-1+#polys) to (cube select 2 else 1) by -1 do
        p1 := (pol*precs[4]) mod precs[1];
        if divcase then
           p2 := (pol-p1*precs[2]) div precs[1];
        else
           p2 := (pol*precs[3]) mod precs[2];
        end if;
        pol := ((k le n0) select PolR!0 else PolR!(polys[k-n0])) +
                   p2 + PrecDiv(2*Derivative(p1),2*k-1);
    end for;
    return pol;

end function;

function UpperCoeffs(M,g,ppow,e_val)

    // Calcs the sequence of the upper coefficients (x^g,..x^(2g))
    // of CharPoly(M) using the trace method. The magma intrinsic
    // could be used but this is slightly more efficient as only
    // Trace(M),Trace(M^2),..Trace(M^g) need be calculated rather
    // than Trace(M),..,Trace(M^(2g)).
    // If Nrows(M) = 2*g+1 then the extra eigenvalue of M, e_val,
    // is discarded. (e_val = q (1) if cube = false (true)).
    // The sequence [r(v)] is returned where, for a p-adic int v,
    // r(v) is the integer n s.t.|n|<ppow/2 and n=v mod ppow.

    pAd := pAdicField(Parent(M[1,1]));
    N := M;
    UCs := [pAd!1];   // coeffs (highest to lowest) of CharPoly(M)
    TrPows := [pAd|]; // [Trace(M),Trace(M^2),...]
    for i in [1..g] do
        Append(~TrPows,Eltseq(Trace(N))[1]);
        Append(~UCs, (- &+[TrPows[j]*UCs[i+1-j] : j in [1..i]])/i);
        if i lt g then N := N*M; end if;
    end for;
    if Nrows(M) ne 2*g then // original Q(x) of even degree
        for i in [1..g] do
            UCs[i+1] := UCs[i+1]+e_val*UCs[i];
        end for;
    end if;
    return [((2*uc gt ppow) select uc-ppow else uc) where uc is
              (IntegerRing()!x) mod ppow : x in UCs];

end function;

/*
    if p eq 0 or n eq 0 then
       R := Parent(poly);
       k := BaseRing(R);
       p := Characteristic(k);
       n := Degree(k);
    end if;
*/

function Kedlaya(poly, p, n)
    // Main function for the (odd char) Kedlaya algorithm.
    //  Computes the numerator of the zeta function for
    //       y^2 = poly (defined over over GF(p^n)).
    //  The cube boolean variable indicates which differential
    // basis is used for cohomological calculations -
    //   (dx/y), x(dx/y), x^2(dx/y), ... if cube = false
    //   (dx/y^3), x(dx/y^3), ...        if cube = true
    //  The 1st is the natural basis, but leads to a non-integral
    // transformation matrix if p is small compared to the genus.
    // The strategy is to use the first if this is not the case
    // unless the ALWAYS_CUBE value is true

    d := Degree(poly)-1;
    g := d div 2;
    cube := true;
    if not ALWAYS_CUBE then
        if (IsEven(d) and p ge d) or  // deg=2*g+1
           (IsOdd(d) and p gt g) then // deg=2*g+2
            cube := false;
        end if;
    end if;

    N1 := Ceiling((n*g/2)+Log(p,2*Binomial(2*g,g)));
    N := N1 + Floor(Log(p,2*N1))+1;
    K := UnramifiedExtension(pAdicField(p,N),n);
    R1 := quo<Integers(K)|UniformizingElement(K)^N>;
    Embed(BaseRing(Parent(poly)),ResidueClassField(R1));
    S := GetRing(R1,poly,N);
    x := S.1;
    //R<T> := LaurentSeriesRing(R1);
    //S1<x> := PolynomialRing(R);
    precs := [PolynomialRing(K)!poly];
    //S<x> := quo<S1|S1!poly-T^-1>;
    Append(~precs,Derivative(precs[1]));
    A,B := myXGCD(precs[1],precs[2]);
    // A,B satisfy A*Q+B*Q'=1 where Q is the lift of poly to char 0.
    Append(~precs,A);
    Append(~precs,B);

    //Stage 1 - compute p*x^(p-1)*(y^Frob)^-1[3]
    vprintf JacHypCnt, 2:
      "Computing (y^Frobenius)^-%o ...\n",cube select 3 else 1;
    tyme :=Cputime();
    x_pow := x^(p-1);
    difl := FrobYInv(PolynomialRing(R1)!poly,p,N,x_pow,S,cube)*x_pow;
    x_pow := x_pow*x;
    vprintf JacHypCnt, 2: "Expansion time: %o\n",Cputime(tyme);

    //Stage 2 - loop to calculate the rows of the Frobenius matrix
    vprint JacHypCnt, 2: "Reducing differentials ...";
    R1 := quo<Integers(K)|UniformizingElement(K)^N1>;
    M := RMatrixSpace(R1,d,d)!0;
    i := 1;
    boundu := p*N+(p div 2)-1;
    S1 := PolynomialRing(BaseRing(S));
    vtime JacHypCnt, 2:
    while true do
        boundl := (p div 2) - Floor((i*p-1)/(d+1));
        polys,bot := Convert(S1!difl,boundu,boundl,K);
        diffr := ReduceA([polys[k] : k in [1..Min(1-bot,#polys)]],precs,-bot)+
	  ReduceB([polys[k] : k in [Max(2-bot,1)..#polys]],precs,Max(bot,1),cube);
	M[i] := RSpace(R1,d)![R1!Coefficient(diffr,k) : k in [0..(d-1)]];
	if i eq d then break; end if;
	i +:= 1;
	difl := difl * x_pow;
    end while;
    print "M = ", M;
    print "CharacteristicPolynomial of M = ", CharacteristicPolynomial(M);

    vprint JacHypCnt, 2: "Computing Frobenius matrix by powering ...";
    vtime JacHypCnt, 2:
    M := SigmaPowMat(M,d,n);
    // Now change the precision to N1+Val(p,g!). The Val(p.. is needed
    // to add random noise as the characteristic polynomial calculation
    // uses the Trace method which involves dividing by 1,2,..g for the
    // top terms (later terms aren't needed) with a corresponding loss
    // in precision.
    N2 := N1 + Valuation(Factorial(g),p);
    if N2 gt N then ChangePrecision(~K,N2); end if;
    M := Matrix(K,M);
    if N2 gt N1 then
      M := Matrix(K,d,[ChangePrecision(K!x,N2-Valuation(x)) : x in Eltseq(M)]);
    end if;
    tyme := Cputime();
    q := #ResidueClassField(Integers(K));
    UCoeffs := UpperCoeffs(M,g,p^N1,cube select 1 else q);
    CharP := PolynomialRing(IntegerRing())!
                  ([UCoeffs[i]*q^(g+1-i): i in [1..g]] cat Reverse(UCoeffs));
    vprintf JacHypCnt,3: "Characteristic polynomial time: %o\n",Cputime(tyme);
    return Parent(CharP)!Reverse(Coefficients(CharP));

end function;

intrinsic kedlaya(poly::., p::., n::.) -> .
{}
   return Kedlaya(poly, p, n);
end intrinsic;


function KedCharPolyOdd(C)

   Fld := BaseField(C);
   if Type(Fld) ne FldFin then
       error "The curve must be defined over a finite field!";
   end if;
   p := Characteristic(Fld);
   if p eq 2 then
       error "Sorry! The kedlaya method can't currently handle char 2";
   end if;
   n := Degree(Fld);
   C1 := SimplifiedModel(C);
   Q := HyperellipticPolynomials(C1);
   twist := false;
   lc := LeadingCoefficient(Q);
   if lc ne Fld!1 then
       if IsOdd(Degree(Q)) then
           Q := Evaluate(Q,lc*Parent(Q).1);
	   Q := Q/(lc^(Degree(Q)+1));
       else
           Q := Q/lc;
           if not IsSquare(lc) then twist := true; end if;
       end if;
   end if;
   vprint JacHypCnt, 1: "Applying Kedlaya's algorithm";
   tyme := Cputime();
   cpol := Kedlaya(Q,p,n);
   vprintf JacHypCnt, 1: "Total time: %o\n", Cputime(tyme);
   if twist then
        cpol := Evaluate(cpol,-(Parent(cpol).1));
   end if;
   return cpol;

end function;

intrinsic kedlayacharpoly(C::.) -> .
{}
  return  KedCharPolyOdd(C);
end intrinsic;


/*
Date: Thu, 8 Jul 2004 13:53:13 -0400
From: Nick Katz <nmk@Math.Princeton.EDU>
Subject: Re: convergence of the Eisenstein series of weight two
To: mazur@math.harvard.edu, nmkatz@Math.Princeton.EDU
Cc: tate@math.utexas.edu, was@math.harvard.edu

It seems to me you want to use the interpretation of P as the
"direction of the unit root subspace", that should make it fast to
compute. Concretely, suppose we have a pair (E, \omega) over Z_p, and
to fix ideas p is not 2 or 3.  Then we write a Weierstrass eqn for E,
y^2 = 4x^3 - g_2x - g_3, so that \omega is dx/y, and we denote by \eta
the differential xdx/y. Then \omega and \eta form a Z_p basis of
H^1_DR = H^1_cris, and the key step is to compute the matrix of
absolute Frobenius (here Z_p linear, the advantage of working over
Z_p: otherwise if over Witt vectors of an F_q, only \sigma-linear).
[This calculation goes fast, because the matrix of Frobenius lives
over the entire p-adic moduli space, and we are back in the glory days
of Washnitzer- Monsky cohomology (of the open curve E - {origin}.]

        Okay, now suppose we have computed the matrix of Frob in the
basis \omega, \eta. The unit root subspace is a direct factor, call it
U, of the H^1, and we know that a complimentary direct factor is Fil^1
:= the Z_p span of \omega. We also know that Frob(\omega) lies in
pH^1, and this tells us that, mod p^n, U is the span of
Frob^n(\eta). What this means concretely is that if we write, for each
n,

   Frob^n(\eta) = a_n\omega + b_n\eta,

then b_n is a unit (cong mod p to the n'th power of the hasse
invariant) and that P is something like the ratio a_n/b_n (up to a
sign and a factor 12 which i don't recall offhand but which is in my
Antwerp appendix and also in my "p-adic interp. of real
anal. Eis. series" paper.

        So in terms of speed of convergence, ONCE you have Frob, you
have to iterate it n times to calculate P mod p^n. Best, Nick
*/

intrinsic padicprint(x::.) ->.
{}
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




/*
*/

intrinsic TateCM_Examples()
{}
   c4 := 2^9*11;
   c6 := 2^9*7*11^2;
   a4 := - c4/(2^4*3);
   a6 := - c6/(2^5*3^3);
   E := EllipticCurve([a4,a6]);
   print "sqrt(-11): mod5    ", (Integers()!E2(E,5,30)) mod 5;
   print "sqrt(-11): mod5^5  ", (Integers()!E2(E,5,30)) mod 5^5;

   c4 := 2^4*3*11;
   c6 := 2^6*3^3*11;
   a4 := - c4/(2^4*3);
   a6 := - c6/(2^5*3^3);
   E := EllipticCurve([a4,a6]);
   print "j=1: mod5    ", (Integers()!E2(E,5,30)) mod 5;
   print "j=1: mod5^5  ", (Integers()!E2(E,5,30)) mod 5^5;

   c4 := 2^5*11;
   c6 := -2^3*7*11^2;
   a4 := - c4/(2^4*3);
   a6 := - c6/(2^5*3^3);
   E := EllipticCurve([a4,a6]);
   print "j=2: mod5    ", (Integers()!E2(E,5,30)) mod 5;
   print "j=2: mod5^5  ", (Integers()!E2(E,5,30)) mod 5^5;

   c4 := -2^4*3;
   c6 := 0;
   a4 := - c4/(2^4*3);
   a6 := - c6/(2^5*3^3);
   E := EllipticCurve([a4,a6]);
   print "j=3: mod5    ", (Integers()!E2(E,5,30)) mod 5;
   print "j=3: mod5^5  ", (Integers()!E2(E,5,30)) mod 5^5;

   c4 := 2^5*3*19;
   c6 := -2^3*3^3*19^2;
   a4 := - c4/(2^4*3);
   a6 := - c6/(2^5*3^3);
   E := EllipticCurve([a4,a6]);
   print "j=4: mod5    ", (Integers()!E2(E,5,30)) mod 5;
   print "j=4: mod5^5  ", (Integers()!E2(E,5,30)) mod 5^5;

end intrinsic;

intrinsic Lagrange(points::SeqEnum) -> RngUPolElt
{
Input:
   points -- a sequence of length n>=1 of pairs (xi,yi=f(xi)), with the xi distinct

Output:
   a polynomial f of degree <= n-1 that passes through those points.
}
   n := #points;
   R<x> := PolynomialRing(FieldOfFractions(Parent(points[1][1])));
   function P(j)
      return points[j][2]* &*[(x-points[k][1])/(points[j][1]-points[k][1]) : k in [1..n] | k ne j];
   end function;
   return &+[P(j) : j in [1..n]];

end intrinsic;

intrinsic IsGoodOrdinary(E::CrvEll, p::RngIntElt) -> BoolElt
{}
   return Conductor(E) mod p ne 0 and
           TraceOfFrobenius(ChangeRing(E,GF(p))) ne 0;
end intrinsic;

/*
Family:

 y^2 = x^3 + x*t + t;

over Q(t).
*/

intrinsic E2oft(t::RngIntElt, p::RngIntElt, prec::RngIntElt) -> FldPadElt
{}
   E := EllipticCurve([t,t]);
   require IsGoodOrdinary(E,p) : "Fiber curve E_t must be good ordinary at argument 2.";
   return E2(E, p, prec);
end intrinsic;

intrinsic E2_def(p::RngIntElt, radius::RngIntElt, num_points::RngIntElt) -> .
{}
   points := [];
   n := 1;
   while #points lt num_points do
      //t := p*Random(2,20) + p^radius*n;
      t := p^radius*n;
      n +:= Random(1,20);
      E := EllipticCurve([-1296+t,11664+t]);
      if IsGoodOrdinary(E,p) then
         e := E2(E,p,20);
         v := [t,e] cat [Integers()!e mod p^a : a in [1..radius+5]];
         Append(~points, v);
      end if;
   end while;
   f := Lagrange(points);
   return f, points, [Valuation(c) : c in Eltseq(f)];
end intrinsic;
