freeze;

// William Stein -- I changed the ClassGroup call on line 191,
// so when bound 0 is passed in the result will be provably
// correct, at least **according to the ClassGroup section of the
// MAGMA documentation**.

// ALWAYS DOCUMENT YOUR CHANGES (what/when/how)
// in the header of this file
// these package files are updated concurrently by
// ourselves (magma) AND M. Stoll
// Thank you!


/*******************************************
 * Hyperelliptic/3descent.m                *
 *                                         *
 * 3-descent for elliptic curves           *
 *                                         *
 * Michael Stoll                           *
 *                                         *
 * started 07-Aug-1999                     *
 *                                         *
 * Needs selmer.m                          *
 *******************************************/

 /*
13/7/01	nicole	comments as to obsoleteness of first two functions




   2001-07: Paulette
   Scheme merge: new types now (CrvHyp, SetPtHyp & PtHyp)
   (there was nothing)

  ------------------------------------------------------------------*/

// Nicely formatted output of sequences

procedure PrintBasis(B, indent)
  for b in B do
    printf " "^indent*"%o\n", b;
  end for;
end procedure;

function Denom(x)
  // x::FldNumElt
  den1 := LCM([Denominator(c) : c in Eltseq(x)]);
  O := Integers(Parent(x));
  x1 := O!(den1*x);
  num := GCD(den1, GCD(ChangeUniverse(Eltseq(x1), Integers())));
  return ExactQuotient(den1, num), O![ExactQuotient(c, num) : c in Eltseq(x1)];
end function;

function MyResidueClassField(pid)
  // Given a prime ideal, returns the residue field Order(pid)/pid
  // together with the reduction map from the ring of pid-integral
  // elements in the corresponding number field to the residue field.
  // This map has an inverse.

  // Type checking on two possible argument types:

  if Type(pid) eq RngOrdIdl then
    // return ResideuClassField(Order(pid), pid);
    O := Order(pid);
    K := FieldOfFractions(O);
    F, m := ResidueClassField(O, pid);
    p := Minimum(pid);
    e := ChineseRemainderTheorem(pid, ideal<O | O!p>/pid^RamificationIndex(pid),
                                 O!1, O!0);
    f := function(x)
           if x eq 0 then return F!0; end if;
           v := Valuation(x, pid);
           error if v lt 0,
                 "Error in reduction modulo a prime ideal:",
                  x, "is not integral with respect to", pid;
           if v gt 0 then return F!0; end if;
           den := Denom(x);
           v := Valuation(den, p);
           y := x * e^v; // y is now p-integral and = x mod pid
           den := Denom(y);
           error if den mod p eq 0,
                 "Something´s wrong in MyResidueClassField in selmer.m!";
           return m(O!(den*y))/F!den;
         end function;
    return F, map< K -> F | x :-> f(x), y :-> y @@ m >;
  elif Type(pid) eq RngInt then
    O := Integers();
    K := FieldOfFractions(O);
    p := Generator(pid);
    // require IsPrime(p): "The ideal must be a prime ideal.";
    F, m := ResidueClassField(O, pid);
    f := function(x)
           if x eq 0 then return F!0; end if;
           v := Valuation(x, p);
           error if v lt 0,
                 "Error in reduction modulo a prime ideal:",
                  x, "is not integral with respect to", pid;
           if v gt 0 then return F!0; end if;
           return (F!Numerator(x))/F!Denominator(x);
         end function;
    return F, map< K -> F | x :-> f(x), y :-> y @@ m >;
  end if;
  error if false, "Unrecognized argument in MyResidueClassField";
end function;

// This fills in a `much needed gap' in Magma.
// This seems not to be a gap anymore? (&* will work on order ideals)
intrinsic mul(S::[RngOrdIdl]) -> RngOrdIdl
{}
    // How can I find the order from an empty sequence of RngOrdIdl's ???
    require not IsEmpty(S): "&* needs a non-empty sequence of ideals.";
    res := S[1];
    for i := 2 to #S do res := res*S[i]; end for;
    return res;
end intrinsic;

// Refers to old p-adic elements which may one day not be exported
// IsSquare exists for locals
/*
intrinsic IsSquare(x::RngPadElt) -> BoolElt, RngPadElt
{ Tests if x is a square. If so, return a square root as second value. }
    if x eq 0 then return true, Parent(x)!0; end if;
    v := Valuation(x);
    prec := AbsolutePrecision(x); // this breaks for exact x
    if v eq prec then return true, Sqrt(x); end if;
    if IsOdd(v) then return false, _; end if;
    v1 := v div 2;
    p := Prime(Parent(x));
    x1 := ExactQuotient(Integers()!x, p^v);
    if p eq 2 then
      case prec - v:
        when 1: return true, Sqrt(x);
        when 2: if x1 eq 1 then
                  return true, Sqrt(x);
                else
                  return false, _;
                end if;
        else if x1 mod 8 eq 1 then
            return true, Sqrt(x);
          else return false, _;
          end if;
      end case;
    else
      if IsSquare(GF(p)!x1) then
        return true, Sqrt(x);
      else
        return false, _;
      end if;
    end if;
end intrinsic;
*/

// First we need a means to construct groups
//  K(S, 3) = { a in K^*/(K^*)^3 | (a) in I_K^3 I_S }
// where K is a number field, I_K is the ideal group of K,
// S is a set of prime ideals of K, and I_S is the subgroup of I_K
// generated by S.


// a function that tries to siplify x in K mod cubes
function red(x)
  den, x1 := Denom(x);
  fd := Factorization(den);
  den1 := (&*[ Integers() | a[1]^((a[2]+2) div 3) : a in fd ])^3;
  x1 := (den1 div den)*x1;
  num := GCD(ChangeUniverse(Eltseq(x1), Integers()));
  fn := Factorization(num);
  num1 := (&*[ Integers() | a[1]^(a[2] div 3) : a in fn ])^3;
  return Parent(x1)![ ExactQuotient(c, num1) : c in Eltseq(x1) ];
end function;


intrinsic KS3(K::FldNum, S::[RngOrdIdl], h::Map : Bound := 0) -> ModTupFld, Map
{ Takes a number field K and a finite sequence of prime ideals of K
  and an involution of K and returns the GF(3)-vector space K(S,3) of
  elements mod cubes that generate an ideal supported in S (mod cubes
  of ideals), restricted to the kernel of x |-> x*h(x).   The
  second value is a map from the abstract vector space to elements
  of K representing the corresponding elements of K(S,3).
  The parameter Bound is passed on to the class group computation.}
    // The first step is the case where S is empty.
    // This involves the 3-torsion of the class group and units mod cubes.
    S1 := [ <Minimum(I), RamificationIndex(I), Degree(I)> : I in S ];
    vprintf JacHypSelmer, 2: "KS3 with\n  K = %o\n  S = %o\n", K, S1;
    O := Integers(K);
    // Set up map on ideals
    hI := func< I | ideal< O | [h(gen) : gen in Generators(I)]> >;
    vprintf JacHypSelmer, 2: " Computing class group...\n";
    Write("polys", DefiningPolynomial(K));
    if Bound eq 0 then
       C, mC := ClassGroup(O);
    else
       C, mC := ClassGroup(O : Bound := Bound);
    end if;
    Cinv := Invariants(C);
    vprintf JacHypSelmer, 2: "  Class group has invariants %o\n", Cinv;
    // Set up norm map
    if #C eq 1 then
      hC := hom< C -> C | x :-> C!0 >;
      CN := sub< C | >;
    else
      hC := hom< C -> C | [C.i + (hI(mC(C.i)) @@ mC) : i in [1..#Cinv]] >;
      CN := Kernel(hC);
    end if;
    CNinv := Invariants(CN);
    vprintf JacHypSelmer, 2: "  Kernel of norm has invariants %o\n", CNinv;
    c := 0;
    cgs0 := [ Codomain(mC) | ];
    cgs := [ Codomain(mC) | ];
    cis := [ Integers() | ];
    for i := 1 to #CNinv do
      if IsDivisibleBy(CNinv[i], 3) then
        c +:= 1;
        Append(~cis, i);
        I := mC((CNinv[i] div 3)*CN.i);
        Append(~cgs0, I*hI(I));
        Append(~cgs, I^3);
      end if;
    end for;
    vprintf JacHypSelmer, 2: "  3-torsion subgroup has rank %o\n", c;
    vprintf JacHypSelmer, 3: "  Generators are\n%o\n", cgs;
    vprintf JacHypSelmer, 2: " Computing unit group...\n";
    U, mU := UnitGroup(O);
    Uinv := Invariants(U);
    vprintf JacHypSelmer, 2: "  Unit group has invariants %o\n", Uinv;
    // Set up norm map
    vprintf JacHypSelmer, 2: "  Setting up norm map on units...\n";
    mat := Matrix([[Log(a) : a in AbsoluteValues(mU(U.i))] : i in [2..#Uinv]]);
    function mUinv(u)
      mm := Matrix([mat[i] : i in [1..#Uinv-1]]
                     cat [Vector([Log(a) : a in AbsoluteValues(u)])]);
      bas, trans := LLL(mm);
      min, pos := Min([&+[bas[i,j]^2 : j in [1..#Uinv]] : i in [1..#Uinv]]);
      error if Abs(trans[pos,#Uinv]) ne 1, "Failure 1 of mUinv!";
      u1 := -trans[pos,#Uinv]*&+[trans[pos,i]*U.(i+1) : i in [1..#Uinv-1]];
      uu1 := mU(u1);
      e := Order(U.1);
      z := mU(U.1);
      error if (uu1/u)^e ne 1, "Failure 2 of mUinv!";
      a := 0; while uu1 ne u do a +:= 1; uu1 := uu1*z; end while;
      return u1 + a*U.1;
    end function;
    // hU := hom< U -> U | [U.i + (h(mU(U.i)) @@ mU) : i in [1..#Uinv]] >;
    hU := hom< U -> U | [U.i + mUinv(h(mU(U.i))) : i in [1..#Uinv]] >;
    vprintf JacHypSelmer, 2: "  Computing kernel of norm map...\n";
    UN := Kernel(hU);
    UNinv := Invariants(UN);
    vprintf JacHypSelmer, 2: "  Kernel of norm has invariants %o\n", UNinv;
    if IsDivisibleBy(UNinv[1], 3) then // cannot happen
      u := #UNinv; ugs := [ K | mU(UN.i) : i in [1..#UNinv] ];
    else
      u := #UNinv - 1; ugs := [ K | mU(UN.i) : i in [2..#UNinv] ];
    end if;
    vprintf JacHypSelmer, 2: "  Ker N on U/U^3 has rank %o\n", u;
    vprintf JacHypSelmer, 3: "  Generators are\n%o\n", ugs;
    // We have an exact sequence
    //  0 --> U/U^3 --> K(empty, 3) --> C[3] --> 0 .
    // An element [I] of C[3] is lifted to K(empty, 3) by taking a
    // generator of I^3.
    dim0 := c + u;
    V0 := KSpace(GF(3), dim0);
    cgls := [ K | gen where _, gen := IsPrincipal(I) : I in cgs ];
    U3, U3q := quo< U | [3*U.i : i in [1..#Uinv]] >;
    hU3 := hom< U3 -> U3 | [U3q(hU(U3.i @@ U3q)) : i in [1..#Generators(U3)]] >;
    // Make sure to get elements with norms that are cubes
    // cgls := [ x*mU(-(U3q((x*h(x)/gen^3) @@ mU) @@ hU3) @@ U3q)
    cgls := [ x*mU(-(U3q(mUinv(x*h(x)/gen^3)) @@ hU3) @@ U3q)
               where _, gen := IsPrincipal(cgs0[i])
               where x := cgls[i]
               : i in [1..#cgls] ];
    require forall{x : x in cgls | IsPower(x*h(x), 3)}:
            "!!! Cl[3] generators are not in kernel of norm !!!";
    gs := ugs cat cgls;
    map0 := map< V0 -> K | v :-> &*[ gs[i]^(Integers()!vs[i]) : i in [1..dim0] ]
                                 where vs := Eltseq(v) >;
    // Now V0, map0 are the firt two values if S is empty.
    vprintf JacHypSelmer, 2: " Ker N on K(empty, 3) has dimension %o\n", dim0;
    if IsEmpty(S) then return V0, map0; end if;
    // We now have to extend K(empty, 3) to K(S, 3).
    // There is an exact sequence
    //  0 --> K(empty, 3) --> K(S, 3) -f-> GF(3)^S -g-> C/C^3 ,
    // where f is given by taking the ideal generated by a representative
    // element and looking at its S-component (mod 3), and g is induced
    // from the usual map I_S --> I_K --> C.
    // Hence we have to find the kernel of g and lift a basis back under f.
    V1 := KSpace(GF(3), #S);
    // Set up norm map
    indices := [ Position(S, hI(pid)) : pid in S ];
    require Seqset(indices) eq {1..#S}:
            "!!! Involution does not permute S !!!";
    hS0 := hom< V1 -> V1 | [V1.(indices[i]) : i in [1..#S]] >;
    hS := hom< V1 -> V1 | [V1.i + V1.(indices[i]) : i in [1..#S]] >;
    V1N := Kernel(hS);
    vprintf JacHypSelmer, 2: " Dimension of kernel of N on F^S = %o\n",
                             Dimension(V1N);
    V2 := KSpace(GF(3), c);
    h12 := hom< V1 -> V2 | [ V2 | Eltseq(CN!(I @@ mC))[cis] : I in S ] >;
    Bker := Basis(Kernel(h12) meet V1N);
    VK := KSpace(GF(3), #Bker);
    homK := hom< VK -> V1 | Bker >;
    s := #Bker;
    vprintf JacHypSelmer, 2:
       " The group of principal ideals mod cubes supported in S has rank %o\n",
       s;
    // To lift back an element I of the kernel, we have to find an ideal
    // J such that I*J^3 = (a) is principal; then a represents a preimage.
    dimS := dim0 + s;
    vprintf JacHypSelmer, 2: " ker N on K(S, 3) has dimension %o\n", dimS;
    if s eq 0 then return V0, map0; end if;
    if #Cinv eq 0 then
      sgs := [ <gen, gen*h(gen)/gen1^3>
                   where _, gen1 := IsPrincipal(I1)
                   where I1 := mul([S[i]^((bs[i]+bs1[i]) div 3) : i in [1..#S]])
                   where bs1 := ChangeUniverse(Eltseq(hS0(b)), Integers())
                   where _, gen := IsPrincipal(I)
                   where I := mul([ S[i]^bs[i] : i in [1..#S] ])
                   where bs := ChangeUniverse(Eltseq(b), Integers())
                  : b in Bker ];
      // now make sure to get elements with cube norms
      // sgs := [ a[1]*mU(-(U3q(a[2] @@ mU) @@ hU3) @@ U3q) : a in sgs ];
      sgs := [ a[1]*mU(-(U3q(mUinv(a[2])) @@ hU3) @@ U3q) : a in sgs ];
    else
      h3 := hom< C -> C | [3*C.i : i in [1..#Cinv]] >;
      sgs := [ <gen, I1*J1>
                   where J1 := hI(J)*J
                   where I1 := mul([S[i]^((bs[i]+bs1[i]) div 3) : i in [1..#S]])
                   where bs1 := ChangeUniverse(Eltseq(hS0(b)), Integers())
                   where _, gen := IsPrincipal(I*J^3)
                   where J := mC(-(I @@ mC) @@ h3)
                   where I := mul([ S[i]^bs[i] : i in [1..#S] ])
                   where bs := ChangeUniverse(Eltseq(b), Integers())
                  : b in Bker ];
      // now make sure to get elements with cube norms
      // first deal with obstruction in Cl[3], then with obstruction in units
      C3t := sub< C | [(Cinv[i] div 3) * C.i : i in cis] >;
      hC3t := hom<C3t -> C3t | [C3t!(hC(C3t.i)) : i in [1..#Generators(C3t)]]>;
      // sgs := [ gen1*mU(-(U3q(u @@ mU) @@ hU3) @@ U3q)
      sgs := [ gen1*mU(-(U3q(mUinv(u)) @@ hU3) @@ U3q)
               where u := gen1*h(gen1)/gen2^3  // a unit
               where _, gen2 := IsPrincipal(a[2]*I*hI(I))
               where gen1 := a[1]*gen
               where _, gen := IsPrincipal(I^3)
               where I := mC(-(C3t!(a[2] @@ mC)) @@ hC3t)
               : a in sgs ];
    end if;
    vprintf JacHypSelmer, 3: " Additional generators are\n%o\n", sgs;
    VS := KSpace(GF(3), dimS);
    gSs := gs cat sgs;
    require forall(g){g : g in gSs | IsPower(g*h(g), 3)} :
            "!!! Generator", g, "is not in the kernel of the norm !!!";
    mapS := map< VS -> K | v :-> &*[gSs[i]^(Integers()!vs[i]) : i in [1..dimS]]
                                 where vs := Eltseq(v) >;
    return VS, mapS;
end intrinsic;

// We represent a p-adic field by the following data:
// + a number field K;
// + a prime ideal in O, the integers of K, lying above p.
// From these, we deduce
// + a uniformizer;
// + a homomorphism K^* -->> finite elementary abelian 3-group with
//   kernel the elements that are cubes in the completion.

// Return the map  K^* -> a GF(3)-vector space
intrinsic MakeModCubes(K::FldNum, pid::RngOrdIdl) -> ModTupFld, Map
{}
  // (K::FldNum, pid::prime ideal in K) -> ModTupFld, Map
  O := Integers(K);
  p := Minimum(pid);
  e := RamificationIndex(pid, p);
  f := Degree(pid);
  _, pi := TwoElementNormal(pid);
  F, m := MyResidueClassField(pid);
  m0 := m;
  vprintf JacHypSelmer, 3: "MakeModCubes: p = %o, e = %o, f = %o\n", p, e, f;
  case #F mod 3:
    when 2:
      V := KSpace(GF(3), 1);
      h := map< K -> V | x :-> V![GF(3) | Valuation(x, pid)] >;
    when 1:
      V := KSpace(GF(3), 2); // the codomain of our homomorphism
      zeta := Roots(PolynomialRing(F)![F | 1,1,1])[1,1];
        // a cube root of unity in F
      if e*f eq 1 then
        pi := K!p;
        if Integers()!zeta gt Integers()!(zeta^2) then zeta := zeta^2; end if;
      end if;
      h := map< K -> V | x :-> V![GF(3) | v, z eq F!1 select 0
                                             else z eq zeta select 1 else 2]
                               where z := m(x2/(pi^v))^((#F-1) div 3)
                               where x2 := O![ c mod ex
                                               : c in ChangeUniverse(Eltseq(O!x1), Integers()) ]
                               where ex := p^(Floor(v/e) + 1)
                               where v := Valuation(x1, pid)
                               where x1 := red(x) >;
    when 0:
      // p = 3

      // c is a pid-adic cube and a pid-unit, but lies in all other
      // prime ideals above 3 in O.
      c := ChineseRemainderTheorem(pid^(Floor(3/2*e)+1), ideal<O | O!3>/pid^e,
                                   O!1, O!0);
      if e*f eq 1 then
        // Special case: completion is Q_3.
        // Fix the basis for Q_3^*/(Q_3^*)^3 to be 3, 1+3.
        // This is to ensure compatibility of bases between various
        //  degree 1 components of A_3; this is needed in the definition
        //  of the local descent map.
        R := quo< O | pid^2 >;
        sc := function(y) // y in K
                vprintf User2: "sc(<%o,%o,%o>, %o):\n", p, e, f, y;
                den, yd := Denom(y);
                vden := Valuation(den, 3);
                denp := ExactQuotient(den, 3^vden);
                _, denpi := XGCD(denp, 9); denpi := denpi mod 9;
                y := denpi*yd*3^((-vden) mod 3);
                vprintf User2: "  Red = %o\n", y;
                v := Valuation(y, pid);
                vprintf User2: "  Val = %o, mod 3 = %o\n", v, v mod 3;
                r := [GF(3) | v];
                y := y*3^(-v);
                while not IsIntegral(y) do y *:= c; end while;
                y := (R!y)^4;
                case y:
                  when R!1: Append(~r, GF(3)!0);
                  when R!4: Append(~r, GF(3)!1);
                  when R!7: Append(~r, GF(3)!2);
                  else error "!!! Something's wrong in sc(Q_3) !!!";
                end case;
                return r;
        end function;
        V := KSpace(GF(3), 2);
        h := map< K -> V | x :-> V!sc(x) >;
        return V, h;
      end if;

      // Our elementary abelian 3-group K_pid^*/(K_pid^*)^3
      // has rank (1 or 2) + e*f. 2 <==> there is a primitive cube root
      // of unity in K_pid.
      // We can check this as follows. First of all, e most be even.
      // Then, write 3 = pi^e*u with a pid-adic unit u. Consider the
      // additive map phi : F --> F, x |-> x*u + x^3. We have a cube root of 1
      // iff this map is not surjective.
      u := m((K!3)/(pi^e));
      vprintf JacHypSelmer, 3: "  u = %o\n", u;
      VF := KSpace(GF(3), f);
      FtoVF := map< F -> VF | x :-> VF!Eltseq(x) >;
      VFtoF := map< VF -> F | x :-> F!Eltseq(x) >;
      phi := hom< VF -> VF | [ FtoVF(x*u + x^3) where x := VFtoF(VF.i)
                               : i in [1..f] ] >;
      vprintf JacHypSelmer, 3: "  Image(phi) = %o\n", Basis(Image(phi));
      flag := IsEven(e) and Image(phi) ne VF;
      vprintf JacHypSelmer, 3:
          "  Cube root of unity is "*(flag select "" else "not ")*"present\n";
      // if flag then
        VQ, epi := quo< VF | Image(phi) >;
      // end if;
      dim := 1 + e*f + (flag select 1 else 0);
      V := KSpace(GF(3), dim);
      // A pid-unit is a cube in K_pid iff it is a cube in R.
      R := quo<O | pid^(Floor(3/2*e)+1)>;
      // reps is a lift to O of an F_3-basis of F.
      reps := [ R!((F![ i eq j select 1 else 0 : i in [1..f] ]) @@ m0)
                 : j in [1..f] ];
      // A basis of pid-units modulo cubes is given by
      //  [ 1 + r*pi^i : r in reps, i in s ] (cat [ unr ]) ,
      // where s := [ i : i in [1..Floor(3/2*e)] | not IsDivisibleBy(i, 3) ]
      // and where in case flag is set, unr = 1 + a*pi^(3/2*e)
      // such that the image of a in F is not in the image of phi above.
      // Together with pi itself, we get a basis of pid-adics modulo cubes.
      sc := function(y) // y in K
              vprintf User2: "sc(<%o,%o,%o>, %o):\n", p, e, f, y;
              den, yd := Denom(y);
              vden := Valuation(den, 3);
              denp := ExactQuotient(den, 3^vden);
              _, denpi := XGCD(denp, 9); denpi := denpi mod 9;
              y := denpi*yd*3^((-vden) mod 3);
              vprintf User2: "  Red = %o\n", y;
              v := Valuation(y, pid);
              vprintf User2: "  Val = %o, mod 3 = %o\n", v, v mod 3;
              r := [GF(3) | v];
              ex := 3^Ceiling((v + 1 + Floor(3*e/2))/e);
              y := O![ (Integers()!c) mod ex : c in Eltseq(y) ];
              vprintf User2: "  Redmod = %o\n", y;
              w := v;
              for j := 1 to v do
                y := (K!y)/pi;
                dy := Denom(y);
                if IsDivisibleBy(dy, 3) then
                  y *:= c^Valuation(dy, 3);
                  dy := Denom(y);
                  error if IsDivisibleBy(dy, 3),
                        "Something's wrong in MakeModCubes in 3descent.m!";
                end if;
                _, dyi := XGCD(dy, 9); // inverse of dy mod 3-adic cubes
                dyi := dyi mod 9;
                y := dyi*O!(dy*y);
                w -:= 1;
                ex := 3^Ceiling((w + 1 + Floor(3*e/2))/e);
                y := O![ (Integers()!c) mod ex : c in Eltseq(y) ];
              end for;
              vprintf User2: "  Integral, pid-unit y = %o\n", y;
              z := (R!y)^(3^f-1);
              vprintf User2: "  z = %o\n", z;
                // put it into 1 + pid; changes class to negative
              for i := 1 to Floor((3*e-1)/2) do
                // check:
                error if z ne 1 and Valuation(O!(z-1), pid) lt i,
                         "In modCubes:sc: i =", i, ", y =", y, ", z =", z,
                         "Valuation(z-1) =", Valuation(O!(z-1), pid);
                z1 := m((K!O!(z - 1))/pi^i);
                if IsDivisibleBy(i, 3) then
                  // Determine contribution of (1 + ?*pi^(i/3))^3
                  _, z2 := IsPower(z1, 3);
                  z *:= (1 - (R!(z2 @@ m0))*(R!pi)^(i div 3))^3;
                else
                  // Determine contribution of (1 + ?*pi^i)
                  seq := Eltseq(z1);
                  r cat:= seq;
                  for j := 1 to f do
                    if seq[j] ne 0 then
                      z *:= (1 - reps[j]*(R!pi)^i)^(Integers()!seq[j]);
                    end if;
                  end for;
                end if;
              end for;
              if flag then
                // Determine unramified contribution
                z1 := m((K!O!(z - 1))/pi^(3*e div 2));
                r cat:= Eltseq(epi(FtoVF(z1)));
              end if;
              vprintf User2: "  --> %o\n", r;
              return r;
            end function;
      h := map< K -> V | x :-> V!sc(x) >;
  end case;
  return V, h;
end intrinsic;

// This is an alternative way of producing a map K --> K(S,3), valid
// on elements of K^* representing something in K(S,3).
// We use a finite set of prime ideals and the local maps to K_p(3).
intrinsic MakeReverseMap(K::FldNum, KS3::ModTupFld, mKS3::Map) -> Map
{}
  vprintf JacHypSelmer, 3: "MakeReverseMap: K = %o\n", K;
  Vs := [];
  maps := [];
  Vms := [];
  O := Integers(K);
  prime := 3;
  repeat
    prime := NextPrime(prime);
    // while prime mod 3 ne 1 do prime := NextPrime(prime); end while;
    vprintf JacHypSelmer, 3: " Using prime %o\n", prime;
    pd := Decomposition(O, prime);
    pd := [ d[1] : d in pd ]; // The ideals
    vprintf JacHypSelmer, 3: " Decomposition: %o\n",
            [ <RamificationIndex(pid), Degree(pid)> : pid in pd ];
    for pid in pd do
      pV, pmap := MakeModCubes(K, pid);
      Append(~Vs, pV);
      Append(~maps, pmap);
      Append(~Vms, <pV, pmap>);
    end for;
    V := KSpace(GF(3), &+[Dimension(Vp) : Vp in Vs]);
    vprintf JacHypSelmer, 3: " Dimension of total space = %o\n", Dimension(V);
    mm := map< K -> V | x :-> V!(&cat[ Eltseq(m(x)) : m in maps ]) >;
    h := hom< KS3 -> V | [ mm(mKS3(v)) : v in Basis(KS3) ] >;
    vprintf JacHypSelmer, 3: " Kernel has dimension %o\n", Dimension(Kernel(h));
  until Dimension(Kernel(h)) eq 0;
  // Find small subset
  kers := [ Kernel(hom<KS3 -> t[1] | [t[2](mKS3(v)) : v in Basis(KS3)]>)
            : t in Vms ];
  _, pos := Min([ Dimension(k) : k in kers ]);
  maps := [ Vms[pos, 2] ];
  Vs :=[ Vms[pos, 1] ];
  V := KSpace(GF(3), &+[Dimension(Vp) : Vp in Vs]);
  mm := map< K -> V | x :-> V!(&cat[ Eltseq(m(x)) : m in maps ]) >;
  h := hom< KS3 -> V | [ mm(mKS3(v)) : v in Basis(KS3) ] >;
  ker := Kernel(h);
  while Dimension(ker) gt 0 do
    _, pos := Min([ Dimension(ker meet k) : k in kers ]);
    Append(~maps, Vms[pos, 2]);
    Append(~Vs, Vms[pos, 1]);
    V := KSpace(GF(3), &+[Dimension(Vp) : Vp in Vs]);
    mm := map< K -> V | x :-> V!(&cat[ Eltseq(m(x)) : m in maps ]) >;
    h := hom< KS3 -> V | [ mm(mKS3(v)) : v in Basis(KS3) ] >;
    ker := Kernel(h);
  end while;
  rmap := map< K -> KS3 | x :-> mm(x) @@ h >;
  return rmap;
end intrinsic;

// Determination of local images
// =============================
//
// For each bad prime p (i.e. p = 3 or 3|c_p), we have to find the image
// of E(Q_p)/3E(Q_p) in A_p^*/(A_p^*)^3. The map is given by
//  (x, y) |-> -2t(t-y) + (3s^2+a)(s-x)
// where  s^4 + 2as^2 + 4bs - a^2/3  = 0
// and    t^2 = s^3 + as + b .
// (A is the Q-algebra generated by s and t.)
//
// If (x,y) is a 3-torsion point, then A_p has two components Q_p with
// (s,t) = (x,y) and (s,t) = (x,-y), resp. The product of the images in
// these two components must be a cube. In the second component, the
// image is -(2t)^2 ~ (2t)^2, so we can take 2t as the image in the first
// component.
//
// If p is not 3, we can find the image from the Q_p-rational 3-power
// torsion.

function LocalImage(V, mAS3, f, pspec, sigma, tau, h, E)
  // + Elliptic curve is  E  over Q
  // + V is a GF(3)-vector space, mAS3 : V -> A a map giving representatives
  //   mod cubes for elements of V subset A(S, 3)
  // + pspec = <p, [ pid1, ..., pidk ]> a prime together with the ideals
  //   above it in A
  // + f : (x,y) --> A, the descent map
  p := pspec[1]; // the prime
  vprintf JacHypSelmer, 1: " Finding local image for p = %o ...\n", p;
  pids := pspec[2]; // the prime ideals
  // Set up map A^* --> prod_pid A_pid^*/(A_pid^*)^3
  A := Codomain(mAS3); O := Integers(A);
  // Find Q_p-rational 3-torsion
  pids0 := [ pid : pid in pids | RamificationIndex(pid) eq 1
                                 and Degree(pid) eq 1 ];
  hI := func< I | ideal< O | [h(gen) : gen in Generators(I)]> >;
  indices := [ Position(pids0, hI(pid)) : pid in pids0 ];
  error if Seqset(indices) ne {1..#pids0},
     "!!! Sonething is seriously wrong in LocalImage in file 3descent.m !!!";
  case #pids0:
    when 0:
      if p ne 3 then
        vprintf JacHypSelmer, 1:
                " No Q_%o-rational 3-torsion --> no restriction.\n", p;
        return V;
      else
        vprintf JacHypSelmer, 1: " dim E(Q_%o)[3] = 0\n", p;
        dim := 0;
        basis := [];
      end if;
    when 2:
      vprintf JacHypSelmer, 1: " dim E(Q_%o)[3] = 1\n", p;
      dim := 1;
      basis := [pids0[1]];
    when 8:
      vprintf JacHypSelmer, 1: " dim E(Q_%o)[3] = 2\n", p;
      dim := 2;
      basis := [pids0[1], pids0[indices[1] eq 2 select 3 else 2]];
    else error "Something is seriously wrong in LocalImage in file 3descent.m!";
  end case;
  // Find p-adic coordinates of points
  prec := 5 + Valuation(Discriminant(MinimalPolynomial(A.1)), p);
  Qp := pAdicField(p, prec); Rp := Integers(Qp);
  Rp1 := pAdicRing(p, Min(prec, 6)); // for output purposes
  PQp<X> := PolynomialRing(Qp : Global := false);
  PRp := PolynomialRing(Rp); PI := PolynomialRing(Integers());
  polp := PQp!PRp!PI!MinimalPolynomial(A.1);
  fact := [ fa[1] : fa in Factorization(polp) ];
  if GetVerbose("JacHypSelmer") ge 2 then
    printf " Defining polynomial of A factors ";
    printf "into degrees %o\n", [Degree(fa) : fa in fact];
  end if;
  roots := [ Rp!(-Coefficient(fa, 0)) : fa in fact | Degree(fa) eq 1 ];
  error if #fact ne #pids,
           "!!! There must be as many factors as there are prime ideals !!!";
  error if #roots ne #pids0,
      "!!! There must be as many zeroes as there are ideals with e*f = 1 !!!";
  // find root corresponding to pid
  broots := [ roots[pos]
              where _, pos := Max([Valuation(O!Integers()!root, pid)
                                    : root in roots])
              : pid in basis ];
  tors3 := [ <Rp!hpid(sigma), Rp!hpid(tau)>
             where hpid := hom< A -> Qp | Qp!root >
              : root in broots ];
  if GetVerbose("JacHypSelmer") ge 3 and dim gt 0 then
    printf " Q_%o-rational 3-torsion has basis\n", p;
    PrintBasis([<Rp1!t[1], Rp1!t[2]> : t in tors3], 3);
  end if;
  mcs := [ <W, mW> where W, mW := MakeModCubes(A, pid) : pid in pids ];
  if GetVerbose("JacHypSelmer") ge 3 then
    printf "  Dimensions of local components:\n";
    for i := 1 to #pids do
      printf "    pid<%o, %o> : Dim = %o\n",
             RamificationIndex(pids[i]), Degree(pids[i]), Dimension(mcs[i,1]);
    end for;
  end if;
  Vp := KSpace(GF(3), &+[ Dimension(mc[1]) : mc in mcs ]);
  maps := [ mc[2] : mc in mcs ];
  for i := 1 to #pids do
    pidi := pids[i];
    if RamificationIndex(pidi)*Degree(pidi) eq 1 then
      j := Position(pids, pids0[indices[Position(pids0, pidi)]]);
      pidj := pids[j];
      mapi := maps[i];
      mapj := maps[j];
      Vi := mcs[i,1];
      maps[i] := map< A -> Vi | x :-> Valuation(x, pidi) gt Valuation(x, pidj)
                                      select Vi!Eltseq(-mapj(x))
                                      else mapi(x) >;
    end if;
  end for;
  mAp := map< A -> Vp | x :-> Vp!(&cat[ Eltseq(m(x)) : m in maps ]) >;
  ff := func< x, y | mAp(f(x, y)) >;
  // Find image of 3-torsion in Vp
  Wp := sub< Vp | [ ff(Rationals()!b[1], Rationals()!b[2]) : b in tors3 ] >;
  if GetVerbose("JacHypSelmer") ge 2 then
    printf "  Local image of 3-torsion is";
    if Dimension(Wp) eq 0 then
      printf " trivial\n";
    else
      printf "\n";
      PrintBasis(Basis(Wp), 4);
    end if;
  end if;
  dim1 := p eq 3 select dim + 1 else dim;
  vprintf JacHypSelmer, 2: " Total local image has dimension %o\n", dim1;
  if Dimension(Wp) lt dim1 then
    vprintf JacHypSelmer, 2: " Need additional generator\n";
  end if;
  a1, a2, a3, a4, a6 := Explode(ChangeUniverse(aInvariants(E), Integers()));
  while Dimension(Wp) lt dim1 do
    x := Random(Rp);
    /*
    pol := PRp.1^2 + a1*x*PRp.1 + a3*PRp.1 - x^3 - a2*x^2 - a4*x - a6;
    roots := [ -Coefficient(fa[1], 0) : fa in Factorization(pol)
                                      | Degree(fa[1]) eq 1 ];
    if not IsEmpty(roots) then // have found a point in E(Q_p)
      y := Random(roots);
    */
    polb := a1*x + a3; polc := -(x^3 + a2*x^2 + a4*x + a6);
    flag, sqrt := IsSquare(polb^2 - 4*polc);
    if flag then
      y := Rp!((-polb + sqrt)/2);
      vprintf JacHypSelmer, 3: "  Found point (%o, %o)\n", Rp1!x, Rp1!y;
      // vprintf JacHypSelmer, 3: "    Check: v_%o(eqn(x,y)) = %o\n", p,
      //         Valuation(y^2 + polb*y + polc);
      im := ff(Rationals()!x, Rationals()!y);
      vprintf JacHypSelmer, 3: "   Image = %o\n", im;
      Wp := sub< Vp | Wp, im >;
      vprintf JacHypSelmer, 3: "   Dimension of image so far = %o\n",
                               Dimension(Wp);
    end if;
  end while;
  // Now we have found the local image
  if GetVerbose("JacHypSelmer") ge 2 then
    printf " Basis of local image:\n";
    PrintBasis(Basis(Wp), 3);
  end if;
  // Find inverse image in V
  homV := hom< V -> Vp | [ mAp(mAS3(V.i)) : i in [1..Dimension(V)] ] >;
  images := [ homV(V.i) : i in [1..Dimension(V)] ];
  if GetVerbose("JacHypSelmer") ge 2 then
    printf " Images of basis elements of V:\n";
    PrintBasis(images, 3);
  end if;
  VpQ, epi := quo< Vp | Wp >;
  homVV := hom< V -> VpQ | [ epi(im) : im in images ] >;
  V1 := Kernel(homVV);
  if GetVerbose("JacHypSelmer") ge 2 then
    printf " Basis of inverse image in V:\n";
    PrintBasis(Basis(V1), 3);
  end if;
  return V1;
end function;


function FindSubspace(V, pred)
  // Given a GF(3)-vector space V and some predicate pred on V such that
  // the subset of elements satisfying pred is a subspace, find this subspace.

  // Divide the underlying set of V into three subsets:
  // - the elements known to satisfy pred
  // - the elements known not to satisfy pred
  // - the rest
  vprintf JacHypSelmer, 3: "  FindSubspace, dim(V) = %o\n", Dimension(V);
  Syes := { V!0 };
  basis := [ V | ]; // a basis of Syes
  Sno := { V | };
  Srest := { V | v : v in V | v ne V!0 };
  while not IsEmpty(Srest) do
    // take some element of Srest and check if it satisfies pred
    vprintf JacHypSelmer, 3: "   #Syes = %o, #Sno = %o, #Srest = %o\n",
                             #Syes, #Sno, #Srest;
    v := Random(Srest);
    vprintf JacHypSelmer, 3: "   Chose random element %o\n", v;
    if pred(v) then
      // found a new basis element
      vprintf JacHypSelmer, 3: "    v is in subspace\n";
      Append(~basis, v);
      Syesnew := { w+v : w in Syes } join { w-v : w in Syes };
      Syes join:= Syesnew;
      Snonew := { v1+v2 : v1 in Sno, v2 in Syesnew };
      Sno join:= Snonew;
      Srest diff:= Syesnew;
      Srest diff:= Snonew;
    else
      // can eliminate v and -v
      vprintf JacHypSelmer, 3: "    v is not in subspace\n";
      S0 := {v, -v};
      Snonew := { v1+v2 : v1 in Syes, v2 in S0 };
      Sno join:= Snonew;
      Srest diff:= Snonew;
    end if;
  end while;
  vprintf JacHypSelmer, 3: "  Subspace has dimension %o\n", #basis;
  return sub< V | basis >;
end function;


// Now the main function...

intrinsic ThreeSelmerGroup(E::CrvEll : Bound := 0, Method := 2)
   -> RngIntElt, RngIntElt
{ Given an elliptic curve over the  rationals of the form y^2 = x^3 + ax + b
  with integers a and b, this computes the dimension of its 3-Selmer group.
  The first value returned is the Selmer rank (i.e., the bound for the
  Mordell-Weil rank one gets from the Selmer group), the second value is
  the GF(3)-dimension of the Selmer group itself.
  The parameter Bound is passed on to the class group computation.
  The parameter Method specifies how to deal with the global restriction
  involving the algebra B. If Method = 0, use class and unit groups and
  ideal factorisation. If Method = 1, use class and unit groups and
  reductions mod primes. If Method = 2, use test on cubes. }
    // Some checks
    require CoefficientRing(E) cmpeq Rationals():
            "E must be defined over the rationals.";
    aInvs := aInvariants(E);
    require forall{a : a in aInvs | IsIntegral(a)}:
            "E must have integral coefficients.";

    vprintf JacHypSelmer, 1: "Compute 3-Selmer group of %o\n", E;
    // Find the algebras A and B
    P<x> := PolynomialRing(Rationals() : Global := false);
    phi3 := DivisionPolynomial(E, 3)/3;
    while not forall{c : c in Coefficients(phi3) | IsIntegral(c)} do
      vprintf JacHypSelmer, 1:
              "3-division polynomial is not integral; scaling equation of E.\n";
      aInvs := [3*aInvs[1], 3^2*aInvs[2], 3^3*aInvs[3], 3^4*aInvs[4],
                3^6*aInvs[5]];
      E := EllipticCurve(aInvs);
      phi3 := DivisionPolynomial(E, 3)/3;
    end while;
    vprintf JacHypSelmer, 1: "3-division polynomial = %o\n", phi3;
    fact1 := [ f[1] : f in Factorization(phi3) ];
    if GetVerbose("JacHypSelmer") ge 2 then
      printf " 3-division polynomial factors into\n";
      PrintBasis(fact1, 3);
    end if;
    require #fact1 eq 1: "Currently, only the case with irreducible phi3",
                         "is implemented.";
    // Now phi3 is irreducible, monic and integral
    tors3 := 0; // Rank of rational 3-torsion
    Aplus := NumberField(phi3);
    sigma0 := Aplus.1;
    a1, a2, a3, a4, a6 := Explode(aInvs);
    // AVOID RELATIVE FIELDS !!!
    // PAplus<Y> := PolynomialRing(Aplus : Global := false);
    // A := ext< Aplus | Y^2 + a1*sigma0*Y + a3*Y
    //                    - sigma0^3 - a2*sigma0^2 - a4*sigma0 - a6 >;
    // sigma1 := A!sigma0;
    // tau1 := A.1;
    // autA := hom< A -> A | -tau1 - a1*sigma0 - a3 >; // pt |-> -pt
    // AA := AbsoluteField(A);
    // sigma := AA!sigma1;
    // tau := AA!tau1;
    // autAA := hom< AA -> AA | AA!autA(A!AA.1) >;
    P2<X> := PolynomialRing(P : Global := false); Y := P2!P.1;
    polA := Resultant(Evaluate(phi3, X),
                      Y^2 + a1*X*Y + a3*Y - X^3 - a2*X^2 - a4*X - a6);
    vprintf JacHypSelmer, 2: "  A is defined by\n   %o\n", polA;
    AA1 := NumberField(polA);
    tau1 := AA1.1;
    PA := PolynomialRing(AA1); XX := PA.1;
    pols := GCD(PA!phi3,
                tau1^2 + a1*XX*tau1 + a3*tau1 - XX^3 - a2*XX^2 - a4*XX - a6);
    require Degree(pols) eq 1: "!!! sigma not uniquely determined by tau !!!";
    sigma1 := -Coefficient(pols, 0);
    autAA1 := hom< AA1 -> AA1 | -tau1 - a1*sigma1 - a3 >;
    AA := OptimizedRepresentation(AA1);
    vprintf JacHypSelmer, 2: "  Using better representation of A, given by\n   %o\n",
            MinimalPolynomial(AA.1);
    flag, isoA := IsIsomorphic(AA1, AA);
    require flag: "!!! New field is not isomorphic to the old one !!!";
    _, isoA1 := IsIsomorphic(AA, AA1);
    require isoA1(isoA(AA1.1)) eq AA1.1 and isoA(isoA1(AA.1)) eq AA.1 :
            "!!! Isomorphisms are not inverse to each other !!!";
    sigma := isoA(sigma1);
    tau := isoA(tau1);
    autAA := hom< AA -> AA | isoA(autAA1(isoA1(AA.1))) >;
    OA := Integers(AA);
    // Define map f
    fx := a1*tau - 3*sigma^2 - 2*a2*sigma - a4;
    fy := 2*tau + a1*sigma + a3;
    fc := -fx*sigma - fy*tau;
    f0 := AA!red(fc^2 - a3*fc*fy - a6*fy^2);
    vprintf JacHypSelmer, 3:
            " f(x,y) = (%o)*((%o)*x + (%o)*y + %o)\n", f0, fx, fy, fc;
    f := func< x, y | f0*(fx*x + fy*y + fc) >;
    require IsPower(f0^2*((fx+fc)^2 - (a1+a3)*(fx+fc)*fy - (1+a2+a4+a6)*fy^2),
                    3): "!!! Something is wrong with f !!!";
    // Find the bad primes
    bad0 := BadPrimes(E);
    bad := [ 3 ];
    for p in bad0 do
      cp := TamagawaNumber(E, p);
      vprintf JacHypSelmer, 2: "  Bad prime %o --> c_%o = %o\n", p, p, cp;
      if IsDivisibleBy(cp, 3) and p ne 3 then Append(~bad, p); end if;
    end for;
    vprintf JacHypSelmer, 1: " Bad primes for 3-descent: %o\n", bad;
    badIdeals := [ <p, [ id[1] : id in Decomposition(OA, p) ]> : p in bad ];
    S := &cat[ bi[2] : bi in badIdeals ];
    // compute A(S,3)
    AS3, mAS3 := KS3(AA, S, autAA : Bound := Bound);
    vprintf JacHypSelmer, 1: " ker(N) on A(S, 3) has dimension %o\n",
                             Dimension(AS3);
    V := AS3;
    if Dimension(V) eq tors3 then return 0, tors3, V, mAS3; end if;
    for pspec in badIdeals do
      // Find local restrictions
      V := LocalImage(V, mAS3, f, pspec, sigma, tau, autAA, E);
      vprintf JacHypSelmer, 1:
              "  Local restriction at %o -> dimension = %o\n",
              pspec[1], Dimension(V);
      if Dimension(V) eq tors3 then return 0, tors3, V, mAS3; end if;
    end for;
    // Find second global restriction
    // Define field B
    // AVOID RELATIVE FIELDS !!!
    // B := ext< Aplus | Y^2 + a1*Y - a2 + Coefficient(phi3, 3) + sigma0 >;
    // autB := hom< B -> B | -B.1 >;
    // BA := AbsoluteField(B);
    // autBA := hom< BA -> BA | BA!autB(B!BA.1) >;
    // slope := BA!B.1;
    // sigmaB := BA!B!sigma0;
    polB := Resultant(Evaluate(phi3, X),
                      Y^2 + a1*Y - a2 + Coefficient(phi3, 3) + X);
    BA := NumberField(polB);
    slope := BA.1;
    sigmaB := -(slope^2 + a1*slope - a2 + Coefficient(phi3, 3));
    autBA := hom< BA -> BA | -slope - a1 >;
    // Define map to B
    linet := -1/(2*slope + a1)
              * (Coefficient(phi3, 2)
                  + sigmaB*(-slope^2 - a1*slope + a2) + a3*slope - a4);
    PB := PolynomialRing(BA);
    phi3red := ExactQuotient(PB!phi3, PB.1 - sigmaB);
    sx1 := -Coefficient(phi3red, 2);
    sx2 := Coefficient(phi3red, 1);
    sx3 := -Coefficient(phi3red, 0);
    sy1 := slope*sx1 + 3*linet;
    sy2 := slope^2*sx2 + 2*slope*linet*sx1 + 3*linet^2;
    sy3 := slope^3*sx3 + slope^2*linet*sx2 + slope*linet^2*sx1 + linet^3;
    mat := MatrixAlgebra(BA, 3)![ BA | 0,1,0, 0,0,1, sy3,-sy2,sy1 ];
    require AA1.1 eq tau1: "!!! AA1.1 and tau1 are distinct !!!";
    AtoB := func< x | Determinant(Evaluate(P!Eltseq(isoA1(x)), mat)) >;
    if Method ne 2 then
      OB := Integers(BA);
      SB := &cat[ [s[1] : s in Decomposition(OB, p)] : p in bad ];
      BS3, mBS3 := KS3(BA, SB, autBA : Bound := Bound);
      mrBS3 := MakeReverseMap(BA, BS3, mBS3);
      vprintf JacHypSelmer, 3: " Images under A->B of basis of V:\n";
      homGlob2 := hom< V -> BS3 | [mrBS3(AtoB(mAS3(v))) : v in Basis(V)] >;
      if GetVerbose("JacHypSelmer") ge 3 then
        for v in Basis(V) do
          printf "  v = %o |--> %o\n", v, homGlob2(v);
        end for;
      end if;
      V := Kernel(homGlob2);
    else
      V := FindSubspace(V, func< v | IsPower(AtoB(mAS3(v)), 3) >);
    end if;
    vprintf JacHypSelmer, 1:
          " Restriction to kernel of map to B^*/(B^*)^3 -> dimension = %o\n",
            Dimension(V);
    return Dimension(V) - tors3, Dimension(V), V, mAS3;
end intrinsic;
