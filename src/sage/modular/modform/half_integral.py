r"""
Compute spaces of half-integral weight modular forms.

AUTHORS:
    -- William Stein, 2007-06-03

from sage.rings.all import Integer

def theta2_qexp(prec=10, var='q'):
    prec = Integer(prec)
    if prec <= 0:
        raise ValueError, "prec must be positive"
    v = [Integer(0)] * prec
    one = Integer(1)
    for m in range(1,floor(sqrt(prec))+2,2):
        v[m*m] = one
    return ZZ[[var]](v)

def theta3_qexp(prec=10, var='q'):
function Theta3(prec)
   q := PowerSeriesRing(Integers()).1;
   return 1 + 2*&+[q^(m^2) : m in [1..Floor(Sqrt(prec))+1]] + O(q^prec);
end function;


intrinsic HalfIntegralWeightForms(chi::GrpDrchElt,
                                    k::RngIntElt,
                                 prec::RngIntElt) -> SeqEnum
{A basis for the space of weight k/2 forms with character chi.
The level of chi must be divisible by 16 and k must be odd and >1.
I don't know what the minimal allowable prec value is, because
I haven't implemented Cohen's algorithm for the dimension of
the space!}

/* Basmaji gives the following algorithm on page 55 of his thesis.

   Let S = S(eps, (k+1)), where eps = chi*psi^((k+1)/2), where
   psi is the nontrivial mod-4 Dirichlet character.
   Let U be the subspace of S x S of elements (a,b) such that
   Theta2*a = Theta3*b.
   Then U is isomorphic to S(chi, k/2) via the map
          (a,b) |----> a/Theta2.
*/

   require Modulus(chi) mod 16 eq 0 :
       "The modulus of argument 1 must be divisible by 16.";
   require IsOdd(k) : "Argument 2 must be odd.";
   psi := DirichletGroup(4).1;
   eps := chi*psi^((k+1) div 2);
   M := BaseExtend(ModularForms([eps], (k+1) div 2),Rationals());
   S := [PowerSeries(f,prec) : f in Basis(CuspidalSubspace(M))];
   T2 := Theta2(prec);
   T3 := Theta3(prec);
   A := RMatrixSpace(BaseRing(eps),2*#S,prec)!0;
   for i in [1..#S] do
      T2f := T2*S[i];
      T3f := T3*S[i];
      for j in [0..prec-1] do
         A[i,j+1] := Coefficient(T2f,j);
         A[#S+i,j+1] := -Coefficient(T3f,j);
      end for;
   end for;
   B := Basis(Kernel(A));
   avec := [&+[b[i]*S[i] : i in [1..#S]] : b in B];  // sequence of a's
   if #avec eq 0 then
      return [];
   end if;
   R<q> := FieldOfFractions(Parent(avec[1]));
   return [(R!a)/(R!T2) : a in avec];
end intrinsic;

/* Example:

> time B := HalfIntegralWeightForms(DirichletGroup(592)!1,3,200);
Time: 18.920
> #B;
34

*/
"""
