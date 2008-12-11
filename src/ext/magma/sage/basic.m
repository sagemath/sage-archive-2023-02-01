// -*- magma -*-

function PreparseElts(R)
    if Type(R) eq RngInt then
       return false;
    end if;
    return true;
end function;

intrinsic Sage(X::.) -> MonStgElt, BoolElt
{Default way to convert a Magma object to Sage if we haven't
written anything better.}
    return Sprintf("%o", X), true;
end intrinsic;

intrinsic Sage(X::SetEnum) -> MonStgElt, BoolElt
{Convert an enumerated set to Sage.}
    Y := [Sage(z) : z in X];
    return Sprintf("Set(%o)", Y), true;
end intrinsic;

intrinsic Sage(X::SetIndx) -> MonStgElt, BoolElt
{Convert an indexed set to Sage.
 WARNING: Sage does not have an analogue of indexed sets (yet!),
 so we just return a Python list.}
    Y := [Sage(z) : z in X];
    return Sprintf("%o", Y), true;
end intrinsic;

intrinsic Sage(X::SetMulti) -> MonStgElt, BoolElt
{Convert a multiset to Sage.
 WARNING: Sage does not have an analogue of multisets yet, so we return a Python list.}
    Y := [Sage(z) : z in X];
    return Sprintf("%o", Y), true;
end intrinsic;

intrinsic Sage(X::RngInt) -> MonStgElt, BoolElt
{Conver the ring of integers to Sage.}
    return "ZZ", false;
end intrinsic;

intrinsic Sage(X::FldRat) -> MonStgElt, BoolElt
{}
    return "QQ", false;
end intrinsic;

intrinsic Sage(X::RngIntElt) -> MonStgElt, BoolElt
{}
    return Sprintf("Integer('%h')", X), false;
end intrinsic;

/* Matrices */

function convert_matrix(X, preparse_entries)
    if preparse_entries then
        return Sprintf("matrix(%o, %o, %o, %o)", Sage(BaseRing(X)),
                    Nrows(X), Ncols(X), [Sage(y) : y in Eltseq(X)]);
    else
        return Sprintf("matrix(%o, %o, %o, %o)", Sage(BaseRing(X)),
                    Nrows(X), Ncols(X), Eltseq(X));
    end if;
end function;

intrinsic Sage(X::AlgMatElt) -> MonStgElt, BoolElt
{}
    pp := PreparseElts(BaseRing(X));
    return convert_matrix(X, pp), pp;
end intrinsic;

intrinsic Sage(X::ModMatRngElt) -> MonStgElt, BoolElt
{}
    pp := PreparseElts(BaseRing(X));
    return convert_matrix(X, pp), pp;
end intrinsic;


intrinsic SageCreateWithNames(X::., names::.) -> .
{Assign the given names to the object X, then return X.}
    AssignNames(~X, names);
    return X;
end intrinsic;

/* Finite fields */

intrinsic Sage(X::FldFin) -> MonStgElt, BoolElt
{}
  if IsPrimeField(X) then
    return Sprintf("GF(%o)", Characteristic(X)), false;
  else
    return Sprintf("GF(%o, '%o'.replace('$.', 'x').replace('.', ''), modulus=%o)", #X, X.1, Sage(DefiningPolynomial(X))), false;
  end if;
end intrinsic;

intrinsic Sage(X::FldFinElt) -> MonStgElt, BoolElt
{}
  return Sprintf("%o(%o)", Sage(Parent(X)), Sage(Polynomial(Eltseq(X)))), false;
end intrinsic;

/* Approximate fields */

intrinsic Sage(X::FldRe) -> MonStgElt, BoolElt
{}
    return Sprintf("RealField(%o)", Precision(X : Bits := true)), false;
end intrinsic;

intrinsic Sage(X::FldCom) -> MonStgElt, BoolElt
{}
    return Sprintf("ComplexField(%o)", Precision(X : Bits := true)), false;
end intrinsic;

intrinsic Sage(X::FldReElt) -> MonStgElt, BoolElt
{}
    return Sprintf("%o(%o)", Sage(Parent(X)), X), false;
end intrinsic;

intrinsic Sage(X::FldComElt) -> MonStgElt, BoolElt
{}
  return Sprintf("%o([%o, %o])", Sage(Parent(X)), Sage(Real(X)), Sage(Imaginary(X))), false;
end intrinsic;

/* Polynomials */

intrinsic SageNamesHelper(X::.) -> MonStgElt
{}
  /* XXX */
  i := NumberOfNames(X);
  if i ge 2 then
      return (&* [ Sprintf("%o, ", X.j) : j in [ 1..i-1 ] ]) * Sprintf("%o", X.i);
  else
      return Sprintf("%o", X.i);
  end if;
end intrinsic;

intrinsic Sage(X::RngUPol) -> MonStgElt, BoolElt
{}
  return Sprintf("%o['%o'.replace('$.', 'x').replace('.', '')]", Sage(BaseRing(X)), SageNamesHelper(X)), false;
end intrinsic;

intrinsic Sage(X::RngUPolElt) -> MonStgElt, BoolElt
{}
  return Sprintf("%o(%o)", Sage(Parent(X)), Sage(Coefficients(X))), false;
end intrinsic;

intrinsic Sage(X::RngMPol) -> MonStgElt, BoolElt
{}
  return Sprintf("%o['%o'.replace('$.', 'x').replace('.', '')]", Sage(BaseRing(X)), SageNamesHelper(X)), false;
end intrinsic;

/* intrinsic Sage(X::RngMPolElt) -> MonStgElt, BoolElt */
/*   {} */
/*   /\* XXX - this doesn't work quite yet *\/ */
/*   return Sprintf("%o(%o)", Sage(Parent(X)), Sage(Coefficients(X))), false; */
/* end intrinsic; */

/* Elliptic curves */

intrinsic Sage(X::CrvEll) -> MonStgElt, BoolElt
{}
  as := aInvariants(X);
  return Sprintf("EllipticCurve(%o)", Sage(as)), true;
end intrinsic;

/* Hyperelliptic curves */

intrinsic Sage(X::CrvHyp) -> MonStgElt, BoolElt
{}
  f, g := HyperellipticPolynomials(X);
  return Sprintf("HyperellipticCurve(%o, %o)", Sage(f), Sage(g)), true;
end intrinsic;

/* Modules and vector spaces */

intrinsic Sage(X::ModTupRng) -> MonStgElt, BoolElt
{}
  if IsIdentity(InnerProductMatrix(X)) then
      return Sprintf("FreeModule(%o, %o)", Sage(BaseRing(X)), Sage(Rank(X))), true;
  else
      return Sprintf("FreeModule(%o, %o, inner_product_matrix=%o)", Sage(BaseRing(X)), Sage(Rank(X)), Sage(InnerProductMatrix(X))), true;
  end if;
end intrinsic;

intrinsic Sage(X::ModTupRngElt) -> MonStgElt, BoolElt
{}
    return Sprintf("%o(%o)", Sage(Parent(X)), Sage(ElementToSequence(X))), true;
end intrinsic;
