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
