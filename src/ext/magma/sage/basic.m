function PreparseElts(R)
    if Type(R) eq RngInt then
       return false;
    end if;
    return true;
end function;

intrinsic Sage(X::.) -> MonStgElt
{}
    return Sprintf("%o", X), true;
end intrinsic;

intrinsic Sage(X::SetEnum) -> MonStgElt
{}
    Y := [Sage(z) : z in X];
    return Sprintf("Set(%o)", Y), true;
end intrinsic;

intrinsic Sage(X::SetIndx) -> MonStgElt
{WARNING: Sage does not have an analogue of indexed sets.}
    Y := [z : z in X];
    return Sprintf("%o", Y), true;
end intrinsic;

intrinsic Sage(X::SetMulti) -> MonStgElt, BoolElt
{WARNING: Sage does not have an analogue of multisets.}
    Y := [z : z in X];
    return Sprintf("%o", Y), true;
end intrinsic;

intrinsic Sage(X::RngInt) -> MonStgElt, BoolElt
{}
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
