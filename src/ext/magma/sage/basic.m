intrinsic Sage(X::.) -> MonStgElt
{}
    return Sprintf("%o", X);
end intrinsic;

intrinsic Sage(X::SetEnum) -> MonStgElt
{}
    Y := [Sage(z) : z in X];
    return Sprintf("Set(%o)", Y);
end intrinsic;

intrinsic Sage(X::SetIndx) -> MonStgElt
{WARNING: Sage does not have an analogue of indexed sets.}
    Y := [z : z in X];
    return Sprintf("%o", Y);
end intrinsic;

intrinsic Sage(X::SetMulti) -> MonStgElt
{WARNING: Sage does not have an analogue of multisets.}
    Y := [z : z in X];
    return Sprintf("%o", Y);
end intrinsic;

intrinsic Sage(X::RngInt) -> MonStgElt
{}
    return "ZZ";
end intrinsic;

function convert_matrix(X)
    return Sprintf("matrix(%o, %o, %o, %o)", Sage(BaseRing(X)),
                    Nrows(X), Ncols(X), [Sage(y) : y in Eltseq(X)]);
end function;

intrinsic Sage(X::AlgMatElt) -> MonStgElt
{}
    return convert_matrix(X);
end intrinsic;

intrinsic Sage(X::ModMatRngElt) -> MonStgElt
{}
    return convert_matrix(X);
end intrinsic;
