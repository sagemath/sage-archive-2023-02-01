from sage.structure.proof.all import polynomial as proof_polynomial
from sage.misc.sage_eval import sage_eval
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

def groebner_basis_libgiac(gens, epsilon=None, prot=False, *args, **kwds):
    from giacpy import libgiac, giacsettings
    F = libgiac(gens)
    if epsilon is None:
        if proof_polynomial():
            giacsettings.proba_epsilon = 0
        else:
            giacsettings.proba_epsilon = 1e-15
    else:
        giacsettings.proba_epsilon = epsilon

    P = iter(gens).next().parent()
    # TODO: check permitted base rings
    if P.term_order() != "degrevlex":
        raise NotImplementedError("Only degrevlex term orderings are supported by Giac's Gröbner basis engine.")

    gb_giac = F.gbasis([P.gens()], "revlex")
    # TODO: we shouldn't use string parsing here
    gens_dict = P.gens_dict()
    gb = []
    for f in gb_giac:
        gb.append(sage_eval(str(f), gens_dict))

    return PolynomialSequence(P, gb)
