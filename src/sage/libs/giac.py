from sage.structure.proof.all import polynomial as proof_polynomial
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence


def groebner_basis_libgiac(gens, epsilon=None, prot=False, *args, **kwds):
    """
    """
    try:
        from giacpy import libgiac, giacsettings
    except ImportError:
        raise ImportError("""One of the optional packages giac or giacpy is missing""")

    proba_epsilon_orig = giacsettings.proba_epsilon  # save the original giac config

    # get the ring from gens
    P = iter(gens).next().parent()
    K = P.base_ring()
    p = K.characteristic()

    if K.is_prime_field() and p == 0:
        F = libgiac(gens)
    elif K.is_prime_field() and p < 2**31:
        F = (libgiac(gens) % p)
    else:
        raise NotImplementedError("only prime fields of cardinal < 2^31 are implemented in giac for Groebner Basis")

    if P.term_order() != "degrevlex":
        raise NotImplementedError("Only degrevlex term orderings are supported by Giac's Groebner basis engine.")

    # proof or probabilistic reconstruction
    if epsilon is None:
        if proof_polynomial():
            giacsettings.proba_epsilon = 0
        else:
            giacsettings.proba_epsilon = 1e-15
    else:
        giacsettings.proba_epsilon = epsilon

    # compute de groebner basis with giac
    gb_giac = F.gbasis([P.gens()], "revlex")

    return PolynomialSequence(gb_giac, P, immutable=True)
