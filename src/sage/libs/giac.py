from sage.structure.proof.all import polynomial as proof_polynomial
from sage.misc.sage_eval import sage_eval
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

def groebner_basis_libgiac(gens, epsilon=None, prot=False, *args, **kwds):

    if not (isinstance(gens,list)):
        raise ValueError("Fisrt argument must be a list")
    try:
       from giacpy import libgiac,giacsettings

    except ImportError:
            raise ImportError("""One of the optional packages giac or giacpy is missing""")

    # save the original giac config
    proba_epsilon_orig=giacsettings.proba_epsilon
    debug_orig=libgiac('debug_infolevel()')

    # get the ring from gens
    P = iter(gens).next().parent()
    K = P.base_ring()
    p = K.characteristic()

    if (  K.is_prime_field()  and  p == 0 ):

          F = libgiac(gens)

    elif ( K.is_prime_field()  and ( p < 2**31 ) ):

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

    # prot
    if (prot==True):
        libgiac('debug_infolevel(2)')

    # compute de groebner basis with giac
    gb_giac = F.gbasis([P.gens()], "revlex")

    # restore giacsettings
    libgiac('debug_infolevel(0)')
    giacsettings.proba_epsilon = proba_epsilon_orig


    # TODO: we shouldn't use string parsing here
    #gens_dict = P.gens_dict()
    #gb = []
    #for f in gb_giac:
    #    gb.append(sage_eval(str(f), gens_dict))
    #
    #return PolynomialSequence( P, gb)
    # ?why not (it seems faster was your method better for ram?):
    return PolynomialSequence( gb_giac, P ) # should it be , immutable=True ?
