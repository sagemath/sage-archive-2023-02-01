"""
Common combinatorial tools


REFERENCES:

.. [NCSF] Gelfand, Krob, Lascoux, Leclerc, Retakh, Thibon,
   *Noncommutative Symmetric Functions*, Adv. Math. 112 (1995), no. 2, 218-348.

"""

from sage.misc.misc_c import prod
from sage.functions.other import factorial

# The following might call for defining a morphism from ``structure
# coefficients'' / matrix using something like:
# Complete.module_morphism( coeff = coeff_pi, codomain=Psi, triangularity="finer" )
# the difficulty is how to best describe the support of the output.

def coeff_pi(J,I):
    r"""
    Returns the coefficient `\pi_{J,I}` as defined in [NCSF]_.

    INPUT:

    - ``J`` -- a composition
    - ``I`` -- a composition refining ``J``

    OUTPUT:

    - integer

    EXAMPLES::

        sage: from sage.combinat.ncsf_qsym.combinatorics import coeff_pi
        sage: coeff_pi(Composition([1,1,1]), Composition([2,1]))
        2
        sage: coeff_pi(Composition([2,1]), Composition([3]))
        6
    """
    return prod(prod(K.partial_sums()) for K in J.refinement_splitting(I))

def coeff_lp(J,I):
    r"""
    Returns the coefficient `lp_{J,I}` as defined in [NCSF]_.

    INPUT:

    - ``J`` -- a composition
    - ``I`` -- a composition refining ``J``

    OUTPUT:

    - integer

    EXAMPLES::

        sage: from sage.combinat.ncsf_qsym.combinatorics import coeff_lp
        sage: coeff_lp(Composition([1,1,1]), Composition([2,1]))
        1
        sage: coeff_lp(Composition([2,1]), Composition([3]))
        1
    """
    return prod(K[-1] for K in J.refinement_splitting(I))

def coeff_ell(J,I):
    r"""
    Returns the coefficient `\ell_{J,I}` as defined in [NCSF]_.

    INPUT:

    - ``J`` -- a composition
    - ``I`` -- a composition refining ``J``

    OUTPUT:

    - integer

    EXAMPLES::

        sage: from sage.combinat.ncsf_qsym.combinatorics import coeff_ell
        sage: coeff_ell(Composition([1,1,1]), Composition([2,1]))
        2
        sage: coeff_ell(Composition([2,1]), Composition([3]))
        2
    """
    return prod(map(len, J.refinement_splitting(I)))

def coeff_sp(J,I):
    r"""
    Returns the coefficient `sp_{J,I}` as defined in [NCSF]_.

    INPUT:

    - ``J`` -- a composition
    - ``I`` -- a composition refining ``J``

    OUTPUT:

    - integer

    EXAMPLES::

        sage: from sage.combinat.ncsf_qsym.combinatorics import coeff_sp
        sage: coeff_sp(Composition([1,1,1]), Composition([2,1]))
        2
        sage: coeff_sp(Composition([2,1]), Composition([3]))
        4
    """
    return prod(factorial(len(K))*prod(K) for K in J.refinement_splitting(I))

