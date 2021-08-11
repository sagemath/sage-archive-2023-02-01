from sage.modular.quasimodform.ring import QuasiModularForms

from sage.rings.all import Integer

def Hn(n):
    r"""
    Return the `n`-th `H_n` function defined by the recurrence relation
    `4n(n-1)H_n = 8D(H_n(n-2)) + PH_n(n-2)`.

    INPUT:

    - `n` - (int, Integer) a nonnegative integer.

    EXAMPLES::

        sage: from sage.combinat.shifted_symmetric_functions.q_bracket import Hn
        sage: Hn(0)
        1
        sage: Hn(1)
        0
        sage: Hn(2)
        1/24 - q - 3*q^2 - 4*q^3 - 7*q^4 - 6*q^5 + O(q^6)
        sage: Hn(2).to_polynomial('P, Q, R')
        1/24*P

    OUTPUT: A quasimodular form defined by the recurrence relation above.

    TESTS::

        sage: Hn('n')
        Traceback (most recent call last):
        ...
        ValueError: n (=n) must be a nonnegative integer
        sage: Hn(-1)
        Traceback (most recent call last):
        ...
        ValueError: n (=-1) must be a nonnegative integer
    """
    if not isinstance(n, (int, Integer)) or n < 0:
        raise ValueError("n (=%s) must be a nonnegative integer" % (n))
    QM = QuasiModularForms(1)
    if n == 0:
        H = QM.one()
    elif n%2 or n == 1:
        H = QM.zero()
    else:
        P, Q, R = QM.gens()
        RR = QM.base_ring()
        u = RR(4 * n * (n + 1)).inverse_of_unit()
        H = u * (RR(8) * Hn(n - 2).derivative() + P * Hn(n - 2))
    return H

# def q_bracket(polynomial, output_as_polynomial=True, names='P, Q, R'):
