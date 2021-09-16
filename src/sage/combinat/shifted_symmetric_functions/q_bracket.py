from sage.modular.quasimodform.ring import QuasiModularForms

from sage.rings.all import Integer
from sage.rings.polynomial.multi_polynomial_element import MPolynomial

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

def _compute_derivative(polynomial):
    r"""
    Compute the shifted symmetric derivative of the given polynomial.

    INPUT:

    - ``polynomial`` -- a multivariate polynomial over `\QQ`.

    TESTS::

        sage: from sage.combinat.shifted_symmetric_functions.q_bracket import _compute_derivative
        sage: P.<Q1, Q2, Q3, Q4, Q5, Q6> = QQ[]
        sage: _compute_derivative(Q6)
        Q5
        sage: _compute_derivative(Q5)
        Q4
        sage: _compute_derivative(Q1)
        1
        sage: _compute_derivative(Q1 + Q3 + Q5)
        Q2 + Q4 + 1
        sage: _compute_derivative(Q3 * Q5)
        Q3*Q4 + Q2*Q5
    """
    if not isinstance(polynomial, MPolynomial):
        raise ValueError("the given polynomial must be a multivariate polynomial")
    poly_parent = polynomial.parent()
    if polynomial in poly_parent.base_ring():
        return poly_parent.zero()
    poly_gens = poly_parent.gens()
    mon_list = polynomial.monomials()
    der = poly_parent.zero()

    # compute the derivative for every monomials
    for m in mon_list:
        coeff = polynomial.monomial_coefficient(m)
        degrees_list = list(m.degrees())
        der_m = poly_parent.zero()
        for i in range(poly_parent.ngens()):
            if degrees_list[i] != 0:
                product_gen = poly_parent.one()
                for j in range(poly_parent.ngens()):
                    d = degrees_list[j]
                    if d != 0:
                        if j == 0:
                            product_gen *= d * poly_gens[j] ** (d - 1)
                        elif j == i:
                            product_gen *= d * poly_gens[j] ** (d - 1) * poly_gens[j - 1]
                        else:
                            product_gen *= poly_gens[j] ** d
                der_m += product_gen
        der += coeff * der_m
    return der

def shifted_symmetric_derivative(polynomial, order=1):
    r"""
    Return the shifted symmetric derivative of the given polynomial.

    INPUT:

    - ``polynomial`` -- a multivariate polynomial over `\QQ`.
    - ``order`` (Integer, default: 1) -- the order of the derivative.

    EXAMPLES::

        sage: from sage.combinat.shifted_symmetric_functions.q_bracket import shifted_symmetric_derivative
        sage: P.<Q1, Q2, Q3, Q4, Q5, Q6> = QQ[]
        sage: shifted_symmetric_derivative(Q2)
        Q1
        sage: shifted_symmetric_derivative(Q2 * Q4)
        Q2*Q3 + Q1*Q4
        sage: shifted_symmetric_derivative(Q4^5)
        5*Q3*Q4^4
        sage: shifted_symmetric_derivative(Q2^3, order=2)
        6*Q1*Q2 + 3*Q2^2
    """
    der_n = polynomial
    for n in range(order):
        der_n = _compute_derivative(der_n)
    return der_n


def _compute_q_bracket(polynomial):
    r"""
    Compute the `q`-bracket of the given (shifted symmetric) polynomial.

    TESTS::

        sage: from sage.combinat.shifted_symmetric_functions.q_bracket import _compute_q_bracket
        sage: P.<Q1, Q2, Q3, Q4, Q5, Q6> = QQ[]
        sage: _compute_q_bracket(Q2)
        -1/24 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
        sage: _compute_q_bracket(Q2 ** 2)
        7/2880 + 1/12*q + 9/4*q^2 + 31/3*q^3 + 343/12*q^4 + 117/2*q^5 + O(q^6)
        sage: _compute_q_bracket(Q4)
        7/5760 + 1/24*q + 9/8*q^2 + 31/6*q^3 + 343/24*q^4 + 117/4*q^5 + O(q^6)
    """
    if not isinstance(polynomial, MPolynomial):
        raise ValueError("polynomial must be a multivariate polynomial")
    QM = QuasiModularForms(1)
    poly_parent = polynomial.parent()
    if polynomial.is_constant():
        return QM(polynomial.base_ring()(polynomial))
    if polynomial % poly_parent.gen(0) == 0:
        return QM.zero()
    weighted_deg = polynomial.weighted_degree(list(range(1, poly_parent.ngens() + 1)))
    bracket = QM.zero()
    for n in range(2, weighted_deg + 1, 2):
        der = shifted_symmetric_derivative(polynomial, order=n)
        bracket += sum(-Hn(n) * der.monomial_coefficient(m) * _compute_q_bracket(m) for m in der.monomials())
    return bracket

def q_bracket(polynomial, output_as_polynomial=True, names='P, Q, R'):
    r"""
    Return the `q`-bracket of the given (shifted symmetric) polynomial.

    If ``output_as_polynomial`` is set to ``True``, then the output will be a
    multivariate polynomial in the generators of the ring of quasimodular forms
    for `\mathrm{SL}_2(\ZZ)`.

    INPUT:

    - ``polynomial`` -- a multivariate polynomial over `\QQ`
    - ``output_as_polynomial`` (default: True) -- a boolean
    - ``names`` (default: 'P, Q, R') -- a string of names representing the
      the generators of the quasimodular forms ring.

    EXAMPLES::

        sage: from sage.combinat.shifted_symmetric_functions.q_bracket import q_bracket
        sage: P.<Q1, Q2, Q3, Q4, Q5, Q6> = QQ[]
        sage: q_bracket(Q2)
        -1/24*P
        sage: q_bracket(Q3)
        0
        sage: q_bracket(Q4)
        1/1152*P^2 + 1/2880*Q
        sage: q_bracket(Q2 ** 2)
        1/576*P^2 + 1/1440*Q
        sage: q_bracket(Q2, output_as_polynomial=False)
        -1/24 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
    """
    return _compute_q_bracket(polynomial).to_polynomial(names) if output_as_polynomial else _compute_q_bracket(polynomial)
