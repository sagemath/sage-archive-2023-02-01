"""
Hypergeometric motives

EXAMPLES::

    sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
    sage: H = Hyp(cyclotomic=([30], [1,2,3,5]))
    sage: H.alpha_beta()
    ([1/30, 7/30, 11/30, 13/30, 17/30, 19/30, 23/30, 29/30],
    [0, 1/5, 1/3, 2/5, 1/2, 3/5, 2/3, 4/5])
    sage: H.M_value() == 30**30 / (15**15 * 10**10 * 6**6)
    True
"""
from collections import defaultdict

from sage.arith.misc import divisors, gcd, euler_phi, moebius, is_prime
from sage.arith.misc import gauss_sum
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.functions.generalized import sgn
from sage.functions.other import floor
from sage.misc.functional import cyclotomic_polynomial
from sage.misc.misc_c import prod
from sage.rings.fraction_field import FractionField
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Zp
from sage.rings.padics.misc import gauss_sum as padic_gauss_sum
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.schemes.generic.spec import Spec
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField

def characteristic_polynomial_from_traces(traces, d, q, i):
    r"""
    Given a sequence of traces `t_1, \dots, t_k`, return the
    corresponding characteristic polynomial with Weil numbers as roots.

    The characteristic polynomial is defined by the generating series

    .. MATH::

        P(T) = \exp\left(- \sum_{k\geq 1} t_k \frac{T^k}{k}\right)

    and should have the property that reciprocals of all roots have
    absolute value `q^{i/2}`.

    There can be two possible signs for the leading coefficient.  If
    the degree ``d`` is even, one need at least ``floor(d/2)`` traces
    and at least one more for odd ``d``. In case the correct sign for
    the leading coefficient cannot be guessed from the given traces,
    the output is a pair (f, dictionary). The number `f` is such that the
    trace for `p^f` allows to fix the ambiguity. The dictionary maps
    the two possible traces for `p^f` to the corresponding Euler
    factors.

    INPUT:

    - ``traces`` -- a list of integers `t_1, \dots, t_k`

    - ``d`` -- the degree of the characteristic polynomial

    - ``q`` -- power of a prime number

    - ``i`` -- integer, the weight in the motivic sense

    OUTPUT:

    a polynomial or a pair (integer, dictionary)

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import characteristic_polynomial_from_traces
        sage: characteristic_polynomial_from_traces([1, 1], 1, 3, 0)
        -T + 1
        sage: characteristic_polynomial_from_traces([1], 1, 3, 0)
        -T + 1

        sage: characteristic_polynomial_from_traces([25,625], 1, 5, 4)
        -25*T + 1
        sage: characteristic_polynomial_from_traces([25], 1, 5, 4)
        -25*T + 1

        sage: characteristic_polynomial_from_traces([3,-1,-18,-49], 2, 5, 1)
        5*T^2 - 3*T + 1
        sage: characteristic_polynomial_from_traces([3], 2, 5, 1)
        5*T^2 - 3*T + 1

        sage: characteristic_polynomial_from_traces([-4,276], 4, 5, 3)
        15625*T^4 + 500*T^3 - 130*T^2 + 4*T + 1
        sage: characteristic_polynomial_from_traces([4,-276], 4, 5, 3)
        15625*T^4 - 500*T^3 + 146*T^2 - 4*T + 1

        sage: characteristic_polynomial_from_traces([1,-13,-20,71], 2, 7, 1)
        7*T^2 - T + 1
        sage: characteristic_polynomial_from_traces([1], 2, 7, 1)
        7*T^2 - T + 1

        sage: characteristic_polynomial_from_traces([36,7620], 4, 17, 3)
        24137569*T^4 - 176868*T^3 - 3162*T^2 - 36*T + 1

    TESTS::

        sage: characteristic_polynomial_from_traces([0],2,17,1)
        (2, {-34: 17*T^2 + 1, 34: -17*T^2 + 1})

        sage: characteristic_polynomial_from_traces([-36], 4, 17, 3)
        Traceback (most recent call last):
        ...
        ValueError: not enough traces were given
    """
    if len(traces) < d // 2 + d % 2:
        raise ValueError('not enough traces were given')
    t = PowerSeriesRing(QQ, 't').gen()
    ring = PolynomialRing(ZZ, 'T')

    series = sum(- api * t**(i + 1) / (i + 1) for i, api in enumerate(traces))
    N = min(len(traces), d)  # never need more than d traces
    series = series.O(N + 1).exp()
    coeffs = list(series)

    rev_coeffs = {d - k: coeffs[k] * q**(-k * i + d * i // 2)
                  for k in range(len(coeffs)) if k <= d}
    intersection = [k for k in range(len(coeffs))
                    if k in rev_coeffs and coeffs[k]]

    def poly(sign):
        data = [0 for _ in range(d + 1)]
        for k in range(len(coeffs)):
            data[k] = coeffs[k]
        for k in rev_coeffs:
            data[k] = sign * rev_coeffs[k]
        return ring(data)

    if intersection:
        idx = intersection[0]
        sign = 1 if coeffs[idx] == rev_coeffs[idx] else -1
        return poly(sign)
    else:
        p1 = poly(1)
        p2 = poly(-1)
        s1 = (p1(t).O(t ** (d + 1))).log()
        s2 = (p2(t).O(t ** (d + 1))).log()
        index = (s1 - s2).valuation()
        return (index, {-s1[index] * index: p1, -s2[index] * index: p2})
        # IDEA: also return the first index of trace
        # that would allow sign recognition
        # this means finding the first coefficient that differs in s1 and s2


def possible_hypergeometric_data(d, weight=None):
    """
    Return the list of possible parameters of hypergeometric motives.

    INPUT:

    - ``d`` -- the degree

    - ``weight`` -- optional integer, to specify the motivic weight

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import possible_hypergeometric_data
        sage: P = possible_hypergeometric_data
        sage: [len(P(i,weight=2)) for i in range(1, 7)]
        [0, 0, 10, 30, 93, 234]
    """
    bound = 2 * d * d  # to make sure that phi(n) <= d
    possible = [(i, euler_phi(i)) for i in range(1, bound + 1)
                if euler_phi(i) <= d]
    poids = [z[1] for z in possible]
    N = len(poids)
    vectors = list(WeightedIntegerVectors(d, poids))

    def iterator():
        for i, u in enumerate(vectors):
            supp_u = [j for j in range(N) if u[j]]
            for v in vectors[i + 1:]:
                if not any(v[j] for j in supp_u):
                    yield (u, v)

    def formule(u):
        return [possible[j][0] for j in range(N) for _ in range(u[j])]

    data = [HypergeometricMotive(cyclotomic=(formule(a), formule(b)))
            for a, b in iterator()]
    if weight is None:
        return data
    else:
        return [H for H in data if H.weight() == weight]


def cyclotomic_to_alpha(cyclo):
    """
    Convert a list of indices of cyclotomic polynomials
    to a list of rational numbers.

    The input represents a product of cyclotomic polynomials.

    The output is the list of arguments of the roots of the
    given product of cyclotomic polynomials.

    This is the inverse of :func:`alpha_to_cyclotomic`.

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import cyclotomic_to_alpha
        sage: cyclotomic_to_alpha([1])
        [0]
        sage: cyclotomic_to_alpha([2])
        [1/2]
        sage: cyclotomic_to_alpha([5])
        [1/5, 2/5, 3/5, 4/5]
        sage: cyclotomic_to_alpha([1,2,3,6])
        [0, 1/6, 1/3, 1/2, 2/3, 5/6]
        sage: cyclotomic_to_alpha([2,3])
        [1/3, 1/2, 2/3]
    """
    alpha = []
    for d in cyclo:
        if d == 1:
            alpha.append(QQ.zero())
        else:
            for k in range(1, d):
                if gcd(k, d) == 1:
                    alpha.append(QQ((k, d)))
    return sorted(alpha)


def alpha_to_cyclotomic(alpha):
    """
    Convert from a list of rationals arguments to a list of integers.

    The input represents arguments of some roots of unity.

    The output represent a product of cyclotomic polynomials with exactly
    the given roots.

    This is the inverse of :func:`cyclotomic_to_alpha`.

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import alpha_to_cyclotomic
        sage: alpha_to_cyclotomic([0])
        [1]
        sage: alpha_to_cyclotomic([1/2])
        [2]
        sage: alpha_to_cyclotomic([1/5,2/5,3/5,4/5])
        [5]
        sage: alpha_to_cyclotomic([1/6, 1/3, 1/2, 2/3, 5/6, 1])
        [1, 2, 3, 6]
        sage: alpha_to_cyclotomic([1/3,2/3,1/2])
        [2, 3]
    """
    cyclo = []
    Alpha = list(alpha)
    while Alpha:
        q = QQ(Alpha.pop())
        d = q.denominator()
        for k in range(1, d):
            if gcd(k, d) == 1 and QQ((k, d)) != q:
                Alpha.remove(QQ((k, d)))
        cyclo.append(d)
    return sorted(cyclo)


def capital_M(n):
    """
    Auxiliary function, used to describe the canonical scheme.

    INPUT:

    - ``n`` -- an integer

    OUTPUT:

    a rational

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import capital_M
        sage: [capital_M(i) for i in range(1,8)]
        [1, 4, 27, 64, 3125, 432, 823543]
    """
    n = ZZ(n)
    return QQ.prod(d ** (d * moebius(n / d)) for d in divisors(n))


def cyclotomic_to_gamma(cyclo_up, cyclo_down):
    """
    Convert a quotient of products of cyclotomic polynomials
    to a quotient of products of polynomials `x^n - 1`.

    INPUT:

    - ``cyclo_up`` -- list of indices of cyclotomic polynomials
    - ``cyclo_down`` -- list of indices of cyclotomic polynomials

    OUTPUT:

    a dictionary mapping an integer `n` to the power of `x^n - 1` that
    appears in the given product

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import cyclotomic_to_gamma
        sage: cyclotomic_to_gamma([6], [1])
        {2: -1, 3: -1, 6: 1}
    """
    dico = defaultdict(lambda: 0)
    for d in cyclo_up:
        dico[d] += 1
    for d in cyclo_down:
        dico[d] -= 1

    resu = defaultdict(lambda: 0)
    for n in dico:
        for d in divisors(n):
            resu[d] += moebius(n / d) * dico[n]

    return {d: resu[d] for d in resu if resu[d]}


def gamma_list_to_cyclotomic(galist):
    r"""
    Convert a quotient of products of polynomials `x^n - 1`
    to a quotient of products of cyclotomic polynomials.

    INPUT:

    - ``galist`` -- a list of integers, where an integer `n` represents
      the power `(x^{|n|} - 1)^{\operatorname{sgn}(n)}`

    OUTPUT:

    a pair of list of integers, where `k` represents the cyclotomic
    polynomial `\Phi_k`

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import gamma_list_to_cyclotomic
        sage: gamma_list_to_cyclotomic([-1, -1, 2])
        ([2], [1])

        sage: gamma_list_to_cyclotomic([-1, -1, -1, -3, 6])
        ([2, 6], [1, 1, 1])

        sage: gamma_list_to_cyclotomic([-1, 2, 3, -4])
        ([3], [4])

        sage: gamma_list_to_cyclotomic([8,2,2,2,-6,-4,-3,-1])
        ([2, 2, 8], [3, 3, 6])
    """
    resu = defaultdict(lambda: 0)
    for n in galist:
        eps = sgn(n)
        for d in divisors(abs(n)):
            resu[d] += eps

    return (sorted(d for d in resu for k in range(resu[d])),
            sorted(d for d in resu for k in range(-resu[d])))


class HypergeometricMotive(object):
    def __init__(self, cyclotomic=None, alpha_beta=None, gamma_list=None):
        """
        Creation of hypergeometric motives.

        INPUT:

        three possibilities are offered:

        - ``cyclotomic`` -- a pair of lists of nonnegative integers,
          each integer `k` represents a cyclotomic polynomial `\Phi_k`

        - ``alpha_beta`` -- a pair of lists of rationals,
          each rational represents a root of unity

        - ``gamma_list`` -- a pair of list of nonnegative integers,
          each integer `n` represents a polynomial `x^n - 1`

        In the last case, it is also allowed to send just one list of signed
        integers where signs indicate to which part the integer belongs to.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(cyclotomic=([2],[1]))
            Hypergeometric motive for [1/2] and [0]

            sage: Hyp(alpha_beta=([1/2],[0]))
            Hypergeometric motive for [1/2] and [0]
            sage: Hyp(alpha_beta=([1/5,2/5,3/5,4/5],[0,0,0,0]))
            Hypergeometric motive for [1/5, 2/5, 3/5, 4/5] and [0, 0, 0, 0]

            sage: Hyp(gamma_list=([5],[1,1,1,1,1]))
            Hypergeometric motive for [1/5, 2/5, 3/5, 4/5] and [0, 0, 0, 0]
        """
        if gamma_list is not None:
            if isinstance(gamma_list[0], (list, tuple)):
                pos, neg = gamma_list
                gamma_list = pos + [-u for u in neg]
            cyclotomic = gamma_list_to_cyclotomic(gamma_list)
        if cyclotomic is not None:
            cyclo_up, cyclo_down = cyclotomic
            if any(x in cyclo_up for x in cyclo_down):
                raise ValueError('must be prime')
            deg = sum(euler_phi(x) for x in cyclo_down)
            up_deg = sum(euler_phi(x) for x in cyclo_up)
            if up_deg != deg:
                msg = 'not the same degree: {} != {}'.format(up_deg, deg)
                raise ValueError(msg)
            cyclo_up.sort()
            cyclo_down.sort()
            alpha = cyclotomic_to_alpha(cyclo_up)
            beta = cyclotomic_to_alpha(cyclo_down)
        elif alpha_beta is not None:
            alpha, beta = alpha_beta
            if len(alpha) != len(beta):
                raise ValueError('alpha and beta not of the same length')
            alpha = sorted(u - floor(u) for u in alpha)
            beta = sorted(u - floor(u) for u in beta)
            cyclo_up = alpha_to_cyclotomic(alpha)
            cyclo_down = alpha_to_cyclotomic(beta)
            deg = sum(euler_phi(x) for x in cyclo_down)

        self._cyclo_up = cyclo_up
        self._cyclo_down = cyclo_down
        self._alpha = alpha
        self._beta = beta
        self._deg = deg

    def __repr__(self):
        """
        Return the string representation.

        This displays the rational arguments of the roots of unity.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/2],[0]))
            Hypergeometric motive for [1/2] and [0]
        """
        txt = "Hypergeometric motive for {} and {}"
        return txt.format(self._alpha, self._beta)

    def twist(self):
        r"""
        Return the twist of ``self``.

        This is defined by adding `1/2` to each rational in `\alpha`
        and `\beta`.

        This is an involution.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(alpha_beta=([1/2],[0]))
            sage: H.twist()
            Hypergeometric motive for [0] and [1/2]

            sage: Hyp(cyclotomic=([6],[1,2])).twist().cyclotomic_data()
            ([3], [1, 2])
        """
        alpha = [x + QQ((1, 2)) for x in self._alpha]
        beta = [x + QQ((1, 2)) for x in self._beta]
        return HypergeometricMotive(alpha_beta=(alpha, beta))

    def swap_alpha_beta(self):
        """
        Return the hypergeometric motive with ``alpha`` and ``beta`` exchanged.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(alpha_beta=([1/2],[0]))
            sage: H.swap_alpha_beta()
            Hypergeometric motive for [0] and [1/2]
        """
        alpha, beta = self.alpha_beta()
        return HypergeometricMotive(alpha_beta=(beta, alpha))

    def primitive_data(self):
        """
        Return a primitive version.

        .. SEEALSO::

            :meth:`is_primitive`, :meth:`primitive_index`,

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(cyclotomic=([3],[4]))
            sage: H2 = Hyp(gamma_list=[-2, 4, 6, -8])
            sage: H2.primitive_data() == H
            True
        """
        g = self.gamma_list()
        d = gcd(g)
        return HypergeometricMotive(gamma_list=[x / d for x in g])

    def weight(self):
        """
        Return the motivic weight of ``self``.

        EXAMPLES:

        With rational inputs::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).weight()
            0
            sage: Hyp(alpha_beta=([1/4,3/4],[0,0])).weight()
            1
            sage: Hyp(alpha_beta=([1/6,1/3,2/3,5/6],[0,0,1/4,3/4])).weight()
            1
            sage: H = Hyp(alpha_beta=([1/6,1/3,2/3,5/6],[1/8,3/8,5/8,7/8]))
            sage: H.weight()
            1

        With cyclotomic inputs::

            sage: Hyp(cyclotomic=([6,2],[1,1,1])).weight()
            2
            sage: Hyp(cyclotomic=([6],[1,2])).weight()
            0
            sage: Hyp(cyclotomic=([8],[1,2,3])).weight()
            0
            sage: Hyp(cyclotomic=([5],[1,1,1,1])).weight()
            3
            sage: Hyp(cyclotomic=([5,6],[1,1,2,2,3])).weight()
            1
            sage: Hyp(cyclotomic=([3,8],[1,1,1,2,6])).weight()
            2
            sage: Hyp(cyclotomic=([3,3],[2,2,4])).weight()
            1

        With gamma list input::

            sage: Hyp(gamma_list=([8,2,2,2],[6,4,3,1])).weight()
            3
        """
        alpha = self._alpha
        beta = self._beta
        D = [sum(1 for a in alpha if a <= x) -
             sum(1 for b in beta if b <= x)
             for x in alpha + beta]
        return ZZ(max(D) - min(D) - 1)

    def degree(self):
        """
        Return the degree.

        This is the sum of the Hodge numbers.

        .. SEEALSO::

            :meth:`hodge_numbers`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).degree()
            1
            sage: Hyp(gamma_list=([2,2,4],[8])).degree()
            4
            sage: Hyp(cyclotomic=([5,6],[1,1,2,2,3])).degree()
            6
            sage: Hyp(cyclotomic=([3,8],[1,1,1,2,6])).degree()
            6
            sage: Hyp(cyclotomic=([3,3],[2,2,4])).degree()
            4
        """
        return self._deg

    def defining_polynomials(self):
        """
        Return the pair of products of cyclotomic polynomials.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/4,3/4],[0,0])).defining_polynomials()
            (x^2 + 1, x^2 - 2*x + 1)
        """
        up = prod(cyclotomic_polynomial(d) for d in self._cyclo_up)
        down = prod(cyclotomic_polynomial(d) for d in self._cyclo_down)
        return (up, down)

    def cyclotomic_data(self):
        """
        Return the pair of lists of indices of cyclotomic polynomials.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).cyclotomic_data()
            ([2], [1])
        """
        return (self._cyclo_up, self._cyclo_down)

    def alpha_beta(self):
        """
        Return the pair of lists of rational arguments.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).alpha_beta()
            ([1/2], [0])
        """
        return (self._alpha, self._beta)

    def M_value(self):
        """
        Return the `M` coefficient that appears in the equations.

        OUTPUT:

        a rational

        .. SEEALSO:: :meth:`canonical_scheme`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(alpha_beta=([1/6,1/3,2/3,5/6],[1/8,3/8,5/8,7/8]))
            sage: H.M_value()
            729/4096
            sage: Hyp(alpha_beta=(([1/2,1/2,1/2,1/2],[0,0,0,0]))).M_value()
            256
            sage: Hyp(cyclotomic=([5],[1,1,1,1])).M_value()
            3125
        """
        up = QQ.prod(capital_M(d) for d in self._cyclo_up)
        down = QQ.prod(capital_M(d) for d in self._cyclo_down)
        return up / down

    def gamma_array(self):
        r"""
        Return the dictionary `\{v: \gamma_v\}` for the expression

        .. MATH::

            \prod_v (T^v - 1)^{\gamma_v}

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).gamma_array()
            {1: -2, 2: 1}
            sage: Hyp(cyclotomic=([6,2],[1,1,1])).gamma_array()
            {1: -3, 3: -1, 6: 1}
        """
        return cyclotomic_to_gamma(self._cyclo_up, self._cyclo_down)

    def gamma_list(self):
        r"""
        Return a list of integers describing the `x^n - 1` factors.

        Each integer `n` stands for `(x^{|n|} - 1)^{\operatorname{sgn}(n)}`.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).gamma_list()
            [-1, -1, 2]

            sage: Hyp(cyclotomic=([6,2],[1,1,1])).gamma_list()
            [-1, -1, -1, -3, 6]

            sage: Hyp(cyclotomic=([3],[4])).gamma_list()
            [-1, 2, 3, -4]
        """
        gamma = self.gamma_array()
        resu = []
        for v in gamma:
            resu += [sgn(gamma[v]) * v] * abs(gamma[v])
        return resu

    def __eq__(self, other):
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H1 = Hyp(alpha_beta=([1/2],[0]))
            sage: H2 = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: H1 == H1
            True
            sage: H1 == H2
            False
        """
        return (self._alpha == other._alpha and
                self._beta == other._beta)

    def __ne__(self, other):
        """
        Return whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H1 = Hyp(alpha_beta=([1/2],[0]))
            sage: H2 = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: H1 != H1
            False
            sage: H1 != H2
            True
        """
        return not (self == other)

    def is_primitive(self):
        """
        Return whether ``self`` is primitive.

        .. SEEALSO::

            :meth:`primitive_index`, :meth:`primitive_data`,

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(cyclotomic=([3],[4])).is_primitive()
            True
            sage: Hyp(gamma_list=[-2, 4, 6, -8]).is_primitive()
            False
            sage: Hyp(gamma_list=[-3, 6, 9, -12]).is_primitive()
            False
        """
        return self.primitive_index() == 1  # ?

    def primitive_index(self):
        """
        Return the primitive index.

        .. SEEALSO::

            :meth:`is_primitive`, :meth:`primitive_data`,

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(cyclotomic=([3],[4])).primitive_index()
            1
            sage: Hyp(gamma_list=[-2, 4, 6, -8]).primitive_index()
            2
            sage: Hyp(gamma_list=[-3, 6, 9, -12]).primitive_index()
            3
        """
        return gcd(self.gamma_list())

    def has_symmetry_at_one(self):
        """
        If ``True``, the motive H(t=1) is a direct sum of two motives.

        Note that simultaneous exchange of (t,1/t) and (alpha,beta)
        always gives the same motive.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: Hyp(alpha_beta=[[1/2]*16,[0]*16]).has_symmetry_at_one()
            True

        REFERENCE:

        - https://www.matrix-inst.org.au/wp_Matrix2016/wp-content/uploads/2016/04/Roberts-2.pdf
        """
        _, beta_twist = self.twist().alpha_beta()
        return self.degree() % 2 == 0 and self._alpha == beta_twist

    def hodge_numbers(self):
        """
        Return the Hodge numbers.

        .. SEEALSO::

            :meth:`degree`, :meth:`hodge_polynomial`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(cyclotomic=([3],[6]))
            sage: H.hodge_numbers()
            [1, 1]

            sage: H = Hyp(cyclotomic=([4],[1,2]))
            sage: H.hodge_numbers()
            [2]

            sage: H = Hyp(gamma_list=([8,2,2,2],[6,4,3,1]))
            sage: H.hodge_numbers()
            [1, 2, 2, 1]

            sage: H = Hyp(gamma_list=([5],[1,1,1,1,1]))
            sage: H.hodge_numbers()
            [1, 1, 1, 1]

            sage: H = Hyp(gamma_list=[6,1,-4,-3])
            sage: H.hodge_numbers()
            [1, 1]

            sage: H = Hyp(gamma_list=[-3]*4 + [1]*12)
            sage: H.hodge_numbers()
            [1, 1, 1, 1, 1, 1, 1, 1]
        """
        alpha = [(x, 'a') for x in self._alpha]
        beta = [(x, 'b') for x in self._beta]
        height = 0
        hodge = defaultdict(lambda: 0)
        for x, letter in sorted(alpha + beta):
            if letter == 'a':
                hodge[height] += 1
                height += 1
            else:
                height -= 1
        return [hodge[i] for i in sorted(hodge)]

    def hodge_polynomial(self):
        """
        Return the Hodge polynomial.

        .. SEEALSO::

            :meth:`hodge_numbers`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(cyclotomic=([6,10],[3,12]))
            sage: H.hodge_polynomial()
            (T^3 + 2*T^2 + 2*T + 1)/T^2
            sage: H = Hyp(cyclotomic=([2,2,2,2,3,3,3,6,6],[1,1,4,5,9]))
            sage: H.hodge_polynomial()
            (T^5 + 3*T^4 + 3*T^3 + 3*T^2 + 3*T + 1)/T^2
        """
        alpha = self._alpha
        beta = self._beta

        def D(x):
            return (sum(1 for a in alpha if a <= x) -
                    sum(1 for b in beta if 1 - b <= x))

        def z(x):
            return alpha.count(x)

        T = polygen(ZZ, 'T')
        return sum(T ** (D(a) - z(a)) * (T**z(a) - 1) // (T - 1)
                   for a in set(alpha))

    def padic_H_value(self, p, f, t, prec=20):
        """
        Return the `p`-adic trace of the Frobenius.

        INPUT:

        - `p` -- a prime number

        - `f` -- an integer such that `q = p^f`

        - `t` -- a rational parameter

        - ``prec`` -- precision (optional, default 20)

        OUTPUT:

        an integer

        .. WARNING::

            This is not yet working correctly.

        This is denoted by `U_q(t)` in the reference below.

        EXAMPLES:

        From Benasque report, page 8::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(alpha_beta=([1/2]*4,[0]*4))
            sage: [H.padic_H_value(3,i,-1) for i in range(1,3)]
            [0, -12]
            sage: [H.padic_H_value(5,i,-1) for i in range(1,3)]
            [-4, 276]
            sage: [H.padic_H_value(7,i,-1) for i in range(1,3)]
            [0, -476]
            sage: [H.padic_H_value(11,i,-1) for i in range(1,3)]
            [0, -4972]

        From slides::

            sage: H = Hyp(gamma_list=[-6,-1,4,3])
            sage: t = 189/125
            sage: H.padic_H_value(13,1,t)
            0

        REFERENCE:

        - http://magma.maths.usyd.edu.au/~watkins/papers/HGM-chapter.pdf
        """
        beta = self._beta
        t = QQ(t)
        gamma = self.gamma_array()
        q = p ** f

        m = {r: beta.count(QQ((r, q - 1))) for r in range(q - 1)}
        M = self.M_value()
        D = (self.weight() + 1 - m[0]) // 2

        gauss_table = [padic_gauss_sum(r, p, f, prec, factored=True) for r in range(q - 1)]

        p_ring = Zp(p, prec=prec)
        teich = p_ring.teichmuller(t * M)
        sigma = sum(q**(D + m[0] - m[r]) *
                    (-p)**(sum(gauss_table[(v * r) % (q - 1)][0] * gv
                             for v, gv in gamma.items())//(p-1)) *
                    prod(gauss_table[(v * r) % (q - 1)][1] ** gv
                         for v, gv in gamma.items()) *
                    teich ** r
                    for r in range(q - 1))
        resu = ZZ(-1) ** m[0] / (1 - q) * sigma
        return IntegerModRing(p**prec)(resu).lift_centered()

    def H_value(self, p, f, t, ring=None):
        """
        Return the trace of the Frobenius.

        INPUT:

        - `p` -- a prime number

        - `f` -- an integer such that `q = p^f`

        - `t` -- a rational parameter

        - ``ring`` -- optional (default ``UniversalCyclotomicfield``)

        The ring could be also ``ComplexField(n)`` or ``QQbar``.

        OUTPUT:

        an integer

        .. WARNING::

            This is apparently working correctly as can be tested
            using ComplexField(70) as value ring.

            Using instead UniversalCyclotomicfield, this is much
            slower than the `p`-adic version :meth:`padic_H_value`.

        EXAMPLES:

        With values in the UniversalCyclotomicField (slow)::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(alpha_beta=([1/2]*4,[0]*4))
            sage: [H.H_value(3,i,-1) for i in range(1,3)]
            [0, -12]
            sage: [H.H_value(5,i,-1) for i in range(1,3)]
            [-4, 276]
            sage: [H.H_value(7,i,-1) for i in range(1,3)]  # not tested
            [0, -476]
            sage: [H.H_value(11,i,-1) for i in range(1,3)]  # not tested
            [0, -4972]
            sage: [H.H_value(13,i,-1) for i in range(1,3)]  # not tested
            [-84, -1420]

        With values in ComplexField::

            sage: [H.H_value(5,i,-1, ComplexField(60)) for i in range(1,3)]
            [-4, 276]

        REFERENCES:

        - https://arxiv.org/pdf/1505.02900.pdf, Theorem 1.3

        - http://users.ictp.it/~villegas/hgm/benasque-2009-report.pdf
        """
        if ring is None:
            ring = UniversalCyclotomicField()
        beta = self._beta
        t = QQ(t)
        gamma = self.gamma_array()
        q = p ** f

        m = {r: beta.count(QQ((r, q - 1))) for r in range(q - 1)}
        D = (self.weight() + 1 - m[0]) // 2
        M = self.M_value()

        Fq = GF(q)
        gen = Fq.multiplicative_generator()
        zeta_q = ring.zeta(q - 1)

        tM = Fq(t * M)
        for k in range(q - 1):
            if gen ** k == tM:
                teich = zeta_q ** k
                break

        gauss_table = [gauss_sum(zeta_q ** r, Fq) for r in range(q - 1)]

        sigma = sum(q**(D + m[0] - m[r]) *
                    prod(gauss_table[(-v * r) % (q - 1)] ** gv
                         for v, gv in gamma.items()) *
                    teich ** r
                    for r in range(q - 1))
        resu = ZZ(-1) ** m[0] / (1 - q) * sigma
        if not ring.is_exact():
            resu = resu.real_part().round()
        return resu

    def euler_factor(self, t, p, degree=0):
        """
        Return the Euler factor of the motive `H_t` at prime `p`.

        INPUT:

        - `t` -- rational number, not 0 or 1

        - `p` -- prime number

        - ``degree`` -- optional integer (default 0)

        OUTPUT:

        a polynomial

        See http://users.ictp.it/~villegas/hgm/benasque-2009-report.pdf
        for explicit examples of Euler factors.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(alpha_beta=([1/2]*4,[0]*4))
            sage: H.euler_factor(-1, 5)
            15625*T^4 + 500*T^3 - 130*T^2 + 4*T + 1

            sage: [Hyp(cyclotomic=([6,2],[1,1,1])).euler_factor(4,p)
            ....:  for p in [5,7,11,13,17,19]]
            [125*T^3 + 20*T^2 + 4*T + 1,
             343*T^3 - 42*T^2 - 6*T + 1,
             -1331*T^3 - 22*T^2 + 2*T + 1,
             -2197*T^3 - 156*T^2 + 12*T + 1,
             4913*T^3 + 323*T^2 + 19*T + 1,
             6859*T^3 - 57*T^2 - 3*T + 1]

            sage: H = Hyp(gamma_list=[-6,-1,4,3])
            sage: t = 189/125
            sage: H.euler_factor(t,11)
            11*T^2 + 4*T + 1
            sage: H.euler_factor(t,13)
            13*T^2 + 1
            sage: H.euler_factor(t,17)
            17*T^2 + 1
            sage: H.euler_factor(t,19)
            19*T^2 + 1
            sage: H.euler_factor(t,23)
            23*T^2 + 8*T + 1
            sage: H.euler_factor(t,29)
            29*T^2 + 2*T + 1

        REFERENCE:

        - https://icerm.brown.edu/materials/Slides/sp-f15-offweeks/Hypergeomteric_Motives,_I_]_David_Roberts,_University_of_Minnesota_-_Morris.pdf
        """
        if t not in QQ or t in [0, 1]:
            raise ValueError('wrong t')
        if not is_prime(p):
            raise ValueError('p not prime')
        if not all(x.denominator() % p for x in self._alpha + self._beta):
            raise ValueError('p is wild')
        if (t.valuation(p) or (t - 1).valuation(p) > 0):
            raise ValueError('p is tame')
        # now p is good
        if degree == 0:
            d = self.degree()
        bound = d // 2 + d % 2
        traces = [self.padic_H_value(p, i + 1, t) for i in range(bound)]

        w = self.weight()

        # One has to handle the cases of sign ambiguity
        resu = characteristic_polynomial_from_traces(traces, d, p, w)
        if isinstance(resu, tuple):
            f, dico = resu
            return dico[self.padic_H_value(p, f, t)]
        else:
            return resu

    def canonical_scheme(self, t=None):
        """
        Return the canonical scheme.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricMotive as Hyp
            sage: H = Hyp(cyclotomic=([3],[4]))
            sage: H.gamma_list()
            [-1, 2, 3, -4]
            sage: H.canonical_scheme()
            Spectrum of Quotient of Multivariate Polynomial Ring
            in X0, X1, Y0, Y1 over Fraction Field of Univariate Polynomial Ring
            in t over Rational Field by the ideal
            (X0 + X1 - 1, Y0 + Y1 - 1, -X0^2*X1^3 + 27/64*t*Y0*Y1^4)

            sage: H = Hyp(gamma_list=[-2, 3, 4, -5])
            sage: H.canonical_scheme()
            Spectrum of Quotient of Multivariate Polynomial Ring
            in X0, X1, Y0, Y1 over Fraction Field of Univariate Polynomial Ring
            in t over Rational Field by the ideal
            (X0 + X1 - 1, Y0 + Y1 - 1, -X0^3*X1^4 + 1728/3125*t*Y0^2*Y1^5)
        """
        if t is None:
            t = FractionField(QQ['t']).gen()
        basering = t.parent()
        gamma_pos = [u for u in self.gamma_list() if u > 0]
        gamma_neg = [u for u in self.gamma_list() if u < 0]
        N_pos = len(gamma_pos)
        N_neg = len(gamma_neg)
        varX = ['X{}'.format(i) for i in range(N_pos)]
        varY = ['Y{}'.format(i) for i in range(N_neg)]
        ring = PolynomialRing(basering, varX + varY)
        gens = ring.gens()
        X = gens[:N_pos]
        Y = gens[N_pos:]
        eq0 = ring.sum(X) - 1
        eq1 = ring.sum(Y) - 1
        eq2_pos = ring.prod(X[i] ** gamma_pos[i] for i in range(N_pos))
        eq2_neg = ring.prod(Y[j] ** -gamma_neg[j] for j in range(N_neg))

        ideal = ring.ideal([eq0, eq1, self.M_value() * t * eq2_neg - eq2_pos])
        return Spec(ring.quotient(ideal))
