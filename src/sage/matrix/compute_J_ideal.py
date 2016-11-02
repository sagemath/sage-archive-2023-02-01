"""
Calculate nu
"""

from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject
from time import time


def lifting(p, t, A, G):
    r"""
    Compute generators of `\{f \in D[X]^d \mid Af \equiv 0 \pmod{p^{t}}\}` given
    generators of `\{f\in D[X]^d \mid Af \equiv 0\pmod{p^{t-1}}\}`.

    INPUT:

    - ``p`` -- a prime element of some ring `D`

    - ``t`` -- a non-negative integer

    - ``A`` -- a `c\times d` matrix  over `D[X]`

    - ``G`` -- a matrix over `D[X]`. The columns of `(p^{t-1}I G)` are generators of
      `\{ f\in D[X]^d \mid Af \equiv 0\pmod{p^{t-1}}\}`. Can be set to ``None`` if ``t`` is zero.

    OUTPUT:

    A matrix ``F`` over `D[X]` such that the columns of `(p^tI F pG)` are generators of
    `\{ f\in D[X]^d \mid Af \equiv 0\pmod{p^{t-1}}\}`.

    EXAMPLES::

        sage: X = polygen(ZZ, 'X')
        sage: A = matrix([[1, X], [2*X, X^2]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        []
        sage: assert (A*G1 % 5).is_zero()
        sage: A = matrix([[1, X, X^2], [2*X, X^2, 3*X^3]])
        sage: G0 = lifting(5, 0, A, None)
        sage: G1 = lifting(5, 1, A, G0); G1
        [3*X^2]
        [    X]
        [    1]
        sage: assert (A*G1 % 5).is_zero()
        sage: G2 = lifting(5, 2, A, G1); G2
        [15*X^2 23*X^2]
        [   5*X      X]
        [     5      1]
        sage: assert (A*G2 % 25).is_zero()
        sage: lifting(5, 10, A, G1)
        Traceback (most recent call last):
        ...
        AssertionError: A*G is not zero mod 5^9
    """
    if t == 0:
        return matrix(A.parent().base(), A.ncols(), 0)

    assert (A*G % p**(t-1)).is_zero(),\
        "A*G is not zero mod %s^%s" % (str(p), str(t-1))

    P = A.parent()
    ZZX = P.base()
    (X,) = ZZX.gens()
    d = A.ncols()

    R = A*G/p**(t-1)
    R.change_ring(ZZX)

    AR = matrix.block([[A, R]])
    Fp = GF(p)
    FpX = PolynomialRing(Fp, name=X)

    ARb = AR.change_ring(FpX)
    #starting = time()
    (Db, Sb, Tb) = ARb.smith_form()
    #print "(SNF: %s sec)" % str(time()-starting)
    assert Sb * ARb * Tb == Db
    assert all(i == j or Db[i, j].is_zero()
               for i in range(Db.nrows())
               for j in range(Db.ncols()))

    r = Db.rank()
    T = Tb.change_ring(ZZX)

    F1 = matrix.block([[p**(t-1) * matrix.identity(d), G]])*T
    F = F1.matrix_from_columns(range(r, F1.ncols()))
    assert (A*F % (p**t)).is_zero(), "A*F=%s" % str(A*F)

    return F


def p_part(f, p):
    r"""
    Compute the `p`-part of a polynomial.

    INPUT:

    - ``f`` -- a polynomial over `D`

    - ``p`` -- a prime in `D`

    OUTPUT:

    A polynomial ``g`` such that `\deg g \le \deg f` and
    all non-zero coefficients of `f - p g` are not divisible by `p`.

    EXAMPLES::

        sage: X = polygen(ZZ, 'X')
        sage: f = X^3 + 5*X+25
        sage: g = p_part(f, 5); g
        X + 5
        sage: f - 5*g
        X^3
    """
    ZZX = f.parent()
    (X,) = ZZX.gens()
    p_prt = 0
    coefficients = f.list()
    for i in range(len(coefficients)):
        coeff = coefficients[i]
        if coeff%p == 0:
            p_prt = p_prt + coeff//p*X**i

    return p_prt


class ComputeMinimalPolynomials(SageObject):
    r"""
    Create an object for computing `(p^t)`-minimal polynomials and `J`-ideals.

    INPUT:

    - ``B`` -- a square matrix over a principal ideal domain `D`.

    OUTPUT:

    An object which allows to call ``p_minimal_polynomials`` and ``null_ideal``.

    The characteristic polynomial of `B` is denoted by `\chi_B`; `n` is the size of `B`.

    EXAMPLES::

        sage: from calculate_nu import ComputeMinimalPolynomials # not tested
        sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
        sage: C = ComputeMinimalPolynomials(B)
        sage: for t in range(3):
        ....:     print C.null_ideal(2**t)
        Ideal (1, x^3 + x^2 - 12*x - 20) of Univariate Polynomial Ring in x over Integer Ring
        Ideal (2, x^3 + x^2 - 12*x - 20, x^2 + x, x^2 + x) of Univariate Polynomial Ring in x over Integer Ring
        Ideal (4, x^3 + x^2 - 12*x - 20, x^2 + 3*x + 2, x^2 + 3*x + 2) of Univariate Polynomial Ring in x over Integer Ring
        sage: C.p_minimal_polynomials(2)
        [x^3 + 7*x^2 + 6*x]
        ([2], {2: x^2 + 3*x + 2})

    .. TODO:: There should not be several polynomials of the same degree.

    .. TODO:: Test composite ``b`` for ``null_ideal``

    .. TODO:: Implement and test ``b=0`` for ``null_ideal``

    """
    def __init__(self, B):
        r"""
        Initialize the ComputeMinimalPolynomials class.

        INPUT:

        - ``B`` -- a square matrix.

        TESTS::

            sage: ComputeMinimalPolynomials(matrix([[1, 2]]))
            Traceback (most recent call last):
            ...
            TypeError: square matrix required.
        """
        super(ComputeMinimalPolynomials, self).__init__()
        if not B.is_square():
            raise TypeError("square matrix required.")

        self._B = B
        X = polygen(B.base_ring())
        adjunct = (X-B).adjoint()
        d = B.nrows()**2
        b = matrix(d, 1, adjunct.list())
        chi_B = B.charpoly(X)
        self._A = matrix.block([[b , -chi_B*matrix.identity(d)]])
        self._A.set_immutable()
        self._ZX = X.parent()


    # ersetzt Polynome im p^t-Ideal durch normierte (Lemma 5.4)
    def find_monic_replacements(self, p, t, poly_set, prev_nu):
        r"""
        Replace possibly non-monic generators of `N_{p^t}(B)` by monic generators

        INPUT:

        - ``p`` -- a prime element of `D`

        - ``t`` -- a non-negative integer

        - ``poly_set`` -- a list of polynomials over `D[X]`. Together with
          `pN_{p^{t-1}}(B)`, they generate `N_{p^t}(B)`.

        - ``prev_nu`` -- a `p^{t-1}`-minimal polynomial of `B`.

        OUTPUT:

        A list of monic polynomials. Together with `pN_{p^{t-1}}(B)`,
        they generate `N_{p^t}(B)`.

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: x = polygen(ZZ, 'x')
            sage: nu_2 = x^2 + x
            sage: generators_4 = [2*x^2 + 2*x, x^2 + 3*x + 2]
            sage: C.find_monic_replacements(2, 2, generators_4, nu_2)
            [x^2 + 3*x + 2]

        TESTS::

            sage: C.find_monic_replacements(2, 3, generators_4, nu_2)
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert all((f(self._B) % p**t).is_zero()
                   for f in poly_set)

        (X,) = self._ZX.gens()

        replacements = []
        for f in poly_set:
            g = self._ZX(f)
            nu = self._ZX(prev_nu)
            p_prt = self._ZX(p_part(g, p))

            while g != p*p_prt:
                r = p_prt.quo_rem(nu)[1]
                g2 = g - p*p_prt
                d,u,v = xgcd(g2.leading_coefficient(), p)
                tmp_h = p*r + g2
                h = u*tmp_h + v*p*prev_nu*X**(tmp_h.degree()-prev_nu.degree())
                replacements.append(h % p**t)
                #reduce coefficients mod p^t to keep coefficients small
                g = g.quo_rem(h)[1]
                p_prt = self._ZX(p_part(g,p))

        replacements = list(set(replacements))
        assert all( g.is_monic() for g in replacements),\
            "Something went wrong in find_monic_replacements"
        return replacements


    # Algorithm 2
    def current_nu(self, p, t, pt_generators, prev_nu):
        r"""
        Compute `p^t`-minimal polynomial of `B`.

        INPUT:

        - ``p`` -- a prime element of `D`

        - ``t`` -- a non-negative integer

        - ``pt_generators`` -- a list of polynomials over `D[X]`. Together with
          `pN_{p^{t-1}}(B)`, they generate `N_{p^t}(B)`.

        - ``prev_nu`` -- a `p^{t-1}`-minimal polynomial of `B`.

        OUTPUT:

        A `p^t`-minimal polynomial of `B`.

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: x = polygen(ZZ, 'x')
            sage: nu_2 = x^2 + x
            sage: generators_4 = [2*x^2 + 2*x, x^2 + 3*x + 2]
            sage: C.current_nu(2, 2, generators_4, nu_2)
            x^2 + 3*x + 2

        TESTS::

            sage: C.current_nu(2, 3, generators_4, nu_2)
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert all((g(self._B) % p**t).is_zero()
                   for g in pt_generators)

        generators = self.find_monic_replacements(p, t, pt_generators, prev_nu)

        verbose("------------------------------------------")
        verbose(pt_generators)
        verbose("Generators with (p^t)-generating property:")
        verbose(generators)


        # find poly of minimal degree
        g = generators[0]
        for f in generators:
            if f.degree() < g.degree():
                g=f

        # find nu
        while len(generators) > 1:
            f = list(set(generators) - set([g]))[0]
            #take first element in generators not equal g
            generators.remove(f)
            r = (f.quo_rem(g)[1]) % p**t
            generators = generators + self.find_monic_replacements(p, t, [r], prev_nu)
            print generators

            if generators[-1].degree() < g.degree():
                g=generators[-1]

        return generators[0]


    def mccoy_column(self, p, t, nu_t):
        r"""
        INPUT:

        - ``p`` -- a prime element in `D`

        - ``t`` -- a positive integer

        - ``nu_t`` -- a `p^t`-minimal polynomial of `B`

        OUTPUT:

        An `(n^2 + 1) \times 1` matrix `g` with first entry ``nu_t`` such that
        `(b -\chi_B I)g \equiv 0\pmod{p^t}` where `b` consists
        of the entries of `\operatorname{adj}(X-B)`.

        EXAMPLES::

           sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
           sage: C = ComputeMinimalPolynomials(B)
           sage: x = polygen(ZZ, 'x')
           sage: nu_4 = x^2 + 3*x + 2
           sage: g = C.mccoy_column(2, 2, nu_4)
           sage: b = matrix(9, 1, (x-B).adjoint().list())
           sage: M = matrix.block([[b , -B.charpoly(x)*matrix.identity(9)]])
           sage: (M*g % 4).is_zero()
           True

        TESTS::

           sage: nu_2 = x^2 + x
           sage: C.mccoy_column(2, 2, nu_2)
           Traceback (most recent call last):
           ...
           AssertionError: x^2 + x not in (2^2)-ideal

        """
        assert (nu_t(self._B) % p**t).is_zero(),\
            "%s not in (%s^%s)-ideal" % (str(nu_t), str(p), str(t))
        chi_B = self._ZX(self._B.characteristic_polynomial())
        column = [nu_t]
        #print nu_t
        for b in self._A[:, 0].list():
            q,r = (nu_t*b).quo_rem(chi_B)

            column.append(q)

        assert (self._A * matrix(self._ZX,
                                 (self._A).ncols(), 1, column) % p**t).is_zero(),\
                                 "McCoy column is not correct"
        return  matrix(self._ZX, self._A.ncols(), 1, column)


    def p_minimal_polynomials(self, p, upto=None):
        r"""
        Return index set `\mathcal{S}` and monic polynomials
        `\nu_s` for `s\in \mathcal{S}` such that `N_{p^t}(B) = \mu_B
        \mathbb{Z}[X] + p^t\mathbb{Z}[X] + \sum_{s\in \mathcal{S}}
        p^{t-s}\nu_s \mathbb{Z}[X]`

        INPUT:

        - ``p`` -- an integer prime

        - ``upto`` -- a nonnegative integer Default is ``None``: Returns
          `\mathcal{S}` such that
          `N_{p^t}(B) = \mu_B \mathbb{Z}[X] + p^t\mathbb{Z}[X] + \sum_{s\in\mathcal{S}} p^{t-s}\nu_s \mathbb{Z}[X]`
          holds for all `t \ge \max\{s\in \mathcal{S}\}`.

        OUTPUT:

        A list (index set `\mathcal{S}`) together with a dictionary
        (keys=indices in `\mathcal{S}`, values=polynomials `\nu_s` )

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: C.p_minimal_polynomials(2)
            [x^3 + 7*x^2 + 6*x]
            ([2], {2: x^2 + 3*x + 2})
            sage: set_verbose(1)
            sage: C.p_minimal_polynomials(2)
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) ------------------------------------------
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) p=2, t=1:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) Result of lifting:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) F=
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) [x^2 + x]
            [      x]
            [      0]
            [      1]
            [      1]
            [  x + 1]
            [      1]
            [      0]
            [      0]
            [  x + 1]
            verbose 1 (...: calculate_nu.py, current_nu) ------------------------------------------
            verbose 1 (...: calculate_nu.py, current_nu) [x^2 + x]
            verbose 1 (...: calculate_nu.py, current_nu) Generators with (p^t)-generating property:
            verbose 1 (...: calculate_nu.py, current_nu) [x^2 + x]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) nu=x^2 + x
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) corresponding columns for G
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) [x^2 + x]
            [  x + 2]
            [      0]
            [      1]
            [      1]
            [  x - 1]
            [     -1]
            [     10]
            [      0]
            [  x + 1]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) ------------------------------------------
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) p=2, t=2:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) Result of lifting:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) F=
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) [  2*x^2 + 2*x x^2 + 3*x + 2]
            [          2*x         x + 4]
            [            0             0]
            [            2             1]
            [            2             1]
            [      2*x + 2         x + 1]
            [            2            -1]
            [            0            10]
            [            0             0]
            [      2*x + 2         x + 3]
            verbose 1 (...: calculate_nu.py, current_nu) ------------------------------------------
            verbose 1 (...: calculate_nu.py, current_nu) [2*x^2 + 2*x, x^2 + 3*x + 2]
            verbose 1 (...: calculate_nu.py, current_nu) Generators with (p^t)-generating property:
            verbose 1 (...: calculate_nu.py, current_nu) [x^2 + 3*x + 2]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) nu=x^2 + 3*x + 2
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) corresponding columns for G
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) [x^2 + 3*x + 2]
            [        x + 4]
            [            0]
            [            1]
            [            1]
            [        x + 1]
            [           -1]
            [           10]
            [            0]
            [        x + 3]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) ------------------------------------------
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) p=2, t=3:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) Result of lifting:
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) F=
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) [x^3 + 7*x^2 + 6*x x^3 + 3*x^2 + 2*x]
            [        x^2 + 8*x         x^2 + 4*x]
            [                0                 0]
            [                x             x + 4]
            [            x + 4                 x]
            [    x^2 + 5*x + 4           x^2 + x]
            [           -x + 4                -x]
            [             10*x              10*x]
            [                0                 0]
            [        x^2 + 7*x     x^2 + 3*x + 4]
            verbose 1 (...: calculate_nu.py, current_nu) ------------------------------------------
            verbose 1 (...: calculate_nu.py, current_nu) [x^3 + 7*x^2 + 6*x, x^3 + 3*x^2 + 2*x]
            verbose 1 (...: calculate_nu.py, current_nu) Generators with (p^t)-generating property:
            verbose 1 (...: calculate_nu.py, current_nu) [x^3 + 7*x^2 + 6*x, x^3 + 3*x^2 + 2*x]
            [x^3 + 7*x^2 + 6*x]
            verbose 1 (...: calculate_nu.py, p_minimal_polynomials) nu=x^3 + 7*x^2 + 6*x
            ([2], {2: x^2 + 3*x + 2})
            sage: set_verbose(0)
            sage: C.p_minimal_polynomials(2, upto=1)
            ([1], {1: x^2 + x})
        """

        mu_B = self._B.minimal_polynomial()
        deg_mu = mu_B.degree()

        t = 0
        calS = []
        p_min_polys = {}
        nu = self._ZX(1)
        d=self._A.ncols()
        G = matrix(self._ZX,d,0)


        while True:
            deg_prev_nu = nu.degree()
            t = t + 1
            verbose("------------------------------------------")
            verbose("p=%s, t=%s:" % (str(p), str(t)))

            verbose("Result of lifting:")
            verbose("F=")
            verbose(lifting(p, t, self._A, G))

            nu = self.current_nu(p,t, list(lifting(p, t, self._A, G)[0]), nu)

            verbose("nu=%s" % str(nu))
            if nu.degree() >= deg_mu:
                return calS, p_min_polys


            if nu.degree() == deg_prev_nu:
                calS.remove(t-1)
                G = G.matrix_from_columns(range(G.ncols()-1))
                del p_min_polys[t-1]

            verbose("corresponding columns for G")
            verbose(self.mccoy_column(p,t,nu))

            G = matrix.block([[p * G, self.mccoy_column(p, t, nu)]])
            calS.append(t)
            p_min_polys[t] = nu

            # allow early stopping for small t
            if t == upto:
                return calS, p_min_polys


    def null_ideal(self, b=0):
        r"""
        Return the ideal `N_{b}(B)=\{ f\in \mathbb{Z}[X] \mid \exists
        M\in\mathbb{Z}^{n\times n}\colon f \operatorname{adj}(X-B)
        \equiv \chi_B M \pmod{b}\}`.

        INPUT:

        - ``b`` -- an integer (Default value is 0)

        OUTPUT:

        An ideal in `\mathbb{Z}[X]`.

        EXAMPLES::

            sage: from calculate_nu import compute_nu # not tested
            sage: B = matrix([[1, 2], [3, 4]])
            sage: C = ComputeMinimalPolynomials(B)
            sage: for t in range(3):
            ....:     print C.null_ideal(3**t)
            Ideal (1, x^2 - 5*x - 2) of Univariate Polynomial Ring in x over Integer Ring
            Ideal (3, x^2 - 5*x - 2) of Univariate Polynomial Ring in x over Integer Ring
            Ideal (9, x^2 - 5*x - 2) of Univariate Polynomial Ring in x over Integer Ring
        """
        factorization = list(factor(b))
        generators = [self._ZX(b), self._ZX((self._B).minimal_polynomial())]
        for (p,t) in factorization:
            #print (p,t)
            cofactor = b // p**t
            calS, p_polys = self.p_minimal_polynomials(p,t)
            #print "++++"
            #print calS
            #print p_polys

            #print p_polys
            #print "Generators before: %s" % str(generators)
            for s in calS+ [t]:
                #print s
                #print p_polys[s]
                generators = generators + \
                             [self._ZX(cofactor*p**(t-s)*p_polys[s]) for s in calS]
                #print "Generators after: %s" % str(generators)


        assert all((g(self._B) % b).is_zero() for g in generators), \
            "Polynomials not in %s-ideal" % str(b)

        return self._ZX.ideal(generators)


    # which primes we need to checke
    def prime_candidates(self):
        r"""
        Determine those primes `p` where `\mu_B` might not be a `p`-minimal polynomial.

        OUTPUT:

        A list of primes.

        EXAMPLES::

             sage: from calculate_nu import compute_nu # not tested
             sage: B = matrix([[1, 2], [3, 4]])
             sage: C = ComputeMinimalPolynomials(B)
             sage: C.prime_candidates()
             [3]
             sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
             sage: C = ComputeMinimalPolynomials(B)
             sage: C.prime_candidates()
             [2, 3, 5]
        """
        F, T = (self._B).frobenius(2)
        factorization = list(factor(T.det()))

        primes = []
        for (p, t) in factorization:
            primes.append(p)

        return primes

