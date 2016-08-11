r"""
Cyclic Code

Let `F` be a field. A `[n ,k]` code `C` over `F` is called cyclic if every cyclic shift
of a codeword in `C` is also a codeword [R06]_:

    .. MATH::

        \forall c \in C,
        c = (c_{0}, c_{1}, \dots , c_{n-1}) \in C
        \Rightarrow (c_{n-1}, c_{0}, \dots , c_{n-2})

Let `c = (c_0, c_2, \dots, c_{n-1})` be a codeword of `C`.
This codeword can be seen as a polynomial over `F_n[x]`, as follows:
`\Sigma_{i=0}^{n-1} c_i \times x^i`.
There is a unique, monic polynomial `g(x)` such that for every
`c(x) \in F_n[x]`, `c(x) \in C \Leftrightarrow g(x) | c(x)`.
This polynomial is called the generator polynomial of `C`.


REFERENCES:

.. [HT] C. Hartmann and K.K. Tzeng. Generalizations of the BCH Bound.
   Information and Control, 20(5):489-498, June 1972

TESTS:

This class uses the following experimental feature:
:class:`sage.coding.relative_finite_field_extension.RelativeFiniteFieldExtension`.
This test block is here only to trigger the experimental warning so it does not
interferes with doctests::

    sage: from sage.coding.relative_finite_field_extension import *
    sage: Fqm.<aa> = GF(16)
    sage: Fq.<a> = GF(4)
    sage: RelativeFiniteFieldExtension(Fqm, Fq)
    doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
    See http://trac.sagemath.org/20284 for details.
    Relative field extension between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from linear_code import (AbstractLinearCode,
                         LinearCodeSyndromeDecoder,
                         LinearCodeNearestNeighborDecoder)
from encoder import Encoder
from decoder import Decoder
from sage.rings.integer import Integer
from sage.arith.all import gcd
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.functions.log import log
from sage.categories.homset import Hom
from copy import copy
from sage.groups.generic import discrete_log
from sage.misc.functional import multiplicative_order
from relative_finite_field_extension import RelativeFiniteFieldExtension

def find_generator_polynomial(code, probabilistic = False):
    r"""
    Returns a possible generator polynomial for ``code``.

    It works as follows: as the generator polynomial of a cyclic code `C`
    divides every word in `C`, we consider the generator matrix of `C` as
    its rows form a basis of `C`. We denote by `n` [resp: `k`] the length
    [resp: dimension] of `C`.

    Then, we take the first row as temporary generator polynomial, and
    compute a gcd between this temporary generator polynomial and the next row.
    The result of this gcd is now our temporary generator polynomial.
    We iterate over all rows, and if the the last temporary generator polynomial
    divides `x^{n} - 1` and has degree `n-k` it is the generator polynoial of `C`.

    The probabilistic version of this algorithm works as above, but instead of
    considering every row of a generator matrix we take random codewords.
    If a gcd of degree less than `n-k` is found, we take it as our generator polynomial.
    If it fails at fidning such a polynomial in less than `k` attemps, the algorithm stops.

    INPUT:

    - ``code`` -- a linear code

    - ``probabilistic`` -- (default : ``False``) enables the probabilistic version of
      the algorithm

    OUTPUT:

    - a polynomial

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(2^3, 'a').list()[1:2^3], 2^2)
        sage: sage.coding.cyclic_code.find_generator_polynomial(C)
        x^3 + (a^2 + 1)*x^2 + a*x + a^2 + 1
    """
    G = code.generator_matrix()
    n = code.length()
    k = code.dimension()
    F = code.base_ring()
    R = F['x']
    g = R(G.row(0).list())

    if probabilistic:
        i = 0
        while i < k or g.degree() > n - k:
            p = R(code.random_element().list())
            while(p == R.zero()):
                p = R(code.random_element().list())
            g = gcd(g, p)
            i += 1
    else:
        i = 1
        while(g.degree() > k - 1 and i < k):
            p = R(G.row(i).list())
            g = gcd(g, p)
            i += 1
    if g.degree() != n - k :
        raise ValueError("Impossible to find a generator polynomial")
    if not g.divides(R.gen() ** n - 1):
         raise ValueError("Impossible to find a generator polynomial")
    return g

def _cyclotomic_coset(n, r, q):
    r"""
    Returns the ``q``-cyclotomic coset of ``r`` modulo ``n``.

    INPUT:

    - ``n``, ``r``, ``q`` -- integers

    OUTPUT:

    - a list of integers

    AUTHORS:

    This function is taken from codinglib [Nielsen]_

    EXAMPLES::

        sage.coding.cyclic_code._cyclotomic_coset(11, 7, 3)
        [8, 10, 2, 6, 7]
    """
    r = r%n
    cyc = set([r])
    rem = (r*q) % n
    while not rem in cyc:
        cyc.add(rem)
        rem = (rem*q) % n
    return list(cyc)

def _complete_list(l, length):
    r"""
    Returns ``l`` with as many zeros appened to it as necessary to make
    it a list of size ``length``.

    INPUT:

    - ``l`` -- a list

    - ``length`` -- an integer

    OUTPUT:

    - a list formed by ``l`` and its completion to ``length``

    EXAMPLES::

        sage: l = [1, 2, 3, 4]
        sage: sage.coding.cyclic_code._complete_list(l, 6)
        [1, 2, 3, 4, 0, 0]
    """
    F = l[0].base_ring()
    l_len = len(l)
    if(l_len != length):
        l = l + ([F.zero()] * (length - l_len))
    return l

def _build_chain_dictionary(D, n):
    r"""
    Returns the dictionnary containing the length of the arithmetic chain for each couple
    ``(d, delta)`` where ``d`` is in ``D`` and ``delta`` is the step of the chain.

    Let `D` be a list of integers, `n` a positive integer. For each `d` in `D`, for each `\delta`
    smaller than `n` and coprime with `n`, we fill the following dictionary, called `\Phi`, such
    that: `\Phi[d, \delta] = l`, where `l` is the biggest integer such that
    `\{d+\delta, d+2\times\delta, \dots, d+l\times\delta\} \subseteq D`.

    INPUT:

    - ``D`` -- a list of integers

    - ``n`` -- an integer

    OUTPUT:

    - a dictionnary

    .. NOTE::

        This is a helper function, used only in :meth:`CyclicCode.bch_bound` and
        :meth:`CyclicCode.hartmann_tzeng_bound`.

    EXAMPLES::

        sage: D = [0, 1, 2, 4]
        sage: n = 5
        sage: sage.coding.cyclic_code._build_chain_dictionary(D, n)
        {(0, 1): 3,
         (0, 2): 4,
         (0, 3): 1,
         (0, 4): 2,
         (1, 1): 2,
         (1, 2): 1,
         (1, 3): 4,
         (1, 4): 3,
         (2, 1): 1,
         (2, 2): 3,
         (2, 3): 2,
         (2, 4): 4,
         (4, 1): 4,
         (4, 2): 2,
         (4, 3): 3,
         (4, 4): 1}
    """
    phi = {}
    for delta in range(1, n):
        if gcd(delta, n) == 1:
            for i in D:
                phi[i, delta] = _fill_chain_dictionary(phi, D, i, delta, n)
    return phi

def _fill_chain_dictionary(phi, D, i, delta, n):
    r"""
    Returns the associated value of ``(i, delta)`` for ``phi``.

    See :meth:`build_chain_dictionary` for details.

    INPUT:

    - ``phi`` -- a dictionary

    - ``D`` -- a list of integers

    - ``i``, ``delta``, ``n`` -- integers

    OUTPUT:

    - an integer

    .. NOTE::

        This is a helper method, for internal use only

    EXAMPLES::

        sage: phi = {}
        sage: D = [0, 1, 2, 4]
        sage: i, delta, n = 0, 1, 5
        sage: sage.coding.cyclic_code._fill_chain_dictionary(phi, D, i, delta, n)
        3
    """
    if not i % n in D:
        return 0
    else:
        try:
            return 1 + phi[i + delta, delta]
        except KeyError:
            return 1 + _fill_chain_dictionary(phi, D, i + delta, delta, n)

def bch_bound(n, D, arithmetic = False, bch_parameters = False):
    r"""
    Returns the BCH bound obtained for a cyclic code of length ``n`` and defining set ``D``.

    Considering a cyclic code `C`, with defining set `D`, length `n`, and minimum
    distance `d`. We have the following bound, called BCH bound, on `d`:
    `d \geq \delta + 1`, where `\delta` is the length of the longest chain of
    consecutive elements of `D` modulo `n`.

    We can also see the BCH bound as an arithmetic sequence: with the same parameters as above,
    if `\exists c, \gcd(c,n)=1` such that `\{l, l+c, \dots, l+\delta\times c\} \subseteq D`,
    `l \in D`, then `d \geq \delta + 1` [1]

    .. NOTE::

        As this is a specific use case of the BCH bound, it is *not* available if the global namespace.
        Call it by using ``sage.coding.cyclic_code.bch_bound``. You can also load it into the global
        namespace by typing ``from sage.coding.cyclic_code import bch_bound``.

    INPUT:

    - ``n`` -- an integer

    - ``D`` -- a list of integers

    - ``arithmetic`` -- (default: ``False``), if it is set to ``True`` it computes the BCH bound using the
      longest arithmetic sequence definition

    EXAMPLES::

        sage: n = 15
        sage: D = [14,1,2,11,12]
        sage: sage.coding.cyclic_code.bch_bound(n, D)
        3

        sage: n = 15
        sage: D = [14,1,2,11,12]
        sage: sage.coding.cyclic_code.bch_bound(n, D, True)
        4
    """
    if arithmetic == True:
        phi = _build_chain_dictionary(D, n)
        val = phi.values()
        longest = 0
        for i in val:
            if longest < i:
                longest = i

        if bch_parameters == True:
            for i in phi.items():
                if i[1] == longest:
                    b = i[0][0]
                    l = i[0][1]
                    break
            bch_params = (b, l)

    else:
        longest = 1
        length = len(D)
        stop = False
        j = 0
        while not stop:
            current_len = 1
            incr = True
            while (D[j % length] + 1) % n == D[(j+1) % length]:
                current_len += 1
                j += 1
                incr = False
            if current_len > longest:
                longest = current_len
            if incr:
                j += 1
            if j >= length and incr == True:
                stop = True

    if bch_parameters == True and arithmetic == True:
        return longest + 1, bch_params
    return longest + 1









class CyclicCode(AbstractLinearCode):
    r"""
    Representation of a cyclic code.

    We propose three different ways to create a new CyclicCode, either by providing:

    - the generator polynomial and the length (1)
    - an existing linear code. In that case, a generator polynomial will be computed
       from the provided linear code's parameters (2)
    - the defining set of the cyclic code (3)

    Depending on which behaviour you want, you need to specify the names of the arguments to
    `CyclicCode`. See EXAMPLES section below for details.

    INPUT:

    - ``generator_pol`` -- (default: ``None``) the unique monic polynomial which divides every
      codeword of ``self``.

    - ``length`` -- (default: ``None``) the length of ``self``. It has to be bigger than the degree
      of ``generator_pol``.

    - ``code`` -- (default: ``None``) a linear code.

    - ``probabilistic`` -- (default: ``False``) enables probabilistic version of the
      generator polynomial from any linear code lookup algorithm.
      See :meth:`sage.find_generator_polynomial` for details.

    - ``D`` -- (default: ``None``) the set of powers of ``primitive element``. It can be partial.

    - ``field`` -- (default: ``None``) the base field of ``self``.

    - ``primitive_element`` -- (default: ``None``) the primitive element
      to use when creating the set of roots for the generating polynomial
      over the splitting field. It has to be of multiplicative order ``length`` over this
      field. If the splitting field is not ``field``, it also have to be a polynomial in ``zx``,
      where ``x`` is the degree of the extension field. For instance,
      over ``GF(16)``, it has to be a polynomial in ``z4``.

    EXAMPLES:

    We can construct a Cyclic Code using three different methods.
    First (1), we provide a generator_polynomial and a length for the code::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol = g, length = n)
        sage: C
        [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial

    We can also provide a code (2). In that case, the program will try to extract a generator
    polynomial for the provided code (see :meth:`find_generator_polynomial` for details)::

        sage: C = codes.GeneralizedReedSolomonCode(GF(2 ** 3, 'a').list()[1:2 ** 3], 2 ** 2)
        sage: Cc = codes.CyclicCode(code = C)
        sage: Cc
        [7, 4] Cyclic Code over Finite Field in a of size 2^3 with x^3 + (a^2 + 1)*x^2 + a*x + a^2 + 1 as generator polynomial

    In that case, it is possible to choose between a probabilistic algorithm and a deterministic
    one to extract a generator polynomial, using the keyword ``probabilistic``::

        sage: C = codes.GeneralizedReedSolomonCode(GF(2 ** 3, 'a').list()[1:2 ** 3], 2 ** 2)
        sage: Cc = codes.CyclicCode(code = C, probabilistic = True)
        sage: Cc
        [7, 4] Cyclic Code over Finite Field in a of size 2^3 with x^3 + (a^2 + 1)*x^2 + a*x + a^2 + 1 as generator polynomial

    We can also provide a defining set for the code (3). In that case, the generator polynomial
    will be computed::

        sage: F = GF(16, 'a')
        sage: n = 15
        sage: Cc = codes.CyclicCode(length = n, field = F, D = [1,2])
        sage: Cc
        [15, 13] Cyclic Code over Finite Field in a of size 2^4 with x^2 + (a^2 + a)*x + a^3 as generator polynomial
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, generator_pol=None, length=None, code=None, probabilistic = False, D = None, field = None, primitive_element = None):
        r"""
        TESTS:

        If one provides a polynomial and the dimension, we check that provided length
        is bigger than provided polynomial's degree::

            sage: F.<x> = GF(2)[]
            sage: n = 2
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            Traceback (most recent call last):
            ...
            ValueError: Length must be bigger than generator polynomial's degree

        We check that provided polynomial is defined over a finite field::

            sage: F.<x> = RR[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            Traceback (most recent call last):
            ...
            ValueError: Generator polynomial must be defined over a finite field

        And we check that the provided polynomial divides `x^{n} - 1`, where `n` is provided length::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 2 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            Traceback (most recent call last):
            ...
            ValueError: Provided polynomial must divide x^n - 1, where n is the provided length

        In the case of a code is passed as argument, if it's not possible to extract a
        generator polynomial, an exception is raised::

            sage: C = codes.GeneralizedReedSolomonCode(GF(2 ** 4, 'a').list()[1:2 ** 4], 2 ** 4 - 1)
            sage: Cc = codes.CyclicCode(code = C)
            Traceback (most recent call last):
            ...
            ValueError: Impossible to find a generator polynomial
        """
        if generator_pol is not None and length is not None:
            F = generator_pol.base_ring()
            R = F[generator_pol.variable_name()]
            deg = generator_pol.degree()
            if not isinstance(length, Integer):
                length = Integer(length)
            if length <= deg:
                raise ValueError("Length must be bigger than generator polynomial's degree")
            if not generator_pol.base_ring().is_finite():
                raise ValueError("Generator polynomial must be defined over a finite field")
            if not generator_pol.divides(R.gen() ** length - 1):
                raise ValueError("Provided polynomial must divide x^n - 1, where n is the provided length")
            self._dimension = length - deg
            self._generator_polynomial = generator_pol
            super(CyclicCode, self).__init__(F, length, "Vector", "Syndrome")
        elif code is not None:
            try:
                g = find_generator_polynomial(code, probabilistic)
            except ValueError, e:
                raise ValueError(e)
            if not g.is_monic():
                g = g.monic()
            self._generator_polynomial = g
            self._dimension = code.dimension()
            super(CyclicCode, self).__init__(code.base_ring(), code.length(), "Vector", "Syndrome")
        elif D is not None and length is not None and field is not None:
            F = field
            n = length
            q = F.cardinality()
            #basic checks over the input
            if not gcd(n, q) == 1:
                raise ValueError("n and q must be coprimes")
            #creation of the extension field (the splitting field)
            R = F['x']
            x = R.gen()
            s = 1
            while not n.divides(q ** s - 1):
                s += 1

            #If s equals to 1, then Fsplit is the same field that F, so we do not need to build
            #the extension field to compute the generator polynomial
            if s == 1:
                pows = set()
                for r in D:
                    if not r in pows:
                        pows = pows.union(_cyclotomic_coset(n, r, q))

                if primitive_element is not None and (primitive_element not in F or
                        multiplicative_order(primitive_element) != n):
                    raise ValueError("primitive_element has to be an element of multplicative order n in the extension field used to compute the generator polynomial")
                elif primitive_element is not None:
                    alpha = primitive_element
                else:
                    alpha = F.zeta(n)
                self._primitive_element = alpha
                alpha = F.gen() ** ((q - 1) // n)
                g = R(prod(x - alpha ** p for p in pows))

            else:
                Fsplit, F_to_Fsplit = F.extension(Integer(s), map = True)
                FE = RelativeFiniteFieldExtension(Fsplit, F, embedding = F_to_Fsplit)
                if primitive_element is not None and (primitive_element not in Fsplit or
                        multiplicative_order(primitive_element) != n):
                    raise ValueError("primitive_element has to be an element of multplicative order n in the extension field used to compute the generator polynomial")
                elif primitive_element is not None:
                    beta = primitive_element
                else:
                    beta = Fsplit.zeta(n)
                self._primitive_element = beta
                Rsplit = Fsplit['xx']
                xx = Rsplit.gen()

                #we compute the generator polynomial over the large field (Fsplit)
                #for this, we need the powers of the roots, so we build the cyclotomic cosets
                pows = set()
                gsplit = Rsplit(1)
                cosets = []
                for r in D:
                    if not r in pows:
                        current = _cyclotomic_coset(n, r, q)
                        cosets.append(current)
                        pows = pows.union(current)

                # then we compute the minimal polynomial to each coset
                min_pols = []
                for i in cosets:
                    pol = Rsplit.one()
                    for j in i:
                        pol = pol * (xx - beta**j)
                    min_pols.append(pol)

                #finally, we get back to the small field and compute the generator polynomial
                R = F['x']
                pols_coeffs = []
                for i in min_pols:
                    tmp = []
                    for j in i:
                        tmp.append(sum(FE.relative_field_representation(j)))
                    pols_coeffs.append(tmp)
                g = R.one()
                for i in pols_coeffs:
                    g = g * R(i)
            #we set class variables (and store some of the things we computed before)
            self._defining_set = list(pows)
            self._defining_set.sort()
            self._generator_polynomial = g
            self._dimension = n - g.degree()
            super(CyclicCode, self).__init__(F, n, "Vector", "Syndrome")

        else:
            raise AttributeError("You must provide either a code, or a list of powers and the length\
                    and the field, or a generator polynomial and the code length")

    def __contains__(self, word):
        r"""
        Returns ``True`` if ``word`` belongs to ``self``, ``False`` otherwise.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: c in C
            True
        """
        g = self.generator_polynomial()
        R = g.parent()
        return g.divides(R(word.list()))

    def __eq__(self, other):
        r"""
        Tests equality between Cyclic Code objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C1 = codes.CyclicCode(generator_pol = g, length = n)
            sage: C2 = codes.CyclicCode(generator_pol = g, length = n)
            sage: C1 == C2
            True
        """
        return isinstance(other, CyclicCode) \
                and self.base_field() == other.base_field() \
                and self.length() == other.length() \
                and self.generator_polynomial() == other.generator_polynomial() \

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C
            [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
        """
        return "[%s, %s] Cyclic Code over %s with %s as generator polynomial"\
                % (self.length(), self.dimension(),\
                self.base_field(), self.generator_polynomial())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: latex(C)
            [7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2} \textnormal{ with } x^{3} + x + 1 \textnormal{ as generator polynomial}
        """
        return "[%s, %s] \\textnormal{ Cyclic Code over } %s \\textnormal{ with } %s \\textnormal{ as generator polynomial}"\
                % (self.length(), self.dimension(),\
                self.base_field()._latex_(), self.generator_polynomial()._latex_())

    def generator_polynomial(self):
        r"""
        Returns the generator polynomial of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.generator_polynomial()
            x^3 + x + 1
        """
        return self._generator_polynomial

    def primitive_element(self):
        r"""
        Returns the primitive element that was used as a root of
        the generator polynomial over the extension field.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.primitive_element()
            z3

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: a = F.gen()
            sage: Cc = codes.CyclicCode(length = n, field = F, D = [1,2], primitive_element = a^2 + 1)
            sage: Cc.primitive_element()
            a^2 + 1
        """
        if hasattr(self, "_primitive_element"):
            return self._primitive_element
        else:
            F = self.base_field()
            q = F.cardinality()
            n = self.length()
            s = 1
            while not n.divides(q ** s - 1):
                s += 1

            Fsplit = F.extension(Integer(s))
            beta = Fsplit.zeta(n)
            self._primitive_element = beta
            return beta

    def check_polynomial(self):
        r"""
        Returns the check polynomial of ``self``.

        Let `C` be a cyclic code of length `n` and `g` its
        generator polynomial.
        The following: `h = \frac{x^n - 1}{g(x)}` is called `C`'s
        check polynomial.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: h = C.check_polynomial()
            sage: h == (x**n - 1)/C.generator_polynomial()
            True
        """
        if hasattr(self, "_check_polynomial"):
            return self._check_polynomial
        R = self.base_ring()['x']
        n = self.length()
        self._check_polynomial = (R.gen() ** n - 1) // self.generator_polynomial()
        return self._check_polynomial

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

        The parity check matrix of a linear code `C` corresponds to the
        generator matrix of the dual code of `C`.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.parity_check_matrix()
            [1 0 1 1 1 0 0]
            [0 1 0 1 1 1 0]
            [0 0 1 0 1 1 1]
        """
        k = self.dimension()
        n = self.length()
        R = self.base_ring()['x']
        h = self.check_polynomial()
        l = h.coefficients(sparse = False)
        l.reverse()
        l = _complete_list(l, n)
        H = matrix(self.base_ring(), n - k, n)
        for i in range(n - k):
            H.set_row(i, l)
            l = l[n-1:] + l[:n-1]
        return H

    def defining_set(self):
        r"""
        Returns the set of powers of the root of ``self``'s generator polynomial
        over the extension field.

        EXAMPLES:

        We provide a defining set at construction time::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: C = codes.CyclicCode(length = n, field = F, D = [1,2])
            sage: C.defining_set()
            [1, 2]

        If the defining set was provided by the user, it might have been expanded
        at construction time. In this case, the expanded defining set will be returned::

            sage: C = codes.CyclicCode(length = 13, field = F, D = [1, 2])
            sage: C.defining_set()
            [1, 2, 3, 5, 6, 9]

        If a generator polynomial was passed at construction time,
        the defining set is computed using this polynomial::

            sage: F.<x> = GF(8, 'a')[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.defining_set()
            [1, 2, 4]

        Both operations give the same result::

            sage: C1 = codes.CyclicCode(length = n, field = GF(8, 'a'), D = [1, 2, 4])
            sage: C1.generator_polynomial() == g
            True

        Another one, in the revert order::

            sage: F = GF(16, 'a')
            sage: n = 13
            sage: C1 = codes.CyclicCode(length = n, field = F, D = [1, 2])
            sage: g = C1.generator_polynomial()
            sage: C2 = codes.CyclicCode(generator_pol = g, length = n)
            sage: C1.defining_set() == C2.defining_set()
            True
        """
        if hasattr(self, "_defining_set"):
            return self._defining_set
        else:
            D = []
            F = self.base_field()
            n = self.length()
            q = F.cardinality()
            g = self.generator_polynomial()
            #creation of the extension field (the splitting field)
            #and embeddings
            s = 1
            while not n.divides(q ** s - 1):
                s += 1

            Fsplit, F_to_Fsplit = F.extension(Integer(s), 'b', map = True)
            beta = Fsplit.zeta(n)

            #computation of the roots of the generator polynomial
            #over Fsplit
            gsplit = []
            Rsplit = Fsplit['xx']
            for i in g.coefficients(sparse = False):
                gsplit.append(F_to_Fsplit(i))
            gsplit = Rsplit(gsplit)
            roots = gsplit.roots(multiplicities = False)

            #recovering defining set
            for i in roots:
                D.append(discrete_log(i, beta))

            D.sort()
            self._defining_set = D
            return D

    def bch_bound(self, arithmetic = False, bch_parameters = False):
        r"""
        Returns the BCH bound of self which is a bound on ``self``'s minimum distance.

        See :meth:`sage.coding.cyclic_code.bch_bound` for details.

        INPUT:

        - ``F`` -- a finite field

        - ``n`` -- an integer

        - ``D`` -- a list of integers

        - ``arithmetic`` -- (default: ``False``), if it is set to ``True`` it computes the BCH bound
          using the longest arithmetic sequence definition

        EXAMPLES::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: D = [14,1,2,11,12]
            sage: C = codes.CyclicCode(field = F, length = n, D = D)
            sage: C.bch_bound()
            3

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: D = [14,1,2,11,12]
            sage: C = codes.CyclicCode(field = F, length = n, D = D)
            sage: C.bch_bound(True)
            4
        """
        return bch_bound(n = self.length(), D = self.defining_set(), arithmetic = arithmetic,\
                bch_parameters = bch_parameters)










class CyclicCodePolynomialEncoder(Encoder):
    r"""
    An encoder encoding polynomials into codewords.

    Let `C` be a cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.

    This encoder encodes any polynomial `p \in F[x]_{<k}`, by doing:
    `\forall p \in F[x]_{<k}, (p \times g) \in C`

    INPUT:

    - ``code`` -- The associated code of this encoder

    EXAMPLES::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol = g, length = n)
        sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
        sage: E
        Polynomial-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E
            Polynomial-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial

        """
        super(CyclicCodePolynomialEncoder, self).__init__(code)
        self._R = code.base_field()['x']

    def __eq__(self, other):
        r"""
        Tests equality between CyclicCodePolynomialEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E1 = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E2 = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E1 == E2
            True
        """
        return isinstance(other, CyclicCodePolynomialEncoder) \
                and self.code() == other.code()

    def __ne__(self, other):
        r"""
        Tests difference between CyclicCodePolynomialEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: g1 = x ** 3 + x + 1
            sage: g2 = (x ** 5 + 1) * (x ** 4 + x + 1)
            sage: C1 = codes.CyclicCode(generator_pol = g1, length = 7)
            sage: C2 = codes.CyclicCode(generator_pol = g2, length = 15)
            sage: E1 = codes.encoders.CyclicCodePolynomialEncoder(C1)
            sage: E2 = codes.encoders.CyclicCodePolynomialEncoder(C2)
            sage: E1 != E2
            True
        """
        return not self.__eq__(other)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E
            Polynomial-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
        """
        return "Polynomial-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: latex(E)
            \textnormal{Polynomial-style encoder for }[7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2} \textnormal{ with } x^{3} + x + 1 \textnormal{ as generator polynomial}
        """
        return "\\textnormal{Polynomial-style encoder for }%s" % self.code()._latex_()

    def encode(self, p):
        r"""
        Transforms ``p`` into an element of the associated code of ``self``.

        INPUT:

        - ``p`` -- A polynomial from ``self`` message space

        OUTPUT:

        - A codeword in associated code of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: m = x ** 2 + 1
            sage: E.encode(m)
            (1, 1, 1, 0, 0, 1, 0)
        """
        k = self.code().dimension()
        n = self.code().length()
        if p.degree() > k:
            raise ValueError("Degree of the message must be at most %s" % k)
        res = (p * self.code().generator_polynomial()).coefficients(sparse = False)
        res = _complete_list(res, n)
        return  vector(self.code().base_field(), res)

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- A vector with the same length as the code

        OUTPUT:

        - An element of the message space

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: E.unencode_nocheck(c)
            x^2 + 1
        """
        R = self.message_space()
        g = self.code().generator_polynomial()
        p = R(c.list())
        return p // g

    def message_space(self):
        r"""
        Returns the message space of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E.message_space()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)
        """
        return self._R










class CyclicCodeVectorEncoder(Encoder):
    r"""
    An encoder which can encode vectors into codewords.

    Let `C` be a cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.

    Let `m = (m_1, m_2, \dots, m_{k-1})` be a vector in `F^{k}`.
    This codeword can be seen as a polynomial over `F_k[x]`, as follows:
    `m_p = \Sigma_{i=0}^{k-1} m_i \times x^i`.

    To encode m, this encoder does the following multiplication: `m_p \times g`.

    INPUT:

    - ``code`` -- The associated code of this encoder

    EXAMPLES::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol = g, length = n)
        sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
        sage: E
        Vector-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
    """

    def __init__(self, code):
        r"""

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E
            Vector-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial

        """
        super(CyclicCodeVectorEncoder, self).__init__(code)

    def __eq__(self, other):
        r"""
        Tests equality between CyclicCodeVectorEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E1 = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E2 = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E1 == E2
            True
        """
        return isinstance(other, CyclicCodeVectorEncoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E
            Vector-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
        """
        return "Vector-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: latex(E)
            \textnormal{Vector-style encoder for }[7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2} \textnormal{ with } x^{3} + x + 1 \textnormal{ as generator polynomial}
        """
        return "\\textnormal{Vector-style encoder for }%s" % self.code()._latex_()

    def encode(self, m):
        r"""
        Transforms ``m`` into an element of the associated code of ``self``.

        INPUT:

        - ``m`` -- A element from ``self``'s message space

        OUTPUT:

        - A codeword in associated code of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: m = vector(GF(2), (1, 0, 1, 0))
            sage: E.encode(m)
            (1, 1, 1, 0, 0, 1, 0)
        """
        if hasattr(self, "_known_generator_matrix"):
            return super(CyclicCodeVectorEncoder, self).encode(m)

        k = self.code().dimension()
        n = self.code().length()
        F  = self.code().base_field()
        R = F['x']
        p = R(m.list())
        if p.degree() > k:
            raise ValueError("Degree of the message must be at most %s" % k)
        res = (p * self.code().generator_polynomial()).coefficients(sparse = False)
        res = _complete_list(res, n)
        return vector(F, res)

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- A vector with the same length as the code

        OUTPUT:

        - An element of the message space

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: E.unencode_nocheck(c)
            (1, 0, 1, 0)
        """

        R = self.message_space().base_ring()['x']
        g = self.code().generator_polynomial()
        p = R(c.list())
        l = (p//g).coefficients(sparse = False)
        l = _complete_list(l, self.message_space().dimension())
        return vector(self.code().base_field(), l)

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E.generator_matrix()
            [1 1 0 1 0 0 0]
            [0 1 1 0 1 0 0]
            [0 0 1 1 0 1 0]
            [0 0 0 1 1 0 1]
        """

        self._known_generator_matrix = True
        C = self.code()
        k = C.dimension()
        n = C.length()
        l = C.generator_polynomial().coefficients(sparse = False)
        l = _complete_list(l, n)
        G = matrix(C.base_ring(), k, C.length())
        for i in range(k):
            G.set_row(i, l)
            l = l[n-1:] + l[:n-1]
        return G

    def message_space(self):
        r"""
        Returns the message space of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E.message_space()
            Vector space of dimension 4 over Finite Field of size 2
        """
        return self.code().base_ring() ** self.code().dimension()


####################### registration ###############################

CyclicCode._registered_encoders["Polynomial"] = CyclicCodePolynomialEncoder
CyclicCode._registered_encoders["Vector"] = CyclicCodeVectorEncoder
CyclicCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
CyclicCode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder
