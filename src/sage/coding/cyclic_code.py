r"""
Cyclic Code

Let `F` be a field. A `[n ,k]` code `C` over `F` is called cyclic if every cyclic shift
of a codeword in `C` is also a codeword [R]:

    .. MATH::
        \begin{aligned}
        \forall c \in C \\
        c = (c_{0}, c_{1}, \dots , c_{n-1}) \in C
        \Rightarrow (c_{n-1}, c_{0}, \dots , c_{n-2})
        \end{aligned}

REFERENCES:

.. [HT] C. Hartmann and K.K. Tzeng. Generalizations of the BCH Bound.
   Information and Control, 20(5):489-498, June 1972

.. [R2] Introduction to Coding Theory, Ron Roth, Cambridge University Press, 2006



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
from field_embedding import FieldEmbedding
from sage.groups.generic import discrete_log

def find_generator_polynomial(code, probabilistic = False):
    r"""
    Returns a possible generator polynomial for ``code``.

    This method is probabilistic. If it does not find a good candidate in $k$
    attempts, where k is ``code``'s dimension, it returns an exception.

    INPUT:

    - ``code`` -- a linear code

    - ``probabilistic`` -- (default : ``False``) enables the probabilistic version of
      this method

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

def cyclotomic_coset(n, r, q):
    r"""
    Returns the q-cyclotomic coset of ``r`` modulo ``n``.

    INPUT:

    - ``n``, ``r``, ``q`` -- integers

    OUTPUT:

    - a list of integers

    AUTHORS:

    This function is taken from codinglib (https://bitbucket.org/jsrn/codinglib/)
    and was written by Johan Nielsen.

    EXAMPLES::

        sage.coding.cyclic_code.cyclotomic_coset(11, 7, 3)
        [8, 10, 2, 6, 7]
    """
    r = r%n
    cyc = set([r])
    rem = (r*q) % n
    while not rem in cyc:
        cyc.add(rem)
        rem = (rem*q) % n
    return list(cyc)

def complete_list(l, length):
    r"""
    Returns ``l`` with as many zeros appened to it as necessary for list to be of size ``length``.

    INPUT:

    - ``l`` -- a list

    - ``length`` -- an integer

    OUTPUT:

    - a list formed by ``l`` and its completion to ``length``

    EXAMPLES::

        sage: l = [1, 2, 3, 4]
        sage: sage.coding.cyclic_code.complete_list(l, 6)
        [1, 2, 3, 4, 0, 0]
    """
    F = l[0].base_ring()
    l_len = len(l)
    if(l_len != length):
        l = l + ([F.zero()] * (length - l_len))
    return l

def build_chain_dictionary(D, n):
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
        sage: sage.coding.cyclic_code.build_chain_dictionary(D, n)
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
                phi[i, delta] = fill_chain_dictionary(phi, D, i, delta, n)
    return phi

def fill_chain_dictionary(phi, D, i, delta, n):
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
        sage: sage.coding.cyclic_code.fill_chain_dictionary(phi, D, i, delta, n)
        3
    """
    if not i % n in D:
        return 0
    else:
        try:
            return 1 + phi[i + delta, delta]
        except KeyError:
            return 1 + fill_chain_dictionary(phi, D, i + delta, delta, n)

def bch_bound(n, D, arithmetic = False, bch_parameters = False):
    r"""
    Returns the BCH bound obtained for a cyclic code of length ``n`` and defining set ``D``.

    Considering a cyclic code `C`, with defining set `D`, length `n`, and minimum
    distance `d`. We have the following bound, called BCH bound on `d`:
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
        phi = build_chain_dictionary(D, n)
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

    ..NOTE::

        We propose three different ways to create a new CyclicCode, either by providing:
        (1) the generator polynomial and the length
        (2) an existing linear code. In that case, a generator polynomial will be computed
          from the provided linear code's parameters
        (3) the defining set of the cyclic code

        Depending on what behaviour you want, you need to specify the names of the arguments to
        CyclicCode. See EXAMPLES section below for details.

    INPUT:

    - ``generator_pol`` -- (default: ``None``) a generator polynomial for this code. It has to
      be monic.

    - ``length`` -- (default: ``None``) the length of the code. It has to be bigger than the degree
      of ``generator_pol``.

    - ``code`` -- (default: ``None``) a linear code.

    - ``probabilistic`` -- (default: ``False``) enables probabilistic version of the
      generator polynomial lookup algorithm.

    - ``D`` -- (default: ``None``) the defining set of the code. It can be partial.

    - ``field`` -- (default: ``None``) the base field for this code.

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

    def __init__(self, generator_pol=None, length=None, code=None, probabilistic = False, D = None, field = None):
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
                        pows = pows.union(cyclotomic_coset(n, r, q))

                alpha = F.gen() ** ((q - 1) // n)
                g = R(prod(x - alpha ** p for p in pows))

            else:
                Fsplit, F_to_Fsplit = F.extension(Integer(s), 'b', map = True)
                FE = FieldEmbedding(Fsplit, F, embedding = F_to_Fsplit)
                beta = Fsplit.zeta(n)
                Rsplit = Fsplit['xx']
                xx = Rsplit.gen()

                #we compute the generator polynomial over the large field (Fsplit)
                #for this, we need the powers of the roots, so we build the cyclotomic cosets
                pows = set()
                gsplit = Rsplit(1)
                cosets = []
                for r in D:
                    if not r in pows:
                        current = cyclotomic_coset(n, r, q)
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
                        tmp.append(FE.small_field_polynomial_representation(j))
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

    def __ne__(self, other):
        r"""
        Tests inequality of Cyclic Code objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: g1 = x ** 3 + x + 1
            sage: g2 = (x ** 5 + 1) * (x ** 4 + x + 1)
            sage: C1 = codes.CyclicCode(generator_pol = g1, length = 7)
            sage: C2 = codes.CyclicCode(generator_pol = g2, length = 15)
            sage: C1 != C2
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

    def check_polynomial(self):
        r"""
        Returns the check polynomial of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.check_polynomial()
            x^4 + x^2 + x + 1
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
        l = complete_list(l, n)
        H = matrix(self.base_ring(), n - k, n)
        for i in range(n - k):
            H.set_row(i, l)
            l = l[n-1:] + l[:n-1]
        return H

    def defining_set(self):
        r"""
        Returns the defining set of ``self``.

        If it was computed at construction time, it returns immediately the computed one,
        else it is computed using the generator polynomial built at construction time.

        EXAMPLES:

        We provide a defining set at construction time::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: C = codes.CyclicCode(length = n, field = F, D = [1,2])
            sage: C.defining_set()
            [1, 2]

        If the defining was expanded while computing cyclotomic classes, the
        expanded defining set will be returned::

            sage: C = codes.CyclicCode(length = 13, field = F, D = [1, 2])
            sage: C.defining_set()
            [1, 2, 3, 5, 6, 9]

        If a generator polynomial was passed at construction time,
        the defining set is computed by this method::

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
            gsplit = Rsplit(gsplit) #WARNING: inverses order of coefficients !!
            roots = gsplit.roots(multiplicities = False)

            #recovering defining set
            for i in roots:
                D.append(discrete_log(beta, i))

            D.sort()
            self._defining_set = D
            return D

    def bch_bound(self, arithmetic = False, bch_parameters = False):
        r"""
        Returns the BCH bound of self. See :meth:`sage.coding.cyclic_code.bch_bound` for details.

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
    An encoder which can encode polynomials into codewords.

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
        res = complete_list(res, n)
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

    def __ne__(self, other):
        r"""
        Tests difference between CyclicCodeVectorEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: g1 = x ** 3 + x + 1
            sage: g2 = (x ** 5 + 1) * (x ** 4 + x + 1)
            sage: C1 = codes.CyclicCode(generator_pol = g1, length = 7)
            sage: C2 = codes.CyclicCode(generator_pol = g2, length = 15)
            sage: E1 = codes.encoders.CyclicCodeVectorEncoder(C1)
            sage: E2 = codes.encoders.CyclicCodeVectorEncoder(C2)
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
        res = complete_list(res, n)
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
        l = complete_list(l, self.message_space().dimension())
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
        l = complete_list(l, n)
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
