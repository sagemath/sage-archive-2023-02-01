# -*- coding: utf-8 -*-
r"""
Dirichlet characters

A :class:`DirichletCharacter` is the extension of a homomorphism

.. math::

    (\ZZ/N\ZZ)^* \to R^*,

for some ring `R`, to the map `\ZZ/N\ZZ \to R` obtained by sending
those `x\in\ZZ/N\ZZ` with `\gcd(N,x)>1` to `0`.

EXAMPLES::

    sage: G = DirichletGroup(35)
    sage: x = G.gens()
    sage: e = x[0]*x[1]^2; e
    Dirichlet character modulo 35 of conductor 35 mapping 22 |--> zeta12^3, 31 |--> zeta12^2 - 1
    sage: e.order()
    12

This illustrates a canonical coercion::

    sage: e = DirichletGroup(5, QQ).0
    sage: f = DirichletGroup(5,CyclotomicField(4)).0
    sage: e*f
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -zeta4

AUTHORS:

-  William Stein (2005-09-02): Fixed bug in comparison of Dirichlet
   characters. It was checking that their values were the same, but
   not checking that they had the same level!

-  William Stein (2006-01-07): added more examples

-  William Stein (2006-05-21): added examples of everything; fix a
   *lot* of tiny bugs and design problem that became clear when
   creating examples.

-  Craig Citro (2008-02-16): speed up __call__ method for
   Dirichlet characters, miscellaneous fixes

-  Julian Rueth (2014-03-06): use UniqueFactory to cache DirichletGroups

"""

#*****************************************************************************
#       Copyright (C) 2004-2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2014 Julian Rueth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.categories.all                  as cat
from sage.misc.all import prod
import sage.misc.prandom                    as random
import sage.modules.free_module             as free_module
import sage.modules.free_module_element     as free_module_element
import sage.rings.all                       as rings
import sage.rings.number_field.number_field as number_field

from sage.categories.map import Map
from sage.rings.rational_field import is_RationalField
from sage.rings.complex_field import is_ComplexField
from sage.rings.qqbar import is_AlgebraicField
from sage.rings.ring import is_Ring

from sage.misc.cachefunc                    import cached_method
from sage.misc.fast_methods                 import WithEqualityById
from sage.structure.element                 import MultiplicativeGroupElement
from sage.structure.gens_py                 import multiplicative_iterator
from sage.structure.parent                  import Parent
from sage.structure.sequence                import Sequence
from sage.structure.factory                 import UniqueFactory
from sage.arith.all import (binomial, bernoulli, kronecker, factor, gcd,
        lcm, fundamental_discriminant, euler_phi, factorial, valuation)

def trivial_character(N, base_ring=rings.RationalField()):
    r"""
    Return the trivial character of the given modulus, with values in the given
    base ring.

    EXAMPLE::

        sage: t = trivial_character(7)
        sage: [t(x) for x in [0..20]]
        [0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1]
        sage: t(1).parent()
        Rational Field
        sage: trivial_character(7, Integers(3))(1).parent()
        Ring of integers modulo 3
    """
    return DirichletGroup(N, base_ring)(1)

TrivialCharacter = trivial_character

def kronecker_character(d):
    """
    Returns the quadratic Dirichlet character (d/.) of minimal
    conductor.

    EXAMPLES::

        sage: kronecker_character(97*389*997^2)
        Dirichlet character modulo 37733 of conductor 37733 mapping 1557 |--> -1, 37346 |--> -1

    ::

        sage: a = kronecker_character(1)
        sage: b = DirichletGroup(2401,QQ)(a)    # NOTE -- over QQ!
        sage: b.modulus()
        2401

    AUTHORS:

    - Jon Hanke (2006-08-06)
    """
    d = rings.Integer(d)
    if d == 0:
        raise ValueError("d must be nonzero")

    D = fundamental_discriminant(d)
    G = DirichletGroup(abs(D), rings.RationalField())
    return G([kronecker(D,u) for u in G.unit_gens()])


def kronecker_character_upside_down(d):
    """
    Returns the quadratic Dirichlet character (./d) of conductor d, for
    d0.

    EXAMPLES::

        sage: kronecker_character_upside_down(97*389*997^2)
        Dirichlet character modulo 37506941597 of conductor 37733 mapping 13533432536 |--> -1, 22369178537 |--> -1, 14266017175 |--> 1

    AUTHORS:

    - Jon Hanke (2006-08-06)
    """
    d = rings.Integer(d)
    if d <= 0:
        raise ValueError("d must be positive")

    G = DirichletGroup(d, rings.RationalField())
    return G([kronecker(u.lift(),d) for u in G.unit_gens()])


def is_DirichletCharacter(x):
    r"""
    Return True if x is of type DirichletCharacter.

    EXAMPLES::

        sage: from sage.modular.dirichlet import is_DirichletCharacter
        sage: is_DirichletCharacter(trivial_character(3))
        True
        sage: is_DirichletCharacter([1])
        False
    """
    return isinstance(x, DirichletCharacter)


class DirichletCharacter(MultiplicativeGroupElement):
    """
    A Dirichlet character.
    """
    def __init__(self, parent, x, check=True):
        r"""
        Create a Dirichlet character with specified values on
        generators of `(\ZZ/n\ZZ)^*`.

        INPUT:

        - ``parent`` -- :class:`DirichletGroup`, a group of Dirichlet
           characters

        - ``x`` -- one of the following:

           - tuple or list of ring elements: the values of the
             Dirichlet character on the standard generators of
             `(\ZZ/N\ZZ)^*` as returned by
             :meth:`sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic.unit_gens`.

           - vector over `\ZZ/e\ZZ`, where `e` is the order of the
             standard root of unity for ``parent``.

           In both cases, the orders of the elements must divide the
           orders of the respective generators of `(\ZZ/N\ZZ)^*`.

        OUTPUT:

        The Dirichlet character defined by `x` (type
        :class:`DirichletCharacter`).

        EXAMPLES::

            sage: G.<e> = DirichletGroup(13)
            sage: G
            Group of Dirichlet characters modulo 13 with values in Cyclotomic Field of order 12 and degree 4
            sage: e
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12
            sage: loads(e.dumps()) == e
            True

        ::

            sage: G, x = DirichletGroup(35).objgens()
            sage: e = x[0]*x[1]; e
            Dirichlet character modulo 35 of conductor 35 mapping 22 |--> zeta12^3, 31 |--> zeta12^2
            sage: e.order()
            12
            sage: loads(e.dumps()) == e
            True

        TESTS::

            sage: G = DirichletGroup(10)
            sage: TestSuite(G[1]).run()

        It is checked that the orders of the elements in `x` are
        admissible (see :trac:`17283`)::

            sage: k.<i> = CyclotomicField(4)
            sage: G = DirichletGroup(192)
            sage: G([i, -1, -1])
            Traceback (most recent call last):
            ...
            ValueError: values (= (zeta16^4, -1, -1)) must have multiplicative orders dividing (2, 16, 2), respectively

            sage: from sage.modular.dirichlet import DirichletCharacter
            sage: M = FreeModule(Zmod(16), 3)
            sage: DirichletCharacter(G, M([4, 8, 8]))
            Traceback (most recent call last):
            ...
            ValueError: values (= (4, 8, 8) modulo 16) must have additive orders dividing (2, 16, 2), respectively
        """
        MultiplicativeGroupElement.__init__(self, parent)
        if check:
            orders = parent.integers_mod().unit_group().gens_orders()
            if len(x) != len(orders):
                raise ValueError("wrong number of values (= {}) on generators (want {})".format(x, len(orders)))
            if free_module_element.is_FreeModuleElement(x):
                x = parent._module(x)
                if any(map(lambda u, v: v*u != 0, x, orders)):
                    raise ValueError("values (= {} modulo {}) must have additive orders dividing {}, respectively"
                                     .format(x, parent.zeta_order(), orders))
                self.element.set_cache(x)
            else:
                R = parent.base_ring()
                x = tuple(map(R, x))
                if R.is_exact() and any(map(lambda u, v: u**v != 1, x, orders)):
                    raise ValueError("values (= {}) must have multiplicative orders dividing {}, respectively"
                                     .format(x, orders))
                self.__values_on_gens = x
        else:
            if free_module_element.is_FreeModuleElement(x):
                self.element.set_cache(x)
            else:
                self.__values_on_gens = x

    @cached_method
    def __eval_at_minus_one(self):
        r"""
        Efficiently evaluate the character at -1 using knowledge of its
        order. This is potentially much more efficient than computing the
        value of -1 directly using dlog and a large power of the image root
        of unity.

        We use the following. Proposition: Suppose eps is a character mod
        `p^n`, where `p` is a prime. Then
        `\varepsilon(-1) = -1` if and only if `p = 2` and
        the factor of eps at 4 is nontrivial or `p > 2` and 2 does
        not divide `\phi(p^n)/\mbox{\rm ord}(\varepsilon)`.

        EXAMPLES::

            sage: chi = DirichletGroup(20).0; chi._DirichletCharacter__eval_at_minus_one()
            -1
        """
        D = self.decomposition()
        val = self.base_ring()(1)
        for e in D:
            if e.modulus() % 2 == 0:
                if e.modulus() % 4 == 0:
                    val *= e.values_on_gens()[0] # first gen is -1 for 2-power modulus
            elif (euler_phi(e.parent().modulus()) / e.order()) % 2 != 0:
                val *= -1
        return val

    def __call__(self, m):
        """
        Return the value of this character at the integer `m`.

        .. warning::

           A table of values of the character is made the first time
           you call this (unless `m` equals -1)

        EXAMPLES::

            sage: G = DirichletGroup(60)
            sage: e = prod(G.gens(), G(1))
            sage: e
            Dirichlet character modulo 60 of conductor 60 mapping 31 |--> -1, 41 |--> -1, 37 |--> zeta4
            sage: e(-1)
            -1
            sage: e(2)
            0
            sage: e(7)
            -zeta4
            sage: Integers(60).unit_gens()
            (31, 41, 37)
            sage: e(31)
            -1
            sage: e(41)
            -1
            sage: e(37)
            zeta4
            sage: e(31*37)
            -zeta4
            sage: parent(e(31*37))
            Cyclotomic Field of order 4 and degree 2
        """
        N = self.modulus()
        m = m % N
        if self.values.is_in_cache() or m != N - 1:
            return self.values()[m]
        else:
            return self.__eval_at_minus_one()

    def change_ring(self, R):
        """
        Return the base extension of ``self`` to ``R``.

        INPUT:

        - ``R`` -- either a ring admitting a conversion map from the
          base ring of ``self``, or a ring homomorphism with the base
          ring of ``self`` as its domain

        EXAMPLE::

            sage: e = DirichletGroup(7, QQ).0
            sage: f = e.change_ring(QuadraticField(3, 'a'))
            sage: f.parent()
            Group of Dirichlet characters modulo 7 with values in Number Field in a with defining polynomial x^2 - 3

        ::

            sage: e = DirichletGroup(13).0
            sage: e.change_ring(QQ)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce zeta12 to a rational

        We test the case where `R` is a map (:trac:`18072`)::

            sage: K.<i> = QuadraticField(-1)
            sage: chi = DirichletGroup(5, K)[1]
            sage: chi(2)
            i
            sage: f = K.complex_embeddings()[0]
            sage: psi = chi.change_ring(f)
            sage: psi(2)
            -1.83697019872103e-16 - 1.00000000000000*I

        """
        if self.base_ring() is R:
            return self
        G = self.parent().change_ring(R)
        return G.element_class(G, [R(x) for x in self.values_on_gens()])

    def _cmp_(self, other):
        """
        Compare ``self`` to ``other``.

        .. NOTE::

            Since there is no coercion between Dirichlet groups
            of different moduli, characters of different moduli
            compare as unequal, even if they define identical
            functions on ``ZZ``.

        EXAMPLES::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: f = e.restrict(8)
            sage: e == e
            True
            sage: f == f
            True
            sage: e == f
            False
            sage: k = DirichletGroup(7)([-1])
            sage: k == e
            False

        """
        return cmp(self.values_on_gens(), other.values_on_gens())

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: hash(e)
            -1497246046           # 32-bit
            -3713082714463545694  # 64-bit
        """
        return hash(self.values_on_gens())

    def __invert__(self):
        """
        Return the multiplicative inverse of self.

        EXAMPLES::

            sage: e = DirichletGroup(13).0
            sage: f = ~e
            sage: f*e
            Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1
        """
        G = self.parent()
        if G.zeta.is_in_cache():
            x = -self.element()
        else:
            x = tuple(~z for z in self.values_on_gens())
        return G.element_class(G, x, check=False)

    def _mul_(self,  other):
        """
        Return the product of self and other.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: b
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4
            sage: a*b # indirect doctest
            Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4

        Multiplying elements whose parents have different zeta orders works::

            sage: a = DirichletGroup(3, QQ, zeta=1, zeta_order=1)(1)
            sage: b = DirichletGroup(3, QQ, zeta=-1, zeta_order=2)([-1])
            sage: a * b # indirect doctest
            Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1
        """
        G = self.parent()
        if G.zeta.is_in_cache():
            x = self.element() + other.element()
        else:
            x = tuple(y * z for y, z in zip(self.values_on_gens(), other.values_on_gens()))
        return G.element_class(G, x, check=False)

    def __copy__(self):
        """
        Return a (shallow) copy of this Dirichlet character.

        EXAMPLE::

            sage: G.<a> = DirichletGroup(11)
            sage: b = copy(a)
            sage: a is b
            False
            sage: a.element() is b.element()
            False
            sage: a.values_on_gens() is b.values_on_gens()
            True
        """
        # This method exists solely because of a bug in the cPickle module --
        # see modsym/manin_symbols.py.
        G = self.parent()
        return G.element_class(G, self.values_on_gens(), check=False)

    def __pow__(self, n):
        """
        Return self raised to the power of n

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a^2
            Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1
            sage: b^2
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -1
        """
        G = self.parent()
        if G.zeta.is_in_cache():
            x = n * self.element()
        else:
            x = tuple(z**n for z in self.values_on_gens())
        return G.element_class(G, x, check=False)

    def _repr_short_(self):
        r"""
        A short string representation of self, often used in string representations of modular forms

        EXAMPLES::

            sage: chi = DirichletGroup(24).0
            sage: chi._repr_short_()
            '[-1, 1, 1]'

        """
        return str(list(self.values_on_gens()))

    def _repr_(self):
        """
        String representation of self.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: repr(a) # indirect doctest
            'Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1'

        TESTS:

        Dirichlet characters modulo 1 and 2 are printed correctly (see
        :trac:`17338`)::

            sage: DirichletGroup(1)[0]
            Dirichlet character modulo 1 of conductor 1
            sage: DirichletGroup(2)[0]
            Dirichlet character modulo 2 of conductor 1
        """
        s = 'Dirichlet character modulo %s of conductor %s' % (self.modulus(), self.conductor())
        r = len(self.values_on_gens())
        if r != 0:
            s += ' mapping '
        for i in range(r):
            if i != 0:
                s += ', '
            s += str(self.parent().unit_gens()[i]) + ' |--> ' + str(self.values_on_gens()[i])
        return s

    def _latex_(self):
        r"""
        LaTeX representation of self.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(16)
            sage: latex(b)  # indirect doctest
            \hbox{Dirichlet character modulo } 16 \hbox{ of conductor } 16 \hbox{ mapping } 15 \mapsto 1,\ 5 \mapsto \zeta_{4}

        TESTS:

        Dirichlet characters modulo 1 and 2 are printed correctly (see
        :trac:`17338`)::

            sage: latex(DirichletGroup(1)[0])
            \hbox{Dirichlet character modulo } 1 \hbox{ of conductor } 1
            sage: latex(DirichletGroup(2)[0])
            \hbox{Dirichlet character modulo } 2 \hbox{ of conductor } 1
        """
        s = r'\hbox{Dirichlet character modulo } %s \hbox{ of conductor } %s' % (self.modulus(), self.conductor())
        r = len(self.values_on_gens())
        if r != 0:
            s += r' \hbox{ mapping } '
        for i in range(r):
            if i != 0:
                s += r',\ '
            s += self.parent().unit_gens()[i]._latex_() + r' \mapsto ' + self.values_on_gens()[i]._latex_()
        return s

    def base_ring(self):
        """
        Returns the base ring of this Dirichlet character.

        EXAMPLES::

            sage: G = DirichletGroup(11)
            sage: G.gen(0).base_ring()
            Cyclotomic Field of order 10 and degree 4
            sage: G = DirichletGroup(11, RationalField())
            sage: G.gen(0).base_ring()
            Rational Field
        """
        return self.parent().base_ring()

    def bar(self):
        """
        Return the complex conjugate of this Dirichlet character.

        EXAMPLES::

            sage: e = DirichletGroup(5).0
            sage: e
            Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4
            sage: e.bar()
            Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -zeta4
        """
        return ~self

    def bernoulli(self, k, algorithm='recurrence', cache=True, **opts):
        r"""
        Returns the generalized Bernoulli number `B_{k,eps}`.

        INPUT:

        - ``k`` -- a non-negative integer

        - ``algorithm`` -- either ``'recurrence'`` (default) or
          ``'definition'``

        - ``cache`` -- if True, cache answers

        - ``**opts`` -- optional arguments; not used directly, but
          passed to the :func:`bernoulli` function if this is called

        OUTPUT:

        Let `\varepsilon` be a (not necessarily primitive) character
        of modulus `N`.  This function returns the generalized
        Bernoulli number `B_{k,\varepsilon}`, as defined by the
        following identity of power series (see for example
        [Diamond-Im]_, Section 2.2):

        .. math::

            \sum_{a=1}^N \frac{\varepsilon(a) t e^{at}}{e^{Nt}-1}
            = sum_{k=0}^{\infty} \frac{B_{k,\varepsilon}}{k!} t^k.

        ALGORITHM:

        The ``'recurrence'`` algorithm computes generalized Bernoulli
        numbers via classical Bernoulli numbers using the formula in
        [Cohen-II]_, Proposition 9.4.5; this is usually optimal.  The
        ``definition`` algorithm uses the definition directly.

        .. WARNING::

            In the case of the trivial Dirichlet character modulo 1,
            this function returns `B_{1,\varepsilon} = 1/2`, in
            accordance with the above definition, but in contrast to
            the value `B_1 = -1/2` for the classical Bernoulli number.
            Some authors use an alternative definition giving
            `B_{1,\varepsilon} = -1/2`; see the discussion in
            [Cohen-II]_, Section 9.4.1.

        REFERENCES:

        .. [Cohen-II] \H. Cohen, Number Theory and Diophantine
           Equations, Volume II.  Graduate Texts in Mathematics 240.
           Springer, 2007.

        .. [Diamond-Im] \F. Diamond and J. Im, Modular forms and
           modular curves.  In: V. Kumar Murty (ed.), Seminar on
           Fermat's Last Theorem (Toronto, 1993-1994), 39-133.  CMS
           Conference Proceedings 17.  American Mathematical Society,
           1995.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.bernoulli(5)
            7430/13*zeta12^3 - 34750/13*zeta12^2 - 11380/13*zeta12 + 9110/13
            sage: eps = DirichletGroup(9).0
            sage: eps.bernoulli(3)
            10*zeta6 + 4
            sage: eps.bernoulli(3, algorithm="definition")
            10*zeta6 + 4

        TESTS:

        Check that :trac:`17586` is fixed::

            sage: DirichletGroup(1)[0].bernoulli(1)
            1/2

        """
        if cache:
            try:
                self.__bernoulli
            except AttributeError:
                self.__bernoulli = {}
            if k in self.__bernoulli:
                return self.__bernoulli[k]
        N = self.modulus()
        K = self.base_ring()

        if N == 1:
            # By definition, the first Bernoulli number of the trivial
            # character is 1/2, in contrast to the value B_1 = -1/2.
            ber = K.one()/2 if k == 1 else K(bernoulli(k))
        elif self(-1) != K((-1)**k):
            ber = K.zero()
        elif algorithm == "recurrence":
            # The following code is pretty fast, at least compared to
            # the other algorithm below.  That said, I'm sure it could
            # be sped up by a factor of 10 or more in many cases,
            # especially since we end up computing all the bernoulli
            # numbers up to k, which should be done with power series
            # instead of calls to the bernoulli function.  Likewise
            # computing all binomial coefficients can be done much
            # more efficiently.
            v = self.values()
            S = lambda n: sum(v[r] * r**n for r in range(1, N))
            ber = K(sum(binomial(k,j) * bernoulli(j, **opts) *
                        N**(j-1) * S(k-j) for j in range(k+1)))
        elif algorithm == "definition":
            # This is better since it computes the same thing, but requires
            # no arith in a poly ring over a number field.
            prec = k+2
            R = rings.PowerSeriesRing(rings.QQ, 't')
            t = R.gen()
            # g(t) = t/(e^{Nt}-1)
            g = t/((N*t).exp(prec) - 1)
            # h(n) = g(t)*e^{nt}
            h = [0] + [g * ((n*t).exp(prec)) for n in range(1,N+1)]
            ber = sum([self(a)*h[a][k] for a in range(1,N+1)]) * factorial(k)
        else:
            raise ValueError("algorithm = '%s' unknown"%algorithm)

        if cache:
            self.__bernoulli[k] = ber
        return ber

    @cached_method
    def conductor(self):
        """
        Computes and returns the conductor of this character.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.conductor()
            4
            sage: b.conductor()
            5
            sage: (a*b).conductor()
            20

        TESTS::

            sage: G.<a, b> = DirichletGroup(20)
            sage: type(G(1).conductor())
            <type 'sage.rings.integer.Integer'>
        """
        if self.modulus() == 1 or self.is_trivial():
            return rings.Integer(1)
        F = factor(self.modulus())
        if len(F) > 1:
            return prod([d.conductor() for d in self.decomposition()])
        p = F[0][0]
        # When p is odd, and x =/= 1, the conductor is the smallest p**r such that
        #   Order(x) divides EulerPhi(p**r) = p**(r-1)*(p-1).
        # For a given r, whether or not the above divisibility holds
        # depends only on the factor of p**(r-1) on the right hand side.
        # Since p-1 is coprime to p, this smallest r such that the
        # divisibility holds equals Valuation(Order(x),p)+1.
        cond = p**(valuation(self.order(),p) + 1)
        if p == 2 and F[0][1] > 2 and self.values_on_gens()[1].multiplicative_order() != 1:
            cond *= 2
        return rings.Integer(cond)

    @cached_method
    def decomposition(self):
        """
        Return the decomposition of self as a product of Dirichlet
        characters of prime power modulus, where the prime powers exactly
        divide the modulus of this character.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: c = a*b
            sage: d = c.decomposition(); d
            [Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4]
            sage: d[0].parent()
            Group of Dirichlet characters modulo 4 with values in Cyclotomic Field of order 4 and degree 2
            sage: d[1].parent()
            Group of Dirichlet characters modulo 5 with values in Cyclotomic Field of order 4 and degree 2

        We can't multiply directly, since coercion of one element into the
        other parent fails in both cases::

            sage: d[0]*d[1] == c
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Group of Dirichlet characters modulo 4 with values in Cyclotomic Field of order 4 and degree 2' and 'Group of Dirichlet characters modulo 5 with values in Cyclotomic Field of order 4 and degree 2'

        We can multiply if we're explicit about where we want the
        multiplication to take place.

        ::

            sage: G(d[0])*G(d[1]) == c
            True

        Conductors that are divisible by various powers of 2 present
        some problems as the multiplicative group modulo `2^k` is
        trivial for `k = 1` and non-cyclic for `k \ge 3`::

            sage: (DirichletGroup(18).0).decomposition()
            [Dirichlet character modulo 2 of conductor 1, Dirichlet character modulo 9 of conductor 9 mapping 2 |--> zeta6]
            sage: (DirichletGroup(36).0).decomposition()
            [Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1, Dirichlet character modulo 9 of conductor 1 mapping 2 |--> 1]
            sage: (DirichletGroup(72).0).decomposition()
            [Dirichlet character modulo 8 of conductor 4 mapping 7 |--> -1, 5 |--> 1, Dirichlet character modulo 9 of conductor 1 mapping 2 |--> 1]
        """
        D = self.parent().decomposition()
        vals = [[z] for z in self.values_on_gens()]
        if self.modulus() % 8 == 0:   # 2 factors at 2.
            vals[0].append(vals[1][0])
            del vals[1]
        elif self.modulus() % 4 == 2: # 0 factors at 2.
            vals = [1] + vals
        return [D[i](vals[i]) for i in range(len(D))]

    def extend(self, M):
        """
        Returns the extension of this character to a Dirichlet character
        modulo the multiple M of the modulus.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: H.<c> = DirichletGroup(4)
            sage: c.extend(20)
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: a
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: c.extend(20) == a
            True
        """
        if M % self.modulus() != 0:
            raise ArithmeticError("M(=%s) must be a multiple of the modulus(=%s)"%(M,self.modulus()))
        H = DirichletGroup(M, self.base_ring())
        return H(self)

    def galois_orbit(self, sort=True):
        r"""
        Return the orbit of this character under the action of the absolute
        Galois group of the prime subfield of the base ring.

        EXAMPLES::

            sage: G = DirichletGroup(30); e = G.1
            sage: e.galois_orbit()
            [Dirichlet character modulo 30 of conductor 5 mapping 11 |--> 1, 7 |--> -zeta4,
             Dirichlet character modulo 30 of conductor 5 mapping 11 |--> 1, 7 |--> zeta4]

        Another example::

            sage: G = DirichletGroup(13)
            sage: G.galois_orbits()
            [
            [Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1],
            ...,
            [Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -1]
            ]
            sage: e = G.0
            sage: e
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12
            sage: e.galois_orbit()
            [Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12,
             Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta12^3 + zeta12,
             Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12^3 - zeta12,
             Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta12]
            sage: e = G.0^2; e
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12^2
            sage: e.galois_orbit()
            [Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12^2, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta12^2 + 1]

        A non-example::

            sage: chi = DirichletGroup(7, Integers(9), zeta = Integers(9)(2)).0
            sage: chi.galois_orbit()
            Traceback (most recent call last):
            ...
            TypeError: Galois orbits only defined if base ring is an integral domain
        """
        if not self.base_ring().is_integral_domain():
            raise TypeError("Galois orbits only defined if base ring is an integral domain")
        k = self.order()
        if k <= 2:
            return [self]
        P = self.parent()
        z = self.element()
        o = int(z.additive_order())
        Auts = set([m % o for m in P._automorphisms()])
        v = [P.element_class(P, m * z, check=False) for m in Auts]
        if sort:
            v.sort()
        return v

    def gauss_sum(self, a=1):
        r"""
        Return a Gauss sum associated to this Dirichlet character.

        The Gauss sum associated to `\chi` is

        .. math::

            g_a(\chi) = \sum_{r \in \ZZ/m\ZZ} \chi(r)\,\zeta^{ar},

        where `m` is the modulus of `\chi` and `\zeta` is a primitive
        `m^{th}` root of unity.

        FACTS: If the modulus is a prime `p` and the character is
        nontrivial, then the Gauss sum has absolute value `\sqrt{p}`.

        CACHING: Computed Gauss sums are *not* cached with this character.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G([-1])
            sage: e.gauss_sum(1)
            2*zeta6 - 1
            sage: e.gauss_sum(2)
            -2*zeta6 + 1
            sage: norm(e.gauss_sum())
            3

        ::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.gauss_sum()
            -zeta156^46 + zeta156^45 + zeta156^42 + zeta156^41 + 2*zeta156^40 + zeta156^37 - zeta156^36 - zeta156^34 - zeta156^33 - zeta156^31 + 2*zeta156^30 + zeta156^28 - zeta156^24 - zeta156^22 + zeta156^21 + zeta156^20 - zeta156^19 + zeta156^18 - zeta156^16 - zeta156^15 - 2*zeta156^14 - zeta156^10 + zeta156^8 + zeta156^7 + zeta156^6 + zeta156^5 - zeta156^4 - zeta156^2 - 1
            sage: factor(norm(e.gauss_sum()))
            13^24

        TESTS:

        The field of algebraic numbers is supported (:trac:`19056`)::

            sage: G = DirichletGroup(7, QQbar)
            sage: G[1].gauss_sum()
            -2.440133358345538? + 1.022618791871794?*I

        Check that :trac:`19060` is fixed::

            sage: K.<z> = CyclotomicField(8)
            sage: G = DirichletGroup(13, K)
            sage: chi = G([z^2])
            sage: chi.gauss_sum()
            zeta52^22 + zeta52^21 + zeta52^19 - zeta52^16 + zeta52^15 + zeta52^14 + zeta52^12 - zeta52^11 - zeta52^10 - zeta52^7 - zeta52^5 + zeta52^4
        """
        G = self.parent()
        K = G.base_ring()
        chi = self
        m = G.modulus()
        if is_ComplexField(K):
            return self.gauss_sum_numerical(a=a)
        elif is_AlgebraicField(K):
            L = K
            zeta = L.zeta(m)
        elif number_field.is_CyclotomicField(K) or is_RationalField(K):
            chi = chi.minimize_base_ring()
            n = lcm(m, G.zeta_order())
            L = rings.CyclotomicField(n)
            zeta = L.gen(0) ** (n // m)
        else:
            raise NotImplementedError("Gauss sums only currently implemented when the base ring is a cyclotomic field, QQ, QQbar, or a complex field")
        zeta = zeta ** a
        g = L.zero()
        z = L.one()
        for c in chi.values()[1:]:
            z *= zeta
            g += L(c)*z
        return g

    def gauss_sum_numerical(self, prec=53, a=1):
        r"""
        Return a Gauss sum associated to this Dirichlet character as an
        approximate complex number with prec bits of precision.

        INPUT:

        - ``prec`` -- integer (default: 53), *bits* of precision

        - ``a`` -- integer, as for :meth:`gauss_sum`.

        The Gauss sum associated to `\chi` is

        .. math::

            g_a(\chi) = \sum_{r \in \ZZ/m\ZZ} \chi(r)\,\zeta^{ar},

        where `m` is the modulus of `\chi` and `\zeta` is a primitive
        `m^{th}` root of unity.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G.0
            sage: abs(e.gauss_sum_numerical())
            1.7320508075...
            sage: sqrt(3.0)
            1.73205080756888
            sage: e.gauss_sum_numerical(a=2)
            -...e-15 - 1.7320508075...*I
            sage: e.gauss_sum_numerical(a=2, prec=100)
            4.7331654313260708324703713917e-30 - 1.7320508075688772935274463415*I
            sage: G = DirichletGroup(13)
            sage: H = DirichletGroup(13, CC)
            sage: e = G.0
            sage: f = H.0
            sage: e.gauss_sum_numerical()
            -3.07497205... + 1.8826966926...*I
            sage: f.gauss_sum_numerical()
            -3.07497205... + 1.8826966926...*I
            sage: abs(e.gauss_sum_numerical())
            3.60555127546...
            sage: abs(f.gauss_sum_numerical())
            3.60555127546...
            sage: sqrt(13.0)
            3.60555127546399

        TESTS:

        The field of algebraic numbers is supported (:trac:`19056`)::

            sage: G = DirichletGroup(7, QQbar)
            sage: G[1].gauss_sum_numerical()
            -2.44013335834554 + 1.02261879187179*I
        """
        G = self.parent()
        K = G.base_ring()
        if is_ComplexField(K):
            phi = lambda t : t
            CC = K
        elif is_AlgebraicField(K):
            from sage.rings.complex_field import ComplexField
            CC = ComplexField(prec)
            phi = CC.coerce_map_from(K)
        elif number_field.is_CyclotomicField(K) or is_RationalField(K):
            phi = K.complex_embedding(prec)
            CC = phi.codomain()
        else:
            raise NotImplementedError("Gauss sums only currently implemented when the base ring is a cyclotomic field, QQ, QQbar, or a complex field")
        zeta = CC.zeta(G.modulus()) ** a
        g = CC.zero()
        z = CC.one()
        for c in self.values()[1:]:
            z *= zeta
            g += phi(c)*z
        return g

    def jacobi_sum(self, char, check=True):
        """
        Return the Jacobi sum associated to these Dirichlet characters
        (i.e., J(self,char)). This is defined as

        .. math::

            J(\chi, \psi) = \sum_{a \in \ZZ / N\ZZ} \chi(a) \psi(1-a)

        where `\chi` and `\psi` are both characters modulo `N`.

        EXAMPLES::

            sage: D = DirichletGroup(13)
            sage: e = D.0
            sage: f = D[-2]
            sage: e.jacobi_sum(f)
            3*zeta12^2 + 2*zeta12 - 3
            sage: f.jacobi_sum(e)
            3*zeta12^2 + 2*zeta12 - 3
            sage: p = 7
            sage: DP = DirichletGroup(p)
            sage: f = DP.0
            sage: e.jacobi_sum(f)
            Traceback (most recent call last):
            ...
            NotImplementedError: Characters must be from the same Dirichlet Group.

            sage: all_jacobi_sums = [(DP[i].values_on_gens(),DP[j].values_on_gens(),DP[i].jacobi_sum(DP[j])) \
            ...                       for i in range(p-1) for j in range(p-1)[i:]]
            ...
            sage: for s in all_jacobi_sums:
            ...       print s
            ((1,), (1,), 5)
            ((1,), (zeta6,), -1)
            ((1,), (zeta6 - 1,), -1)
            ((1,), (-1,), -1)
            ((1,), (-zeta6,), -1)
            ((1,), (-zeta6 + 1,), -1)
            ((zeta6,), (zeta6,), -zeta6 + 3)
            ((zeta6,), (zeta6 - 1,), 2*zeta6 + 1)
            ((zeta6,), (-1,), -2*zeta6 - 1)
            ((zeta6,), (-zeta6,), zeta6 - 3)
            ((zeta6,), (-zeta6 + 1,), 1)
            ((zeta6 - 1,), (zeta6 - 1,), -3*zeta6 + 2)
            ((zeta6 - 1,), (-1,), 2*zeta6 + 1)
            ((zeta6 - 1,), (-zeta6,), -1)
            ((zeta6 - 1,), (-zeta6 + 1,), -zeta6 - 2)
            ((-1,), (-1,), 1)
            ((-1,), (-zeta6,), -2*zeta6 + 3)
            ((-1,), (-zeta6 + 1,), 2*zeta6 - 3)
            ((-zeta6,), (-zeta6,), 3*zeta6 - 1)
            ((-zeta6,), (-zeta6 + 1,), -2*zeta6 + 3)
            ((-zeta6 + 1,), (-zeta6 + 1,), zeta6 + 2)

        Let's check that trivial sums are being calculated correctly::

            sage: N = 13
            sage: D = DirichletGroup(N)
            sage: g = D(1)
            sage: g.jacobi_sum(g)
            11
            sage: sum([g(x)*g(1-x) for x in IntegerModRing(N)])
            11

        And sums where exactly one character is nontrivial (see :trac:`6393`)::

            sage: G = DirichletGroup(5); X=G.list(); Y=X[0]; Z=X[1]
            sage: Y.jacobi_sum(Z)
            -1
            sage: Z.jacobi_sum(Y)
            -1

        Now let's take a look at a non-prime modulus::

            sage: N = 9
            sage: D = DirichletGroup(N)
            sage: g = D(1)
            sage: g.jacobi_sum(g)
            3

        We consider a sum with values in a finite field::

            sage: g = DirichletGroup(17, GF(9,'a')).0
            sage: g.jacobi_sum(g**2)
            2*a

        TESTS:

        This shows that :trac:`6393` has been fixed::

            sage: G = DirichletGroup(5); X = G.list(); Y = X[0]; Z = X[1]
            sage: # Y is trivial and Z is quartic
            sage: sum([Y(x)*Z(1-x) for x in IntegerModRing(5)])
            -1
            sage: # The value -1 above is the correct value of the Jacobi sum J(Y, Z).
            sage: Y.jacobi_sum(Z); Z.jacobi_sum(Y)
            -1
            -1
        """
        if check:
            if self.parent() != char.parent():
                raise NotImplementedError("Characters must be from the same Dirichlet Group.")

        return sum([self(x) * char(1-x) for x in rings.IntegerModRing(self.modulus())])

    def kloosterman_sum(self, a=1,b=0):
        r"""
        Return the "twisted" Kloosterman sum associated to this Dirichlet character.
        This includes Gauss sums, classical Kloosterman sums, Salie sums, etc.

        The Kloosterman sum associated to `\chi` and the integers a,b is

        .. math::

            K(a,b,\chi) = \sum_{r \in (\ZZ/m\ZZ)^\times} \chi(r)\,\zeta^{ar+br^{-1}},

        where `m` is the modulus of `\chi` and `\zeta` is a primitive
        `m` th root of unity. This reduces to to the Gauss sum if `b=0`.

        This method performs an exact calculation and returns an element of a
        suitable cyclotomic field; see also :meth:`.kloosterman_sum_numerical`,
        which gives an inexact answer (but is generally much quicker).

        CACHING: Computed Kloosterman sums are *not* cached with this
        character.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G([-1])
            sage: e.kloosterman_sum(3,5)
            -2*zeta6 + 1
            sage: G = DirichletGroup(20)
            sage: e = G([1 for  u in G.unit_gens()])
            sage: e.kloosterman_sum(7,17)
            -2*zeta20^6 + 2*zeta20^4 + 4

        """
        G = self.parent()
        K = G.base_ring()
        if not (number_field.is_CyclotomicField(K) or is_RationalField(K)):
            raise NotImplementedError("Kloosterman sums only currently implemented when the base ring is a cyclotomic field or QQ.")
        g = 0
        m = G.modulus()
        L = rings.CyclotomicField(lcm(m,G.zeta_order()))
        zeta = L.gen(0)
        n = zeta.multiplicative_order()
        zeta = zeta ** (n // m)
        for c in range(1,m):
            if gcd(c,m)==1:
                e = rings.Mod(c,m)
                z = zeta ** int(a*e + b*(e**(-1)))
                g += self(c)*z
        return g

    def kloosterman_sum_numerical(self, prec=53, a=1,b=0):
        r"""
        Return the Kloosterman sum associated to this Dirichlet character as
        an approximate complex number with prec bits of precision. See also
        :meth:`.kloosterman_sum`, which calculates the sum exactly (which is
        generally slower).

        INPUT:

        - ``prec`` -- integer (default: 53), *bits* of precision
        - ``a`` -- integer, as for :meth:`.kloosterman_sum`
        - ``b`` -- integer, as for :meth:`.kloosterman_sum`.

        EXAMPLES::

            sage: G = DirichletGroup(3)
            sage: e = G.0

        The real component of the numerical value of e is near zero::

            sage: v=e.kloosterman_sum_numerical()
            sage: v.real() < 1.0e15
            True
            sage: v.imag()
            1.73205080756888
            sage: G = DirichletGroup(20)
            sage: e = G.1
            sage: e.kloosterman_sum_numerical(53,3,11)
            3.80422606518061 - 3.80422606518061*I
        """
        G = self.parent()
        K = G.base_ring()
        if not (number_field.is_CyclotomicField(K) or is_RationalField(K)):
            raise NotImplementedError("Kloosterman sums only currently implemented when the base ring is a cyclotomic field or QQ.")
        phi = K.complex_embedding(prec)
        CC = phi.codomain()
        g = 0
        m = G.modulus()
        zeta = CC.zeta(m)

        for c in range(1,m):
            if gcd(c,m)==1:
                e = rings.Mod(c,m)
                z = zeta ** int(a*e + b*(e**(-1)))
                g += phi(self(c))*z
        return g

    @cached_method
    def is_even(self):
        r"""
        Return ``True`` if and only if `\varepsilon(-1) = 1`.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.is_even()
            False
            sage: e(-1)
            -1
            sage: [e.is_even() for e in G]
            [True, False, True, False, True, False, True, False, True, False, True, False]

            sage: G = DirichletGroup(13, CC)
            sage: e = G.0
            sage: e.is_even()
            False
            sage: e(-1)
            -1.000000...
            sage: [e.is_even() for e in G]
            [True, False, True, False, True, False, True, False, True, False, True, False]

            sage: G = DirichletGroup(100000, CC)
            sage: G.1.is_even()
            True

        Note that ``is_even`` need not be the negation of
        is_odd, e.g., in characteristic 2::

            sage: G.<e> = DirichletGroup(13, GF(4,'a'))
            sage: e.is_even()
            True
            sage: e.is_odd()
            True
        """
        R = self.base_ring()
        # self(-1) is either +1 or -1
        if not R.is_exact():
            return abs(self(-1) - R(1)) < 0.5
        return self(-1) == R(1)

    @cached_method
    def is_odd(self):
        r"""
        Return ``True`` if and only if
        `\varepsilon(-1) = -1`.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.is_odd()
            True
            sage: [e.is_odd() for e in G]
            [False, True, False, True, False, True, False, True, False, True, False, True]

            sage: G = DirichletGroup(13)
            sage: e = G.0
            sage: e.is_odd()
            True
            sage: [e.is_odd() for e in G]
            [False, True, False, True, False, True, False, True, False, True, False, True]

            sage: G = DirichletGroup(100000, CC)
            sage: G.0.is_odd()
            True

        Note that ``is_even`` need not be the negation of
        is_odd, e.g., in characteristic 2::

            sage: G.<e> = DirichletGroup(13, GF(4,'a'))
            sage: e.is_even()
            True
            sage: e.is_odd()
            True
        """
        R = self.base_ring()
        # self(-1) is either +1 or -1
        if not R.is_exact():
            return abs(self(-1) - R(-1)) < 0.5
        return self(-1) == R(-1)

    @cached_method
    def is_primitive(self):
        """
        Return ``True`` if and only if this character is
        primitive, i.e., its conductor equals its modulus.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.is_primitive()
            False
            sage: b.is_primitive()
            False
            sage: (a*b).is_primitive()
            True
            sage: G.<a,b> = DirichletGroup(20, CC)
            sage: a.is_primitive()
            False
            sage: b.is_primitive()
            False
            sage: (a*b).is_primitive()
            True
        """
        return (self.conductor() == self.modulus())

    @cached_method
    def is_trivial(self):
        r"""
        Returns ``True`` if this is the trivial character,
        i.e., has order 1.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.is_trivial()
            False
            sage: (a^2).is_trivial()
            True
        """
        if self.element.is_in_cache():
            return not self.element()
        one = self.base_ring().one()
        return all(x == one for x in self.values_on_gens())

    def kernel(self):
        r"""
        Return the kernel of this character.

        OUTPUT: Currently the kernel is returned as a list. This may
        change.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.kernel()
            [1, 9, 13, 17]
            sage: b.kernel()
            [1, 11]
        """
        one = self.base_ring().one()
        return [x for x in range(self.modulus()) if self(x) == one]

    def maximize_base_ring(self):
        r"""
        Let

        .. math::

           \varepsilon : (\ZZ/N\ZZ)^* \to \QQ(\zeta_n)

        be a Dirichlet character. This function returns an equal Dirichlet
        character

        .. math::

           \chi : (\ZZ/N\ZZ)^* \to \QQ(\zeta_m)

        where `m` is the least common multiple of `n` and
        the exponent of `(\ZZ/N\ZZ)^*`.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20,QQ)
            sage: b.maximize_base_ring()
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -1
            sage: b.maximize_base_ring().base_ring()
            Cyclotomic Field of order 4 and degree 2
            sage: DirichletGroup(20).base_ring()
            Cyclotomic Field of order 4 and degree 2
        """
        g = rings.IntegerModRing(self.modulus()).unit_group_exponent()
        if g == 1:
            g = 2
        z = self.base_ring().zeta()
        n = z.multiplicative_order()
        m = lcm(g,n)
        if n == m:
            return self
        K = rings.CyclotomicField(m)
        return self.change_ring(K)

    def minimize_base_ring(self):
        r"""
        Return a Dirichlet character that equals this one, but over as
        small a subfield (or subring) of the base ring as possible.

        .. note::

           This function is currently only implemented when the base
           ring is a number field. It's the identity function in
           characteristic p.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: e = DirichletGroup(13).0
            sage: e.base_ring()
            Cyclotomic Field of order 12 and degree 4
            sage: e.minimize_base_ring().base_ring()
            Cyclotomic Field of order 12 and degree 4
            sage: (e^2).minimize_base_ring().base_ring()
            Cyclotomic Field of order 6 and degree 2
            sage: (e^3).minimize_base_ring().base_ring()
            Cyclotomic Field of order 4 and degree 2
            sage: (e^12).minimize_base_ring().base_ring()
            Rational Field

        TESTS:

        Check that :trac:`18479` is fixed::

            sage: f = Newforms(Gamma1(25), names='a')[1]
            sage: eps = f.character()
            sage: eps.minimize_base_ring() == eps
            True
        
        A related bug (see :trac:`18086`)::

            sage: K.<a,b>=NumberField([x^2 + 1, x^2 - 3])
            sage: chi = DirichletGroup(7, K).0
            sage: chi.minimize_base_ring()
            Dirichlet character modulo 7 of conductor 7 mapping 3 |--> -1/2*b*a + 1/2
        """
        R = self.base_ring()
        if R.is_prime_field():
            return self
        p = R.characteristic()

        if p:
            K = rings.IntegerModRing(p)
        elif self.order() <= 2:
            K = rings.QQ
        elif (isinstance(R, number_field.NumberField_generic)
              and euler_phi(self.order()) < R.absolute_degree()):
            K = rings.CyclotomicField(self.order())
        else:
            return self

        try:
            return self.change_ring(K)
        except (TypeError, ValueError, ArithmeticError):
            return self

    def modulus(self):
        """
        The modulus of this character.

        EXAMPLES::

            sage: e = DirichletGroup(100, QQ).0
            sage: e.modulus()
            100
            sage: e.conductor()
            4
        """
        return self.parent().modulus()

    def level(self):
        """
        Synonym for modulus.

        EXAMPLE::

            sage: e = DirichletGroup(100, QQ).0
            sage: e.level()
            100
        """
        return self.modulus()

    @cached_method
    def multiplicative_order(self):
        """
        The order of this character.

        EXAMPLES::

            sage: e = DirichletGroup(100).1
            sage: e.order()    # same as multiplicative_order, since group is multiplicative
            20
            sage: e.multiplicative_order()
            20
            sage: e = DirichletGroup(100).0
            sage: e.multiplicative_order()
            2
        """
        if self.parent().zeta.is_in_cache():
            return self.element().additive_order()
        return lcm([z.multiplicative_order() for z in self.values_on_gens()])

    def primitive_character(self):
        """
        Returns the primitive character associated to self.

        EXAMPLES::

            sage: e = DirichletGroup(100).0; e
            Dirichlet character modulo 100 of conductor 4 mapping 51 |--> -1, 77 |--> 1
            sage: e.conductor()
            4
            sage: f = e.primitive_character(); f
            Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1
            sage: f.modulus()
            4
        """
        return self.restrict(self.conductor())

    def restrict(self, M):
        """
        Returns the restriction of this character to a Dirichlet character
        modulo the divisor M of the modulus, which must also be a multiple
        of the conductor of this character.

        EXAMPLES::

            sage: e = DirichletGroup(100).0
            sage: e.modulus()
            100
            sage: e.conductor()
            4
            sage: e.restrict(20)
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: e.restrict(4)
            Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1
            sage: e.restrict(50)
            Traceback (most recent call last):
            ...
            ValueError: conductor(=4) must divide M(=50)
        """
        M = int(M)
        if self.modulus()%M != 0:
            raise ValueError("M(=%s) must divide the modulus(=%s)"%(M,self.modulus()))
        if M%self.conductor() != 0:
            raise ValueError("conductor(=%s) must divide M(=%s)"%(self.conductor(),M))
        H = DirichletGroup(M, self.base_ring())
        return H(self)

    @cached_method
    def values(self):
        """
        Return a list of the values of this character on each integer
        between 0 and the modulus.

        EXAMPLES::

            sage: e = DirichletGroup(20)(1)
            sage: e.values()
            [0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1]
            sage: e = DirichletGroup(20).gen(0)
            sage: print e.values()
            [0, 1, 0, -1, 0, 0, 0, -1, 0, 1, 0, -1, 0, 1, 0, 0, 0, 1, 0, -1]
            sage: e = DirichletGroup(20).gen(1)
            sage: e.values()
            [0, 1, 0, -zeta4, 0, 0, 0, zeta4, 0, -1, 0, 1, 0, -zeta4, 0, 0, 0, zeta4, 0, -1]
            sage: e = DirichletGroup(21).gen(0) ; e.values()
            [0, 1, -1, 0, 1, -1, 0, 0, -1, 0, 1, -1, 0, 1, 0, 0, 1, -1, 0, 1, -1]
            sage: e = DirichletGroup(21, base_ring=GF(37)).gen(0) ; e.values()
            [0, 1, 36, 0, 1, 36, 0, 0, 36, 0, 1, 36, 0, 1, 0, 0, 1, 36, 0, 1, 36]
            sage: e = DirichletGroup(21, base_ring=GF(3)).gen(0) ; e.values()
            [0, 1, 2, 0, 1, 2, 0, 0, 2, 0, 1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 2]

        ::

            sage: chi = DirichletGroup(100151, CyclotomicField(10)).0
            sage: ls = chi.values() ; ls[0:10]
            [0,
            1,
            -zeta10^3,
            -zeta10,
            -zeta10,
            1,
            zeta10^3 - zeta10^2 + zeta10 - 1,
            zeta10,
            zeta10^3 - zeta10^2 + zeta10 - 1,
            zeta10^2]

        TESTS:

        Test that :trac:`11783` and :trac:`14368` are fixed::

            sage: chi = DirichletGroup(1).list()[0]
            sage: chi.values()
            [1]
            sage: chi(1)
            1
        """
        G = self.parent()
        R = G.base_ring()

        mod = self.parent().modulus()
        if mod == 1:
            return [R.one()]
        elif mod == 2:
            return [R.zero(), R.one()]

        result_list = [R.zero()] * mod
        gens = G.unit_gens()
        orders = G.integers_mod().unit_group().gens_orders()

        R_values = G._zeta_powers
        val_on_gen = self.element()

        exponents = [0] * len(orders)
        n = G.integers_mod().one()
        value = val_on_gen.base_ring().zero()

        while True:
            # record character value on n
            result_list[n] = R_values[value]
            # iterate:
            #   increase the exponent vector by 1,
            #   increase n accordingly, and increase value
            i = 0
            while True:
                try:
                    exponents[i] += 1
                except IndexError:  # Done!
                    return result_list
                value += val_on_gen[i]
                n *= gens[i]
                if exponents[i] < orders[i]:
                    break
                exponents[i] = 0
                i += 1

    def values_on_gens(self):
        r"""
        Return a tuple of the values of ``self`` on the standard
        generators of `(\ZZ/N\ZZ)^*`, where `N` is the modulus.

        EXAMPLES::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: e.values_on_gens ()
            (-1, 1)
        """
        try:
            return self.__values_on_gens
        except AttributeError:
            pows = self.parent()._zeta_powers
            v = tuple([pows[i] for i in self.element()])
            self.__values_on_gens = v
            return v

    @cached_method
    def element(self):
        r"""
        Return the underlying `\ZZ/n\ZZ`-module
        vector of exponents.

        .. warning::

           Please do not change the entries of the returned vector;
           this vector is mutable *only* because immutable vectors are
           implemented yet.

        EXAMPLES::

            sage: G.<a,b> = DirichletGroup(20)
            sage: a.element()
            (2, 0)
            sage: b.element()
            (0, 1)
        """
        P = self.parent()
        M = P._module
        if is_ComplexField(P.base_ring()):
            zeta = P.zeta()
            zeta_argument = zeta.argument()
            v = M([int(round(x.argument()/zeta_argument))
                   for x in self.values_on_gens()])
        else:
            dlog = P._zeta_dlog
            v = M([dlog[x] for x in self.values_on_gens()])
        return v


class DirichletGroupFactory(UniqueFactory):
    r"""
    Construct a group of Dirichlet characters modulo `N`.

    INPUT:

    - ``N`` -- positive integer

    - ``base_ring`` -- commutative ring; the value ring for the
      characters in this group (default: the cyclotomic field
      `\QQ(\zeta_n)`, where `n` is the exponent of `(\ZZ/N\ZZ)^*`)

    - ``zeta`` -- (optional) root of unity in ``base_ring``

    - ``zeta_order`` -- (optional) positive integer; this must be the
      order of ``zeta`` if both are specified

    - ``names`` -- ignored (needed so ``G.<...> = DirichletGroup(...)``
      notation works)

    - ``integral`` -- boolean (default: ``False``); whether to replace
      the default cyclotomic field by its rings of integers as the
      base ring.  This is ignored if ``base_ring`` is not ``None``.

    OUTPUT:

    The group of Dirichlet characters modulo `N` with values in a
    subgroup `V` of the multiplicative group `R^*` of ``base_ring``.
    This is the group of homomorphisms `(\ZZ/N\ZZ)^* \to V` with
    pointwise multiplication.  The group `V` is determined as follows:

    - If both ``zeta`` and ``zeta_order`` are omitted, then `V` is
      taken to be `R^*`, or equivalently its `n`-torsion subgroup,
      where `n` is the exponent of `(\ZZ/N\ZZ)^*`.  Many operations,
      such as finding a set of generators for the group, are only
      implemented if `V` is cyclic and a generator for `V` can be
      found.

    - If ``zeta`` is specified, then `V` is taken to be the cyclic
      subgroup of `R^*` generated by ``zeta``.  If ``zeta_order`` is
      also given, it must be the multiplicative order of ``zeta``;
      this is useful if the base ring is not exact or if the order of
      ``zeta`` is very large.

    - If ``zeta`` is not specified but ``zeta_order`` is, then `V` is
      taken to be the group of roots of unity of order dividing
      ``zeta_order`` in `R`.  In this case, `R` must be a domain (so
      `V` is cyclic), and `V` must have order ``zeta_order``.
      Furthermore, a generator ``zeta`` of `V` is computed, and an
      error is raised if such ``zeta`` cannot be found.

    EXAMPLES:

    The default base ring is a cyclotomic field of order the exponent
    of `(\ZZ/N\ZZ)^*`::

        sage: DirichletGroup(20)
        Group of Dirichlet characters modulo 20 with values in Cyclotomic Field of order 4 and degree 2

    We create the group of Dirichlet character mod 20 with values in
    the rational numbers::

        sage: G = DirichletGroup(20, QQ); G
        Group of Dirichlet characters modulo 20 with values in Rational Field
        sage: G.order()
        4
        sage: G.base_ring()
        Rational Field

    The elements of G print as lists giving the values of the character
    on the generators of `(Z/NZ)^*`::

        sage: list(G)
        [Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1, Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1, Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -1, Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -1]

    Next we construct the group of Dirichlet character mod 20, but with
    values in `\QQ(\zeta_n)`::

        sage: G = DirichletGroup(20)
        sage: G.1
        Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4

    We next compute several invariants of ``G``::

        sage: G.gens()
        (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1, Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4)
        sage: G.unit_gens()
        (11, 17)
        sage: G.zeta()
        zeta4
        sage: G.zeta_order()
        4

    In this example we create a Dirichlet group with values in a
    number field::

        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(x^4 + 1)
        sage: DirichletGroup(5, K)
        Group of Dirichlet characters modulo 5 with values in Number Field in a with defining polynomial x^4 + 1

    An example where we give ``zeta``, but not its order::

        sage: G = DirichletGroup(5, K, a); G
        Group of Dirichlet characters modulo 5 with values in the group of order 8 generated by a in Number Field in a with defining polynomial x^4 + 1
        sage: G.list()
        [Dirichlet character modulo 5 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> a^2, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -a^2]

    We can also restrict the order of the characters, either with or
    without specifying a root of unity::

        sage: DirichletGroup(5, K, zeta=-1, zeta_order=2)
        Group of Dirichlet characters modulo 5 with values in the group of order 2 generated by -1 in Number Field in a with defining polynomial x^4 + 1
        sage: DirichletGroup(5, K, zeta_order=2)
        Group of Dirichlet characters modulo 5 with values in the group of order 2 generated by -1 in Number Field in a with defining polynomial x^4 + 1

    ::

        sage: G.<e> = DirichletGroup(13)
        sage: loads(G.dumps()) == G
        True

    ::

        sage: G = DirichletGroup(19, GF(5))
        sage: loads(G.dumps()) == G
        True

    We compute a Dirichlet group over a large prime field::

        sage: p = next_prime(10^40)
        sage: g = DirichletGroup(19, GF(p)); g
        Group of Dirichlet characters modulo 19 with values in Finite Field of size 10000000000000000000000000000000000000121

    Note that the root of unity has small order, i.e., it is not the
    largest order root of unity in the field::

        sage: g.zeta_order()
        2

    ::

        sage: r4 = CyclotomicField(4).ring_of_integers()
        sage: G = DirichletGroup(60, r4)
        sage: G.gens()
        (Dirichlet character modulo 60 of conductor 4 mapping 31 |--> -1, 41 |--> 1, 37 |--> 1, Dirichlet character modulo 60 of conductor 3 mapping 31 |--> 1, 41 |--> -1, 37 |--> 1, Dirichlet character modulo 60 of conductor 5 mapping 31 |--> 1, 41 |--> 1, 37 |--> zeta4)
        sage: val = G.gens()[2].values_on_gens()[2] ; val
        zeta4
        sage: parent(val)
        Gaussian Integers in Cyclotomic Field of order 4 and degree 2
        sage: r4.residue_field(r4.ideal(29).factor()[0][0])(val)
        17
        sage: r4.residue_field(r4.ideal(29).factor()[0][0])(val) * GF(29)(3)
        22
        sage: r4.residue_field(r4.ideal(29).factor()[0][0])(G.gens()[2].values_on_gens()[2]) * 3
        22
        sage: parent(r4.residue_field(r4.ideal(29).factor()[0][0])(G.gens()[2].values_on_gens()[2]) * 3)
        Residue field of Fractional ideal (-2*zeta4 + 5)

    ::

        sage: DirichletGroup(60, integral=True)
        Group of Dirichlet characters modulo 60 with values in Gaussian Integers in Cyclotomic Field of order 4 and degree 2
        sage: parent(DirichletGroup(60, integral=True).gens()[2].values_on_gens()[2])
        Gaussian Integers in Cyclotomic Field of order 4 and degree 2

    If the order of ``zeta`` cannot be determined automatically, we
    can specify it using ``zeta_order``::

        sage: DirichletGroup(7, CC, zeta=exp(2*pi*I/6))
        Traceback (most recent call last):
        ...
        NotImplementedError: order of element not known

        sage: DirichletGroup(7, CC, zeta=exp(2*pi*I/6), zeta_order=6)
        Group of Dirichlet characters modulo 7 with values in the group of order 6 generated by 0.500000000000000 + 0.866025403784439*I in Complex Field with 53 bits of precision

    If the base ring is not a domain (in which case the group of roots
    of unity is not necessarily cyclic), some operations still work,
    such as creation of elements::

        sage: G = DirichletGroup(5, Zmod(15)); G
        Group of Dirichlet characters modulo 5 with values in Ring of integers modulo 15
        sage: chi = G([13]); chi
        Dirichlet character modulo 5 of conductor 5 mapping 2 |--> 13
        sage: chi^2
        Dirichlet character modulo 5 of conductor 5 mapping 2 |--> 4
        sage: chi.multiplicative_order()
        4

    Other operations only work if ``zeta`` is specified::

        sage: G.gens()
        Traceback (most recent call last):
        ...
        NotImplementedError: factorization of polynomials over rings with composite characteristic is not implemented
        sage: G = DirichletGroup(5, Zmod(15), zeta=2); G
        Group of Dirichlet characters modulo 5 with values in the group of order 4 generated by 2 in Ring of integers modulo 15
        sage: G.gens()
        (Dirichlet character modulo 5 of conductor 5 mapping 2 |--> 2,)

    TESTS:

    Dirichlet groups are cached, creating two groups with the same parameters
    yields the same object::

        sage: DirichletGroup(60) is DirichletGroup(60)
        True

    """
    def create_key(self, N, base_ring=None, zeta=None, zeta_order=None, names=None, integral=False):
        """
        Create a key that uniquely determines a Dirichlet group.

        TESTS::

            sage: DirichletGroup.create_key(60)
            (Cyclotomic Field of order 4 and degree 2, 60, None, None)

        An example to illustrate that ``base_ring`` is a part of the key::

            sage: k = DirichletGroup.create_key(2, base_ring=QQ); k
            (Rational Field, 2, None, None)
            sage: l = DirichletGroup.create_key(2, base_ring=CC); l
            (Complex Field with 53 bits of precision, 2, None, None)
            sage: k == l
            False
            sage: G = DirichletGroup.create_object(None, k); G
            Group of Dirichlet characters modulo 2 with values in Rational Field
            sage: H = DirichletGroup.create_object(None, l); H
            Group of Dirichlet characters modulo 2 with values in Complex Field with 53 bits of precision
            sage: G == H
            False

        If ``base_ring`` was not be a part of the key, the keys would compare
        equal and the caching would be broken::

            sage: k = k[1:]; k
            (2, None, None)
            sage: l = l[1:]; l
            (2, None, None)
            sage: k == l
            True
            sage: DirichletGroup(2, base_ring=QQ) is DirichletGroup(2, base_ring=CC)
            False

        If the base ring is not an integral domain, an error will be
        raised if only ``zeta_order`` is specified::

            sage: DirichletGroup(17, Integers(15))
            Group of Dirichlet characters modulo 17 with values in Ring of integers modulo 15
            sage: DirichletGroup(17, Integers(15), zeta_order=4)
            Traceback (most recent call last):
            ...
            ValueError: base ring (= Ring of integers modulo 15) must be an integral domain if only zeta_order is specified
            sage: G = DirichletGroup(17, Integers(15), zeta=7); G
            Group of Dirichlet characters modulo 17 with values in the group of order 4 generated by 7 in Ring of integers modulo 15
            sage: G.order()
            4

        """
        modulus = rings.Integer(N)

        if base_ring is None:
            if not (zeta is None and zeta_order is None):
                raise ValueError("zeta and zeta_order must be None if base_ring not specified")
            e = rings.IntegerModRing(modulus).unit_group_exponent()
            base_ring = rings.CyclotomicField(e)
            if integral:
                base_ring = base_ring.ring_of_integers()

        if not is_Ring(base_ring):
            raise TypeError("base_ring (= %s) must be a ring" % base_ring)

        # If either zeta or zeta_order is given, compute the other.
        if zeta is not None:
            zeta = base_ring(zeta)
            if zeta_order is None:
                zeta_order = zeta.multiplicative_order()
        elif zeta_order is not None:
            if not base_ring.is_integral_domain():
                raise ValueError("base ring (= %s) must be an integral domain if only zeta_order is specified"
                                 % base_ring)
            zeta_order = rings.Integer(zeta_order)
            zeta = base_ring.zeta(zeta_order)

        return (base_ring, modulus, zeta, zeta_order)

    def create_object(self, version, key, **extra_args):
        """
        Create the object from the key (extra arguments are ignored). This is
        only called if the object was not found in the cache.

        TESTS::

            sage: K = CyclotomicField(4)
            sage: DirichletGroup.create_object(None, (K, 60, K.gen(), 4))
            Group of Dirichlet characters modulo 60 with values in the group of order 4 generated by zeta4 in Cyclotomic Field of order 4 and degree 2

        """
        base_ring, modulus, zeta, zeta_order = key
        return DirichletGroup_class(base_ring, modulus, zeta, zeta_order)

DirichletGroup = DirichletGroupFactory("DirichletGroup")

def is_DirichletGroup(x):
    """
    Returns True if x is a Dirichlet group.

    EXAMPLES::

        sage: from sage.modular.dirichlet import is_DirichletGroup
        sage: is_DirichletGroup(DirichletGroup(11))
        True
        sage: is_DirichletGroup(11)
        False
        sage: is_DirichletGroup(DirichletGroup(11).0)
        False
    """
    return isinstance(x, DirichletGroup_class)


class DirichletGroup_class(WithEqualityById, Parent):
    """
    Group of Dirichlet characters modulo `N` with values in a ring `R`.
    """

    Element = DirichletCharacter

    def __init__(self, base_ring, modulus, zeta, zeta_order):
        """
        Create a Dirichlet group.

        Not to be called directly (use the factory function ``DirichletGroup``).

        The ``DirichletGroup`` factory ensures that either both
        ``zeta`` and ``zeta_order`` are specified, or that both are
        ``None``.  In the former case, it also ensures that ``zeta``
        is an element of ``base_ring`` and that ``zeta_order`` is an
        element of ``ZZ``.

        TESTS::

            sage: G = DirichletGroup(7, base_ring=Integers(9), zeta=2)  # indirect doctest
            sage: TestSuite(G).run()
            sage: G.base()  # check that Parent.__init__ has been called
            Ring of integers modulo 9

            sage: DirichletGroup(13) == DirichletGroup(13)
            True
            sage: DirichletGroup(13) == DirichletGroup(13, QQ)
            False

        """
        from sage.categories.groups import Groups
        category = Groups().Commutative()
        if base_ring.is_integral_domain() or base_ring.is_finite():
            # The group of n-th roots of unity in the base ring is
            # finite, and hence this Dirichlet group is finite too.
            # In particular, it is finitely generated; the added
            # FinitelyGenerated() here means that the group has a
            # distinguished set of generators.
            category = category.Finite().FinitelyGenerated()
        Parent.__init__(self, base_ring, category=category)
        self._zeta = zeta
        self._zeta_order = zeta_order
        self._modulus = modulus
        self._integers = rings.IntegerModRing(modulus)

    def __setstate__(self, state):
        """
        Used for unpickling old instances.

        TESTS::

            sage: G = DirichletGroup(9)
            sage: loads(dumps(G)) is G
            True
        """
        self._set_element_constructor()
        if '_zeta_order' in state:
            state['_zeta_order'] = rings.Integer(state['_zeta_order'])
        super(DirichletGroup_class, self).__setstate__(state)

    @property
    def _module(self):
        """
        Return the free module used to represent Dirichlet characters.

        TESTS::

            sage: DirichletGroup(12)._module
            Vector space of dimension 2 over Ring of integers modulo 2
        """
        return free_module.FreeModule(rings.IntegerModRing(self.zeta_order()),
                                      len(self.unit_gens()))

    @property
    def _zeta_powers(self):
        """
        Return a list of powers of the distinguished root of unity.

        TESTS::

            sage: DirichletGroup(5)._zeta_powers
            [1, zeta4, -1, -zeta4]
        """
        R = self.base_ring()
        a = R.one()
        w = [a]
        zeta = self.zeta()
        zeta_order = self.zeta_order()
        if is_ComplexField(R):
            for i in range(1, zeta_order):
                a = a * zeta
                a._set_multiplicative_order(zeta_order/gcd(zeta_order, i))
                w.append(a)
        else:
            for i in range(1, zeta_order):
                a = a * zeta
                w.append(a)
        return w

    @property
    def _zeta_dlog(self):
        """
        Return a dictionary that can be used to compute discrete
        logarithms in the value group of this Dirichlet group.

        TESTS::

            sage: DirichletGroup(5)._zeta_dlog
            {-1: 2, -zeta4: 3, zeta4: 1, 1: 0}
        """
        return {z: i for i, z in enumerate(self._zeta_powers)}

    def change_ring(self, R, zeta=None, zeta_order=None):
        """
        Return the base extension of ``self`` to ``R``.

        INPUT:

        - ``R`` -- either a ring admitting a conversion map from the
          base ring of ``self``, or a ring homomorphism with the base
          ring of ``self`` as its domain

        - ``zeta`` -- (optional) root of unity in ``R``

        - ``zeta_order`` -- (optional) order of ``zeta``

        EXAMPLES::

            sage: G = DirichletGroup(7,QQ); G
            Group of Dirichlet characters modulo 7 with values in Rational Field
            sage: G.change_ring(CyclotomicField(6))
            Group of Dirichlet characters modulo 7 with values in Cyclotomic Field of order 6 and degree 2

        TESTS:

        We test the case where `R` is a map (:trac:`18072`)::

            sage: K.<i> = QuadraticField(-1)
            sage: f = K.complex_embeddings()[0]
            sage: D = DirichletGroup(5, K)
            sage: D.change_ring(f)
            Group of Dirichlet characters modulo 5 with values in Complex Field with 53 bits of precision

        """
        if zeta is None and self._zeta is not None:
            # A root of unity was explicitly given; we use it over the
            # new base ring as well.
            zeta = self._zeta
            if zeta_order is None:
                # We reuse _zeta_order if we know that it stays the
                # same; otherwise it will be recomputed as the order
                # of R(zeta) by the DirichletGroup factory.
                p = R.characteristic()
                if p == 0 or p.gcd(self._zeta_order) == 1:
                    zeta_order = self._zeta_order
        else:
            # No root of unity specified; use the same zeta_order
            # (which may still be None).
            zeta_order = self._zeta_order
        # Map zeta to the new parent
        if zeta is not None:
            zeta = R(zeta)
        if isinstance(R, Map):
            R = R.codomain()
        return DirichletGroup(self.modulus(), R,
                              zeta=zeta,
                              zeta_order=zeta_order)

    def base_extend(self, R):
        """
        Return the base extension of ``self`` to ``R``.

        INPUT:

        - ``R`` -- either a ring admitting a *coercion* map from the
          base ring of ``self``, or a ring homomorphism with the base
          ring of ``self`` as its domain

        EXAMPLES::

            sage: G = DirichletGroup(7,QQ); G
            Group of Dirichlet characters modulo 7 with values in Rational Field
            sage: H = G.base_extend(CyclotomicField(6)); H
            Group of Dirichlet characters modulo 7 with values in Cyclotomic Field of order 6 and degree 2

        Note that the root of unity can change::

            sage: H.zeta()
            zeta6

        This method (in contrast to :meth:`change_ring`) requires a
        coercion map to exist::

            sage: G.base_extend(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no coercion map from Rational Field to Integer Ring is defined

        Base-extended Dirichlet groups do not silently get roots of
        unity with smaller order than expected (:trac:`6018`)::

            sage: G = DirichletGroup(10, QQ).base_extend(CyclotomicField(4))
            sage: H = DirichletGroup(10, CyclotomicField(4))
            sage: G is H
            True

            sage: G3 = DirichletGroup(31, CyclotomicField(3))
            sage: G5 = DirichletGroup(31, CyclotomicField(5))
            sage: K30 = CyclotomicField(30)
            sage: G3.gen(0).base_extend(K30) * G5.gen(0).base_extend(K30)
            Dirichlet character modulo 31 of conductor 31 mapping 3 |--> -zeta30^7 + zeta30^5 + zeta30^4 + zeta30^3 - zeta30 - 1

        When a root of unity is specified, base extension still works
        if the new base ring is not an integral domain::

            sage: f = DirichletGroup(17, ZZ, zeta=-1).0
            sage: g = f.base_extend(Integers(15))
            sage: g(3)
            14
            sage: g.parent().zeta()
            14

        """
        if not (isinstance(R, Map)
                or R.has_coerce_map_from(self.base_ring())):
            raise TypeError("no coercion map from %s to %s is defined"
                            % (self.base_ring(), R))
        return self.change_ring(R)

    def _element_constructor_(self, x):
        """
        Construct a Dirichlet character from `x`.

        EXAMPLES::

            sage: G = DirichletGroup(13)
            sage: K = G.base_ring()
            sage: G(1)
            Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1
            sage: G([-1])
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -1
            sage: G([K.0])
            Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta12
            sage: G(0)
            Traceback (most recent call last):
            ...
            TypeError: cannot convert 0 to an element of Group of Dirichlet characters modulo 13 with values in Cyclotomic Field of order 12 and degree 4

            sage: G = DirichletGroup(6)
            sage: G(DirichletGroup(3).0)
            Dirichlet character modulo 6 of conductor 3 mapping 5 |--> -1
            sage: G(DirichletGroup(15).0)
            Dirichlet character modulo 6 of conductor 3 mapping 5 |--> -1
            sage: G(DirichletGroup(15).1)
            Traceback (most recent call last):
            ...
            TypeError: conductor must divide modulus
            sage: H = DirichletGroup(16, QQ); H(DirichletGroup(16).1)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce zeta4 to a rational
        """
        R = self.base_ring()
        try:
            if x == R.one():
                x = [R.one()] * len(self.unit_gens())
        except (TypeError, ValueError, ArithmeticError):
            pass
        if isinstance(x, list):  # list of values on each unit generator
            return self.element_class(self, x)
        elif not isinstance(x, DirichletCharacter):
            raise TypeError("cannot convert %s to an element of %s" % (x, self))
        elif not x.conductor().divides(self.modulus()):
            raise TypeError("conductor must divide modulus")
        a = []
        for u in self.unit_gens():
            v = u.lift()
            # have to do this, since e.g., unit gens mod 11 are not units mod 22.
            while x.modulus().gcd(v) != 1:
                v += self.modulus()
            a.append(R(x(v)))
        return self.element_class(self, a)

    def _coerce_map_from_(self, X):
        """
        Decide whether there is a coercion map from `X`.

        There is conversion between Dirichlet groups of different
        moduli, but no coercion.  This implies that Dirichlet
        characters of different moduli do not compare as equal.

        TESTS::

            sage: trivial_character(6) == trivial_character(3) # indirect doctest
            False
            sage: trivial_character(3) == trivial_character(9)
            False
            sage: trivial_character(3) == DirichletGroup(3, QQ).0^2
            True
        """
        return (isinstance(X, DirichletGroup_class) and
                self.modulus() == X.modulus() and
                self.base_ring().has_coerce_map_from(X.base_ring()) and
                (self._zeta is None or
                 (X._zeta is not None and
                  self.base_ring()(X._zeta) in self._zeta_powers)))

    def __len__(self):
        """
        Return the number of elements of this Dirichlet group. This is the
        same as self.order().

        EXAMPLES::

            sage: len(DirichletGroup(20))
            8
            sage: len(DirichletGroup(20, QQ))
            4
            sage: len(DirichletGroup(20, GF(5)))
            8
            sage: len(DirichletGroup(20, GF(2)))
            1
            sage: len(DirichletGroup(20, GF(3)))
            4
        """
        return self.order()

    def _repr_(self):
        """
        Return a print representation of this group, which can be renamed.

        EXAMPLES::

            sage: G = DirichletGroup(11)
            sage: repr(G) # indirect doctest
            'Group of Dirichlet characters modulo 11 with values in Cyclotomic Field of order 10 and degree 4'
            sage: G.rename('Dir(11)')
            sage: G
            Dir(11)
        """
        s = "Group of Dirichlet characters modulo %s with values in " % self.modulus()
        if self._zeta is not None:
            s += "the group of order %s generated by %s in " % (self._zeta_order, self._zeta)
        s += str(self.base_ring())
        return s

    @cached_method
    def decomposition(self):
        r"""
        Returns the Dirichlet groups of prime power modulus corresponding
        to primes dividing modulus.

        (Note that if the modulus is 2 mod 4, there will be a "factor" of
        `(\ZZ/2\ZZ)^*`, which is the trivial group.)

        EXAMPLES::

            sage: DirichletGroup(20).decomposition()
            [
            Group of Dirichlet characters modulo 4 with values in Cyclotomic Field of order 4 and degree 2,
            Group of Dirichlet characters modulo 5 with values in Cyclotomic Field of order 4 and degree 2
            ]
            sage: DirichletGroup(20,GF(5)).decomposition()
            [
            Group of Dirichlet characters modulo 4 with values in Finite Field of size 5,
            Group of Dirichlet characters modulo 5 with values in Finite Field of size 5
            ]
        """
        R = self.base_ring()
        return Sequence([DirichletGroup(p**r,R) for p, r \
                           in factor(self.modulus())],
                                cr=True,
                                universe = cat.Objects())

    def exponent(self):
        """
        Return the exponent of this group.

        EXAMPLES::

            sage: DirichletGroup(20).exponent()
            4
            sage: DirichletGroup(20,GF(3)).exponent()
            2
            sage: DirichletGroup(20,GF(2)).exponent()
            1
            sage: DirichletGroup(37).exponent()
            36
        """
        return self.zeta_order()

    @cached_method
    def _automorphisms(self):
        """
        Compute the automorphisms of self. These are always given by raising to
        a power, so the return value is a list of integers.

        At present this is only implemented if the base ring has characteristic 0 or a prime.

        EXAMPLES::

            sage: DirichletGroup(17)._automorphisms()
            [1, 3, 5, 7, 9, 11, 13, 15]
            sage: DirichletGroup(17, GF(11^4, 'a'))._automorphisms()
            [1, 11, 121, 1331]
            sage: DirichletGroup(17, Integers(6), zeta=Integers(6)(5))._automorphisms()
            Traceback (most recent call last):
            ...
            NotImplementedError: Automorphisms for finite non-field base rings not implemented
            sage: DirichletGroup(17, Integers(9), zeta=Integers(9)(2))._automorphisms()
            Traceback (most recent call last):
            ...
            NotImplementedError: Automorphisms for finite non-field base rings not implemented
        """
        n = self.zeta_order()
        R = self.base_ring()
        p = R.characteristic()
        if p == 0:
            Auts = [e for e in xrange(1,n) if gcd(e,n) == 1]
        else:
            if not rings.ZZ(p).is_prime():
                raise NotImplementedError("Automorphisms for finite non-field base rings not implemented")
            # The automorphisms in characteristic p are
            # k-th powering for
            #         k = 1, p, p^2, ..., p^(r-1),
            # where p^r = 1 (mod n), so r is the mult order of p modulo n.
            r = rings.IntegerModRing(n)(p).multiplicative_order()
            Auts = [p**m for m in xrange(0,r)]
        return Auts

    def galois_orbits(self, v=None, reps_only=False, sort=True, check=True):
        """
        Return a list of the Galois orbits of Dirichlet characters in self,
        or in v if v is not None.

        INPUT:

        -  ``v`` - (optional) list of elements of self

        -  ``reps_only`` - (optional: default False) if True
           only returns representatives for the orbits.

        -  ``sort`` - (optional: default True) whether to sort
           the list of orbits and the orbits themselves (slightly faster if
           False).

        -  ``check`` - (optional, default: True) whether or not
           to explicitly coerce each element of v into self.

        The Galois group is the absolute Galois group of the prime subfield
        of Frac(R). If R is not a domain, an error will be raised.

        EXAMPLES::

            sage: DirichletGroup(20).galois_orbits()
            [
            [Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -1],
            ...,
            [Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1]
            ]
            sage: DirichletGroup(17, Integers(6), zeta=Integers(6)(5)).galois_orbits()
            Traceback (most recent call last):
            ...
            TypeError: Galois orbits only defined if base ring is an integral domain
            sage: DirichletGroup(17, Integers(9), zeta=Integers(9)(2)).galois_orbits()
            Traceback (most recent call last):
            ...
            TypeError: Galois orbits only defined if base ring is an integral domain
        """
        if v is None:
            v = self.list()
        else:
            if check:
                v = [self(x) for x in v]

        G = []
        seen_so_far = set([])
        for x in v:
            z = x.element()
            e = tuple(z)   # change when there are immutable vectors (and below)
            if e in seen_so_far:
                continue
            orbit = x.galois_orbit(sort=sort)
            if reps_only:
                G.append(x)
            else:
                G.append(orbit)
            for z in orbit:
                seen_so_far.add(tuple(z.element()))
        G = Sequence(G, cr=True)
        if sort:
            G.sort()
        return G

    def gen(self, n=0):
        """
        Return the n-th generator of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.gen(0)
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: G.gen(1)
            Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4
            sage: G.gen(2)
            Traceback (most recent call last):
            ...
            IndexError: n(=2) must be between 0 and 1

        ::

            sage: G.gen(-1)
            Traceback (most recent call last):
            ...
            IndexError: n(=-1) must be between 0 and 1
        """
        n = int(n)
        g = self.gens()
        if n<0 or n>=len(g):
            raise IndexError("n(=%s) must be between 0 and %s"%(n,len(g)-1))
        return g[n]

    @cached_method
    def gens(self):
        """
        Returns generators of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.gens()
            (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1, Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4)
        """
        g = []
        ord = self.zeta_order()
        M = self._module
        zero = M(0)
        orders = self.integers_mod().unit_group().gens_orders()
        for i in range(len(self.unit_gens())):
            z = zero.__copy__()
            z[i] = ord//gcd(ord, orders[i])
            g.append(self.element_class(self, z, check=False))
        return tuple(g)

    def integers_mod(self):
        r"""
        Returns the group of integers `\ZZ/N\ZZ`
        where `N` is the modulus of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.integers_mod()
            Ring of integers modulo 20
        """
        return self._integers

    __iter__ = multiplicative_iterator

    def list(self):
        """
        Return a list of the Dirichlet characters in this group.

        EXAMPLES::

            sage: DirichletGroup(5).list()
            [Dirichlet character modulo 5 of conductor 1 mapping 2 |--> 1,
             Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4,
             Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1,
             Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -zeta4]
        """
        return self._list_from_iterator_cached()

    def modulus(self):
        """
        Returns the modulus of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.modulus()
            20
        """
        return self._modulus

    def ngens(self):
        """
        Returns the number of generators of self.

        EXAMPLES::

            sage: G = DirichletGroup(20)
            sage: G.ngens()
            2
        """
        return len(self.gens())

    @cached_method
    def order(self):
        """
        Return the number of elements of self. This is the same as
        len(self).

        EXAMPLES::

            sage: DirichletGroup(20).order()
            8
            sage: DirichletGroup(37).order()
            36
        """
        ord = rings.Integer(1)
        for g in self.gens():
            ord *= int(g.order())
        return ord

    def random_element(self):
        """
        Return a random element of self.

        The element is computed by multiplying a random power of each
        generator together, where the power is between 0 and the order of
        the generator minus 1, inclusive.

        EXAMPLES::

            sage: DirichletGroup(37).random_element()
            Dirichlet character modulo 37 of conductor 37 mapping 2 |--> zeta36^4
            sage: DirichletGroup(20).random_element()
            Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1
            sage: DirichletGroup(60).random_element()
            Dirichlet character modulo 60 of conductor 3 mapping 31 |--> 1, 41 |--> -1, 37 |--> 1
        """
        e = self(1)
        for i in range(self.ngens()):
            g = self.gen(i)
            n = random.randrange(g.order())
            e *= g**n
        return e

    def unit_gens(self):
        r"""
        Returns the minimal generators for the units of
        `(\ZZ/N\ZZ)^*`, where `N` is the
        modulus of self.

        EXAMPLES::

            sage: DirichletGroup(37).unit_gens()
            (2,)
            sage: DirichletGroup(20).unit_gens()
            (11, 17)
            sage: DirichletGroup(60).unit_gens()
            (31, 41, 37)
            sage: DirichletGroup(20,QQ).unit_gens()
            (11, 17)
        """
        return self._integers.unit_gens()

    @cached_method
    def zeta(self):
        """
        Return the chosen root of unity in the base ring.

        EXAMPLES::

            sage: DirichletGroup(37).zeta()
            zeta36
            sage: DirichletGroup(20).zeta()
            zeta4
            sage: DirichletGroup(60).zeta()
            zeta4
            sage: DirichletGroup(60,QQ).zeta()
            -1
            sage: DirichletGroup(60, GF(25,'a')).zeta()
            2
        """
        zeta = self._zeta
        if zeta is None:
            R = self.base_ring()
            e = self._integers.unit_group_exponent()
            for d in reversed(e.divisors()):
                try:
                    zeta = R.zeta(d)
                    break
                except ValueError:
                    pass
            self.zeta_order.set_cache(d)
        return zeta

    @cached_method
    def zeta_order(self):
        """
        Return the order of the chosen root of unity in the base ring.

        EXAMPLES::

            sage: DirichletGroup(20).zeta_order()
            4
            sage: DirichletGroup(60).zeta_order()
            4
            sage: DirichletGroup(60, GF(25,'a')).zeta_order()
            4
            sage: DirichletGroup(19).zeta_order()
            18
        """
        order = self._zeta_order
        if order is None:
            order = self.zeta().multiplicative_order()
        return order
