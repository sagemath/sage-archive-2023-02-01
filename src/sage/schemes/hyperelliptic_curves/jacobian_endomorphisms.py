# -*- coding: utf-8 -*-
r"""
Endomorphism rings of Jacobians of hyperelliptic curves.

Given the Jacobian `J` of a hyperelliptic curve over a number field `K`,
the class ``EndomorphismRing`` constructs the ring of endomorphisms of J
which are defined over `K`.

Currently the only functionality is for genus 2 curves over the rational
numbers `\QQ`, for which there are two methods associated to an instance
`E` of ``EndomorphismRing``:

- ``is_absolutely_field(E)``
- ``is_absolutely_trivial(E)``

``is_absolutely_field(E)`` determines whether the endomorphism algebra (which
is the endomorphism ring tensored with `\QQ`) is a field.

``is_absolutely_trivial(E)`` determines whether the endomorphism ring is
equal to the integer ring `\ZZ`. If this is the case, one says that `J` is
``generic``.

Both of these are important attributes of the Jacobian J.

The algorithms of these two methods are implementations of Algorithms 4.10 and
4.15 from Lombardo's paper [Lom2019]_.

The following examples have been verified with the corresponding LMFDB entries.

EXAMPLES::

Here is an example of a generic Jacobian; the LMFDB label of the curve is
249.a.249.1::

    sage: f = x^6 + 2*x^3 + 4*x^2 + 4*x + 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: E = A.endomorphism_ring()
    sage: E.is_absolutely_trivial()
    True
    sage: E.is_absolutely_simple()
    True

Here is an example of a Jacobian whose endomorphism algebra is a field but not
the rational number field; the LMFDB label of the curve is 529.a.529.1::

    sage: f = y^2 = x^6 - 4*x^5 + 2*x^4 + 2*x^3 + x^2 + 2*x + 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: E = A.endomorphism_ring()
    sage: E.is_absolutely_trivial()
    False
    sage: E.is_absolutely_simple()
    True

Here is an example of a Jacobian whose endomorphism algebra is not a field;
the LMFDB label of the curve is 169.a.169.1::

    sage: f = x^6 + 4*x^5 + 6*x^4 + 2*x^3 + x^2 + 2*x + 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: E = A.endomorphism_ring()
    sage: E.is_absolutely_trivial()
    False
    sage: E.is_absolutely_simple()
    False

.. WARNING::

    There is a very small chance that the algorithms return ``False`` for the
    two methods described above when in fact one or both of them are ``True``.
    In this case, as explained in the discussion immediately preceding
    Algorithm 4.15 of [Lom2019]_, this can be established by increasing the
    optional `B` parameter of the two methods. Mathematically, the algorithms
    given the correct answer only in the limit as `B \to \infty`, although in
    practice `B = 200` was sufficient to correctly verify every single entry
    in the LMFDB.

AUTHORS:

- Barinder S. Banwait and Davide Lombardo (2021-06-09): initial version

"""

# ****************************************************************************
#       Copyright (C) 2021 Barinder S. Banwait and Davide Lombardo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.all import QQ, ZZ, PolynomialRing, FiniteField, NumberField
from sage.misc.lazy_import import lazy_import
from sage.rings.fast_arith import prime_range
from sage.arith.all import gcd
lazy_import('sage.interfaces.genus2reduction', ['genus2reduction', 'Genus2reduction'])


def satisfies_coefficient_condition(g,p):
    """This is the coefficient condition in the definition of Omega_K'
    on page 912 of the published version of paper"""

    if g[0] != p^2:
        return False
    if g[3]*p != g[1]:
        return False
    if g[2]%p == 0:
        return False
    return True


def get_is_geom_irred(f,C,bad_primes,B=200):
    """Implements Algorithm 4.10 from the paper.

    An additional optimisation comes from Part (2) of Theorem 4.8, from
    which we can conclude that the endomorphism ring is absolutely simple.

    Args:
        f (polynomial): Polynomial defining the curve
        C ([HyperellipticCrv]): the curve
        bad_primes (list): List of odd primes of bad reduction
        B (int, optional): Bound in the algorithm. Defaults to 200.

    Returns:
        (bool1, bool2): Pair of booleans. `bool1` indicates if the
        endomorphism algebra is a field; `bool2` indicates if the
        endomorphism algebra is the field of rational numbers.
    """

    if C.has_odd_degree_model():
        C_odd = C.odd_degree_model()
        f_odd, h_odd = C_odd.hyperelliptic_polynomials()
        # if f was odd to begin with, then f_odd = f
        assert f_odd.degree() == 5

        if (4*f_odd + h_odd**2).degree() == 5:
            f_new = 4*f_odd + h_odd**2
            if f_new.is_irreducible():
                # i.e. the Jacobian is absolutely simple
                f_disc_odd_prime_exponents = [v for _,v in f_new.discriminant().prime_to_S_part([2]).factor()]
                if 1 in f_disc_odd_prime_exponents:
                    return (True, True)  # Theorem 4.8 (2)
                return (True, False)  # Algorithm 4.10 Step 1

    if f.is_irreducible():
        assert f.degree() == 6  # else we should already have exited by now
        G = f.galois_group()
        if G.order() in [360, 720]:
            return (True, True)  # Algorithm 4.10 Step 2

    R = PolynomialRing(ZZ,2,"xv")
    x,v = R.gens()
    T = PolynomialRing(QQ,'v')
    g = v - x**12

    for p in prime_range(3,B):
        if p not in bad_primes:
            fp = C.change_ring(FiniteField(p)).frobenius_polynomial()

            # This defines the polynomial f_v**[12] from the paper
            fp12 = T(R(fp).resultant(g))

            if fp12.is_irreducible():
                # i.e. the Jacobian is absolutely simpl
                f_disc_odd_prime_exponents = [v for _,v in f.discriminant().prime_to_S_part([ZZ(2)]).factor()]
                if 1 in f_disc_odd_prime_exponents :
                    return (True, True)  # Theorem 4.8 (2)
                return (True, False)  # Algorithm 4.10 Step 3
    return (False, False)


def is_geom_trivial_when_field(C, bad_primes, B=200):
    """Codes up Lombardo's Algorithm 4.15"""

    running_gcd = 0
    R = PolynomialRing(ZZ,2,"xv")
    x,v = R.gens()
    T = PolynomialRing(QQ,'v')
    g = v - x**4

    for p in prime_range(3,B):
        if p not in bad_primes:
            Cp = C.change_ring(FiniteField(p))
            fp = Cp.frobenius_polynomial()
            if satisfies_coefficient_condition(fp, p):
                # This defines the polynomial f_v**[4] from the paper
                fp4 = T(R(fp).resultant(g))
                if fp4.is_irreducible():
                    running_gcd = gcd(running_gcd, NumberField(fp,'a').discriminant())
                    if running_gcd <= 24:
                        return True
    return False


class EndomorphismRing(SageObject):
    r"""
    The ring of endomorphisms of this Jacobian which are defined over the
    base field.

    More text here.

    EXAMPLES::

        sage: R.<x>=QQ[]
        sage: f = x**5 + 17
        sage: C = HyperellipticCurve(f)
        sage: J = C.jacobian()
        sage: E = J.endomorphism_ring()
        sage: E
        Endomorphism ring of Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + 17
    """

    def __init__(self, A):
        r"""
        see ``EndomorphismRing`` for documentation
        """
        self._A = A
        self._have_established_absolutely_field = False
        self._have_established_absolutely_trivial = False

    def __repr__(self):
        r"""
        string representation of the class

        EXAMPLES::

            sage: rho = EllipticCurve([0,1]).galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y**2 = x**3 + 1 over Rational Field

        """
        return "Endomorphism ring of " + repr(self._A)

    def is_absolutely_field(self, B=200):
        if self._have_established_absolutely_field:
            return True
        C = self._A.curve()
        if C.genus() != 2:
            raise NotImplementedError("Current implementation requires the curve to be of genus 2")
        if C.base_ring() != QQ:
            raise NotImplementedError("Current implementation requires the curve to be defined over the rationals")
        f, h = C.hyperelliptic_polynomials()
        if h != 0:
            raise NotImplementedError("Current implementation requires the curve to be in the form y^2 = f(x)")
        red_data = genus2reduction(0,f)
        cond_C = red_data.conductor  # WARNING: this is only the prime_to_2 conductor.
        bad_primes = cond_C.prime_divisors()
        self.bad_primes = bad_primes

        is_abs_simp, is_def_geom_trivial = get_is_geom_irred(f, C, bad_primes, B)

        if is_def_geom_trivial:
            self._have_established_absolutely_trivial = True
        if is_abs_simp:
            self._have_established_absolutely_field = True
            return True
        return False

    def is_absolutely_trivial(self, B=200):

        if self._have_established_absolutely_trivial:
            return True
        is_abs_simple = self.is_absolutely_field()
        if self._have_established_absolutely_trivial:
            return True
        if is_abs_simple and is_geom_trivial_when_field(self._A.curve(), self.bad_primes):
            return True
        return False
