"""
Local Generic

Superclass for `p`-adic and power series rings.

AUTHORS:

- David Roe
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from sage.rings.ring import CommutativeRing
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings, CompleteDiscreteValuationFields
from sage.structure.category_object import check_default_category
from sage.structure.parent import Parent
from sage.rings.integer import Integer

class LocalGeneric(CommutativeRing):
    def __init__(self, base, prec, names, element_class, category=None):
        """
        Initializes self.

        EXAMPLES::

            sage: R = Zp(5) #indirect doctest
            sage: R.precision_cap()
            20

        In :trac:`14084`, the category framework has been implemented for p-adic rings::

            sage: TestSuite(R).run()
            sage: K = Qp(7)
            sage: TestSuite(K).run()

        TESTS::

            sage: R = Zp(5, 5, 'fixed-mod')
            sage: R._repr_option('element_is_atomic')
            False
        """
        self._prec = prec
        self.Element = element_class
        default_category = getattr(self, '_default_category', None)
        if self.is_field():
            category = CompleteDiscreteValuationFields()
        else:
            category = CompleteDiscreteValuationRings()
        category = category.Metric().Complete()
        if default_category is not None:
            category = check_default_category(default_category, category)
        Parent.__init__(self, base, names=(names,), normalize=False, category=category, element_constructor=element_class)

    def is_capped_relative(self):
        """
        Returns whether this `p`-adic ring bounds precision in a capped
        relative fashion.

        The relative precision of an element is the power of `p`
        modulo which the unit part of that element is defined.  In a
        capped relative ring, the relative precision of elements are
        bounded by a constant depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_relative()
            False
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_relative()
            True
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return False

    def is_capped_absolute(self):
        """
        Returns whether this `p`-adic ring bounds precision in a
        capped absolute fashion.

        The absolute precision of an element is the power of `p`
        modulo which that element is defined.  In a capped absolute
        ring, the absolute precision of elements are bounded by a
        constant depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_absolute()
            True
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_absolute()
            False
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return False

    def is_fixed_mod(self):
        """
        Returns whether this `p`-adic ring bounds precision in a fixed
        modulus fashion.

        The absolute precision of an element is the power of `p`
        modulo which that element is defined.  In a fixed modulus
        ring, the absolute precision of every element is defined to be
        the precision cap of the parent.  This means that some
        operations, such as division by `p`, don't return a well defined
        answer.

        EXAMPLES::

            sage: R = ZpFM(5,15)
            sage: R.is_fixed_mod()
            True
            sage: R(5^7,absprec=9)
            5^7 + O(5^15)
            sage: S = ZpCA(5, 15)
            sage: S.is_fixed_mod()
            False
            sage: S(5^7,absprec=9)
            5^7 + O(5^9)
        """
        return False

    def is_lazy(self):
        """
        Returns whether this `p`-adic ring bounds precision in a lazy
        fashion.

        In a lazy ring, elements have mechanisms for computing
        themselves to greater precision.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.is_lazy()
            False
        """
        return False

    def _latex_(self):
        r"""
        Latex.

        EXAMPLES::

            sage: latex(Zq(27,names='a')) #indirect doctest
            \mathbf{Z}_{3^{3}}
        """
        return self._repr_(do_latex = True)

    def change(self, **kwds):
        """
        Return a new ring with changed attributes.

        INPUT:

        Keyword arguments are passed on to the :mod:`constructors <sage.rings.padics.factory>`
        via the :meth:`construction` functor.

        NOTES:

        For extension rings, some keywords take effect only on the extension
        (``names``, ``var_name``, ``res_name``, ``unram_name``, ``ram_name``)
        and others on both the extension and, recursively, the base ring
        (``print_mode``, ``halt``, ``print_pos``, ``print_sep``,
         ``print_alphabet``, ``print_max_ram_terms``, ``check``).

        If the precision is increased on an extension ring,
        the precision on the base is increased as necessary.
        If the precision is decreased, the precision of the base is unchanged.

        EXAMPLES:

        We can use this method to change the precision::

            sage: Zp(5).change(prec=40)
            5-adic Ring with capped relative precision 40

        or the precision type::

            sage: Zp(5).change(type="capped-abs")
            5-adic Ring with capped absolute precision 20

        or even the prime::

            sage: ZpCA(3).change(p=17)
            17-adic Ring with capped absolute precision 20

        You can switch between the ring of integers and its fraction field::

            sage: ZpCA(3).change(field=True)
            3-adic Field with capped relative precision 20

        You can also change print modes::

            sage: R = Zp(5).change(prec=5, print_mode='digits')
            sage: ~R(17)
            ...13403

        You can change extensions::

            sage: K.<a> = QqFP(125, prec=4)
            sage: K.change(q=64)
            Unramified Extension of 2-adic Field with floating precision 4 in a defined by x^6 + x^4 + x^3 + x + 1
            sage: R.<x> = QQ[]
            sage: K.change(modulus = x^2 - x + 2)
            Unramified Extension of 5-adic Field with floating precision 4 in a defined by x^2 - x + 2

        and variable names::

            sage: K.change(names='b')
            Unramified Extension of 5-adic Field with floating precision 4 in b defined by x^3 + 3*x + 3

        and precision::

            sage: Kup = K.change(prec=8); Kup
            Unramified Extension of 5-adic Field with floating precision 8 in a defined by x^3 + 3*x + 3
            sage: Kup.base_ring()
            5-adic Field with floating precision 8

        If you decrease the precision, the precision of the base stays the same::

            sage: Kdown = K.change(prec=2); Kdown
            Unramified Extension of 5-adic Field with floating precision 2 in a defined by x^3 + 3*x + 3
            sage: Kdown.base_ring()
            5-adic Field with floating precision 4
        """
        functor, ring = self.construction()
        functor = copy(functor)
        # There are two kinds of functors possible:
        # CompletionFunctor and AlgebraicExtensionFunctor
        # We distinguish them by the presence of ``prec``,
        if hasattr(functor, "prec"):
            functor.extras = copy(functor.extras)
            if 'type' in kwds and kwds['type'] not in functor._dvr_types:
                raise ValueError("completion type must be one of %s"%(", ".join(functor._dvr_types[1:])))
            if 'field' in kwds:
                field = kwds.pop('field')
                if field:
                    ring = ring.fraction_field()
                else:
                    ring = ring.ring_of_integers()
            for atr in ('p', 'prec', 'type'):
                if atr in kwds:
                    setattr(functor, atr, kwds.pop(atr))
            for atr in ('print_mode', 'halt', 'names', 'ram_name', 'print_pos', 'print_sep', 'print_alphabet', 'print_max_terms', 'check'):
                if atr in kwds:
                    functor.extras[atr] = kwds.pop(atr)
            if kwds:
                raise ValueError("Extra arguments received: %s"%(", ".join(kwds.keys())))
        else:
            functor.kwds = copy(functor.kwds)
            if 'prec' in kwds:
                prec = kwds.pop('prec')
                baseprec = (prec - 1) // self.e() + 1
                if baseprec > self.base_ring().precision_cap():
                    kwds['prec'] = baseprec
                functor.kwds['prec'] = prec
            from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
            n = None
            if 'q' in kwds and isinstance(ring, pAdicBaseGeneric):
                q = kwds.pop('q')
                if not isinstance(q, Integer):
                    raise TypeError("q must be an integer")
                p, n = q.is_prime_power(get_data=True)
                if n == 0:
                    raise ValueError("q must be a prime power")
                if 'p' in kwds and kwds['p'] != p:
                    raise ValueError("q does not match p")
                kwds['p'] = p
            if 'modulus' in kwds:
                modulus = kwds.pop('modulus')
                if n is not None and modulus.degree() != n:
                    raise ValueError("modulus must have degree matching q")
                functor.polys = [modulus]
            elif n is not None:
                functor.polys = [n]
            for atr in ('names', 'var_name', 'res_name', 'unram_name', 'ram_name'):
                functor.kwds[atr] = kwds.pop(atr)
            for atr in ('print_mode', 'halt', 'print_pos', 'print_sep', 'print_alphabet', 'print_max_terms', 'check'):
                functor.kwds[atr] = kwds[atr]
            try:
                ring = ring.change(**kwds)
            except AttributeError:
                raise NotImplementedError
        return functor(ring)

    def precision_cap(self):
        r"""
        Returns the precision cap for ``self``.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - integer -- ``self``'s precision cap

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.precision_cap()
            10
            sage: R = Zp(3, 10,'capped-rel'); R.precision_cap()
            10
            sage: R = Zp(3, 10,'capped-abs'); R.precision_cap()
            10

        NOTES::

            This will have different meanings depending on the type of
            local ring.  For fixed modulus rings, all elements are
            considered modulo ``self.prime()^self.precision_cap()``.
            For rings with an absolute cap (i.e. the class
            ``pAdicRingCappedAbsolute``), each element has a precision
            that is tracked and is bounded above by
            ``self.precision_cap()``.  Rings with relative caps
            (e.g. the class ``pAdicRingCappedRelative``) are the same
            except that the precision is the precision of the unit
            part of each element.  For lazy rings, this gives the
            initial precision to which elements are computed.
        """
        return self._prec

    def is_exact(self):
        r"""
        Returns whether this p-adic ring is exact, i.e. False.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            boolean -- whether self is exact, i.e. False.

        EXAMPLES:
            #sage: R = Zp(5, 3, 'lazy'); R.is_exact()
            #False
            sage: R = Zp(5, 3, 'fixed-mod'); R.is_exact()
            False
        """
        return False

    def residue_characteristic(self):
        r"""
        Returns the characteristic of ``self``'s residue field.

        INPUT:

        - ``self`` -- a p-adic ring.

        OUTPUT:

        - integer -- the characteristic of the residue field.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.residue_characteristic()
            3
        """
        return self.residue_class_field().characteristic()

    def defining_polynomial(self, var = 'x'):
        r"""
        Returns the defining polynomial of this local ring, i.e. just ``x``.

        INPUT:

        - ``self`` -- a local ring
        - ``var`` -- string (default: ``'x'``) the name of the variable

        OUTPUT:

        - polynomial -- the defining polynomial of this ring as an extension over its ground ring

        EXAMPLES::

            sage: R = Zp(3, 3, 'fixed-mod'); R.defining_polynomial('foo')
            (1 + O(3^3))*foo + (O(3^3))
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(self, var).gen()

    def ground_ring(self):
        r"""
        Returns ``self``.

        Will be overridden by extensions.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - the ground ring of ``self``, i.e., itself

        EXAMPLES::

            sage: R = Zp(3, 5, 'fixed-mod')
            sage: S = Zp(3, 4, 'fixed-mod')
            sage: R.ground_ring() is R
            True
            sage: S.ground_ring() is R
            False
        """
        return self

    def ground_ring_of_tower(self):
        r"""
        Returns ``self``.

        Will be overridden by extensions.

        INPUT:

        - ``self`` -- a `p`-adic ring

        OUTPUT:

        - the ground ring of the tower for ``self``, i.e., itself

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.ground_ring_of_tower()
            5-adic Ring with capped relative precision 20
        """
        return self

    def degree(self):
        r"""
        Returns the degree of ``self`` over the ground ring, i.e. 1.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - integer -- the degree of this ring, i.e., 1

        EXAMPLES::

            sage: R = Zp(3, 10, 'capped-rel'); R.degree()
            1
        """
        return Integer(1)

    def ramification_index(self, K = None):
        r"""
        Returns the ramification index over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - integer -- the ramification index of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.ramification_index()
            1
        """
        if K is None or K is self:
            return Integer(1)
        else:
            raise ValueError("K should be a subring of self")

    def e(self, K = None):
        r"""
        Returns the ramification index over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring of ``self`` (default ``None``)

        OUTPUT:

        - integer -- the ramification index of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.e()
            1
        """
        return self.ramification_index(K)

    def inertia_degree(self, K=None):
        r"""
        Returns the inertia degree over ``K`` (defaults to the ground ring): 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring of ``self`` (default None)

        OUTPUT:

        - integer -- the inertia degree of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.inertia_degree()
            1
        """
        return Integer(1)

    def residue_class_degree(self, K=None):
        r"""
        Returns the inertia degree over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring (default ``None``)

        OUTPUT:

        - integer -- the inertia degree of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.residue_class_degree()
            1
        """
        return self.inertia_degree(K)

    def f(self, K=None):
        r"""
        Returns the inertia degree over the ground ring: 1 unless overridden.

        INPUT:

        - ``self`` -- a local ring
        - ``K`` -- a subring (default ``None``)

        OUTPUT:

        - integer -- the inertia degree of this ring: 1 unless overridden.

        EXAMPLES::

            sage: R = Zp(3, 5, 'capped-rel'); R.f()
            1
        """
        return self.inertia_degree(K)

    def inertia_subring(self):
        r"""
        Returns the inertia subring, i.e. ``self``.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - the inertia subring of self, i.e., itself

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.inertia_subring()
            5-adic Ring with capped relative precision 20
        """
        return self

    def maximal_unramified_subextension(self):
        r"""
        Returns the maximal unramified subextension.

        INPUT:

        - ``self`` -- a local ring

        OUTPUT:

        - the maximal unramified subextension of ``self``

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.maximal_unramified_subextension()
            5-adic Ring with capped relative precision 20
        """
        return self.inertia_subring()

#    def get_extension(self):
#        r"""
#        Returns the trivial extension of self.
#        """
#        raise NotImplementedError

    def uniformiser(self):
        """
        Returns a uniformiser for ``self``, ie a generator for the unique maximal ideal.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.uniformiser()
            5 + O(5^21)
            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.uniformiser()
            t + O(t^21)
        """
        return self.uniformizer()

    def uniformiser_pow(self, n):
        """
        Returns the `n`th power of the uniformiser of ``self`` (as an element of ``self``).

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.uniformiser_pow(5)
            5^5 + O(5^25)
        """
        return self.uniformizer_pow(n)

    def is_finite(self):
        r"""
        Returns whether this ring is finite, i.e. ``False``.

        INPUT:

        - ``self`` -- a `p`-adic ring

        OUTPUT:

        - boolean -- whether self is finite, i.e., ``False``

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.is_finite()
            False
        """
        return False

    def ext(self, *args, **kwds):
        """
        Constructs an extension of self.  See ``extension`` for more details.

        EXAMPLES::

            sage: A = Zp(7,10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2+7)
            sage: B.uniformiser()
            t + O(t^21)
        """
        return self.extension(*args, **kwds)

    def _test_residue(self, **options):
        r"""
        Perform some tests on the residue field of this ring.

        EXAMPLES::

            sage: R = Zp(2)
            sage: R._test_residue()

        """
        tester = self._tester(**options)
        tester.assertEqual(self.residue_field().characteristic(), self.residue_characteristic())

        for x in tester.some_elements():
            errors = []
            if x.precision_absolute() <= 0:
                from .precision_error import PrecisionError
                errors.append(PrecisionError)
            if x.valuation() < 0:
                errors.append(ValueError)
            if errors:
                with tester.assertRaises(tuple(errors)):
                    x.residue()
                continue
            y = x.residue()
            # residue() is in `Z/pZ` which is not identical to the residue field `F_p`
            tester.assertEqual(y.parent().cardinality(), self.residue_field().cardinality())
            z = self(y)
            tester.assertGreater((x-z).valuation(), 0)

        for x in self.residue_field().some_elements():
            y = self(x)
            if x.is_zero():
                tester.assertGreater(y.valuation(), 0)
            else:
                tester.assertEqual(y.valuation(), 0)
            z = y.residue()
            tester.assertEqual(x, z)
