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

    def is_floating_point(self):
        """
        Returns whether this `p`-adic ring bounds precision in a floating
        point fashion.

        The relative precision of an element is the power of `p`
        modulo which the unit part of that element is defined.  In a
        floating point ring, elements do not store precision, but arithmetic
        operations truncate to a relative precision depending on the ring.

        EXAMPLES::

            sage: R = ZpCR(5, 15)
            sage: R.is_floating_point()
            False
            sage: R(5^7)
            5^7 + O(5^22)
            sage: S = ZpFP(5, 15)
            sage: S.is_floating_point()
            True
            sage: S(5^7)
            5^7
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

        EXAMPLES::

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

    def _matrix_smith_form(self, M, transformation, integral):
        r"""
        Return the Smith normal form of the matrix ``M``.

        This method gets called by
        :meth:`sage.matrix.matrix2.Matrix.smith_form` to compute the Smith
        normal form over local rings and fields.

        The entries of the Smith normal form are normalized such that non-zero
        entries of the diagonal are powers of the distinguished uniformizer.

        INPUT:

        - ``M`` -- a matrix over this ring

        - ``transformation`` -- a boolean; whether the transformation matrices
        are returned

        - ``integral`` -- a subring of the base ring or ``True``; the entries
        of the transformation matrices are in this ring.  If ``True``, the
        entries are in the ring of integers of the base ring.

        EXAMPLES::

            sage: A = Zp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

            sage: S, L, R = M.smith_form()  # indirect doctest
            sage: S
            [ ...1     0]
            [    0 ...10]
            sage: L
            [...222222223          ...]
            [...444444444         ...2]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        If not needed, it is possible to avoid the computation of
        the transformations matrices `L` and `R`::

            sage: M.smith_form(transformation=False)  # indirect doctest
            [ ...1     0]
            [    0 ...10]

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: S, L, R = M.smith_form()  # indirect doctest
            sage: S
            [ ...1     0]
            [    0 ...10]
            [    0     0]
            sage: L
            [...222222223          ...          ...]
            [...444444444         ...2          ...]
            [...444444443         ...1         ...1]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        An error is raised if the precision on the entries is
        not enough to determine the Smith normal form::

            sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
            sage: M.smith_form()  # indirect doctest
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to compute Smith normal form

        """
        if integral is None:
            integral = self
        if integral is True:
            integral = self.integer_ring()
        if integral is not self and integral is not self.integer_ring():
            raise NotImplementedError("Smith normal form over this subring")

        n = M.nrows()
        m = M.ncols()
        S = M.parent()(M.list())
        smith = M.parent()(0)
        R = M.base_ring()
        if R.tracks_precision():
            precM = min([ x.precision_absolute() for x in M.list() ])

        if transformation:
            from sage.matrix.special import identity_matrix
            left = identity_matrix(R,n)
            right = identity_matrix(R,m)

        from sage.rings.all import Infinity
        val = -Infinity
        for piv in range(min(n,m)):
            curval = Infinity
            pivi = pivj = piv
            for i in range(piv,n):
                for j in range(piv,m):
                    v = S[i,j].valuation()
                    if v < curval:
                        pivi = i; pivj = j
                        curval = v
                        if v == val: break
                else:
                    continue
                break
            val = curval

            if R.tracks_precision() and precM is not Infinity and val >= precM:
                from .precision_error import PrecisionError
                raise PrecisionError("not enough precision to compute Smith normal form")

            if val is Infinity:
                break

            S.swap_rows(pivi,piv)
            S.swap_columns(pivj,piv)
            if transformation:
                left.swap_rows(pivi,piv)
                right.swap_columns(pivj,piv)

            smith[piv,piv] = R(1) << val
            inv = ~(S[piv,piv] >> val)
            for i in range(piv+1,n):
                scalar = -inv * (S[i,piv] >> val)
                if R.tracks_precision():
                    scalar = scalar.lift_to_precision()
                S.add_multiple_of_row(i,piv,scalar,piv+1)
                if transformation:
                    left.add_multiple_of_row(i,piv,scalar)
            if transformation:
                left.rescale_row(piv,inv)
                for j in range(piv+1,m):
                    scalar = -inv * (S[piv,j] >> val)
                    if R.tracks_precision():
                        scalar = scalar.lift_to_precision()
                    right.add_multiple_of_column(j,piv,scalar)

        if transformation:
            if R.tracks_precision() and precM is not Infinity:
                left = left.apply_map(lambda x: x.add_bigoh(precM-val))
            left = left.change_ring(integral)
            right = right.change_ring(integral)
            return smith, left, right
        else:
            return smith

    def _test_matrix_smith(self, **options):
        r"""
        Test that :meth:`_matrix_smith_form` works correctly.

        EXAMPLES::

            sage: ZpCA(5, 15)._test_matrix_smith()

        """
        tester = self._tester(**options)
        tester.assertEqual(self.residue_field().characteristic(), self.residue_characteristic())

        from itertools import chain
        from sage.all import MatrixSpace
        from .precision_error import PrecisionError
        matrices = chain(*[MatrixSpace(self, n, m).some_elements() for n in (1,3,7) for m in (1,4,7)])
        for M in tester.some_elements(matrices):
            bases = [self]
            if self is not self.integer_ring():
                bases.append(self.integer_ring())
            for base in bases:
                try:
                    S,U,V = M.smith_form(integral=base)
                except PrecisionError:
                    continue

                if self.is_exact() or self.tracks_precision():
                    tester.assertEqual(U*M*V, S)

                tester.assertEqual(U.nrows(), U.ncols())
                tester.assertEqual(U.base_ring(), base)

                tester.assertEqual(V.nrows(), V.ncols())
                tester.assertEqual(V.base_ring(), base)

                for d in S.diagonal():
                    if not d.is_zero():
                        tester.assertTrue(d.unit_part().is_one())

                for (d,dd) in zip(S.diagonal(), S.diagonal()[1:]):
                    tester.assertTrue(d.divides(dd))
