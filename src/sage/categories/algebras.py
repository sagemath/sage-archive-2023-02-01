r"""
Algebras

AUTHORS:

 - David Kohel & William Stein (2005): initial revision
 - Nicolas M. Thiery (2008): rewrote for new category framework
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.homset import Hom
from sage.categories.cartesian_product import CartesianProductsCategory, cartesian_product
from sage.categories.dual import DualObjectsCategory
from sage.categories.tensor import TensorProductsCategory, tensor
from sage.categories.morphism import SetMorphism
from sage.categories.rings import Rings
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import have_same_parent

class Algebras(Category_over_base_ring):
    """
    The category of algebras over a given base ring.

    An algebra over a ring `R` is a module over `R` which is itself a ring.

    TODO: should `R` be a commutative ring?

    EXAMPLES::

        sage: Algebras(ZZ)
        Category of algebras over Integer Ring
        sage: Algebras(ZZ).super_categories()
        [Category of rings, Category of modules over Integer Ring]

    TESTS::

        sage: TestSuite(Algebras(ZZ)).run()
    """

    # For backward compatibility?
    def __contains__(self, x):
        """
        Membership testing

        EXAMPLES::

            sage: QQ[x] in Algebras(QQ)
            True

            sage: QQ^3 in Algebras(QQ)
            False
            sage: QQ[x] in Algebras(CDF)
            False
        """
        if super(Algebras, self).__contains__(x):
            return True
        from sage.rings.ring import Algebra
        return isinstance(x, Algebra) and x.base_ring() == self.base_ring()

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(ZZ).super_categories()
            [Category of rings, Category of modules over Integer Ring]
        """
        R = self.base_ring()
        from sage.categories.rings import Rings
        from sage.categories.modules import Modules
        return [Rings(), Modules(R)]

    class ParentMethods: # (Algebra):  # Eventually, the content of Algebra should be moved here
        def from_base_ring(self, r):
            """
            Canonical embedding from base ring

            INPUT:

             - ``r`` -- an element of ``self.base_ring()``

            Returns the canonical embedding of `r` into self.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: A.from_base_ring(1)
                B[word: ]
            """
            return self.one()._lmul_(r)

        def __init_extra__(self):
            """
            Declares the canonical coercion from ``self.base_ring()`` to ``self``,
            if there has been none before.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: coercion_model = sage.structure.element.get_coercion_model()
                sage: coercion_model.discover_coercion(QQ, A)
                (Generic morphism:
                  From: Rational Field
                  To:   An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field, None)
                sage: A(1)          # indirect doctest
                B[word: ]

            """
            # If self has an attribute _no_generic_basering_coercion
            # set to True, then this declaration is skipped.
            # This trick, introduced in #11900, is used in
            # sage.matrix.matrix_space.py and
            # sage.rings.polynomial.polynomial_ring.
            # It will hopefully be refactored into something more
            # conceptual later on.
            if getattr(self,'_no_generic_basering_coercion',False):
                return

            base_ring = self.base_ring()
            if base_ring is self:
                # There are rings that are their own base rings. No need to register that.
                return
            if self.is_coercion_cached(base_ring):
                # We will not use any generic stuff, since a (presumably) better conversion
                # has already been registered.
                return
            mor = None
            # This could be a morphism of Algebras(self.base_ring()); however, e.g., QQ is not in Algebras(QQ)
            H = Hom(base_ring, self, Rings())

            # Idea: There is a generic method "from_base_ring", that just does multiplication with
            # the multiplicative unit. However, the unit is constructed repeatedly, which is slow.
            # Hence, if the unit is available *now*, then we store it.
            #
            # However, if there is a specialised from_base_ring method, then it should be used!
            try:
                has_custom_conversion = self.category().parent_class.from_base_ring.__func__ is not self.from_base_ring.__func__
            except AttributeError:
                # Sometimes from_base_ring is a lazy attribute
                has_custom_conversion = True
            if has_custom_conversion:
                mor = SetMorphism(function = self.from_base_ring, parent = H)
                try:
                    self.register_coercion(mor)
                except AssertionError:
                    pass
                return

            try:
                one = self.one()
            except (NotImplementedError, AttributeError, TypeError):
                # The unit is not available, yet. But there are cases
                # in which it will be available later. Hence:
                mor = SetMorphism(function = self.from_base_ring, parent = H)
            # try sanity of one._lmul_
            if mor is None:
                try:
                    if one._lmul_(base_ring.an_element()) is None:
                        # There are cases in which lmul returns None, believe it or not.
                        # One example: Hecke algebras.
                        # In that case, the generic implementation of from_base_ring would
                        # fail as well. Hence, unless it is overruled, we will not use it.
                        #mor = SetMorphism(function = self.from_base_ring, parent = H)
                        return
                except (NotImplementedError, AttributeError, TypeError):
                    # it is possible that an_element or lmul are not implemented.
                    return
                    #mor = SetMorphism(function = self.from_base_ring, parent = H)
                mor = SetMorphism(function = one._lmul_, parent = H)
            try:
                self.register_coercion(mor)
            except AssertionError:
                pass

    class ElementMethods:
        # TODO: move the content of AlgebraElement here

        # Workaround: this sets back Semigroups.Element.__mul__, which is currently overriden by Modules.Element.__mul__
        # What does this mean in terms of inheritance order?
        # Could we do a variant like __mul__ = Semigroups.Element.__mul__
        def __mul__(self, right):
            """
            EXAMPLES::

                sage: s = SFASchur(QQ)
                sage: a = s([2])
                sage: a._mul_(a) #indirect doctest
                s[2, 2] + s[3, 1] + s[4]

            Todo: use AlgebrasWithBasis(QQ).example()
            """
            if have_same_parent(self, right) and hasattr(self, "_mul_"):
                return self._mul_(right)
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(self, right, operator.mul)

#        __imul__ = __mul__

        # Parents in this category should implement _lmul_ and _rmul_

        def _div_(self, y):
            """
            Division by invertible elements

            # TODO: move in Monoids

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ).example()
                sage: x = C(2); x
                2*B[word: ]
                sage: y = C.algebra_generators().first(); y
                B[word: a]

                sage: y._div_(x)
                1/2*B[word: a]
                sage: x._div_(y)
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= B[word: a])
            """
            return self.parent().product(self, ~y)

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of algebras constructed as cartesian products of algebras

        This construction gives the direct product of algebras. See
        discussion on:

         - http://groups.google.fr/group/sage-devel/browse_thread/thread/35a72b1d0a2fc77a/348f42ae77a66d16#348f42ae77a66d16
         - http://en.wikipedia.org/wiki/Direct_product
        """
        def extra_super_categories(self):
            """
            A cartesian product of algebras is endowed with a natural
            algebra structure.

            EXAMPLES::

                sage: Algebras(QQ).CartesianProducts().extra_super_categories()
                [Category of algebras over Rational Field]
                sage: Algebras(QQ).CartesianProducts().super_categories()
                [Category of algebras over Rational Field, Category of Cartesian products of monoids]
            """
            return [self.base_category()]


    class TensorProducts(TensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Algebras(QQ).TensorProducts().extra_super_categories()
                [Category of algebras over Rational Field]
                sage: Algebras(QQ).TensorProducts().super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of algebras is an algebra
            """
            return [self.base_category()]

        class ParentMethods:
            #def coproduct(self):
            #    tensor products of morphisms are not yet implemented
            #    return tensor(module.coproduct for module in self.modules)
            pass

        class ElementMethods:
            pass

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Returns the dual category

            EXAMPLES:

            The category of algebras over the Rational Field is dual
            to the category of coalgebras over the same field::

                sage: C = Algebras(QQ)
                sage: C.dual()
                Category of duals of algebras over Rational Field
                sage: C.dual().extra_super_categories()
                [Category of coalgebras over Rational Field]
            """
            from sage.categories.coalgebras import Coalgebras
            return [Coalgebras(self.base_category().base_ring())]
