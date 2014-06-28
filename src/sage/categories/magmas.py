r"""
Magmas
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.misc.abstract_method import abstract_method
from sage.categories.subquotients import SubquotientsCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_singleton
from sage.categories.sets_cat import Sets
from sage.categories.realizations import RealizationsCategory
from sage.structure.sage_object import have_same_parent

class Magmas(Category_singleton):
    """
    The category of (multiplicative) magmas.

    A magma is a set with a binary operation `*`.

    EXAMPLES::

        sage: Magmas()
        Category of magmas
        sage: Magmas().super_categories()
        [Category of sets]
        sage: Magmas().all_super_categories()
        [Category of magmas, Category of sets,
         Category of sets with partial maps, Category of objects]

    The following axioms are defined by this category::

        sage: Magmas().Associative()
        Category of semigroups
        sage: Magmas().Unital()
        Category of unital magmas
        sage: Magmas().Commutative()
        Category of commutative magmas
        sage: Magmas().Unital().Inverse()
        Category of inverse unital magmas
        sage: Magmas().Associative()
        Category of semigroups
        sage: Magmas().Associative().Unital()
        Category of monoids
        sage: Magmas().Associative().Unital().Inverse()
        Category of groups

    TESTS::

        sage: C = Magmas()
        sage: TestSuite(C).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: Magmas().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class SubcategoryMethods:

        @cached_method
        def Associative(self):
            """
            Return the full subcategory of the associative objects
            of ``self``.

            A (multiplicative) :class:`magma Magmas` `M` is
            *associative* if, for all `x,y,z\in M`,

            .. MATH:: x * (y * z) = (x * y) * z

            .. SEEALSO:: :wikipedia:`Associative_property`

            EXAMPLES::

                sage: Magmas().Associative()
                Category of semigroups

            TESTS::

                sage: TestSuite(Magmas().Associative()).run()
                sage: Rings().Associative.__module__
                'sage.categories.magmas'
            """
            return self._with_axiom('Associative')

        @cached_method
        def Commutative(self):
            """
            Return the full subcategory of the commutative objects
            of ``self``.

            A (multiplicative) :class:`magma Magmas` `M` is
            *commutative* if, for all `x,y\in M`,

            .. MATH:: x * y = y * x

            .. SEEALSO:: :wikipedia:`Commutative_property`

            EXAMPLES::

                sage: Magmas().Commutative()
                Category of commutative magmas
                sage: Monoids().Commutative()
                Category of commutative monoids

            TESTS::

                sage: TestSuite(Magmas().Commutative()).run()
                sage: Rings().Commutative.__module__
                'sage.categories.magmas'
            """
            return self._with_axiom('Commutative')

        @cached_method
        def Unital(self):
            r"""
            Return the subcategory of the unital objects of ``self``.

            A (multiplicative) :class:`magma Magmas` `M` is *unital*
            if it admits an element `1`, called *unit*, such that for
            all `x\in M`,

            .. MATH:: 1 * x = x * 1 = x

            This element is necessarily unique, and should be provided
            as ``M.one()``.

            .. SEEALSO:: :wikipedia:`Unital_magma#unital`

            EXAMPLES::

                sage: Magmas().Unital()
                Category of unital magmas
                sage: Semigroups().Unital()
                Category of monoids
                sage: Monoids().Unital()
                Category of monoids
                sage: from sage.categories.associative_algebras import AssociativeAlgebras
                sage: AssociativeAlgebras(QQ).Unital()
                Category of algebras over Rational Field

            TESTS::

                sage: TestSuite(Magmas().Unital()).run()
                sage: Semigroups().Unital.__module__
                'sage.categories.magmas'
            """
            return self._with_axiom("Unital")

        @cached_method
        def Distributive(self):
            """
            Return the full subcategory of the objects of ``self``
            where `*` is distributive on `+`.

            INPUT:

            - ``self`` -- a subcategory of :class:`Magmas`
              and :class:`AdditiveMagmas`

            Given that Sage does not yet know that the category
            :class:`MagmasAndAdditiveMagmas` is the intersection of
            the categories :class:`Magmas` and
            :class:`AdditiveMagmas`, the method
            :meth:`MagmasAndAdditiveMagmas.SubcategoryMethods.Distributive`
            is not available, as would be desirable, for this intersection.

            This method is a workaround. It checks that ``self`` is a
            subcategory of both :class:`Magmas` and
            :class:`AdditiveMagmas` and upgrades it to a subcategory
            of :class:`MagmasAndAdditiveMagmas` before applying the
            axiom. It complains overwise, since the ``Distributive``
            axiom does not make sense for a plain magma.

            EXAMPLES::

                sage: (Magmas() & AdditiveMagmas()).Distributive()
                Category of distributive magmas and additive magmas
                sage: (Monoids() & CommutativeAdditiveGroups()).Distributive()
                Category of rings

                sage: Magmas().Distributive()
                Traceback (most recent call last):
                ...
                ValueError: The distributive axiom only makes sense on a magma which is simultaneously an additive magma
                sage: Semigroups().Distributive()
                Traceback (most recent call last):
                ...
                ValueError: The distributive axiom only makes sense on a magma which is simultaneously an additive magma

            TESTS::

                sage: Semigroups().Distributive.__module__
                'sage.categories.magmas'
                sage: Rings().Distributive.__module__
                'sage.categories.magmas_and_additive_magmas'
            """
            from additive_magmas import AdditiveMagmas
            if not self.is_subcategory(AdditiveMagmas()):
                raise ValueError("The distributive axiom only makes sense on a magma which is simultaneously an additive magma")
            from magmas_and_additive_magmas import MagmasAndAdditiveMagmas
            return (self & MagmasAndAdditiveMagmas()).Distributive()

    Associative = LazyImport('sage.categories.semigroups', 'Semigroups', at_startup=True)

    class Algebras(AlgebrasCategory):

        def extra_super_categories(self):
            """
            EXAMPLES:

                sage: Magmas().Commutative().Algebras(QQ).extra_super_categories()
                [Category of commutative magmas]

            This implements the fact that the algebra of a commutative
            magma is commutative::

                sage: Magmas().Commutative().Algebras(QQ).super_categories()
                [Category of magma algebras over Rational Field, Category of commutative magmas]

            In particular, commutative monoid algebras are
            commutative algebras::

                sage: Monoids().Commutative().Algebras(QQ).is_subcategory(Algebras(QQ).Commutative())
                True
            """
            from sage.categories.magmatic_algebras import MagmaticAlgebras
            return [MagmaticAlgebras(self.base_ring())]

    class Commutative(CategoryWithAxiom):

        class ParentMethods:
            def is_commutative(self):
                """
                Return ``True``, since commutative magmas are commutative.

                EXAMPLES::

                    sage: Parent(QQ,category=CommutativeRings()).is_commutative()
                    True
                """
                return True

        class Algebras(AlgebrasCategory):

            def extra_super_categories(self):
                """
                EXAMPLES:

                    sage: Magmas().Commutative().Algebras(QQ).extra_super_categories()
                    [Category of commutative magmas]

                This implements the fact that the algebra of a commutative
                magma is commutative::

                    sage: Magmas().Commutative().Algebras(QQ).super_categories()
                    [Category of magma algebras over Rational Field,
                     Category of commutative magmas]

                In particular, commutative monoid algebras are
                commutative algebras::

                    sage: Monoids().Commutative().Algebras(QQ).is_subcategory(Algebras(QQ).Commutative())
                    True
                """
                return [Magmas().Commutative()]

    class Unital(CategoryWithAxiom):

        class ParentMethods:
            @cached_method
            def one(self):
                r"""
                Return the unit of the monoid, that is the unique neutral
                element for `*`.

                .. NOTE::

                   The default implementation is to coerce `1` into ``self``.
                   It is recommended to override this method because the
                   coercion from the integers:

                    - is not always meaningful (except for `1`);
                    - often uses ``self.one()``.

                EXAMPLES::

                    sage: M = Monoids().example(); M
                    An example of a monoid: the free monoid generated by ('a', 'b', 'c', 'd')
                    sage: M.one()
                    ''
                """
                return self(1)

            def _test_one(self, **options):
                r"""
                Test that ``self.one()`` is an element of ``self`` and is
                neutral for the operation ``*``.

                INPUT:

                - ``options`` -- any keyword arguments accepted by :meth:`_tester`

                EXAMPLES:

                By default, this method tests only the elements returned by
                ``self.some_elements()``::

                    sage: S = Monoids().example()
                    sage: S._test_one()

                However, the elements tested can be customized with the
                ``elements`` keyword argument::

                    sage: S._test_one(elements = (S('a'), S('b')))

                See the documentation for :class:`TestSuite` for more information.
                """
                tester = self._tester(**options)
                one = self.one()
                tester.assert_(self.is_parent_of(one))
                for x in tester.some_elements():
                    tester.assert_(x * one == x)
                    tester.assert_(one * x == x)
                # Check that one is immutable by asking its hash;
                tester.assertEqual(type(one.__hash__()), int)
                tester.assertEqual(one.__hash__(), one.__hash__())


        class SubcategoryMethods:

            @cached_method
            def Inverse(self):
                r"""
                Return the full subcategory of the inverse objects of ``self``.

                An inverse :class:` (multiplicative) magma <Magmas>`
                is a :class:`unital magma <Magmas.Unital>` such that
                every element admits both an inverse on the left and
                on the right. Such a magma is also called a *loop*.

                .. SEEALSO::

                    :wikipedia:`Inverse_element`, :wikipedia:`Quasigroup`

                EXAMPLES::

                    sage: Magmas().Unital().Inverse()
                    Category of inverse unital magmas
                    sage: Monoids().Inverse()
                    Category of groups

                TESTS::

                    sage: TestSuite(Magmas().Unital().Inverse()).run()
                    sage: Algebras(QQ).Inverse.__module__
                    'sage.categories.magmas'
                """
                return self._with_axiom("Inverse")

        class Inverse(CategoryWithAxiom):
            class CartesianProducts(CartesianProductsCategory):
                def extra_super_categories(self):
                    """
                    Implement the fact that a cartesian product of magmas with
                    inverses is a magma with inverse.

                    EXAMPLES::

                        sage: C = Magmas().Unital().Inverse().CartesianProducts()
                        sage: C.extra_super_categories();
                        [Category of inverse unital magmas]
                        sage: sorted(C.axioms())
                        ['Inverse', 'Unital']
                    """
                    return [Magmas().Unital().Inverse()]

        class CartesianProducts(CartesianProductsCategory):
            def extra_super_categories(self):
                """
                Implement the fact that a cartesian product of unital magmas is
                a unital magma

                EXAMPLES::

                    sage: C = Magmas().Unital().CartesianProducts()
                    sage: C.extra_super_categories();
                    [Category of unital magmas]
                    sage: C.axioms()
                    frozenset(['Unital'])

                    sage: Monoids().CartesianProducts().is_subcategory(Monoids())
                    True
                """
                return [Magmas().Unital()]

            class ParentMethods:

                @cached_method
                def one(self):
                    """
                    Return the unit of this cartesian product.

                    It is built from the units for the cartesian factors of ``self``.

                    EXAMPLES::

                        sage: cartesian_product([QQ, ZZ, RR]).one()
                        (1, 1, 1.00000000000000)
                    """
                    return self._cartesian_product_of_elements(
                        _.one() for _ in self.cartesian_factors())

            class ElementMethods:
                def __invert__(self):
                    r"""
                    Return the inverse of ``self``, if it exists.

                    The inverse is computed by inverting each
                    cartesian factor and attempting to convert the
                    result back to the original parent.

                    For example, if one of the cartesian factor is an
                    element ``x`` of `\ZZ`, the result of ``~x`` is in
                    `\QQ`. So we need to convert it back to `\ZZ`. As
                    a side effect, this checks that ``x`` is indeed
                    invertible in `\ZZ`.

                    If needed an optimized version without this
                    conversion could be implemented in
                    :class:`Magmas.Unital.Inverse.ElementMethods`.

                    EXAMPLES::

                        sage: C = cartesian_product([QQ, ZZ, RR, GF(5)])
                        sage: c = C([2,-1,2,2]); c
                        (2, -1, 2.00000000000000, 2)
                        sage: ~c
                        (1/2, -1, 0.500000000000000, 3)

                    This fails as soon as one of the entries is not
                    invertible::

                        sage: ~C([0,2,2,2])
                        Traceback (most recent call last):
                        ...
                        ZeroDivisionError: rational division by zero

                        sage: ~C([2,2,2,2])
                        Traceback (most recent call last):
                        ...
                        TypeError: no conversion of this rational to integer
                    """
                    # variant without coercion:
                    # return self.parent()._cartesian_product_of_elements(
                    return self.parent()(
                        ~x for x in self.cartesian_factors())

        class Algebras(AlgebrasCategory):

            def extra_super_categories(self):
                """
                EXAMPLES:

                    sage: Magmas().Commutative().Algebras(QQ).extra_super_categories()
                    [Category of commutative magmas]

                This implements the fact that the algebra of a
                commutative magma is commutative::

                    sage: Magmas().Commutative().Algebras(QQ).super_categories()
                    [Category of magma algebras over Rational Field,
                     Category of commutative magmas]

                In particular, commutative monoid algebras are
                commutative algebras::

                    sage: Monoids().Commutative().Algebras(QQ).is_subcategory(Algebras(QQ).Commutative())
                    True
                """
                return [Magmas().Unital()]

    class ParentMethods:

        def product(self, x, y):
            """
            The binary multiplication of the magma.

            INPUT:

            - ``x``, ``y`` -- elements of this magma

            OUTPUT:

            - an element of the magma (the product of ``x`` and ``y``)

            EXAMPLES::

                sage: S = Semigroups().example("free")
                sage: x = S('a'); y = S('b')
                sage: S.product(x, y)
                'ab'

            A parent in ``Magmas()`` must either implement
            :meth:`.product` in the parent class or ``_mul_`` in the
            element class. By default, the addition method on elements
            ``x._mul_(y)`` calls ``S.product(x,y)``, and reciprocally.

            As a bonus, ``S.product`` models the binary function from
            ``S`` to ``S``::

                sage: bin = S.product
                sage: bin(x,y)
                'ab'

            Currently, ``S.product`` is just a bound method::

                sage: bin
                <bound method FreeSemigroup_with_category.product of An example of a semigroup: the free semigroup generated by ('a', 'b', 'c', 'd')>

            When Sage will support multivariate morphisms, it will be
            possible, and in fact recommended, to enrich ``S.product``
            with extra mathematical structure. This will typically be
            implemented using lazy attributes.::

                sage: bin                 # todo: not implemented
                Generic binary morphism:
                From: (S x S)
                To:   S
            """
            return x._mul_(y)

        product_from_element_class_mul = product

        def __init_extra__(self):
            """
                sage: S = Semigroups().example("free")
                sage: S('a') * S('b') # indirect doctest
                'ab'
                sage: S('a').__class__._mul_ == S('a').__class__._mul_parent
                True
            """
            # This should instead register the multiplication to the coercion model
            # But this is not yet implemented in the coercion model
            #
            # Trac ticket #11900: The following used to test whether
            # self.product != self.product_from_element_class_mul. But
            # that is, of course, a bug. Namely otherwise, if the parent
            # has an optimized `product` then its elements will *always* use
            # a slow generic `_mul_`.
            #
            # So, in addition, it should be tested whether the element class exists
            # *and* has a custom _mul_, because in this case it must not be overridden.
            if (self.product.__func__ == self.product_from_element_class_mul.__func__):
                return
            if not (hasattr(self, "element_class") and hasattr(self.element_class, "_mul_parent")):
                return
            E = self.element_class
            if hasattr(E._mul_,'__func__'):
                try:
                    el_class_mul = self.category().element_class._mul_.__func__
                except AttributeError: # abstract method
                    return
                if E._mul_.__func__ is el_class_mul:
                    # self.product is custom, thus, we rely on it
                    E._mul_ = E._mul_parent
            else: # E._mul_ has so far been abstract
                E._mul_ = E._mul_parent

        def multiplication_table(self, names='letters', elements=None):
            r"""
            Returns a table describing the multiplication operation.

            .. note:: The order of the elements in the row and column
              headings is equal to the order given by the table's
              :meth:`~sage.matrix.operation_table.OperationTable.list`
              method.  The association can also be retrieved with the
              :meth:`~sage.matrix.operation_table.OperationTable.dict`
              method.

            INPUTS:

            - ``names`` - the type of names used

              * ``'letters'`` - lowercase ASCII letters are used
                for a base 26 representation of the elements'
                positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading 'a's.
              * ``'digits'`` - base 10 representation of the
                elements' positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading zeros.
              * ``'elements'`` - the string representations
                of the elements themselves.
              * a list - a list of strings, where the length
                of the list equals the number of elements.
            - ``elements`` - default = ``None``.  A list of
              elements of the magma, in forms that can be
              coerced into the structure, eg. their string
              representations. This may be used to impose an
              alternate ordering on the elements, perhaps when
              this is used in the context of a particular structure.
              The default is to use whatever ordering the ``S.list``
              method returns. Or the ``elements`` can be a subset
              which is closed under the operation. In particular,
              this can be used when the base set is infinite.

            OUTPUT:
            The multiplication table as an object of the class
            :class:`~sage.matrix.operation_table.OperationTable`
            which defines several methods for manipulating and
            displaying the table.  See the documentation there
            for full details to supplement the documentation
            here.

            EXAMPLES:

            The default is to represent elements as lowercase
            ASCII letters.  ::

                sage: G=CyclicPermutationGroup(5)
                sage: G.multiplication_table()
                *  a b c d e
                 +----------
                a| a b c d e
                b| b c d e a
                c| c d e a b
                d| d e a b c
                e| e a b c d

            All that is required is that an algebraic structure
            has a multiplication defined.  A
            :class:`~sage.categories.examples.finite_semigroups.LeftRegularBand`
            is an example of a finite semigroup.  The ``names`` argument allows
            displaying the elements in different ways.  ::

                sage: from sage.categories.examples.finite_semigroups import LeftRegularBand
                sage: L=LeftRegularBand(('a','b'))
                sage: T=L.multiplication_table(names='digits')
                sage: T.column_keys()
                ('a', 'b', 'ab', 'ba')
                sage: T
                *  0 1 2 3
                 +--------
                0| 0 2 2 2
                1| 3 1 3 3
                2| 2 2 2 2
                3| 3 3 3 3

            Specifying the elements in an alternative order can provide
            more insight into how the operation behaves.  ::

                sage: L=LeftRegularBand(('a','b','c'))
                sage: elts = sorted(L.list())
                sage: L.multiplication_table(elements=elts)
                *  a b c d e f g h i j k l m n o
                 +------------------------------
                a| a b c d e b b c c c d d e e e
                b| b b c c c b b c c c c c c c c
                c| c c c c c c c c c c c c c c c
                d| d e e d e e e e e e d d e e e
                e| e e e e e e e e e e e e e e e
                f| g g h h h f g h i j i j j i j
                g| g g h h h g g h h h h h h h h
                h| h h h h h h h h h h h h h h h
                i| j j j j j i j j i j i j j i j
                j| j j j j j j j j j j j j j j j
                k| l m m l m n o o n o k l m n o
                l| l m m l m m m m m m l l m m m
                m| m m m m m m m m m m m m m m m
                n| o o o o o n o o n o n o o n o
                o| o o o o o o o o o o o o o o o

            The ``elements`` argument can be used to provide
            a subset of the elements of the structure.  The subset
            must be closed under the operation.  Elements need only
            be in a form that can be coerced into the set.  The
            ``names`` argument can also be used to request that
            the elements be represented with their usual string
            representation.  ::

                sage: L=LeftRegularBand(('a','b','c'))
                sage: elts=['a', 'c', 'ac', 'ca']
                sage: L.multiplication_table(names='elements', elements=elts)
                   *   'a'  'c' 'ac' 'ca'
                    +--------------------
                 'a'|  'a' 'ac' 'ac' 'ac'
                 'c'| 'ca'  'c' 'ca' 'ca'
                'ac'| 'ac' 'ac' 'ac' 'ac'
                'ca'| 'ca' 'ca' 'ca' 'ca'

            The table returned can be manipulated in various ways.  See
            the documentation for
            :class:`~sage.matrix.operation_table.OperationTable` for more
            comprehensive documentation. ::

                sage: G=AlternatingGroup(3)
                sage: T=G.multiplication_table()
                sage: T.column_keys()
                ((), (1,2,3), (1,3,2))
                sage: sorted(T.translation().items())
                [('a', ()), ('b', (1,2,3)), ('c', (1,3,2))]
                sage: T.change_names(['x', 'y', 'z'])
                sage: sorted(T.translation().items())
                [('x', ()), ('y', (1,2,3)), ('z', (1,3,2))]
                sage: T
                *  x y z
                 +------
                x| x y z
                y| y z x
                z| z x y
            """
            from sage.matrix.operation_table import OperationTable
            import operator
            return OperationTable(self, operation=operator.mul, names=names, elements=elements)

    class ElementMethods:

        def __mul__(self, right):
            r"""
            Product of two elements

            INPUT::

             - ``self``, ``right`` -- two elements

            This calls the `_mul_` method of ``self``, if it is
            available and the two elements have the same parent.

            Otherwise, the job is delegated to the coercion model.

            Do not override; instead implement a ``_mul_`` method in the
            element class or a ``product`` method in the parent class.

            EXAMPLES::

                sage: S = Semigroups().example("free")
                sage: x = S('a'); y = S('b')
                sage: x * y
                'ab'
            """
            if have_same_parent(self, right) and hasattr(self, "_mul_"):
                return self._mul_(right)
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(self, right, operator.mul)

        __imul__ = __mul__

        @abstract_method(optional = True)
        def _mul_(self, right):
            """
            Product of two elements

            INPUT::

             - ``self``, ``right`` -- two elements with the same parent

            OUTPUT::

             - an element of the same parent

            EXAMPLES::

                sage: S = Semigroups().example("free")
                sage: x = S('a'); y = S('b')
                sage: x._mul_(y)
                'ab'
            """

        def _mul_parent(self, other):
            r"""
            Returns the product of the two elements, calculated using
            the ``product`` method of the parent.

            This is the default implementation of _mul_ if
            ``product`` is implemented in the parent.

            INPUT::

             - ``other`` -- an element of the parent of ``self``

            OUTPUT::

             - an element of the parent of ``self``

            EXAMPLES::

                sage: S = Semigroups().example("free")
                sage: x = S('a'); y = S('b')
                sage: x._mul_parent(y)
                'ab'

            """
            return self.parent().product(self, other)

        def is_idempotent(self):
            r"""
            Test whether ``self`` is idempotent.

            EXAMPLES::

                sage: S = Semigroups().example("free"); S
                An example of a semigroup: the free semigroup generated by ('a', 'b', 'c', 'd')
                sage: a = S('a')
                sage: a^2
                'aa'
                sage: a.is_idempotent()
                False

            ::

                sage: L = Semigroups().example("leftzero"); L
                An example of a semigroup: the left zero semigroup
                sage: x = L('x')
                sage: x^2
                'x'
                sage: x.is_idempotent()
                True

            """
            return self * self == self

    class CartesianProducts(CartesianProductsCategory):

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Semigroups().CartesianProducts().extra_super_categories()
                [Category of semigroups]
                sage: Semigroups().CartesianProducts().super_categories()
                [Category of semigroups, Category of Cartesian products of magmas]
            """
            return [Magmas()]

        def example(self):
            """
            Returns an example of cartesian product of magmas

            EXAMPLES::

                sage: C = Magmas().CartesianProducts().example(); C
                The cartesian product of (Rational Field, Integer Ring, Integer Ring)
                sage: C.category()
                Join of Category of rings ...
                sage: sorted(C.category().axioms())
                ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital', 'Associative', 'Distributive', 'Unital']

                sage: TestSuite(C).run()
            """
            from cartesian_product import cartesian_product
            from sage.rings.integer_ring import ZZ
            from sage.rings.rational_field import QQ
            return cartesian_product([QQ, ZZ, ZZ])

        class ParentMethods:

            def product(self, left, right):
                """
                EXAMPLES::

                    sage: C = Magmas().CartesianProducts().example(); C
                    The cartesian product of (Rational Field, Integer Ring, Integer Ring)
                    sage: x = C.an_element(); x
                    (1/2, 1, 1)
                    sage: x * x
                    (1/4, 1, 1)

                    sage: A = SymmetricGroupAlgebra(QQ, 3);
                    sage: x = cartesian_product([A([1,3,2]), A([2,3,1])])
                    sage: y = cartesian_product([A([1,3,2]), A([2,3,1])])
                    sage: cartesian_product([A,A]).product(x,y)
                    B[(0, [1, 2, 3])] + B[(1, [3, 1, 2])]
                    sage: x*y
                    B[(0, [1, 2, 3])] + B[(1, [3, 1, 2])]
                """
                return self._cartesian_product_of_elements([(a*b) for (a,b) in zip(left.cartesian_factors(), right.cartesian_factors())])

    class Subquotients(SubquotientsCategory):
        r"""
        The category of subquotient magmas.

        See :meth:`Sets.SubcategoryMethods.Subquotients` for the
        general setup for subquotients. In the case of a subquotient
        magma `S` of a magma `G`, the condition that `r` be a
        morphism in ``As`` can be rewritten as follows:

         - for any two `a,b \in S` the identity
           `a \times_S b = r(l(a) \times_G l(b))` holds.

        This is used by this category to implement the product
        `\times_S` of `S` from `l` and `r` and the product of `G`.

        EXAMPLES::

            sage: Semigroups().Subquotients().all_super_categories()
            [Category of subquotients of semigroups, Category of semigroups,
             Category of subquotients of magmas, Category of magmas,
             Category of subquotients of sets, Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """

        class ParentMethods:

            def product(self, x, y):
                """
                Return the product of two elements of ``self``.

                EXAMPLES::

                    sage: S = Semigroups().Subquotients().example()
                    sage: S.product(S(19), S(3))
                    19
                """
                assert(x in self)
                assert(y in self)
                return self.retract(self.lift(x) * self.lift(y))

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def product_by_coercion(self, left, right):
                r"""
                Default implementation of product for realizations.

                This method coerces to the realization specified by
                ``self.realization_of().a_realization()``, computes
                the product in that realization, and then coerces
                back.

                EXAMPLES::

                    sage: Out = Sets().WithRealizations().example().Out(); Out
                    The subset algebra of {1, 2, 3} over Rational Field in the Out basis
                    sage: Out.product
                    <bound method SubsetAlgebra.Out_with_category.product_by_coercion of The subset algebra of {1, 2, 3} over Rational Field in the Out basis>
                    sage: Out.product.__module__
                    'sage.categories.magmas'
                    sage: x = Out.an_element()
                    sage: y = Out.an_element()
                    sage: Out.product(x, y)
                    Out[{}] + 4*Out[{1}] + 9*Out[{2}] + Out[{1, 2}]

                """
                R = self.realization_of().a_realization()
                return self(R(left) * R(right))
