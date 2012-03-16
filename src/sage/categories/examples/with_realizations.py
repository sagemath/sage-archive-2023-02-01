r"""
Examples of parents endowed with multiple realizations
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.all import Rings, Algebras, AlgebrasWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.sets.set import Set
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.subset import Subsets

class SubsetAlgebra(UniqueRepresentation, Parent):
    r"""
    An example of parent endowed with several realizations

    We consider an algebra `A(S)` whose bases are indexed by the
    subsets `s` of a given set `S`. We consider three natural basis of
    this algebra: ``F``, ``In``, and ``Out``. In the first basis, the
    product is given by the union of the indexing sets. That is, for any
    `s, t\subset S`

    .. MATH::

        F_s F_t  = F_{s\cup t}

    The ``In`` basis and ``Out`` basis are defined respectively by:

    .. MATH::

        In_s  = \sum_{t\subset s} F_t
        \qquad\text{and}\qquad
        F_s   = \sum_{t\supset s} Out_t

    Each such basis gives a realization of `A`, where the elements are
    represented by their expansion in this basis.

    This parent, and its code, demonstrate how to implement this
    algebra and its three realizations, with coercions and mixed
    arithmetic between them.

    .. SEEALSO::

       - :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`

    EXAMPLES::

        sage: A = Sets().WithRealizations().example(); A
        The subset algebra of {1, 2, 3} over Rational Field

    The three bases of ``A``::

        sage: F   = A.F()  ; F
        The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
        sage: In  = A.In() ; In
        The subset algebra of {1, 2, 3} over Rational Field on the in basis
        sage: Out = A.Out(); Out
        The subset algebra of {1, 2, 3} over Rational Field on the out basis

        sage: F.basis().list()
        [F[{}], F[{1}], F[{2}], F[{3}], F[{1, 2}], F[{1, 3}], F[{2, 3}], F[{1, 2, 3}]]

    Some conversions::

        sage: F(In.from_set(2,3))
        F[{}] + F[{2}] + F[{3}] + F[{2, 3}]
        sage: In(F.from_set(2,3))
        In[{}] - In[{2}] - In[{3}] + In[{2, 3}]

        sage: Out(F.from_set(3))
        Out[{3}] + Out[{1, 3}] + Out[{2, 3}] + Out[{1, 2, 3}]
        sage: F(Out.from_set(3))
        F[{3}] - F[{1, 3}] - F[{2, 3}] + F[{1, 2, 3}]

        sage: Out(In.from_set(2,3))
        Out[{}] + Out[{1}] + 2*Out[{2}] + 2*Out[{3}] + 2*Out[{1, 2}] + 2*Out[{1, 3}] + 4*Out[{2, 3}] + 4*Out[{1, 2, 3}]

    We can now mix expressions::

        sage: (1 + Out.from_set(1)) * In.from_set(2,3)
        Out[{}] + 2*Out[{1}] + 2*Out[{2}] + 2*Out[{3}] + 2*Out[{1, 2}] + 2*Out[{1, 3}] + 4*Out[{2, 3}] + 4*Out[{1, 2, 3}]
    """

    def __init__(self, R, S):
        r"""
        EXAMPLES::

            sage: from sage.categories.examples.with_realizations import SubsetAlgebra
            sage: A = SubsetAlgebra(QQ, Set((1,2,3))); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: Sets().WithRealizations().example() is A
            True
            sage: TestSuite(A).run()
        """
        assert(R in Rings())
        self._base = R # Won't be needed when CategoryObject won't override anymore base_ring
        self._S = S
        Parent.__init__(self, category = Algebras(R).WithRealizations())

    # Could possibly go in Monoids.WithRealizations.ParentMethods
    def one(self):
        r"""
        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.one()
            F[{}]
            sage: A.one() is A.F().one()
            True
        """
        return self.F().one()

    # Could possibly go in CommutativeAdditiveMonoids.WithRealizations.ParentMethods
    def zero(self):
        r"""
        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.zero()
            0
            sage: A.zero() is A.F().zero()
            True
        """
        return self.F().zero()

    # Could be inherited from ParentWithBase?
    def base_ring(self):
        r"""
        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.base_ring()
            Rational Field
        """
        return self._base

    def base_set(self):
        r"""
        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.base_set()
            {1, 2, 3}
        """
        return self._S

    def indices(self):
        r"""
        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.indices()
            Subsets of {1, 2, 3}
        """
        return Subsets(self._S)

    def indices_cmp(self, x, y):
        r"""
        A comparison function on sets which gives a linear extension
        of the inclusion order.

        INPUT:

        - ``x``, ``y`` -- sets

        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: sorted(A.indices(), A.indices_cmp)
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
        """
        s = cmp(len(x), len(y))
        if s != 0:
            return s
        return cmp(list(x), list(y))

    def supsets(self, set):
        r"""
        INPUT:

        - ``set`` -- a subset of the base set `S` of ``self``

        Returns all the subsets of `S` containing ``set``

        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.supsets(Set((2,)))
            [{1, 2, 3}, {2, 3}, {1, 2}, {2}]
        """
        S = self.base_set()
        return list(S.difference(s) for s in Subsets(S.difference(set)))

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: Sets().WithRealizations().example()   # indirect doctest
            The subset algebra of {1, 2, 3} over Rational Field
        """
        return "The subset algebra of %s over %s"%(self.base_set(), self.base_ring())

    # Eventually it will be possible to put the class directly here if desired
    def F(self):
        """
        Returns the fundamental basis of ``self``

        EXAMPLES::

            sage: Sets().WithRealizations().example().F()
            The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
        """
        return Fundamental(self)

    def In(self):
        """
        Returns the in basis of ``self``

        EXAMPLES::

            sage: Sets().WithRealizations().example().In()
            The subset algebra of {1, 2, 3} over Rational Field on the in basis
        """
        return In(self)

    def Out(self):
        """
        Returns the out basis of ``self``

        EXAMPLES::

            sage: Sets().WithRealizations().example().Out()
            The subset algebra of {1, 2, 3} over Rational Field on the out basis
        """
        return Out(self)

    def __init_extra__(self):
        r"""
        Initializes the bases and change of bases of ``self``

        TESTS::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: F, In, Out = A.realizations()
            sage: type(F.coerce_map_from(In))
            <class 'sage.categories.modules_with_basis.TriangularModuleMorphism'>
            sage: type(In.coerce_map_from(F))
            <class 'sage.categories.modules_with_basis.TriangularModuleMorphism'>
            sage: type(F.coerce_map_from(Out))
            <class 'sage.categories.modules_with_basis.TriangularModuleMorphism'>
            sage: type(Out.coerce_map_from(F))
            <class 'sage.categories.modules_with_basis.TriangularModuleMorphism'>
            sage: In.coerce_map_from(Out)
            Composite map:
              From: The subset algebra of {1, 2, 3} over Rational Field on the out basis
              To:   The subset algebra of {1, 2, 3} over Rational Field on the in basis
              Defn:   Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field on the out basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
                    then
                      Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field on the in basis
            sage: Out.coerce_map_from(In)
            Composite map:
              From: The subset algebra of {1, 2, 3} over Rational Field on the in basis
              To:   The subset algebra of {1, 2, 3} over Rational Field on the out basis
              Defn:   Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field on the in basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
                    then
                      Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field on the out basis
        """
        category   = self.Realizations()
        F   = self.F()
        In  = self.In()
        Out = self.Out()

        In_to_F = In.module_morphism(F.sum_of_monomials * Subsets,
                                     codomain = F, category = category,
                                     triangular = 'upper', unitriangular = True,
                                     cmp = self.indices_cmp)
        In_to_F   .register_as_coercion()
        (~In_to_F).register_as_coercion()

        F_to_Out = F.module_morphism(Out.sum_of_monomials * self.supsets,
                                     codomain = Out, category = category,
                                     triangular = 'lower', unitriangular = True,
                                     cmp = self.indices_cmp)
        F_to_Out   .register_as_coercion()
        (~F_to_Out).register_as_coercion()

    # Alternatively, this category can be defined elsewhere, with just a link
    class Realizations(Category_realization_of_parent):
        r"""
        The category of the realizations of the subset algebra
        """

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: A = Sets().WithRealizations().example(); A
                The subset algebra of {1, 2, 3} over Rational Field
                sage: C = A.Realizations(); C
                The category of realizations of The subset algebra of {1, 2, 3} over Rational Field
                sage: C.super_categories()
                [Join of Category of algebras over Rational Field and Category of realizations of sets, Category of algebras with basis over Rational Field]
            """
            R = self.base().base_ring()
            return [Algebras(R).Realizations(), AlgebrasWithBasis(R)]

        class ParentMethods:

            def from_set(self, *args):
                r"""
                Construct the monomial indexed by the set containing the
                elements passed as arguments.

                EXAMPLES::

                    sage: In = Sets().WithRealizations().example().In(); In
                    The subset algebra of {1, 2, 3} over Rational Field on the in basis
                    sage: In.from_set(2,3)
                    In[{2, 3}]

                .. TODO:: rename to __getitem__ once ``Parent.__getitem__``
                   won't be anymore in the way::

                        sage: In[2,3]     # todo: not implemented
                        In[{2, 3}]
                """
                return self.monomial(Set(args))

            # This could go in some super category VectorSpaces().Realizations()
            def _repr_(self):
                r"""
                EXAMPLES::

                    sage: Sets().WithRealizations().example().In()  # indirect doctest
                    The subset algebra of {1, 2, 3} over Rational Field on the in basis
                """
                return "%s %s"%(self.realization_of(), "on the %s basis"%(self._realization_name()))

            @cached_method
            def one(self):
                r"""
                Returns the unit of this algebra.

                This default implementation takes the unit in the
                fundamental basis, and coerces it in ``self``.

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: In = A.In(); Out = A.Out()
                    sage: In.one()
                    In[{}]
                    sage: Out.one()
                    Out[{}] + Out[{1}] + Out[{2}] + Out[{3}] + Out[{1, 2}] + Out[{1, 3}] + Out[{2, 3}] + Out[{1, 2, 3}]
                """
                return self(self.realization_of().F().one())

class Fundamental(CombinatorialFreeModule):
    r"""
    The class for the parent modeling the fundamental basis of the
    Subset Algebra

    INPUT:

    - ``A`` -- a parent with realization in :class:`SubsetAlgebra`
    """

    def __init__(self, A):
        r"""
        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: F = A.F(); F
            The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
            sage: TestSuite(F).run()
        """
        CombinatorialFreeModule.__init__(self,
            A.base_ring(), A.indices(),
            category=A.Realizations(), prefix='F', monomial_cmp=A.indices_cmp)

    def product_on_basis(self, left, right):
        r"""
        Product of basis elements, as per :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis`.

        INPUT:

        - ``left``, ``right`` -- sets indexing basis elements

        EXAMPLES::

            sage: F = Sets().WithRealizations().example().F(); F
            The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
            sage: S = F.basis().keys(); S
            Subsets of {1, 2, 3}
            sage: F.product_on_basis(S([]), S([]))
            F[{}]
            sage: F.product_on_basis(S({1}), S({3}))
            F[{1, 3}]
            sage: F.product_on_basis(S({1,2}), S({2,3}))
            F[{1, 2, 3}]
        """
        return self.monomial( left.union(right)  )

    def one_basis(self):
        r"""
        Returns the index of the basis element which is equal to '1'.

        EXAMPLES::

            sage: F = Sets().WithRealizations().example().F(); F
            The subset algebra of {1, 2, 3} over Rational Field on the fundamental basis
            sage: F.one_basis()
            {}
            sage: F.one()
            F[{}]
        """
        return Set([])

    one = AlgebrasWithBasis.ParentMethods.one

class In(CombinatorialFreeModule):
    r"""
    The class for the parent modeling the ``In`` basis of the Subset
    Algebra

    INPUT:

    - ``A`` -- a parent with realization in :class:`SubsetAlgebra`
    """

    def __init__(self, A):
        r"""
        EXAMPLES::

            sage: In = Sets().WithRealizations().example().In(); In
            The subset algebra of {1, 2, 3} over Rational Field on the in basis
            sage: TestSuite(In).run()
        """
        CombinatorialFreeModule.__init__(self,
            A.base_ring(), A.indices(),
            category=A.Realizations(), prefix='In', monomial_cmp=A.indices_cmp)

    # This won't be needed any more once #8878 will be closed
    def product(self, x, y):
        r"""
        Returns the product of ``x`` by ``y``

        .. SEEALSO:: :meth:`Magmas.ParentMethods.product`

        EXAMPLES::

            sage: In = Sets().WithRealizations().example().In(); In
            The subset algebra of {1, 2, 3} over Rational Field on the in basis
            sage: x = In.an_element()
            sage: y = In.an_element()
            sage: In.product(x, y)
            -21*In[{}] - 2*In[{1}] + 19*In[{2}] + 53*In[{1, 2}]

        .. TODO:: this method won't be needed once :trac:`8878` will be closed
        """
        F = self.realization_of().F()
        return self(F(x) * F(y))

class Out(CombinatorialFreeModule):
    r"""
    The class for the parent modeling the `Out` basis of the Subset
    Algebra

    INPUT:

    - ``A`` -- a parent with realization in :class:`SubsetAlgebra`
    """

    def __init__(self, A):
        r"""
        EXAMPLES::

            sage: Out = Sets().WithRealizations().example().Out(); Out
            The subset algebra of {1, 2, 3} over Rational Field on the out basis
            sage: TestSuite(Out).run()
        """
        CombinatorialFreeModule.__init__(self,
            A.base_ring(), A.indices(),
            category=A.Realizations(), prefix='Out', monomial_cmp=A.indices_cmp)

    def product(self, x, y):
        r"""
        Returns the product of ``x`` by ``y``

        .. SEEALSO:: :meth:`Magmas.ParentMethods.product`

        EXAMPLES::

            sage: Out = Sets().WithRealizations().example().Out(); Out
            The subset algebra of {1, 2, 3} over Rational Field on the out basis
            sage: x = Out.an_element()
            sage: y = Out.an_element()
            sage: Out.product(x, y)
            Out[{}] + 4*Out[{1}] + 9*Out[{2}] + Out[{1, 2}]

        .. TODO:: this method won't be needed once :trac:`8878` will be closed
        """
        F = self.realization_of().F()
        return self(F(x) * F(y))
