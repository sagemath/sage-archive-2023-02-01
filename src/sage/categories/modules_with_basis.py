r"""
Modules With Basis

AUTHORS:

- Nicolas M. Thiery (2008-2014): initial revision, axiomatization
- Jason Bandlow and Florent Hivert (2010): Triangular Morphisms
- Christian Stump (2010): :trac:`9648` module_morphism's to a wider class
  of codomains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method
from sage.misc.misc import attrcall
from sage.misc.superseded import deprecated_function_alias
from sage.misc.sage_itertools import max_cmp, min_cmp
from sage.categories.category import HomCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.tensor import tensor, TensorProductsCategory
from sage.categories.dual import DualObjectsCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.morphism import SetMorphism, Morphism
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
from sage.categories.commutative_additive_semigroups import CommutativeAdditiveSemigroups
from sage.categories.modules import Modules
from sage.structure.element import Element, parent

class ModulesWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of modules with a distinguished basis.

    The elements are represented by expanding them in the distinguished basis.
    The morphisms are not required to respect the distinguished basis.

    EXAMPLES::

        sage: ModulesWithBasis(ZZ)
        Category of modules with basis over Integer Ring
        sage: ModulesWithBasis(ZZ).super_categories()
        [Category of modules over Integer Ring]

    If the base ring is actually a field, this constructs instead the
    category of vector spaces with basis::

        sage: ModulesWithBasis(QQ)
        Category of vector spaces with basis over Rational Field

        sage: ModulesWithBasis(QQ).super_categories()
        [Category of modules with basis over Rational Field,
         Category of vector spaces over Rational Field]

    Let `X` and `Y` be two modules with basis. We can build `Hom(X,Y)`::

        sage: X = CombinatorialFreeModule(QQ, [1,2]); X.__custom_name = "X"
        sage: Y = CombinatorialFreeModule(QQ, [3,4]); Y.__custom_name = "Y"
        sage: H = Hom(X, Y); H
        Set of Morphisms from X to Y in Category of vector spaces with basis over Rational Field

    The simplest morphism is the zero map::

        sage: H.zero()         # todo: move this test into module once we have an example
        Generic morphism:
          From: X
          To:   Y

    which we can apply to elements of `X`::

        sage: x = X.monomial(1) + 3 * X.monomial(2)
        sage: H.zero()(x)
        0

    TESTS::

        sage: f = H.zero().on_basis()
        sage: f(1)
        0
        sage: f(2)
        0

    EXAMPLES:

    We now construct a more interesting morphism by extending a
    function by linearity::

        sage: phi = H(on_basis = lambda i: Y.monomial(i+2)); phi
        Generic morphism:
          From: X
          To:   Y
        sage: phi(x)
        B[3] + 3*B[4]

    We can retrieve the function acting on indices of the basis::

        sage: f = phi.on_basis()
        sage: f(1), f(2)
        (B[3], B[4])

    `Hom(X,Y)` has a natural module structure (except for the zero,
    the operations are not yet implemented though). However since the
    dimension is not necessarily finite, it is not a module with
    basis; but see :class:`FiniteDimensionalModulesWithBasis` and
    :class:`GradedModulesWithBasis`::

        sage: H in ModulesWithBasis(QQ), H in Modules(QQ)
        (False, True)

    Some more playing around with categories and higher order homsets::

        sage: H.category()
        Category of hom sets in Category of modules with basis over Rational Field
        sage: Hom(H, H).category()
        Category of hom sets in Category of modules over Rational Field

    .. TODO:: ``End(X)`` is an algebra.

    TESTS::

        sage: TestSuite(ModulesWithBasis(ZZ)).run()
    """

    def _call_(self, x):
        """
        Construct a module with basis (resp. vector space) from the data in ``x``.

        EXAMPLES::

            sage: CZ = ModulesWithBasis(ZZ); CZ
            Category of modules with basis over Integer Ring
            sage: CQ = ModulesWithBasis(QQ); CQ
            Category of vector spaces with basis over Rational Field

        ``x`` is returned unchanged if it is already in this category::

            sage: CZ(CombinatorialFreeModule(ZZ, ('a','b','c')))
            Free module generated by {'a', 'b', 'c'} over Integer Ring
            sage: CZ(ZZ^3)
            Ambient free module of rank 3 over the principal ideal domain Integer Ring

        If needed (and possible) the base ring is changed appropriately::

            sage: CQ(ZZ^3)                       # indirect doctest
            Vector space of dimension 3 over Rational Field

        If ``x`` itself is not a module with basis, but there is a
        canonical one associated to it, the later is returned::

            sage: CQ(AbelianVariety(Gamma0(37))) # indirect doctest
            Vector space of dimension 4 over Rational Field
        """
        try:
            M = x.free_module()
            if M.base_ring() != self.base_ring():
                M = M.change_ring(self.base_ring())
        except (TypeError, AttributeError) as msg:
            raise TypeError("%s\nunable to coerce x (=%s) into %s"%(msg,x,self))
        return M

    def is_abelian(self):
        """
        Returns whether this category is abelian.

        This is the case if and only if the base ring is a field.

        EXAMPLES::

            sage: ModulesWithBasis(QQ).is_abelian()
            True
            sage: ModulesWithBasis(ZZ).is_abelian()
            False
        """
        return self.base_ring().is_field()

    FiniteDimensional = LazyImport('sage.categories.finite_dimensional_modules_with_basis', 'FiniteDimensionalModulesWithBasis')
    Graded = LazyImport('sage.categories.graded_modules_with_basis', 'GradedModulesWithBasis')

    class ParentMethods:
        @cached_method
        def basis(self):
            """
            Return the basis of ``self``.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
                sage: F.basis()
                Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}

            ::

                sage: QS3 = SymmetricGroupAlgebra(QQ,3)
                sage: list(QS3.basis())
                [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
            """
            from sage.combinat.family import Family
            return Family(self._basis_keys, self.monomial)

        def module_morphism(self, on_basis = None, diagonal = None, triangular = None, **keywords):
            r"""
            Construct a morphism from ``self`` to ``codomain`` by
            linearity from its restriction ``on_basis`` to the basis of
            ``self``.

            Let ``self`` be the module `X` with a basis indexed by `I`.
            This constructs a morphism `f: X \to Y` by linearity from
            a map `I \to Y` which is to be its restriction to the
            basis `(x_i)_{i \in I}` of `X`.

            INPUT:

            - ``codomain`` -- the codomain `Y` of `f`: defaults to
              ``f.codomain()`` if the latter is defined
            - ``zero`` -- the zero of the codomain; defaults to
              ``codomain.zero()``; can be used (with care) to define affine maps
            - ``position`` -- a non-negative integer; defaults to 0
            - ``on_basis`` -- a function `f` which accepts elements of `I`
              (the indexing set of the basis of `X`) as ``position``-th argument
              and returns elements of `Y`
            - ``diagonal`` -- a function `d` from `I` to `R` (the base ring
              of ``self`` and ``codomain``)
            - ``triangular`` --  (default: ``None``) ``"upper"`` or
              ``"lower"`` or ``None``:

              * ``"upper"`` - if the :meth:`leading_support()` of the image
                of the basis vector `x_i` is `i`, or
              * ``"lower"`` - if the :meth:`trailing_support()` of the image
                of the basis vector `x_i` is `i`

            - ``category`` -- a category; by default, this is
              ``ModulesWithBasis(R)`` if `Y` is in this category, and
              otherwise this lets `Hom(X,Y)` decide

            Exactly one of ``on_basis`` and ``diagonal`` options should
            be specified.

            With the ``on_basis`` option, this returns a function `g`
            obtained by extending `f` by linearity on the ``position``-th
            positional argument. For example, for ``position == 1`` and a
            ternary function `f`, one has:

            .. MATH::

                g\left( a,\ \sum_i \lambda_i x_i,\ c \right)
                = \sum_i \lambda_i f(a, i, c).

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1,2,3]);   X.rename("X")
                sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4]); Y.rename("Y")
                sage: phi = X.module_morphism(lambda i: Y.monomial(i) + 2*Y.monomial(i+1), codomain = Y)
                sage: phi
                Generic morphism:
                From: X
                To:   Y
                sage: phi.category_for()
                Category of vector spaces with basis over Rational Field
                sage: x = X.basis(); y = Y.basis()
                sage: phi(x[1] + x[3])
                B[1] + 2*B[2] + B[3] + 2*B[4]

            With the ``zero`` argument, one can define affine morphisms::

                sage: phi = X.module_morphism(lambda i: Y.monomial(i) + 2*Y.monomial(i+1), codomain = Y, zero = 10*y[1])
                sage: phi(x[1] + x[3])
                11*B[1] + 2*B[2] + B[3] + 2*B[4]
                sage: phi.category_for()
                Category of sets

            One can construct morphisms with the base ring as codomain::

                sage: X = CombinatorialFreeModule(ZZ,[1,-1])
                sage: phi = X.module_morphism( on_basis=lambda i: i, codomain=ZZ )
                sage: phi( 2 * X.monomial(1) + 3 * X.monomial(-1) )
                -1
                sage: phi.category_for()
                Category of commutative additive semigroups
                sage: phi.category_for() # todo: not implemented (ZZ is currently not in Modules(ZZ))
                Category of modules over Integer Ring

            Or more generaly any ring admitting a coercion map from the base ring::

                sage: phi = X.module_morphism(on_basis= lambda i: i, codomain=RR )
                sage: phi( 2 * X.monomial(1) + 3 * X.monomial(-1) )
                -1.00000000000000
                sage: phi.category_for()
                Category of commutative additive semigroups
                sage: phi.category_for() # todo: not implemented (RR is currently not in Modules(ZZ))
                Category of modules over Integer Ring

                sage: phi = X.module_morphism(on_basis= lambda i: i, codomain=Zmod(4) )
                sage: phi( 2 * X.monomial(1) + 3 * X.monomial(-1) )
                3

                sage: phi = Y.module_morphism(on_basis= lambda i: i, codomain=Zmod(4) )
                Traceback (most recent call last):
                ...
                ValueError: codomain(=Ring of integers modulo 4) should be a module over the base ring of the domain(=Y)

            On can also define module morphisms between free modules
            over different base rings; here we implement the natural
            map from `X = \RR^2` to `Y = \CC`::

                sage: X = CombinatorialFreeModule(RR,['x','y'])
                sage: Y = CombinatorialFreeModule(CC,['z'])
                sage: x = X.monomial('x')
                sage: y = X.monomial('y')
                sage: z = Y.monomial('z')
                sage: def on_basis( a ):
                ....:     if a == 'x':
                ....:         return CC(1) * z
                ....:     elif a == 'y':
                ....:         return CC(I) * z
                sage: phi = X.module_morphism( on_basis=on_basis, codomain=Y )
                sage: v = 3 * x + 2 * y; v
                3.00000000000000*B['x'] + 2.00000000000000*B['y']
                sage: phi(v)
                (3.00000000000000+2.00000000000000*I)*B['z']
                sage: phi.category_for()
                Category of commutative additive semigroups
                sage: phi.category_for() # todo: not implemented (CC is currently not in Modules(RR)!)
                Category of vector spaces over Real Field with 53 bits of precision

                sage: Y = CombinatorialFreeModule(CC['q'],['z'])
                sage: z = Y.monomial('z')
                sage: phi = X.module_morphism( on_basis=on_basis, codomain=Y )
                sage: phi(v)
                (3.00000000000000+2.00000000000000*I)*B['z']

            Of course, there should be a coercion between the
            respective base rings of the domain and the codomain for
            this to be meaningful::

                sage: Y = CombinatorialFreeModule(QQ,['z'])
                sage: phi = X.module_morphism( on_basis=on_basis, codomain=Y )
                Traceback (most recent call last):
                ...
                ValueError: codomain(=Free module generated by {'z'} over Rational Field) should be a module over the base ring of the domain(=Free module generated by {'x', 'y'} over Real Field with 53 bits of precision)

                sage: Y = CombinatorialFreeModule(RR['q'],['z'])
                sage: phi = Y.module_morphism( on_basis=on_basis, codomain=X )
                Traceback (most recent call last):
                ...
                ValueError: codomain(=Free module generated by {'x', 'y'} over Real Field with 53 bits of precision) should be a module over the base ring of the domain(=Free module generated by {'z'} over Univariate Polynomial Ring in q over Real Field with 53 bits of precision)


            With the ``diagonal`` argument, this returns the module
            morphism `g` such that:

                `g(x_i) = d(i) y_i`

            This assumes that the respective bases `x` and `y` of `X`
            and `Y` have the same index set `I`.

            With ``triangular = upper``, the constructed module
            morphism is assumed to be upper triangular; that is its
            matrix in the distinguished basis of `X` and `Y` would be
            upper triangular with invertible elements on its
            diagonal. This is used to compute preimages and
            inverting the morphism::

                sage: I = range(1,200)
                sage: X = CombinatorialFreeModule(QQ, I); X.rename("X"); x = X.basis()
                sage: Y = CombinatorialFreeModule(QQ, I); Y.rename("Y"); y = Y.basis()
                sage: f = Y.sum_of_monomials * divisors
                sage: phi = X.module_morphism(f, triangular="upper", codomain = Y)
                sage: phi(x[2])
                B[1] + B[2]
                sage: phi(x[6])
                B[1] + B[2] + B[3] + B[6]
                sage: phi(x[30])
                B[1] + B[2] + B[3] + B[5] + B[6] + B[10] + B[15] + B[30]
                sage: phi.preimage(y[2])
                -B[1] + B[2]
                sage: phi.preimage(y[6])
                B[1] - B[2] - B[3] + B[6]
                sage: phi.preimage(y[30])
                -B[1] + B[2] + B[3] + B[5] - B[6] - B[10] - B[15] + B[30]
                sage: (phi^-1)(y[30])
                -B[1] + B[2] + B[3] + B[5] - B[6] - B[10] - B[15] + B[30]

            For details and further optional arguments, see
            :class:`sage.categories.modules_with_basis.TriangularModuleMorphism`.


            Caveat: the returned element is in ``Hom(codomain, domain,
            category``). This is only correct for unary functions.

            .. TODO::

                Should codomain be ``self`` by default in the
                diagonal and triangular cases?
            """
            if diagonal is not None:
                return DiagonalModuleMorphism(diagonal = diagonal, domain = self, **keywords)
            elif on_basis is not None:
                if triangular is not None:
                    return TriangularModuleMorphism(on_basis, domain = self, triangular = triangular, **keywords)
                return ModuleMorphismByLinearity(on_basis = on_basis, domain = self, **keywords)
            raise ValueError("module morphism requires either on_basis or diagonal argument")

        _module_morphism = module_morphism

        def _repr_(self):
            """
            EXAMPLES::

                sage: class FooBar(CombinatorialFreeModule): pass
                sage: C = FooBar(QQ, (1,2,3)); C # indirect doctest
                Free module generated by {1, 2, 3} over Rational Field

                sage: C._name = "foobar"; C
                foobar over Rational Field

                sage: C.rename("barfoo"); C
                barfoo

                sage: class FooBar(Parent):
                ....:     def basis(self): return Family({1:"foo", 2:"bar"})
                ....:     def base_ring(self): return QQ
                sage: FooBar(category = ModulesWithBasis(QQ))
                Free module generated by [1, 2] over Rational Field
            """
            if hasattr(self, "_name"):
                name = self._name
            else:
                name = "Free module generated by {}".format(self.basis().keys())
            return name + " over {}".format(self.base_ring())


        def tensor(*parents):
            """
            Return the tensor product of the parents.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example(); A.rename("A")
                sage: A.tensor(A,A)
                A # A # A
                sage: A.rename(None)
            """
            return parents[0].__class__.Tensor(parents, category = tensor.category_from_parents(parents))


    class ElementMethods:
        # TODO: Define the appropriate element methods here (instead of in
        # subclasses).  These methods should be consistent with those on
        # polynomials.

#         def _neg_(self):
#             """
#             Default implementation of negation by trying to multiply by -1.
#             TODO: doctest
#             """
#             return self._lmul_(-self.parent().base_ring().one(), self)


        def support_of_term(self):
            """
            Return the support of ``self``, where ``self`` is a monomial
            (possibly with coefficient).

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1,2,3,4]); X.rename("X")
                sage: X.monomial(2).support_of_term()
                2
                sage: X.term(3, 2).support_of_term()
                3

            An exception is raised if ``self`` has more than one term::

                sage: (X.monomial(2) + X.monomial(3)).support_of_term()
                Traceback (most recent call last):
                ...
                ValueError: B[2] + B[3] is not a single term
            """
            if len(self) == 1:
                return self.support()[0]
            else:
                raise ValueError("%s is not a single term"%(self))

        def leading_support(self, cmp=None):
            r"""
            Return the maximal element of the support of ``self``. Note
            that this may not be the term which actually appears first when
            ``self`` is printed.

            If the default ordering of the basis elements is not what is
            desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, ``0`` if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + 4*X.monomial(3)
                sage: x.leading_support()
                3
                sage: def cmp(x,y): return y-x
                sage: x.leading_support(cmp=cmp)
                1

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_support()
                [3]
            """
            return max_cmp(self.support(), cmp)


        def leading_item(self, cmp=None):
            r"""
            Return the pair ``(k, c)`` where

            .. MATH::

                c \cdot (\mbox{the basis element indexed by } k)

            is the leading term of ``self``.

            Here 'leading term' means that the corresponding basis element is
            maximal.  Note that this may not be the term which actually appears
            first when ``self`` is printed.  If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be
            provided.  This should return a negative value if ``x < y``, ``0``
            if ``x == y`` and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + 4*X.monomial(3)
                sage: x.leading_item()
                (3, 4)
                sage: def cmp(x,y): return y-x
                sage: x.leading_item(cmp=cmp)
                (1, 3)

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_item()
                ([3], -5)
            """
            k = self.leading_support(cmp=cmp)
            return k, self[k]

        def leading_monomial(self, cmp=None):
            r"""
            Return the leading monomial of ``self``.

            This is the monomial whose corresponding basis element is
            maximal. Note that this may not be the term which actually appears
            first when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, ``0`` if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.leading_monomial()
                B[3]
                sage: def cmp(x,y): return y-x
                sage: x.leading_monomial(cmp=cmp)
                B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_monomial()
                s[3]
            """
            return self.parent().monomial( self.leading_support(cmp=cmp) )

        def leading_coefficient(self, cmp=None):
            r"""
            Returns the leading coefficient of ``self``.

            This is the coefficient of the term whose corresponding basis element is
            maximal. Note that this may not be the term which actually appears
            first when ``self`` is printed.  If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X")
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.leading_coefficient()
                1
                sage: def cmp(x,y): return y-x
                sage: x.leading_coefficient(cmp=cmp)
                3

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_coefficient()
                -5
            """
            return self.leading_item(cmp=cmp)[1]

        def leading_term(self, cmp=None):
            r"""
            Return the leading term of ``self``.

            This is the term whose corresponding basis element is
            maximal. Note that this may not be the term which actually appears
            first when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.leading_term()
                B[3]
                sage: def cmp(x,y): return y-x
                sage: x.leading_term(cmp=cmp)
                3*B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.leading_term()
                -5*s[3]
            """
            return self.parent().term(*self.leading_item(cmp=cmp))

        def trailing_support(self, cmp=None):
            r"""
            Return the minimal element of the support of ``self``. Note
            that this may not be the term which actually appears last when
            ``self`` is printed.

            If the default ordering of the basis elements is not what is
            desired, a comparison function, ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, `0` if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + 4*X.monomial(3)
                sage: x.trailing_support()
                1
                sage: def cmp(x,y): return y-x
                sage: x.trailing_support(cmp=cmp)
                3

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_support()
                [1]
            """
            return min_cmp(self.support(), cmp)

        def trailing_item(self, cmp=None):
            r"""
            Returns the pair ``(c, k)`` where ``c*self.parent().monomial(k)``
            is the trailing term of ``self``.

            This is the monomial whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed.  If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_item()
                (1, 3)
                sage: def cmp(x,y): return y-x
                sage: x.trailing_item(cmp=cmp)
                (3, 1)

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_item()
                ([1], 2)
            """
            k = self.trailing_support(cmp=cmp)
            return k, self[k]

        def trailing_monomial(self, cmp=None):
            r"""
            Return the trailing monomial of ``self``.

            This is the monomial whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_monomial()
                B[1]
                sage: def cmp(x,y): return y-x
                sage: x.trailing_monomial(cmp=cmp)
                B[3]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_monomial()
                s[1]
            """
            return self.parent().monomial( self.trailing_support(cmp=cmp) )

        def trailing_coefficient(self, cmp=None):
            r"""
            Return the trailing coefficient of ``self``.

            This is the coefficient of the monomial whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_coefficient()
                3
                sage: def cmp(x,y): return y-x
                sage: x.trailing_coefficient(cmp=cmp)
                1

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_coefficient()
                2
            """

            return self.trailing_item(cmp=cmp)[1]

        def trailing_term(self, cmp=None):
            r"""
            Return the trailing term of ``self``.

            This is the term whose corresponding basis element is
            minimal. Note that this may not be the term which actually appears
            last when ``self`` is printed. If the default term ordering is not
            what is desired, a comparison function ``cmp(x,y)``, can be provided.
            This should return a negative value if ``x < y``, 0 if ``x == y``
            and a positive value if ``x > y``.

            EXAMPLES::

                sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
                sage: x = 3*X.monomial(1) + 2*X.monomial(2) + X.monomial(3)
                sage: x.trailing_term()
                3*B[1]
                sage: def cmp(x,y): return y-x
                sage: x.trailing_term(cmp=cmp)
                B[3]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = 2*s[1] + 3*s[2,1] - 5*s[3]
                sage: f.trailing_term()
                2*s[1]
            """
            return self.parent().term( *self.trailing_item(cmp=cmp) )

        def map_coefficients(self, f):
            """
            Mapping a function on coefficients.

            INPUT:

            - ``f`` -- an endofunction on the coefficient ring of the
              free module

            Return a new element of ``self.parent()`` obtained by applying the
            function ``f`` to all of the coefficients of ``self``.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
                sage: B = F.basis()
                sage: f = B['a'] - 3*B['c']
                sage: f.map_coefficients(lambda x: x+5)
                6*B['a'] + 2*B['c']

            Killed coefficients are handled properly::

                sage: f.map_coefficients(lambda x: 0)
                0
                sage: list(f.map_coefficients(lambda x: 0))
                []

            ::

                sage: s = SymmetricFunctions(QQ).schur()
                sage: a = s([2,1])+2*s([3,2])
                sage: a.map_coefficients(lambda x: x*2)
                2*s[2, 1] + 4*s[3, 2]
            """
            return self.parent().sum_of_terms( (m, f(c)) for m,c in self )

        def map_support(self, f):
            """
            Mapping a function on the support.

            INPUT:

            - ``f`` -- an endofunction on the indices of the free module

            Return a new element of ``self.parent()`` obtained by
            applying the function ``f`` to all of the objects indexing
            the basis elements.

            EXAMPLES::

                sage: B = CombinatorialFreeModule(ZZ, [-1, 0, 1])
                sage: x = B.an_element(); x
                2*B[-1] + 2*B[0] + 3*B[1]
                sage: x.map_support(lambda i: -i)
                3*B[-1] + 2*B[0] + 2*B[1]

            ``f`` needs not be injective::

                sage: x.map_support(lambda i: 1)
                7*B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: a = s([2,1])+2*s([3,2])
                sage: a.map_support(lambda x: x.conjugate())
                s[2, 1] + 2*s[2, 2, 1]

            TESTS::

                sage: B.zero()      # This actually failed at some point!!! See #8890
                0

                sage: y = B.zero().map_support(lambda i: i/0); y
                0
                sage: y.parent() is B
                True
            """
            return self.parent().sum_of_terms( (f(m), c) for m,c in self )

        def map_support_skip_none(self, f):
            """
            Mapping a function on the support.

            INPUT:

            - ``f`` -- an endofunction on the indices of the free module

            Returns a new element of ``self.parent()`` obtained by
            applying the function `f` to all of the objects indexing
            the basis elements.

            EXAMPLES::

                sage: B = CombinatorialFreeModule(ZZ, [-1, 0, 1])
                sage: x = B.an_element(); x
                2*B[-1] + 2*B[0] + 3*B[1]
                sage: x.map_support_skip_none(lambda i: -i if i else None)
                3*B[-1] + 2*B[1]

            ``f`` needs not be injective::

                sage: x.map_support_skip_none(lambda i: 1 if i else None)
                5*B[1]

            TESTS::

                sage: y = x.map_support_skip_none(lambda i: None); y
                0
                sage: y.parent() is B
                True
            """
            return self.parent().sum_of_terms( (fm,c) for (fm,c) in ((f(m), c) for m,c in self) if fm is not None)

        def map_item(self, f):
            """
            Mapping a function on items.

            INPUT:

            - ``f`` -- a function mapping pairs ``(index, coeff)`` to
              other such pairs

            Return a new element of ``self.parent()`` obtained by
            applying the function `f` to all items ``(index, coeff)`` of
            ``self``.

            EXAMPLES::

                sage: B = CombinatorialFreeModule(ZZ, [-1, 0, 1])
                sage: x = B.an_element(); x
                2*B[-1] + 2*B[0] + 3*B[1]
                sage: x.map_item(lambda i, c: (-i, 2*c))
                6*B[-1] + 4*B[0] + 4*B[1]

            ``f`` needs not be injective::

                sage: x.map_item(lambda i, c: (1, 2*c))
                14*B[1]

                sage: s = SymmetricFunctions(QQ).schur()
                sage: f = lambda m,c: (m.conjugate(), 2*c)
                sage: a = s([2,1]) + s([1,1,1])
                sage: a.map_item(f)
                2*s[2, 1] + 2*s[3]

            The following methods are deprecated::

                sage: a.map_term(f)
                doctest:...: DeprecationWarning: map_term is deprecated. Please use map_item instead.
                See http://trac.sagemath.org/8890 for details.
                2*s[2, 1] + 2*s[3]
                sage: a.map_mc(f)
                doctest:...: DeprecationWarning: map_mc is deprecated. Please use map_item instead.
                See http://trac.sagemath.org/8890 for details.
                2*s[2, 1] + 2*s[3]
            """
            return self.parent().sum_of_terms( f(m,c) for m,c in self )

        map_term = deprecated_function_alias(8890, map_item)
        map_mc   = deprecated_function_alias(8890, map_item)

        def tensor(*elements):
            """
            Return the tensor product of its arguments, as an element of
            the tensor product of the parents of those elements.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example()
                sage: (a,b,c) = A.algebra_generators()
                sage: a.tensor(b, c)
                B[word: a] # B[word: b] # B[word: c]

            FIXME: is this a policy that we want to enforce on all parents?
            """
            assert(all(isinstance(element, Element) for element in elements))
            parents = [parent(element) for element in elements]
            return tensor(parents)._tensor_of_elements(elements) # good name???

    class HomCategory(HomCategory):
        """
        The category of homomorphisms sets `Hom(X,Y)` for `X`, `Y`
        modules with basis.
        """
        class ParentMethods:
            def __call_on_basis__(self, **options):
                """
                Construct a morphism in this homset from a function defined
                on the basis.

                INPUT:

                - ``on_basis`` -- a function from the indices of the
                  basis of the domain of ``self`` to the codomain of
                  ``self``

                This method simply delegates the work to
                :meth:`ModulesWithBasis.ParentMethods.module_morphism`. It
                is used by :meth:`Homset.__call__` to handle the
                ``on_basis`` argument, and will disapear as soon as
                the logic will be generalized.

                EXAMPLES::

                    sage: X = CombinatorialFreeModule(QQ, [1,2,3]);   X.rename("X")
                    sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4]); Y.rename("Y")
                    sage: H = Hom(X, Y)
                    sage: x = X.basis()

                    sage: phi = H(on_basis = lambda i: Y.monomial(i) + 2*Y.monomial(i+1)) # indirect doctest
                    sage: phi
                    Generic morphism:
                    From: X
                    To:   Y
                    sage: phi(x[1] + x[3])
                    B[1] + 2*B[2] + B[3] + 2*B[4]

                Diagonal functions can be constructed using the ``diagonal`` option::

                    sage: X = CombinatorialFreeModule(QQ, [1,2,3,4]); X.rename("X")
                    sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4], key="Y"); Y.rename("Y")
                    sage: H = Hom(X, Y)
                    sage: x = X.basis()
                    sage: phi = H(diagonal = lambda x: x^2)
                    sage: phi(x[1] + x[2] + x[3])
                    B[1] + 4*B[2] + 9*B[3]

                TESTS::

                As for usual homsets, the argument can be a Python function::

                    sage: phi = H(lambda x: Y.zero())
                    sage: phi
                    Generic morphism:
                      From: X
                      To:   Y
                    sage: phi(x[1] + x[3])
                    0

               We check that the homset category is properly set up::

                    sage: category = FiniteDimensionalModulesWithBasis(QQ)
                    sage: X = CombinatorialFreeModule(QQ, [1,2,3], category = category);   X.rename("X")
                    sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4], category = category); Y.rename("Y")
                    sage: H = Hom(X, Y)
                    sage: H.zero().category_for()
                    Category of finite dimensional vector spaces with basis over Rational Field
                """
                return self.domain().module_morphism(codomain = self.codomain(),
                                                     **options)

        class ElementMethods:
            """
            Abstract class for morphisms of modules with basis.
            """
            @cached_method
            def on_basis(self):
                """
                Return the action of this morphism on basis elements.

                OUTPUT:

                - a function from the indices of the basis of the domain to
                  the codomain

                EXAMPLES::

                    sage: X = CombinatorialFreeModule(QQ, [1,2,3]);   X.rename("X")
                    sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4]); Y.rename("Y")
                    sage: H = Hom(X, Y)
                    sage: x = X.basis()

                    sage: f = H(lambda x: Y.zero()).on_basis()
                    sage: f(2)
                    0

                    sage: f = lambda i: Y.monomial(i) + 2*Y.monomial(i+1)
                    sage: g = H(on_basis = f).on_basis()
                    sage: g(2)
                    B[2] + 2*B[3]
                    sage: g == f
                    True
                """
                monomial = self.domain().monomial
                return lambda t: self(monomial(t))

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of modules with basis constructed by cartesian products
        of modules with basis.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: ModulesWithBasis(QQ).CartesianProducts().extra_super_categories()
                [Category of modules with basis over Rational Field]
                sage: ModulesWithBasis(QQ).CartesianProducts().super_categories()
                [Category of modules with basis over Rational Field,
                 Category of Cartesian products of commutative additive groups]
            """
            return [self.base_category()]

        class ParentMethods:

            def _an_element_(self):
                """
                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: B = HopfAlgebrasWithBasis(QQ).example(); B
                    An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                    sage: A.an_element(), B.an_element()
                    (2*B[word: ] + 2*B[word: a] + 3*B[word: b], B[()] + 2*B[(2,3)] + 3*B[(1,2)] + B[(1,2,3)])
                    sage: cartesian_product((A, B, A)).an_element()           # indirect doctest
                    2*B[(0, word: )] + 2*B[(0, word: a)] + 3*B[(0, word: b)]
                """
                from cartesian_product import cartesian_product
                return cartesian_product([module.an_element() for module in self.modules])

    class TensorProducts(TensorProductsCategory):
        """
        The category of modules with basis constructed by tensor product of
        modules with basis.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: ModulesWithBasis(QQ).TensorProducts().extra_super_categories()
                [Category of modules with basis over Rational Field]
                sage: ModulesWithBasis(QQ).TensorProducts().super_categories()
                [Category of modules with basis over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            Implements operations on tensor products of modules with basis.
            """
            pass

        class ElementMethods:
            """
            Implements operations on elements of tensor products of modules
            with basis.
            """

            def apply_multilinear_morphism(self, f, codomain = None):
                r"""
                Return the result of applying the morphism induced by ``f``
                to ``self``.

                INPUT:

                - ``f`` -- a multilinear morphism from the component
                  modules of the parent tensor product to any module

                - ``codomain`` -- the codomain of ``f`` (optional)

                By the universal property of the tensor product, ``f``
                induces a linear morphism from `self.parent()` to the
                target module. Returns the result of applying that
                morphism to ``self``.

                The codomain is used for optimizations purposes
                only. If it's not provided, it's recovered by calling
                ``f`` on the zero input.

                EXAMPLES:

                We start with simple (admittedly not so interesting)
                examples, with two modules `A` and `B`::

                    sage: A = CombinatorialFreeModule(ZZ, [1,2], prefix="A"); A.rename("A")
                    sage: B = CombinatorialFreeModule(ZZ, [3,4], prefix="B"); B.rename("B")

                and `f` the bilinear morphism `(a,b) \mapsto b \otimes a`
                from `A \times B` to `B \otimes A`::

                    sage: def f(a,b):
                    ....:     return tensor([b,a])

                Now, calling applying `f` on `a \otimes b` returns the same
                as `f(a,b)`::

                    sage: a = A.monomial(1) + 2 * A.monomial(2); a
                    A[1] + 2*A[2]
                    sage: b = B.monomial(3) - 2 * B.monomial(4); b
                    B[3] - 2*B[4]
                    sage: f(a,b)
                    B[3] # A[1] + 2*B[3] # A[2] - 2*B[4] # A[1] - 4*B[4] # A[2]
                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    B[3] # A[1] + 2*B[3] # A[2] - 2*B[4] # A[1] - 4*B[4] # A[2]

                `f` may be a bilinear morphism to any module over the
                base ring of `A` and `B`. Here the codomain is `\ZZ`::

                    sage: def f(a,b):
                    ....:     return sum(a.coefficients(), 0) * sum(b.coefficients(), 0)
                    sage: f(a,b)
                    -3
                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    -3

                Mind the `0` in the sums above; otherwise `f` would
                not return `0` in `\ZZ`::

                    sage: def f(a,b):
                    ....:     return sum(a.coefficients()) * sum(b.coefficients())
                    sage: type(f(A.zero(), B.zero()))
                    <type 'int'>

                Which would be wrong and break this method::

                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    Traceback (most recent call last):
                    ...
                    AttributeError: 'int' object has no attribute 'parent'

                Here we consider an example where the codomain is a
                module with basis with a different base ring::

                    sage: C = CombinatorialFreeModule(QQ, [(1,3),(2,4)], prefix="C"); C.rename("C")
                    sage: def f(a,b):
                    ....:     return C.sum_of_terms( [((1,3), QQ(a[1]*b[3])), ((2,4), QQ(a[2]*b[4]))] )
                    sage: f(a,b)
                    C[(1, 3)] - 4*C[(2, 4)]
                    sage: tensor([a,b]).apply_multilinear_morphism(f)
                    C[(1, 3)] - 4*C[(2, 4)]

                 We conclude with a real life application, where we
                 check that the antipode of the Hopf algebra of
                 Symmetric functions on the Schur basis satisfies its
                 defining formula::

                    sage: Sym = SymmetricFunctions(QQ)
                    sage: s = Sym.schur()
                    sage: def f(a,b): return a*b.antipode()
                    sage: x = 4*s.an_element(); x
                    8*s[] + 8*s[1] + 12*s[2]
                    sage: x.coproduct().apply_multilinear_morphism(f)
                    8*s[]
                    sage: x.coproduct().apply_multilinear_morphism(f) == x.counit()
                    True

                We recover the constant term of `x`, as desired.

                .. TODO::

                    Extract a method to linearize a multilinear
                    morphism, and delegate the work there.
                """
                K = self.parent().base_ring()
                modules = self.parent()._sets
                if codomain is None:
                    try:
                        codomain = f.codomain()
                    except AttributeError:
                        codomain = f(*[module.zero() for module in modules]).parent()
                if codomain in ModulesWithBasis(K):
                    return codomain.linear_combination((f(*[module.monomial(t) for (module,t) in zip(modules, m)]), c)
                                                       for m,c in self)
                else:
                    return sum((c * f(*[module.monomial(t) for (module,t) in zip(modules, m)])
                                for m,c in self),
                               codomain.zero())

    class DualObjects(DualObjectsCategory):

        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: ModulesWithBasis(ZZ).DualObjects().extra_super_categories()
                [Category of modules over Integer Ring]
                sage: ModulesWithBasis(QQ).DualObjects().super_categories()
                [Category of duals of vector spaces over Rational Field, Category of duals of modules with basis over Rational Field]
            """
            return [Modules(self.base_category().base_ring())]


class ModuleMorphismByLinearity(Morphism):
    """
    A class for module morphisms obtained by extending a function by linearity.
    """
    def __init__(self, domain, on_basis = None, position = 0, zero = None, codomain = None, category = None):
        """
        Construct a module morphism by linearity.

        INPUT:

        - ``domain`` -- a parent in ``ModulesWithBasis(...)``
        - ``codomain`` -- a parent in ``Modules(...)``; defaults to
          ``f.codomain()`` if the latter is defined
        - ``position`` -- a non-negative integer; defaults to 0
        - ``on_basis`` -- a function which accepts indices of the basis of
          ``domain`` as ``position``-th argument (optional)
        - ``zero`` -- the zero of the codomain; defaults to ``codomain.zero()``

        ``on_basis`` may alternatively be provided in derived classes by
        implementing or setting ``_on_basis``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(ZZ, [-2, -1, 1, 2])
            sage: Y = CombinatorialFreeModule(ZZ, [1, 2])
            sage: phi = sage.categories.modules_with_basis.ModuleMorphismByLinearity(X, on_basis = Y.monomial * abs)

        TESTS::

            sage: TestSuite(phi).run() # known issue
            Failure in _test_category:
            ...
            The following tests failed: _test_category

        Analysis: ``phi`` does not inherit from the element class of
        the category of its parent::

            sage: isinstance(phi, phi.parent().category().element_class)
            False

        To be fixed in the general morphism overhaul (#....), possibly
        by making sure to create ``phi`` through its parent.
        """
        # Might want to assert that domain is a module with basis
        base_ring = domain.base_ring()

        if codomain is None and hasattr(on_basis, 'codomain'):
            codomain = on_basis.codomain()
        if not hasattr( codomain, 'base_ring' ):
            raise ValueError("codomain(=%s) needs to have a base_ring attribute"%(codomain))
        # codomain should be a module over base_ring
        # The natural test would be ``codomains in Modules(base_ring)``
        # But this is not properly implemented yet:
        #     sage: CC in Modules(QQ)
        #     False
        #     sage: QQ in Modules(QQ)
        #     False
        #     sage: CC[x] in Modules(QQ)
        #     False
        # The test below is a bit more restrictive
        if (not codomain.base_ring().has_coerce_map_from(base_ring)) \
           and (not codomain.has_coerce_map_from(base_ring)):
            raise ValueError("codomain(=%s) should be a module over the base ring of the domain(=%s)"%(codomain, domain))

        if zero is None:
            zero = codomain.zero()
        self._zero = zero

        self._is_module_with_basis_over_same_base_ring = \
            codomain in ModulesWithBasis( base_ring ) and zero == codomain.zero()

        if category is None:
            if self._is_module_with_basis_over_same_base_ring:
                category = ModulesWithBasis(base_ring)
            elif zero == codomain.zero():
                if codomain in Modules(base_ring):
                    category = Modules(base_ring)
                else:
                    # QQ is not in Modules(QQ)!
                    category = CommutativeAdditiveSemigroups()
            else:
                category = Sets()

        Morphism.__init__(self, Hom(domain, codomain, category = category))


        self._position = position
        if on_basis is not None:
            self._on_basis = on_basis

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: X = CombinatorialFreeModule(ZZ, [-2, -1, 1, 2])
            sage: Y = CombinatorialFreeModule(ZZ, [1, 2])
            sage: f  = X.module_morphism(on_basis = Y.monomial * abs)
            sage: g  = X.module_morphism(on_basis = Y.monomial * abs)
            sage: h1 = X.module_morphism(on_basis = X.monomial * abs)
            sage: h2 = X.module_morphism(on_basis = X.monomial * factorial)
            sage: h3 = X.module_morphism(on_basis = Y.monomial * abs, category = Modules(ZZ))
            sage: f == g, f == h1, f == h2, f == h3, f == 1, 1 == f
            (True, False, False, False, False, False)
        """
        return self.__class__ is other.__class__ and parent(self) == parent(other) and self.__dict__ == other.__dict__


    def on_basis(self):
        """
        Return the action of this morphism on basis elements, as per
        :meth:`ModulesWithBasis.HomCategory.ElementMethods.on_basis`.

        OUTPUT:

        - a function from the indices of the basis of the domain to the
          codomain

        EXAMPLES::

            sage: X = CombinatorialFreeModule(ZZ, [-2, -1, 1, 2])
            sage: Y = CombinatorialFreeModule(ZZ, [1, 2])
            sage: phi_on_basis = Y.monomial * abs
            sage: phi = sage.categories.modules_with_basis.ModuleMorphismByLinearity(X, on_basis = phi_on_basis, codomain = Y)
            sage: x = X.basis()
            sage: phi.on_basis()(-2)
            B[2]
            sage: phi.on_basis() == phi_on_basis
            True
        """
        return self._on_basis

    def __call__(self, *args):
        r"""
        Apply this morphism to ``*args``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(ZZ, [-2, -1, 1, 2])
            sage: Y = CombinatorialFreeModule(ZZ, [1, 2])
            sage: def phi_on_basis(i): return Y.monomial(abs(i))
            sage: phi = sage.categories.modules_with_basis.ModuleMorphismByLinearity(X, on_basis = Y.monomial * abs, codomain = Y)
            sage: x = X.basis()
            sage: phi(x[1]), phi(x[-2]), phi(x[1] + 3 * x[-2])
            (B[1], B[2], B[1] + 3*B[2])

        .. TODO::

            Add more tests for multi-parameter module morphisms.
        """
        before = args[0:self._position]
        after = args[self._position+1:len(args)]
        x = args[self._position]
        assert(x.parent() is self.domain())

        if self._is_module_with_basis_over_same_base_ring:
            return self.codomain().linear_combination( (self._on_basis(*(before+(index,)+after)), coeff ) for (index, coeff) in args[self._position] )
        else:
            return sum(( coeff * self._on_basis(*(before+(index,)+after)) for (index, coeff) in args[self._position]), self._zero)

    # As per the specs of Map, we should in fact implement _call_.
    # However we currently need to abuse Map.__call__ (which strict
    # type checking) for multi-parameter module morphisms
    # To be cleaned up
    _call_ = __call__

class TriangularModuleMorphism(ModuleMorphismByLinearity):
    r"""
    A class for triangular module morphisms; that is, module morphisms
    from `X` to `Y` whose representing matrix in the distinguished
    bases of `X` and `Y` is upper triangular with invertible elements
    on its diagonal.

    See :meth:`ModulesWithBasis.ParentMethods.module_morphism`

    INPUT:

    - ``domain`` -- a module `X` with basis `F`
    - ``codomain`` -- a module `Y` with basis `G` (defaults to `X`)
    - ``on_basis`` -- a function from the index set of the basis `F`
      to the module `Y` which determines the morphism by linearity
    - ``unitriangular`` -- boolean (default: ``False``)
    - ``triangular`` -- (default: ``"upper"``) ``"upper"`` or ``"lower"``:

      * ``"upper"`` - if the :meth:`leading_support()` of the image of
        `F(i)` is `i`, or
      * ``"lower"`` - if the :meth:`trailing_support()` of the image of
        `F(i)` is `i`

    - ``cmp`` -- an optional comparison function on the index set `J` of
      the basis `G` of the codomain.
    - ``invertible`` -- boolean or ``None`` (default: ``None``); should
      be set to ``True`` if Sage is to compute an inverse for ``self``.
      Automatically set to ``True`` if the domain and codomain share the
      same indexing set and to ``False`` otherwise.
    - ``inverse_on_support`` - compute the inverse on the support if the
      codomain and domain have different index sets. See assumptions
      below.

    Assumptions:

    - `X` and `Y` have the same base ring `R`.

    - Let `I` and `J` be the respective index sets of the bases `F` and
      `G`. Either `I = J`, or ``inverse_on_support`` is a
      function `r : J \to I` with the following property: for any `j \in J`,
      `r(j)` should return an `i \in I` such that the leading term (or
      trailing term, if ``triangular`` is set to ``"lower"``) of
      ``on_basis(i)`` (with respect to the comparison ``cmp``, if the
      latter is set, or just the default comparison otherwise) is `j` if
      there exists such an `i`, or ``None`` if not.

    OUTPUT:

    The triangular module morphism from `X` to `Y` which maps `F(i)`
    to ``on_basis(i)`` and is extended by linearity.

    EXAMPLES:

    We construct and invert an upper unitriangular module morphism between
    two free `\QQ`-modules::

        sage: I = range(1,200)
        sage: X = CombinatorialFreeModule(QQ, I); X.rename("X"); x = X.basis()
        sage: Y = CombinatorialFreeModule(QQ, I); Y.rename("Y"); y = Y.basis()
        sage: f = Y.sum_of_monomials * divisors   # This * is map composition.
        sage: phi = X.module_morphism(f, triangular="upper", unitriangular = True, codomain = Y)
        sage: phi(x[2])
        B[1] + B[2]
        sage: phi(x[6])
        B[1] + B[2] + B[3] + B[6]
        sage: phi(x[30])
        B[1] + B[2] + B[3] + B[5] + B[6] + B[10] + B[15] + B[30]
        sage: phi.preimage(y[2])
        -B[1] + B[2]
        sage: phi.preimage(y[6])
        B[1] - B[2] - B[3] + B[6]
        sage: phi.preimage(y[30])
        -B[1] + B[2] + B[3] + B[5] - B[6] - B[10] - B[15] + B[30]
        sage: (phi^-1)(y[30])
        -B[1] + B[2] + B[3] + B[5] - B[6] - B[10] - B[15] + B[30]

    A lower triangular (but not unitriangular) morphism::

        sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
        sage: def ut(i): return sum(j*x[j] for j in range(i,4))
        sage: phi = X.module_morphism(ut, triangular="lower", codomain = X)
        sage: phi(x[2])
        2*B[2] + 3*B[3]
        sage: phi.preimage(x[2])
        1/2*B[2] - 1/2*B[3]
        sage: phi(phi.preimage(x[2]))
        B[2]

    Using the ``cmp`` keyword, we can use triangularity even if
    the map becomes triangular only after a permutation of the basis::

        sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
        sage: def vt(i): return (x[1] + x[2] if i == 1 else x[2] + (x[3] if i == 3 else 0))
        sage: perm = [0, 2, 1, 3]
        sage: phi = X.module_morphism(vt, triangular="upper", codomain = X,
        ....:                         cmp=lambda a, b: cmp(perm[a], perm[b]))
        sage: [phi(x[i]) for i in range(1, 4)]
        [B[1] + B[2], B[2], B[2] + B[3]]
        sage: [phi.preimage(x[i]) for i in range(1, 4)]
        [B[1] - B[2], B[2], -B[2] + B[3]]

    The same works in the lower-triangular case::

        sage: def wt(i): return (x[1] + x[2] + x[3] if i == 2 else x[i])
        sage: phi = X.module_morphism(wt, triangular="lower", codomain = X,
        ....:                         cmp=lambda a, b: cmp(perm[a], perm[b]))
        sage: [phi(x[i]) for i in range(1, 4)]
        [B[1], B[1] + B[2] + B[3], B[3]]
        sage: [phi.preimage(x[i]) for i in range(1, 4)]
        [B[1], -B[1] + B[2] - B[3], B[3]]

    An injective but not surjective morphism cannot be inverted,
    but the ``inverse_on_support`` keyword allows Sage to find a
    partial inverse::

        sage: X = CombinatorialFreeModule(QQ, [1,2,3]); x = X.basis()
        sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4,5]); y = Y.basis()
        sage: uut = lambda i: sum(  y[j] for j in range(i+1,6)  )
        sage: phi = X.module_morphism(uut, codomain = Y,
        ....:        triangular=True, unitriangular=True,
        ....:        inverse_on_support=lambda i: i-1 if i in [2,3,4] else None)
        sage: phi(x[2])
        B[3] + B[4] + B[5]
        sage: phi.preimage(y[3])
        B[2] - B[3]

    The ``inverse_on_support`` keyword can also be used if the
    bases of the domain and the codomain are identical but one of
    them has to be permuted in order to render the morphism
    triangular. For example::

        sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
        sage: def zt(i):
        ....:     return (x[3] if i == 1 else x[1] if i == 2
        ....:             else x[1] + x[2])
        sage: def perm(i):
        ....:     return (2 if i == 1 else 3 if i == 2 else 1)
        sage: phi = X.module_morphism(zt, triangular="upper", codomain = X,
        ....:                         inverse_on_support=perm)
        sage: [phi(x[i]) for i in range(1, 4)]
        [B[3], B[1], B[1] + B[2]]
        sage: [phi.preimage(x[i]) for i in range(1, 4)]
        [B[2], -B[2] + B[3], B[1]]

    The same works if the permutation induces lower triangularity::

        sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
        sage: def zt(i):
        ....:     return (x[3] if i == 1 else x[2] if i == 2
        ....:             else x[1] + x[2])
        sage: def perm(i):
        ....:     return 4 - i
        sage: phi = X.module_morphism(zt, triangular="lower", codomain = X,
        ....:                         inverse_on_support=perm)
        sage: [phi(x[i]) for i in range(1, 4)]
        [B[3], B[2], B[1] + B[2]]
        sage: [phi.preimage(x[i]) for i in range(1, 4)]
        [-B[2] + B[3], B[2], B[1]]

    The ``inverse_on_basis`` and ``cmp`` keywords can be combined::

        sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
        sage: def zt(i):
        ....:     return (2*x[2] + 3*x[3] if i == 1
        ....:             else x[1] + x[2] + x[3] if i == 2
        ....:             else 4*x[2])
        sage: def perm(i):
        ....:     return (2 if i == 1 else 3 if i == 2 else 1)
        sage: perverse_cmp = lambda a, b: cmp((a-2) % 3, (b-2) % 3)
        sage: phi = X.module_morphism(zt, triangular="upper", codomain = X,
        ....:                         inverse_on_support=perm, cmp=perverse_cmp)
        sage: [phi(x[i]) for i in range(1, 4)]
        [2*B[2] + 3*B[3], B[1] + B[2] + B[3], 4*B[2]]
        sage: [phi.preimage(x[i]) for i in range(1, 4)]
        [-1/3*B[1] + B[2] - 1/12*B[3], 1/4*B[3], 1/3*B[1] - 1/6*B[3]]
    """
    def __init__(self, on_basis, domain, triangular = "upper", unitriangular=False,
                 codomain = None, category = None, cmp = None,
                 inverse = None, inverse_on_support = None, invertible = None):
        """
        TESTS::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
            sage: def ut(i): return sum(j*x[j] for j in range(i,4))
            sage: import __main__; __main__.ut = ut
            sage: phi = X.module_morphism(ut, triangular="lower", codomain = X)
            sage: phi.__class__
            <class 'sage.categories.modules_with_basis.TriangularModuleMorphism'>
            sage: TestSuite(phi).run() # known issue; see ModuleMorphism above
            Failure in _test_category:
            ...
            The following tests failed: _test_category
        """
        ModuleMorphismByLinearity.__init__(self, domain = domain, codomain = codomain, category = category)
        if triangular == "upper":
            self._dominant_item = attrcall("leading_item",  cmp)
        else:
            self._dominant_item = attrcall("trailing_item", cmp)
        # We store those two just be able to pass them down to the inverse function
        self._triangular = triangular
        self._cmp = cmp

        self._unitriangular = unitriangular
        self._inverse = inverse
        self.on_basis = on_basis # should this be called on_basis (or _on_basis)?
        self._inverse_on_support = inverse_on_support
        if invertible is not None:
            self._invertible = invertible
        else:
            self._invertible = (domain.basis().keys() == codomain.basis().keys())

    def _test_triangular(self, **options):
        """
        Test that ``self`` is actually triangular

        See also: :class:`sage.misc.sage_unittest.TestSuite`.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); y = Y.basis()
            sage: f = lambda i: sum(  y[j] for j in range(i,4)  )
            sage: phi = X.module_morphism(f, triangular="lower", codomain = Y)
            sage: phi._test_triangular()

            sage: fw = lambda i: sum(  y[j] for j in range(i+1,4)  )
            sage: phi = X.module_morphism(fw, triangular="lower", codomain = Y)
            sage: phi._test_triangular()
            Traceback (most recent call last):
            ...
            AssertionError: morphism is not triangular on 1

            sage: X = CombinatorialFreeModule(QQ, [1,2,3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4,5]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i+1,6)  )
            sage: phi = X.module_morphism(uut, codomain = Y,
            ....:      triangular=True, unitriangular=True,
            ....:      inverse_on_support=lambda i: i-1 if i in [2,3,4] else None)
            sage: phi._test_triangular()

            sage: uutw = lambda i: sum(  2*y[j] for j in range(i+1,6)  ) # uni-upper
            sage: phi = X.module_morphism(uutw, codomain = Y,
            ....:      triangular=True, unitriangular=True,
            ....:      inverse_on_support=lambda i: i-1 if i in [2,3,4] else None)
            sage: phi._test_triangular()
            Traceback (most recent call last):
            ...
            AssertionError: morphism is not untriangular on 1
        """
        from sage.misc.lazy_format import LazyFormat
        tester = self._tester(**options)
        for x in self.domain().basis().keys():
            # or should it be self.domain().basis().some_elements() # ?
            bs, co = self._dominant_item(self._on_basis(x))
            if self._unitriangular:
                tester.assertEqual(co, self.domain().base_ring().one(),
                    LazyFormat("morphism is not untriangular on %s")%(x))
            if self._inverse_on_support is not None:
                xback = self._inverse_on_support(bs)
            else:
                xback = bs
            tester.assertEqual(x, xback,
                LazyFormat("morphism is not triangular on %s")%(x))


    def _on_basis(self, i):
        """
        Return the image, by ``self``, of the basis element indexed by ``i``.

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); y = Y.basis()
            sage: f = lambda i: sum(  y[j] for j in range(i,4)  )
            sage: phi = X.module_morphism(f, triangular="lower", codomain = Y)
            sage: phi._on_basis(2)
            B[2] + B[3]
        """
        return self.on_basis(i)

    def __invert__(self):
        """
        Return the triangular morphism which is the inverse of ``self``.

        Raises an error if ``self`` is not invertible.

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(1,i+1)) # uni-upper
            sage: ult = lambda i: sum(  y[j] for j in range(i,4)  ) # uni-lower
            sage: ut =  lambda i: sum(j*y[j] for j in range(1,i+1)) # upper
            sage: lt =  lambda i: sum(j*y[j] for j in range(i,4  )) # lower
            sage: f_uut = X.module_morphism(uut, triangular="upper", unitriangular=True,  codomain = Y)
            sage: f_ult = X.module_morphism(ult, triangular="lower", unitriangular=True,  codomain = Y)
            sage: f_ut =  X.module_morphism(ut,  triangular="upper",                      codomain = Y)
            sage: f_lt =  X.module_morphism(lt,  triangular="lower",                      codomain = Y)
            sage: (~f_uut)(y[2])
            -B[1] + B[2]
            sage: (~f_ult)(y[2])
            B[2] - B[3]
            sage: (~f_ut)(y[2])
            -1/2*B[1] + 1/2*B[2]
            sage: (~f_lt)(y[2])
            1/2*B[2] - 1/2*B[3]
        """
        if not self._invertible:
            raise ValueError("Non invertible morphism")
        else:
            return self.section()

    def section(self):
        """
        Return the section (partial inverse) of ``self``.

        Return a partial triangular morphism which is a section of
        ``self``. The section morphism raise a ``ValueError`` if asked to
        apply on an element which is not in the image of ``self``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, [1,2,3]); x = X.basis()
            sage: X.rename('X')
            sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4,5]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i+1,6)  ) # uni-upper
            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y,
            ....:      inverse_on_support=lambda i: i-1 if i in [2,3,4] else None)
            sage: ~phi
            Traceback (most recent call last):
            ...
            ValueError: Non invertible morphism
            sage: phiinv = phi.section()
            sage: map(phiinv*phi, X.basis().list()) == X.basis().list()
            True
            sage: phiinv(Y.basis()[1])
            Traceback (most recent call last):
            ...
            ValueError: B[1] is not in the image
        """
        if self._inverse is not None:
            return self._inverse
        if self._inverse_on_support is None:
            retract_dom = None
        else:
            def retract_dom(i):
                self._dominant_item(self._on_basis(i))[0]

        if self._invertible:
            return self.__class__( self._invert_on_basis,
                domain = self.codomain(),               codomain = self.domain(),
                unitriangular = self._unitriangular,  triangular = self._triangular,
                cmp = self._cmp,
                inverse = self,                       category = self.category_for(),
                inverse_on_support=retract_dom, invertible = self._invertible)
        else:
            return SetMorphism(Hom(self.codomain(), self.domain(),
                                   SetsWithPartialMaps()),
                               self.preimage)

    # This should be removed and optimized (the inverse should not be computed
    # on the basis
    def _invert_on_basis(self, i):
        r"""
        Return the image, by the inverse of ``self``, of the basis element
        indexed by ``i``.

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i,4)  ) # uni-upper
            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y)
            sage: phi._invert_on_basis(2)
            B[2] - B[3]
        """
        return self.preimage( self.codomain().monomial(i) )

    def preimage(self, f):
        """
        Return the preimage of `f` under ``self``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i,4)  ) # uni-upper
            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y)
            sage: phi.preimage(y[1] + y[2])
            B[1] - B[3]

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3, 4]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i,5)  ) # uni-upper
            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y)
            sage: phi.preimage(y[1] + y[2])
            B[1] - B[3]

            sage: X = CombinatorialFreeModule(QQ, [1,2,3]); x = X.basis()
            sage: X.rename("X")
            sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4,5]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i+1,6)  ) # uni-upper
            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y,
            ....:         inverse_on_support=lambda i: i-1 if i in [2,3,4] else None)
            sage: phi.preimage(y[2] + y[3])
            B[1] - B[3]
            sage: phi(phi.preimage(y[2] + y[3])) == y[2] + y[3]
            True
            sage: el = x[1] + 3*x[2] + 2*x[3]
            sage: phi.preimage(phi(el)) == el
            True

            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y,
            ....:         inverse_on_support=lambda i: i-1 if i in [2,3,4] else None)
            sage: phi.preimage(y[1])
            Traceback (most recent call last):
            ...
            ValueError: B[1] is not in the image
        """
        F = self.domain()
        G = self.codomain()
        basis_map = self._on_basis
        if not f in G:
            raise ValueError("f(={}) must be in the codomain of the morphism to have a preimage under the latter".format(f))

        remainder = f

        out = F.zero()
        while not remainder.is_zero():
            (j,c) = self._dominant_item(remainder)

            if self._inverse_on_support is None:
                j_preimage = j
            else:
                j_preimage = self._inverse_on_support(j)
                if j_preimage is None:
                    raise ValueError("{} is not in the image".format(f))
            s = basis_map(j_preimage)
            if not j == self._dominant_item(s)[0]:
                raise ValueError("The morphism (={}) is not triangular at {}, and therefore a preimage cannot be computed".format(f, s))

            if not self._unitriangular:
                c /= s[j]

            remainder -= s._lmul_(c)
            out += F.term(j_preimage, c)

        return out

    def co_reduced(self, y):
        """
        Reduce element `y` of codomain of ``self`` w.r.t. the image of
        ``self``.

        Suppose that ``self`` is a morphism from `X` to `Y`. Then for any
        `y \in Y`, the call ``self.co_reduced(y)`` returns a normal form for
        `y` in the quotient `Y / I` where `I` is the image of ``self``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i,4)  ) # uni-upper
            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y)
            sage: phi.co_reduced(y[1] + y[2])
            0
        """
        G = self.codomain()
        basis_map = self._on_basis
        assert y in G

        result    = G.zero()
        remainder = y

        while not remainder.is_zero():
            (j,c) = self._dominant_item(remainder)
            if self._inverse_on_support is None:
                j_preimage = j
            else:
                j_preimage = self._inverse_on_support(j)
            if j_preimage is None:
                dom_term = G.term(j,c)
                remainder -= dom_term
                result += dom_term
            else:
                s = basis_map(j_preimage)
                assert j == self._dominant_item(s)[0]
                if not self._unitriangular:
                    c /= s[j]
                    remainder -= s._lmul_(c)
        return result

    def co_kernel_projection(self, category = None):
        """
        Return a projection on the co-kernel of ``self``.

        INPUT:

        - ``category`` -- the category of the result

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, [1,2,3]); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1,2,3,4,5]); y = Y.basis()
            sage: uut = lambda i: sum(  y[j] for j in range(i+1,6)  ) # uni-upper
            sage: phi = X.module_morphism(uut, triangular=True, codomain = Y,
            ....:      inverse_on_support=lambda i: i-1 if i in [2,3,4] else None)
            sage: phipro = phi.co_kernel_projection()
            sage: phipro(y[1] + y[2])
            B[1]
            sage: all(phipro(phi(x)).is_zero() for x in X.basis())
            True
            sage: phipro(y[1])
            B[1]
            sage: phipro(y[4])
            -B[5]
            sage: phipro(y[5])
            B[5]
        """
        category = ModulesWithBasis(self.codomain().base_ring()).or_subcategory(category)
        return SetMorphism(Hom(self.codomain(), self.codomain(), category),
                           self.co_reduced)

class DiagonalModuleMorphism(ModuleMorphismByLinearity):
    r"""
    A class for diagonal module morphisms.

    See :meth:`ModulesWithBasis.ParentMethods.module_morphism`.

    INPUT:

    - ``domain``, ``codomain`` -- two modules with basis `F` and `G`,
      respectively
    - ``diagonal`` -- a function `d`

    Assumptions:

    - ``domain`` and ``codomain`` have the same base ring `R`,
    - their respective bases `F` and `G` have the same index set `I`,
    - `d` is a function `I \to R`.

    Return the diagonal module morphism from ``domain`` to ``codomain``
    sending `F(i) \mapsto d(i) G(i)` for all `i \in I`.

    By default, ``codomain`` is currently assumed to be ``domain``.
    (Todo: make a consistent choice with ``*ModuleMorphism``.)

    .. TODO::

        - Implement an optimized ``_call_()`` function.
        - Generalize to a mapcoeffs.
        - Generalize to a mapterms.

    EXAMPLES::

        sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X")
        sage: phi = X.module_morphism(diagonal = factorial, codomain = X)
        sage: x = X.basis()
        sage: phi(x[1]), phi(x[2]), phi(x[3])
        (B[1], 2*B[2], 6*B[3])
    """
    def __init__(self, diagonal, domain, codomain = None, category = None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X")
            sage: phi = X.module_morphism(diagonal = factorial, codomain = X)
            sage: phi.__class__
            <class 'sage.categories.modules_with_basis.DiagonalModuleMorphism'>
            sage: TestSuite(phi).run() # known issue; see ModuleMorphismByLinearity.__init__
            Failure in _test_category:
            ...
            The following tests failed: _test_category
        """
        assert codomain is not None
        assert domain.basis().keys() == codomain.basis().keys()
        assert domain.base_ring()    == codomain.base_ring()
        if category is None:
            category = ModulesWithBasis(domain.base_ring())
        ModuleMorphismByLinearity.__init__(self, domain = domain, codomain = codomain, category = category)
        self._diagonal = diagonal

    def _on_basis(self, i):
        """
        Return the image by ``self`` of the basis element indexed by ``i``.

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); Y.rename("Y"); y = Y.basis()
            sage: phi = X.module_morphism(diagonal = factorial, codomain = X)
            sage: phi._on_basis(3)
            6*B[3]
        """
        return self.codomain().term(i, self._diagonal(i))

    def __invert__(self):
        """
        Return the inverse diagonal morphism.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, [1, 2, 3]); X.rename("X"); x = X.basis()
            sage: Y = CombinatorialFreeModule(QQ, [1, 2, 3]); Y.rename("Y"); y = Y.basis()
            sage: phi = X.module_morphism(diagonal = factorial, codomain = X)
            sage: phi_inv = ~phi
            sage: phi_inv
            Generic endomorphism of Y
            sage: phi_inv(y[3])
            1/6*B[3]

        Caveat: this inverse morphism is only well defined if
        `d(\lambda)` is always invertible in the base ring. This is
        condition is *not* tested for, so using an ill defined inverse
        morphism will trigger arithmetic errors.
        """
        return self.__class__(
            pointwise_inverse_function(self._diagonal),
            domain = self.codomain(), codomain = self.domain(), category = self.category_for())


def pointwise_inverse_function(f):
    r"""
    Return the function `x \mapsto 1 / f(x)`.

    INPUT:

    - ``f`` -- a function

    EXAMPLES::

        sage: from sage.categories.modules_with_basis import pointwise_inverse_function
        sage: def f(x): return x
        ....:
        sage: g = pointwise_inverse_function(f)
        sage: g(1), g(2), g(3)
        (1, 1/2, 1/3)

    :func:`pointwise_inverse_function` is an involution::

        sage: f is pointwise_inverse_function(g)
        True

    .. TODO::

        This has nothing to do here!!! Should there be a library for
        pointwise operations on functions somewhere in Sage?
    """
    if hasattr(f, "pointwise_inverse"):
        return f.pointwise_inverse()
    return PointwiseInverseFunction(f)

from sage.structure.sage_object import SageObject
class PointwiseInverseFunction(SageObject):
    r"""
    A class for pointwise inverse functions.

    The pointwise inverse function of a function `f` is the function
    sending every `x` to `1 / f(x)`.

    EXAMPLES::

        sage: from sage.categories.modules_with_basis import PointwiseInverseFunction
        sage: f = PointwiseInverseFunction(factorial)
        sage: f(0), f(1), f(2), f(3)
        (1, 1, 1/2, 1/6)
    """

    def __eq__(self, other):
        """
        TESTS::

            sage: from sage.categories.modules_with_basis import PointwiseInverseFunction
            sage: f = PointwiseInverseFunction(factorial)
            sage: g = PointwiseInverseFunction(factorial)
            sage: f is g
            False
            sage: f == g
            True
        """
        return self.__class__ is other.__class__ and self.__dict__ == other.__dict__

    def __init__(self, f):
        """
        TESTS::

            sage: from sage.categories.modules_with_basis import PointwiseInverseFunction
            sage: f = PointwiseInverseFunction(factorial)
            sage: f(0), f(1), f(2), f(3)
            (1, 1, 1/2, 1/6)
            sage: TestSuite(f).run()
        """
        self._pointwise_inverse = f

    def __call__(self, *args):
        """
        TESTS::

            sage: from sage.categories.modules_with_basis import PointwiseInverseFunction
            sage: g = PointwiseInverseFunction(operator.mul)
            sage: g(5,7)
            1/35
        """
        return ~(self._pointwise_inverse(*args))

    def pointwise_inverse(self):
        """
        TESTS::

            sage: from sage.categories.modules_with_basis import PointwiseInverseFunction
            sage: g = PointwiseInverseFunction(operator.mul)
            sage: g.pointwise_inverse() is operator.mul
            True
        """
        return self._pointwise_inverse

