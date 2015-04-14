r"""
Filtered Modules

A *filtered module* over a ring `R` with a totally ordered
indexing set `I` (typically `I = \NN`) is an `R`-module `M` equipped
with a family `(F_i)_{i \in I}` of `R`-submodules satisfying
`F_i \subseteq F_j` for all `i,j \in I` having `i \leq j`, and
`M = \bigcup_{i \in I} F_i`. This family is called a *filtration*
of the given module `M`.

.. TODO::

    Implement a notion for decreasing filtrations: where `F_j \subseteq F_i`
    when `i \leq j`.

.. TODO::

    Implement filtrations for all concrete categories.

.. TODO::

    Implement `\operatorname{gr}` as a functor.
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory

class FilteredModulesCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = Algebras(QQ).Filtered()
            sage: C
            Category of filtered algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of algebras over Rational Field,
             Category of filtered modules over Rational Field]

            sage: AlgebrasWithBasis(QQ).Filtered().base_ring()
            Rational Field
            sage: HopfAlgebrasWithBasis(QQ).Filtered().base_ring()
            Rational Field
        """
        super(FilteredModulesCategory, self).__init__(base_category, base_category.base_ring())

    _functor_category = "Filtered"

    @lazy_class_attribute
    def _base_category_class(cls):
        """
        Recover the class of the base category.

        OUTPUT:

        A *tuple* whose first entry is the base category class.

        .. WARNING::

            This is only used for filtered categories that are not
            implemented as nested classes, and won't work otherwise.

        .. SEEALSO:: :meth:`__classcall__`

        EXAMPLES::

            sage: from sage.categories.filtered_modules import FilteredModules
            sage: FilteredModules._base_category_class
            (<class 'sage.categories.modules.Modules'>,)
            sage: from sage.categories.filtered_algebras_with_basis import FilteredAlgebrasWithBasis
            sage: FilteredAlgebrasWithBasis._base_category_class
            (<class 'sage.categories.algebras_with_basis.AlgebrasWithBasis'>,)

        The reason for wrapping the base category class in a tuple is
        that, often, the base category class implements a
        :meth:`__classget__` method which would get in the way upon
        attribute access::

            sage: F = FilteredAlgebrasWithBasis
            sage: F._foo = F._base_category_class[0]
            sage: F._foo
            Traceback (most recent call last):
            ...
            AssertionError: base category class for <...AlgebrasWithBasis'> mismatch;
            expected <...Algebras'>, got <...FilteredAlgebrasWithBasis'>
        """
        module_name = cls.__module__.replace("filtered_","")
        import sys
        name   = cls.__name__.replace("Filtered","")
        __import__(module_name)
        module = sys.modules[module_name]
        return (module.__dict__[name],)

    @staticmethod
    def __classcall__(cls, category, *args):
        """
        Magic support for putting Filtered categories in their own file.

        EXAMPLES::

            sage: from sage.categories.filtered_modules import FilteredModules
            sage: FilteredModules(ZZ)   # indirect doctest
            Category of filtered modules over Integer Ring
            sage: Modules(ZZ).Filtered()
            Category of filtered modules over Integer Ring
            sage: FilteredModules(ZZ) is Modules(ZZ).Filtered()
            True

        .. TODO::

            Generalize this support for all other functorial
            constructions if at some point we have a category ``Blah`` for
            which we want to implement the construction ``Blah.Foo`` in a
            separate file like we do for e.g. :class:`FilteredModules`,
            :class:`FilteredAlgebras`, ...

        .. SEEALSO:: :meth:`_base_category_class`
        """
        base_category_class = cls._base_category_class[0]
        if isinstance(category, base_category_class):
            return super(FilteredModulesCategory, cls).__classcall__(cls, category, *args)
        else:
            return base_category_class(category, *args).Filtered()

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Filtered()  # indirect doctest
            Category of filtered algebras with basis over Rational Field
        """
        return "filtered {}".format(self.base_category()._repr_object_names())

class FilteredModules(FilteredModulesCategory):
    r"""
    The category of filtered modules over a given ring `R`.

    A *filtered module* over a ring `R` with a totally ordered
    indexing set `I` (typically `I = \NN`) is an `R`-module `M` equipped
    with a family `(F_i)_{i \in I}` of `R`-submodules satisfying
    `F_i \subseteq F_j` for all `i,j \in I` having `i \leq j`, and
    `M = \bigcup_{i \in I} F_i`. This family is called a *filtration*
    of the given module `M`.

    EXAMPLES::

        sage: Modules(ZZ).Filtered()
        Category of filtered modules over Integer Ring
        sage: Modules(ZZ).Filtered().super_categories()
        [Category of modules over Integer Ring]

    TESTS::

        sage: TestSuite(Modules(ZZ).Filtered()).run()

    REFERENCES:

    - :wikipedia:`Filtration_(mathematics)`
    """

    def extra_super_categories(self):
        r"""
        Add :class:`VectorSpaces` to the super categories of ``self`` if
        the base ring is a field.

        EXAMPLES::

            sage: Modules(QQ).Filtered().extra_super_categories()
            [Category of vector spaces over Rational Field]
            sage: Modules(ZZ).Filtered().extra_super_categories()
            []

        This makes sure that ``Modules(QQ).Filtered()`` returns an
        instance of :class:`FilteredModules` and not a join category of
        an instance of this class and of ``VectorSpaces(QQ)``::

            sage: type(Modules(QQ).Filtered())
            <class 'sage.categories.filtered_modules.FilteredModules_with_category'>

        .. TODO::

            Get rid of this workaround once there is a more systematic
            approach for the alias ``Modules(QQ)`` -> ``VectorSpaces(QQ)``.
            Probably the latter should be a category with axiom, and
            covariant constructions should play well with axioms.
        """
        from sage.categories.modules import Modules
        from sage.categories.fields import Fields
        base_ring = self.base_ring()
        if base_ring in Fields:
            return [Modules(base_ring)]
        else:
            return []

    class SubcategoryMethods:

        @cached_method
        def Connected(self):
            r"""
            Return the full subcategory of the connected objects of ``self``.

            A filtered `R`-module `M` with filtration
            `(F_0, F_1, F_2, \ldots)` (indexed by `\NN`)
            is said to be *connected* if `F_0` is isomorphic
            to `R`.

            EXAMPLES::

                sage: Modules(ZZ).Filtered().Connected()
                Category of filtered connected modules over Integer Ring
                sage: Coalgebras(QQ).Filtered().Connected()
                Join of Category of filtered connected modules over Rational Field
                    and Category of coalgebras over Rational Field
                sage: AlgebrasWithBasis(QQ).Filtered().Connected()
                Category of filtered connected algebras with basis over Rational Field

            TESTS::

                sage: TestSuite(Modules(ZZ).Filtered().Connected()).run()
                sage: Coalgebras(QQ).Filtered().Connected.__module__
                'sage.categories.filtered_modules'
            """
            return self._with_axiom("Connected")

    class Connected(CategoryWithAxiom_over_base_ring):
        pass

    class ParentMethods:
        pass

    class ElementMethods:
        pass

