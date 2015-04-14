r"""
Super modules
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory
from sage.categories.modules import Modules

class SuperModulesCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = Algebras(QQ).Super()
            sage: C
            Category of super algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of graded algebras over Rational Field,
             Category of super modules over Rational Field]

            sage: AlgebrasWithBasis(QQ).Super().base_ring()
            Rational Field
            sage: HopfAlgebrasWithBasis(QQ).Super().base_ring()
            Rational Field
        """
        super(SuperModulesCategory, self).__init__(base_category, base_category.base_ring())

    _functor_category = "Super"

    @lazy_class_attribute
    def _base_category_class(cls):
        """
        Recover the class of the base category.

        OUTPUT:

        A *tuple* whose first entry is the base category class.

        .. WARNING::

            This is only used for super categories that are not
            implemented as nested classes, and won't work otherwise.

        .. SEEALSO:: :meth:`__classcall__`

        EXAMPLES::

            sage: from sage.categories.super_modules import SuperModules
            sage: from sage.categories.super_algebras_with_basis import SuperAlgebrasWithBasis
            sage: SuperModules._base_category_class
            (<class 'sage.categories.modules.Modules'>,)
            sage: SuperAlgebrasWithBasis._base_category_class
            (<class 'sage.categories.algebras_with_basis.AlgebrasWithBasis'>,)

        The reason for wrapping the base category class in a tuple is
        that, often, the base category class implements a
        :meth:`__classget__` method which would get in the way upon
        attribute access::

            sage: F = SuperAlgebrasWithBasis
            sage: F._foo = F._base_category_class[0]
            sage: F._foo
            Traceback (most recent call last):
            ...
            AssertionError: base category class for <...AlgebrasWithBasis'> mismatch;
            expected <...Algebras'>, got <...SuperAlgebrasWithBasis'>
        """
        module_name = cls.__module__.replace("super_","")
        import sys
        name   = cls.__name__.replace("Super","")
        __import__(module_name)
        module = sys.modules[module_name]
        return (module.__dict__[name],)

    @staticmethod
    def __classcall__(cls, category, *args):
        """
        Magic support for putting Super categories in their own file.

        EXAMPLES::

            sage: from sage.categories.super_modules import SuperModules
            sage: SuperModules(ZZ)   # indirect doctest
            Category of super modules over Integer Ring
            sage: Modules(ZZ).Super()
            Category of super modules over Integer Ring
            sage: SuperModules(ZZ) is Modules(ZZ).Super()
            True

        .. TODO::

            Generalize this support for all other functorial
            constructions if at some point we have a category ``Blah`` for
            which we want to implement the construction ``Blah.Foo`` in a
            separate file like we do for e.g. :class:`SuperModules`,
            :class:`SuperAlgebras`, ...

        .. SEEALSO:: :meth:`_base_category_class`
        """
        base_category_class = cls._base_category_class[0]
        if isinstance(category, base_category_class):
            return super(SuperModulesCategory, cls).__classcall__(cls, category, *args)
        else:
            return base_category_class(category, *args).Super()

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Super()  # indirect doctest
            Category of super algebras with basis over Rational Field
        """
        return "super {}".format(self.base_category()._repr_object_names())

    @classmethod
    def default_super_categories(cls, category, *args):
        """
        Returns the default super categories of ``category.Subobjects()``

        Mathematical meaning: if `A` is a super version of `B`,
        then `A` is also a graded version of `B`.

        INPUT:

         - ``cls`` -- the class ``SubobjectsCategory``
         - ``category`` -- a category `Cat`

        OUTPUT: a (join) category

        In practice, this returns ``category.Subquotients()``, joined
        together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories() <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category`` and ``cat.Subobjects()`` for
        each ``cat`` in the super categories of ``category``).

        EXAMPLES::
        """
        return Category.join([category.Graded(), super(SuperModulesCategory, cls).default_super_categories(category, *args)])

class SuperModules(SuperModulesCategory):
    """
    The category of super modules.

    EXAMPLES::

        sage: Modules(ZZ).Super()
        Category of super modules over Integer Ring
        sage: Modules(ZZ).Super().super_categories()
        [Category of graded modules over Integer Ring]

    The category of super modules defines the super structure which
    shall be preserved by morphisms::

        sage: Modules(ZZ).Super().additional_structure()
        Category of super modules over Integer Ring

    TESTS::

        sage: TestSuite(Modules(ZZ).Super()).run()
    """
    def extra_super_categories(self):
        r"""
        Adds :class:`VectorSpaces` to the super categories of ``self`` if
        the base ring is a field.

        EXAMPLES::

            sage: Modules(QQ).Super().extra_super_categories()
            [Category of vector spaces over Rational Field]
            sage: Modules(ZZ).Super().extra_super_categories()
            []

        This makes sure that ``Modules(QQ).Super()`` returns an
        instance of :class:`GradedModules` and not a join category of
        an instance of this class and of ``VectorSpaces(QQ)``::

            sage: type(Modules(QQ).Super())
            <class 'sage.categories.super_modules.SuperModules_with_category'>

        .. TODO::

            Get rid of this workaround once there is a more systematic
            approach for the alias ``Modules(QQ)`` -> ``VectorSpaces(QQ)``.
            Probably the later should be a category with axiom, and
            covariant constructions should play well with axioms.
        """
        from sage.categories.modules import Modules
        from sage.categories.fields import Fields
        base_ring = self.base_ring()
        if base_ring in Fields:
            return [Modules(base_ring)]
        else:
            return []

    class ParentMethods:
        pass

    class ElementMethods:
        def is_even_odd(self):
            """
            Return ``0`` if ``self`` is an even element or ``1``
            if an odd element.

            EXAMPLES::

                sage: cat = Algebras(QQ).WithBasis().Super()
                sage: C = CombinatorialFreeModule(QQ, Partitions(), category=cat)
                sage: C.degree_on_basis = sum
                sage: C.basis()[2,2,1].is_even_odd()
                1
                sage: C.basis()[2,2].is_even_odd()
                0
            """
            return self.degree() % 2

        def is_even(self):
            """
            Return if ``self`` is an even element.

            EXAMPLES::

                sage: cat = Algebras(QQ).WithBasis().Super()
                sage: C = CombinatorialFreeModule(QQ, Partitions(), category=cat)
                sage: C.degree_on_basis = sum
                sage: C.basis()[2,2,1].is_even()
                False
                sage: C.basis()[2,2].is_even()
                True
            """
            return self.is_even_odd() == 0

        def is_odd(self):
            """
            Return if ``self`` is an odd element.

            EXAMPLES::

                sage: cat = Algebras(QQ).WithBasis().Super()
                sage: C = CombinatorialFreeModule(QQ, Partitions(), category=cat)
                sage: C.degree_on_basis = sum
                sage: C.basis()[2,2,1].is_odd()
                True
                sage: C.basis()[2,2].is_odd()
                False
            """
            return self.is_even_odd() == 1

