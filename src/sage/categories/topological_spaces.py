r"""
Topological Spaces
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
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory
from sage.categories.sets_cat import Sets

class TopologicalSpacesCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "Topological"

    # Currently no use case for this
    @lazy_class_attribute
    def _base_category_class(cls):
        """
        Recover the class of the base category.

        OUTPUT:

        A *tuple* whose first entry is the base category class.

        .. WARNING::

            This is only used for graded categories that are not
            implemented as nested classes, and won't work otherwise.

        .. SEEALSO:: :meth:`__classcall__`

        EXAMPLES::

            sage: GradedModules._base_category_class
            (<class 'sage.categories.modules.Modules'>,)
            sage: GradedAlgebrasWithBasis._base_category_class
            (<class 'sage.categories.algebras_with_basis.AlgebrasWithBasis'>,)

        The reason for wrapping the base category class in a tuple is
        that, often, the base category class implements a
        :meth:`__classget__` method which would get in the way upon
        attribute access::

            sage: F = GradedAlgebrasWithBasis
            sage: F._foo = F._base_category_class[0]
            sage: F._foo
            Traceback (most recent call last):
            ...
            AssertionError: base category class for <...AlgebrasWithBasis'> mismatch;
            expected <...Algebras'>, got <...GradedAlgebrasWithBasis'>
        """
        module_name = cls.__module__.replace("topological_","")
        import sys
        name   = cls.__name__.replace("Topological","")
        __import__(module_name)
        module = sys.modules[module_name]
        return (module.__dict__[name],)

    @staticmethod
    def __classcall__(cls, category=None, *args):
        """
        Magic support for putting topological categories in their own file.

        EXAMPLES::

            sage: TopologicalSpaces()   # indirect doctest
            Category of graded modules over Integer Ring
            sage: Sets().Topological()
            Category of graded modules over Integer Ring
            sage: TopologicalSpaces() is Sets().Topological()
            True

        .. TODO::

            Generalize this support for all other functorial
            constructions if at some point we have a category ``Blah`` for
            which we want to implement the construction ``Blah.Foo`` in a
            separate file like we do for e.g. :class:`GradedModules`,
            :class:`GradedAlgebras`, ...

        .. SEEALSO:: :meth:`_base_category_class`
        """
        base_category_class = cls._base_category_class[0]
        if isinstance(category, base_category_class):
            return super(TopologicalSpacesCategory, cls).__classcall__(cls, category, *args)
        else:
            return base_category_class(category, *args).Topological()

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Groups().Topological()  # indirect doctest
            Category of graded algebras with basis over Rational Field
        """
        return "topological {}".format(self.base_category()._repr_object_names())

class TopologicalSpaces(TopologicalSpacesCategory):
    """
    The category of topological spaces.

    EXAMPLES::

        sage: Sets().Topological()
        Category of topological spaces
        sage: Sets().Topological().super_categories()
        [Category of modules over Integer Ring]

    The category of topological spaces defines the topological structure,
    which shall be preserved by morphisms::

        sage: Sets().Topological().additional_structure()
        Category of topological spaces

    TESTS::

        sage: TestSuite(Sets().Topological()).run()
    """
    # We must override the general object because the names don't match
    _base_category_class = (Sets,)

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Sets().Topological()  # indirect doctest
            Category of topological spaces
        """
        return "topological spaces"

