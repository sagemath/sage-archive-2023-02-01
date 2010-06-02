"""
Dual functorial construction

AUTHORS:

 - Nicolas M. Thiery (2009-2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from category import Category

# could do SelfDualCategory

from sage.categories.category import Category
from sage.categories.covariant_functorial_construction import CovariantFunctorialConstruction, CovariantConstructionCategory
from sage.categories.category_types import Category_over_base_ring
import sage.structure.parent
import sage.structure.element

# This is Category.DualObjects
def DualObjects(self):
    """
    Returns the category of duals of objects of ``self``.

    INPUT:

     - ``self`` -- a subcategory of vector spaces over some base ring

    The dual of a vector space `V` is the space consisting of all
    linear functionals on `V` (http://en.wikipedia.org/wiki/Dual_space).
    Additional structure on `V` can endow its dual with additional
    structure; e.g. if `V` is an algebra, then its dual is a
    coalgebra.

    This returns the category of dual of spaces in ``self`` endowed
    with the appropriate additional structure.

    See also
    :class:`~sage.categories.covariant_functorial_construction.CovariantFunctorialConstruction`.

    TODO: add support for graded duals.

    EXAMPLES::

        sage: VectorSpaces(QQ).DualObjects()
        Category of duals of vector spaces over Rational Field

    The dual of a vector space is a vector space::

        sage: VectorSpaces(QQ).DualObjects().super_categories()
        [Category of vector spaces over Rational Field]

    The dual of an algebra space is a coalgebra::

        sage: Algebras(QQ).DualObjects().super_categories()
        [Category of coalgebras over Rational Field, Category of duals of vector spaces over Rational Field]

    The dual of a coalgebra space is an algebra::

        sage: Coalgebras(QQ).DualObjects().super_categories()
        [Category of algebras over Rational Field, Category of duals of vector spaces over Rational Field]

    As a shorthand, this category can be accessed with the
    :meth:`dual` method::

        sage: VectorSpaces(QQ).dual()
        Category of duals of vector spaces over Rational Field

    TESTS::

        sage: C = VectorSpaces(QQ).DualObjects()
        sage: C.base_category()
        Category of vector spaces over Rational Field
        sage: C.super_categories()
        [Category of vector spaces over Rational Field]
        sage: latex(C)
        \mathbf{DualObjects}(\mathbf{VectorSpaces}_{\Bold{Q}})
        sage: TestSuite(C).run()
    """
    return DualObjectsCategory.category_of(self)

Category.DualObjects = DualObjects
Category.dual = DualObjects

class DualFunctor(CovariantFunctorialConstruction):
    """
    A singleton class for the dual functor
    """
    _functor_name = "dual"
    _functor_category = "DualObjects"
    symbol = "^*"

class DualObjectsCategory(CovariantConstructionCategory):

    _functor_category = "DualObjects"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ).DualObjects() # indirect doctest
            Category of duals of vector spaces over Rational Field
        """
        # Just to remove the `objects`
        return "duals of %s"%(self.base_category()._repr_object_names())
