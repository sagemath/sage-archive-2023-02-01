"""
Subquotient Functorial Construction

AUTHORS:

 - Nicolas M. Thiery (2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category import Category
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory

# This is Category.Subquotients
def Subquotients(self):
    r"""
    INPUT:
     - ``self`` -- a concrete category

    Given a concrete category ``self == As()`` (i.e. a subcategory of
    ``Sets()``), ``As().Subquotients()`` returns the category of objects
    of ``As()`` endowed with a distinguished description as
    subquotient of some other object of ``As()``.

    EXAMPLES::

        sage: Monoids().Subquotients()
        Category of subquotients of monoids

    A parent `A` in ``As()`` is further in ``As().Subquotients()`` if
    there is a distinguished parent `B` in ``As()``, called the
    *ambient space*, a subspace `B'` of `B` and a pair of *structure
    preserving* maps:

    .. math::

        l: A \mapsto B'  \text{ and }  r: B' \mapsto A

    called respectively the *lifting map* and *retract map* such that
    `r \circ l` is the identity of `A`. What exactly *structure
    preserving* means is explicited in each category; this typically
    states that, for each operation `op` of the category, there is a
    commutative diagram such that:

        for all `e\in A`, one has `op_A(e) = r(op_B(l(e)))`

    This allows for deriving the operations on `A` from those on `B`.

    Note: this is a slightly weaker definition than that found on
    http://en.wikipedia.org/wiki/Subquotient: B' is not necessarily
    required to be a subobject of B.

    Assumptions:

    - For any category ``As()``, ``As().Subquotients()`` is a
      subcategory of ``As()``.

      Example: a subquotient of a group is a group (e.g. a left or right
      quotients of a group by a non normal subgroup is not in this category).

    - This construction is covariant: if ``As()`` is a subcategory of
      ``Bs()``, then ``As().Subquotients()`` is a subcategory of
      ``Bs().Subquotients()``

      Example: if `A` is a distinguished subquotient of `B` in the category of
      groups, then is is also a subquotient of `B` in the category of monoids.

    - If the user (or a program) calls ``As().Subquotients()``, then it is
      assumed that subquotients are well defined in this category. This is not
      checked, and probably never will. Note that, if a category `As()` does
      not specify anything about its subquotients, then it's subquotient
      category looks like this::

         sage: EuclideanDomains().Subquotients()
         Join of Category of euclidean domains and Category of subquotients of monoids

    Interface: the ambient space of `B` is given by ``B.ambient()``. The
    lifting and retract map are implemented respectively as methods
    ``B.lift(b)`` and ``B.retract(a)``. As a shorthand, one can use
    alternatively ``b.lift()``::

        sage: S = Semigroups().Subquotients().example(); S
        An example of a (sub)quotient semigroup: a quotient of the left zero semigroup
        sage: S.ambient()
        An example of a semigroup: the left zero semigroup
        sage: S(3).lift().parent()
        An example of a semigroup: the left zero semigroup
        sage: S(3) * S(1) == S.retract( S(3).lift() * S(1).lift() )
        True

    See ``S?`` for more.

    TODO: use a more interesting example, like `\ZZ/n\ZZ`.

    The two most common use cases are:

     - *quotients*, when `A'=A` and `r` is a morphism; then `r` is a
       canonical quotient map from `A` to `B`)
     - *subobjects* (when `l` is an embedding from `B` into `A`).

    See respectively :class:`~sage.categories.quotients.Quotients` and
    :class:`~sage.categories.subobjects.Subobjects`.

    TESTS::

        sage: TestSuite(Sets().Subquotients()).run()

    """
    return SubquotientsCategory.category_of(self)

Category.Subquotients = Subquotients

class SubquotientsCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "Subquotients"
