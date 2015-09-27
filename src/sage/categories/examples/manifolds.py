"""
Examples of manifolds
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.manifolds import Manifolds
from sage.rings.all import QQ

class Plane(UniqueRepresentation, Parent):
    r"""
    An example of a manifold: the `n`-dimensional plane.

    This class illustrates a minimal implementation of a manifold.

    EXAMPLES::

        sage: from sage.categories.manifolds import Manifolds
        sage: M = Manifolds().example(); M
        An example of a manifold: the 3-dimensional plane

        sage: M.category()
        Category of manifolds

    We conclude by running systematic tests on this manifold::

        sage: TestSuite(M).run()
    """

    def __init__(self, n=3):
        r"""
        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: M = Manifolds().example(6); M
            An example of a manifold: the 6-dimensional plane

        TESTS::

            sage: TestSuite(M).run()
        """
        self._n = n
        Parent.__init__(self, category=Manifolds())

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds().example()
            An example of a manifold: the 3-dimensional plane
        """
        return "An example of a manifold: the {}-dimensional plane".format(self._n)

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: M = Manifolds().example()
            sage: M.dimension()
            3
        """
        return self._n

    def an_element(self):
        r"""
        Return an element of the manifold, as per
        :meth:`Sets.ParentMethods.an_element`.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: M = Manifolds().example()
            sage: M.an_element()
            (0, 0, 0)
        """
        zero = QQ.zero()
        return self(tuple([zero]*self._n))

    Element = ElementWrapper

Example = Plane

