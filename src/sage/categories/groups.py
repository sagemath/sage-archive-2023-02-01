r"""
Groups
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

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.monoids import Monoids

class Groups(Category):
    """
    The category of (multiplicative) groups, i.e. monoids with
    inverses.

    EXAMPLES::

        sage: Groups()
        Category of groups
        sage: Groups().super_categories()
        [Category of monoids]

    TESTS::

        sage: TestSuite(Groups()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Groups().super_categories()
            [Category of monoids]
        """
        return [Monoids()]

    class ParentMethods:

        def group_generators(self):
            """
            Returns group generators for self.

            This default implementation calls :meth:`.gens`, for
            backward compatibility.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.group_generators()
                [(1,2,3), (2,3,4)]
            """
            return self.gens()

        def _test_inverse(self, **options):
            """
            Run generic tests on the method :meth:`.__invert__`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: G = SymmetricGroup(3)
                sage: G._test_inverse()

            """
            tester = self._tester(**options)
            for x in tester.some_elements():
                tester.assertEquals(x * ~x, self.one())
                tester.assertEquals(~x * x, self.one())

    class ElementMethods:
        ## inv(x), x/y
        pass
