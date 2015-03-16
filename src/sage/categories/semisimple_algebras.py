r"""
Semisimple Algebras
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from algebras import Algebras
from sage.misc.cachefunc import cached_method
import operator

class SemisimpleAlgebras(Category_over_base_ring):
    """
    The category of semisimple algebras over a given base ring.

    EXAMPLES::

        from sage.categories.semisimple_algebras import SemisimpleAlgebras
        sage: C = SemisimpleAlgebras(QQ); C
        Category of semisimple algebras over Rational Field

    This category is best constructed as::

        sage: D = Algebras(QQ).Semisimple(); D
        Category of semisimple algebras over Rational Field
        sage: D is C
        True

        sage: C.super_categories()
        [Category of algebras over Rational Field]

    Typically, finite group algebras are semisimple::

        sage: DihedralGroup(5).algebra(QQ) in SemisimpleAlgebras
        True

    Unless the characteristic of the field divides the order of the group::

        sage: DihedralGroup(5).algebra(IntegerModRing(5)) in SemisimpleAlgebras
        False

        sage: DihedralGroup(5).algebra(IntegerModRing(7)) in SemisimpleAlgebras
        True

    .. SEEALSO:: `<http://en.wikipedia.org/wiki/Semisimple_algebra>`_

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: SemisimpleAlgebras(QQ).super_categories()
            [Category of algebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R)]
