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


class SemisimpleAlgebras(Category_over_base_ring):
    """
    The category of semisimple algebras over a given base ring.

    EXAMPLES::

        sage: SemisimpleAlgebras(QQ)
        Category of semisimple algebras over Rational Field
        sage: SemisimpleAlgebras(QQ).super_categories()
        [Category of algebras over Rational Field]

    Typically, finite group algebras are semisimple::

        sage: DihedralGroup(5).algebra(QQ) in SemisimpleAlgebras
        True

    Unless the characteristic of the field divides the order of the group::

        sage: DihedralGroup(5).algebra(IntegerModRing(5)) in SemisimpleAlgebras
        False

        sage: DihedralGroup(5).algebra(IntegerModRing(7)) in SemisimpleAlgebras # todo: not implemented
        True

    .. seealso:: `<http://en.wikipedia.org/wiki/Semisimple_algebra>`_

    TESTS::

        sage: TestSuite(SemisimpleAlgebras(QQ)).run()
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
