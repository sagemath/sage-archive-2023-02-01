r"""
Ribbons
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from ribbon_shaped_tableau import RibbonShapedTableau, StandardRibbonShapedTableaux

def Ribbon(r):
    """
    Deprecated in :trac:`14101`. Use :class:`RibbonShapedTableau` instead.

    EXAMPLES::

        sage: sage.combinat.ribbon.Ribbon([[2,3],[1,4,5]])
        doctest:...: DeprecationWarning: this class is deprecated. Use RibbonShapedTableau instead
        See http://trac.sagemath.org/14101 for details.
        [[None, None, 2, 3], [1, 4, 5]]
    """
    from sage.misc.superseded import deprecation
    deprecation(14101, 'this class is deprecated. Use RibbonShapedTableau instead')
    return RibbonShapedTableau(r)

def StandardRibbons(shape=None):
    """
    Deprecated in :trac:`14101`. Use :class:`RibbonShapedTableaux` instead.

    EXAMPLES::

        sage: sage.combinat.ribbon.StandardRibbons([3,3,1])
        doctest:...: DeprecationWarning: this class is deprecated. Use StandardRibbonShapedTableaux instead
        See http://trac.sagemath.org/14101 for details.
        Standard ribbon tableaux of shape [3, 3, 1]
    """
    from sage.misc.superseded import deprecation
    deprecation(14101, 'this class is deprecated. Use StandardRibbonShapedTableaux instead')
    return StandardRibbonShapedTableaux(shape)

def from_shape_and_word(shape, word):
    """
    This is deprecated in :trac:`14101`. Use instead
    :meth:`StandardRibbonShapedTableaux.from_shape_and_word()`.

    EXAMPLES::

        sage: sage.combinat.ribbon.from_shape_and_word([2,3],[1,2,3,4,5])
        doctest:...: DeprecationWarning: this function is deprecated. Use StandardRibbonShapedTableaux().from_shape_and_word instead
        See http://trac.sagemath.org/14101 for details.
        [[None, None, 1, 2], [3, 4, 5]]
    """
    from sage.misc.superseded import deprecation
    deprecation(14101,'this function is deprecated. Use StandardRibbonShapedTableaux().from_shape_and_word instead')
    return StandardRibbonShapedTableaux().from_shape_and_word(shape, word)

def from_permutation(p):
    """
    This is deprecated in :trac:`14101`. Use instead
    :meth:`StandardRibbonShapedTableaux.from_shape_and_word()`.

    EXAMPLES::

        sage: sage.combinat.ribbon.from_permutation(Permutation([1, 2, 3]))
        doctest:...: DeprecationWarning: this function is deprecated. Use StandardRibbonShapedTableaux().from_permutation instead
        See http://trac.sagemath.org/14101 for details.
        [[1, 2, 3]]
    """
    from sage.misc.superseded import deprecation
    deprecation(14101,'this function is deprecated. Use StandardRibbonShapedTableaux().from_permutation instead')
    return StandardRibbonShapedTableaux().from_permutation(p)

def StandardRibbons_shape(shape):
    """
    This is deprecated in :trac:`14101`. Use instead
    :class:`StandardRibbonShapedTableaux`.

    EXAMPLES::

        sage: sage.combinat.ribbon.StandardRibbons_shape([3,3,1])
        doctest:...: DeprecationWarning: this class is deprecated. Use StandardRibbonShapedTableaux instead
        See http://trac.sagemath.org/14101 for details.
        Standard ribbon tableaux of shape [3, 3, 1]
    """
    from sage.misc.superseded import deprecation
    deprecation(14101,'this class is deprecated. Use StandardRibbonShapedTableaux instead')
    return StandardRibbonShapedTableaux(shape)

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.ribbon', 'Ribbon_class', RibbonShapedTableau)
register_unpickle_override('sage.combinat.ribbon', 'StandardRibbons_shape', StandardRibbonShapedTableaux)

