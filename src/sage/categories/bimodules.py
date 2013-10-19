r"""
Bimodules
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

from sage.categories.category import Category, CategoryWithParameters
from sage.categories.left_modules import LeftModules
from sage.categories.right_modules import RightModules

from sage.categories.rings import Rings
_Rings = Rings()

#?class Bimodules(Category_over_base_rng, Category_over_base_rng):
class Bimodules(CategoryWithParameters):
    """
    The category of `(R,S)`-bimodules

    For `R` and `S` rings, a `(R,S)`-bimodule `X` is a left `R`-module
    and right `S`-module such that the left and right actions commute:
    `r*(x*s) = (r*x)*s`.

    EXAMPLES::

        sage: Bimodules(QQ, ZZ)
        Category of bimodules over Rational Field on the left and Integer Ring on the right
        sage: Bimodules(QQ, ZZ).super_categories()
        [Category of left modules over Rational Field, Category of right modules over Integer Ring]
    """

    def __init__(self, left_base, right_base, name=None):
        """
        EXAMPLES::

            sage: C = Bimodules(QQ, ZZ)
            sage: TestSuite(C).run()
        """
        assert left_base  in _Rings, "The left base must be a ring"
        assert right_base in _Rings, "The right base must be a ring"
        self._left_base_ring = left_base
        self._right_base_ring = right_base
        Category.__init__(self, name)

    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... classes depend on.

        Since :trac:`11935`, the element and parent classes of a
        bimodule only depend on the categories of the left and right
        base ring.

        .. SEEALSO::

            - :meth:`CategoryWithParameters`
            - :meth:`CategoryWithParameters._make_named_class_key`

        EXAMPLES::

            sage: Bimodules(QQ,ZZ)._make_named_class_key('parent_class')
            (Category of quotient fields, Category of euclidean domains)
        """
        return (self._left_base_ring.category(), self._right_base_ring.category())

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class

        EXAMPLES::

            sage: Bimodules.an_instance()
            Category of bimodules over Rational Field on the left and Real Field with 53 bits of precision on the right
        """
        from sage.rings.all import QQ, RR
        return cls(QQ, RR)

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Bimodules(QQ, ZZ) # indirect doctest
            Category of bimodules over Rational Field on the left and Integer Ring on the right
        """
        return "bimodules over %s on the left and %s on the right" \
            %(self._left_base_ring, self._right_base_ring)

    def left_base_ring(self):
        """
        Returns the left base ring over which elements of this category are
        defined.

        EXAMPLES::

            sage: Bimodules(QQ, ZZ).left_base_ring()
            Rational Field
        """
        return self._left_base_ring

    def right_base_ring(self):
        """
        Returns the right base ring over which elements of this category are
        defined.

        EXAMPLES::

            sage: Bimodules(QQ, ZZ).right_base_ring()
            Integer Ring
        """
        return self._right_base_ring

    def _latex_(self):
        """
        EXAMPLES::

            sage: print Bimodules(QQ, ZZ)._latex_()
            {\mathbf{Bimodules}}_{\Bold{Q}}_{\Bold{Z}}
        """
        from sage.misc.latex import latex
        return "{%s}_{%s}_{%s}"%(Category._latex_(self), latex(self._left_base_ring), latex(self._right_base_ring))

    def super_categories(self):
        """
        EXAMPLES::

            sage: Bimodules(QQ, ZZ).super_categories()
            [Category of left modules over Rational Field, Category of right modules over Integer Ring]
        """
        R = self.left_base_ring()
        S = self.right_base_ring()
        return [LeftModules(R), RightModules(S)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
