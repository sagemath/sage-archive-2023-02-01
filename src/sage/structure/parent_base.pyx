###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "../ext/stdsage.pxi"

## def make_parent_with_base_v0(_class, _dict, base, has_coerce_map_from):
##     """
##     This should work for any Python class deriving from this, as long
##     as it doesn't implement some screwy __new__() method.
##     """
##     new_object = _class.__new__(_class)
##     if base is None:
##         (<ParentWithBase>new_object)._base = new_object
##     else:
##         (<ParentWithBase>new_object)._base = base
##     (<ParentWithBase>new_object)._has_coerce_map_from = has_coerce_map_from
##     if not _dict is None:
##         new_object.__dict__ = _dict
##     return new_object

def is_ParentWithBase(x):
    """
    Return True if x is a parent object with base.
    """
    return bool(PY_TYPE_CHECK(x, ParentWithBase))

cdef class ParentWithBase(parent.Parent):
    def __init__(self, base):
        self._base = base
        self._has_coerce_map_from = {}

##     def x__reduce__(self):
##         if HAS_DICTIONARY(self):
##             _dict = self.__dict__
##         else:
##             _dict = None
##         if self._base is self:
##             return (make_parent_with_base_v0, (self.__class__, _dict, None, self._has_coerce_map_from))
##         else:
##             return (make_parent_with_base_v0, (self.__class__, _dict, self._base, self._has_coerce_map_from))

    def base_ring(self):
        return self._base

    # Derived class *must* define base_extend.
    def base_extend(self, X):
        raise TypeError, "BUG: the base_extend method must be defined for '%s' (class '%s')"%(
            self, type(self))

    def base(self):
        return self._base


    ############################################################################
    # Homomorphism --
    ############################################################################
    def Hom(self, codomain, cat=None):
        r"""
        self.Hom(codomain, cat=None):

        Return the homspace \code{Hom(self, codomain, cat)} of all
        homomorphisms from self to codomain in the category cat.  The
        default category is \code{self.category()}.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: R.Hom(QQ)
            Set of Homomorphisms from Polynomial Ring in x, y over Rational Field to Rational Field

        Homspaces are defined for very general \sage objects, even elements of familiar rings.
            sage: n = 5; Hom(n,7)
            Set of Morphisms from 5 to 7 in Category of elements of Integer Ring
            sage: z=(2/3); Hom(z,8/1)
            Set of Morphisms from 2/3 to 8 in Category of elements of Rational Field

        This example illustrates the optional third argument:
            sage: QQ.Hom(ZZ, Sets())
            Set of Morphisms from Rational Field to Integer Ring in Category of sets
        """
        try:
            return self._Hom_(codomain, cat)
        except (TypeError, AttributeError):
            pass
        from sage.categories.all import Hom
        return Hom(self, codomain, cat)
