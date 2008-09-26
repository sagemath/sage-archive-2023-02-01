"""
Functors

AUTHORS:
    - David Kohel and William Stein
    - David Joyner (2005-12-17): examples
    - Robert Bradshaw (2007-06-23): Pyrexify
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import category

cdef class Functor(SageObject):
    """
    EXAMPLES:
        sage: rings  = Rings()
        sage: abgrps = AbelianGroups()
        sage: F = ForgetfulFunctor(rings, abgrps)
        sage: F.domain()
        Category of rings
        sage: F.codomain()
        Category of abelian groups
        sage: from sage.categories.functor import is_Functor
        sage: is_Functor(F)
        True
        sage: I = IdentityFunctor(abgrps)
        sage: I
        The identity functor on AbelianGroups
        sage: I.domain()
        Category of abelian groups
        sage: is_Functor(I)
        True
    """
    def __init__(self, domain, codomain):
        if not category.is_Category(domain):
            raise TypeError, "domain (=%s) must be a category"%domain
        if not category.is_Category(codomain):
            raise TypeError, "codomain (=%s) must be a category"%codomain
        self.__domain = domain
        self.__codomain = codomain

    def __repr__(self):
        return "Functor from %s to %s"%(self.__domain.short_name(), self.__codomain.short_name())

    def __call__(self, x):
        if not (x in  self.__domain):
            try:
                x = self.__domain(x)
            except TypeError:
                raise TypeError, "x (=%s) must be coercible to an object of %s"%(x, self.__domain)
        y = self._apply_functor(x)
        if not (y in self.__codomain):
            raise TypeError, "The functor %s is ill-defined, since it sends x (=%s) to something that is not an object of %s"%(x, self.__codomain)

    def domain(self):
        return self.__domain

    def codomain(self):
        return self.__codomain


def is_Functor(x):
    return isinstance(x, Functor)


###########################################
# The natural functors in SAGE
###########################################

class ForgetfulFunctor_generic(Functor):
    def __reduce__(self):
        """
        EXAMPLES:
            sage: F = ForgetfulFunctor(Groups(), Sets())
            sage: loads(F.dumps()) == F
            True
        """
        return ForgetfulFunctor, (self.domain(), self.codomain())

    def __repr__(self):
        return "The forgetful functor from %s to %s"%(
            self.domain().short_name(), self.codomain().short_name())

    def _apply_functor(self, x):
        return self.codomain()(x)

    def __cmp__(self, other):
        if self.domain() == other.domain() and \
           self.codomain() == other.codomain():
            return 0
        return -1


class IdentityFunctor_generic(ForgetfulFunctor_generic):
    def __init__(self, C):
        ForgetfulFunctor_generic.__init__(self, C, C)

    def __reduce__(self):
        """
        EXAMPLES:
            sage: F = IdentityFunctor(Groups())
            sage: loads(F.dumps()) == F
            True
        """
        return IdentityFunctor, (self.domain(), )

    def __repr__(self):
        return "The identity functor on %s"%(self.domain().short_name())

    def _apply_functor(self, x):
        return x

    def __call__(self, x):  # no
        if not x in self.domain():
            raise TypeError, "x (=%s) must be in the category %s"%(x,self.domain())
        return x

def IdentityFunctor(C):
    return IdentityFunctor_generic(C)

def ForgetfulFunctor(domain, codomain):
    """
    Construct the forgetful function from one category to another.

    EXAMPLES:
        sage: rings = Rings()
        sage: abgrps = AbelianGroups()
        sage: F = ForgetfulFunctor(rings, abgrps)
        sage: F
        The forgetful functor from Rings to AbelianGroups
    """
    from category_types import category_hierarchy
    if domain == codomain:
        return IdentityFunctor(domain, codomain)
    if not domain.is_subcategory(codomain):
        raise ValueError, "Forgetful functor not supported for domain %s"%domain
    return ForgetfulFunctor_generic(domain, codomain)


