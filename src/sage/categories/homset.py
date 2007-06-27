"""
Homsets

AUTHORS:
    -- David Kohel and William Stein
    -- David Joyner (2005-12-17): added examples
    -- William Stein (2006-01-14): Changed from Homspace to Homset.
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>, William Stein <wstein@gmail.com>
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

import weakref

import category
import morphism
from sage.structure.parent import Set_generic
from sage.structure.parent_base import ParentWithBase

_cache = {}
def Hom(X, Y, cat=None):
    """
    Create the space of homomorphisms from X to Y in the category cat.

    INPUT:
        X -- anything
        Y -- anything
        cat -- (optional) category in which the morphisms must be

    OUTPUT:
        a homset in cat

    EXAMPLES:
        sage: V = VectorSpace(QQ,3)
        sage: Hom(V, V)
        Set of Morphisms from Vector space of dimension 3 over Rational
        Field to Vector space of dimension 3 over Rational Field in
        Category of vector spaces over Rational Field
        sage: G = SymmetricGroup(3)
        sage: Hom(G, G)
        Set of Morphisms from SymmetricGroup(3) to SymmetricGroup(3) in Category of groups
        sage: Hom(ZZ, QQ, Sets())
        Set of Morphisms from Integer Ring to Rational Field in Category of sets
    """
    import category_types
    global _cache
    key = (X,Y,cat)
    if _cache.has_key(key):
        H = _cache[key]()
        if H: return H

    if cat is None or (cat is X.category() and cat is Y.category()):
        try:
            H = X._Hom_(Y)
        except AttributeError:
            pass

    cat_X = X.category()
    cat_Y = Y.category()
    if cat is None:
        if cat_X.is_subcategory(cat_Y):
            cat = cat_Y
        elif cat_Y.is_subcategory(cat_X):
            if not (cat is None) and not (cat_X is cat_Y):
                raise ValueError, "No unambiguous category found for Hom from %s to %s."%(X,Y)
            cat = cat_X
        else:
            raise TypeError, "No suitable category found for Hom from %s to %s."%(X,Y)

    elif isinstance(cat, category.Category):
        if not isinstance(cat, category.Category):
            raise TypeError, "Argument cat (= %s) must be a category."%cat
        if not cat_X.is_subcategory(cat) \
               or not cat_Y.is_subcategory(cat):
            raise TypeError, \
                  "Argument cat (= %s) is incompatible with %s and %s."%(cat, X, Y)
    else:
        raise TypeError, "Argument cat (= %s) must be a category."%cat


    # coercing would be incredibly annoying, since the domain and codomain
    # are totally different objects
    #X = cat(X); Y = cat(Y)

    # construct H
    if cat._is_subclass(category_types.HeckeModules):

        from sage.modular.hecke.homspace import HeckeModuleHomspace
        H = HeckeModuleHomspace(X, Y)

    elif cat._is_subclass(category_types.FreeModules):

        from sage.modules.free_module_homspace import FreeModuleHomspace
        H = FreeModuleHomspace(X, Y, cat)

    elif cat._is_subclass(category_types.Rings):

        from sage.rings.homset import RingHomset
        H = RingHomset(X, Y)

    elif cat._is_subclass(category_types.Schemes) or cat._is_subclass(category_types.Schemes_over_base):

        from sage.schemes.generic.homset import SchemeHomset
        H = SchemeHomset(X, Y)

    else:  # default
        if isinstance(X, ParentWithBase):
            H = HomsetWithBase(X, Y, cat)
        else:
            H = Homset(X, Y, cat)

    ##_cache[key] = weakref.ref(H)
    _cache[(X, Y, cat)] = weakref.ref(H)

    return H

def hom(X, Y, f):
    """
    Return Hom(X,Y)(f), where f is data that defines an element of Hom(X,Y).

    EXAMPLES:
        sage: R, x = PolynomialRing(QQ,'x').objgen()
        sage: phi = hom(R, QQ, [2])
        sage: phi(x^2 + 3)
        7
    """
    return Hom(X,Y)(f)

def End(X, cat=None):
    r"""
    Create the set of endomorphisms of X in the category cat.

    INPUT:
        X -- anything
        cat -- (optional) category in which to coerce X

    OUTPUT:
        a set of endomorphisms in cat

    EXAMPLES:
        sage: V = VectorSpace(QQ, 3)
        sage: End(V)
        Set of Morphisms from Vector space of dimension 3 over Rational
        Field to Vector space of dimension 3 over Rational Field in
        Category of vector spaces over Rational Field

        sage: G = SymmetricGroup(3)
        sage: S = End(G); S
        Set of Morphisms from SymmetricGroup(3) to SymmetricGroup(3) in Category of groups
        sage: is_Endset(S)
        True
        sage: S.domain()
        Symmetric group of order 3! as a permutation group

    Homsets are \emph{not} objects in their category.  They are
    currently sets.
        sage: S.category()
        Category of sets
        sage: S.domain().category()
        Category of groups
    """
    return Hom(X,X, cat)

def end(X, f):
    """
    Return End(X)(f), where f is data that defines an element of End(X).

    EXAMPLES:
        sage: R, x = PolynomialRing(QQ,'x').objgen()
        sage: phi = end(R, [x + 1])
        sage: phi
        Ring endomorphism of Univariate Polynomial Ring in x over Rational Field
          Defn: x |--> x + 1
        sage: phi(x^2 + 5)
        x^2 + 2*x + 6
    """
    return End(X)(f)

class Homset(Set_generic):
    """
    The class for collections of morphisms in a category.

    EXAMPLES:
        sage: H = Hom(QQ^2, QQ^3)
        sage: loads(H.dumps()) == H
        True
        sage: E = End(AffineSpace(2, names='x,y'))
        sage: loads(E.dumps()) == E
        True
    """
    def __init__(self, X, Y, cat=None, check=True):
        self.__domain = X
        self.__codomain = Y
        if cat is None:
            cat = X.category()
        self.__category = cat
        if check:
            if not isinstance(cat, category.Category):
                raise TypeError, "cat (=%s) must be a category"%cat
            #if not X in cat:
            #    raise TypeError, "X (=%s) must be in cat (=%s)"%(X, cat)
            #if not Y in cat:
            #    raise TypeError, "Y (=%s) must be in cat (=%s)"%(Y, cat)

    def _repr_(self):
        return "Set of Morphisms from %s to %s in %s"%(
            self.__domain, self.__codomain, self.__category)

    def __call__(self, x, y=None):
        """
        Construct a morphism in this homset from x if possible.
        """
        if x in self:
            return x
        raise TypeError, "Unable to coerce x (=%s) to a morphism in %s"%(x,self)

    def __cmp__(self, other):
        if not isinstance(other, Homset):
            return cmp(type(self), type(other))
        if self.__domain == other.__domain and self.__codomain == other.__codomain \
               and self.__category == other.__category:
            return 0
        return cmp(self.__domain, other.__domain)

    def __contains__(self, x):
        try:
            return x.parent() == self
        except AttributeError:
            pass
        return False

    def natural_map(self):
        return morphism.FormalCoercionMorphism(self)   # good default in many cases

    def identity(self):
        if self.is_endomorphism_set():
            return morphism.IdentityMorphism(self)
        else:
            raise TypeError, "Identity map only defined for endomorphisms. Try natural_map() instead."

    def domain(self):
        return self.__domain

    def codomain(self):
        return self.__codomain

    def is_endomorphism_set(self):
        """
        Return True if the domain and codomain of self are the
        same object.
        """
        return self.__domain is self.__codomain

    def reversed(self):
        """
        Return the corresponding homset, but with the domain and codomain
        reversed.
        """
        return Homset(self.__codomain, self.__domain, self.__category)

class HomsetWithBase(ParentWithBase, Homset):
    def __init__(self, X, Y, cat=None, check=True, base=None):
        Homset.__init__(self, X, Y, cat, check)
        if base is None:
            ParentWithBase.__init__(self, X.base_ring())
        else:
            ParentWithBase.__init__(self, base)

def is_Homset(x):
    """
    Return True if x is a set of homomorphisms in a category.
    """
    return isinstance(x, Homset)

def is_Endset(x):
    """
    Return True if x is a set of endomorphisms in a category.
    """
    return isinstance(x, Homset) and x.is_endomorphism_set()

