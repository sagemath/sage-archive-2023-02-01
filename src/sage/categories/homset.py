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
from sage.structure.parent import Parent, Set_generic

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
    if hasattr(X, '_Hom_'):
        return X._Hom_(Y, cat)
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
            # try hard to find a suitable base category
            subcats_X = category_types.category_hierarchy[cat_X.__class__]
            subcats_Y = category_types.category_hierarchy[cat_Y.__class__]
            cats = set(subcats_X).intersection(set(subcats_Y))
            params = tuple(set(cat_X._parameters()).intersection(cat_Y._parameters()))

            cat = None
            size = -1

            for c in cats:
                if len(category_types.category_hierarchy[c]) > size:
                    try:
                        cat = c(*params)
                        size = len(category_types.category_hierarchy[c])
                    except TypeError:
                        pass

            if cat is None:
                for c in cats:
                    if len(category_types.category_hierarchy[c]) > size:
                        try:
                            cat = c()
                            size = len(category_types.category_hierarchy[c])
                        except TypeError:
                            pass

            if cat is None:
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
        if hasattr(X, '_base') and X._base is not X and X._base is not None:
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
        sage: from sage.categories.homset import is_Endset
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
        self._domain = X
        self._codomain = Y
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
            self._domain, self._codomain, self.__category)

    def __nonzero__(self):
        return True

    def _generic_convert_map(self, S):
        if self._element_constructor is None:
            from sage.categories.morphism import CallMorphism
            from sage.categories.homset import Hom
            return CallMorphism(Hom(S, self))
        else:
            return parent.Parent._generic_convert_map(self, S)

    def homset_category(self):
        """
        Return the category that this is a Hom in, i.e., this is
        typically the category of the domain or codomain object.

        EXAMPLES:
            sage: H = Hom(SymmetricGroup(4), SymmetricGroup(7))
            sage: H.homset_category()
            Category of groups
        """
        return self.__category

    def __call__(self, x, y=None, check=True):
        """
        Construct a morphism in this homset from x if possible.

        EXAMPLES:
            sage: H = Hom(SymmetricGroup(4), SymmetricGroup(7))
            sage: phi = Hom(SymmetricGroup(5), SymmetricGroup(6)).natural_map()
            sage: phi
            Coercion morphism:
              From: SymmetricGroup(5)
              To:   SymmetricGroup(6)
            sage: H(phi)
            Composite map:
              From: SymmetricGroup(4)
              To:   SymmetricGroup(7)
              Defn:   Composite map:
                      From: SymmetricGroup(4)
                      To:   SymmetricGroup(6)
                      Defn:   Call morphism:
                              From: SymmetricGroup(4)
                              To:   SymmetricGroup(5)
                            then
                              Coercion morphism:
                              From: SymmetricGroup(5)
                              To:   SymmetricGroup(6)
                    then
                      Call morphism:
                      From: SymmetricGroup(6)
                      To:   SymmetricGroup(7)

        AUTHOR:
            -- Robert Bradshaw
        """
        if isinstance(x, morphism.Morphism):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                x._set_parent(self) # needed due to non-uniqueness of homsets
                return x
            else:
                if x.domain() != self.domain():
                    mor = x.domain().coerce_map_from(self.domain())
                    if mor is None:
                        raise TypeError, "Incompatible domains: x (=%s) cannot be an element of %s"%(x,self)
                    x = x * mor
                if x.codomain() != self.codomain():
                    mor = self.codomain().coerce_map_from(x.codomain())
                    if mor is None:
                        raise TypeError, "Incompatible codomains: x (=%s) cannot be an element of %s"%(x,self)
                    x = mor * x
                return x
        raise TypeError, "Unable to coerce x (=%s) to a morphism in %s"%(x,self)

    def __cmp__(self, other):
        if not isinstance(other, Homset):
            return cmp(type(self), type(other))
        if self._domain == other._domain:
            if self._codomain == other._codomain:
                if self.__category == other.__category:
                    return 0
                else: return cmp(self.__category, other.__category)
            else: return cmp(self._codomain, other._codomain)
        else: return cmp(self._domain, other._domain)

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
        return self._domain

    def codomain(self):
        return self._codomain

    def is_endomorphism_set(self):
        """
        Return True if the domain and codomain of self are the
        same object.
        """
        return self._domain is self._codomain

    def reversed(self):
        """
        Return the corresponding homset, but with the domain and codomain
        reversed.
        """
        return Homset(self._codomain, self._domain, self.__category)

    ############### For compatibility with old coercion model #######################

    def get_action_c(self, R, op, self_on_left):
        return None

    def coerce_map_from_c(self, R):
        return None

class HomsetWithBase(Homset):
    def __init__(self, X, Y, cat=None, check=True, base=None):
        Homset.__init__(self, X, Y, cat, check)
        if base is None:
            Parent.__init__(self, X.base_ring())
        else:
            Parent.__init__(self, base)

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

