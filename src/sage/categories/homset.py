"""
Homsets

AUTHORS:

- David Kohel and William Stein

- David Joyner (2005-12-17): added examples

- William Stein (2006-01-14): Changed from Homspace to Homset.

- Nicolas M. Thiery (2008-12-): Updated for the new category framework
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

from sage.categories.category import Category
import morphism
from sage.structure.parent import Parent, Set_generic
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_function
import types

_cache = {}
def Hom(X, Y, category=None):
    """
    Create the space of homomorphisms from X to Y in the category category.

    INPUT:


    -  ``X`` - anything

    -  ``Y`` - anything

    -  ``category`` - (optional) category in which the morphisms
       must be


    OUTPUT: a homset in category

    EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: Hom(V, V)
        Set of Morphisms from Vector space of dimension 3 over Rational
        Field to Vector space of dimension 3 over Rational Field in
        Category of vector spaces over Rational Field
        sage: G = AlternatingGroup(3)
        sage: Hom(G, G)
        Set of Morphisms from AlternatingGroup(3) to AlternatingGroup(3) in Category of finite permutation groups
        sage: Hom(ZZ, QQ, Sets())
        Set of Morphisms from Integer Ring to Rational Field in Category of sets

        sage: Hom(FreeModule(ZZ,1), FreeModule(QQ,1))
        Set of Morphisms from Ambient free module of rank 1 over the principal ideal domain Integer Ring to Vector space of dimension 1 over Rational Field in Category of modules with basis over Integer Ring
        sage: Hom(FreeModule(QQ,1), FreeModule(ZZ,1))
        Set of Morphisms from Vector space of dimension 1 over Rational Field to Ambient free module of rank 1 over the principal ideal domain Integer Ring in Category of vector spaces over Rational Field


    TODO: design decision: how much of the homset comes from the
    category of X and Y, and how much from the specific X and Y.  In
    particular, do we need several parent classes depending on X and
    Y, or does the difference only lie in the elements (i.e. the
    morphism), and of course how the parent calls their constructors.

    """


    # This should use cache_function instead
    # However it breaks somehow the coercion (see e.g. sage -t sage.rings.real_mpfr)
    # To be investigated.
    global _cache
    key = (X,Y,category)
    if _cache.has_key(key):
        H = _cache[key]()
        # What is this test for? Why does the cache ever contain a 0 value?
        # This actually occurs: see e.g. sage -t  sage-combinat/sage/categories/modules_with_basis.py
        if H:
            return H

    try:
        return X._Hom_(Y, category)
    except (AttributeError, TypeError):
        pass

#    if category is None or (category is X.category() and category is Y.category()):
#        try:
#            H = X._Hom_(Y)
#        except AttributeError:
#            pass
#    import category_types

    cat_X = X.category()
    cat_Y = Y.category()
    if category is None:
        if cat_X.is_subcategory(cat_Y):
            category = cat_Y
        elif cat_Y.is_subcategory(cat_X):
            # NT: this "category is None" test is useless and could be removed
            # NT: why is there an assymmetry between X and Y?
            if not (category is None) and not (cat_X is cat_Y):
                raise ValueError, "No unambiguous category found for Hom from %s to %s."%(X,Y)
            category = cat_X
        else:
            # Search for the lowest common super category
            subcats_X = cat_X.all_super_categories(proper = True)
            subcats_Y = set(cat_Y.all_super_categories(proper = True))
            category = None
            for c in subcats_X:
                if c in subcats_Y:
                    category = c
                    break

            if category is None:
                raise TypeError, "No suitable category found for Hom from %s to %s."%(X,Y)

    elif isinstance(category, Category):
        if not isinstance(category, Category):
            raise TypeError, "Argument category (= %s) must be a category."%category
        if not cat_X.is_subcategory(category) \
               or not cat_Y.is_subcategory(category):
            raise TypeError, \
                  "Argument category (= %s) is incompatible with %s: %s and %s: %s."%(category, X, cat_X, Y, cat_Y)
    else:
        raise TypeError, "Argument category (= %s) must be a category."%category


    # coercing would be incredibly annoying, since the domain and codomain
    # are totally different objects
    #X = category(X); Y = category(Y)

    # construct H
    # Design question: should the Homset classes get the category or the homset category?
    # For the moment, this is the category, for compatibility with the current implementations
    # of Homset in rings, schemes, ...
    H = category.hom_category().parent_class(X, Y, category = category)

    ##_cache[key] = weakref.ref(H)
    _cache[(X, Y, category)] = weakref.ref(H)

    return H

def hom(X, Y, f):
    """
    Return Hom(X,Y)(f), where f is data that defines an element of
    Hom(X,Y).

    EXAMPLES::

        sage: R, x = PolynomialRing(QQ,'x').objgen()
        sage: phi = hom(R, QQ, [2])
        sage: phi(x^2 + 3)
        7
    """
    return Hom(X,Y)(f)

def End(X, category=None):
    r"""
    Create the set of endomorphisms of X in the category category.

    INPUT:


    -  ``X`` - anything

    -  ``category`` - (optional) category in which to coerce X


    OUTPUT: a set of endomorphisms in category

    EXAMPLES::

        sage: V = VectorSpace(QQ, 3)
        sage: End(V)
        Set of Morphisms from Vector space of dimension 3 over Rational
        Field to Vector space of dimension 3 over Rational Field in
        Category of vector spaces over Rational Field

    ::

        sage: G = AlternatingGroup(3)
        sage: S = End(G); S
        Set of Morphisms from AlternatingGroup(3) to AlternatingGroup(3) in Category of finite permutation groups
        sage: from sage.categories.homset import is_Endset
        sage: is_Endset(S)
        True
        sage: S.domain()
        Alternating group of order 3!/2 as a permutation group

    To avoid creating superfluous categories, homsets are in the
    homset category of the lowest category which currently says
    something specific about its homsets. For example, S is not
    in the category of hom sets of the category of groups::

        sage: S.category()
        Category of hom sets in Category of sets
        sage: End(QQ).category()
        Category of hom sets in Category of rings
    """
    return Hom(X,X, category)

def end(X, f):
    """
    Return End(X)(f), where f is data that defines an element of
    End(X).

    EXAMPLES::

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

    EXAMPLES::

        sage: H = Hom(QQ^2, QQ^3)
        sage: loads(H.dumps()) == H
        True
        sage: E = End(AffineSpace(2, names='x,y'))
        sage: loads(E.dumps()) == E
        True
    """
    def __init__(self, X, Y, category=None, base = None, check=True):
        """
        TESTS::
            sage: X = ZZ['x']; X.rename("X")
            sage: Y = ZZ['y']; Y.rename("Y")
            sage: class MyHomset(Homset):
            ...       def my_function(self, x):
            ...           return Y(x[0])
            ...       def _an_element_(self):
            ...           return sage.categories.morphism.SetMorphism(self, self.my_function)
            ...
            sage: import __main__; __main__.MyHomset = MyHomset # fakes MyHomset being defined in a Python module
            sage: H = MyHomset(X, Y, category=Monoids(), base = ZZ)
            sage: H
            Set of Morphisms from X to Y in Category of monoids
            sage: TestSuite(H).run()

            sage: H = MyHomset(X, Y, category=1, base = ZZ)
            Traceback (most recent call last):
            ...
            TypeError: category (=1) must be a category

            sage: H
            Set of Morphisms from X to Y in Category of monoids
            sage: TestSuite(H).run()
            sage: H = MyHomset(X, Y, category=1, base = ZZ, check = False)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.integer.Integer' object has no attribute 'hom_category'
        """

        self._domain = X
        self._codomain = Y
        if category is None:
            category = X.category()
        self.__category = category
        if check:
            if not isinstance(category, Category):
                raise TypeError, "category (=%s) must be a category"%category
            #if not X in category:
            #    raise TypeError, "X (=%s) must be in category (=%s)"%(X, category)
            #if not Y in category:
            #    raise TypeError, "Y (=%s) must be in category (=%s)"%(Y, category)

        Parent.__init__(self, base = base, category = category.hom_category())

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
            return Parent._generic_convert_map(self, S)

    def homset_category(self):
        """
        Return the category that this is a Hom in, i.e., this is typically
        the category of the domain or codomain object.

        EXAMPLES::

            sage: H = Hom(AlternatingGroup(4), AlternatingGroup(7))
            sage: H.homset_category()
            Category of finite permutation groups
        """
        return self.__category

    def __call__(self, x=None, y=None, check=True, on_basis=None):
        """
        Construct a morphism in this homset from x if possible.

        EXAMPLES::

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

            sage: H = Hom(ZZ, ZZ, Sets())
            sage: f = H( lambda x: x + 1 )
            sage: f.parent()
            Set of Morphisms from Integer Ring to Integer Ring in Category of sets
            sage: f.domain()
            Integer Ring
            sage: f.codomain()
            Integer Ring
            sage: f(1), f(2), f(3)
            (2, 3, 4)

            sage: H = Hom(Set([1,2,3]), Set([1,2,3]))
            sage: f = H( lambda x: 4-x )
            sage: f.parent()
            Set of Morphisms from {1, 2, 3} to {1, 2, 3} in Category of sets
            sage: f(1), f(2), f(3) # todo: not implemented


        - Robert Bradshaw, with changes by Nicolas M. Thiery
        """
        # Temporary workaround: currently, HomCategory.ParentMethods's cannot override
        # this __call__ method because of the class inheritance order
        # This dispatches back the call there
        if on_basis is not None:
            return self.__call_on_basis__(on_basis = on_basis)
        assert x is not None

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

        if isinstance(x, types.FunctionType) or isinstance(x, types.MethodType):
            return self.element_class_set_morphism(self, x)

        raise TypeError, "Unable to coerce x (=%s) to a morphism in %s"%(x,self)

    @lazy_attribute
    def element_class_set_morphism(self):
        """
        A base class for elements of this homset which are
        also ``SetMorphism``, i.e. implemented by mean of a
        Python function.

        This is currently plain ``SetMorphism``, without inheritance
        from categories.

        Todo: refactor during the upcoming homset cleanup.

        EXAMPLES::

            sage: H = Hom(ZZ, ZZ)
            sage: H.element_class_set_morphism
            <type 'sage.categories.morphism.SetMorphism'>
        """
        return self.__make_element_class__(morphism.SetMorphism)

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
        Return True if the domain and codomain of self are the same
        object.
        """
        return self._domain is self._codomain

    def reversed(self):
        """
        Return the corresponding homset, but with the domain and codomain
        reversed.

        EXAMPLES::

            sage: H = Hom(ZZ^2, ZZ^3); H
            Set of Morphisms from Ambient free module of rank 2 over the principal ideal domain Integer Ring to Ambient free module of rank 3 over the principal ideal domain Integer Ring in Category of modules with basis over Integer Ring
            sage: type(H)
            <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
            sage: H.reversed()
            Set of Morphisms from Ambient free module of rank 3 over the principal ideal domain Integer Ring to Ambient free module of rank 2 over the principal ideal domain Integer Ring in Category of hom sets in Category of modules with basis over Integer Ring
            sage: type(H.reversed())
            <class 'sage.modules.free_module_homspace.FreeModuleHomspace_with_category'>
        """
        return Hom(self.codomain(), self.domain(), category = self.category())

    ############### For compatibility with old coercion model #######################

    def get_action_c(self, R, op, self_on_left):
        return None

    def coerce_map_from_c(self, R):
        return None

# Really needed???
class HomsetWithBase(Homset):
    def __init__(self, X, Y, category=None, check=True, base=None):
        """
        TESTS::
            sage: X = ZZ['x']; X.rename("X")
            sage: Y = ZZ['y']; Y.rename("Y")
            sage: class MyHomset(HomsetWithBase):
            ...       def my_function(self, x):
            ...           return Y(x[0])
            ...       def _an_element_(self):
            ...           return sage.categories.morphism.SetMorphism(self, self.my_function)
            ...
            sage: import __main__; __main__.MyHomset = MyHomset # fakes MyHomset being defined in a Python module
            sage: H = MyHomset(X, Y, category=Monoids())
            sage: H
            Set of Morphisms from X to Y in Category of monoids
            sage: H.base()
            Integer Ring
            sage: TestSuite(H).run()
        """
        if base is None:
            base = X.base_ring()
        Homset.__init__(self, X, Y, check=check, category=category, base = base)

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

