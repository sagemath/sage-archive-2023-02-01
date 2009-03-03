"""
Categories

AUTHORS:

- David Kohel and William Stein

Every Sage object lies in a category. Categories in Sage are
modeled on the mathematical idea of category, and are distinct from
Python classes, which are a programming construct.

In most cases, typing ``x.category()`` returns the
category to which `x` belongs. If `C` is a category
and `x` is any object, `C(x)` tries to make an
object in `C` from `x`.

EXAMPLES: We create a couple of categories.

::

    sage: Sets()
    Category of sets
    sage: GSets(AbelianGroup([2,4,9]))
    Category of G-sets for Multiplicative Abelian Group isomorphic to C2 x C4 x C9
    sage: Semigroups()
    Category of semigroups
    sage: VectorSpaces(FiniteField(11))
    Category of vector spaces over Finite Field of size 11
    sage: Ideals(IntegerRing())
    Category of ring ideals in Integer Ring

The default category for elements `x` of an objects
`O` is the category of all objects of `O`. For
example,

::

    sage: V = VectorSpace(RationalField(), 3)
    sage: x = V.gen(1)
    sage: x.category()
    Category of elements of Vector space of dimension 3 over Rational Field
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

__doc_exclude = ['k', 'uniq', 'uniq1', 'latex', \
                 'Category', 'Category_ideal', 'Category_in_ambient', \
                 'Category_module', 'Category_over_base', \
                 'Category_over_base_ring','SageObject']

import weakref
from sage.misc.latex import latex
from sage.structure.sage_object import SageObject


class Category(SageObject):
    """
    The base class for all categories.
    """
    def __init__(self, s=None):
        if s is None:  # figure out from the type name!
            t = str(type(self))
            t = t[t.rfind('.')+1:]
            s = t[:t.rfind("'")]
            self.__label = s
            i = -1
            while i < len(s)-1:
                for i in range(len(s)):
                    if s[i].isupper():
                        s = s[:i] + " " + s[i].lower() + s[i+1:]
                        break
            s = s.lstrip()
        elif isinstance(s, str):
            self.__label = s
        else:
            raise TypeError, "Argument string must be a string."
        self.__category = s

    def __call__(self, x):
        if x in self:
            return x
        return self._call_(x)

    def _call_(self, x):
        raise NotImplementedError

    def _repr_(self):
        return "Category of %s"%self.__category

    def _latex_(self):
        return "\\mathbf{%s}"%self.__label

    def __hash__(self):
        return hash(self.__category)

    def short_name(self):
        return ''.join([x.capitalize() for x in self.__category.split()])

    def __contains__(self, x):
        try:
            c = x.category()
        except AttributeError:
            return False
        return c.is_subcategory(self)

    def is_abelian(self):
        return False

    def is_subcategory(self, c):
        """
        Returns True if self is naturally embedded as a subcategory of c.

        EXAMPLES::

            sage: Rings  = Rings()
            sage: AbGrps = AbelianGroups()
            sage: Rings.is_subcategory(AbGrps)
            True
            sage: AbGrps.is_subcategory(Rings)
            False

        The ``is_subcategory`` function takes into account the
        base.

        ::

            sage: M3 = VectorSpaces(FiniteField(3))
            sage: M9 = VectorSpaces(FiniteField(9, 'a'))
            sage: M3.is_subcategory(M9)
            False
        """
        if not isinstance(c, Category):
            raise TypeError, "Argument c (= %s, type = %s) must be a category"%(c, type(c))
        if self == c: return True
        from category_types import category_hierarchy
        if category_hierarchy.has_key(self.__class__):
            S = category_hierarchy[self.__class__]
            if not c.__class__ in S:
                return False
            return c._parameters().issubset(self._parameters())
        return False

    def _is_subclass(self, c,):
        from category_types import category_hierarchy
        if isinstance(c, Category):
            return self.is_subcategory(c)
        if self.__class__ == c:
            return True
        return c in category_hierarchy[self.__class__]

    def __eq__(self, c):
        if self is c:
            return True
        return False

    def _parameters(self):
        return set([])

    def category(self):
        return Objects()

def is_Category(x):
    """
    Returns True if x is a category.
    """
    return isinstance(x, Category)


#############################################################
# ...
#############################################################
cache = {}
class uniq(object):
    def __new__(cls):
        global cache
        if cache.has_key(cls):
            return cache[cls]
        O = object.__new__(cls)
        cache[cls] = O
        return O


cache1 = {}
class uniq1(object):
    def __new__(cls, arg1):
        global cache1
        key = (cls, arg1)
        if cache1.has_key(key):
            X = cache1[key]()
            if X: return X
        O = object.__new__(cls)
        cache1[key] = weakref.ref(O)
        return O


#############################################################
# Unique category
#############################################################
class Category_uniq(uniq, Category):
    pass

#############################################################
# Sets
#############################################################

class Sets(Category_uniq):
    """
    The category of sets.

    EXAMPLES:
        sage: Sets()
        Category of sets
    """
    def __call__(self, X):
        """
        EXAMPLES:
            sage: Sets()(ZZ)
            Set of elements of Integer Ring
        """
        import sage.sets.all
        return sage.sets.all.Set(X)

    def __reduce__(self):
        """
        For pickling.

        TESTS:
            sage: loads(dumps(Sets())) is Sets()
            True
        """
        return Sets, tuple([])
