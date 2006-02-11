"""
Generic Set object
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

from sage.structure.element import Element
from sage.structure.all import SageObject
import sage.categories.all
from sage.misc.latex import latex

def Set(X):
    if is_Set(X):
        return X

    if isinstance(X, Element):
        raise TypeError, "Element has no defined underlying set"
    try:
        if isinstance(X, (list, tuple, set)) or X.is_finite():
            return Set_object_enumerated(X)
    except AttributeError:
        pass
    return Set_object(X)

def EnumeratedSet(X):
    return Set_object_enumerated(X)

def is_Set(x):
    return isinstance(x, Set_generic)

class Set_generic(SageObject):
    def category(self):
        return sage.categories.all.Sets()

class Set_object(Set_generic):
    """
    EXAMPLES:
        sage: K = GF(19)
        sage: Set(K)
        {11, 10, 13, 12, 15, 14, 17, 16, 18, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8}
        sage: S = Set(K)

        sage: print latex(S)
        \left\{11, 10, 13, 12, 15, 14, 17, 16, 18, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8\right\}
        sage: loads(S.dumps()) == S
        True

        sage: print latex(Set(ZZ))
        \mbox{\bf{}Z}
    """
    def __init__(self, X):
        self.__object = X

    def _latex_(self):
        return latex(self.__object)

    def _repr_(self):
        return "Set of elements of %s"%self.__object

    def __iter__(self):
        return self.__object.__iter__()

    def __contains__(self, x):
        return x in self.__object

    def __cmp__(self, right):
        if not is_Set(right):
            return -1
        return cmp(self.__object, right.__object)

    def union(self, X):
        if is_Set(X):
            if self == X:
                return self
            return Set_object_union(self.__object, X.__object)
        raise TypeError, "X must be a Set"

    def intersection(self, X):
        if is_Set(X):
            if self == X:
                return self
            return Set_object_intersection(self.__object, X.__object)
        raise TypeError, "X must be a Set"


    def __len__(self):
        try:
            return len(self.__object)
        except TypeError:
            raise NotImplementedError, "computation of order of %s not yet implemented"%self.__object

    def order(self):
        return len(self)

    def object(self):
        return self.__object

class Set_object_enumerated(Set_object):
    def __init__(self, X):
        """
        EXAMPLES:
            sage: S = EnumeratedSet(GF(19)); S
            {11, 10, 13, 12, 15, 14, 17, 16, 18, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8}
            sage: print latex(S)
            \left\{11, 10, 13, 12, 15, 14, 17, 16, 18, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8\right\}
            sage: loads(S.dumps()) == S
            True
        """
        Set_object.__init__(self, X)

    def __iter__(self):
        for x in self.set():
            yield x

    def _latex_(self):
        return '\\left\\{' + ', '.join([latex(x) for x in self.set()])  + '\\right\\}'

    def _repr_(self):
        s = str(self.set())
        return "{" + s[5:-2] + "}"
        #    return "Finite set of elements of %s"%self.__object

    def __contains__(self, x):
        return x in self.set()

    def set(self):
        try:
            return self.__set
        except AttributeError:
            self.__set = set(self.object())
            return self.__set

    def __cmp__(self, other):
        if isinstance(other, Set_object_enumerated):
            if self.set() == other.set():
                return 0
            else:
                return -1
        else:
            return Set_object.__cmp__(self, other)

    def union(self, other):
        if not isinstance(other, Set_object_enumerated):
            return Set_object.union(self, other)
        return Set_object_enumerated(self.set().union(other.set()))

    def intersection(self, other):
        if not isinstance(other, Set_object_enumerated):
            return Set_object.intersection(self, other)
        return Set_object_enumerated(self.set().intersection(other.set()))


class Set_object_union(Set_object):
    def __init__(self, X, Y):
        """
        EXAMPLES:
            sage: S = Set(QQ^2)
            sage: T = Set(ZZ)
            sage: X = S.union(T); X
            Set-theoretic union of Vector space of dimension 2 over Rational Field and Integer Ring

            sage: print latex(X)
            \mbox{\bf{}Q}^{2} \cup \mbox{\bf{}Z}

            sage: loads(X.dumps()) == X
            True
        """
        self.__X = X
        self.__Y = Y
        Set_object.__init__(self, self)

    def __cmp__(self, right):
        """
        TODO: Comparison is basically not implemented!  I don't
        even know how one could implement it.
        """
        if not is_Set(right):
            return -1
        if not isinstance(right, Set_object_union):
            raise NotImplementedError
        if self.__X == right.__X and self.__Y == right.__Y:
            return 0
        raise NotImplementedError

    def _repr_(self):
        return "Set-theoretic union of %s and %s"%(self.__X, self.__Y)

    def _latex_(self):
        return '%s \\cup %s'%(latex(self.__X), latex(self.__Y))

    def __iter__(self):
        for x in self.__X:
            yield x
        for y in self.__Y:
            yield y

    def __contains__(self, x):
        return x in self.__X or x in self.__Y

    def order(self):
        return self.__X.order() + self.__Y.order()

class Set_object_intersection(Set_object):
    def __init__(self, X, Y):
        """
        EXAMPLES:
            sage: S = Set(QQ^2)
            sage: T = Set(ZZ)
            sage: X = S.intersection(T); X
            Set-theoretic intersection of Vector space of dimension 2 over Rational Field and Integer Ring
            sage: print latex(X)
            \mbox{\bf{}Q}^{2} \cap \mbox{\bf{}Z}

            sage: loads(X.dumps()) == X
            True
        """
        self.__X = X
        self.__Y = Y
        Set_object.__init__(self, self)


    def __cmp__(self, right):
        if not is_Set(right):
            return -1
        if not isinstance(right, Set_object_intersection):
            raise NotImplementedError
        if self.__X == right.__X and self.__Y == right.__Y:
            return 0
        raise NotImplementedError

    def _repr_(self):
        return "Set-theoretic intersection of %s and %s"%(self.__X, self.__Y)

    def _latex_(self):
        return '%s \\cap %s'%(latex(self.__X), latex(self.__Y))

    def __iter__(self):
        for x in self.__X:
            if x in self.__Y:
                yield x

    def __contains__(self, x):
        return x in self.__X and x in self.__Y

    def order(self):
        raise NotImplementedError

