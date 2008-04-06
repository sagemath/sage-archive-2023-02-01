r"""
Combinatorial Algebras

A combinatorial algebra is an algebra whose basis elements are indexed
by a combinatorial class.  Some examples of combinatorial algebras are
the symmetric group algebra of order n (indexed by permutations of size n)
and the algebra of symmetric functions (indexed by integer partitions).

The CombinatorialAlgebra base class makes it easy to define and work
with new combinatorial algebras in Sage.  For example, the following
code constructs an algebra which models the power-sum symmetric
functions.

sage: class PowerSums(CombinatorialAlgebra):
...     def __init__(self, R):
...         self._combinatorial_class = Partitions()
...         self._one = Partition([])
...         self._name = 'Power-sum symmetric functions'
...         self._prefix = 'p'
...         CombinatorialAlgebra.__init__(self, R)
...     def _multiply_basis(self, a, b):
...         l = list(a)+list(b)
...         l.sort(reverse=True)
...         return Partition(l)
...

sage: ps = PowerSums(QQ); ps
Power-sum symmetric functions over Rational Field
sage: ps([2,1])^2
p[2, 2, 1, 1]
sage: ps([2,1])+2*ps([1,1,1])
2*p[1, 1, 1] + p[2, 1]
sage: ps(2)
2*p[]

The important things to define are ._combinatorial_class which specifies
the combinatorial class that indexes the basis elements, ._one which
specifies the identity element in the algebra, ._name which specifies
the name of the algebra, ._prefix which is the string put in front of
each basis element, and finally a _multiply or _multiply basis method
that defines the multiplication in the algebra.

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

from sage.rings.all import Ring, Integer
from sage.misc.misc import repr_lincomb
from sage.algebras.algebra import Algebra
from sage.algebras.algebra_element import AlgebraElement
import sage.structure.parent_base
import sage.combinat.partition
from sage.modules.free_module_element import vector
from sage.matrix.all import MatrixSpace

class CombinatorialAlgebraElement(AlgebraElement):
    def __init__(self, A, x):
        """
        Create a combinatorial algebra element x.  This should never
        be called directly, but only through the parent combinatorial
        algebra's __call__ method.

        TESTS:
            sage: s = SFASchur(QQ)
            sage: a = s._element_class(s, {Partition([2,1]):QQ(2)}); a
            2*s[2, 1]
            sage: a == loads(dumps(a))
            True
        """
        AlgebraElement.__init__(self, A)
        self._monomial_coefficients = x

    def __iter__(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: [i for i in sorted(a)]
            [([2, 1], 1), ([3], 1)]

        """
        return self._monomial_coefficients.iteritems()

    def __contains__(self, x):
        """
        Returns whether or not a combinatorial object x indexing a basis
        element is in the support of self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: Partition([2,1]) in a
            True
            sage: Partition([1,1,1]) in a
            False
        """
        return x in self._monomial_coefficients and self._monomial_coefficients[x] != 0

    def monomial_coefficients(self):
        """
        Return the internal dictionary which has the combinatorial
        objects indexing the basis as keys and their corresponding
        coefficients as values.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: d = a.monomial_coefficients()
            sage: type(d)
            <type 'dict'>
            sage: d[ Partition([2,1]) ]
            1
            sage: d[ Partition([3,2]) ]
            2
        """
        return self._monomial_coefficients

    def __repr__(self):
        """
        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: a = 2 + QS3([2,1,3])
            sage: print a.__repr__()
            2*[1, 2, 3] + [2, 1, 3]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        prefix = self.parent().prefix()
        mons = [ prefix + repr(m) for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def _latex_(self):
        """
        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: a = 2 + QS3([2,1,3])
            sage: latex(a) #indirect doctest
            2[1,2,3] + [2,1,3]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        prefix = self.parent().prefix()
        if prefix == "":
            mons = [ prefix + '[' + ",".join(map(str, m)) + ']' for (m, _) in v ]
        else:
            mons = [ prefix + '_{' + ",".join(map(str, m)) + '}' for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs, is_latex=True).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def __cmp__(left, right):
        """
        The ordering is the one on the underlying sorted list of
        (monomial,coefficients) pairs.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: cmp(a,b) #indirect doctest
            1
        """
        nonzero = lambda mc: mc[1] != 0
        v = filter(nonzero, left._monomial_coefficients.items())
        v.sort()
        w = filter(nonzero, right._monomial_coefficients.items())
        w.sort()
        return cmp(v, w)

    def _add_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]) + s([5,4]) # indirect doctest
            s[2, 1] + s[5, 4]
            sage: a = s([2,1]) + 0
            sage: len(a.monomial_coefficients())
            1
        """
        A = self.parent()
        BR = A.base_ring()
        z_elt = dict(self._monomial_coefficients)
        for m, c in y._monomial_coefficients.iteritems():
            if z_elt.has_key(m):
                cm = z_elt[m] + c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = c


        #Remove all entries that are equal to 0
        del_list = []
        zero = BR(0)
        for m, c in z_elt.iteritems():
            if c == zero:
                del_list.append(m)
        for m in del_list:
            del z_elt[m]

        return A._from_dict(z_elt)


    def _neg_(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: -s([2,1]) # indirect doctest
            -s[2, 1]
        """
        return self.map_coefficients(lambda c: -c)


    def _sub_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]) - s([5,4]) # indirect doctest
            s[2, 1] - s[5, 4]
        """
        A = self.parent()
        BR = A.base_ring()
        z_elt = dict(self._monomial_coefficients)
        for m, c in y._monomial_coefficients.iteritems():
            if z_elt.has_key(m):
                cm = z_elt[m] - c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = -c

        #Remove all entries that are equal to 0
        zero = BR(0)
        del_list = []
        for m, c in z_elt.iteritems():
            if c == zero:
                del_list.append(m)
        for m in del_list:
            del z_elt[m]

        return A._from_dict(z_elt)

    def _mul_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2])
            sage: a._mul_(a) #indirect doctest
            s[2, 2] + s[3, 1] + s[4]
        """
        return self.parent().multiply(self, y)

    def _div_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2])
            sage: a._div_(s(2))
            1/2*s[2]
            sage: a._div_(a)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert self (= s[2])

        """
        return self.parent().multiply(self, ~y)

    def __invert__(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: ~s(2)
            1/2*s[]
            sage: ~s([2,1])
            Traceback (most recent call last):
            ...
            ValueError: cannot invert self (= s[2, 1])
        """
        mcs = self._monomial_coefficients
        one = self.parent()._one
        if len(mcs) == 1 and one in mcs:
            return self.parent()( ~mcs[ one ] )
        else:
            raise ValueError, "cannot invert self (= %s)"%self

    def __pow__(self, n):
        """
        Returns self to the $n^{th}$ power.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2])^2
            s[2, 2] + s[3, 1] + s[4]


        TESTS:
            sage: s = SFASchur(QQ)
            sage: z = s([2,1])
            sage: z
            s[2, 1]
            sage: z^2
            s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

            sage: e = SFAElementary(QQ)
            sage: y = e([1])
            sage: y^2
            e[1, 1]
            sage: y^4
            e[1, 1, 1, 1]
        """
        if not isinstance(n, (int, Integer)):
            raise TypeError, "n must be an integer"
        A = self.parent()
        z = A(Integer(1))
        for _ in range(n):
            z *= self
        return z

    def _coefficient_fast(self, m, default=None):
        """
        Returns the coefficient of m in self, where m is key
        in self._monomial_coefficients.

        EXAMPLES:
            sage: p = Partition([2,1])
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: a._coefficient_fast([2,1])
            Traceback (most recent call last):
            ...
            TypeError: list objects are unhashable
            sage: a._coefficient_fast(p)
            1
            sage: a._coefficient_fast(p, 2)
            1
        """
        if default is None:
            default = self.base_ring()(0)
        return self._monomial_coefficients.get(m, default)

    def coefficient(self, m):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) - 2*s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficient([4])
            1
            sage: z.coefficient([2,1])
            -2
        """
        p = self.parent()
        if isinstance(m, p._combinatorial_class.object_class):
            return self._monomial_coefficients.get(m, p.base_ring().zero_element())
        if m in p._combinatorial_class:
            return self._monomial_coefficients.get(p._combinatorial_class.object_class(m), p.base_ring().zero_element())
        else:
            raise TypeError, "you must specify an element of %s"%p._combinatorial_class


    def is_zero(self):
        """
        Returns True if and only self == 0.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]).is_zero()
            False
            sage: s(0).is_zero()
            True
            sage: (s([2,1]) - s([2,1])).is_zero()
            True
        """
        BR = self.parent().base_ring()
        for v in self._monomial_coefficients.values():
            if v != BR(0):
                return False
        return True

    def __len__(self):
        """
        Returns the number of basis elements of self with
        nonzero coefficients.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: len(z)
            4
        """
        return self.length()

    def length(self):
        """
        Returns the number of basis elements of self with
        nonzero coefficients.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.length()
            4
        """
        return len( filter(lambda x: self._monomial_coefficients[x] != 0, self._monomial_coefficients) )

    def support(self):
        """
        Returns a pair [mons, cffs] of lists of the monomials
        of self (mons) and their respective coefficients (cffs).

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.support()
            [[[1], [1, 1, 1], [2, 1], [4]], [1, 1, 1, 1]]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        mons = [ m for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        return [mons, cffs]

    def monomials(self):
        """
        Returns a list of the combinatorial objects indexing
        the basis elements of self which non-zero coefficients.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.monomials()
            [[1], [1, 1, 1], [2, 1], [4]]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        mons = [ m for (m, _) in v ]
        return mons

    def coefficients(self):
        """
        Returns a list of the coefficents appearing on the
        basiselements in self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficients()
            [1, 1, 1, 1]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        cffs = [ c for (_, c) in v ]
        return cffs

    def _vector_(self, new_BR=None):
        """
        Returns a vector version of self. If new_BR is specified,
        then in returns a vector over new_BR.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a._vector_()
            (2, 0, 0, 0, 0, 4)
            sage: vector(a)
            (2, 0, 0, 0, 0, 4)
            sage: a._vector_(RDF)
            (2.0, 0.0, 0.0, 0.0, 0.0, 4.0)
        """
        parent = self.parent()
        if parent.get_order() is None:
            cc = parent._combinatorial_class
        else:
            cc = parent.get_order()

        if new_BR is None:
            new_BR = parent.base_ring()

        return vector(new_BR, [new_BR(self._monomial_coefficients.get(m, 0)) for m in cc])

    def to_vector(self):
        """
        Returns a vector version of self.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a.to_vector()
            (2, 0, 0, 0, 0, 4)
        """
        return self._vector_()

    def _matrix_(self, new_BR = None):
        """
        Returns a matrix version of self obtained by
        the action of self on the left.  If new_BR
        is specified, then the matrix will be over
        new_BR.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = QS3([2,1,3])
            sage: a._matrix_()
            [0 0 1 0 0 0]
            [0 0 0 0 1 0]
            [1 0 0 0 0 0]
            [0 0 0 0 0 1]
            [0 1 0 0 0 0]
            [0 0 0 1 0 0]
            sage: a._matrix_(RDF)
            [0.0 0.0 1.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 1.0 0.0]
            [1.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 1.0]
            [0.0 1.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 1.0 0.0 0.0]

        """
        parent = self.parent()

        if parent.get_order() is None:
            cc = parent._combinatorial_class
        else:
            cc = parent.get_order()

        BR = parent.base_ring()
        if new_BR is None:
            new_BR = BR

        MS = MatrixSpace(new_BR, parent._dim, parent._dim)
        l = [ (self*parent(m)).to_vector() for m in cc ]
        return MS(l).transpose()

    def to_matrix(self):
        """
        Returns a matrix version of self obtained by
        the action of self on the left.  If new_BR
        is specified, then the matrix will be over
        new_BR.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = QS3([2,1,3])
            sage: a._matrix_() # indirect doctest
            [0 0 1 0 0 0]
            [0 0 0 0 1 0]
            [1 0 0 0 0 0]
            [0 0 0 0 0 1]
            [0 1 0 0 0 0]
            [0 0 0 1 0 0]
            sage: a._matrix_(RDF)
            [0.0 0.0 1.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 1.0 0.0]
            [1.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 1.0]
            [0.0 1.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 1.0 0.0 0.0]

        """
        return self._matrix_()

    def map_coefficients(self, f):
        """
        Returns a new element of self.parent() obtained
        by applying the function f to all of the coefficients
        of self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: a.map_coefficients(lambda x: x*2)
            2*s[2, 1] + 4*s[3, 2]
        """
        res = self.parent()(0)
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            z_elt[m] = f(c)
        res._monomial_coefficients = z_elt
        return res

    def map_basis(self, f):
        """
        Returns a new element of self.parent() obtained
        by applying the function f to all of the combinatorial
        objects indexing the basis elements.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: a.map_basis(lambda x: x.conjugate())
            s[2, 1] + 2*s[2, 2, 1]
        """
        res = self.parent()(0)
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            z_elt[f(m)] = c
        res._monomial_coefficients = z_elt
        return res

    def map_mc(self, f):
        """
        Returns a new element of self.parent() obtained
        by applying the function f to a monomial coefficient
        (m,c) pair.  f returns a (new_m, new_c) pair.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: f = lambda m,c: (m.conjugate(), 2*c)
            sage: a = s([2,1]) + s([1,1,1])
            sage: a.map_mc(f)
            2*s[2, 1] + 2*s[3]
        """
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            new_m, new_c = f(m,c)
            z_elt[new_m] = new_c
        return self.parent()._from_dict(z_elt)

class CombinatorialAlgebra(Algebra):
    def __init__(self, R, element_class=None):
        """
        TESTS:
            sage: s = SFASchur(QQ)
            sage: s == loads(dumps(s))
            True
        """
        #Make sure R is a ring with unit element
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        try:
            z = R(Integer(1))
        except:
            raise ValueError, "R must have a unit element"

        #Check to make sure that the user defines the necessary
        #attributes / methods to make the combinatorial algebra
        #work
        required = ['_combinatorial_class','_one',]
        for r in required:
            if not hasattr(self, r):
                raise ValueError, "%s is required"%r
        if not hasattr(self, '_multiply') and not hasattr(self, '_multiply_basis'):
            raise ValueError, "either _multiply or _multiply_basis is required"

        #Create a class for the elements of this combinatorial algebra
        #We need to do this so to distinguish between element of different
        #combinatorial algebras
        if element_class is None:
            if not hasattr(self, '_element_class'):
                class CAElement(CombinatorialAlgebraElement):
                    pass
                self._element_class = CAElement
        else:
            self._element_class = element_class

        #Set the dimension
        try:
            self._dim = self._combinatorial_class.count()
        except (ValueError, NotImplementedError):
            self._dim = None

        self._order = None

        #Initialize the base structure
        sage.structure.parent_base.ParentWithBase.__init__(self, R)

    _prefix = ""
    _name   = "CombinatorialAlgebra -- change me"

    def basis(self):
        """
        Returns a list of the basis elements of self.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3.basis()
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

        """
        return [self(x) for x in self._combinatorial_class]

    def __call__(self, x):
        """
        Coerce x into self.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3(2)
            2*[1, 2, 3]
            sage: QS3([2,3,1])
            [2, 3, 1]
        """
        R = self.base_ring()
        eclass = self._element_class

        #Coerce ints to Integers
        if isinstance(x, int):
            x = Integer(x)


        if hasattr(self, '_coerce_start'):
            try:
                return self._coerce_start(x)
            except TypeError:
                pass

        #x is an element of the same type of combinatorial algebra
        if hasattr(x, 'parent') and x.parent().__class__ is self.__class__:
            P = x.parent()
            #same base ring
            if P is self:
                return x
            #different base ring -- coerce the coefficients from into R
            else:
                return eclass(self, dict([ (e1,R(e2)) for e1,e2 in x._monomial_coefficients.items()]))
        #x is an element of the basis combinatorial class
        elif isinstance(x, self._combinatorial_class.object_class):
            return eclass(self, {x:R(1)})
        elif x in self._combinatorial_class:
            return eclass(self, {self._combinatorial_class.object_class(x):R(1)})
        #Coerce elements of the base ring
        elif hasattr(x, 'parent') and x.parent() is R:
            if x == R(0):
                return eclass(self, {})
            else:
                return eclass(self, {self._one:x})
        #Coerce things that coerce into the base ring
        elif R.has_coerce_map_from(x.parent()):
            rx = R(x)
            if rx == R(0):
                return eclass(self, {})
            else:
                return eclass(self, {self._one:R(x)})
        else:
            if hasattr(self, '_coerce_end'):
                try:
                    return self._coerce_start(x)
                except TypeError:
                    pass
            raise TypeError, "do not know how to make x (= %s) an element of self (=%s)"%(x,self)


    def _an_element_impl(self):
        """
        Returns an element of self, namely the unit element.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s._an_element_impl()
            s[]
            sage: _.parent() is s
            True
        """
        return self._element_class(self, {self._one:self.base_ring()(1)})

    def __repr__(self):
        """
        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: print QS3.__repr__()
            Symmetric group algebra of order 3 over Rational Field

        """
        return self._name + " over %s"%self.base_ring()

    def combinatorial_class(self):
        """
        Returns the combinatorial class that indexes the basis
        elements.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.combinatorial_class()
            Partitions
        """
        return self._combinatorial_class

    def _coerce_impl(self, x):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s._coerce_impl(2)
            2*s[]
        """
        try:
            R = x.parent()
            if R.__class__ is self.__class__:
                #Only perform the coercion if we can go from the base
                #ring of x to the base ring of self
                if self.base_ring().has_coerce_map_from( R.base_ring() ):
                    return self(x)
        except AttributeError:
            pass

        # any ring that coerces to the base ring
        return self._coerce_try(x, [self.base_ring()])

    def dimension(self):
        """
        Returns the dimension of the combinatorial algebra (which is given
        by the number of elements in the associated combinatorial class).

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.dimension()
            +Infinity
        """
        return self._combinatorial_class.count()

    def set_order(self, order):
        """
        Sets the order of the elements of the combinatorial class.

        If .set_order() has not been called, then the ordering is
        the one used in the generation of the elements of self's
        associated combinatorial class.

        EXAMPLES:
            sage: QS2 = SymmetricGroupAlgebra(QQ,2)
            sage: b = QS2.basis()
            sage: b.reverse()
            sage: QS2.set_order(b)
            sage: QS2.get_order()
            [[2, 1], [1, 2]]

        """
        self._order = order

    def get_order(self):
        """
        Returns the order of the elements in the basis.

        EXAMPLES:
            sage: QS2 = SymmetricGroupAlgebra(QQ,2)
            sage: QS2.get_order()
            [[1, 2], [2, 1]]
        """
        if self._order is None:
            self._order = self.combinatorial_class().list()
        return self._order

    def prefix(self):
        """
        Returns the prefix used when displaying elements of self.

        EXAMPLES:
            sage: X = SchubertPolynomialRing(QQ)
            sage: X.prefix()
            'X'
        """
        return self._prefix

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: XQ = SchubertPolynomialRing(QQ)
            sage: XZ = SchubertPolynomialRing(ZZ)
            sage: XQ == XZ #indirect doctest
            False
            sage: XQ == XQ
            True
        """
        if not isinstance(other, self.__class__):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        return 0

    def _apply_module_morphism(self, x, f):
        """
        Returns the image of x under the module morphism defined by
        extending f by linearity.

        INPUT:
            -- x : a element of self
            -- f : a function that takes in a combinatorial object
                   indexing a basis element and returns an element
                   of the target domain

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([3]) + s([2,1]) + s([1,1,1])
            sage: b = 2*a
            sage: f = lambda part: len(part)
            sage: s._apply_module_morphism(a, f) #1+2+3
            6
            sage: s._apply_module_morphism(b, f) #2*(1+2+3)
            12


        """
        res = 0
        for m, c in x._monomial_coefficients.iteritems():
            res += c*f(m)
        return res


    def _apply_module_endomorphism(self, a, f):
        """
        This takes in a function from the basis elements
        to the elements of self and applies it linearly
        to a. Note that _apply_module_endomorphism does not
        require multiplication on self to be defined.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: f = lambda part: 2*s(part.conjugate())
            sage: s._apply_module_endomorphism( s([2,1]) + s([1,1,1]), f)
            2*s[2, 1] + 2*s[3]

        """
        mcs = a.monomial_coefficients()
        base_ring = self.base_ring()
        zero = base_ring(0)

        z_elt = {}
        for basis_element in mcs:
            f_mcs = f(basis_element).monomial_coefficients()
            for f_basis_element in f_mcs:
                z_elt[ f_basis_element ] = z_elt.get(f_basis_element, zero) + mcs[basis_element]*f_mcs[f_basis_element]

        return self._from_dict(z_elt)


    def multiply(self,left,right):
        """
        Returns left*right where left and right are elements of self.
        multiply() uses either _multiply or _multiply basis to carry
        out the actual multiplication.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2])
            sage: s.multiply(a,a)
            s[2, 2] + s[3, 1] + s[4]

        """
        A = left.parent()
        BR = A.base_ring()
        z_elt = {}

        #Do the case where the user specifies how to multiply basis
        #elements
        if hasattr(self, '_multiply_basis'):
            for (left_m, left_c) in left._monomial_coefficients.iteritems():
                for (right_m, right_c) in right._monomial_coefficients.iteritems():
                    res = self._multiply_basis(left_m, right_m)
                    #Handle the case where the user returns a dictionary
                    #where the keys are the monomials and the values are
                    #the coefficients.  If res is not a dictionary, then
                    #it is assumed to be an element of self
                    if not isinstance(res, dict):
                        if isinstance(res, self._element_class):
                            res = res._monomial_coefficients
                        else:
                            res = {res: BR(1)}
                    for m in res:
                        if m in z_elt:
                            z_elt[ m ] = z_elt[m] + left_c * right_c * res[m]
                        else:
                            z_elt[ m ] = left_c * right_c * res[m]

        #We assume that the user handles the multiplication correctly on
        #his or her own, and returns a dict with monomials as keys and
        #coefficients as values
        else:
            m = self._multiply(left, right)
            if isinstance(m, self._element_class):
                return m
            if not isinstance(m, dict):
                z_elt = m.monomial_coefficients()
            else:
                z_elt = m

        #Remove all entries that are equal to 0
        BR = self.base_ring()
        zero = BR(0)
        del_list = []
        for m, c in z_elt.iteritems():
            if c == zero:
                del_list.append(m)
        for m in del_list:
            del z_elt[m]

        return self._from_dict(z_elt)

    def _from_dict(self, d):
        """
        Given a monomial coefficient dictionary d, return the element
        of self with the dictionary.

        EXAMPLES:
            sage: e = SFAElementary(QQ)
            sage: s = SFASchur(QQ)
            sage: a = e([2,1]) + e([1,1,1]); a
            e[1, 1, 1] + e[2, 1]
            sage: s._from_dict(a.monomial_coefficients())
            s[1, 1, 1] + s[2, 1]

        """
        return self._element_class(self, d)
