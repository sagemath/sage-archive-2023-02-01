r"""
Combinatorial Algebras
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

from sage.interfaces.all import gap, maxima
from sage.rings.all import QQ, RR, ZZ
from sage.rings.arith import binomial
from sage.misc.sage_eval import sage_eval
from sage.libs.all import pari
from sage.rings.all import Ring, factorial, Integer
from random import randint
from sage.misc.misc import prod, repr_lincomb
from sage.structure.sage_object import SageObject
import __builtin__
from sage.algebras.algebra import Algebra
from sage.algebras.algebra_element import AlgebraElement
import sage.structure.parent_base
import sage.combinat.partition
from sage.modules.free_module_element import vector
from sage.matrix.all import matrix, MatrixSpace

class CombinatorialAlgebraElement(AlgebraElement):
    def __init__(self, A, x):
        """
        Create a combinatorial algebra element x.  This should never
        be called directly, but only through the parent combinatorial
        algebra's __call__ method.
        """
        AlgebraElement.__init__(self, A)
        self._monomial_coefficients = x


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

    def _repr_(self):
        v = self._monomial_coefficients.items()
        v.sort()
        prefix = self.parent().prefix()
        mons = [ prefix + str(m) for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def _latex_(self):
        v = self._monomial_coefficients.items()
        v.sort()
        prefix = self.parent().prefix()
        mons = [ prefix + '_{' + ",".join(map(str, m)) + '}' for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs, is_latex=True).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def _eq_(self, right):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: b = s([2,1])
            sage: c = s([2,1]) + 0
            sage: d = s([2,1]) + 3
            sage: a == b
            True
            sage: b == c
            True
            sage: c == d
            False
        """
        for b in self.parent()._combinatorial_class:
            if self.coefficient(b) != right.coefficient(b):
                return False
        return True

    def __cmp__(left, right):
        """
        The ordering is the one on the underlying sorted list of (monomial,coefficients) pairs.
        """
        v = left._monomial_coefficients.items()
        v.sort()
        w = right._monomial_coefficients.items()
        w.sort()
        return cmp(v, w)

    def _add_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]) + s([5,4])
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

        z = A(Integer(0))
        z._monomial_coefficients = z_elt
        return z

    def _neg_(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: -s([2,1])
            -s[2, 1]
        """
        return self.map_coefficients(lambda c: -c)


    def _sub_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]) - s([5,4])
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

        z = A(Integer(0))
        z._monomial_coefficients = z_elt
        return z

    def _mul_(self, y):
        return self.parent().multiply(self, y)

    def _div_(self, y):
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
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2])^2
            s[2, 2] + s[3, 1] + s[4]

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
        for i in range(n):
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
        if m in p._combinatorial_class:
            return self._monomial_coefficients.get(p._combinatorial_class.object_class(m), p.base_ring()(0))
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
        """
        res = self.parent()(0)
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            new_m, new_c = f(m,c)
            z_elt[new_m] = new_c
        res._monomial_coefficients = z_elt
        return res

class CombinatorialAlgebra(Algebra):
    """
    """
    def __init__(self, R, element_class=None):
        """
        INPUT:
            R -- ring
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
        except:
            self._dim = None

        self._order = None

        #Initialize the base structure
        sage.structure.parent_base.ParentWithBase.__init__(self, R)

    _prefix = ""
    _name   = "CombinatorialAlgebra -- change me"

    def basis(self):
        """
        Returns a list of the basis elements of self.
        """
        return map(self, self._combinatorial_class.list())

    def __call__(self, x):
        """
        Coerce x into self.
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
        elif x in self._combinatorial_class:
            return eclass(self, {self._combinatorial_class.object_class(x):R(1)})
        #Coerce elements of the base ring
        elif hasattr(x, 'parent') and x.parent() is R:
            if x == R(0):
                return eclass(self, {})
            else:
                return eclass(self, {self._combinatorial_class.object_class(self._one):x})
        #Coerce things that coerce into the base ring
        elif R.has_coerce_map_from(x.parent()):
            rx = R(x)
            if rx == R(0):
                return eclass(self, {})
            else:
                return eclass(self, {self._combinatorial_class.object_class(self._one):R(x)})
        else:
            if hasattr(self, '_coerce_end'):
                try:
                    return self._coerce_start(x)
                except TypeError:
                    pass
            raise TypeError, "do not know how to make x (= %s) an element of self (=%s)"%(x,self)


    def _an_element_impl(self):
        return self._element_class(self, {self._combinatorial_class.object_class(self._one):self.base_ring()(1)})

    def _repr_(self):
        return self._name + " over %s"%self.base_ring()

    def combintorial_class(self):
        return self._combinatorial_class

    def _coerce_impl(self, x):
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
        """
        return self._dim

    def set_order(self, order):
        """
        Sets the order of the elements of the combinatorial class.
        If .set_order() has not been called, then the ordering is
        the one used in the generation of the elements of self's
        associated combinatorial class.
        """

        #TODO: add checks

        self._order = order

    def get_order(self):
        return self._order

    def prefix(self):
        return self._prefix


    def _linearize(self, f, a):
        """
        This takes in a function from the basis elements
        to the elements of self.
        """
        mcs = a.monomial_coefficients()
        base_ring = self.base_ring()
        zero = base_ring(0)

        z_elt = {}
        for basis_element in mcs:
            f_mcs = f(basis_element).monomial_coefficients()
            print f_mcs
            for f_basis_element in f_mcs:
                z_elt[ f_basis_element ] = z_elt.get(f_basis_element, zero) + mcs[basis_element]*f_mcs[f_basis_element]

        res = self(0)
        print z_elt
        res._monomial_coefficients = z_elt
        print res.monomial_coefficients()
        return res

    def multiply(self,left,right):
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

        z = A(Integer(0))
        z._monomial_coefficients = z_elt
        return z


#class TestCA(CombinatorialAlgebra):
#    _combinatorial_class = partition.Partitions()
#    _name = "Test Combinatorial Algebra"
#    _one = partition.Partition([])
#    def _multiply_basis(self, left_basis, right_basis):
#        #Append the two lists
#        m = list(left_basis) + list(right_basis)
#        #Sort the new lists
#        m.sort(reverse=True)
#        return partition.Partition(m)
