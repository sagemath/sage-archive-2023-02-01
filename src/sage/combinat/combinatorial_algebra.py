r"""
Combinatorial Algebras

A combinatorial algebra is an algebra whose basis elements are
indexed by a combinatorial class. Some examples of combinatorial
algebras are the symmetric group algebra of order n (indexed by
permutations of size n) and the algebra of symmetric functions
(indexed by integer partitions).

The CombinatorialAlgebra base class makes it easy to define and
work with new combinatorial algebras in Sage. For example, the
following code constructs an algebra which models the power-sum
symmetric functions.

::

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

::

    sage: ps = PowerSums(QQ); ps
    Power-sum symmetric functions over Rational Field
    sage: ps([2,1])^2
    p[2, 2, 1, 1]
    sage: ps([2,1])+2*ps([1,1,1])
    2*p[1, 1, 1] + p[2, 1]
    sage: ps(2)
    2*p[]

The important things to define are ._combinatorial_class which
specifies the combinatorial class that indexes the basis elements,
._one which specifies the identity element in the algebra, ._name
which specifies the name of the algebra, ._prefix which is the
string put in front of each basis element, and finally a _multiply
or _multiply basis method that defines the multiplication in the
algebra.
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

from sage.rings.all import Integer
from sage.algebras.algebra import Algebra
from sage.rings.ring import Ring
from sage.algebras.algebra_element import AlgebraElement
from sage.matrix.all import MatrixSpace
from sage.combinat.free_module import CombinatorialFreeModuleElement, CombinatorialFreeModule, CombinatorialFreeModuleInterface
from sage.misc.misc import repr_lincomb

class CombinatorialAlgebraElement(AlgebraElement, CombinatorialFreeModuleElement):
    def __init__(self, A, x):
        """
        Create a combinatorial algebra element x. This should never be
        called directly, but only through the parent combinatorial
        algebra's __call__ method.

        TESTS::

            sage: s = SFASchur(QQ)
            sage: a = s._element_class(s, {Partition([2,1]):QQ(2)}); a
            2*s[2, 1]
            sage: a == loads(dumps(a))
            True
        """
        AlgebraElement.__init__(self, A)
        self._monomial_coefficients = x


    def _mul_(self, y):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2])
            sage: a._mul_(a) #indirect doctest
            s[2, 2] + s[3, 1] + s[4]
        """
        return self.parent().multiply(self, y)

    def _div_(self, y):
        """
        EXAMPLES::

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
        EXAMPLES::

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
        Returns self to the `n^{th}` power.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s([2])^2
            s[2, 2] + s[3, 1] + s[4]

        TESTS::

            sage: s = SFASchur(QQ)
            sage: z = s([2,1])
            sage: z
            s[2, 1]
            sage: z^2
            s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

        ::

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

    def _matrix_(self, new_BR = None):
        """
        Returns a matrix version of self obtained by the action of self on
        the left. If new_BR is specified, then the matrix will be over
        new_BR.

        EXAMPLES::

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

        MS = MatrixSpace(new_BR, parent.dimension(), parent.dimension())
        l = [ (self*parent(m)).to_vector() for m in cc ]
        return MS(l).transpose()

    def to_matrix(self):
        """
        Returns a matrix version of self obtained by the action of self on
        the left. If new_BR is specified, then the matrix will be over
        new_BR.

        EXAMPLES::

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


    def __repr__(self):
        """
        EXAMPLES::

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


class CombinatorialAlgebra(CombinatorialFreeModuleInterface, Algebra):
    def __init__(self, R, element_class = None):
        """
        TESTS::

            sage: s = SFASchur(QQ)
            sage: s == loads(dumps(s))
            True
        """
        #Check to make sure that the user defines the necessary
        #attributes / methods to make the combinatorial module
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

        CombinatorialFreeModuleInterface.__init__(self, R, self._element_class)

    def __call__(self, x):
        try:
            return CombinatorialFreeModuleInterface.__call__(self, x)
        except TypeError:
            pass

        R = self.base_ring()
        eclass = self._element_class
        #Coerce elements of the base ring
        if not hasattr(x, 'parent'):
            x = R(x)
        if x.parent() == R:
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

        raise TypeError, "do not know how to make x (= %s) an element of self (=%s)"%(x,self)

    def _an_element_impl(self):
        """
        Returns an element of self, namely the unit element.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s._an_element_impl()
            s[]
            sage: _.parent() is s
            True
        """
        return self._element_class(self, {self._one:self.base_ring()(1)})

    def _coerce_impl(self, x):
        """
        EXAMPLES::

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

    def multiply(self,left,right):
        """
        Returns left\*right where left and right are elements of self.
        multiply() uses either _multiply or _multiply basis to carry out
        the actual multiplication.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2])
            sage: s.multiply(a,a)
            s[2, 2] + s[3, 1] + s[4]
            sage: ZS3 = SymmetricGroupAlgebra(ZZ, 3)
            sage: a = 2 + ZS3([2,1,3])
            sage: a*a
            5*[1, 2, 3] + 4*[2, 1, 3]
            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ,3)
            sage: j2 = H3.jucys_murphy(2)
            sage: j2*j2
            (q^3-q^2+q)*T[1, 2, 3] + (q^3-q^2+q-1)*T[2, 1, 3]
            sage: X = SchubertPolynomialRing(ZZ)
            sage: X([1,2,3])*X([2,1,3])
            X[2, 1]
        """
        A = left.parent()
        ABR = A.base_ring()
        ABRzero = ABR(0)
        z_elt = {}

        #Do the case where the user specifies how to multiply basis
        #elements
        if hasattr(self, '_multiply_basis'):
            for (left_m, left_c) in left._monomial_coefficients.iteritems():
                for (right_m, right_c) in right._monomial_coefficients.iteritems():
                    res = self._multiply_basis(left_m, right_m)
                    coeffprod = left_c * right_c
                    #Handle the case where the user returns a dictionary
                    #where the keys are the monomials and the values are
                    #the coefficients.  If res is not a dictionary, then
                    #it is assumed to be an element of self
                    if not isinstance(res, dict):
                        if isinstance(res, self._element_class):
                            res = res._monomial_coefficients
                        else:
                            z_elt[res] = z_elt.get(res, ABRzero) + coeffprod
                            continue
                    for m, c in res.iteritems():
                        z_elt[m] = z_elt.get(m, ABRzero) + coeffprod * c

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
