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


class CombinatorialAlgebraElement(AlgebraElement):
    def __init__(self, A, x):
        """
        """
        AlgebraElement.__init__(self, A)
        self._monomial_coefficients = x


    def monomial_coefficients(self):
        return self._monomial_coefficients

    def _repr_(self):
        v = self._monomial_coefficients.items()

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
        y = self.parent()(0)
        y_elt = {}
        for m, c in self._monomial_coefficients.iteritems():
            y_elt[m] = -c
        y._monomial_coefficients = y_elt
        return y

    def _sub_(self, y):
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

    def __pow__(self, n):
        """

        """
        if not isinstance(n, (int, Integer)):
            raise TypeError, "n must be an integer"
        A = self.parent()
        z = A(Integer(1))
        for i in range(n):
            z *= self
        return z

    def _coefficient_fast(self, m, default):
        return self._monomial_coefficients.get(m, default)

    def coefficient(self, m):
        """

        """
        p = self.parent()
        if m in p._combinatorial_class:
            return self._monomial_coefficients.get(p._combinatorial_class.object_class(m), p.base_ring()(0))
        else:
            raise TypeError, "you must specify an element of %s"%p._combinatorial_class

    def is_zero(self):
        BR = self.parent().base_ring()
        for v in self._monomial_coefficients.values():
            if v != BR(0):
                return False
        return True

    def __len__(self):
        return self.length()

    def length(self):
        """
        EXAMPLES:
        """
        return len( filter(lambda x: self._monomial_coefficients[x] != 0, self._monomial_coefficients) )

    def support(self):
        """
        EXAMPLES:
        """
        v = self._monomial_coefficients.items()
        mons = [ m for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        return [mons, cffs]

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
            except:
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
        elif x.parent() is R:
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
                except:
                    pass
            raise TypeError, "do not know how to make x (= %s) an element of self"%(x)


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

    def prefix(self):
        return self._prefix

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
                    #the coefficients
                    if isinstance(res, dict):
                        for m in res:
                            if m in z_elt:
                                z_elt[ m ] = z_elt[m] + left_c * right_c * res[m]
                            else:
                                z_elt[ m ] = left_c * right_c * res[m]
                    #Otherwise, res is assumed to be an element of the
                    #object class
                    else:
                        m = res
                        if m  in z_elt:
                            z_elt[ m ] = z_elt[m] + left_c * right_c
                        else:
                            z_elt[ m ] = left_c * right_c
        #We assume that the user handles the multiplication correctly on
        #his or her own, and returns a dict with monomials as keys and
        #coefficients as values
        else:
            z_elt = self._multiply(left, right)

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
