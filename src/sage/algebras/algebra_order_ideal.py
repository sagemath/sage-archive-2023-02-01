"""
Algebra order ideals

AUTHOR: David Kohel, 2005-09
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
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

from sage.modules.free_module import FreeModule, VectorSpace
from algebra import Algebra
from algebra_order import AlgebraOrder_generic
from algebra_ideal import AlgebraLeftIdeal, AlgebraRightIdeal
from sage.algebras.algebra_element import AlgebraElement
from algebra_order_ideal_element import AlgebraOrderIdealElement

class AlgebraOrderIdeal_generic(object):
    """
    A algebra order ideal.
    """
    def __init__(self, A, R, gens = [], left_order=None, right_order=None):
        if not isinstance(A, Algebra):
            raise TypeError, "Argument A (= %s) must be an algebra."%A
        self.__ambient_algebra = A
        self.__base_ring = R
        self.__gens = gens
        self.__left_order = left_order
        self.__right_order = right_order

    def __repr__(self):
        NotImplementedError

    def __call__(self, x, check=True):
        if isinstance(x, AlgebraOrderIdealElement) and x.parent() is self:
            return x
        return AlgebraOrderIdealElement(self, x, check=check)

    def __contains__(self, x):
        if isinstance(x, AlgebraOrderIdealElement) and x.parent() is self:
	    return True
        elif isinstance(x, AlgebraElement) \
               and x.parent() is self.ambient_algebra():
            return NotImplementedError
        elif isinstance(x, RingElement) and x in self.base_ring():
            return True
        else:
            return False

    def ambient_algebra(self):
        return self.__ambient_algebra

    def base_ring(self):
        return self.__base_ring()

    def basis(self):
        try:
            return tuple([ self(x) for x in self.__algebra_basis_elements ])
        except AttributeError:
            if isinstance(self,AlgebraOrderTwoSidedIdeal):
                self.__compute_ideal_basis(left=True, right=True)
            elif isinstance(self,AlgebraOrderLeftIdeal):
                self.__compute_ideal_basis(left=True)
            elif isinstance(self,AlgebraOrderRightIdeal):
                self.__compute_ideal_basis(right=True)
            else:
                raise RuntimeError, \
                      "Bug: unassigned left, right or two-sided order action"
        return tuple([ self(x, check=False) for x in self.__algebra_basis_elements ])


    def __compute_ideal_basis(self, left=False, right=False):
        """
        Compute a basis as a left ideal.
        """
        A = self.__ambient_algebra
        if left:
            O = self.__left_order
        else:
            O = self.__right_order
        R = self.__base_ring
        V = A.vector_space()
        M = FreeModule(R,O.rank())
        N = M.submodule([]) # syntax should change
        try:
            order_basis = O._AlgebraOrder_generic__algebra_basis_elements
        except AttributeError:
            O._AlgebraOrder_generic__compute_order_basis()
            order_basis = O._AlgebraOrder_generic__algebra_basis_elements
        ideal_basis = self.__gens
        saturated = False
        while not saturated:
            saturated = True
            for x in ideal_basis:
                for y in order_basis:
                    if left:
                        v = (y*x).vector()
                        if not v in N:
                            saturated = False
                            N += M.submodule([ v ])
                    if right:
                        v = (x*y).vector()
                        if not v in N:
                            saturated = False
                            N += M.submodule([ v ])
            ideal_basis = [ A(v.list()) for v in N.basis() ]
            ideal_basis.reverse()
        self.__algebra_basis_elements = ideal_basis

    def gen(self,i):
        try:
            return self.__gens[i]
        except IndexError:
            raise IndexError, \
                  "Argument i (= %s) must be between 0 and %s."%(i, ngens(self)-1)

    def gens(self):
        return tuple([ self(x) for x in self.__gens ])

    def module(self):
        try:
            return self.__module
        except AttributeError:
            M = FreeModule(self.__base_ring,self.__rank)
            vecs = [ x.vector() for x in self.__algebra_basis_elements ]
            # This syntax should change to V.submodule(R,vecs)
            self.__module = M.submodule(vecs)
        return self.__module

    def ngens(self):
        return len(self.__gens)

    def left_order(self):
        return self.__left_order

    def right_order(self):
        return self.__right_order

    def rank(self):
        return self.__rank


class AlgebraOrderLeftIdeal(AlgebraOrderIdeal_generic):
    """
    A left ideal over an order in an algebra.
    """
    def __init__(self, O, gens = []):
        if not isinstance(O, AlgebraOrder_generic):
            raise TypeError, \
                  "Argument O (= %s) must be an order in an order in an algebra."%O
        A = O.ambient_algebra()
        R = O.base_ring()
        AlgebraOrderIdeal_generic.__init__(self, A, R, gens, left_order=O)

    def __repr__(self):
        return "Left ideal generated by %s over %s"%(
            self._AlgebraOrderIdeal_generic__gens, self.left_order())

class AlgebraOrderRightIdeal(AlgebraOrderIdeal_generic):
    """
    A right ideal over an order in an algebra.
    """
    def __init__(self, O, gens = []):
        if not isinstance(O, AlgebraOrder_generic):
            raise TypeError, "Argument O (= %s) must be an order in an algebra."%O
        A = O.ambient_algebra()
        R = O.base_ring()
        AlgebraOrderIdeal_generic.__init__(self, A, R, gens, right_order=O)

    def __repr__(self):
        return "Right ideal generated by %s over %s"%(
            self._AlgebraOrderIdeal_generic__gens, self.right_order())


class AlgebraOrderTwoSidedIdeal(AlgebraOrderIdeal_generic):
    """
    A right ideal over an order in an algebra.
    """
    def __init__(self, O, gens = []):
        if not isinstance(O, AlgebraOrder_generic):
            raise TypeError, "Argument O (= %s) must be an order in an algebra."%O
        A = O.ambient_algebra()
        R = O.base_ring()
        AlgebraOrderIdeal_generic.__init__(
            self, A, R, gens, left_order=O, right_order=O)

    def __repr__(self):
        return "Two-sided ideal generated by %s over %s"%(
            self._AlgebraOrderIdeal_generic__gens, self.right_order())

