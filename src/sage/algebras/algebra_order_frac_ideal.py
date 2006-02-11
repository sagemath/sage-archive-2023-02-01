"""
Fractional ideals over orders in an algebra

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

class AlgebraOrderFractionalIdeal_generic(object):
    """
    Generic left, right, or two-sided algebra order ideals.
    """
    def __init__(self, A, O1=None, O2=None, R=None, gens = []):
        if not isinstance(R, AlgebraOrder):
            raise TypeError, "Argument R must be an algebra."
        self.__algebra = A
        self.__left_algebra = left
        self.__right_algebra = right
        self.__base_ring = R
        self.__gens = gens

    def __repr__(self):
        NotImplementedError

    def __call__(self, x):
        raise NotImplementedError

    def __contains__(self, x):
        raise NotImplementedError

    def left_order(self):
        return self.__left_order

    def right_order(self):
        return self.__right_order

    def base_ring(self):
        return self.__base_ring

    def gen(self,i):
        n = ngens(self)
        if i < 0 or not i < n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        return self.__gens[i]

    def ngens(self):
        return len(self.__gens)

class AlgebraOrderLeftFractionalIdeal(AlgebraOrderFractionalIdeal_generic):
    """
    A fractional left ideal in an algebra order.
    """
    def __init__(self, A, O, gens = []):
        if not isinstance(O, AlgebraOrder):
            raise TypeError, "Argument O must be an algebra order."
        R = O.base_ring()
        AlgebraOrderFractionalIdeal_generic.__init__(self, A, O, None, R, gens)

    def __repr__(self):
        return "Left fractional ideal with generators %s over %s"%(
            self.gens, self.right_order)

class AlgebraOrderRightFractionalIdeal(AlgebraOrderFractionalIdeal_generic):
    """
    A fractional right ideal in an algebra order.
    """
    def __init__(self, A, O, gens = []):
        if not isinstance(O, AlgebraOrder):
            raise TypeError, "Argument O must be an algebra order."
        R = O.base_ring()
        AlgebraOrderFractionalIdeal_generic.__init__(self, A, None, O, R, gens)

    def __repr__(self, maximal=False):
        if maximal:
            return "Right fractional ideal with generators %s over %s"%(
                self.gens(), self.right_order())
        else:
            return "Right fractional ideal with generators %s"%(self.gens())

class AlgebraOrderFractionalBiIdeal(AlgebraOrderFractionalIdeal_generic):
    """
    A fractional left and right ideal over orders in an algebra.
    """
    def __init__(self, A, O1, O2, gens = []):
        if not isinstance(O1, AlgebraOrder) and isinstance(O2, AlgebraOrder):
            raise TypeError, "Argument O must be an algebra order."
        if O1.ambient_algebra() is O2.ambient_algebra():
            raise ValueError, \
                  "Arguments O1 and O2 must be orders in a common algebra."
        R = O.base_ring()
        AlgebraOrderFractionalIdeal_generic.__init__(self, A, O1, O2, R, gens)

    def __repr__(self, maximal=False):
        if maximal:
            return "Fractional bi-ideal with generators %s over left order %s and right order %s"%(
            self.gens(), self.left_order(), self.right_order())
        else:
            return "Fractional bi-ideal with generators %s"%(self.gens())

class AlgebraOrderFractionalTwoSidedIdeal(AlgebraOrderFractionalIdeal_generic):
    """
    A fractional two-sided ideal for an order in an algebra.
    """
    def __init__(self, A, O, gens = []):
        if not isinstance(O1, AlgebraOrder) and isinstance(O2, AlgebraOrder):
            raise TypeError, "Argument O must be an algebra order."
        R = O.base_ring()
        AlgebraOrderFractionalIdeal_generic.__init__(self, A, O, O, R, gens)

    def __repr__(self, maximal=False):
        if maximal:
            return "Fractional two-sided ideal with generators %s over %s"%(
                self.gens(), self.left_order())
        else:
            return "Fractional two-sided ideal with generators %s"%(self.gens())

