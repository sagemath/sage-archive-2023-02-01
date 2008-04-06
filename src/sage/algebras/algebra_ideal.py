"""
Algebra left, right, and two-sided ideals

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

from sage.algebras.algebra import Algebra

class AlgebraIdeal(object):
    """
    Generic left, right, and two-sided algebra ideals.

    classes AlgebraLeftIdeal, AlgebraRightIdeal and
    AlgebraTwosidedIdeal are derived from this class
    """
    def __init__(self, A, gens = []):
        """
        Create an ideal in algebra A with given gens
        """
        if not isinstance(A, Algebra): raise TypeError, "Argument A must be an algebra."
        self.__algebra = A
        self.__gens = gens

    def __repr__(self):
        NotImplementedError

    def __call__(self, x):
        raise NotImplementedError

    def __contains__(self, x):
        raise NotImplementedError

    def base_ring(self):
        """
        Return the base ring of this algebra ideal

        EXAMPLES:
        """
        return self.__algebra.base_ring()

    def gen(self,i):
        """
        Return the i'th generator of this algebra ideal

        EXAMPLES:
        """
        if i < 0 or not i < ngens(self):
            raise IndexError, \
                  "Argument i (= %s) must be between 0 and %s."%(i, ngens(self)-1)
        return self.__gens[i]

    def ngens(self):
        """
        Return the number of generators of this algebra ideal

        EXAMPLES:
        """
        return len(self.__gens)

class AlgebraLeftIdeal(AlgebraIdeal):
    """
    A left ideal in an algebra.
    """
    def __init__(self, A, gens = []):
        AlgebraIdeal.__init__(self, A, gens)

    def __repr__(self):
        return "Left ideal on generators %s over %s"%(
            self.gens(), self.algebra())

    def left_algebra(self):
        """
        Return the left algebra of this algebra left-ideal

        EXAMPLES:
        """
        return self.algebra()

class AlgebraRightIdeal(AlgebraIdeal):
    """
    A right ideal in an algebra.
    """
    def __init__(self, A, gens = []):
        AlgebraIdeal.__init__(self, A, gens)

    def __repr__(self):
        return "Right ideal on generators %s over %s"%(
            self.gens(), self.algebra())

    def right_algebra(self):
        """
        Return the right algebra of this algebra right-ideal

        EXAMPLES:
        """
        return self.algebra()

class AlgebraTwoSidedIdeal(AlgebraIdeal):
    """
    A two-sided ideal in an algebra, such that the left and right orders coincide.
    """
    def __init__(self, A, gens = []):
        AlgebraIdeal.__init__(self, A, gens)

    def __repr__(self):
        return "Two-sided ideal on generators %s over %s"%(
            self.gens(), self.algebra())

    def __call__(self, x):
        raise NotImplementedError

    def __contains__(self, x):
        raise NotImplementedError

    def right_algebra(self):
        """
        Return the right algebra of this algebra two-sided-ideal

        EXAMPLES:
        """
        return self.algebra()

    def left_algebra(self):
        """
        Return the right algebra of this algebra two-sided-ideal

        EXAMPLES:
        """
        return self.algebra()
