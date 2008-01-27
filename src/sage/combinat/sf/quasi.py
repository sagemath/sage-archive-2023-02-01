r"""
Quasisymmetric Functions
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

from combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
import word
from permutation import Permutations
from composition import Compositions, Composition


#Generic classes
class QuasiSymmetricFunctions_generic(CombinatorialAlgebra):
    pass
class QuasiSymmetricFunction_generic(CombinatorialAlgebraElement):
    pass

#Monomials -- dual of NCSF::Complete
class QuasiSymmetricFunctions_monomial(QuasiSymmetricFunctions_generic):
    def __init__(self, R):
        self._combinatorial_class = Compositions()
        self._name = "Quasisymmetric functions in the monomial basis"
        self._one = Composition([])
        self._prefix = "M"
        self._element_class = QuasiSymmetricFunction_monomial
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        z_elt = {}
        one = self.base_ring()(1)
        zero = self.base_ring()(0)
        for c in word.ShuffleProduct(left, right, overlapping=True):
            cc = Composition(c)
            z_elt[ cc ] = z_elt.get(cc, zero) + one
        return z_elt

class QuasiSymmetricFunction_monomial(QuasiSymmetricFunction_generic):
    pass
##     def expand(self, n, alphabet='x'):
##         resPR = PolynomialRing(self.parent().base_ring(), n, alphabet)
##         res = resPR(0)
##         self_mc = self._monomial_coefficients
##         for c in self_mc:
##             #Expand c
##             for
##             res += self_mc[c] * resPR(e(part, n, alphabet))
##         return res


class QuasiSymmetricFunctions_ribbon(QuasiSymmetricFunctions_generic):
    def __init__(self, R):
        self._combinatorial_class = Compositions()
        self._name = "Quasisymmetric functions in the monomial basis"
        self._one = Composition([])
        self._prefix = "F"
        self._element_class = QuasiSymmetricFunction_ribbon
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        p1 = Permutations(descents=left).first()
        p2 = Permutations(descents=right).first()
        k = len(p1)
        p2 = [k+i for i in p2]

        z_elt = {}
        one = self.base_ring()(1)
        zero = self.base_ring()(0)
        for c in word.ShuffleProduct(p1, p2):
            cc = Permutation(c).descents_composition()
            z_elt[ cc ] = z_elt.get(cc, zero) + one
        return z_elt

class QuasiSymmetricFunction_ribbon(QuasiSymmetricFunction_generic):
    pass


#Psi -- dual of NCSF::Psi
class QuasiSymmetricFunctions_psi(QuasiSymmetricFunctions_generic):
    def __init__(self, R):
        self._combinatorial_class = Compositions()
        self._name = "Quasisymmetric functions in the psi basis"
        self._one = Composition([])
        self._prefix = "psi"
        self._element_class = QuasiSymmetricFunction_psi
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        z_elt = {}
        one = self.base_ring()(1)
        zero = self.base_ring()(0)
        for c in word.ShuffleProduct(left, right):
            cc = Composition(c)
            z_elt[ cc ] = z_elt.get(cc, zero) + one
        return z_elt

class QuasiSymmetricFunction_psi(QuasiSymmetricFunction_generic):
    pass


#Phi -- dual of NCSF::Phi
class QuasiSymmetricFunctions_phi(QuasiSymmetricFunctions_generic):
    def __init__(self, R):
        self._combinatorial_class = Compositions()
        self._name = "Quasisymmetric functions in the phi basis"
        self._one = Composition([])
        self._prefix = "phi"
        self._element_class = QuasiSymmetricFunction_phi
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        z_elt = {}
        one = self.base_ring()(1)
        zero = self.base_ring()(0)
        for c in word.ShuffleProduct(left, right):
            cc = Composition(c)
            z_elt[ cc ] = z_elt.get(cc, zero) + one
        return z_elt

class QuasiSymmetricFunction_phi(QuasiSymmetricFunction_generic):
    pass


#G -- Hall-Littlewood Aanalogue
class QuasiSymmetricFunctions_g(QuasiSymmetricFunctions_generic):
    def __init__(self, R):
        self._combinatorial_class = Compositions()
        self._name = "Quasisymmetric functions in the g basis"
        self._one = Composition([])
        self._prefix = "G"
        self._element_class = QuasiSymmetricFunction_g
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        raise NotImplementedError

class QuasiSymmetricFunction_g(QuasiSymmetricFunction_generic):
    pass

#Q -- Gessel's quasisymmetric functions
class QuasiSymmetricFunctions_g(QuasiSymmetricFunctions_generic):
    def __init__(self, R):
        self._combinatorial_class = Compositions()
        self._name = "Quasisymmetric functions in the g basis"
        self._one = Composition([])
        self._prefix = "Q"
        self._element_class = QuasiSymmetricFunction_g
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        raise NotImplementedError

    def _coerce_start(self):
        pass

class QuasiSymmetricFunction_g(QuasiSymmetricFunction_generic):
    pass
