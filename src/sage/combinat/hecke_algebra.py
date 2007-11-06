r"""
Hecke Algebra
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

from combinatorial_algebra import CombinatorialAlgebra
import permutation
from sage.rings.all import PolynomialRing

class HeckeAlgebra_generic(CombinatorialAlgebra):
    def __init__(self, R, n, q=None):
        """
        TESTS:
            #sage: QH3 = HeckeAlgebra(QQ, 3)
            #sage: QH3 == loads(dumps(QH3))
            #True
        """
        self.n = n
        self._combinatorial_class = permutation.Permutations(n)
        self._name = "Hecke Algebra of order %s"%self.n
        self._one = permutation.Permutation(range(1,n+1))
        self._prefix = ""

        if q is None:
            q = PolynomialRing(R, 'q').fraction_field().gen()
            R = q.parent()
        else:
            if q not in R:
                raise ValueError, "q must be in R (= %s)"%R

        self._q = q

        CombinatorialAlgebra.__init__(self, R)

    def q(self):
        return self._q

    def q_inverse(self):
        return ~self._q

    def _coerce_start(self, x):
        ###################################################
        # Coerce permutations of size smaller that self.n #
        ###################################################
        if x == []:
            return self( self._one )
        if len(x) < self.n and x in permutation.Permutations():
            return self( list(x) + range(len(x)+1, self.n+1) )
        raise TypeError



class HeckeAlgebra_t(HeckeAlgebra_generic):
    def __init__(self, R, n, q=None):
        HeckeAlgebra_generic.__init__(self, R, n, q)
        self._prefix = "T"
        self._name += "on the T_sigma basis"

    def t_action_on_basis(self, perm, i):
        if i not in range(1, self.n):
            raise ValueError, "i must be between 1 and n (= %s)"%self.n

        t_i = permutation.Permutation( (i, i+1) )

        perm_i = t_i * perm

        print self(perm_i)

        if perm[i-1] < perm[i]:
            return self(perm_i)
        else:
            #Ti^2 = (q + q^(-1))*Ti - q1*q2
            q = self.q()
            qi = self.q_inverse()

            z_elt = {perm_i:-q*qi, perm:q+qi}
            res = self(0)
            res._monomial_coefficients = z_elt
            return res
            #return (-q*qi)*self(perm_i) + (q+qi)*self(perm)

    def t_action(self, a, i):
        t_i = lambda x: self.t_action_on_basis(x, i)
        return self._linearize(t_i, a)


    def _multiply_basis(self, perm1, perm2):
        res = self(perm1)
        for i in perm2.reduced_word():
            res = self.t_action(res, i)
        return res

    def t(self, i):
        if i not in range(1, self.n):
            raise ValueError, "i must be between 1 and n-1 (= %s)"%(self.n-1)

        return self( permutation.Permutation( (i, i+1) ) )

    def algebra_generators(self):
        """
        Return the generators of the algebra.
        """
        return map(self.t, range(1, self.n))

    def jucys_murphy(self, k):
        """
        Returns the Jucys-Murphy element J_k of the Hecke algebra.
        The Jucys-Murphy elements generate the maximal commutative
        sub-algebra of the Hecke algebra.
        """
        if k not in range(2, self.n+1):
            raise ValueError, "k must be between 2 and n (= %s)"%self.n

        left = 1
        right = 1
        for j in range(1, k):
            left *= self.t(k-j)
            right *= self.t(j)

        return left*right

    def _t_to_y(self, sigm, parameters=None):
        """
        Compute the change of basis between the usual basis
        and the Yang-Baxter basis.
        """

        if parameters is None:
            parameters = PolynomialRing(self.base_ring(), self.n, 'x').gens()

        w = sigma.reduced_word()
        res = self(1)
        perm = permutation.Permutation_class(range(1, self.n+1))
        q = self.q()
        qi = self.q_inverse()
        for i in range(len(w)):
            res *= (self.t(w[i])-(q+qi))/(1-parameters[perm[w[i]]]/parameters[perm[w[i]]-1])

        return res

