"""
Power-sum symmetric functions
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
import sfa, multiplicative, classical

class SymmetricFunctionAlgebra_power(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, R):
        """
        TESTS::

            sage: p = SFAPower(QQ)
            sage: p == loads(dumps(p))
            True
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, R, "power", SymmetricFunctionAlgebraElement_power, 'p')




class SymmetricFunctionAlgebraElement_power(classical.SymmetricFunctionAlgebraElement_classical):
    def omega(self):
        """
        Returns the image of self under the Frobenius / omega
        automorphism.

        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: a = p([2,1]); a
            p[2, 1]
            sage: a.omega()
            -p[2, 1]
        """
        f = lambda part, coeff: (part, (-1)**(sum(part)-len(part))*coeff)
        return self.map_mc(f)

    def scalar(self, x):
        r"""
        Returns the standard scalar product of self and x.

        Note that the power-sum symmetric functions are orthogonal under
        this scalar product. The value of `\langle p_\lambda, p_\lambda \rangle`
        is given by the size of the centralizer in `S_n` of a
        permutation of cycle type `\lambda`.

        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: p4 = Partitions(4)
            sage: matrix([ [p(a).scalar(p(b)) for a in p4] for b in p4])
            [ 4  0  0  0  0]
            [ 0  3  0  0  0]
            [ 0  0  8  0  0]
            [ 0  0  0  4  0]
            [ 0  0  0  0 24]
        """
        parent = self.parent()
        R = parent.base_ring()
        x = parent(x)
        f = lambda part1, part2:  sfa.zee(part1)
        return parent._apply_multi_module_morphism(self, x, f, orthogonal=True)

    def _derivative(self, part):
        """
        Returns the 'derivative' of p([part]) with respect to p([1]).

        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: a = p([2,1])
            sage: a._derivative(Partition([2,1]))
            p[2]
            sage: a._derivative(Partition([1,1,1]))
            3*p[1, 1]
        """
        p = self.parent()
        if 1 not in part:
            return 0
        else:
            return len([i for i in part if i == 1])*p(part[:-1])

    def _derivative_with_respect_to_p1(self):
        """
        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: a = p([1,1,1])
            sage: a._derivative_with_respect_to_p1()
            3*p[1, 1]
            sage: a = p([3,2])
            sage: a._derivative_with_respect_to_p1()
            0
        """
        p = self.parent()
        return p._apply_module_morphism(self, self._derivative)

    def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n
        variables.

        EXAMPLES::

            sage: p = SFAPower(QQ)
            sage: a = p([2])
            sage: a.expand(2)
            x0^2 + x1^2
            sage: a.expand(3, alphabet=['a','b','c'])
            a^2 + b^2 + c^2
            sage: p([2,1,1]).expand(2)
            x0^4 + 2*x0^3*x1 + 2*x0^2*x1^2 + 2*x0*x1^3 + x1^4
            sage: p([7]).expand(4)
            x0^7 + x1^7 + x2^7 + x3^7
            sage: p([7]).expand(4,alphabet='t')
            t0^7 + t1^7 + t2^7 + t3^7
            sage: p([7]).expand(4,alphabet='x,y,z,t')
            x^7 + y^7 + z^7 + t^7
        """
        condition = lambda part: False
        return self._expand(condition, n, alphabet)
