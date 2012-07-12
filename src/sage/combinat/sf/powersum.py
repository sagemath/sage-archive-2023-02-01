"""
Power sum symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#                     2012 Anne Schilling <anne at math.ucdavis.edu>
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
from sage.combinat.partition import Partition

class SymmetricFunctionAlgebra_power(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, Sym):
        """
        A class for methods associated to the power sum basis of the symmetric functions

        INPUT:

        - ``self`` -- the power sum basis of the symmetric functions
        - ``Sym`` -- an instance of the ring of symmetric functions

        TESTS::

            sage: p = SymmetricFunctions(QQ).p()
            sage: p == loads(dumps(p))
            True
            sage: TestSuite(p).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(p).run(elements = [p[1,1]+p[2], p[1]+2*p[1,1]])
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, Sym, "power", 'p')

    def coproduct_on_generators(self, i):
        r"""
        Returns coproduct on generators for power sums `p_i`.
        The elements `p_i` are primitive elements.

        INPUT:

        - ``self`` -- the power sum basis of the symmetric functions
        - ``i`` -- a positive integer

        OUTPUT:

        - the result of the coproduct on the generator `p(i)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.powersum()
            sage: p.coproduct_on_generators(2)
            p[] # p[2] + p[2] # p[]
        """
        def P(k): return Partition([k]) if k else Partition([])
        T = self.tensor_square()
        return T.sum_of_monomials( [(P(i), P(0)), (P(0), P(i))] )

    def antipode_on_basis(self, partition):
        r"""
        Returns the antipode of ``self[partition]``.
        The antipode on the generator `p_i` is `-p_i` and the
        antipode on `p_\mu` is `(-1)^{length(\mu)} p_\mu`.

        INPUT:

        - ``self`` -- the power sum basis of the symmetric functions
        - ``partition`` -- a partition

        OUTPUT:

        - the result of the antipode on ``self(partition)``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: p = Sym.p()
            sage: p.antipode_on_basis([2])
            -p[2]
            sage: p.antipode_on_basis([3])
            -p[3]
            sage: p.antipode_on_basis([2,2])
            p[2, 2]
        """
        return (-1)**len(partition) * self[partition]

    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        def omega(self):
            """
            Returns the image of ``self`` under the Frobenius / omega automorphism.

            INPUT:

            - ``self`` -- an element of the symmetric functions in the power sum basis

            OUTPUT:

            - the image of ``self`` under the omega automorphism

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: a = p([2,1]); a
                p[2, 1]
                sage: a.omega()
                -p[2, 1]
                sage: p([]).omega()
                p[]
                sage: p(0).omega()
                0
            """
            f = lambda part, coeff: (part, (-1)**(sum(part)-len(part))*coeff)
            return self.map_item(f)

        def scalar(self, x, zee=None):
            r"""
            Returns the standard scalar product of ``self`` and ``x``.

            INPUT:

            - ``self`` -- an element of the symmetric functions in the power sum basis
            - ``x`` -- an power sum symmetric function
            - ``zee`` -- optional input specifying the scalar product on the power sum basis
              with normalization `<p_\mu, p_\mu> = zee(\mu)`. ``zee`` should be
              a function on partitions. (default: uses standard ``zee`` function)

            Note that the power-sum symmetric functions are orthogonal under
            this scalar product. The value of `\langle p_\lambda, p_\lambda \rangle`
            is given by the size of the centralizer in `S_n` of a
            permutation of cycle type `\lambda`.

            OUTPUT:

            - the standard scalar product between ``self`` and ``x`` or if the optional
              parameter ``zee`` is specified then the scalar product with respect to the
              normalization `<p_\mu, p_\mu> = zee(\mu)` with the power sum bases
              elements are orthogonal

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: p4 = Partitions(4)
                sage: matrix([ [p(a).scalar(p(b)) for a in p4] for b in p4])
                [ 4  0  0  0  0]
                [ 0  3  0  0  0]
                [ 0  0  8  0  0]
                [ 0  0  0  4  0]
                [ 0  0  0  0 24]
                sage: p(0).scalar(p(1))
                0
                sage: p(1).scalar(p(2))
                2

                sage: zee = lambda x : 1
                sage: matrix( [[p[la].scalar(p[mu], zee) for la in Partitions(3)] for mu in Partitions(3)])
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            parent = self.parent()
            x = parent(x)
            if zee is None:
                f = lambda part1, part2:  sfa.zee(part1)
            else:
                f = lambda part1, part2:  zee(part1)
            return parent._apply_multi_module_morphism(self, x, f, orthogonal=True)


        def _derivative(self, part):
            """
            Returns the 'derivative' of `p([part])` with respect to `p([1])`.

            INPUT:

            - ``self`` -- an element of the symmetric functions in the power sum basis
            - ``part`` -- a partition

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: a = p([2,1])
                sage: a._derivative(Partition([2,1]))
                p[2]
                sage: a._derivative(Partition([1,1,1]))
                3*p[1, 1]
            """
            p = self.parent()
            if 1 not in part:
                return p(0)
            else:
                return len([i for i in part if i == 1])*p(part[:-1])

        def _derivative_with_respect_to_p1(self):
            """
            Returns the `derivative` of a symmetric function in the power sum basis with respective to `p([1])`.
            On the Frobenius image of an `S_n` module the resulting character is
            the Frobenius image of the restriction of this module to `S_{n-1}`.

            INPUT:

            - ``self`` -- an element of the symmetric functions in the power sum basis

            OUPUT:

            - a power symmetric function of degree one smaller corresponding
              differentiation by `p_1`.

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
                sage: a = p([2,1,1,1])
                sage: a._derivative_with_respect_to_p1()
                3*p[2, 1, 1]
                sage: a = p([3,2])
                sage: a._derivative_with_respect_to_p1()
                0
                sage: p(0)._derivative_with_respect_to_p1()
                0
                sage: p(1)._derivative_with_respect_to_p1()
                0
                sage: p([1])._derivative_with_respect_to_p1()
                p[]
                sage: f = p[1] + p[2,1]
                sage: f._derivative_with_respect_to_p1()
                p[] + p[2]
            """
            p = self.parent()
            if self == p.zero():
                return self
            return p._apply_module_morphism(self, self._derivative)

        def expand(self, n, alphabet='x'):
            """
            Expands the symmetric function as a symmetric polynomial in `n` variables.

            INPUT:

            - ``self`` -- an element of the symmetric functions in the power sum basis
            - ``n`` -- a positive integer
            - ``alphabet`` -- a variable for the expansion (default: `x`)

            OUTPUT:

            - a polynomial expansion of an instance of ``self`` in `n` variables

            EXAMPLES::

                sage: p = SymmetricFunctions(QQ).p()
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
                sage: p(1).expand(4)
                1
                sage: p(0).expand(4)
                0
                sage: (p([]) + 2*p([1])).expand(3)
                2*x0 + 2*x1 + 2*x2 + 1
            """
            condition = lambda part: False
            return self._expand(condition, n, alphabet)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.powersum', 'SymmetricFunctionAlgebraElement_power',  SymmetricFunctionAlgebra_power.Element)
