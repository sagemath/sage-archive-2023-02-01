r"""
Symmetric Functions

AUTHOR: Mike Hansen, 2007-06-15

sage: s = SymmetricFunctionAlgebra(QQ, basis='schur')
sage: e = SymmetricFunctionAlgebra(QQ, basis='elementary')
sage: f1 = s([2,1])
sage: f1
s[2, 1]
sage: f2 = e(f1)
sage: f2
e[2, 1] - e[3]
sage: f1 == f2
False
sage: f1.expand(3, alphabet=['x','y','z'])
x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2
sage: f2.expand(3, alphabet=['x','y','z'])
x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2


sage: m = SFAMonomial(QQ)
sage: m([3,1])
m[3, 1]
sage: m(4)
4*m[]
sage: m([4])
m[4]
sage: 3*m([3,1])-1/2*m([4])
3*m[3, 1] - 1/2*m[4]


Code needs to be added to coerce symmetric polynomials into symmetric functions.

sage: p = SFAPower(QQ)
sage: m = p(3)
sage: m
3*p[]
sage: m.parent()
Symmetric Algebra over Rational Field, Power symmetric functions as basis
sage: m + p([3,2])
3*p[] + p[3, 2]


sage: s = SFASchur(QQ)
sage: h = SFAHomogeneous(QQ)
sage: P = SFAPower(QQ)
sage: e = SFAElementary(QQ)
sage: m = SFAMonomial(QQ)
sage: a = s([3,1])
sage: s(a)
s[3, 1]
sage: h(a)
h[3, 1] - h[4]
sage: p(a)
1/8*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] - 1/8*p[2, 2] - 1/4*p[4]
sage: e(a)
e[2, 1, 1] - e[2, 2] - e[3, 1] + e[4]
sage: m(a)
3*m[1, 1, 1, 1] + 2*m[2, 1, 1] + m[2, 2] + m[3, 1]
sage: a.expand(4)
x0^3*x1 + x0^2*x1^2 + x0*x1^3 + x0^3*x2 + 2*x0^2*x1*x2 + 2*x0*x1^2*x2 + x1^3*x2 + x0^2*x2^2 + 2*x0*x1*x2^2 + x1^2*x2^2 + x0*x2^3 + x1*x2^3 + x0^3*x3 + 2*x0^2*x1*x3 + 2*x0*x1^2*x3 + x1^3*x3 + 2*x0^2*x2*x3 + 3*x0*x1*x2*x3 + 2*x1^2*x2*x3 + 2*x0*x2^2*x3 + 2*x1*x2^2*x3 + x2^3*x3 + x0^2*x3^2 + 2*x0*x1*x3^2 + x1^2*x3^2 + 2*x0*x2*x3^2 + 2*x1*x2*x3^2 + x2^2*x3^2 + x0*x3^3 + x1*x3^3 + x2*x3^3


sage: h(m([1]))
h[1]
sage: h( m([2]) +m([1,1]) )
h[2]
sage: h( m([3]) + m([2,1]) + m([1,1,1]) )
h[3]
sage: h( m([4]) + m([3,1]) + m([2,2]) + m([2,1,1]) + m([1,1,1,1]) )
h[4]
sage: k = 5
sage: h( sum([ m(part) for part in Partitions(k)]) )
h[5]
sage: k = 10
sage: h( sum([ m(part) for part in Partitions(k)]) )
h[10]

sage: P3 = Partitions(3)
sage: P3.list()
[[3], [2, 1], [1, 1, 1]]
sage: m = SFAMonomial(QQ)
sage: f = sum([m(p) for p in P3])
sage: m.get_print_style()
'lex'
sage: f
m[1, 1, 1] + m[2, 1] + m[3]
sage: m.set_print_style('length')
sage: f
m[3] + m[2, 1] + m[1, 1, 1]
sage: m.set_print_style('maximal_part')
sage: f
m[1, 1, 1] + m[2, 1] + m[3]


sage: s = SFASchur(QQ)
sage: m = SFAMonomial(QQ)
sage: m([3])*s([2,1])
2*m[3, 1, 1, 1] + m[3, 2, 1] + m[4, 2] + 2*m[4, 1, 1] + m[5, 1]
sage: s(m([3])*s([2,1]))
s[2, 1, 1, 1, 1] - s[2, 2, 2] - s[3, 3] + s[5, 1]
sage: s(s([2,1])*m([3]))
s[2, 1, 1, 1, 1] - s[2, 2, 2] - s[3, 3] + s[5, 1]
sage: e = SFAElementary(QQ)
sage: e([4])*e([3])*e([1])
e[4, 3, 1]


sage: s = SFASchur(QQ)
sage: z = s([2,1]) + s([1,1,1])
sage: z.coefficient([2,1])
1
sage: z.length()
2
sage: z.support()
[[[1, 1, 1], [2, 1]], [1, 1]]
sage: z.degree()
3




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

from sage.rings.ring import Ring
from sage.rings.integer import Integer

from sage.algebras.algebra import Algebra

import sage.combinat.partition
import sage.combinat.skew_partition
import sage.structure.parent_gens
import sage.libs.symmetrica.all as symmetrica
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
from sage.matrix.constructor import matrix

from sage.rings.integer_ring import IntegerRing
from sage.rings.rational_field import RationalField

from sage.misc.misc import repr_lincomb, prod
from sage.algebras.algebra_element import AlgebraElement

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

import operator


ZZ = IntegerRing()
QQ = RationalField()



def SymmetricFunctionAlgebra(R, basis="schur"):
    """
    Return the free algebra over the ring $R$ on $n$ generators with
    given names.

    INPUT:
        R -- ring with identity
        basis

    OUTPUT:
        A SymmetricFunctionAlgebra

    EXAMPLES:
        sage: SymmetricFunctionAlgebra(QQ)
        Symmetric Algebra over Rational Field, Schur symmetric functions as basis


        sage: SymmetricFunctionAlgebra(QQ, basis='m')
        Symmetric Algebra over Rational Field, Monomial symmetric functions as basis


        sage: SymmetricFunctionAlgebra(QQ, basis='power')
        Symmetric Algebra over Rational Field, Power symmetric functions as basis

    """
    if basis == 'schur' or basis == 's':
        return cache_s(R)
    elif basis == "elementary" or  basis ==  'e':
        return cache_e(R)
    elif basis == "homogeneous" or basis ==  'h':
        return cache_h(R)
    elif basis == 'power' or basis ==  'p':
        return cache_p(R)
    elif basis == 'monomial' or basis ==  'm':
        return cache_m(R)
    else:
        raise ValueError, "unknown basis (= %s)"%basis

def SFAPower(R):
    """
    Returns the symmetric function algebra over R with the power-sum symmetric
    functions as the basis.

    EXAMPLES:
        sage: SFAPower(QQ)
        Symmetric Algebra over Rational Field, Power symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='power')


def SFAElementary(R):
    """
    Returns the symmetric function algebra over R with the elementary symmetric
    functions as the basis.

    EXAMPLES:
        sage: SFAElementary(QQ)
        Symmetric Algebra over Rational Field, Elementary symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='elementary')


def SFAHomogeneous(R):
    """
    Returns the symmetric function algebra over R with the Homogeneous symmetric
    functions as the basis.

    EXAMPLES:
        sage: SFAHomogeneous(QQ)
        Symmetric Algebra over Rational Field, Homogeneous symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='homogeneous')


def SFASchur(R):
    """
    Returns the symmetric function algebra over R with the Schur symmetric
    functions as the basis.

    EXAMPLES:
        sage: SFASchur(QQ)
        Symmetric Algebra over Rational Field, Schur symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='schur')


def SFAMonomial(R):
    """
    Returns the symmetric function algebra over R with the monomial symmetric
    functions as the basis.

    EXAMPLES:
        sage: SFAMonomial(QQ)
        Symmetric Algebra over Rational Field, Monomial symmetric functions as basis
    """
    return SymmetricFunctionAlgebra(R, basis='monomial')

def is_SymmetricFunctionAlgebra(x):
    """
    Return True if x is a symmetric function algebra; otherwise, return False.

    EXAMPLES:
        sage: sage.combinat.sf.sfa.is_SymmetricFunctionAlgebra(5)
        False
        sage: sage.combinat.sf.sfa.is_SymmetricFunctionAlgebra(ZZ)
        False
        sage: sage.combinat.sf.sfa.is_SymmetricFunctionAlgebra(SymmetricFunctionAlgebra(ZZ,'schur'))
        True
    """
    return isinstance(x, SymmetricFunctionAlgebra_generic)



def zee(part):
    p = sage.combinat.partition.Partition_class(part)
    return p.centralizer_size()


def is_SymmetricFunction(x):
    return isinstance(x, SymmetricFunctionAlgebraElement_generic)

class SymmetricFunctionAlgebra_generic(CombinatorialAlgebra):
    def basis_name(self):
        return None

    def _change_by_proportionality(self, x, function):
        """
        INPUT:
            x: a symmetric function
            function: a function which takes in a partition and returns
                      a scalar

        OUTPUT:
            a symmetric function in self which is a scaled version
            of x
        """
        BR = self.base_ring()
        z_elt = {}
        for m, c in x._monomial_coefficients.iteritems():
            coeff = function(m)
            z_elt[m] = BR( c*coeff )
        z = self(0)
        z._monomial_coefficients = z_elt
        return z

    def _apply_multi_module_morphism(self, x, y, f, orthogonal=False):
        """
        INPUT:
            -- x : a element of self
            -- y : a element of self
            -- f : a function that takes in two partitions (basis elements)
                   and returns an element of the target domain
            -- orthogonal: if orthogonal is set to True, then f(part1, part2)
                           is assumed to be 0 if part1 != part2.
        """
        res = 0
        if orthogonal:
            for mx, cx in x._monomial_coefficients.iteritems():
                if mx not in y._monomial_coefficients:
                    continue
                res += cx*y._monomial_coefficients[mx]*f(mx, mx)
            return res
        else:
            for mx, cx in x._monomial_coefficients.iteritems():
                for my, cy in y._monomial_coefficients.iteritems():
                    res += cx*cy*f(mx,my)
                return res

    def _apply_module_morphism(self, x, f):
        """
        Returns the image of x under the module morphism defined by
        extending f by linearity.

        INPUT:
            -- x : a element of self
            -- f : a function that takes in a partition (basis element)
                   and returns an element of the target domain

        """
        res = 0
        for m, c in x._monomial_coefficients.iteritems():
            res += c*f(m)
        return res


    def _from_dict(self, d):
        """
        Given a monomial coefficient dictionary d, return the element
        of self with the dictionary.

        """
        res = self(0)
        res._monomial_coefficients = d
        return res

    def _from_element(self, x):
        """
        Return the element of self with the same 'internal structure' as
        x.
        """
        return self._from_dict(x.monomial_coefficients())


class SymmetricFunctionAlgebraElement_generic(CombinatorialAlgebraElement):
    def __call__(self, x):
        """
        Plethysm.

        This is inefficient right now as it not only does it term by term, it also converts both arguments
        into the schur basis before computing the plethysm.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: h = SFAHomogeneous(QQ)
            sage: s ( h([3])( h([2]) ) )
            s[2, 2, 2] + s[4, 2] + s[6]
            sage: p = SFAPower(QQ)
            sage: p([3])( s([2,1]) )
            1/3*p[3, 3, 3] - 1/3*p[9]
            sage: e = SFAElementary(QQ)
            sage: e([3])( e([2]) )
            e[3, 3] + e[4, 1, 1] - 2*e[4, 2] - e[5, 1] + e[6]


            sage: R.<t> = QQ[]
            sage: s = SFASchur(R);
            sage: a = s([3])
            sage: f = t*s([2])
            sage: a(f)
            t^3*s[2, 2, 2] + t^3*s[4, 2] + t^3*s[6]
            sage: f(a)
            t*s[4, 2] + t*s[6]

        """

        if not is_SymmetricFunction(x):
            raise TypeError, "only know how to compute plethysms between symmetric functions"

        R = self.parent().base_ring()
        zero = R(0)
        gens = R.gens()
        s = SFASchur(R)

        self_schur = s(self)
        x_schur = s(x)

        self_mcs = self_schur.monomial_coefficients()
        x_mcs    = x_schur.monomial_coefficients()
        import pdb
        z_elt = {}
        for self_part in self_mcs:
            for x_part in x_mcs:
                #Make sure that the non-one generators for
                #degree one element for the plethysm
                coeff = x_mcs[x_part]
                for gen in gens:
                    if gen != 1:
                        #This is due to bad behavior with subs /
                        #substitute.  A bug report will be filed,
                        #and the bug should be fixed in a future
                        #version
                        if len(gens) > 1:
                            coeff = coeff.subs(gen=gen**sum(self_part))
                        else:
                            coeff = coeff.substitute(gen**sum(self_part))
                coeff = self_mcs[self_part] * coeff

                plet = symmetrica.schur_schur_plet(self_part, x_part)
                plet_mcs = plet.monomial_coefficients()
                for plet_part in plet_mcs:
                    z_elt[ plet_part ] = z_elt.get( plet_part, zero ) + plet_mcs[plet_part]*coeff

        z = s(0)
        z._monomial_coefficients = z_elt

        return self.parent()(z)

    def frobenius(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.
        """
        raise NotImplementedError

    def omega(self):
        """
        An alias for self.frobenius() .
        """
        return self.frobenius()

    def theta(self,a):
        """
        Returns the image of self under the theta automorphism
        which sends $p[k]$ to $a*p[k]$.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]).theta(2)
            2*s[1, 1, 1] + 6*s[2, 1] + 2*s[3]
            sage: p = SFAPower(QQ)
            sage: p([2]).theta(2)
            2*p[2]
        """
        p = SFAPower(self.parent().base_ring())
        p_self = p(self)
        res = p_self.map_mc(lambda m,c: (m, c*a**len(m)))
        return self.parent()(res)

    def theta_qt(self,q,t):
        """
        Returns the image of self under the theta automorphism
        which sends $p[k]$ to $(1-q^k)/(1-t^k)*p[k]$.

        EXAMPLES:
            sage: QQqt = QQ['q,t'].fraction_field()
            sage: q,t = QQqt.gens()
            sage: p = SFAPower(QQqt)
            sage: p([2]).theta_qt(q,t)
            ((-q^2+1)/(-t^2+1))*p[2]
            sage: p([2,1]).theta_qt(q,t)
            ((q^3-q^2-q+1)/(t^3-t^2-t+1))*p[2, 1]

        """
        BR = self.parent().base_ring()
        p = SFAPower(BR)
        p_self = p(self)
        res = p_self.map_mc(lambda m,c: (m, BR(prod([(1-q**k)/(1-t**k) for k in m])*c)))
        return self.parent()(res)


    def itensor(self, x):
        """
        Returns the inner tensor product of self and x
        in the basis of self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: b = s([3])
            sage: a.itensor(b)
            s[2, 1]
            sage: c = s([3,2,1])
            sage: c.itensor(c)
            s[1, 1, 1, 1, 1, 1] + 2*s[2, 1, 1, 1, 1] + 3*s[2, 2, 1, 1] + 2*s[2, 2, 2] + 4*s[3, 1, 1, 1] + 5*s[3, 2, 1] + 2*s[3, 3] + 4*s[4, 1, 1] + 3*s[4, 2] + 2*s[5, 1] + s[6]

        """
        #Convert both self and x to the p basis
        p = SFAPower(self.parent().base_ring())
        p_self = p(self)
        p_x    = p(x)

        #Determine whether p_self or p_x has
        #the smallest support
        if len(p_self) <= len(p_x):
            g = p_self
            f = p_x
        else:
            g = p_x
            f = p_self

        #Get the internal data structures for
        #efficientcy
        g_mcs = g._monomial_coefficients
        f_mcs = f._monomial_coefficients


        #Compute the inner tensor product
        res = 0
        for part in g_mcs:
            if part in f_mcs:
                c = zee(part)
                res += c*g_mcs[part]*f_mcs[part]*p(part)

        #Return the result in terms of self's basis
        return self.parent()(res)


    internal_product = itensor
    kronecker_product = itensor
    inner_tensor = itensor

#############
#   Cache   #
#############
from sage.misc.cache import Cache
import schur, monomial, powersum, elementary, homogeneous
cache_s = Cache(schur.SymmetricFunctionAlgebra_schur)
cache_m = Cache(monomial.SymmetricFunctionAlgebra_monomial)
cache_p = Cache(powersum.SymmetricFunctionAlgebra_power)
cache_e = Cache(elementary.SymmetricFunctionAlgebra_elementary)
cache_h = Cache(homogeneous.SymmetricFunctionAlgebra_homogeneous)
