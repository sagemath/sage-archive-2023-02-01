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

#Print style
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
2*m[3, 1, 1, 1] + m[3, 2, 1] + 2*m[4, 1, 1] + m[4, 2] + m[5, 1]
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

import partition
import skew_partition
import sage.structure.parent_gens
import sage.libs.symmetrica.all as symmetrica
from sage.matrix.constructor import matrix

from sage.rings.integer_ring import IntegerRing
from sage.rings.rational_field import RationalField

from sage.misc.misc import repr_lincomb
from sage.algebras.algebra_element import AlgebraElement

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

import operator
import hall_littlewood

ZZ = IntegerRing()
QQ = RationalField()

translate = {'monomial':'MONOMIAL', 'homogeneous':'HOMSYM', 'power':'POWSYM', 'elementary':'ELMSYM', 'schur':'SCHUR'}

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
    return SymmetricFunctionAlgebra(R, basis='power')
def SFAElementary(R):
    return SymmetricFunctionAlgebra(R, basis='elementary')
def SFAHomogeneous(R):
    return SymmetricFunctionAlgebra(R, basis='homogeneous')
def SFASchur(R):
    return SymmetricFunctionAlgebra(R, basis='schur')
def SFAMonomial(R):
    return SymmetricFunctionAlgebra(R, basis='monomial')

def is_SymmetricFunctionAlgebra(x):
    """
    Return True if x is a symmetric function algebra; otherwise, return False.

    EXAMPLES:
        sage: sage.combinat.sfa.is_SymmetricFunctionAlgebra(5)
        False
        sage: sage.combinat.sfa.is_SymmetricFunctionAlgebra(ZZ)
        False
        sage: sage.combinat.sfa.is_SymmetricFunctionAlgebra(SymmetricFunctionAlgebra(ZZ,'schur'))
        True
    """
    return isinstance(x, SymmetricFunctionAlgebra_generic)


class SymmetricFunctionAlgebra_generic(Algebra):
    """
    EXAMPLES:
    """
    def __init__(self, R, basis, element_class, prefix):
        """
        INPUT:
            R -- ring
        """
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        try:
            z = R(Integer(1))
        except:
            raise ValueError, "R must have a unit element"

        self.__basis = basis
        self.__element_class = element_class
        self.__prefix = prefix
        self.__print_style = 'lex'
        sage.structure.parent_gens.ParentWithGens.__init__(self, R, None)

    def __call__(self, x):
        """
        Coerce x into self.
        """
        R = self.base_ring()
        eclass = self.__element_class
        if isinstance(x, int):
            x = Integer(x)

        #same basis
        if x in partition.Partitions():
            return eclass(self, {partition.Partition(x):R(1)})
        elif isinstance(x, eclass):
            P = x.parent()
            #same base ring
            if P is self:
                return x
            #different base ring
            else:
                return eclass(self, dict([ (e1,R(e2)) for e1,e2 in x._monomial_coefficients.items()]))

        #symmetric function, but different basis
        elif is_SymmetricFunction(x):
            R = self.base_ring()
            xP = x.parent()
            xm = x.monomial_coefficients()

            #determine the conversion function.
            try:
                t = eval('symmetrica.t_' + translate[xP.basis()] + '_' +  translate[self.basis()])
            except AttributeError:
                raise TypeError, "do not know how to convert from %s to %s"%(xP.basis(), self.basis())

            z_elt = {}

            for part in xm:
                if xm[part] == R(0):
                    continue
                xmprime = t( {part:Integer(1)} ).monomial_coefficients()
                for part2 in xmprime:
                    z_elt[part2] = z_elt.get(part2, R(0)) + xmprime[part2]*R(xm[part])
            z = self(Integer(0))
            z._monomial_coefficients = z_elt
            return z
        #############################
        #Hall-Littlewood Polynomials#
        #############################
        #Qp: Convert to Schur basis and then convert to self
        elif isinstance(x, hall_littlewood.HallLittlewoodElement_qp):
            Qp = x.parent()
            BR = self.base_ring()
            zero = BR(0)
            s = SFASchur(BR)
            QpBR = Qp.base_ring()
            if not BR.has_coerce_map_from(QpBR):
                raise TypeError, "no coerce map from x's parent's base ring (= %s) to self's base ring (= %s)"%(self.QpBR, self.BR)

            z_elt = {}
            for m, c in x._monomial_coefficients.iteritems():
                n = sum(m)
                Qp._s_cache(n)
                for part in Qp._qp_to_s_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + BR(c*Qp._qp_to_s_cache[n][m][part])

            z = s(0)
            z._monomial_coefficients = z_elt
            return self(z)
        #P: Convert to monomials and then convert to self
        elif isinstance(x, hall_littlewood.HallLittlewoodElement_p):
            res = 0
            for m, c in x._monomial_coefficients.iteritems():
                n = sum(m)
                P._m_cache(n)
                res += c*P._p_to_m_cache[m]
            return self(z)


        elif x in skew_partition.SkewPartitions():
            skewschur = symmetrica.part_part_skewschur(x[0], x[1])
            return self(skewschur)
        elif isinstance(x, list):
            if len(x) == 2 and isinstance(x[0], list) and isinstance(x[1], list):
                skewschur = symmetrica.part_part_skewschur(x[0], x[1])
                return self(skewschur)
            else:
                return eclass(self, {partition.Partition(x):R(1)})
        elif x.parent() is R:
            return eclass(self, {partition.Partition([]):x})
        elif R.has_coerce_map_from(x.parent()):
            return eclass(self, {partition.Partition([]):R(x)})
        else:
            try:
                return eclass(self, {Partition([]):self.base_ring()(x)})
            except:
                raise TypeError, "do not know how to make x (= %s) an element of self"%(x)

    def basis(self):
        return self.__basis

    def is_field(self):
        if self.__ngens == 0:
            return self.base_ring().is_field()
        return False

    def is_commutative(self):
        """
        Return True if this symmetric function algebra is commutative.

        EXAMPLES:

        """
        return self.__ngens <= 1 and self.base_ring().is_commutative()

    def _an_element_impl(self):
        return self.__element_class(self, {partition.Partition([]):self.base_ring()(0)})

    def __cmp__(self, other):
        """
        Two free algebras are considered the same if they have the
        same base ring, number of generators and variable names.

        EXAMPLES:

        """
        if not isinstance(other, SymmetricFunctionAlgebra_generic):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        c = cmp(self.__basis, other.__basis)
        if c: return c
        return 0

    def get_print_style(self):
        return self.__print_style

    def set_print_style(self, ps):
        styles = ['lex', 'length', 'maximal_part']
        if ps not in styles:
            raise ValueError, "the print style must be one of ", styles
        self.__print_style = ps

    def _repr_(self):
        """
        Text representation of this symmetric function algebra.

        EXAMPLES:
        """
        return "Symmetric Algebra over %s, %s symmetric functions as basis"%(
            self.base_ring(), self.__basis.capitalize())

    #def _corece_impl(self, x):
    #     raise NotImplementedError
    def _coerce_impl(self, x):
        try:
            R = x.parent()
            #Coerce other symmetric functions in
            if is_SymmetricFunctionAlgebra(R):
                #Only perform the coercion if we can go from the base
                #ring of x to the base ring of self
                if self.base_ring().has_coerce_map_from( R.base_ring() ):
                    return self(x)
        except AttributeError:
            pass

        # any ring that coerces to the base ring of this free algebra.
        return self._coerce_try(x, [self.base_ring()])

    def ngens(self):
        """
        The number of generators of the algebra.

        EXAMPLES:
        """
        return infinity

    def prefix(self):
        return self.__prefix

    def transition_matrix(self, basis, n):
        """
        Returns the transitions matrix.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: m = SFAMonomial(QQ)
            sage: s.transition_matrix(m,5)
            [1 1 1 1 1 1 1]
            [0 1 1 2 2 3 4]
            [0 0 1 1 2 3 5]
            [0 0 0 1 1 3 6]
            [0 0 0 0 1 2 5]
            [0 0 0 0 0 1 4]
            [0 0 0 0 0 0 1]

            sage: p = SFAPower(QQ)
            sage: s.transition_matrix(p, 4)
            [ 1/4  1/3  1/8  1/4 1/24]
            [-1/4    0 -1/8  1/4  1/8]
            [   0 -1/3  1/4    0 1/12]
            [ 1/4    0 -1/8 -1/4  1/8]
            [-1/4  1/3  1/8 -1/4 1/24]
            sage: StoP = s.transition_matrix(p,4)
            sage: a = s([3,1])+5*s([1,1,1,1])-s([4])
            sage: a
            5*s[1, 1, 1, 1] + s[3, 1] - s[4]
            sage: mon, coef = a.support()
            sage: coef
            [5, -1, 1]
            sage: mon
            [[1, 1, 1, 1], [4], [3, 1]]
            sage: cm = matrix([[-1,1,0,0,5]])
            sage: cm * StoP
            [-7/4  4/3  3/8 -5/4 7/24]
            sage: p(a)
            7/24*p[1, 1, 1, 1] - 5/4*p[2, 1, 1] + 3/8*p[2, 2] + 4/3*p[3, 1] - 7/4*p[4]


            sage: h = SFAHomogeneous(QQ)
            sage: e = SFAElementary(QQ)
            sage: s.transition_matrix(m,7) == h.transition_matrix(s,7).transpose()
            True

            sage: h.transition_matrix(m, 7) == h.transition_matrix(m, 7).transpose()
            True

            sage: h.transition_matrix(e, 7) == e.transition_matrix(h, 7)
            True


            sage: p.transition_matrix(s, 5)
            [ 1 -1  0  1  0 -1  1]
            [ 1  0 -1  0  1  0 -1]
            [ 1 -1  1  0 -1  1 -1]
            [ 1  1 -1  0 -1  1  1]
            [ 1  0  1 -2  1  0  1]
            [ 1  2  1  0 -1 -2 -1]
            [ 1  4  5  6  5  4  1]

            sage: e.transition_matrix(m,7) == e.transition_matrix(m,7).transpose()
            True


        """
        P = partition.Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( map( lambda col_part: z.coefficient(col_part), Plist ) )
        return matrix(m)


class SymmetricFunctionAlgebra_schur(SymmetricFunctionAlgebra_generic):
    def __init__(self, R):
        SymmetricFunctionAlgebra_generic.__init__(self, R, "schur", SymmetricFunctionAlgebraElement_schur, 's')
    def is_schur_basis(self):
        return True



class SymmetricFunctionAlgebra_monomial(SymmetricFunctionAlgebra_generic):
    def __init__(self, R):
        SymmetricFunctionAlgebra_generic.__init__(self, R, "monomial", SymmetricFunctionAlgebraElement_monomial, 'm')

class SymmetricFunctionAlgebra_elementary(SymmetricFunctionAlgebra_generic):
    def __init__(self, R):
        SymmetricFunctionAlgebra_generic.__init__(self, R, "elementary", SymmetricFunctionAlgebraElement_elementary, 'e')

class SymmetricFunctionAlgebra_power(SymmetricFunctionAlgebra_generic):
    def __init__(self, R):
        SymmetricFunctionAlgebra_generic.__init__(self, R, "power", SymmetricFunctionAlgebraElement_power, 'p')

class SymmetricFunctionAlgebra_homogeneous(SymmetricFunctionAlgebra_generic):
    def __init__(self, R):
        SymmetricFunctionAlgebra_generic.__init__(self, R, "homogeneous", SymmetricFunctionAlgebraElement_homogeneous, 'h')

from sage.misc.cache import Cache
cache_s = Cache(SymmetricFunctionAlgebra_schur)
cache_m = Cache(SymmetricFunctionAlgebra_monomial)
cache_p = Cache(SymmetricFunctionAlgebra_power)
cache_e = Cache(SymmetricFunctionAlgebra_elementary)
cache_h = Cache(SymmetricFunctionAlgebra_homogeneous)



############
# Elements #
############

def zee(part):
    p = partition.Partition_class(part)
    return p.centralizer_size()


def is_SymmetricFunction(x):
    return isinstance(x, SymmetricFunctionAlgebraElement_generic)


class SymmetricFunctionAlgebraElement_generic(AlgebraElement):
    """
    A symmetric function element.
    """
    def __init__(self, A, x):
        """
        Create a symmetric function x.  This should never be called directly, but only through the
        symmetric function algebra's __call__ method.
        """
        AlgebraElement.__init__(self, A)
        self._monomial_coefficients = x

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

        """

        if not is_SymmetricFunction(x):
            raise TypeError, "only know how to compute plethysms between symmetric functions"

        R = self.parent().base_ring()
        s = SFASchur(R)

        if self.parent().basis() == 'schur':
            self_schur = self
        else:
            self_schur = s(self)

        if x.parent().basis() == 'schur':
            x_schur = x
        else:
            x_schur = s(x)


        self_mcs = self_schur.monomial_coefficients()
        x_mcs    = x_schur.monomial_coefficients()

        z_elt = {}


        for self_part in self_mcs:
            for x_part in x_mcs:
                scalar = self_mcs[self_part]*x_mcs[x_part]
                plet = symmetrica.schur_schur_plet(self_part, x_part)
                plet_mcs = plet.monomial_coefficients()
                for plet_part in plet_mcs:
                    z_elt[ plet_part ] = z_elt.get( plet_part, R(0)) + plet_mcs[plet_part]*scalar

        z = s(0)
        z._monomial_coefficients = z_elt

        if self.parent().basis() == 'schur':
            return z
        else:
            return self.parent()(z)

    def monomial_coefficients(self):
        return self._monomial_coefficients

    def _repr_(self):
        v = self._monomial_coefficients.items()

        ps = self.parent().get_print_style()
        if ps == 'lex':
            v.sort(key=lambda x: x[0])
        if ps == 'length':
            v.sort(key=lambda x: len(x[0]))
        if ps == 'maximal_part':
            def lmax(x):
                if x[0] == []:
                    return 0
                else:
                    return max(x[0])
            v.sort(key=lmax)

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

        ps = self.parent().get_print_style()
        if ps == 'lex':
            v.sort(key=lambda x: x[0])
        if ps == 'length':
            v.sort(key=lambda x: len(x[0]))
        if ps == 'maximal_part':
            def lmax(x):
                if x[0] == []:
                    return 0
                else:
                    return max(x[0])
            v.sort(key=lmax)

        prefix = self.parent().prefix()
        mons = [ prefix + '_{' + ",".join(map(str, m)) + '}' for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs, is_latex=True).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def __cmp__(left, right):
        """
        Compare two symmetric functions with the same parents.

        The ordering is the one on the underlying sorted list of (monomial,coefficients) pairs.

        EXAMPLES:
        """
        v = left._monomial_coefficients.items()
        v.sort()
        w = right._monomial_coefficients.items()
        w.sort()
        return cmp(v, w)

    def _add_(self, y):
        A = self.parent()
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
        z = A(Integer(0))
        z._monomial_coefficients = z_elt
        return z

    def _mul_(self, y):
        raise NotImplementedError

    def __pow__(self, n):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([2,1])
            sage: z
            s[2, 1]
            sage: z^2
            s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

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
        for i in range(n):
            z *= self
        return z

    def _coefficient_fast(self, m, default):
        return self._monomial_coefficients.get(m, default)

    def coefficient(self, m):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) - 2*s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficient([4])
            1
            sage: z.coefficient([2,1])
            -2
        """
        if isinstance(m, partition.Partition_class):
            return self._monomial_coefficients.get(m, self.parent().base_ring()(0))
        elif m in partition.Partitions():
            return self._monomial_coefficients[partition.Partition(m)]
        else:
            raise TypeError, "you must specify a partition"

    def __len__(self):
        return self.length()

    def length(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.length()
            4
        """
        return len( filter(lambda x: self._monomial_coefficients[x] != 0, self._monomial_coefficients) )

    def support(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.support()
            [[[4], [1, 1, 1], [1], [2, 1]], [1, 1, 1, 1]]
        """
        v = self._monomial_coefficients.items()
        mons = [ m for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        return [mons, cffs]

    def degree(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1]) + 3
            sage: z.degree()
            4
        """
        return max( map( sum, self._monomial_coefficients ) )

    def restrict_degree(self, d, exact=True):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_degree(2)
            0
            sage: z.restrict_degree(1)
            s[1]
            sage: z.restrict_degree(3)
            s[1, 1, 1] + s[2, 1]
            sage: z.restrict_degree(3, exact=False)
            s[1] + s[1, 1, 1] + s[2, 1]
        """
        res = self.parent()(0)
        if exact:
            res._monomial_coefficients = dict( filter( lambda x: sum(x[0]) == d, self._monomial_coefficients.items()) )
        else:
            res._monomial_coefficients = dict( filter( lambda x: sum(x[0]) <= d, self._monomial_coefficients.items()) )
        return res

    def restrict_partition_lengths(self, l, exact=True):
        """

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_partition_lengths(2)
            s[2, 1]
            sage: z.restrict_partition_lengths(2, exact=False)
            s[1] + s[2, 1] + s[4]
        """
        res = self.parent()(0)
        if exact:
            res._monomial_coefficients = dict( filter( lambda x: len(x[0]) == l, self._monomial_coefficients.items()) )
        else:
            res._monomial_coefficients = dict( filter( lambda x: len(x[0]) <= l, self._monomial_coefficients.items()) )
        return res

    def restrict_parts(self, n):
        """

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_parts(2)
            s[1] + s[1, 1, 1] + s[2, 1]
            sage: z.restrict_parts(1)
            s[1] + s[1, 1, 1]
        """
        def lmax(x):
            """Return 0 as the max of an empty list."""
            if x[0] == []:
                return 0
            else:
                return max(x[0])
        res = self.parent()(0)
        res._monomial_coefficients = dict( filter( lambda x: lmax(x) <= n, self._monomial_coefficients.items()) )
        return res

    def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n variables.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: a.expand(2)
            x0^2*x1 + x0*x1^2
            sage: a.expand(3)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
            sage: a.expand(4)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + x0^2*x3 + 2*x0*x1*x3 + x1^2*x3 + 2*x0*x2*x3 + 2*x1*x2*x3 + x2^2*x3 + x0*x3^2 + x1*x3^2 + x2*x3^2
            sage: a.expand(2, alphabet='y')
            y0^2*y1 + y0*y1^2
            sage: a.expand(2, alphabet=['a','b'])
            a^2*b + a*b^2

            sage: p = SFAPower(QQ)
            sage: a = p([2])
            sage: a.expand(2)
            x0^2 + x1^2
            sage: a.expand(3, alphabet=['a','b','c'])
            a^2 + b^2 + c^2
            sage: p([2,1,1]).expand(2)
            x0^4 + 2*x0^3*x1 + 2*x0^2*x1^2 + 2*x0*x1^3 + x1^4
        """
        e = eval('symmetrica.compute_' + str(translate[self.parent().basis()]).lower() + '_with_alphabet')
        resPR = PolynomialRing(self.parent().base_ring(), n, alphabet)
        res = resPR(0)
        self_mc = self._monomial_coefficients
        for part in self_mc:
            #if len(part) > n:
            #    continue
            res += self_mc[part] * resPR(e(part, n, alphabet))
        return res

    def frobenius(self):
        raise NotImplementedError

    def omega(self):
        return self.frobenius()

    def scalar(self, x):
        """
        Returns standard scalar product between self and s.

        This is the default implementation that converts both self
        and x into Schur functions and performs the scalar product
        that basis.

        EXAMPLES:
            sage: e = SFAElementary(QQ)
            sage: h = SFAHomogeneous(QQ)
            sage: m = SFAMonomial(QQ)
            sage: p4 = Partitions(4)
            sage: matrix([ [e(a).scalar(h(b)) for a in p4] for b in p4])
            [ 0  0  0  0  1]
            [ 0  0  0  1  4]
            [ 0  0  1  2  6]
            [ 0  1  2  5 12]
            [ 1  4  6 12 24]
            sage: matrix([ [h(a).scalar(e(b)) for a in p4] for b in p4])
            [ 0  0  0  0  1]
            [ 0  0  0  1  4]
            [ 0  0  1  2  6]
            [ 0  1  2  5 12]
            [ 1  4  6 12 24]
            sage: matrix([ [m(a).scalar(e(b)) for a in p4] for b in p4])
            [-1  2  1 -3  1]
            [ 0  1  0 -2  1]
            [ 0  0  1 -2  1]
            [ 0  0  0 -1  1]
            [ 0  0  0  0  1]
            sage: matrix([ [m(a).scalar(h(b)) for a in p4] for b in p4])
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]

        """
        sp = self.parent()
        xp = x.parent()
        BR = sp.base_ring()

        s = SFASchur(BR)
        s_self = s(self)
        s_x = s(x)
        return s_self.scalar(s_x)


    def scalar_hl(self, x, t=None):
        R = self.parent().base_ring()
        p = SFAPower(R)

        p_self = p(self)
        p_x    = p(x)

        if len(p_self) < len(p_x):
            smaller = p_self
            greater = p_x
        else:
            smaller = p_x
            greater = p_self

        res = R(0)
        if t is None:
            Zt = ZZ['t'].fraction_field()
            t = Zt.gen()
        smcs = smaller._monomial_coefficients
        gmcs = greater._monomial_coefficients
        for s_part in smcs :
            if s_part in gmcs:
                res += smcs[s_part]*gmcs[s_part]*s_part.centralizer_size(t=t)

        return res


class SymmetricFunctionAlgebraElement_ehp(SymmetricFunctionAlgebraElement_generic):
    def _mul_(left, right):
        A = left.parent()
        R = A.base_ring()
        z_elt = {}
        for (left_m, left_c) in left._monomial_coefficients.iteritems():
            for (right_m, right_c) in right._monomial_coefficients.iteritems():
                m = list(left_m)+list(right_m)
                m.sort(reverse=True)
                m = partition.Partition(m)
                if m in z_elt:
                    z_elt[ m ] = z_elt[m] + left_c * right_c
                else:
                    z_elt[ m ] = left_c * right_c
        z = A(Integer(0))
        z._monomial_coefficients = z_elt
        return z



class SymmetricFunctionAlgebraElement_schur(SymmetricFunctionAlgebraElement_generic):
    def _mul_(left, right):
        #Use symmetrica to do the multiplication
        A = left.parent()
        R = A.base_ring()

        if  R is ZZ or R is QQ:
            return symmetrica.mult_schur_schur(left, right)

        z_elt = {}
        for (left_m, left_c) in left._monomial_coefficients.iteritems():
            for (right_m, right_c) in right._monomial_coefficients.iteritems():
                d = symmetrica.mult_schur_schur({left_m:Integer(1)}, {right_m:Integer(1)})._monomial_coefficients
                for m in d:
                    if m in z_elt:
                        z_elt[ m ] = z_elt[m] + left_c * right_c * d[m]
                    else:
                        z_elt[ m ] = left_c * right_c * d[m]
        z = A(Integer(0))
        z._monomial_coefficients = z_elt
        return z

    def frobenius(self):
        parent = self.parent()
        z = {}
        mcs = self.monomial_coefficients()
        for part in mcs:
            z[part.conjugate()] = mcs[part]
        res = parent(0)
        res._monomial_coefficients = z
        return res


    def scalar(self, x):
        """
        Returns the standard scalar product between self and x.

        Note that the Schur functions are self-dual with respect
        to this scalar product. They are also lower-triangularly
        related to the monomial symmetric functions with respect
        to this scalar product.

        EXAMPLES:
            sage: s = SFASchur(ZZ)
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: c = 2*s([1,1,1])
            sage: d = a + b
            sage: a.scalar(a)
            1
            sage: b.scalar(b)
            1
            sage: b.scalar(a)
            0
            sage: b.scalar(c)
            2
            sage: c.scalar(c)
            4
            sage: d.scalar(a)
            1
            sage: d.scalar(b)
            1
            sage: d.scalar(c)
            2

            sage: m = SFAMonomial(ZZ)
            sage: p4 = Partitions(4)
            sage: l = [ [s(p).scalar(m(q)) for q in p4] for p in p4]
            sage: matrix(l)
            [ 1  0  0  0  0]
            [-1  1  0  0  0]
            [ 0 -1  1  0  0]
            [ 1 -1 -1  1  0]
            [-1  2  1 -3  1]
        """
        R = self.parent().base_ring()

        if self.parent() != x.parent():
            try:
                x = self.parent()( x )
            except:
                raise TypeError, "cannot compute the scalar product of self and x (= %s)"%x

        if len(self) < len(x):
            smaller = self
            greater = x
        else:
            smaller = x
            greater = self

        res = R(0)
        smcs = smaller._monomial_coefficients
        gmcs = greater._monomial_coefficients
        for s_part in smcs :
            if s_part in gmcs:
                res += smcs[s_part]*gmcs[s_part]

        return res




class SymmetricFunctionAlgebraElement_monomial(SymmetricFunctionAlgebraElement_generic):
    def _mul_(left, right):
        #Use symmetrica to do the multiplication
        A = left.parent()
        R = A.base_ring()

        if  R is ZZ or R is QQ:
            return symmetrica.mult_monomial_monomial(left, right)

        z_elt = {}
        for (left_m, left_c) in left._monomial_coefficients.iteritems():
            for (right_m, right_c) in right._monomial_coefficients.iteritems():
                d = symmetrica.mult_monomial_monomial({left_m:Integer(1)}, {right_m:Integer(1)})
                for m in d:
                    if m in z_elt:
                        z_elt[ m ] = z_elt[m] + left_c * right_c * d[m]
                    else:
                        z_elt[ m ] = left_c * right_c * d[m]
        z = A(Integer(0))
        z._monomial_coefficients = z_elt
        return z

    def frobenius(self):
        parent = self.parent()
        s = SFASchur(parent.base_ring())
        return parent(s(self).frobenius())


class SymmetricFunctionAlgebraElement_elementary(SymmetricFunctionAlgebraElement_ehp):
    def frobenius(self):
        base_ring = self.parent().base_ring()
        h = SFAHomogeneous(base_ring)
        mcs = self.monomial_coefficients()
        res = h(0)
        res._monomial_coefficients = mcs
        return res

class SymmetricFunctionAlgebraElement_homogeneous(SymmetricFunctionAlgebraElement_ehp):
    def frobenius(self):
        base_ring = self.parent().base_ring()
        e = SFAElementary(base_ring)
        mcs = self.monomial_coefficients()
        res = e(0)
        res._monomial_coefficients = mcs
        return res


class SymmetricFunctionAlgebraElement_power(SymmetricFunctionAlgebraElement_ehp):
    def frobenius(self):
        parent = self.parent()
        base_ring = parent.base_ring()
        z = {}
        mcs = self.monomial_coefficients()
        for part in mcs:
            z[part] = (-1)**(sum(part)-len(part))*mcs[part]
        res = parent(0)
        res._monomial_coefficients = z
        return res

    def scalar(self, x):
        """
        Returns the standard scalar product of self and x.

        Note that the power-sum symmetric functions are orthogonal
        under this scalar product.  The value of <p_lambda, p_lambda>
        is given by the size of the centralizer in S_n of a permutation
        of cycle type lambda.

        EXAMPLES:
            sage: p = SFAPower(QQ)
            sage: p4 = Partitions(4)
            sage: matrix([ [p(a).scalar(p(b)) for a in p4] for b in p4])
            [ 4  0  0  0  0]
            [ 0  3  0  0  0]
            [ 0  0  8  0  0]
            [ 0  0  0  4  0]
            [ 0  0  0  0 24]
        """

        R = self.parent().base_ring()

        if self.parent() != x.parent():
            try:
                x = self.parent()( x )
            except:
                raise TypeError, "cannot compute the scalar product of self and x (= %s)"%x

        if len(self) < len(x):
            smaller = self
            greater = x
        else:
            smaller = x
            greater = self

        res = R(0)
        smcs = smaller._monomial_coefficients
        gmcs = greater._monomial_coefficients
        for s_part in smcs :
            if s_part in gmcs:
                res += smcs[s_part]*gmcs[s_part]*zee(s_part)

        return res
