"""
Classical symmetric functions.
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

from sage.misc.misc import repr_lincomb
from sage.algebras.algebra_element import AlgebraElement

from sage.rings.all import is_Polynomial, PolynomialRing, is_MPolynomial
import operator
import hall_littlewood
import sfa

ZZ = IntegerRing()
QQ = RationalField()

translate = {'monomial':'MONOMIAL', 'homogeneous':'HOMSYM', 'power':'POWSYM', 'elementary':'ELMSYM', 'schur':'SCHUR'}

conversion_functions = {}

def init():
    #Set up the conversion functions
    for other_basis in translate:
        for basis in translate:
            try:
                conversion_functions[(other_basis, basis)] = eval('symmetrica.t_' + translate[other_basis] + '_' +  translate[basis])
            except AttributeError:
                pass


init()


###################################
#                                 #
#  Classical Symmetric Functions  #
#                                 #
###################################
class SymmetricFunctionAlgebra_classical(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, basis, element_class, prefix):
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        try:
            z = R(Integer(1))
        except:
            raise ValueError, "R must have a unit element"

        self.__basis = basis
        self._element_class = element_class
        self._prefix = prefix
        self.__print_style = 'lex'
        self._combinatorial_class = sage.combinat.partition.Partitions_all()
        self._one = sage.combinat.partition.Partition_class([])


        CombinatorialAlgebra.__init__(self, R)



    def __call__(self, x):
        """
        Coerce x into self.
        """
        R = self.base_ring()
        eclass = self._element_class
        if isinstance(x, int):
            x = Integer(x)


        ##############
        # Partitions #
        ##############
        if x in sage.combinat.partition.Partitions_all():
            return eclass(self, {sage.combinat.partition.Partition_class(filter(lambda x: x!=0, x)):R(1)})

        ##############
        # Dual bases #
        ##############
        elif sfa.is_SymmetricFunction(x) and hasattr(x, 'dual'):
            #Check to see if it is the dual of some other basis
            #If it is, try to coerce its corresponding element
            #in the other basis
            return self(x.dual())

        ########################################
        # Symmetric Functions, different basis #
        ########################################
        elif isinstance(x, eclass):
            P = x.parent()
            #same base ring
            if P is self:
                return x
            #different base ring
            else:
                return eclass(self, dict([ (e1,R(e2)) for e1,e2 in x._monomial_coefficients.items()]))

        ##################################################
        # Classical Symmetric Functions, different basis #
        ##################################################
        elif isinstance(x, SymmetricFunctionAlgebraElement_classical):


            R = self.base_ring()
            xP = x.parent()
            xm = x.monomial_coefficients()

            #determine the conversion function.
            try:
                t = conversion_functions[(xP.basis_name(),self.basis_name())]
            except AttributeError:
                raise TypeError, "do not know how to convert from %s to %s"%(xP.basis_name(), self.basis_name())

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

        ###############################
        # Hall-Littlewood Polynomials #
        ###############################
        #
        #Qp: Convert to Schur basis and then convert to self
        #
        elif isinstance(x, hall_littlewood.HallLittlewoodElement_qp):
            Qp = x.parent()
            BR = self.base_ring()
            zero = BR(0)
            s = sfa.SFASchur(BR)
            QpBR = Qp.base_ring()
            if not BR.has_coerce_map_from(QpBR):
                raise TypeError, "no coerce map from x's parent's base ring (= %s) to self's base ring (= %s)"%(QpBR, self.base_ring())

            z_elt = {}
            for m, c in x._monomial_coefficients.iteritems():
                n = sum(m)
                Qp._s_cache(n)
                for part in Qp._qp_to_s_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + BR(c*Qp._qp_to_s_cache[n][m][part].subs(t=Qp.t))

            z = s(0)
            z._monomial_coefficients = z_elt
            return self(z)
        #
        #P: Convert to Schur basis and then convert to self
        #
        elif isinstance(x, hall_littlewood.HallLittlewoodElement_p):
            P = x.parent()
            BR = self.base_ring()
            zero = BR(0)
            s = sfa.SFASchur(BR)
            PBR = P.base_ring()
            if not BR.has_coerce_map_from(PBR):
                raise TypeError, "no coerce map from x's parent's base ring (= %s) to self's base ring (= %s)"%(PBR, self.base_ring())

            z_elt = {}
            for m, c in x._monomial_coefficients.iteritems():
                n = sum(m)
                P._s_cache(n)
                for part in P._p_to_s_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + BR(c*P._p_to_s_cache[n][m][part].subs(t=P.t))

            z = s(0)
            z._monomial_coefficients = z_elt
            return self(z)
        #
        #Q: Convert to P basis and then convert to self
        #
        elif isinstance(x, hall_littlewood.HallLittlewoodElement_q):
            return self( x.parent()._P( x ) )

        ###################
        # Skew Partitions #
        ###################
        elif x in sage.combinat.skew_partition.SkewPartitions():
            skewschur = symmetrica.part_part_skewschur(x[0], x[1])
            return self(skewschur)

        #####################
        # Untyped lists --  #
        #####################
        #elif isinstance(x, list):
        #    if len(x) == 2 and isinstance(x[0], list) and isinstance(x[1], list):
        #        skewschur = symmetrica.part_part_skewschur(x[0], x[1])
        #        return self(skewschur)
        #    else:
        #        return eclass(self, {partition.Partition_class(x):R(1)})

        #############################
        # Elements of the base ring #
        #############################
        elif x.parent() is R:
            return eclass(self, {sage.combinat.partition.Partition_class([]):x})

        ###########################################
        # Elements that coerce into the base ring #
        ###########################################
        elif R.has_coerce_map_from(x.parent()):
            return eclass(self, {sage.combinat.partition.Partition_class([]):R(x)})

        #################################
        # Last shot -- try calling R(x) #
        #################################
        else:
            try:
                return eclass(self, {Partition_class([]):R(x)})
            except:
                raise TypeError, "do not know how to make x (= %s) an element of self"%(x)

    def basis_name(self):
        """
        Returns the name of the basis of self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.basis_name()
            'schur'
            sage: p = SFAPower(QQ)
            sage: p.basis_name()
            'power'
            sage: h = SFAHomogeneous(QQ)
            sage: h.basis_name()
            'homogeneous'
            sage: e = SFAElementary(QQ)
            sage: e.basis_name()
            'elementary'
            sage: m = SFAMonomial(QQ)
            sage: m.basis_name()
            'monomial'
        """
        return self.__basis

    def is_field(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.is_field()
            False
        """
        return False

    def is_commutative(self):
        """
        Return True if this symmetric function algebra is commutative.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.is_commutative()
            True
        """
        return self.base_ring().is_commutative()

    def _an_element_impl(self):
        """
        Returns 0 in self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s._an_element_impl(); a
            0
            sage: a.parent() is s
            True
        """
        return self._element_class(self, {sage.combinat.partition.Partition_class([]):self.base_ring()(0)})

    def __cmp__(self, other):
        """
        EXAMPLES:

        """
        if not isinstance(other, SymmetricFunctionAlgebra_classical):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        c = cmp(self.__basis, other.__basis)
        if c: return c
        return 0

    def get_print_style(self):
        """
        Returns the value of the current print style for self.
        """
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
        return "Symmetric Algebra over %s, %s symmetric functions as basis"%(self.base_ring(), self.__basis.capitalize())

    def _coerce_impl(self, x):
        try:
            R = x.parent()
            #Coerce other symmetric functions in
            if sfa.is_SymmetricFunctionAlgebra(R):
                #Only perform the coercion if we can go from the base
                #ring of x to the base ring of self
                if self.base_ring().has_coerce_map_from( R.base_ring() ):
                    return self(x)
        except AttributeError:
            pass

        # any ring that coerces to the base ring of this free algebra.
        return self._coerce_try(x, [self.base_ring()])


    def prefix(self):
        """
        Returns the prefix on the elements of self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.prefix()
            's'
        """
        return self._prefix

    def transition_matrix(self, basis, n):
        """
        Returns the transitions matrix between self and basis
        for the homogenous component of degree n.

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
            [5, 1, -1]
            sage: mon
            [[1, 1, 1, 1], [3, 1], [4]]
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
        P = sage.combinat.partition.Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( map( lambda col_part: z.coefficient(col_part), Plist ) )
        return matrix(m)


    def dual_basis(self, scalar=None):
        """
        Returns the basis dual to self with respect to the
        scalar product scalar.  If scalar is not specified,
        then it returns the basis dual to the standard scalar
        product for the classical symmetric functions.
        """
        raise NotImplementedError

class SymmetricFunctionAlgebraElement_classical(sfa.SymmetricFunctionAlgebraElement_generic):
    """
    A symmetric function.
    """

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


    def __pow__(self, n):
        """
        EXAMPLES:

        """
        if not isinstance(n, (int, Integer)):
            raise TypeError, "n must be an integer"
        A = self.parent()
        z = A(Integer(1))
        for i in range(n):
            z *= self
        return z



    def degree(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1]) + 3
            sage: z.degree()
            4
        """
        return max( map( sum, self._monomial_coefficients ) + [0] )

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
            sage: h = SFAHomogeneous(QQ)
            sage: h([3]).expand(2)
            x0^3 + x0^2*x1 + x0*x1^2 + x1^3
            sage: h([1,1,1]).expand(2)
            x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3
            sage: h([2,1]).expand(3)
            x0^3 + 2*x0^2*x1 + 2*x0*x1^2 + x1^3 + 2*x0^2*x2 + 3*x0*x1*x2 + 2*x1^2*x2 + 2*x0*x2^2 + 2*x1*x2^2 + x2^3

            sage: p = SFAPower(QQ)
            sage: a = p([2])
            sage: a.expand(2)
            x0^2 + x1^2
            sage: a.expand(3, alphabet=['a','b','c'])
            a^2 + b^2 + c^2
            sage: p([2,1,1]).expand(2)
            x0^4 + 2*x0^3*x1 + 2*x0^2*x1^2 + 2*x0*x1^3 + x1^4
        """
        e = eval('symmetrica.compute_' + str(translate[self.parent().basis_name()]).lower() + '_with_alphabet')
        resPR = PolynomialRing(self.parent().base_ring(), n, alphabet)
        res = resPR(0)
        self_mc = self._monomial_coefficients
        for part in self_mc:
            res += self_mc[part] * resPR(e(part, n, alphabet))
        return res


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

        s = sfa.SFASchur(BR)
        s_self = s(self)
        s_x = s(x)
        return s_self.scalar(s_x)


    def scalar_hl(self, x, t=None):
        """
        Returns the standard Hall-Littlewood scalar product of self and x.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: sp = a.scalar_hl(a); sp
            (-t^2 - 1)/(t^5 - 2*t^4 + t^3 - t^2 + 2*t - 1)
            sage: sp.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        R = self.parent().base_ring()
        p = sfa.SFAPower(R)

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

    def is_schur_positive(self):
        """
        Returns True if and only if self is Schur positive. If s is the space
        of Schur functions over self's base ring, then this is the same
        as self._is_positive(s).

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: a.is_schur_positive()
            True
            sage: a = s([2,1]) - s([3])
            sage: a.is_schur_positive()
            False

            sage: QQx = QQ['x']
            sage: s = SFASchur(QQx)
            sage: x = QQx.gen()
            sage: a = (1+x)*s([2,1])
            sage: a.is_schur_positive()
            True
            sage: a = (1-x)*s([2,1])
            sage: a.is_schur_positive()
            False
        """
        return self._is_positive( sfa.SFASchur(self.parent().base_ring()) )


    def _is_positive(self, s):
        """
        Returns True if and only if self has positive coefficients in the
        basis s.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: a._is_positive(s)
            True
            sage: a = s([2,1]) - s([3])
            sage: a._is_positive(s)
            False
        """
        s_self = s(self)
        def positive_coefficients(x):
            if is_Polynomial(x) or is_MPolynomial(x):
                return all([ c > 0 for c in x.coeffs() ])
            else:
                return x > 0

        return all([ positive_coefficients(c) for c in s_self.coefficients()])
