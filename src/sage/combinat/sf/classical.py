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


import operator
import hall_littlewood
import sfa
import llt
import macdonald
import jack
import orthotriang
import kschur

ZZ = IntegerRing()
QQ = RationalField()

translate = {'monomial':'MONOMIAL', 'homogeneous':'HOMSYM', 'power':'POWSYM', 'elementary':'ELMSYM', 'schur':'SCHUR'}

conversion_functions = {}

def init():
    """
    Set up the conversion functions between the classical bases.

    EXAMPLES:
        sage: from sage.combinat.sf.classical import init
        sage: sage.combinat.sf.classical.conversion_functions = {}
        sage: init()
        sage: sage.combinat.sf.classical.conversion_functions[('schur', 'power')]
        <built-in function t_SCHUR_POWSYM_symmetrica>
    """
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
        """
        TESTS:
            sage: from sage.combinat.sf.classical import SymmetricFunctionAlgebra_classical
            sage: s = SFASchur(QQ)
            sage: isinstance(s, SymmetricFunctionAlgebra_classical)
            True
        """
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        try:
            z = R(Integer(1))
        except:
            raise ValueError, "R must have a unit element"

        self._basis = basis
        self._element_class = element_class
        self._prefix = prefix
        self._combinatorial_class = sage.combinat.partition.Partitions_all()
        self._one = sage.combinat.partition.Partition_class([])
        CombinatorialAlgebra.__init__(self, R)



    def __call__(self, x):
        """
        Coerce x into self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s(2)
            2*s[]
            sage: s([2,1])
            s[2, 1]

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

            if R == QQ and xP.base_ring() == QQ:
                if xm:
                    return self._from_dict(t(xm)._monomial_coefficients, coerce=True)
                else:
                    return self(0)
            else:
                f = lambda part: self._from_dict(t( {part:Integer(1)} )._monomial_coefficients)
                return self._apply_module_endomorphism(x, f)


        ###############################
        # Hall-Littlewood Polynomials #
        ###############################
        elif isinstance(x, hall_littlewood.HallLittlewoodElement_generic):
            #
            #Qp: Convert to Schur basis and then convert to self
            #
            if isinstance(x, hall_littlewood.HallLittlewoodElement_qp):
                Qp = x.parent()
                sx = Qp._s._from_cache(x, Qp._s_cache, Qp._self_to_s_cache, t=Qp.t)
                return self(sx)
            #
            #P: Convert to Schur basis and then convert to self
            #
            elif isinstance(x, hall_littlewood.HallLittlewoodElement_p):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, t=P.t)
                return self(sx)
            #
            #Q: Convert to P basis and then convert to self
            #
            elif isinstance(x, hall_littlewood.HallLittlewoodElement_q):
                    return self( x.parent()._P( x ) )

        #######
        # LLT #
        #######
        #Convert to m and then to self.
        elif isinstance(x, llt.LLTElement_generic):
            P = x.parent()
            BR = self.base_ring()
            zero = BR(0)
            PBR = P.base_ring()
            if not BR.has_coerce_map_from(PBR):
                raise TypeError, "no coerce map from x's parent's base ring (= %s) to self's base ring (= %s)"%(PBR, self.base_ring())

            z_elt = {}
            for m, c in x._monomial_coefficients.iteritems():
                n = sum(m)
                P._m_cache(n)
                for part in P._self_to_m_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + BR(c*P._self_to_m_cache[n][m][part].subs(t=P.t))

            m = sfa.SFAMonomial(BR)
            return self( m._from_dict(z_elt) )

        #########################
        # Macdonald Polynomials #
        #########################
        elif isinstance(x, macdonald.MacdonaldPolynomial_generic):
            if isinstance(x, macdonald.MacdonaldPolynomial_j):
                J = x.parent()
                sx = J._s._from_cache(x, J._s_cache, J._j_to_s_cache, q=J.q, t=J.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomial_q, macdonald.MacdonaldPolynomial_p)):
                J = x.parent()._J
                jx = J(x)
                sx = J._s._from_cache(jx, J._s_cache, J._j_to_s_cache, q=J.q, t=J.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomial_h,macdonald.MacdonaldPolynomial_ht)):
                H = x.parent()
                sx = H._s._from_cache(x, H._s_cache, H._self_to_s_cache, q=H.q, t=H.t)
                return self(sx)
            elif isinstance(x, macdonald.MacdonaldPolynomial_s):
                S = x.parent()
                sx = S._s._from_cache(x, S._s_cache, S._self_to_s_cache, q=S.q, t=S.t)
                return self(sx)
            else:
                raise TypeError

        ####################
        # Jack Polynomials #
        ####################
        elif isinstance(x, jack.JackPolynomial_generic):
            if isinstance(x, jack.JackPolynomial_p):
                P = x.parent()
                mx = P._m._from_cache(x, P._m_cache, P._self_to_m_cache, t=P.t)
                return self(mx)
            if isinstance(x, (jack.JackPolynomial_j, jack.JackPolynomial_q)):
                return self( x.parent()._P(x) )
            else:
                raise TypeError

        #####################
        # k-Schur Functions #
        #####################
        if isinstance(x, kschur.kSchurFunction_generic):
            if isinstance(x, kschur.kSchurFunction_t):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, t=P.t)
                return self(sx)
            else:
                raise TypeError

        ####################################################
        # Bases defined by orthogonality and triangularity #
        ####################################################
        elif isinstance(x, orthotriang.SymmetricFunctionAlgebraElement_orthotriang):
            #Convert to its base and then to self
            xp = x.parent()
            if self is xp._sf_base:
                return xp._sf_base._from_cache(x, xp._base_cache, xp._self_to_base_cache)
            else:
                return self( xp._sf_base(x) )

        ###################
        # Skew Partitions #
        ###################
        elif x in sage.combinat.skew_partition.SkewPartitions():
            skewschur = symmetrica.part_part_skewschur(x[0], x[1])
            return self(skewschur)


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
                return eclass(self, {sage.combinat.partition.Partition_class([]):R(x)})
            except:
                raise TypeError, "do not know how to make x (= %s) an element of self"%(x)


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


    def __repr__(self):
        """
        Text representation of this symmetric function algebra.

        EXAMPLES:
            sage: SFASchur(QQ).__repr__()
            'Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis'

        """
        return "Symmetric Function Algebra over %s, %s symmetric functions as basis"%(self.base_ring(), self._basis.capitalize())


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


class SymmetricFunctionAlgebraElement_classical(sfa.SymmetricFunctionAlgebraElement_generic):
    """
    A symmetric function.
    """


##     def __pow__(self, n):
##         """
##         EXAMPLES:

##         """
##         if not isinstance(n, (int, Integer)):
##             raise TypeError, "n must be an integer"
##         A = self.parent()
##         z = A(Integer(1))
##         for i in range(n):
##             z *= self
##         return z



