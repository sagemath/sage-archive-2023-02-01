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

from sage.rings.integer import Integer

import sage.combinat.skew_partition
import sage.libs.symmetrica.all as symmetrica

from sage.rings.integer_ring import IntegerRing
from sage.rings.rational_field import RationalField



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

    EXAMPLES::

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
    """
    TESTS::

        sage: TestSuite(SymmetricFunctions(QQ).s()).run()
    """

    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        """
        A symmetric function.
        """
        pass

    def _element_constructor_(self, x):
        """
        Convert x into self, if coercion failed

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s(2)
            2*s[]
            sage: s([2,1]) # indirect doctest
            s[2, 1]
            sage: s([[2,1],[1]])
            s[1, 1] + s[2]
            sage: s([[],[]])
            s[]

            sage: McdJ = MacdonaldPolynomialsJ(QQ)
            sage: s = SymmetricFunctions(McdJ.base_ring()).s()
            sage: s._element_constructor(McdJ(s[2,1]))
            s[2, 1]
        """
        R = self.base_ring()

        eclass = self.element_class
        if isinstance(x, int):
            x = Integer(x)


        ##############
        # Partitions #
        ##############
        if x in sage.combinat.partition.Partitions_all():
            return eclass(self, {sage.combinat.partition.Partition_class(filter(lambda x: x!=0, x)):R(1)})

        # Todo: discard all of this which is taken care by Sage's coercion
        # (up to changes of base ring)

        ##############
        # Dual bases #
        ##############
        elif sfa.is_SymmetricFunction(x) and hasattr(x, 'dual'):
            #Check to see if it is the dual of some other basis
            #If it is, try to coerce its corresponding element
            #in the other basis
            return self(x.dual())

        ##################################################################
        # Symmetric Functions, same basis, possibly different coeff ring #
        ##################################################################

        # self.Element is used below to test if another symmetric
        # function is expressed in the same basis but in another
        # ground ring.  This idiom is fragile and depends on the
        # internal (unstable) specifications of parents and categories
        #
        # TODO: find the right idiom
        #
        # One cannot use anymore self.element_class: it is build by
        # the category mecanism, and depends on the coeff ring.

        elif isinstance(x, self.Element):
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
        elif isinstance(x, SymmetricFunctionAlgebra_classical.Element):


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
        elif isinstance(x, hall_littlewood.HallLittlewood_generic.Element):
            #
            #Qp: Convert to Schur basis and then convert to self
            #
            if isinstance(x, hall_littlewood.HallLittlewood_qp.Element):
                Qp = x.parent()
                sx = Qp._s._from_cache(x, Qp._s_cache, Qp._self_to_s_cache, t=Qp.t)
                return self(sx)
            #
            #P: Convert to Schur basis and then convert to self
            #
            elif isinstance(x, hall_littlewood.HallLittlewood_p.Element):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, t=P.t)
                return self(sx)
            #
            #Q: Convert to P basis and then convert to self
            #
            elif isinstance(x, hall_littlewood.HallLittlewood_q.Element):
                    return self( x.parent()._P( x ) )

        #######
        # LLT #
        #######
        #Convert to m and then to self.
        elif isinstance(x, llt.LLT_generic.Element):
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
        elif isinstance(x, macdonald.MacdonaldPolynomials_generic.Element):
            if isinstance(x, macdonald.MacdonaldPolynomials_j.Element):
                J = x.parent()
                sx = J._s._from_cache(x, J._s_cache, J._self_to_s_cache, q=J.q, t=J.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomials_q.Element, macdonald.MacdonaldPolynomials_p.Element)):
                J = x.parent()._J
                jx = J(x)
                sx = J._s._from_cache(jx, J._s_cache, J._self_to_s_cache, q=J.q, t=J.t)
                return self(sx)
            elif isinstance(x, (macdonald.MacdonaldPolynomials_h.Element,macdonald.MacdonaldPolynomials_ht.Element)):
                H = x.parent()
                sx = H._s._from_cache(x, H._s_cache, H._self_to_s_cache, q=H.q, t=H.t)
                return self(sx)
            elif isinstance(x, macdonald.MacdonaldPolynomials_s.Element):
                S = x.parent()
                sx = S._s._from_cache(x, S._s_cache, S._self_to_s_cache, q=S.q, t=S.t)
                return self(sx)
            else:
                raise TypeError

        ####################
        # Jack Polynomials #
        ####################
        elif isinstance(x, jack.JackPolynomials_generic.Element):
            if isinstance(x, jack.JackPolynomials_p.Element):
                P = x.parent()
                mx = P._m._from_cache(x, P._m_cache, P._self_to_m_cache, t=P.t)
                return self(mx)
            if isinstance(x, (jack.JackPolynomials_j.Element, jack.JackPolynomials_q.Element)):
                return self( x.parent()._P(x) )
            else:
                raise TypeError

        #####################
        # k-Schur Functions #
        #####################
        if isinstance(x, kschur.kSchurFunctions_generic.Element):
            if isinstance(x, kschur.kSchurFunctions_t.Element):
                P = x.parent()
                sx = P._s._from_cache(x, P._s_cache, P._self_to_s_cache, t=P.t)
                return self(sx)
            else:
                raise TypeError

        ####################################################
        # Bases defined by orthogonality and triangularity #
        ####################################################
        elif isinstance(x, orthotriang.SymmetricFunctionAlgebra_orthotriang.Element):
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
            import sage.libs.lrcalc.lrcalc as lrcalc
            skewschur = lrcalc.skew(x[0], x[1])
            return self._from_dict(skewschur)


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


    def is_field(self, proof = True):
        """
        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s.is_field()
            False
        """
        return False

    def is_commutative(self):
        """
        Return True if this symmetric function algebra is commutative.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: s.is_commutative()
            True
        """
        return self.base_ring().is_commutative()


    def _repr_(self):
        """
        Text representation of this symmetric function algebra.

        EXAMPLES::

            sage: SFASchur(QQ)._repr_()
            'Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis'
        """
        return "Symmetric Function Algebra over %s, %s symmetric functions as basis"%(self.base_ring(), self._basis.capitalize())




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



