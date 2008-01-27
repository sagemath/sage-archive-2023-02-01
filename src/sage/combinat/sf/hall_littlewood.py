"""
Hall-Littlewood Polynomials
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

from sage.combinat.combinat import CombinatorialClass
from sage.libs.symmetrica.all import hall_littlewood
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
import sfa
import sage.combinat.partition
import kfpoly
from sage.matrix.all import matrix, MatrixSpace
from sage.rings.all import ZZ, QQ


##################################
#Still under major development!!!#
##################################

def HallLittlewoodP(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood
    P basis.  This is the same as the HL basis in John Stembridge's
    SF examples file.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: HallLittlewoodP(QQ)
        Hall-Littlewood polynomials in the P basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodP(QQ, t=-1)
        Hall-Littlewood polynomials in the P basis with t=-1 over Rational Field
        sage: HLP = HallLittlewoodP(QQ)
        sage: s = SFASchur(HLP.base_ring())
        sage: s(HLP([2,1]))
        (-t^2-t)*s[1, 1, 1] + s[2, 1]


    """
    return cache_p(R,t)

def HallLittlewoodQ(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood
    Q basis.  This is the same as the Q basis in John Stembridge's
    SF examples file.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: HallLittlewoodQ(QQ)
        Hall-Littlewood polynomials in the Q basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodQ(QQ, t=-1)
        Hall-Littlewood polynomials in the Q basis with t=-1 over Rational Field

    """
    return cache_q(R,t)
def HallLittlewoodQp(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood
    Qp basis.  This is dual to the Hall-Littlewood P basis with
    respect to the standard scalar product.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: HallLittlewoodQp(QQ)
        Hall-Littlewood polynomials in the Qp basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodQp(QQ, t=-1)
        Hall-Littlewood polynomials in the Qp basis with t=-1 over Rational Field

    """
    return cache_qp(R,t)


##################################

class HallLittlewood_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, t=None):
        if t is None:
            R = R['t'].fraction_field()
            self.t = R.gen()
        elif t not in R:
            raise ValueError, "t (=%s) must be in R (=%s)"%(t,R)
        else:
            self.t = R(t)
            self._name += " with t=%s"%self.t

        self._combinatorial_class = sage.combinat.partition.Partitions()
        self._one = sage.combinat.partition.Partition([])

        CombinatorialAlgebra.__init__(self, R)


    def transition_matrix(self, basis, n):
        """
        Returns the transitions matrix between self and basis
        for the homogenous component of degree n.

        EXAMPLES:
            sage: HLP = HallLittlewoodP(QQ)
            sage: s   = SFASchur(HLP.base_ring())
            sage: HLP.transition_matrix(s, 4)
            [             1             -t              0            t^2           -t^3]
            [             0              1             -t             -t      t^3 + t^2]
            [             0              0              1             -t            t^3]
            [             0              0              0              1 -t^3 - t^2 - t]
            [             0              0              0              0              1]
            sage: HLQ = HallLittlewoodQ(QQ)
            sage: HLQ.transition_matrix(s,3)
            [                        -t + 1                        t^2 - t                     -t^3 + t^2]
            [                             0                  t^2 - 2*t + 1           -t^4 + t^3 + t^2 - t]
            [                             0                              0 -t^6 + t^5 + t^4 - t^2 - t + 1]
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLQp.transition_matrix(s,3)
            [      1       0       0]
            [      t       1       0]
            [    t^3 t^2 + t       1]


        """
        P = sage.combinat.partition.Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( map( lambda col_part: z.coefficient(col_part), Plist ) )
        return matrix(m)

class HallLittlewoodElement_generic(sfa.SymmetricFunctionAlgebraElement_generic):
    def frobenius(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLP([2,1]).frobenius()
            (-t^5-t^4+t^2+t)*P[1, 1, 1] + (-t^3-t^2+1)*P[2, 1] + (-t^2-t)*P[3]
            sage: HLP([2]).frobenius()
            (-t^2+1)*P[1, 1] - t*P[2]
            sage: HLQp([2]).frobenius()
            Qp[1, 1] - t*Qp[2]
        """
        sp = self.parent()
        BR = sp.base_ring()
        s = sfa.SFASchur(BR)
        return sp( s(self).frobenius() )

    def omega(self):
        """
        An alias for self.frobenius() .
        """
        return self.frobenius()

    def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n variables.

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLP([2]).expand(2)
            x0^2 + (-t + 1)*x0*x1 + x1^2
            sage: HLQ([2]).expand(2)
            (-t + 1)*x0^2 + (t^2 - 2*t + 1)*x0*x1 + (-t + 1)*x1^2
            sage: HLQp([2]).expand(2)
            x0^2 + x0*x1 + x1^2
        """
        sp = self.parent()
        BR = sp.base_ring()
        s = sfa.SFASchur(BR)
        return s(self).expand(n, alphabet=alphabet)

    def scalar(self, x):
        """
        Returns standard scalar product between self and s.

        This is the default implementation that converts both self
        and x into Schur functions and performs the scalar product
        that basis.

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLP([2]).scalar(HLQp([2]))
            1
            sage: HLP([2]).scalar(HLQp([1,1]))
            0
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
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLP([2]).scalar_hl(HLQ([2]))
            1
            sage: HLP([2]).scalar_hl(HLQ([1,1]))
            0
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
            t = self.parent().t

        smcs = smaller._monomial_coefficients
        gmcs = greater._monomial_coefficients
        for s_part in smcs :
            if s_part in gmcs:
                res += smcs[s_part]*gmcs[s_part]*s_part.centralizer_size(t=t)

        return res
###########
# P basis #
###########
p_to_m_cache = {}
m_to_p_cache = {}
p_to_s_cache = {}
s_to_p_cache = {}

class HallLittlewoodElement_p(HallLittlewoodElement_generic):
    pass

class HallLittlewood_p(HallLittlewood_generic):
    def __init__(self, R, t=None):
        self._name = "Hall-Littlewood polynomials in the P basis"
        self._prefix = "P"
        self._element_class = HallLittlewoodElement_p

        HallLittlewood_generic.__init__(self, R, t=t)

        self._m = sfa.SFAMonomial(self.base_ring())
        self._s = sfa.SFASchur(self.base_ring())

        self._p_to_m_cache = p_to_m_cache
        self._m_to_p_cache = m_to_p_cache

        self._p_to_s_cache = p_to_s_cache
        self._s_to_p_cache = s_to_p_cache



    def _multiply(self, left, right):
        """
        Convert to the Schur basis, do the multiplication there, and
        convert back to the P basis.

        EXAMPLES:
            sage: HLP = HallLittlewoodP(QQ)
            sage: HLP([2])^2
            (t+1)*P[2, 2] + (-t+1)*P[3, 1] + P[4]

        """
        return self( self._s(left) * self._s(right) )


    def _coerce_start(self, x):
        """
        Coerce things into the Hall-Littlewood P basis.

        1) Hall-Littlewood polynomials in the Q basis
        2) Hall-Littlewood polynomials in the Qp basis (via the Schurs)
        3) Classical symmetric functions

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLP(HLQ([2]))
            (-t+1)*P[2]
            sage: HLP(HLQp([2]))
            t*P[1, 1] + P[2]
            sage: HLP(s([2]))
            t*P[1, 1] + P[2]
            sage: HLP(p([2]))
            (t-1)*P[1, 1] + P[2]

        """
        BR = self.base_ring()
        if isinstance(x, HallLittlewoodElement_q):
            t = self.t
            def f(m):
                coeff = (1-t)**len(m)
                for i in m.to_exp():
                    for j in range(1,i+1):
                        coeff *= (1-t**j)/(1-t)
                return coeff
            return self._change_by_proportionality(x, f)
        elif isinstance(x, HallLittlewoodElement_qp):
            return self( self._s(x) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            #Convert x to the Schur basis
            x = self._s(x)
            zero = self.base_ring()(0)
            z_elt = {}
            for part,c in x.monomial_coefficients().iteritems():
                if sum(part) not in self._s_to_p_cache:
                    self._s_cache(sum(part))
                for part2, c2 in self._s_to_p_cache[sum(part)][part].iteritems():
                    z_elt[ part2 ] = z_elt.get(part2, zero) + BR(c*c2.subs(t=self.t))
            res = self(0)
            res._monomial_coefficients = z_elt
            return res

        else:
            raise TypeError

    def _get_hl_matrix(self, n, pn):
        #Univariate polynomial arithmetic is faster
        #over ZZ.  Since that is all we need to compute
        #the transition matrices between S and Qp, we
        #should use that.
        Zt = ZZ['t']
        t = Zt.gen()
        zero = Zt(0)

        #Generate the expansions of the Qp polynomials in terms
        #of the Schur functions.
        hlpn = [ kfpoly.schur_to_hl(p, t) for p in pn ]

        #Since the coefficients returned by hall_littlewood are in ZZ['x'], we
        #need to replace the x's with t's.
        hlpn_m = [[ x.get(p, zero) for p in pn ] for x in hlpn ]

        return hlpn_m

    def _s_cache(self, n):
        """
        Computes the change of basis between the P polynomials and
        the Schur functions for partitions of size n.

        Should use the fact that the transformation matrix is lower-triangular
        in order to obtain the inverse transformation.

        """
        global p_to_s_cache, s_to_p_cache

        #Do nothing if we've already computed the transition matrices
        #for degree n.
        if n in p_to_s_cache:
            return

        #Univariate polynomial arithmetic is faster
        #over ZZ.  Since that is all we need to compute
        #the transition matrices between S and P, we
        #should use that.
        Zt = ZZ['t']
        t = Zt.gen()

        #We have to handle the case where n == 0 seperately
        #since symmetrica.hall_littlewood does not like
        #the empty sage.combinat.partition.
        if n == 0:
            p = sage.combinat.partition.Partition_class([])
            p_to_s_cache[ n ] = {p: {p: Zt(1)}}
            s_to_p_cache[ n ] = {p: {p: Zt(1)}}
            return


        #Make sure we don't need to spend extra time
        #coercing things into Zt
        one = Zt(1)
        zero = Zt(0)
        def delta(i):
            def f(j):
                if i == j:
                    return one
                else:
                    return zero
            return f

        #Get and store the list of partition we'll need
        pn = sage.combinat.partition.Partitions_n(n).list()

        #Get the matrix with all the coefficients of the
        #expansions of the P polynomials in S.
        #Note: The function should probably be written so
        #that we don't need to store all of these in memory
        #at once
        hlpn_m = self._get_hl_matrix(n, pn)


        #Create the initial cache dictionaries
        p2s_n = {}
        s2p_n = {}
        for i in range(len(pn)):
            p2s_part = {}
            s2p_part = {}
            #Since we already have to coefficients from
            #S -> P, we can store them here.
            for j in range(i, len(pn)):
                if hlpn_m[i][j] != zero:
                    s2p_part[ pn[ j ] ] = hlpn_m[i][j]
            p2s_n[ pn[i] ] = p2s_part
            s2p_n[ pn[i] ] = s2p_part


##         #Compute the inverse of hlpn_m by using back
##         #substitution.  We solve a len(pn) systems of
##         #equations hlpn_m*x = b_i for x, where e_i
##         #is the ith standard basis vector
##         for column in range(len(pn)):
##             e = delta(column)
##             x = []
##             for i in range(len(pn)):
##                 value = e(i)
##                 for j in range(len(x)):
##                     value -= hlpn_m[i][j]*x[j]
##                 x.append(value)
##             for j in range(column,len(x)):
##                 if x[j] != zero:
##                     s2p_n[ pn[j] ][ pn[column] ] = x[ j ]


        #Just compute the matrix inverse for now
        inverse = ~matrix(hlpn_m)
        for i in range(len(pn)):
            for j in range(len(pn)):
                if inverse[i,j] != zero:
                    p2s_n[ pn[i] ] [ pn[j] ] = inverse[i,j]

        p_to_s_cache[ n ] = p2s_n
        s_to_p_cache[ n ] = s2p_n



##     def _scalar_hl_part(part1, part2):
##         return self._s(part1).scalar_hl(part2, t=self.t)

##     def _m_cache(self, n):
##         if n in p_to_m_cache:
##             return

##         #All the coefficients stored will be in
##         #Zt even though we need to do the computations
##         #in Ztff
##         Qt = QQ['t']
##         Qtff = Qt.fraction_field()
##         t = Qtff.gen()
##         m = sfa.SFAMonomial(QQ)
##         one = Qt(1)
##         zero = Qtff(0)

##         pn = sage.combinat.partition.Partitions(n)
##         len_pn = pn.count()
##         pn = pn.list()
##         mpn = map(m, pn)

##         p2m_n = {}
##         m2p_n = {}

##         #Store the list of scalar products <P_mu, P_lambda>_t
##         p_scalars = [None]*len_pn

##         ################################
##         #Compute the expansions of HL  #
##         #from [1,...,1] to [n]         #
##         ################################

##         #The coefficient of the all ones partition is 1
##         p2m_n[pn[-1]] = {pn[-1]: one}
##         p_scalars[-1] = mpn[-1].scalar_hl(mpn[-1], t=t)


##         M = sfa.SFAMonomial(Qt)
##         def convert_to_m(i):
##             res = M(0)
##             res._monomial_coefficients = p2m_n[ pn[i] ]
##             return res

##         for i in range(len_pn-2,-1,-1):
##             p2m_n[pn[i]] = {}
##             mi = mpn[ i ]

##             #Coefficient on mu = lambda is 1
##             p2m_n[ pn[i] ][ pn[i] ] = one

##             for j in range(i+1, len_pn):
##                 mj = mpn[ j ]
##                 #Calculate <m_i, P_j>_t and store it in value
##                 Pj = convert_to_m(j)
##                 value = M(pn[i]).scalar_hl(Pj, t=t)
##                 p2m_n[ pn[i] ][ pn[j] ]  = (-value/p_scalars[j])
##                 #p2m_n[ pn[i] ][ pn[j] ] = p2m_n[ pn[i] ][ pn[j] ]

##             Pi = convert_to_m(i)
##             p_scalars[i] = Pi.scalar_hl(Pi, t=t)

##         p_to_m_cache[n] = p2m_n

##         print p_scalars



###########
# Q basis #
###########
class HallLittlewoodElement_q(HallLittlewoodElement_generic):
    pass

class HallLittlewood_q(HallLittlewood_generic):
    def __init__(self, R, t=None):
        self._name = "Hall-Littlewood polynomials in the Q basis"
        self._prefix = "Q"
        self._element_class = HallLittlewoodElement_q

        HallLittlewood_generic.__init__(self, R, t=t)

        self._P = HallLittlewood_p(R, t=t)



    def _multiply(self, left, right):
        """
        Converts to the P basis, does the multiplication there, and
        converts back to the Q basis.

        EXAMPLES:
            sage: HLQ = HallLittlewoodQ(QQ)
            sage: HLQ([2])^2
            Q[2, 2] + (-t+1)*Q[3, 1] + (-t+1)*Q[4]

        """
        return self( self._P(left) * self._P(right) )


    def _coerce_start(self, x):
        """
        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLQ( HLP([2,1]) + HLP([3]) )
            (1/(t^2-2*t+1))*Q[2, 1] + (1/(-t+1))*Q[3]
            sage: HLQ(HLQp([2]))
            (-t/(-t^3+t^2+t-1))*Q[1, 1] + (1/(-t+1))*Q[2]
            sage: HLQ(s([2]))
            (-t/(-t^3+t^2+t-1))*Q[1, 1] + (1/(-t+1))*Q[2]
            sage: HLQ(p([2]))
            (-1/(-t^2+1))*Q[1, 1] + (1/(-t+1))*Q[2]
        """
        if isinstance(x, HallLittlewoodElement_p):
            t = self.t
            def f(m):
                coeff = 1/(1-t)**len(m)
                for i in m.to_exp():
                    for j in range(1,i+1):
                        coeff *= (1-t)/(1-t**j)
                return coeff
            return self._change_by_proportionality(x, f)
        elif isinstance(x, HallLittlewoodElement_qp):
            return self( self._P(x) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            return self( self._P(x) )
        else:
            raise TypeError


############
# Qp basis #
############
qp_to_s_cache = {}
s_to_qp_cache = {}

class HallLittlewoodElement_qp(HallLittlewoodElement_generic):
    pass


class HallLittlewood_qp(HallLittlewood_generic):
    def __init__(self, R, t=None):
        self._name = "Hall-Littlewood polynomials in the Qp basis"
        self._prefix = "Qp"
        self._element_class = HallLittlewoodElement_qp

        HallLittlewood_generic.__init__(self,R, t=t)

        self._s = sfa.SFASchur(self.base_ring())

        self._qp_to_s_cache = qp_to_s_cache
        self._s_to_qp_cache = s_to_qp_cache




    def _coerce_start(self, x):
        """
        Coerce things into the Hall-Littlewood Qp basis.

        1) Hall-Littlewood polynomials in the Q basis (via the Schurs)
        2) Hall-Littlewood polynomials in the P basis (via the Schurs)
        3) Symmetric Functions

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLQp(HLP([2]))
            -t*Qp[1, 1] + (t^2+1)*Qp[2]
            sage: HLQp(HLQ([2]))
            (t^2-t)*Qp[1, 1] + (-t^3+t^2-t+1)*Qp[2]
            sage: HLQp(s([2]))
            Qp[2]
            sage: HLQp(p([2]))
            -Qp[1, 1] + (t+1)*Qp[2]
        """
        if isinstance(x, HallLittlewoodElement_q):
            return self( self._s( x ) )
        elif isinstance(x, HallLittlewoodElement_p):
            return self( self._s( x ) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            sx = self._s( x )
            BR = self.base_ring()
            zero = BR(0)

            z_elt = {}
            for m, c in sx._monomial_coefficients.iteritems():
                n = sum(m)
                self._s_cache(n)
                for part in self._s_to_qp_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + BR(c*self._s_to_qp_cache[n][m][part].subs(t=self.t))

            z = self(0)
            z._monomial_coefficients = z_elt
            return z
        else:
            raise TypeError

    def _multiply(self, left, right):
        """
        Converts the Hall-Littlewood polynomial in the Qp basis
        to a Schur function, performs the multiplication there,
        and converts it back to the Qp basis.

        EXAMPLES:
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLQp([2])^2
            Qp[2, 2] + (-t+1)*Qp[3, 1] + (-t+1)*Qp[4]

        """
        return self( self._s(left)*self._s(right) )._monomial_coefficients

    def _get_hl_matrix(self, n, pn):
        #Univariate polynomial arithmetic is faster
        #over ZZ.  Since that is all we need to compute
        #the transition matrices between S and Qp, we
        #should use that.
        Zt = ZZ['t']
        t = Zt.gen()

        #Generate the expansions of the Qp polynomials in terms
        #of the Schur functions.
        hlqpn = [ hall_littlewood(p) for p in pn ]

        #Since the coefficients returned by hall_littlewood are in ZZ['x'], we
        #need to replace the x's with t's.
        hlqpn_m = [[ x.coefficient(p).subs(x=t) for p in pn ] for x in hlqpn ]

        return hlqpn_m


    def _s_cache_local(self, n):
        #Compute the generic transition matrices
        #if needed
        self._s_cache_generic(n)

        #Convert the generic transition matrices to matrices in
        #our self's base ring
        self._qp_to_s_cache[n] = {}
        self._s_to_qp_cache[n] = {}

        BR = self.base_ring()

        coerce = True
        if BR.base() == ZZ['t']:
            coerce=False

        for part1 in qp_to_s_cache[n]:
            self._qp_to_s_cache[n][part1] = {}
            self._s_to_qp_cache[n][part1] = {}
            for part2 in qp_to_s_cache[n][part1]:
                self._qp_to_s_cache[n][part1][part2] = BR(qp_to_s_cache[n][part1][part2], coerce=coerce)
            for part2 in s_to_qp_cache[n][part1]:
                self._s_to_qp_cache[n][part1][part2] = BR(s_to_qp_cache[n][part1][part2], coerce=coerce)

    def _s_cache(self, n):
        """
        Computes the change of basis between the Qp polynomials and
        the Schur functions for partitions of size n.

        Uses the fact that the transformation matrix is lower-triangular
        in order to obtain the inverse transformation.

        """
        global qp_to_s_cache, s_to_qp_cache

        #Do nothing if we've already computed the transition matrices
        #for degree n.
        if n in qp_to_s_cache:
            return

        #Univariate polynomial arithmetic is faster
        #over ZZ.  Since that is all we need to compute
        #the transition matrices between S and Qp, we
        #should use that.
        Zt = ZZ['t']
        t = Zt.gen()

        #We have to handle the case where n == 0 seperately
        #since symmetrica.hall_littlewood does not like
        #the empty sage.combinat.partition.
        if n == 0:
            p = sage.combinat.partition.Partition_class([])
            qp_to_s_cache[ n ] = {p: {p: Zt(1)}}
            s_to_qp_cache[ n ] = {p: {p: Zt(1)}}
            return


        #Make sure we don't need to spend extra time
        #coercing things into Zt
        one = Zt(1)
        zero = Zt(0)
        def delta(i):
            def f(j):
                if i == j:
                    return one
                else:
                    return zero
            return f

        #Get and store the list of partition we'll need
        pn = sage.combinat.partition.Partitions_n(n).list()

        #Get the matrix with all the coefficients of the
        #expansions of the Qp polynomials in S.
        #Note: The function should probably be written so
        #that we don't need to store all of these in memory
        #at once
        hlqpn_m = self._get_hl_matrix(n, pn)


        #Create the initial cache dictionaries
        qp2s_n = {}
        s2qp_n = {}
        for i in range(len(pn)):
            qp2s_part = {}
            s2qp_part = {}
            #Since we already have to coefficients from
            #Qp -> S, we can store them here.
            for j in range(i+1):
                if hlqpn_m[i][j] != zero:
                    qp2s_part[ pn[ j ] ] = hlqpn_m[i][j]
            qp2s_n[ pn[i] ] = qp2s_part
            s2qp_n[ pn[i] ] = s2qp_part


        #Compute the inverse of hlqpn_m by using forward
        #substitution.  We solve a len(pn) systems of
        #equations hlqpn_m*x = b_i for x, where e_i
        #is the ith standard basis vector
        for column in range(len(pn)):
            e = delta(column)
            x = []
            for i in range(len(pn)):
                value = e(i)
                for j in range(len(x)):
                    value -= hlqpn_m[i][j]*x[j]
                x.append(value)
            for j in range(column,len(x)):
                if x[j] != zero:
                    s2qp_n[ pn[j] ][ pn[column] ] = x[ j ]

        qp_to_s_cache[ n ] = qp2s_n
        s_to_qp_cache[ n ] = s2qp_n



#############
#   Cache   #
#############
from sage.misc.cache import Cache
cache_p = Cache(HallLittlewood_p)
cache_q = Cache(HallLittlewood_q)
cache_qp = Cache(HallLittlewood_qp)

