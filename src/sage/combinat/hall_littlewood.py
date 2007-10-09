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

from combinat import CombinatorialClass
from sage.libs.symmetrica.all import hall_littlewood
from combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
import sfa
import partition
from sage.matrix.all import matrix, MatrixSpace
from sage.rings.all import ZZ, QQ


##################################
#Still under major development!!!#
##################################

class HallLittlewood_generic(CombinatorialAlgebra):
    pass


###########
# P basis #
###########
p_to_m_cache = {}
m_to_p_cache = {}

class HallLittlewoodElement_p(CombinatorialAlgebraElement):
    pass

class HallLittlewood_p(HallLittlewood_generic):
    def __init__(self, R):
        Rt = R['t'].fraction_field()
        self.t = Rt.gen()
        self._combinatorial_class = partition.Partitions()
        self._name = "Hall-Littlewood polynomials in the P basis"
        self._one = partition.Partition([])
        self._prefix = "P"
        self._element_class = HallLittlewoodElement_p

        self._m = sfa.SFAMonomial(Rt)

        self._p_to_m_cache = p_to_m_cache
        self._m_to_p_cache = m_to_p_cache

        CombinatorialAlgebra.__init__(self, Rt)

    def _multiply(self, left, right):
        """
        Convert to the monomials, do the multiplication there, and
        convert back to the P basis.
        """
        return self( self._m(left) * self._m(right) )


    def _coerce_start(self, x):
        """
        Coerce things into the Hall-Littlewood P basis.

        EXAMPLES:
          1) Hall-Littlewood polynomials in the Q basis

          2) Monomial symmetric functions
        """

        if isinstance(x, HallLittlewoodElement_q):
            BR = self.base_ring()
            t = self.t
            z_elt = {}

            #(1-Sym::vHL)^nops(part)*_mult(_mult((1-(Sym::vHL)^j)/(1-Sym::vHL) $j=1..i)
            # $i in combinat::partitions::toExp(part))*HallLittlewood::P(part)
            for m, c in x._monomial_coefficients.iteritems():
                coeff = (1-t)**len(m)
                for i in m.to_exp():
                    for j in range(1,i+1):
                        coeff *= (1-t**j)/(1-t)
                z_elt[m] = BR( c*coeff )

            z = self(0)
            z._monomial_coefficients = z_elt
            return z
        else:
            raise TypeError

    def _scalar_hl_part(part1, part2):
        return self._m(part1).scalar_hl(part2, t=self.t)

    def _m_cache(self, n):
        if n in p_to_m_cache:
            return

        #All the coefficients stored will be in
        #Zt even though we need to do the computations
        #in Ztff
        Qt = QQ['t']
        Qtff = Qt.fraction_field()
        t = Qtff.gen()
        m = sfa.SFAMonomial(QQ)
        one = Qt(1)
        zero = Qtff(0)

        pn = partition.Partitions(n)
        len_pn = pn.count()
        pn = pn.list()
        mpn = map(m, pn)

        p2m_n = {}
        m2p_n = {}

        #Store the list of scalar products <P_mu, P_lambda>_t
        p_scalars = [None]*len_pn

        ################################
        #Compute the expansions of HL  #
        #from [1,...,1] to [n]         #
        ################################

        #The coefficient of the all ones partition is 1
        p2m_n[pn[-1]] = {pn[-1]: one}
        p_scalars[-1] = mpn[-1].scalar_hl(mpn[-1], t=t)


        M = sfa.SFAMonomial(Qt)
        def convert_to_m( i):
            res = M(0)
            res._monomial_coefficients = p2m_n[ pn[i] ]
            return res

        for i in range(len_pn-2,-1,-1):
            p2m_n[pn[i]] = {}
            mi = mpn[ i ]

            #Coefficient on mu = lambda is 1
            p2m_n[ pn[i] ][ pn[i] ] = one

            for j in range(i+1, len_pn):
                mj = mpn[ j ]
                #Calculate <m_i, P_j>_t and store it in value
                Pj = convert_to_m(j)
                value = M(pn[i]).scalar_hl(Pj, t=t)
                p2m_n[ pn[i] ][ pn[j] ]  = (-value/p_scalars[j])
                #p2m_n[ pn[i] ][ pn[j] ] = p2m_n[ pn[i] ][ pn[j] ]

            Pi = convert_to_m(i)
            p_scalars[i] = Pi.scalar_hl(Pi, t=t)

        p_to_m_cache[n] = p2m_n

        print p_scalars



###########
# Q basis #
###########
class HallLittlewoodElement_q(CombinatorialAlgebraElement):
    pass

class HallLittlewood_q(HallLittlewood_generic):
    def __init__(self, R):
        Rt = R['t'].fraction_field()
        self.t = Rt.gen()
        self._combinatorial_class = partition.Partitions()
        self._name = "Hall-Littlewood polynomials in the Q basis"
        self._one = partition.Partition([])
        self._prefix = "Q"
        self._element_class = HallLittlewoodElement_q

        self._P = HallLittlewood_p(R)

        CombinatorialAlgebra.__init__(self,Rt)

    def _multiply(self, left, right):
        """
        Converts to the P basis, does the multiplication there, and
        converts back to the Q basis.

        """
        return self( self._P(left) * self._P(right) )


    def _coerce_start(self, x):
        """
        Coerce things into the Hall-Littlewood Q basis.

        1) Hall-Littlewood polynomials in the P basis
        """

        if isinstance(x, HallLittlewoodElement_p):
            BR = self.base_ring()
            t = self.t
            z_elt = {}

            #(1/(1-Sym::vHL)^nops(part))*_mult(_mult((1-Sym::vHL)/(1-Sym::vHL^j) $j=1..i)
            # $i in combinat::partitions::toExp(part))*HallLittlewood::Q(part)
            for m, c in x._monomial_coefficients.iteritems():
                coeff = 1/(1-t)**len(m)
                for i in m.to_exp():
                    for j in range(1,i+1):
                        coeff *= (1-t)/(1-t**j)
                z_elt[m] = BR( c*coeff )

            z = self(0)
            z._monomial_coefficients = z_elt
            return z
        else:
            raise TypeError


############
# Qp basis #
############
qp_to_s_cache = {}
s_to_qp_cache = {}

class HallLittlewoodElement_qp(CombinatorialAlgebraElement):
    def expand(self, n, alphabet='x'):
        return self.parent()._s(self).expand(n, alphabet=alphabet)

class HallLittlewood_qp(HallLittlewood_generic):
    def __init__(self, R):
        Rt = R['t'].fraction_field()
        self.t = Rt.gen()
        self._combinatorial_class = partition.Partitions()
        self._name = "Hall-Littlewood polynomials in the Qp basis"
        self._one = partition.Partition([])
        self._prefix = "Qp"
        self._element_class = HallLittlewoodElement_qp

        self._s = sfa.SFASchur(Rt)

        self._qp_to_s_cache = qp_to_s_cache
        self._s_to_qp_cache = s_to_qp_cache

        CombinatorialAlgebra.__init__(self,Rt)


    def _coerce_start(self, x):
        """
        Coerce things into the Hall-Littlewood Qp basis.

        1) Symmetric Functions
        """
        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            sx = self._s( x )
            BR = self.base_ring()
            zero = BR(0)

            z_elt = {}
            for m, c in sx._monomial_coefficients.iteritems():
                n = sum(m)
                self._s_cache(n)
                for part in self._s_to_qp_cache[n][m]:
                    z_elt[part] = z_elt.get(part, zero) + BR(c*self._s_to_qp_cache[n][m][part])

            z = self(0)
            z._monomial_coefficients = z_elt
            return z

        raise TypeError

    def _multiply(self, left, right):
        """
        Converts the Hall-Littlewood polynomial in the Qp basis
        to a Schur function, performs the multiplication there,
        and converts it back to the Qp basis.
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
        #the empty partition.
        if n == 0:
            p = partition.Partition_class([])
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
        pn = partition.Partitions_n(n).list()

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
