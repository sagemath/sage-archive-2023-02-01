"""
Macdonald Polynomials -- under development.
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
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
import sfa
import sage.combinat.partition
from sage.matrix.all import matrix, MatrixSpace
from sage.rings.all import ZZ, QQ
from sage.misc.misc import prod



def c1(part, q, t):
    res = 1
    for i in range(part.size()):
        res *= 1-q^(sum(part.arm_lengths(),[])[i]+1)*t^(sum(part.leg_lengths(),[])[i])
    return res

def c2(part, q, t):
    res = 1
    for i in range(part.size()):
        res *= 1-q^(sum(part.arm_lengths(),[])[i])*t^(sum(part.leg_lengths(),[])[i]+1)
    return res


#Generic MacdonaldPolynomials
class Macdonald_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R):
        Rqt = R['q','t'].fraction_field()
        (self.q, self.t) = Rqt.gens()
        self._combinatorial_class = sage.combinat.partition.Partitions()
        self._one = sage.combinat.partition.Partition([])
        CombinatorialAlgebra.__init__(self, Rqt)

class MacdonaldElement_generic(sfa.SymmetricFunctionAlgebraElement_generic):
    def scalar(self, x):
        """
        Returns the Macdonald scalar product of self and x by converting
        both to the power-sum basis.
        """
        P = self.parent()
        p = sfa.SFAPower(self.parent().base_ring())
        p_self = p(self)
        p_x = p(x)

        def f(part1, part2):
            return part1.centralizer_size(t=P.t, q=P.q)

        return P._apply_multi_module_morphism(p_self, p_x, f, orthogonal=True)


    def frobenius(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.
        """
        P = self.parent()
        p = sfa.SFAPower(self.parent().base_ring())
        p_self = p(self)
        def f(part):
            return (-1)**(part.size()-part.length())*prod([(1-P.q**part[i])/(1-P.t**part[i]) for i in part])*p(part)
        return P( P._apply_module_morphism(self, f) )


#P basis
class Macdonald_p(Macdonald_generic):
    def __init__(self, R):
        self._name = "Macdonald polynomials in the P basis"
        self._prefix = "McdP"
        self._element_class = MacdonaldElement_p
        Macdonald_generic.__init__(self, R)

    def _coerce_start(self, x):
        """
        Coerce things into the P basis

        1) Q basis (proportional)
        2) J basis (proportional)
        """
        BR = self.base_ring()

        #Macdonald Q-basis
        if isinstance(x, Macdonald_q):
            def f(m):
                return 1/self(m).scalar(self(m))
            return self._change_by_proportionality(x, f)
        if isinstance(x, Macdonald_j):
            def f(m):
                return c2(m, self._q, self._t)
            return self._change_by_proportionality(x, f)
        else:
            raise TypeError
class MacdonaldElement_p(MacdonaldElement_generic):
    def scalar(self, x):
        if isinstance(x, MacdonaldElement_p):
            P = self.parent()
            p_x = P(x)
            def f(part1, part2):
                return c1(part1, P.q, P.t)/c2(part, P.q, P.t)
            return P._apply_multi_module_morphism(self, self, p_x, f, orthogonal=True)
        elif isinstance(x, MacdonaldElement_q):
            return x.scalar(self)
        else:
            return MacdonaldElement_generic.scalar(self, x)


#Q basis
class Macdonald_q(Macdonald_generic):
    def __init__(self, R):
        self._name = "Macdonald polynomials in the Q basis"
        self._prefix = "McdQ"
        self._element_class = MacdonaldElement_q
        Macdonald_generic.__init__(self, R)

    def _coerce_start(self, x):
        """
        Coerce things into the Q basis

        1) P basis (proportional)
        2) J basis (proportional)
        """
        BR = self.base_ring()

        #Macdonald P-basis
        if isinstance(x, Macdonald_p):
            def f(m):
                return self(m).scalar(self(m))
            return self._change_by_proportionality(x, f)
        elif isinstance(x, Macdonald_j):
            def f(m):
                return self(m).scalar(self(m))*c2(m, self.q, self._t)
            return self._change_by_proportionality(x, f)
        else:
            raise TypeError

class MacdonaldElement_q(MacdonaldElement_generic):
    def scalar(self, x):
        if isinstance(x, MacdonaldElement_p):
            def f(part1, part2):
                return 1
            return P._apply_multi_module_morphism(self, self, x, f, orthogonal=True)
        else:
            return MacdonaldElement_generic.scalar(self, x)


#J basis
j_to_s_cache = {}
s_to_j_cache = {}
class Macdonald_j(Macdonald_generic):
    def __init__(self, R):
        self._name = "Macdonald polynomials in the J basis"
        self._prefix = "McdJ"
        self._element_class = MacdonaldElement_j
        Macdonald_generic.__init__(self, R)

        self._s = sfa.SFASchur(self.base_ring())
        self._j_to_s_cache = j_to_s_cache
        self._s_to_j_cache = s_to_j_cache

    def _coerce_start(self, x):
        """
        Coerce things into the J basis

        1) P basis (proportional)
        2) Q basis (proportional)
        4) Schur functions ( Lascoux, Lapointe and Morse creation operators on Macdonald polynomials )
        3) Everything else through s
        """
        BR = self.base_ring()

        #Macdonald P-basis
        if isinstance(x, Macdonald_p):
            def f(m):
                return 1/c2(m, self.q, self._t)
            return self._change_by_proportionality(x, f)
        #Macdonald Q-basis
        elif isinstance(x, Macdonald_q):
            def f(m):
                return 1/(self(m).scalar(self(m))*c2(m, self.q, self.t))
            return self._change_by_proportionality(x, f)
        #Every
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
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
        else:
            raise TypeError

    def _s_cache(self, n):
        """
        Computes the change of basis between the P polynomials and
        the Schur functions for partitions of size n.

        """
        global j_to_s_cache, s_to_j_cache

        #Do nothing if we've already computed the transition matrices
        #for degree n.
        if n in j_to_s_cache:
            return

        #Univariate polynomial arithmetic is faster
        #over ZZ.  Since that is all we need to compute
        #the transition matrices between S and P, we
        #should use that.
        QQqt = QQ['q','t']
        q,t = QQqt.gen()

        #We have to handle the case where n == 0 seperately
        #since symmetrica.hall_littlewood does not like
        #the empty sage.combinat.partition.
        if n == 0:
            p = sage.combinat.partition.Partition_class([])
            j_to_s_cache[ n ] = {p: {p: QQqt(1)}}
            s_to_j_cache[ n ] = {p: {p: QQqt(1)}}
            return


        #Make sure we don't need to spend extra time
        #coercing things into QQqt
        one = QQqt(1)
        zero = QQqt(0)
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

        j_to_s_cache[ n ] = p2s_n
        s_to_j_cache[ n ] = s2p_n

class MacdonaldElement_j(MacdonaldElement_generic):
    def scalar(self, x):
        if isinstance(x, MacdonaldElement_j):
            J = self.parent()
            j_x = J(x)
            def f(part1, part2):
                return c1(part1, J.q, J.t)*c2(part, J.q, J.t)
            return J._apply_multi_module_morphism(self, self, j_x, f, orthogonal=True)
        else:
            return MacdonaldElement_generic.scalar(self, x)



#H basis
class Macdonald_h(Macdonald_generic):
    pass
class MacdonaldElement_h(MacdonaldElement_generic):
    def nabla(self):
        """
        Returns the value of the nabla operator applied to self.  The
        eigenvectors of the nabla operator are the Macdonald polynomials.
        """
        P = self.parent()
        def f(part):
            return (P.t)**(part.weighted_size())*(P.q)**(part.conjugate().weighted_size())*P(part)
        return P._apply_module_morphism(self, f)

#Ht basis
class Macdonald_ht(Macdonald_generic):
    pass
class MacdonaldElement_ht(MacdonaldElement_generic):
    pass

#S basis
class Macdonald_s(Macdonald_generic):
    def __init__(self, R):
        self._name = "Macdonald polynomials in the S basis"
        self._prefix = "McdS"
        self._element_class = MacdonaldElement_s
        Macdonald_generic.__init__(self, R)

        self._s = sfa.SFASchur(self.base_ring())

    def _multiply(self, left, right):
        """
        The multiplication of the modified Schur functions behaves the
        same as the multiplication of the Schur functions.
        """
        s_left = self._s._from_element(left)
        s_right = self._s._from_element(right)
        product = s_left*s_right
        return self._from_element(product)


class MacdonaldElement_s(MacdonaldElement_generic):
    def frobenius(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.
        """
        S = self.parent()
        def f(part):
            return S._s(part)
        return S( S._apply_module_morphism(self, f) )

    def creation(self, k):
        return self._creation_by_determinant(k)

    def _creation_by_determinant(self, k):
        S = self.parent()
        q,t = S.q, S.t
        def f(part):
            part += [0]*(k-len(part))

            if len(part) > k:
                raise ValueError, "the column to add is too small"

            #Create the matrix over the homogeneous symmetric
            #functions and take its determinant
            MS = MatrixSpace(sfa.SFAHomogeneous(S.base_ring()), k, k)
            h  = MS.base_ring()
            m = []
            for i in range(k):
                row = [0]*max(0, (i+1)-2-part[i])
                for j in range(max(0, (i+1)-1-part[i]),k):
                    row.append( (1-q**(part[i]+j-i+1)*t**(k-(j+1)))*h([part[i]+j-i+1]) )
                m.append(row)
            M = MS(m)
            res = M.det()

            #Convert to the Schurs
            res = S._s( res )

            return S._from_element(res)



        return S._apply_module_morphism(self, f)
