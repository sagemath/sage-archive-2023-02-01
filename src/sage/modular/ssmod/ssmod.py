"""
The module of supersingular points

AUTHORS:
    -- William Stein, David Kohel, and Iftikhar Burhanuddin
"""

WARN=True

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 Iftikhar Burhanuddin <burhanud@usc.edu>
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


import sage.modular.hecke.all as hecke
import sage.rings.all as rings
from sage.matrix.matrix_space import MatrixSpace
from sage.modular.dims import dimension_modular_forms
from sage.modular.congroup import Gamma0, Gamma1
from sage.databases.db_class_polynomials import HilbertClassPolynomialDatabase
from sage.databases.db_modular_polynomials \
     import ClassicalModularPolynomialDatabase


def Phi2(x,j):
    j2 = j**2
    j3 = j2*j
    return x.parent()([j3 - 162000*j2 + 8748000000*j - 157464000000000, \
                       1488*j2 + 40773375*j + 8748000000, \
                       - (j2 - 1488*j + 162000), \
                       1])

def Phi2_quad(J3, ssJ1, ssJ2):
    """
    Given a root (j1,j2) to the polynomial Phi_2(J1,J2), the pairs
    (j2,j3) not equal to (j2,j1) which solve Phi_2(j2,j3) are roots
    of the quadratic equation:

    J3^2 + (-j2^2 + 1488*j2 + (j1 - 162000))*J3 + (-j1 + 1488)*j2^2
    + (1488*j1 + 40773375)*j2 + j1^2 - 162000*j1 + 8748000000

    This will be of use to extend the 2-isogeny graph, once the
    initial three roots are determined for Phi_2(J1,J2).
    """
    ssJ1_square = ssJ1**2
    ssJ2_square = ssJ2**2

    return J3.parent()([(-ssJ1 + 1488)*ssJ2_square+ (1488*ssJ1 + 40773375)*ssJ2 \
			  + ssJ1_square - 162000*ssJ1 + 8748000000, \
			-ssJ2_square + 1488*ssJ2 + (ssJ1 - 162000), \
			1])


def dimension_supersingular_module(prime, level):
    if level == 1:
        return dimension_modular_forms(Gamma0(prime), 2)
    #elif (conductorN in [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 25]): #genus(X_0(N))zero list
        #compute basis
    else:
        raise NotImplementedError



def supersingular_D(prime):
    #Making picking D more intelligent
    D = 0
    while True:
        D = D - 1
        Dmod4 = rings.Mod(D,4)
        if Dmod4 in (0,1) and (rings.kronecker(D,prime) != 1):
            return D

def supersingular_j(FF):
    """
    Find a supersingular j-invariant.

    Example: p = 15073 has no class number one cm_j_invariant
    """
    prime = FF.characteristic()
    if not(rings.Integer(prime).is_prime()):
        raise ValueError, "%s is not a prime"%prime
    if rings.kronecker(-1, prime) != 1:
        j_invss = 1728                 #(2^2 * 3)^3
    elif rings.kronecker(-2, prime) != 1:
        j_invss = 8000                 #(2^2 * 5)^3
    elif rings.kronecker(-3, prime) != 1:
        j_invss = 0                    #0^3
    elif rings.kronecker(-7, prime) != 1:
        j_invss = 16581375             #(3 * 5 * 17)^3
    elif rings.kronecker(-11, prime) != 1:
        j_invss = -32768               #-(2^5)^3
    elif rings.kronecker(-19, prime) != 1:
        j_invss = -884736              #-(2^5 * 3)^3
    elif rings.kronecker(-43, prime) != 1:
        j_invss = -884736000           #-(2^6 * 3 * 5)^3
    elif rings.kronecker(-67, prime) != 1:
        j_invss = -147197952000        #-(2^5 * 3 * 5 * 11)^3
    elif rings.kronecker(-163, prime) != 1:
        j_invss = -262537412640768000  #-(2^6 * 3 * 5 * 23 * 29)^3
    else:
        D = supersingular_D(prime)
        DBCP = HilbertClassPolynomialDatabase()
        hc_poly = rings.PolynomialRing(FF)(DBCP[D])
        root_hc_poly_list = list(hc_poly.roots())
        j_invss = root_hc_poly_list[0][0]
    return FF(j_invss)

class SupersingularModule(hecke.HeckeModule_free_module):
    def __init__(self, prime=2, level=1, base_ring=rings.IntegerRing()):
        if WARN:
            print "Supersingular Module -- work in progress; use at own risk. (2006-08-08)"
        self.__prime = prime
        self.__finite_field = rings.FiniteField(prime**2)
        self.__level = level
        self.__base_ring = base_ring
        self.__hecke_matrices = {}
        hecke.HeckeModule_free_module.__init__(
            self, base_ring, prime*level, weight=2)

    def __repr__(self):
        return "Module of supersingular points on X_0(%s)/F_%s over %s"%(
            self.__level, self.__prime, self.__base_ring)

    def base_ring(self):
        return self.__base_ring

    def dimension(self):
        try:
            return self.__dimension
        except:
            pass
        if self.__level == 1:
            G = Gamma0(self.__prime)
            self.__dimension = dimension_modular_forms(G, 2)
        else:
            raise NotImplementedError
        return self.__dimension

    rank = dimension

    def level(self):
        return self.__level

    def prime(self):
        return self.__prime

    def weight(self):
        return 2

    def supersingular_points(self):
        try:
            return (self._ss_points_dic, self._ss_points)
        except:
            pass
        Fp2 = self.__finite_field
        level = self.__level
        prime = Fp2.characteristic()
        X = rings.PolynomialRing(Fp2).gen()
        jinv = supersingular_j(Fp2)

        dim = dimension_supersingular_module(prime, level)

        pos = int(0)
        #using list to keep track of explored nodes using pos
        ss_points = [jinv]

        #using  to keep track of index of the previous node
        ss_points_pre = [-1]

        #using dictionary for fast j-invariant look-up
        ss_points_dic = {jinv:pos}

        T2_matrix = MatrixSpace(rings.Integers(), dim, sparse=True)(0)

        while pos < len(ss_points):
            if pos == 0:
                neighs = Phi2(X,ss_points[pos]).roots()
            else:
                j_prev = ss_points_pre[pos]
                neighs = Phi2_quad(X, ss_points[j_prev], ss_points[pos]).roots()

            for (xj,ej) in neighs:
                if not ss_points_dic.has_key(xj):
                    j = len(ss_points)
                    ss_points += [xj]
                    ss_points_pre += [pos]
                    ss_points_dic[xj] = j
                else:
                    j = ss_points_dic[xj]
                T2_matrix[pos, j] += ej
                if pos != 0:
                    # also record the root from j_prev
                    T2_matrix[pos, j_prev] += 1
            pos += int(1)

        self.__hecke_matrices[2] = T2_matrix
        return (ss_points, ss_points_dic)

    def hecke_matrix(self,L):
        if self.__hecke_matrices.has_key(L):
            return self.__hecke_matrices[L]
        SS, II = self.supersingular_points()
        DBMP = ClassicalModularPolynomialDatabase()
        phi_L = DBMP[L]
        Fp2 = self.__finite_field
        h = len(SS)
        R = self.__base_ring
        T_L = MatrixSpace(R,h,sparse=True)(0)
        S, X = rings.PolynomialRing(Fp2).objgen()
        for i in range(len(SS)):
            x_i = SS[i]
            phi_L_in_y = phi_L(S(x_i),X)
            rts = phi_L_in_y.roots()
            for r in rts:
                T_L[i,int(II[r[0]])] = r[1]
        self.__hecke_matrices[L] = T_L
        return T_L
