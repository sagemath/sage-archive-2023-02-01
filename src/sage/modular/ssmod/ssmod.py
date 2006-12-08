"""
The module of supersingular points

AUTHORS:
    -- William Stein
    -- David Kohel
    -- Iftikhar Burhanuddin
"""

#*****************************************************************************
#       Copyright (C) 2004,2006 William Stein <wstein@gmail.com>
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
    r"""
    This function returns a certain cubic polynomial in the
    indeterminate x over a finite field.

    The roots of the \emph{modular} polynomial $\Phi_2(x,j)$ are the
    2-isogenous supersingular j-invariants of j.

    INPUT:
       x -- indeterminate of a univariate polynomial ring defined over
            a finite field with p^2 elements, where p is a prime
            number
       j -- supersingular j-invariant over the finite field

    OUTPUT:
       polynomial -- defined over the finite field

    EXAMPLES:
    The following code snippet produces the modular polynomial
    $\Phi_{2}(x,j_{in})$, where $j_{in}$ is a supersingular j-invariant
    defined over the finite field with $7^2$ elements.
        sage: F = GF(7^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: sage.modular.ssmod.ssmod.Phi2(X,j_in)
        x^3 + 3*x^2 + 3*x + 1

    AUTHORS:
       William Stein - wstein@gmail.com
       David Kohel - kohel@maths.usyd.edu.au
       Iftikhar Burhanuddin - burhanud@usc.edu
    """
    j_pow2 = j**2
    j_pow3 = j_pow2*j
    return x.parent()([j_pow3 - 162000*j_pow2 + 8748000000*j -
    157464000000000, 1488*j_pow2 + 40773375*j + 8748000000, - (j_pow2
    - 1488*j + 162000), 1])

def Phi2_quad(J3, ssJ1, ssJ2):
    r"""
    This function returns a certain quadratic polynomial over a finite
    field in indeterminate J3.

    The roots of the polynomial along with ssJ1 are the
    neighboring/2-isogenous supersingular j-invariants of ssJ2.

    INPUT:
       J3 -- indeterminate of a univariate polynomial ring defined
       over a finite field with p^2 elements where p is a prime number
       ssJ2, ssJ2 -- supersingular j-invariants over the finite field

    OUTPUT:
       polynomial -- defined over the finite field

    EXAMPLES:
    The following code snippet produces a factor of the modular polynomial
    $\Phi_{2}(x,j_{in})$, where $j_{in}$ is a supersingular j-invariant
    defined over the finite field with $37^2$ elements.
        sage: F = GF(37^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: poly = sage.modular.ssmod.ssmod.Phi2(X,j_in)
        sage: poly.roots()
        [(8, 1), (27*a + 23, 1), (10*a + 20, 1)]
        sage: sage.modular.ssmod.ssmod.Phi2_quad(X, F(8), j_in)
        x^2 + 31*x + 31

    NOTES:
    Given a root (j1,j2) to the polynomial $Phi_2(J1,J2)$, the pairs
    (j2,j3) not equal to (j2,j1) which solve $Phi_2(j2,j3)$ are roots
    of the quadratic equation:

    \begin{verbatim}
    J3^2 + (-j2^2 + 1488*j2 + (j1 - 162000))*J3 + (-j1 + 1488)*j2^2
    + (1488*j1 + 40773375)*j2 + j1^2 - 162000*j1 + 8748000000
    \end{verbatim}

    This will be of use to extend the 2-isogeny graph, once the
    initial three roots are determined for $Phi_2(J1,J2)$.

    AUTHORS:
       David Kohel - kohel@maths.usyd.edu.au
       Iftikhar Burhanuddin - burhanud@usc.edu
    """
    ssJ1_pow2 = ssJ1**2
    ssJ2_pow2 = ssJ2**2

    return J3.parent()([(-ssJ1 + 1488)*ssJ2_pow2+ (1488*ssJ1 +
    40773375)*ssJ2 + ssJ1_pow2 - 162000*ssJ1 + 8748000000,
    -ssJ2_pow2 + 1488*ssJ2 + (ssJ1 - 162000), 1])

def Phi3(x,j):
    r"""
    This function returns a certain quartic polynomial in the indeterminate x
    over a finite field.

    The roots of the \emph{modular} polynomial $\Phi_3(x,j)$ are the
    3-isogenous supersingular j-invariants of j.

    INPUT:
       x -- indeterminate of a univariate polynomial ring defined over
            a finite field with p^2 elements, where p is a prime
            number
       j -- supersingular j-invariant over the finite field

    OUTPUT:
       polynomial -- defined over the finite field

    EXAMPLES:
    The following code snippet produces the modular polynomial
    $\Phi_{3}(x,j_{in})$, where $j_{in}$ is a supersingular j-invariant
    defined over the finite field with $7^2$ elements.
        sage: F = GF(7^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: sage.modular.ssmod.ssmod.Phi3(X,j_in)
        x^4 + 4*x^3 + 6*x^2 + 4*x + 1

    AUTHORS:
       Iftikhar Burhanuddin - burhanud@usc.eduw
    """
    j_pow2 = j**2
    j_pow3 = j_pow2*j
    j_pow4 = j_pow3*j

    return x.parent()([1855425871872000000000*j +
    452984832000000*j_pow2 + 36864000*j_pow3 + j_pow4,
    1855425871872000000000 - 770845966336000000*j +
    8900222976000*j_pow2 - 1069956*j_pow3, 452984832000000 +
    8900222976000*j + 2587918086*j_pow2 + 2232*j_pow3, 36864000 -
    1069956*j + 2232*j_pow2 - j_pow3, 1])

def Phi5(x,j):
    r"""
    This function returns a particular sextic polynomial in the indeterminate x
    over a finite field.

    The roots of the \emph{modular} polynomial $\Phi_5(x,j)$ are the
    5-isogenous supersingular j-invariants of j.

    INPUT:
       x -- indeterminate of a univariate polynomial ring defined over
            a finite field with p^2 elements, where p is a prime
            number
       j -- supersingular j-invariant over the finite field

    OUTPUT:
       polynomial -- defined over the finite field

    EXAMPLES:
    The following code snippet produces the modular polynomial
    $\Phi_{5}(x,j_{in})$, where $j_{in}$ is a supersingular j-invariant
    defined over the finite field with $7^2$ elements.
        sage: F = GF(7^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: sage.modular.ssmod.ssmod.Phi5(X,j_in)
        x^6 + 6*x^5 + x^4 + 6*x^3 + x^2 + 6*x + 1

    AUTHORS:
       Iftikhar Burhanuddin - burhanud@usc.eduw

    """
    j_pow2 = j**2
    j_pow3 = j_pow2*j
    j_pow4 = j_pow3*j
    j_pow5 = j_pow4*j
    j_pow6 = j_pow5*j

    return x.parent()([
        141359947154721358697753474691071362751004672000 +
        53274330803424425450420160273356509151232000*j +
        6692500042627997708487149415015068467200*j_pow2 +
        280244777828439527804321565297868800*j_pow3 +
        1284733132841424456253440*j_pow4 + 1963211489280*j_pow5 +
        j_pow6, 53274330803424425450420160273356509151232000 -
        264073457076620596259715790247978782949376*j +
        36554736583949629295706472332656640000*j_pow2 -
        192457934618928299655108231168000*j_pow3 +
        128541798906828816384000*j_pow4 - 246683410950*j_pow5,
        6692500042627997708487149415015068467200 +
        36554736583949629295706472332656640000*j +
        5110941777552418083110765199360000*j_pow2 +
        26898488858380731577417728000*j_pow3 +
        383083609779811215375*j_pow4 + 2028551200*j_pow5,
        280244777828439527804321565297868800 -
        192457934618928299655108231168000*j +
        26898488858380731577417728000*j_pow2 -
        441206965512914835246100*j_pow3 + 107878928185336800*j_pow4 -
        4550940*j_pow5, 1284733132841424456253440 +
        128541798906828816384000*j + 383083609779811215375*j_pow2 +
        107878928185336800*j_pow3 + 1665999364600*j_pow4 +
        3720*j_pow5, 1963211489280 - 246683410950*j +
        2028551200*j_pow2 - 4550940*j_pow3 + 3720*j_pow4 - j_pow5, 1])


def dimension_supersingular_module(prime, level=1):
    r"""
    This function returns the dimension Supersingular module, which is
    equal to the dimension of the space of cusp forms of weight $2$
    and conductor equal to prime times level.

    INPUT:
       prime -- integer, prime
       level -- integer, positive

    OUTPUT:
       dimension -- integer, nonnegative

    EXAMPLES:
    The code below illustrates the usage of this function.
        sage: dimension_supersingular_module(7)
        1

        sage: dimension_supersingular_module(15073)
        1256

        sage: dimension_supersingular_module(83401)
        6950

    NOTES:
    The case of level > 1 has not been implemented yet.

    AUTHORS:
       David Kohel - kohel@maths.usyd.edu.au
       Iftikhar Burhanuddin - burhanud@usc.edu
    """
    if not(rings.Integer(prime).is_prime()):
        raise ValueError, "%s is not a prime"%prime

    if level == 1:
        return dimension_modular_forms(Gamma0(prime), 2)

    #list of genus(X_0(level)) equal to zero
    #elif (level in [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 25]):
    #compute basis

    else:
        raise NotImplementedError


def supersingular_D(prime):
    r"""
    This function returns a fundamental discriminant $D$ of an
    imaginary quadratic field, where the given prime does not split
    (see Silverman's Advanced Topics in the Arithmetic of Elliptic
    Curves, page 184, exercie 2.30(d).)

    INPUT:
        prime -- integer, prime

    OUTPUT:
        D -- integer, negative

    EXAMPLES:
    These examples illustrate the functionality of the procedure.
        sage: supersingular_D(7)
        -4

        sage: supersingular_D(15073)
        -15

        sage: supersingular_D(83401)
        -7

    AUTHORS:
       David Kohel - kohel@maths.usyd.edu.au
       Iftikhar Burhanuddin - burhanud@usc.edu
    """
    if not(rings.Integer(prime).is_prime()):
        raise ValueError, "%s is not a prime"%prime

    #Making picking D more intelligent
    D = -1
    while True:
        Dmod4 = rings.Mod(D,4)
        if Dmod4 in (0,1) and (rings.kronecker(D,prime) != 1):
            return D
        D = D - 1

def supersingular_j(FF):
    r"""
    This function returns a supersingular j-invariant over the finite
    field FF.

    INPUT:
       FF  -- finite field with p^2 elements, where p is a prime number

    OUTPUT:
       finite field element -- a supersingular j-invariant
       defined over the finite field FF

    EXAMPLES:
    The following examples showcase the usage of the function.
        sage: supersingular_j(GF(7^2, 'a'))
        6

        sage: supersingular_j(GF(15073^2,'a'))
        4443*a + 13964

        sage: supersingular_j(GF(83401^2, 'a'))
        67977

    AUTHORS:
       David Kohel - kohel@maths.usyd.edu.au
       Iftikhar Burhanuddin - burhanud@usc.edu
    """
    if not(FF.is_field()) or not(FF.is_finite()):
        raise ValueError, "%s is not a finite field"%FF
    prime = FF.characteristic()
    if not(rings.Integer(prime).is_prime()):
        raise ValueError, "%s is not a prime"%prime
    if not(rings.Integer(FF.cardinality())) == rings.Integer(prime**2):
        raise ValueError, "%s is not a quadratic extension"%FF
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
        hc_poly = rings.PolynomialRing(FF, 'x')(DBCP[D])
        root_hc_poly_list = list(hc_poly.roots())

        j_invss = root_hc_poly_list[0][0]
    return FF(j_invss)

class SupersingularModule(hecke.HeckeModule_free_module):
    def __init__(self, prime=2, level=1, base_ring=rings.IntegerRing()):
        self.__prime = prime
        self.__finite_field = rings.FiniteField(prime**2,'a')
        self.__level = level
        self.__hecke_matrices = {}
        hecke.HeckeModule_free_module.__init__(
            self, base_ring, prime*level, weight=2)

    def __repr__(self):
        return "Module of supersingular points on X_0(%s)/F_%s over %s"%(
            self.__level, self.__prime, self.base_ring())

    def dimension(self):
        r"""
        This function returns the dimension of the space of modular
        forms of weight 2 and level equal to the level associated to
        self.

        INPUT:
            self -- SupersingularModule object

        OUTPUT:
            integer -- dimension, nonnegative

        EXAMPLES:
            sage: S = SupersingularModule(7)
            sage: S.dimension()
            1

            sage: S = SupersingularModule(15073)
            sage: S.dimension()
            1256

            sage: S = SupersingularModule(83401)
            sage: S.dimension()
            6950

        NOTES:
           The case of level > 1 has not yet been implemented.

        AUTHORS:
            David Kohel - kohel@maths.usyd.edu.au
            Iftikhar Burhanuddin - burhanud@usc.edu
        """
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
        r"""
        This function returns the level associated to self.

        INPUT:
            self -- SupersingularModule object

        OUTPUT:
            integer -- the level, positive

        EXAMPLES:
            sage: S = SupersingularModule(15073)
            sage: S.level()
            1

        AUTHORS:
            David Kohel - kohel@maths.usyd.edu.au
            Iftikhar Burhanuddin - burhanud@usc.edu
        """
        return self.__level

    def prime(self):
        r"""
        This function returns the characteristic of the finite field
        associated to self.

        INPUT:
            self -- SupersingularModule object

        OUTPUT:
            integer -- charateristic, positive

        EXAMPLES:
        This example shows the usage of this function.
            sage: S = SupersingularModule(19)
            sage: S.prime()
            19

        AUTHORS:
            David Kohel - kohel@maths.usyd.edu.au
            Iftikhar Burhanuddin - burhanud@usc.edu
        """
        return self.__prime

    def weight(self):
        r"""
        This function returns the weight associated to self.

        INPUT:
            self -- SupersingularModule object

        OUTPUT:
            integer -- weight, positive

        EXAMPLES:
        This example shows the usage of this function.
            sage: S = SupersingularModule(19)
            sage: S.weight()
            2

        AUTHORS:
            David Kohel - kohel@maths.usyd.edu.au
            Iftikhar Burhanuddin - burhanud@usc.edu
        """
        return 2

    def supersingular_points(self):
        r"""
        This function computes the supersingular j-invariants over the
        finite field associated to self.

        INPUT:
            self -- SupersingularModule object

        OUTPUT: list_j, dict_j -- list_j is the list of supersingular
            j-invariants, dict_j is a dictionary with these
            j-invariants as keys and their indexes as values. The
            latter is used to speed up j-invariant look-up. The
            indexes are based on the order of their \emph{discovery}.

        EXAMPLES:
        The following examples describe the usage of the function.
            sage: S = SupersingularModule(7)
            sage: S.supersingular_points()
            ([6], {6: 0})

            sage: S = SupersingularModule(11)
            sage: S.supersingular_points()
            ([1, 0], {1: 0, 0: 1})

            sage: S = SupersingularModule(37)
            sage: S.supersingular_points()
            ([8, 27*a + 23, 10*a + 20], {8: 0, 10*a + 20: 2, 27*a + 23: 1})

        AUTHORS:
            David Kohel - kohel@maths.usyd.edu.au
            Iftikhar Burhanuddin - burhanud@usc.edu
        """
        try:
            return (self._ss_points_dic, self._ss_points)
        except:
            pass
        Fp2 = self.__finite_field
        level = self.__level
        prime = Fp2.characteristic()
        X = rings.PolynomialRing(Fp2, 'x').gen()
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
        r"""
        This function returns the $L$th Hecke matrix.

        INPUT:
            self -- SupersingularModule object
            L -- integer, positive

        OUTPUT:
            matrix -- sparse integer matrix

        EXAMPLES:
        This example computes the action of the Hecke operator $T_2$
        on the module of supersingular points on $X_0(1)/F_37$.
            sage: S = SupersingularModule(37)
            sage: M = S.hecke_matrix(2)
            sage: M
            [1 1 1]
            [1 0 2]
            [1 2 0]

        AUTHORS:
            David Kohel - kohel@maths.usyd.edu.au
            Iftikhar Burhanuddin - burhanud@usc.edu
        """
        if self.__hecke_matrices.has_key(L):
            return self.__hecke_matrices[L]
        SS, II = self.supersingular_points()
        DBMP = ClassicalModularPolynomialDatabase()
        phi_L = DBMP[L]
        Fp2 = self.__finite_field
        h = len(SS)
        R = self.base_ring()
        T_L = MatrixSpace(R,h,sparse=True)(0)
        S, X = rings.PolynomialRing(Fp2, 'x').objgen()
        for i in range(len(SS)):
            x_i = SS[i]
            phi_L_in_y = phi_L(S(x_i),X)
            rts = phi_L_in_y.roots()
            for r in rts:
                T_L[i,int(II[r[0]])] = r[1]
        self.__hecke_matrices[L] = T_L
        return T_L
