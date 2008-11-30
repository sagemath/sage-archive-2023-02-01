r"""
Torsion subgroups of elliptic curves over number fields (including Q).

AUTHORS:
    -- Nick Alexander: original implementation over Q
    -- Chris Wuthrich: original implementation over number fields
    -- John Cremona: rewrote p-primary part to use division
       polynomials, added some features, unified NF and Q code.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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


from sage.structure.sage_object import SageObject

import sage.misc.misc as misc
import sage.rings.all as rings
from sage.rings.all import (Integer, RationalField)
import sage.groups.abelian_gps.abelian_group as groups
import sage.groups.generic as generic

class EllipticCurveTorsionSubgroup(groups.AbelianGroup_class):
    r"""
    The torsion subgroup of an elliptic curve over a number field.

    EXAMPLES:

    Examples over Q:

        sage: E = EllipticCurve([-4, 0]); E
        Elliptic Curve defined by y^2  = x^3 - 4*x over Rational Field
        sage: G = E.torsion_subgroup(); G
        Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C2 x C2 associated to the Elliptic Curve defined by y^2  = x^3 - 4*x over Rational Field
        sage: G.order()
        4
        sage: G.gen(0)
        (2 : 0 : 1)
        sage: G.gen(1)
        (0 : 0 : 1)
        sage: G.ngens()
        2

        sage: E = EllipticCurve([17, -120, -60, 0, 0]); E
        Elliptic Curve defined by y^2 + 17*x*y - 60*y = x^3 - 120*x^2 over Rational Field
        sage: G = E.torsion_subgroup(); G
        Torsion Subgroup isomorphic to Trivial Abelian Group associated to the Elliptic Curve defined by y^2 + 17*x*y - 60*y = x^3 - 120*x^2 over Rational Field
        sage: G.gens()
        ()
        sage: e = EllipticCurve([0, 33076156654533652066609946884,0,\
        347897536144342179642120321790729023127716119338758604800,\
        1141128154369274295519023032806804247788154621049857648870032370285851781352816640000])
        sage: e.torsion_order()
        16

        Constructing points from the torsion subgroup (which is an
        abstract abelian group):
            sage: E = EllipticCurve('14a1')
            sage: T = E.torsion_subgroup()
            sage: [E(t) for t in T]
            [(0 : 1 : 0),
            (9 : 23 : 1),
            (2 : 2 : 1),
            (1 : -1 : 1),
            (2 : -5 : 1),
            (9 : -33 : 1)]

        An example where the torsion subgroup is not cyclic:
            sage: E = EllipticCurve([0,0,0,-49,0])
            sage: T = E.torsion_subgroup()
            sage: [E(t) for t in T]
            [(0 : 1 : 0), (0 : 0 : 1), (7 : 0 : 1), (-7 : 0 : 1)]

        Note that the following behaviour for curves with trivial
        torsion is due to a bug in the AbelianGroup class:

            sage: E = EllipticCurve('37a1')
            sage: T = E.torsion_subgroup()
            sage: T
            Torsion Subgroup isomorphic to Trivial Abelian Group associated to the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: [E(t) for t in T]
            []
            sage: T.list()
            []

    Examples over other Number Fields:

        sage: E=EllipticCurve('11a1')
        sage: K.<i>=NumberField(x^2+1)
        sage: EK=E.change_ring(K)
        sage: from sage.schemes.elliptic_curves.ell_torsion import EllipticCurveTorsionSubgroup
        sage: EllipticCurveTorsionSubgroup(EK)
        Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1

    Note: this class is normally constructed indirectly as follows:
        sage: T = EK.torsion_subgroup(); T
        Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1
        sage: type(T)
        <class 'sage.schemes.elliptic_curves.ell_torsion.EllipticCurveTorsionSubgroup'>


    AUTHOR: Nick Alexamder - initial implementation over Q
            Chris Wuthrich - initial implementation over number fields
            John Cremona - additional features and unification
    """
    def __init__(self, E, algorithm=None):
        """
        Initialization function for EllipticCurveTorsionSubgroup class

        INPUT:
            E -- An elliptic curve defined over a number field (including Q)

            algorithm -- (string, default None): If not None, must be one
                     of 'pari', 'doud', 'lutz_nagell'.  For curves
                     defined over Q, pari is then used with the
                     appropriate flag passed to pari's elltors()
                     function; this parameter is ignored for curves
                     whose base field is not Q.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_torsion import EllipticCurveTorsionSubgroup
            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: EllipticCurveTorsionSubgroup(EK)
            Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1

        Note: this class is normally constructed indirectly as follows:
            sage: T = EK.torsion_subgroup(); T
            Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1
            sage: type(T)
            <class 'sage.schemes.elliptic_curves.ell_torsion.EllipticCurveTorsionSubgroup'>
        """
        self.__E = E
        self.__K = E.base_field()

        pari_torsion_algorithms = ["pari","doud","lutz_nagell"]

        if self.__K is RationalField() and algorithm in pari_torsion_algorithms:
            flag = pari_torsion_algorithms.index(algorithm)

            G = None
            loop = 0
            while G is None and loop < 3:
                loop += 1
                try:
                    G = self.__E.pari_curve(prec = 400).elltors(flag) # pari_curve will return the curve of maximum known precision
                except RuntimeError:
                    self.__E.pari_curve(factor = 2) # caches a curve of twice the precision
            if G is not None:
                order = G[0].python()
                structure = G[1].python()
                gens = G[2].python()

                self.__torsion_gens = [ self.__E(P) for P in gens ]
                groups.AbelianGroup_class.__init__(self, order, structure, names='P')
                return

        T1 = E(0) # these will be the two generators
        T2 = E(0)
        k1 = 1    # with their order
        k2 = 1

        # find a multiple of the order of the torsion group
        bound = E._torsion_bound(number_of_places=20)

        # now do prime by prime
        for p,e in bound.factor():
            ptor = E._p_primary_torsion_basis(p)
            # print p,'-primary part is ',ptor
            if len(ptor)>0:
                T1 += ptor[0][0]
                k1 *= p**(ptor[0][1])
            if len(ptor)>1:
                T2 += ptor[1][0]
                k2 *= p**(ptor[1][1])

        order = k1*k2
        if k1 == 1:
            structure = []
            gens = []
        elif k2 == 1:
            structure = [k1]
            gens = [T1]
        else:
            structure = [k1,k2]
            gens = [T1,T2]

        self.__torsion_gens = gens
        groups.AbelianGroup_class.__init__(self, order, structure, names='P')


    def _repr_(self):
        """
        Output an instance of the EllipticCurveTorsionSubgroup class.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: T = EK.torsion_subgroup(); T._repr_()
            'Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1'
        """
        return "Torsion Subgroup isomorphic to %s associated to the %s" % (groups.AbelianGroup_class._repr_(self), self.__E)

    def gen(self, i=0):
        """
        Return the i'th torsion generator.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: T = EK.torsion_subgroup()
            sage: T.gen()
            (16 : 60 : 1)

        """
        return self.__torsion_gens[i]

    def ngens(self):
        """
        Return the number of torsion generators.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: T = EK.torsion_subgroup()
            sage: T.ngens()
            1
        """
        return len(self.__torsion_gens)

    def curve(self):
        """
        Return the curve of this torsion subgroup.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: T = EK.torsion_subgroup()
            sage: T.curve() is EK
            True
        """
        return self.__E

    def points(self):
        """
        Return a list of all the points in this torsion subgroup.
        The list is cached.

        EXAMPLES:
            sage: K.<i>=NumberField(x^2 + 1)
            sage: E = EllipticCurve(K,[0,0,0,1,0])
            sage: tor = E.torsion_subgroup()
            sage: tor.points()
            [(i : 0 : 1), (0 : 0 : 1), (-i : 0 : 1), (0 : 1 : 0)]
        """
        try:
            return self.__points
        except AttributeError:
            pass
        E = self.curve()
        ni = self.invariants()
        r = len(ni)
        if r == 0:
            return [E(0)]
        H0 = list(generic.multiples(self.gen(0),ni[0],operation='+'))
        if r == 1:
            H0.sort()
            self.__points = H0
            return self.__points
        H1 = list(generic.multiples(self.gen(1),ni[1],operation='+'))
        H = [P+Q for P in H0 for Q in H1]
        H.sort()
        self.__points = H
        return self.__points

#
#
#
#     def __p_primary_part(self,p,b):
#         r"""
#         Find the $p$-primary part of the torsion subgroup

#         INPUT :
#             p -- a prime
#             b -- an integer such that it is known that the resulting order divides p^b

#         OUTPUT :
#             [[T1,k1],[T2,k2]]
#             with T1,T2 a basis and k1>=k2 their orders

#         """
#         E = self.__E
#         if b <= 0:
#             return [[E(0),1],[E(0),1]]

#         k1 = 0
#         k2 = 0
#         k = 0  # we are computing the p^k torsion
#         fk = 1
#         while k1+k2 < b :
#             k += 1
#             fk_minus_1 = fk
#             fk = E.division_polynomial(p**k)
#             g = fk.quo_rem(fk_minus_1)[0]
#             xs = g.roots()
#             xs = [a[0] for a in xs if E.is_x_coord(a[0])]
#             # print 'torsion : ',p**k,' has x in ',xs
#             # xs now contains the list of x-coordinates of p^k - torsion points
#             if len(xs) > p:  # there are two copies of Z/p^k
#                 k1 += 1
#                 k2 += 1
#             elif len(xs) > 0: # there is one copy of Z/p^k
#                 k1 += 1
#             else:
#                 b = 0  # breaks
#         # we have determined the structure, now to the gens
#         # print 'order determined to :', [k1,k2]
#         if k1 == 0 :
#             return [[E(0),1],[E(0),1]]
#         else :
#             T1 = E.lift_x(xs[0])  # this is a p^k1 torsion point
#             if k2 == 0:
#                 T2 = E(0)
#             else:
#                 fk_minus_1 = E.division_polynomial(p**(k2-1))
#                 fk = E.division_polynomial(p**k2)
#                 g = fk.quo_rem(fk_minus_1)[0]
#                 xs = g.roots()
#                 xs = [a[0] for a in xs if E.is_x_coord(a[0])]
#                 if p == 2:
#                     li = 2
#                 else:
#                     li = (p+1)//2
#                 dp = [(c*(p**(k1-1))*T1)[0] for c in range(1,li)]
#                 # these are the x-coordinates of p-torsion points that are lin dep of T1
#                 for xx in xs:
#                     T2 = E.lift_x(xx)
#                     if not ((p**(k2-1))*T2)[0] in dp:
#                         break  # as then T2 is lin indep of the correct order
#             return [[T1,k1],[T2,k2]]


