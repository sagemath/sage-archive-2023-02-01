r"""
Torsion subgroups of elliptic curves over rational field.

TODO:
    -- Torsion subgroups over number fields!

AUTHORS:
    -- Nick Alexander: original implementation
"""

from sage.structure.sage_object import SageObject

import sage.misc.misc as misc
import sage.rings.all as rings

import sage.schemes.elliptic_curves.ell_rational_field
import sage.groups.abelian_gps.abelian_group as groups

class EllipticCurveTorsionSubgroup(groups.AbelianGroup_class):
    r"""
    The torsion subgroup of an elliptic curve over the rational field.

    EXAMPLES:
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

    AUTHOR: Nick Alexander
    """
    def __init__(self, E, algorithm="pari"):
        self.__E = E

        flag = 0
        if algorithm=="doud":
            flag = 1
        if algorithm=="lutz_nagell":
            flag = 2
        G = None
        loop = 0
        while G is None and loop < 3:
            loop += 1
            try:
                G = self.__E.pari_curve().elltors(flag) # pari_curve will return the curve of maximum known precision
            except RuntimeError:
                self.__E.pari_curve(factor = 2) # caches a curve of double the precision
        if G is None:
            raise RuntimeError, "Could not compute torsion subgroup"

        order = G[0].python()
        structure = G[1].python()
        gens = G[2].python()
        # print order, structure, gens

        # for consistency, we sort our generators by height
        self.__torsion_gens = [ self.__E(P) for P in gens ]
        self.__torsion_gens.sort(key=lambda P: P.height())
        groups.AbelianGroup_class.__init__(self, order, structure, names='P')

    def _repr_(self):
        return "Torsion Subgroup isomorphic to %s associated to the %s" % (groups.AbelianGroup_class._repr_(self), self.__E)

    # Derived class *must* define gen method.
    def gen(self, i=0):
        return self.__torsion_gens[i]
