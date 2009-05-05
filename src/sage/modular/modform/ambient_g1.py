r"""
Modular Forms for `\Gamma_1(N)` over `\QQ`.

EXAMPLES::

    sage: M = ModularForms(Gamma1(13),2); M
    Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: S = M.cuspidal_submodule(); S
    Cuspidal subspace of dimension 2 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: S.basis()
    [
    q - 4*q^3 - q^4 + 3*q^5 + O(q^6),
    q^2 - 2*q^3 - q^4 + 2*q^5 + O(q^6)
    ]

TESTS::

    sage: m = ModularForms(Gamma1(20),2)
    sage: loads(dumps(m)) == m
    True
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import sage.rings.all as rings

import sage.modular.arithgroup.all as arithgroup

import ambient
import cuspidal_submodule
import eisenstein_submodule
import submodule

class ModularFormsAmbient_g1_Q(ambient.ModularFormsAmbient):
    """
    A space of modular forms for the group `\Gamma_1(N)` over the rational numbers.
    """
    def __init__(self, level, weight):
        r"""
        Create a space of modular forms for `\Gamma_1(N)` of integral weight over the
        rational numbers.

        EXAMPLES::

            sage: m = ModularForms(Gamma1(100),5); m
            Modular Forms space of dimension 1270 for Congruence Subgroup Gamma1(100) of weight 5 over Rational Field
            sage: type(m)
            <class 'sage.modular.modform.ambient_g1.ModularFormsAmbient_g1_Q'>
        """
        ambient.ModularFormsAmbient.__init__(self, arithgroup.Gamma1(level), weight, rings.QQ)

    ####################################################################
    # Computation of Special Submodules
    ####################################################################
    def cuspidal_submodule(self):
        """
        Return the cuspidal submodule of this modular forms space.

        EXAMPLES::

            sage: m = ModularForms(Gamma1(17),2); m
            Modular Forms space of dimension 20 for Congruence Subgroup Gamma1(17) of weight 2 over Rational Field
            sage: m.cuspidal_submodule()
            Cuspidal subspace of dimension 5 of Modular Forms space of dimension 20 for Congruence Subgroup Gamma1(17) of weight 2 over Rational Field
        """
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            if self.level() == 1:
                self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_level1_Q(self)
            else:
                self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_g1_Q(self)
        return self.__cuspidal_submodule

    def eisenstein_submodule(self):
        """
        Return the Eisenstein submodule of this modular forms space.

        EXAMPLES::

            sage: ModularForms(Gamma1(13),2).eisenstein_submodule()
            Eisenstein subspace of dimension 11 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
            sage: ModularForms(Gamma1(13),10).eisenstein_submodule()
            Eisenstein subspace of dimension 12 of Modular Forms space of dimension 69 for Congruence Subgroup Gamma1(13) of weight 10 over Rational Field
        """
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = eisenstein_submodule.EisensteinSubmodule_g1_Q(self)
        return self.__eisenstein_submodule

    def _compute_hecke_matrix(self, n):
        r"""
        Compute the matrix of the Hecke operator T_n acting on this space.

        EXAMPLE::

            sage: ModularForms(Gamma1(7), 4).hecke_matrix(3) # indirect doctest
            [           0          -42          133            0            0            0            0            0            0]
            [           0          -28           91            0            0            0            0            0            0]
            [           1           -8           19            0            0            0            0            0            0]
            [           0            0            0           28            0            0            0            0            0]
            [           0            0            0   -10152/259            0      5222/37    -13230/37    -22295/37     92504/37]
            [           0            0            0    -6087/259            0  312067/4329 1370420/4329   252805/333 3441466/4329]
            [           0            0            0     -729/259            1       485/37      3402/37      5733/37      7973/37]
            [           0            0            0      729/259            0      -189/37     -1404/37     -2366/37     -3348/37]
            [           0            0            0      255/259            0  -18280/4329  -51947/4329   -10192/333 -190855/4329]
        """
        return self.cuspidal_submodule().hecke_matrix(n).block_sum(self.eisenstein_submodule().hecke_matrix(n))
