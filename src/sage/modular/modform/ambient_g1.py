r"""
Modular Forms for $\Gamma_1(N)$ over $\Q$.

EXAMPLES:
    sage: M = ModularForms(Gamma1(13),2); M
    Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: S = M.cuspidal_submodule(); S
    Cuspidal subspace of dimension 2 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: S.basis()
    [
    q - 4*q^3 - q^4 + 3*q^5 + O(q^6),
    q^2 - 2*q^3 - q^4 + 2*q^5 + O(q^6)
    ]
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import sage.rings.all as rings

import sage.modular.congroup as congroup

import ambient
import cuspidal_submodule
import eisenstein_submodule
import submodule

class ModularFormsAmbient_g1_Q(ambient.ModularFormsAmbient):
    def __init__(self, level, weight):
        ambient.ModularFormsAmbient.__init__(self, congroup.Gamma1(level), weight, rings.QQ)

    ####################################################################
    # Computation of Special Submodules
    ####################################################################
    def cuspidal_submodule(self):
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            if self.level() == 1:
                self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_level1_Q(self)
            else:
                self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_g1_Q(self)
        return self.__cuspidal_submodule

    def eisenstein_submodule(self):
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = eisenstein_submodule.EisensteinSubmodule_g1_Q(self)
        return self.__eisenstein_submodule

