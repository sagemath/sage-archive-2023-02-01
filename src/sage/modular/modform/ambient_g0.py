r"""
Modular Forms for $\Gamma_0(N)$ over $\Q$.
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

class ModularFormsAmbient_g0_Q(ambient.ModularFormsAmbient):
    def __init__(self, level, weight):
        ambient.ModularFormsAmbient.__init__(self, congroup.Gamma0(level), weight, rings.QQ)

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
                self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_g0_Q(self)
        return self.__cuspidal_submodule

    def eisenstein_submodule(self):
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = eisenstein_submodule.EisensteinSubmodule_g0_Q(self)
        return self.__eisenstein_submodule

