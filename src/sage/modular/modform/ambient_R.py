"""
Modular Forms over a Non-minimal Base Ring
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import ambient

class ModularFormsAmbient_R(ambient.ModularFormsAmbient):
    def __init__(self, M, base_ring):
        self.__M = M
        ambient.ModularFormsAmbient.__init__(self, M.group(), M.weight(), base_ring)

    def _repr_(self):
        s = str(self.__M)
        i = s.find('over')
        if i != -1:
            s = s[:i]
        return s + 'over %s'%self.base_ring()

    def _compute_q_expansion_basis(self, prec):
        if self.base_ring().characteristic() == 0:
            B = self.__M.q_expansion_basis(prec)
        else:
            B = self.__M.q_expansion_integral_basis(prec)
        R = self._q_expansion_ring()
        return [R(f) for f in B]

