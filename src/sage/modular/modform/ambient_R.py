"""
Modular Forms over a Non-minimal Base Ring
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import ambient
from cuspidal_submodule import CuspidalSubmodule_R

class ModularFormsAmbient_R(ambient.ModularFormsAmbient):
    def __init__(self, M, base_ring):
        """
        Ambient space of modular forms over a ring other than QQ.

        EXAMPLES:
            sage: M = ModularForms(23,2,base_ring=GF(7)) ## indirect doctest
            sage: M
            Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(23) of weight 2 over Finite Field of size 7
            sage: M == loads(dumps(M))
            True
        """
        self.__M = M
        ambient.ModularFormsAmbient.__init__(self, M.group(), M.weight(), base_ring)

    def _repr_(self):
        """
        String representation for self.

        EXAMPLES:
            sage: M = ModularForms(23,2,base_ring=GF(7)) ## indirect doctest
            sage: M._repr_()
            'Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(23) of weight 2 over Finite Field of size 7'
        """
        s = str(self.__M)
        i = s.find('over')
        if i != -1:
            s = s[:i]
        return s + 'over %s'%self.base_ring()

    def _compute_q_expansion_basis(self, prec=None):
        """
        Compute q-expansions for a basis of self to precision prec.

        EXAMPLES:
            sage: M = ModularForms(23,2,base_ring=GF(7))
            sage: M._compute_q_expansion_basis(5)
            [1 + 5*q^3 + 5*q^4 + O(q^5),
            q + 6*q^3 + 6*q^4 + O(q^5),
            q^2 + 5*q^3 + 6*q^4 + O(q^5)]
        """
        if prec == None:
            prec = self.prec()
        if self.base_ring().characteristic() == 0:
            B = self.__M.q_expansion_basis(prec)
        else:
            B = self.__M.q_integral_basis(prec)
        R = self._q_expansion_ring()
        return [R(f) for f in B]

    def cuspidal_submodule(self):
        r"""
        Return the cuspidal subspace of this space.

        EXAMPLE::

            sage: C = CuspForms(7, 4, base_ring=CyclotomicField(5)) # indirect doctest
            sage: type(C)
            <class 'sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_R'>
        """
        return CuspidalSubmodule_R(self)
