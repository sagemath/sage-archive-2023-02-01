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
from sage.rings.all import ZZ

class ModularFormsAmbient_R(ambient.ModularFormsAmbient):
    def __init__(self, M, base_ring):
        """
        Ambient space of modular forms over a ring other than QQ.

        EXAMPLES::

            sage: M = ModularForms(23,2,base_ring=GF(7)) ## indirect doctest
            sage: M
            Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(23) of weight 2 over Finite Field of size 7
            sage: M == loads(dumps(M))
            True
        """
        self.__M = M
        if M.character() is not None:
            self.__R_character = M.character().change_ring(base_ring)
        else:
            self.__R_character = None
        ambient.ModularFormsAmbient.__init__(self, M.group(), M.weight(), base_ring, M.character())

    def modular_symbols(self,sign=0):
        r"""
        Return the space of modular symbols attached to this space, with the given sign (default 0).

        TESTS::

            sage: K.<i> = QuadraticField(-1)
            sage: chi = DirichletGroup(5, base_ring = K).0
            sage: L.<c> = K.extension(x^2 - 402*i)
            sage: M = ModularForms(chi, 7, base_ring = L)
            sage: symbs = M.modular_symbols()
            sage: symbs.character() == chi
            True
            sage: symbs.base_ring() == L
            True
        """
        sign = ZZ(sign)
        try:
            return self.__modular_symbols[sign]
        except AttributeError:
            self.__modular_symbols = {}
        except KeyError:
            pass
        M = self.__M.modular_symbols(sign).change_ring(self.base_ring())
        self.__modular_symbols[sign] = M
        return M

    def _repr_(self):
        """
        String representation for self.

        EXAMPLES::

            sage: M = ModularForms(23,2,base_ring=GF(7)) ## indirect doctest
            sage: M._repr_()
            'Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(23) of weight 2 over Finite Field of size 7'

            sage: chi = DirichletGroup(109).0 ** 36
            sage: ModularForms(chi, 2, base_ring = chi.base_ring())
            Modular Forms space of dimension 9, character [zeta3] and weight 2 over Cyclotomic Field of order 108 and degree 36
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
            <class 'sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_R_with_category'>
        """
        return CuspidalSubmodule_R(self)


    def change_ring(self, R):
        r"""
        Return this modular forms space with the base ring changed to the ring R.

        EXAMPLE::

            sage: chi = DirichletGroup(109, CyclotomicField(3)).0
            sage: M9 = ModularForms(chi, 2, base_ring = CyclotomicField(9))
            sage: M9.change_ring(CyclotomicField(15))
            Modular Forms space of dimension 10, character [zeta3 + 1] and weight 2 over Cyclotomic Field of order 15 and degree 8
            sage: M9.change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: Space cannot be defined over Rational Field
        """
        if not R.has_coerce_map_from(self.__M.base_ring()):
            raise ValueError("Space cannot be defined over %s" % R)
        return ModularFormsAmbient_R(self.__M, R)
