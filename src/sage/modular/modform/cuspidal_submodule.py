"""
The Cuspidal Subspace

EXAMPLES:
    sage: S = CuspForms(SL2Z,12); S
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
    Modular Group SL(2,Z) of weight 12 over Rational Field
    sage: S.basis()
    [
    q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
    ]

    sage: S = CuspForms(Gamma0(33),2); S
    Cuspidal subspace of dimension 3 of Modular Forms space of dimension 6 for
    Congruence Subgroup Gamma0(33) of weight 2 over Rational Field
    sage: S.basis()
    [
    q - q^5 + O(q^6),
    q^2 - q^4 - q^5 + O(q^6),
    q^3 + O(q^6)
    ]

    sage: S = CuspForms(Gamma1(3),6); S
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3 for
    Congruence Subgroup Gamma1(3) of weight 6 over Rational Field
    sage: S.basis()
    [
    q - 6*q^2 + 9*q^3 + 4*q^4 + 6*q^5 + O(q^6)
    ]
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.rings.all import Integer
from sage.misc.all import verbose

import submodule
import vm_basis

class CuspidalSubmodule(submodule.ModularFormsSubmodule):
    """
    The cuspidal submodule of an ambient space of modular forms.
    """
    def __init__(self, ambient_space):
        """
        The cuspidal submodule of an ambient space of modular forms.

        EXAMPLES:
            sage: S = CuspForms(SL2Z,12); S
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
            Modular Group SL(2,Z) of weight 12 over Rational Field
            sage: S.basis()
            [
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
            ]

            sage: S = CuspForms(Gamma0(33),2); S
            Cuspidal subspace of dimension 3 of Modular Forms space of dimension 6 for
            Congruence Subgroup Gamma0(33) of weight 2 over Rational Field
            sage: S.basis()
            [
            q - q^5 + O(q^6),
            q^2 - q^4 - q^5 + O(q^6),
            q^3 + O(q^6)
            ]

            sage: S = CuspForms(Gamma1(3),6); S
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3 for
            Congruence Subgroup Gamma1(3) of weight 6 over Rational Field
            sage: S.basis()
            [
            q - 6*q^2 + 9*q^3 + 4*q^4 + 6*q^5 + O(q^6)
            ]
            sage: S == loads(dumps(S))
            True
        """
        verbose('creating cuspidal submodule of %s'%ambient_space)
        d = ambient_space._dim_cuspidal()
        V = ambient_space.module()
        G = [V.gen(i) for i in range(d)]
        S = V.submodule(G, check=False, already_echelonized=True)
        submodule.ModularFormsSubmodule.__init__(self, ambient_space, S)

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: S = CuspForms(Gamma1(3),6); S._repr_()
            'Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3 for Congruence Subgroup Gamma1(3) of weight 6 over Rational Field'
        """
        return "Cuspidal subspace of dimension %s of %s"%(self.dimension(), self.ambient_module())

    def modular_symbols(self, sign=0):
        """
        Return the corresponding space of modular symbols with the given sign.

        EXAMPLES:
            sage: S = ModularForms(11,2).cuspidal_submodule()
            sage: S.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space
            of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field

            sage: S.modular_symbols(sign=-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space
            of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

            sage: M = S.modular_symbols(sign=1); M
            Modular Symbols subspace of dimension 1 of Modular Symbols space of
            dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: M.sign()
            1

            sage: S = ModularForms(1,12).cuspidal_submodule()
            sage: S.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of
            dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field

            sage: eps = DirichletGroup(13).0
            sage: S = CuspForms(eps^2, 2)

            sage: S.modular_symbols(sign=0)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2

            sage: S.modular_symbols(sign=1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 and level 13, weight 2, character [zeta6], sign 1, over Cyclotomic Field of order 6 and degree 2

            sage: S.modular_symbols(sign=-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 1 and level 13, weight 2, character [zeta6], sign -1, over Cyclotomic Field of order 6 and degree 2
        """
        try:
            return self.__modular_symbols[sign]
        except AttributeError:
            self.__modular_symbols = {}
        except KeyError:
            pass
        A = self.ambient_module()
        S = A.modular_symbols(sign).cuspidal_submodule()
        self.__modular_symbols[sign] = S
        return S

# A cuspidal subspace of a space of modular forms
# that is computed using modular symbols.
class CuspidalSubmodule_modsym_qexp(CuspidalSubmodule):
    """
    Cuspidal submodule with q-expansions calculated via modular symbols.
    """
    def _compute_q_expansion_basis(self, prec=None):
        """
        Compute q-expansions of a basis for self (via modular symbols).

        EXAMPLES:
            sage: sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_modsym_qexp(ModularForms(11,2))._compute_q_expansion_basis()
            [
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
            ]
        """
        if prec == None:
            prec = self.prec()
        else:
            prec = Integer(prec)
        M = self.modular_symbols(sign = 1)
        return M.q_expansion_basis(prec)

class CuspidalSubmodule_level1_Q(CuspidalSubmodule):
    """
    Space of cusp forms of level 1 over Q.
    """
    def _compute_q_expansion_basis(self, prec=None):
        """
        Compute q-expansions of a basis for self.

        EXAMPLES:
            sage: sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_level1_Q(ModularForms(1,12))._compute_q_expansion_basis()
            [
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
            ]
        """
        if prec == None:
            prec = self.prec()
        else:
            prec = Integer(prec)
        return vm_basis.victor_miller_basis(self.weight(), prec,
                                            cusp_only=True, var='q')

class CuspidalSubmodule_g0_Q(CuspidalSubmodule_modsym_qexp):
    """
    Space of cusp forms for Gamma0(N).
    """

class CuspidalSubmodule_g1_Q(CuspidalSubmodule_modsym_qexp):
    """
    Space of cusp forms for Gamma1(N).
    """

class CuspidalSubmodule_eps(CuspidalSubmodule_modsym_qexp):
    """
    Space of cusp forms with given Dirichlet character.

    EXAMPLES:
        sage: S = CuspForms(DirichletGroup(5).0,5); S
        Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3, character [zeta4] and weight 5 over Cyclotomic Field of order 4 and degree 2

        sage: S.basis()
        [
        q + (-zeta4 - 1)*q^2 + (6*zeta4 - 6)*q^3 - 14*zeta4*q^4 + (15*zeta4 + 20)*q^5 + O(q^6)
        ]
        sage: f = S.0
        sage: f.qexp()
        q + (-zeta4 - 1)*q^2 + (6*zeta4 - 6)*q^3 - 14*zeta4*q^4 + (15*zeta4 + 20)*q^5 + O(q^6)
        sage: f.qexp(7)
        q + (-zeta4 - 1)*q^2 + (6*zeta4 - 6)*q^3 - 14*zeta4*q^4 + (15*zeta4 + 20)*q^5 + 12*q^6 + O(q^7)
        sage: f.qexp(3)
        q + (-zeta4 - 1)*q^2 + O(q^3)
        sage: f.qexp(2)
        q + O(q^2)
        sage: f.qexp(1)
        O(q^1)
    """

    #def _repr_(self):
    #    A = self.ambient_module()
    #    return "Cuspidal subspace of dimension %s of Modular Forms space with character %s and weight %s over %s"%(self.dimension(), self.character(), self.weight(), self.base_ring())


