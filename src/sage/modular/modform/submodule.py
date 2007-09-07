"""
Submodules of spaces of modular forms

EXAMPLES:
    sage: M = ModularForms(Gamma1(13),2); M
    Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: M.eisenstein_subspace()
    Eisenstein subspace of dimension 11 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: M.cuspidal_subspace()
    Cuspidal subspace of dimension 2 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
    sage: M.new_subspace()
    Modular Forms subspace of dimension 13 of Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import space

import sage.modular.hecke.submodule

class ModularFormsSubmodule(space.ModularFormsSpace,
                            sage.modular.hecke.submodule.HeckeSubmodule):
    """
    A submodule of an ambient space of modular forms.
    """
    def __init__(self, ambient_module, submodule, dual=None, check=False):
        """
            ambient_module -- ModularFormsSpace
            submodule -- a submodule of the ambient space.
            dual_module -- (default: None) ignored
            check -- (default: False) whether to check that the
                     submodule is Hecke equivariant
        """
        A = ambient_module
        sage.modular.hecke.submodule.HeckeSubmodule.__init__(self, A, submodule, check=check)
        space.ModularFormsSpace.__init__(self, A.group(), A.weight(),
                                         A.character(), A.base_ring())

    def _repr_(self):
        return "Modular Forms subspace of dimension %s of %s"%(self.dimension(), self.ambient_module())

    def change_ring(self, base_ring):
        raise NotImplementedError, "Base change only currently implemented for ambient spaces."

    def _compute_coefficients(self, element, X):
        raise NotImplementedError

    def _compute_q_expansion_basis(self, prec):
        A = self.ambient_module()
        return [A._q_expansion(element = f.element(), prec=prec) for f in self.basis()]


# TODO
class ModularFormsSubmoduleWithBasis(ModularFormsSubmodule):
    pass





