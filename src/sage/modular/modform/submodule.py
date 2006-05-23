"""
Submodules of spaces of modular forms
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@ucsd.edu>
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
    def __init__(self, ambient_space, submodule):
        """
            ambient_space -- ModularFormsSpace
            submodule -- a submodule of the ambient space.
        """
        A = ambient_space
        sage.modular.hecke.submodule.HeckeSubmodule.__init__(self, A, submodule)
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


class ModularFormsSubmoduleWithBasis(ModularFormsSubmodule):
    pass





