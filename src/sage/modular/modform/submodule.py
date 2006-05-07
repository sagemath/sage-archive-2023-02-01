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

class ModularForms_submodule(space.ModularFormsSpace):
    """
    A submodule of an ambient space of modular forms.
    """
    def __init__(self, ambient_space, vector_space):
        """
            ambient_space -- ModularFormsSpace
            submodule -- a vector submodule of the underlying vector space of the ambient space.
        """
        self.__ambient_space = ambient_space
        self.__vector_space = vector_space
        A = ambient_space
        ModularFormsSpace.__init__(self, A.group(), A.weight(), A.character(), A.base_field())

    def __repr__(self):
        return "ModularFormsSubmodule(%s,dim=%s)"%(self.ambient_space(), self.dimension())
        #return "Submodule of dimension %s of modular forms of weight %s on %s with character %s."%(
        #    self.dimension(), self.weight(), self.group(), self.character())

    def is_ambient(self):
        return False

    def ambient_space(self):
        return self.__ambient_space

    def change_ring(self):
        raise NotImplementedError, "Base change only currently implemented for ambient spaces."

    def vector_space(self):
        return self.__vector_space

    def dimension(self):
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.vector_space().dimension()
        return self.__dimension



