"""
Morphisms of free modules
"""

# A matrix morphism is a morphism that is defined by multiplication by a
# matrix.  Elements of domain must either have a method "vector()" that
# returns a vector that the defining matrix can hit from the left, or
# be coercible into vector space of appropriate dimension.

import sage.matrix.all as matrix
import sage.misc.misc as misc
import sage.modules.free_module as free_module
import matrix_morphism

import free_module_homspace

def is_FreeModuleMorphism(x):
    return isinstance(x, FreeModuleMorphism)

class FreeModuleMorphism(matrix_morphism.MatrixMorphism):
    def __init__(self, parent, A):
        """
        INPUT:


        -  ``parent`` - a homspace in a (sub) category of free
           modules

        -  ``A`` - matrix
        """
        if not free_module_homspace.is_FreeModuleHomspace(parent):
            raise TypeError, "parent (=%s) must be a free module hom space"%parent
        if isinstance(A, matrix_morphism.MatrixMorphism):
            A = A.matrix()
        A = parent._matrix_space()(A)
        matrix_morphism.MatrixMorphism.__init__(self, parent, A)

    def __repr__(self):
        if max(self.matrix().nrows(),self.matrix().ncols()) > 5:
            mat = "(not printing %s x %s matrix)"%(self.matrix().nrows(), self.matrix().ncols())
        else:
            mat = str(self.matrix())
        return "Free module morphism defined by the matrix\n%s\nDomain: %s\nCodomain: %s"%(\
            mat, misc.strunc(self.domain()), misc.strunc(self.codomain()))

