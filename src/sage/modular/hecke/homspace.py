r"""
Hom spaces between objects of the category of hecke modules over a
given base ring.
"""

#*****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.categories.homset
import morphism
import module

def is_HeckeModuleHomspace(x):
    return isinstance(x, HeckeModuleHomspace)

class HeckeModuleHomspace(sage.categories.homset.HomsetWithBase):
    def __init__(self, X, Y):
        if not module.is_HeckeModule(X) or not module.is_HeckeModule(Y):
            raise TypeError, "X and Y must be Hecke modules"
        if X.base_ring() != Y.base_ring():
            raise TypeError, "X and Y must have the same base ring"
        sage.categories.homset.HomsetWithBase.__init__(self, X, Y, X.category())

    def __call__(self, A, name=''):
        try:
            A = A.hecke_module_morphism()
            if A.parent() == self:
                return A
            else:
                raise TypeError, "unable to coerce A to self"
        except AttributeError:
            pass
        return morphism.HeckeModuleMorphism_matrix(self, A, name)


