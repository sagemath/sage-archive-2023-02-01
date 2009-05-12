"""
Abstract base class for modules
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import random

cdef class Module(sage.structure.parent_gens.ParentWithAdditiveAbelianGens):
    """
    Generic module class.
    """
    def category(self):
        """
        Return the category to which this module belongs.
        """
        import sage.categories.all
        return sage.categories.all.Modules(self.base_ring())

    def endomorphism_ring(self):
        """
        Return the endomorphism ring of this module in its category.
        """
        try:
            return self.__endomorphism_ring
        except AttributeError:
            import sage.categories.all
            E = sage.categories.all.End(self)
            self.__endomorphism_ring = E
            return E


    def is_atomic_repr(self):
        """
        True if the elements have atomic string representations, in the
        sense that they print if they print at s, then -s means the
        negative of s. For example, integers are atomic but polynomials are
        not.
        """
        return False

    def __hash__(self):
        return hash(self.__repr__())



def is_Module(x):
    """
    Return True if x is a module.

    EXAMPLES::

        sage: from sage.modules.module import is_Module
        sage: M = FreeModule(RationalField(),30)
        sage: is_Module(M)
        True
        sage: is_Module(10)
        False
    """
    return isinstance(x, Module)

def is_VectorSpace(x):
    """
    Return True if x is a vector space.

    EXAMPLES::

        sage: from sage.modules.module import is_Module, is_VectorSpace
        sage: M = FreeModule(RationalField(),30)
        sage: is_VectorSpace(M)
        True
        sage: M = FreeModule(IntegerRing(),30)
        sage: is_Module(M)
        True
        sage: is_VectorSpace(M)
        False
    """
    import sage.modules.free_module
    return isinstance(x, sage.modules.free_module.FreeModule_generic_field)


