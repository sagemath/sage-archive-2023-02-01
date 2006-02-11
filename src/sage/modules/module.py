"""
PYREX: sage.ext.module
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

from sage.ext.module import Module

def is_Module(x):
    """
    Return True if x is a module.

    EXAMPLES:
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

    EXAMPLES:
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


