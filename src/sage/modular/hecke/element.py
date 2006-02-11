"""
Elements of Hecke modules

AUTHOR:
    -- William Stein
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
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

import sage.modules.module_element

def is_HeckeModuleElement(x):
    return isinstance(x, HeckeModuleElement)

class HeckeModuleElement(sage.modules.module_element.ModuleElement):
    """
    Element of a Hecke module.
    """
    def __init__(self, parent, x):
        """
        INPUT:
            parent -- a Hecke module
            x -- element of the free module associated to parent
        """
        sage.modules.module_element.ModuleElement.__init__(self, parent)
        self.__element = x

    def _repr_(self):
        return str(self.__element)

    def element(self):
        return self.__element

