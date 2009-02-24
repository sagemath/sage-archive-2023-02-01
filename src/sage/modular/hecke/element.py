"""
Elements of Hecke modules

AUTHORS:

- William Stein
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
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

import sage.modules.module_element

def is_HeckeModuleElement(x):
    return isinstance(x, HeckeModuleElement)

class HeckeModuleElement(sage.modules.module_element.ModuleElement):
    """
    Element of a Hecke module.
    """
    def __init__(self, parent, x=None):
        """
        INPUT:


        -  ``parent`` - a Hecke module

        -  ``x`` - element of the free module associated to
           parent
        """
        sage.modules.module_element.ModuleElement.__init__(self, parent)
        if not x is None:
            self.__element = x

    def _repr_(self):
        return str(self.element())

    def _compute_element(self):
        raise NotImplementedError

    def element(self):
        try:
            return self.__element
        except AttributeError:
            self.__element = self._compute_element()
        return self.__element

    def ambient_module(self):
        return self.parent().ambient_module()

    def _lmul_(self, x):
        return self.parent()(self.element()*x)

    def _rmul_(self, x):
        return self.parent()(x * self.element())

    def _neg_(self):
        return self.parent()(-self.element())

    def _pos_(self):
        return self

    def _sub_(self, right):
        return self.parent()(self.element() - right.element())
