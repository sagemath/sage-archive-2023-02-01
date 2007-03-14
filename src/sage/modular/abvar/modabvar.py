"""
Modular abelian varieties
"""

########################################################################
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
########################################################################

#from sage.modular.hecke.module import HeckeModule_generic

from sage.structure.parent_base import ParentWithBase

class ModularAbelianVariety(ParentWithBase):
    def __init__(self, base_ring, level):
        """
        INPUT:
            base_ring -- the field over which the abelian
                         variety is defined
            level -- an integer N such that this abelian variety
                     is a quotient of J_1(N).
        """
        ParentWithBase.__init__(self, base_ring)
        self._base_ring = base_ring
        self._level = level

    def base_ring(self):
        return self._base_ring

    def level(self):
        return self._level



