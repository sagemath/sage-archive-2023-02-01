"""
The module of supersingular points
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
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


import sage.modular.hecke.all as hecke
import sage.rings.all as rings

class SupersingularModule(hecke.HeckeModule_free_module):
    def __init__(self, level=1, prime=2, base_ring=rings.RationalField()):
        self.__level = level
        self.__prime = prime
        self.__base_ring = base_ring
        hecke.HeckeModule_free_module.__init__(self, base_ring, level, weight=2)

    def __repr__(self):
        return "Module of supersingular points on X_0(%s)/F_%s over %s",\
               (self.__level, self.__prime, self.__base_ring)

    def base_ring(self):
        return self.__base_ring

    def dimension(self):
        raise NotImplementedError

    def hecke_matrix(self, n):
        """
        Matrix of Hecke operator on module of supersingular points.
        """
        raise NotImplementedError

    def level(self):
        return self.__level

    def prime(self):
        return self.__prime

    def weight(self):
        return 2

