r"""
Finite Combinatorial Classes
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from combinat import CombinatorialClass

def FiniteCombinatorialClass(l):
    """
    Returns the combinatorial class with elements in l.
    """
    if not isinstance(l, list):
        l = list(l)
    return FiniteCombinatorialClass_l(l)

class FiniteCombinatorialClass_l(CombinatorialClass):
    def __init__(self, l):
        self.l = l

    def __repr__(self):
        return "Combinatorial class with elements in %s"%self.l

    def __contains__(self, x):
        return x in self.l

    def list(self):
        return self.l

    def count(self):
        return len(self.l)

