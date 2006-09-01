r"""
Mutability
"""

##########################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
##########################################################################

from sage.structure.mutability_pyx import Mutability

def __reduce__Mutability(x):
    return Mutability(x)
