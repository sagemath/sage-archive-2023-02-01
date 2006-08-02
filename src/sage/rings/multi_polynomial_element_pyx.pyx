"""
Multivariate Polynomials (Pyrex Part)

AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-07-31)
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
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

def mpoly_repr_pyx(self):
    """
    Converts a polynomial over an aribitrary field
    to a string readable by Singular et al.

    This function is much faster than the native str()
    approach.
    """
    cdef int i
    cdef int n
    cdef int zero

    ret = []
    zero = 0
    vars = self.parent().variable_names()
    n = len(vars)

    poly = self.element()
    d = poly._PolyDict__repn
    for e,c in d.iteritems():
        ret = ret + ["(",str(c),")","*"]
        for i from 0 <= i < n:
            if e[i] != zero:
                ret = ret + [vars[i],"^",str(e[i]),"*"]
        ret.pop() #remove last "*"
        ret.append("+")
    try:
        ret.pop() #remove last "+"
    except IndexError:
        return "0"
    return "".join(ret)
