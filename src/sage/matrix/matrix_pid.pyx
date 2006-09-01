"""
Matrices over a PID
"""

########################################################################
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
########################################################################


cimport matrix_domain
import matrix_domain

cdef class Matrix_pid(matrix_domain.Matrix_domain):
    def __init__(self, parent):
        matrix_domain.Matrix_domain.__init__(self, parent)
