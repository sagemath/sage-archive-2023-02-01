##############################################################################
#
#  DSAGE: Distributed SAGE
#
#       Copyright (C) 2006, 2007 Yi Qiang <yqiang@gmail.com>
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
#
##############################################################################
"""
Miscellaneous helper methods.

"""

import random

def random_str(length=500):
    """
    Generates a random string.

    INPUT:
    length -- the length of the string

    """
    s = [chr(i) for i in [random.randint(65, 123) for n in range(length)]]

    return ''.join(s)
