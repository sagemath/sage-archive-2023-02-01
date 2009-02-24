"""
Matrix Group Homsets

AUTHORS:

- William Stein (2006-05-07): initial version
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.groups.group_homset import GroupHomset_generic
from sage.categories.homset import HomsetWithBase
from sage.categories.category_types import Groups
GROUPS = Groups()

import matrix_group_morphism

def is_MatrixGroupHomset(x):
    return isinstance(x, MatrixGroupHomset)

class MatrixGroupHomset(GroupHomset_generic):
    def __init__(self, G, H):
        HomsetWithBase.__init__(self, G, H, GROUPS, G.base_ring())

    def __call__(self, im_gens, check=True):
        """
        EXAMPLES:
        """
        try:
            return matrix_group_morphism.MatrixGroupMorphism_im_gens(self,
                                                            im_gens, check=check)
        except (NotImplementedError, ValueError), err:
            raise TypeError, "images (=%s) do not define a valid homomorphism"%im_gens


