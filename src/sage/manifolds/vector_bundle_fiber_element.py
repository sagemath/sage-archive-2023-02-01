r"""
Vector Bundle Fiber Elements

The class :class:`VectorBundleFiberElement` implements vectors in the fiber of
a vector bundle.

AUTHORS:

- Michael Jung (2019): initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement

class VectorBundleFiberElement(FiniteRankFreeModuleElement):
    r"""

    """
    def __init__(self, parent, name=None, latex_name=None):
        r"""
        Construct a vector in the given fiber of a given vector bundle.

        TESTS::



        """
        FiniteRankFreeModuleElement.__init__(self, parent, name=name,
                                             latex_name=latex_name)
        # Extra data (with respect to FiniteRankFreeModuleElement):
        self._point = parent._point
        self._vector_bundle = parent._vector_bundle

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::



        """
        desc = "Vector {} in the fiber of {} at {}".format(self._name,
                                                     self._vector_bundle._name,
                                                     self._point)
        return desc