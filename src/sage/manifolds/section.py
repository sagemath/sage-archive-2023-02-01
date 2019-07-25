r"""
Sections

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
from sage.structure.element import ModuleElement
from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement


class Section(ModuleElement):
    pass

#******************************************************************************

class SectionParal(FiniteRankFreeModuleElement, Section):
    pass
