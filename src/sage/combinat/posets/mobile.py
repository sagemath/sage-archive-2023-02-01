# ****************************************************************************
#       Copyright (C) 2020 Stefan Grosser <stefan.grosser1@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.posets.d_complete import DCompletePoset
from sage.misc.lazy_attribute import lazy_attribute
#from .linear_extensions import LinearExtensionsOfForest

class MobilePoset(FinitePoset):
    r"""
    Mobile posets are an extension of d-complete posets which permit a determinant
    formula for counting linear extensions. 

    """
    
    #_lin_ext_type = LinearExtensionsOfForest
    _desc = 'Finite mobile poset'
    
    def _compute_mobile_structure(self):
        pass
            
            
    
    
    
