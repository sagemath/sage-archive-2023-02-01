"""
Morphisms between modular abelian varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

import sage.categories.morphism

class Morphism_abstract(sage.categories.morphism.Morphism):
    """
    A morphism between modular abelian varieties.
    """
    def __init__(self, parent):
        """
        Create a morphism between modular abelian varieties.

        INPUT:
             parent -- a homset

        EXAMPLES:
            sage: t = J0(11).hecke_operator(2)
            sage: from sage.modular.abvar.morphism import Morphism_abstract
            sage: isinstance(t, Morphism_abstract)
            True
        """
        sage.categories.morphism.Morphism.__init__(self, parent)

    def matrix(self):
        raise NotImplementedError

class Morphism(Morphism_abstract):
    def __init__(self, parent, matrix):
        self.__matrix = matrix
        Morphism_abstract.__init__(self, parent)

    def matrix(self):
        return self.__matrix



