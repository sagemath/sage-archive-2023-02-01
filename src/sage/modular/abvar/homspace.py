"""
Spaces of homomorphisms between modular abelian varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.categories.homset import HomsetWithBase

import abvar as abelian_variety

class Homspace(HomsetWithBase):
    """
    A space of homomorphisms between two modular abelian varieties.
    """
    def __init__(self, domain, codomain, cat):
        """
        Create a homspace.

        INPUT:
            domain, codomain -- modular abelian varieties
            cat -- category

        EXAMPLES:
            sage: H = Hom(J0(11), J0(22)); H
            Space of homomorphisms from Jacobian of the modular curve associated to the congruence subgroup Gamma0(11) to Jacobian of the modular curve associated to the congruence subgroup Gamma0(22)
            sage: type(H)
            <class 'sage.modular.abvar.homspace.Homspace'>
            sage: H.homset_category()
            Category of modular abelian varieties over Rational Field
        """
        if not isinstance(domain, abelian_variety.ModularAbelianVariety):
            raise TypeError, "domain must be a modular abelian variety"
        if not isinstance(codomain, abelian_variety.ModularAbelianVariety):
            raise TypeError, "codomain must be a modular abelian variety"
        HomsetWithBase.__init__(self, domain, codomain, cat)

    def _repr_(self):
        """
        String representation of a modular abelian variety homspace.

        EXAMPLES:
            sage: J = J0(11)
            sage: End(J)._repr_()
            'Space of homomorphisms from Jacobian of the modular curve associated to the congruence subgroup Gamma0(11) to Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)'
        """
        return "Space of homomorphisms from %s to %s"%\
               (self.domain(), self.codomain())
