"""
TODO
"""


###########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################


from sage.rings.all  import QQ, ZZ

import abvar

class ModularAbelianVariety_newform(abvar.ModularAbelianVariety):
    """
    A modular abelian variety attached to a specific newform.
    """
    def __init__(self, f):
        """
        Create the modular abelian variety $A_f$ attached to the
        newform $f$.

        INPUT:
            f -- a newform

        This class is not used anywhere else yet!

        EXAMPLES:
            sage: from sage.modular.abvar.abvar_newform import ModularAbelianVariety_newform
            sage: f = CuspForms(11).0
            sage: Af = ModularAbelianVariety_newform(f); Af
            Modular abelian variety attached to the newform q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
        """
        # todo: check that f is is newform -- need a newform class.
        self.__f = f
        levels = (f.group(),)
        L = f.modular_symbols(sign=0).free_module()
        lattice = L.intersection(ZZ**L.degree())
        abvar.ModularAbelianVariety(self, levels, lattice, QQ)

    def newform(self):
        """
        Return the newform that this modular abelian variety is attached to.

        EXAMPLES:
            sage: from sage.modular.abvar.abvar_newform import ModularAbelianVariety_newform
            sage: f = CuspForms(11).0
            sage: Af = ModularAbelianVariety_newform(f)
            sage: Af.newform()
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
        """
        return self.__f

    def _repr_(self):
        """
        String representation of this modular abelian variety.

        EXAMPLES:
            sage: from sage.modular.abvar.abvar_newform import ModularAbelianVariety_newform
            sage: f = CuspForms(11).0
            sage: Af = ModularAbelianVariety_newform(f)
            sage: Af._repr_()
            'Modular abelian variety attached to the newform q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)'
        """
        return "Modular abelian variety attached to the newform %s"%self.newform()

