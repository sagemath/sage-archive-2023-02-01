"""
TODO
"""


###########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################


from sage.rings.all  import QQ, ZZ

from sage.modular.modform.element import Newform

from abvar import ModularAbelianVariety_modsym_abstract

class ModularAbelianVariety_newform(ModularAbelianVariety_modsym_abstract):
    """
    A modular abelian variety attached to a specific newform.
    """
    def __init__(self, f):
        """
        Create the modular abelian variety $A_f$ attached to the
        newform $f$.

        INPUT:
            f -- a newform

        EXAMPLES:
            sage: S = CuspForms(37)
            sage: f = sage.modular.modform.element.Newform(S, S.modular_symbols(1)[1])
            sage: f.abelian_variety()
            Modular abelian variety attached to the newform q + q^3 - 2*q^4 + O(q^6)
        """
        if not isinstance(f, Newform):
            raise TypeError, "f must be a newform"
        self.__f = f
        ModularAbelianVariety_modsym_abstract.__init__(self, QQ)

    def _modular_symbols(self):
        return self.__f.modular_symbols()

    def newform(self):
        """
        Return the newform that this modular abelian variety is attached to.

        EXAMPLES:
            sage: S = CuspForms(37)
            sage: f = sage.modular.modform.element.Newform(S, S.modular_symbols(1)[1])
            sage: A = f.abelian_variety()
            sage: A.newform()
            q + q^3 - 2*q^4 + O(q^6)
            sage: A.newform() is f
            True
        """
        return self.__f

    def _repr_(self):
        """
        String representation of this modular abelian variety.

        EXAMPLES:
        """
        return "Modular abelian variety attached to the newform %s"%self.newform()

