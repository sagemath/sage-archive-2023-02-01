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
        self.__f = f

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

