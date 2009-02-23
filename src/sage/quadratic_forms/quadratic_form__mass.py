
######################################################
## Routines to compute the mass of a quadratic form ##
######################################################

## Import all general mass finding routines
from sage.quadratic_forms.quadratic_form__mass__Siegel_densities import *
from sage.quadratic_forms.quadratic_form__mass__Conway_Sloane_masses import *


###################################################


def shimura_mass__maximal(self, p):
    """
    Use Shimuras exact mass formula to compute the mass of a maximal
    quadratic lattice. This works for any totally real number field,
    but has a small technical restriction when n is odd.
    """
    pass


def hanke_mass__maximal(self, p):
    """
    Use Shimuras exact mass formula to compute the mass of a maximal
    quadratic lattice.  This works for any totally real number field,
    but has a small technical restriction when n is odd.
    """
    pass


def GHY_mass_maximal(self, p):
    """
    Use the GHY formula to compute the mass of a (maximal?) quadratic
    lattice. This works for any number field.

    See [GHY, Prop 7.4 and 7.5, p121] and [GY, Thrm 10.20, p25].
    """
    pass
