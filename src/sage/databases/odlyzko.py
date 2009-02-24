"""
Tables of zeros of the Riemann-Zeta function.
"""

#*****************************************************************************
#
#       SAGE: Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

import sage.misc.db as db
import sage.misc.misc as misc

PATH = "%s/data/"%misc.SAGE_ROOT

def zeta_zeros():
    r"""
    List of the imaginary parts of the first 100,000 nontrivial zeros
    of the Riemann zeta function. Andrew Odlyzko computed these to
    precision within `3\cdot 10^{-9}`.

    In order to use ``zeta_zeros()``, you will need to
    install the optional Odlyzko database package: ``sage -i
    database_odlyzko_zeta``. You can see a list of all
    available optional packages with ``sage -optional``.

    REFERENCES:

    - http://www.dtc.umn.edu/~odlyzko/zeta_tables/

    EXAMPLES:

    The following example prints the imaginary part of the 13th
    nontrivial zero of the Riemann zeta function. Note that only the
    first 9 digits after the decimal come from the database. Subsequent
    digits are the result of the inherent imprecision of a binary
    representation of decimal numbers.

    ::

        sage: zz = zeta_zeros() # optional
        sage: zz[12]            # optional
        59.347044003000001
    """
    path = "%s/odlyzko"%PATH
    file = "%s/zeros1"%path
    if os.path.exists(file+".pickle"):
        misc.verbose("Loading Odlyzko database from " + file + ".pickle")
        return db.load(file+".pickle")
    misc.verbose("Creating Odlyzko Database.")
    F = [eval(x) for x in open(file).read().split()]
    db.save(F, file+".pickle")
    return F

