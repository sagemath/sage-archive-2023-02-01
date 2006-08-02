"""
Database of Modular Polynomials
"""

#######################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                          David Kohel
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#######################################################################


import bz2, os
import sage.misc.misc
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.multi_polynomial_ring import MPolynomialRing

DB_HOME = '%s/kohel'%sage.misc.misc.SAGE_DATA

def _dbz_to_integer_list(name):
    file = '%s/%s'%(DB_HOME, name)
    if not os.path.exists(file):
        raise RuntimeError, "Modular polynomial database file %s not available"%file
    data = bz2.decompress(open(file).read())
    data = "[[" + data.replace("\n","],[").replace(" ",",")[:-2] + "]"
    return eval(data)

def _pad_int(s,n):
    return "0"*(n-len(str(s))) + str(s)

class ModularPolynomialDatabase:
    def __init__(self):
        """
        Initialize the database.
        """
        pass

    def __repr__(self):
        return "Modular polynomial database"

    def __getitem__(self,level):
        model = "Cls"
        modpol = "PolMod/%s/pol.%s.dbz"%(model, _pad_int(level,3))
        try:
            coeff_list = _dbz_to_integer_list(modpol)
        except RuntimeError, msg:
            print msg
            raise RuntimeError, \
                  "No database entry for modular polynomial of level %s"%level
        P = MPolynomialRing(IntegerRing(),2)
        (X,Y) = P.gens()
        Phi = P(0)
        for k in range(len(coeff_list)):
            cff = coeff_list[k]
            i = cff[0]
            j = cff[1]
            if i == j:
                mon = X**i*Y**j
            else:
                mon = X**i*Y**j + X**j*Y**i
            Phi += Integer(cff[2])*mon
        return Phi
