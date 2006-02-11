"""
Frank Luebeck's tables of Conway polynomials over finite fields.
"""

#*****************************************************************************
#
#       SAGE: Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

import sage.misc.misc
import sage.databases.db   # very important that this be fully qualified
_CONWAYDATA = "%s/data/src/conway/conway_table.py.bz2"%sage.misc.misc.SAGE_ROOT


class ConwayPolynomials(sage.databases.db.Database):
    def __init__(self, read_only=True):
        """
        Initialize the database.

        INPUT:
            read_only -- bool (default: True), if True, then the
                         database is read_only and changes cannot be
                         commited to disk.
        """
        sage.databases.db.Database.__init__(self,
             name="conway_polynomials", read_only=read_only)

    def _init(self):
        if not os.path.exists(_CONWAYDATA):
            raise RuntimeError, "In order to initialize the database, the file %s must exist."%_CONWAYDATA
        os.system("cp %s ."%_CONWAYDATA)
        os.system("bunzip2 conway_table.py.bz2")
        from conway_table import conway_polynomials
        for X in conway_polynomials:
            (p, n, v) = tuple(X)
            if not self.has_key(p):
                self[p] = {}
            self[p][n] = v
            self[p] = self[p]  # so database knows it was changed
        os.system("bzip2 conway_table.py; rm conway_table.pyc; cp conway_table.py %s"%_CONWAYDATA)
        self.commit()

    def __repr__(self):
        return "Frank Luebeck's database of Conway polynomials"


    def polynomial(self, p, n):
        try:
            return self[int(p)][int(n)]
        except KeyError:
            raise RuntimeError, "Conway polynomial over F_%s of degree %s not in database."%(p,n)

    def has_polynomial(self, p, n):
        return self.has_key(int(p)) and self[int(p)].has_key(int(n))

    def primes(self):
        return self.keys()

    def degrees(self, p):
        p = int(p)
        if not p in self.primes():
            return []
        return self[p].keys()


