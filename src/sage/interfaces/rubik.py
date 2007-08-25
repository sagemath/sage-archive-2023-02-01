r"""nodoctest

Interface to two Rubik's cube solvers.

The first is by Michael Reid, and tries to find an optimal solution given
the cube's state, and may take a long time.
See http://www.math.ucf.edu/~reid/Rubik/optimal_solver.html

The second is by Eric Dietz, and uses a standard (?) algorithm to
solve the cube one level at a time. It is extremly fast, but often
returns a far from optimal solution.
See http://wrongway.org/?rubiksource


TODO:
    I'm sure there's algorithms and/or implementations that provide
    a happy medium between the above two. Find or implement and interface.


AUTHOR:
   -- Optimal was written by Michael Reid <reid@math.ucf.edu> (2004)
   -- Cubex was written by Eric Dietz <root@wrongway.org> (2003)
   -- Initial interface by Robert Bradshaw (2007-08)
"""

########################################################################
#   Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, pexpect
import cleaner

from sage.groups.perm_gps.cubegroup import *



# Can't seem to find consistancy in letter ordering
# between us and them... These are copied from the source.
optimal_solver_tokens = ["UF", "UR", "UB", "UL", \
                        "DF", "DR", "DB", "DL", \
                        "FR", "FL", "BR", "BL", \
                        "FU", "RU", "BU", "LU", \
                        "FD", "RD", "BD", "LD", \
                        "RF", "LF", "RB", "LB", \
                        "UFR", "URB", "UBL", "ULF", \
                        "DRF", "DFL", "DLB", "DBR", \
                        "FRU", "RBU", "BLU", "LFU", \
                        "RFD", "FLD", "LBD", "BRD", \
                        "RUF", "BUR", "LUB", "FUL", \
                        "FDR", "LDF", "BDL", "RDB"]

# The imput format.
optimal_solver_format = "UF UR UB UL DF DR DB DL FR FL BR BL UFR URB UBL ULF DRF DFL DLB DBR"

class SingNot:
    """
    This class is to resolve difference between various Singmaster notation.
    Case is ignored, and the second and third letters may be swapped.

    EXAMPLE:
        sage: SingNot("acb") == SingNot("ACB")
        True
        sage: SingNot("acb") == SingNot("bca")
        False
    """
    def __init__(self, s):
        self.rep = s
        self.cannonical = (s[0] + "".join(sorted(list(s[1:])))).lower()
    def __eq__(self, other):
        return isinstance(other, SingNot) and other.cannonical == self.cannonical
    def __repr__(self):
        return self.rep
    def __hash__(self):
        return hash(self.cannonical)

# This is our list
singmaster_list = [''] + [SingNot(index2singmaster(i+1)) for i in range(48)]; singmaster_list

class Optimal:
    """
    Interface to Michael Reid's optimal Rubik's Cube solver.
    """
    __cmd = "/Users/robert/Desktop/optimal/optimal.out"

    def __init__(self, verbose=False, wait=True):
        self.verbose = verbose
        self.start()
        if wait:
            print "Initializing tables..."
            self.ready()
            print "Done."

    def start(self):
        child = pexpect.spawn(self.__cmd)
        cleaner.cleaner(child.pid, self.__cmd)
        child.timeout = None
        self.child = child
        self._ready = False

    def stop(self):
        if child:
            self.child.sendline(chr(3)) # send ctrl-c
            self.child.sendline(chr(4)) # send ctrl-d
            self.child.close(True)
            self.child = None

    def ready(self):
        if not self._ready:
            self.child.expect('enter cube')
            self._ready = True

    def solve(self, facets):
        """
        TODO: Let it keep searching once it found a solution?
        """
        self.ready()
        self.child.sendline(self.format_cube(facets))
        self.child.expect(r"([LRUDBF' ]+)\s+\((\d+)q\*?, (\d+)f\*?\)")
        self.child.sendline(chr(3)) # send ctrl-c
        return self.child.match.groups()[0].strip()

    def format_cube(self, facets):
        L = []
        optimal_solver_list = [SingNot(x) for x in optimal_solver_tokens]
        for f in optimal_solver_format.split(" "):
            ix = facets[singmaster_list.index(SingNot(f))-1]
            facet = singmaster_list[ix]
            L.append(optimal_solver_list[optimal_solver_list.index(facet)])
        return " ".join([str(f) for f in L])




