"""nodoctest
all.py -- much of sage is imported into this module, so you don't
          have to import everything individually.
"""

from __future__ import with_statement

###############################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
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
###############################################################################

# Error message that matches the SAGE/IPython defaults
quit = "Use Ctrl-D (i.e. EOF), %Exit, or %Quit to exit without confirmation."
exit = quit

import os, sys

if sys.version_info[:2] < (2, 4):
    print >>sys.stderr, "SAGE requires Python 2.4 or newer"
    sys.exit(1)

try:
    _l = '%s/local/lib'%os.environ['SAGE_ROOT']
    if os.environ.has_key('LD_LIBRARY_PATH'):
        if not _l in os.environ['LD_LIBRARY_PATH']:
            raise KeyError
        elif not _l in os.environ['DYLD_LIBRARY_PATH']:
            raise KeyError
    del _l
except KeyError:
     raise RuntimeError, "To use the SAGE libraries, set the environment variable SAGE_ROOT to the SAGE build directory and LD_LIBRARY_PATH to $SAGE_ROOT/local/lib"


###################################################################

from random              import *
from time                import sleep

from sage.interfaces.get_sigs import get_sigs
get_sigs()

from sage.rings.memory import pmem_malloc
pmem_malloc()

from sage.misc.all       import *         # takes a while

from sage.libs.all       import *
from sage.rings.all      import *
from sage.matrix.all     import *

pmem_malloc()

from sage.modules.all    import *
from sage.monoids.all    import *
from sage.algebras.all   import *
from sage.modular.all    import *
from sage.schemes.all    import *
from sage.graphs.all     import *
from sage.groups.all     import *
from sage.databases.all  import *
from sage.structure.all  import *
from sage.categories.all import *
from sage.sets.all       import *
from sage.probability.all import *
from sage.interfaces.all import *
from sage.functions.all  import *
#from sage.calculus.all   import *
from sage.server.all     import *
from sage.dsage.all      import *
import sage.tests.all as tests

from sage.crypto.all     import *

from sage.plot.all       import *
from sage.coding.all     import *
from sage.combinat.all   import *

from sage.lfunctions.all import *

from sage.geometry.all   import *

from sage.quadratic_forms.all import *

from sage.gsl.all        import *

from sage.games.all      import *

from copy import copy


###########################################################
#### WARNING:
# DO *not* import numpy / matplotlib / networkx here!!
# Each takes a surprisingly long time to initialize,
# and that initialization should be done more on-the-fly
# when they are first needed.
###########################################################

###################################################################

# maximize memory resources
try:
    import resource   # unix only...
    resource.setrlimit(resource.RLIMIT_AS, (-1,-1))
except:
    pass

# very useful 2-letter shortcuts
CC = ComplexField()
I = CC.gen(0)
QQ = RationalField()
RR = RealField()  # default real field
ZZ = IntegerRing()
# NOTE: QQ, RR, and ZZ are used by the pre-parser, and should not be
# overwritten by the user, unless they want to change the meaning of
# int and real in the interpreter (which is a potentially valid thing
# to do, and doesn't mess up anything else in the SAGE library).
# E.g., typing "int = ZZ" in the SAGE interpreter makes int literals
# acts as Python ints again.



# Some shorter shortcuts:
# Q = QQ
# Z = ZZ
# C = CC
#i = CC.gen(0)
true = True
false = False

oo = infinity
x = PolynomialRing(QQ,'x').gen()

# grab signal handling back from PARI or other C libraries
get_sigs()

from sage.misc.copying import license
copying = license
copyright = license

#_copyright = str(copyright)
#class __Copyright2:
#    def __repr__(self):
#        return _copyright + '\n\nCopyright (c) 2004-2006 William Stein.\nAll Rights Reserved.'
#
#copyright = __Copyright2()



def save_session(state, name='default_session', verbose=False):
    """
    Save all variables defined during this session (that can be saved)
    to the given filename.  The variables will be saved to a dictionary,
    which can be loaded using load(name).

    From the interpreter use \code{save_session(name)} (i.e., ignore
    that state argument).
    """
    S = set(show_identifiers(state))
    D = {}
    for k, x in state.items():
        try:
            if k in S:
                if type(x) == type(save_session):
                    raise TypeError, '%s is a function'%k
                _ = loads(dumps(x))
                if verbose:
                    print "Saving %s"%k
                D[k] = x
        except (IOError, TypeError), msg:
            if verbose:
                print "Not saving %s: %s"%(k, msg)
            pass
    save(D, name)

def load_session(state, name='default_session', verbose=True):
    """
    Load a saved session.

    From the interpreter use \code{load_session(name)} (i.e., ignore
    that state argument).

    This merges in all variables from a previously saved session.  It
    does not clear out the variables in the current sessions, unless
    they are overwritten.  You can thus merge multiple sessions, and
    don't loose all your current work when you use this command.
    """
    D = load(name)
    for k, x in D.items():
        if verbose:
            print "Loading %s"%k
        state[k] = x

def show_identifiers(state):
    """
    Returns a list of all variables you have defined during this session.

    In general, call using show_identifiers(locals()).

    If you are typing at the prompt you can type show_identifiers()
    alone on a line to get a list of identifiers.
    """
    return [k for k in state.keys() if not (k in iglob) and k[0] != '_']

iglob = set(['quit_sage', 'InteractiveShell', \
             'Out', 'In', 'help', 'view_all', \
             'os', 'iglob', 'variables', '__', '_dh', "_ii", '_i2',
             '_i1', '_iii', '_i', '_ih', '___', '_oh', '_', '__nonzero__', '__IP']).union(set(globals().keys()))

_cpu_time_ = cputime()
_wall_time_ = walltime()

def quit_sage(verbose=True):
    """
    If you use SAGE in library mode, you should call this function
    when your application quits.

    It makes sure any child processes are also killed, etc.
    """
    if verbose:
        t1 = cputime(_cpu_time_)
        t1m = int(t1/60); t1s=t1-t1m*60
        t2 = walltime(_wall_time_)
        t2m = int(t2/60); t2s=t2-t2m*60
        print "Exiting SAGE (CPU time %sm%.2fs, Wall time %sm%.2fs)."%(
               t1m,t1s,t2m,t2s)
    from sage.interfaces.quit import expect_quitall
    expect_quitall(verbose=verbose)
    from sage.misc.misc import delete_tmpfiles
    delete_tmpfiles()

def _quit_sage_(self):
    import sage.misc.preparser_ipython
    if sage.misc.preparser_ipython.interface != None:
        sage.misc.preparser_ipython.switch_interface('sage')
        self.exit_now = False
        return

    from IPython.genutils import ask_yes_no
    if self.rc.confirm_exit:
        if ask_yes_no('Do you really want to exit ([y]/n)?','y'):
            self.exit_now = True
    else:
        self.exit_now = True
    if self.exit_now:
        quit_sage()
        self.exit_now = True
    return self.exit_now

from IPython.iplib import InteractiveShell
InteractiveShell.exit = _quit_sage_

from sage.ext.interactive_constructors_c import inject_on, inject_off

from catalogue.all import new

import sage.ext.sig
sage.ext.sig.get_bad_sigs()
