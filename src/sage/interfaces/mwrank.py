r"""
Interface to mwrank
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import os, weakref
from expect import Expect

instances={}
def Mwrank(options="", server=None, server_tmpdir=None):
    """
    Create and return an mwrank interpreter, with given options.

    INPUT:

    -  ``options`` - string; passed when starting mwrank.
       The format is::

       -h       help            prints this info and quits
       -q       quiet           turns OFF banner display and prompt
       -v n     verbosity       sets verbosity to n (default=1)
       -o       PARI/GP output  turns ON extra PARI/GP short output (default is OFF)
       -p n     precision       sets precision to n decimals (default=15)
       -b n     quartic bound   bound on quartic point search (default=10)
       -x n     n aux           number of aux primes used for sieving (default=6)
       -l       list            turns ON listing of points (default ON unless v=0)
       -s       selmer_only     if set, computes Selmer rank only (default: not set)
       -d       skip_2nd_descent        if set, skips the second descent for curves with 2-torsion (default: not set)
       -S n     sat_bd          upper bound on saturation primes (default=100, -1 for automatic)

    .. warning:

       Do not use the option "-q" which turns off the prompt.

    EXAMPLES::

        sage: M = Mwrank('-v 0 -l')
        sage: print M('0 0 1 -1 0')
        Curve [0,0,1,-1,0] :    Rank = 1
        Generator 1 is [0:-1:1]; height 0.0511114082399688
        Regulator = 0.0511114082399688
    """
    global instances
    try:
        X = instances[options]()
        if X:
            return X
    except KeyError:
        pass
    X = Mwrank_class(options, server=server,server_tmpdir=server_tmpdir)
    instances[options] = weakref.ref(X)
    return X


class Mwrank_class(Expect):
    """
    Interface to the Mwrank interpreter.
    """
    def __init__(self, options="", server=None,server_tmpdir=None):
        """
        INPUT:


        -  ``options`` - string; passed when starting mwrank.
           The format is::

           -h       help            prints this info and quits
           -q       quiet           turns OFF banner display and prompt
           -v n     verbosity       sets verbosity to n (default=1)
           -o       PARI/GP output  turns ON extra PARI/GP short output (default is OFF)
           -p n     precision       sets precision to n decimals (default=15)
           -b n     quartic bound   bound on quartic point search (default=10)
           -x n     n aux           number of aux primes used for sieving (default=6)
           -l       list            turns ON listing of points (default ON unless v=0)
           -s       selmer_only     if set, computes Selmer rank only (default: not set)
           -d       skip_2nd_descent        if set, skips the second descent for curves with 2-torsion (default: not set)
           -S n     sat_bd          upper bound on saturation primes (default=100, -1 for automatic)

    .. warning:

       Do not use the option "-q" which turns off the prompt.


        .. note::

           Normally instances of this class would be created by
           calling the global function :meth:`Mwrank`.

        EXAMPLES::

            sage: from sage.interfaces.mwrank import Mwrank_class
            sage: M = Mwrank_class('-v 0 -l')
            sage: M('0 -1 1 0 0')
            'Curve [0,-1,1,0,0] :   Rank = 0\n\n\nRegulator = 1\n'

            sage: from sage.interfaces.mwrank import Mwrank_class
            sage: TestSuite(Mwrank_class).run()
        """
        Expect.__init__(self,
                        name = 'mwrank',
                        prompt = 'Enter curve: ',
                        command = "mwrank %s"%options,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        maxread = 10000,
                        restart_on_ctrlc = True,
                        verbose_start = False)

    def __getattr__(self, attrname):
        """
        Standard function to return an attribute.

        EXAMPLES::

            sage: mwrank.zzz
            Traceback (most recent call last):
            ...
            AttributeError
        """
        raise AttributeError

    def __reduce__(self):
        """
        EXAMPLES::

            sage: mwrank.__reduce__()
            (<function _reduce_load_mwrank at 0x...>, ())
        """

        return _reduce_load_mwrank, tuple([])

    def __call__(self, cmd):
        """
        Standard call function.

        EXAMPLES::

            sage: s = mwrank("0 0 0 0 1"); print s
            Curve [0,0,0,0,1] :
            ...
            Rank = 0
            ...
        """
        return self.eval(str(cmd))

    def eval(self, *args, **kwds):
        """
        Send a line of input to mwrank, then when it finishes return
        everything that mwrank output.

        .. note::

           If a RuntimeError exception is raised, then the mwrank
           interface is restarted and the command is retried once.

        EXAMPLES::

            sage: mwrank.eval('12 3 4 5 6')
            'Curve [12,3,4,5,6] :...'
        """
        if self._expect is not None and not self._expect.isalive():
            # if mwrank is interrupted twice in rapid succession,
            # then it doesn't restart correctly, and we're left with:
            #   "RuntimeError: [Errno 9] Bad file descriptor"
            # Doing _start again fixes that always. See trac #5157.
            self._start()
        return Expect.eval(self, *args, **kwds).replace('\t','   ')

    def console(self):
        """
        Start the mwrank console.

        EXAMPLE::

            sage: mwrank.console() # not tested: expects console input
            Program mwrank: ...

        """
        mwrank_console()

    def quit(self, verbose=False):
        """
        Quit the mwrank process using kill -9 (so exit doesn't dump core, etc.).

        INPUT:

        - ``verbose`` -- ignored

        EXAMPLES::

            sage: m = Mwrank()
            sage: e = m('1 2 3 4 5')
            sage: m.quit()
        """
        if self._expect is None: return
        try:
            os.kill(self._expect.pid, 9)
        except OSError:
            pass
        self._expect = None


# An instance
mwrank = Mwrank()

def _reduce_load_mwrank():
    """
    Return the standard mwrank instance

    EXAMPLES::

        sage: from sage.interfaces.mwrank import _reduce_load_mwrank
        sage: _reduce_load_mwrank()
        Mwrank
    """
    return mwrank

import os
def mwrank_console():
    """
    Start the mwrank console.

    EXAMPLE::

        sage: mwrank_console() # not tested: expects console input
        Program mwrank: ...
    """
    os.system('mwrank')

