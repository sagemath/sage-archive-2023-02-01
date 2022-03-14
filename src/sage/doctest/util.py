"""
Utility functions

This module contains various utility functions and classes used in doctesting.

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.
"""

# ****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from os import sysconf, times

from sage.misc.misc import walltime


def count_noun(number, noun, plural=None, pad_number=False, pad_noun=False):
    """
    EXAMPLES::

        sage: from sage.doctest.util import count_noun
        sage: count_noun(1, "apple")
        '1 apple'
        sage: count_noun(1, "apple", pad_noun=True)
        '1 apple '
        sage: count_noun(1, "apple", pad_number=3)
        '  1 apple'
        sage: count_noun(2, "orange")
        '2 oranges'
        sage: count_noun(3, "peach", "peaches")
        '3 peaches'
        sage: count_noun(1, "peach", plural="peaches", pad_noun=True)
        '1 peach  '
    """
    if plural is None:
        plural = noun + "s"
    if pad_noun:
        # We assume that the plural is never shorter than the noun....
        pad_noun = " " * (len(plural) - len(noun))
    else:
        pad_noun = ""
    if pad_number:
        number_str = ("%%%sd"%pad_number)%number
    else:
        number_str = "%d"%number
    if number == 1:
        return "%s %s%s"%(number_str, noun, pad_noun)
    else:
        return "%s %s"%(number_str, plural)


def dict_difference(self, other):
    """
    Return a dict with all key-value pairs occurring in ``self`` but not
    in ``other``.

    EXAMPLES::

        sage: from sage.doctest.util import dict_difference
        sage: d1 = {1: 'a', 2: 'b', 3: 'c'}
        sage: d2 = {1: 'a', 2: 'x', 4: 'c'}
        sage: dict_difference(d2, d1)
        {2: 'x', 4: 'c'}

    ::

        sage: from sage.doctest.control import DocTestDefaults
        sage: D1 = DocTestDefaults()
        sage: D2 = DocTestDefaults(foobar="hello", timeout=100)
        sage: dict_difference(D2.__dict__, D1.__dict__)
        {'foobar': 'hello', 'timeout': 100}
    """
    D = dict()
    for k, v in self.items():
        try:
            if other[k] == v:
                continue
        except KeyError:
            pass
        D[k] = v
    return D


class Timer:
    """
    A simple timer.

    EXAMPLES::

        sage: from sage.doctest.util import Timer
        sage: Timer()
        {}
        sage: TestSuite(Timer()).run()
    """

    def _pid_cpu_seconds(self, pid):
        r"""
        Parse the ``/proc`` filesystem to get the cputime of the given
        ``pid``.

        This also includes the times for child processes **that have
        been ``wait()``ed for and terminated**. Specifically, pexpect
        processes DO NOT fall into that category.
        """
        with open(f"/proc/{pid}/stat", "r") as statfile:
            stats = statfile.read().split()

        # man 5 proc (linux)
        cputicks = sum( float(s) for s in stats[13:17] )

        hertz = sysconf("SC_CLK_TCK")
        return (cputicks / hertz)


    def _quick_cputime(self):
        r"""
        A fast and replacement for ``sage.misc.misc.cputime``.

        This is a "reliable" replacement (on Linux/BSD) that takes
        subprocesses (particularly pexpect interfaces) into
        account. The ``cputime()`` function from the ``misc`` module
        can be passed ``subprocesses=True``, but this has a few
        faults; mainly that it relies on each pexpect interface to
        implement its own ``cputime()`` function. And most of our
        pexpect interfaces either don't implement one, or implement
        one in a way that requires the subprocess (being pexpected) to
        be in perfect working condition -- that will often not be the
        case at the end of a doctest line.

        OUTPUT:

        A float measuring the cputime in seconds of the sage process
        and all its subprocesses.
        """
        # Start by using os.times() to get the cputime for sage itself
        # and any subprocesses that have been wait()ed for and that
        # have terminated.
        cputime = sum( times()[:4] )

        # Now try to get the times for any pexpect interfaces, since
        # they do not fall into the category above.
        from sage.interfaces.quit import expect_objects
        for s in expect_objects:
            S = s()
            if S and S.is_running():
                try:
                    cputime += self._pid_cpu_seconds(S.pid())
                except (ArithmeticError, LookupError, OSError,
                        TypeError, ValueError):
                    # This will fail anywhere but linux/BSD, but
                    # there's no good cross-platform way to get the
                    # cputimes from pexpect interfaces without
                    # totally mucking up the doctests.
                    pass

        return cputime


    def start(self):
        """
        Start the timer.

        Can be called multiple times to reset the timer.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer().start()
            {'cputime': ..., 'walltime': ...}
        """
        self.cputime = self._quick_cputime()
        self.walltime = walltime()
        return self

    def stop(self):
        """
        Stops the timer, recording the time that has passed since it
        was started.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: import time
            sage: timer = Timer().start()
            sage: time.sleep(float(0.5))
            sage: timer.stop()
            {'cputime': ..., 'walltime': ...}
        """
        self.cputime = self._quick_cputime() - self.cputime
        self.walltime = walltime(self.walltime)
        return self

    def annotate(self, object):
        """
        Annotates the given object with the cputime and walltime
        stored in this timer.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer().start().annotate(EllipticCurve)
            sage: EllipticCurve.cputime # random
            2.817255
            sage: EllipticCurve.walltime # random
            1332649288.410404
        """
        object.cputime = self.cputime
        object.walltime = self.walltime

    def __repr__(self):
        """
        String representation is from the dictionary.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: repr(Timer().start()) # indirect doctest
            "{'cputime': ..., 'walltime': ...}"
        """
        return str(self)

    def __str__(self):
        """
        String representation is from the dictionary.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: str(Timer().start()) # indirect doctest
            "{'cputime': ..., 'walltime': ...}"
        """
        return str(self.__dict__)

    def __eq__(self, other):
        """
        Comparison.

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer() == Timer()
            True
            sage: t = Timer().start()
            sage: loads(dumps(t)) == t
            True
        """
        if not isinstance(other, Timer):
            return False
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Test for non-equality

        EXAMPLES::

            sage: from sage.doctest.util import Timer
            sage: Timer() == Timer()
            True
            sage: t = Timer().start()
            sage: loads(dumps(t)) != t
            False
        """
        return not (self == other)


# Inheritance rather then delegation as globals() must be a dict
class RecordingDict(dict):
    """
    This dictionary is used for tracking the dependencies of an example.

    This feature allows examples in different doctests to be grouped
    for better timing data.  It's obtained by recording whenever
    anything is set or retrieved from this dictionary.

    EXAMPLES::

        sage: from sage.doctest.util import RecordingDict
        sage: D = RecordingDict(test=17)
        sage: D.got
        set()
        sage: D['test']
        17
        sage: D.got
        {'test'}
        sage: D.set
        set()
        sage: D['a'] = 1
        sage: D['a']
        1
        sage: D.set
        {'a'}
        sage: D.got
        {'test'}

    TESTS::

        sage: TestSuite(D).run()
    """
    def __init__(self, *args, **kwds):
        """
        Initialization arguments are the same as for a normal dictionary.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D.got
            set()
        """
        dict.__init__(self, *args, **kwds)
        self.start()

    def start(self):
        """
        We track which variables have been set or retrieved.
        This function initializes these lists to be empty.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D.set
            set()
            sage: D['a'] = 4
            sage: D.set
            {'a'}
            sage: D.start(); D.set
            set()
        """
        self.set = set([])
        self.got = set([])

    def __getitem__(self, name):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4
            sage: D.got
            set()
            sage: D['a'] # indirect doctest
            4
            sage: D.got
            set()
            sage: D['d']
            42
            sage: D.got
            {'d'}
        """
        if name not in self.set:
            self.got.add(name)
        return dict.__getitem__(self, name)

    def __setitem__(self, name, value):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4 # indirect doctest
            sage: D.set
            {'a'}
        """
        self.set.add(name)
        dict.__setitem__(self, name, value)

    def __delitem__(self, name):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: del D['d'] # indirect doctest
            sage: D.set
            {'d'}
        """
        self.set.add(name)
        dict.__delitem__(self, name)

    def get(self, name, default=None):
        """
        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D.get('d')
            42
            sage: D.got
            {'d'}
            sage: D.get('not_here')
            sage: sorted(list(D.got))
            ['d', 'not_here']
        """
        if name not in self.set:
            self.got.add(name)
        return dict.get(self, name, default)

    def copy(self):
        """
        Note that set and got are not copied.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4
            sage: D.set
            {'a'}
            sage: E = D.copy()
            sage: E.set
            set()
            sage: sorted(E.keys())
            ['a', 'd']
        """
        return RecordingDict(dict.copy(self))

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: from sage.doctest.util import RecordingDict
            sage: D = RecordingDict(d = 42)
            sage: D['a'] = 4
            sage: D.get('not_here')
            sage: E = loads(dumps(D))
            sage: E.got
            {'not_here'}
        """
        return make_recording_dict, (dict(self), self.set, self.got)


def make_recording_dict(D, st, gt):
    """
    Auxiliary function for pickling.

    EXAMPLES::

        sage: from sage.doctest.util import make_recording_dict
        sage: D = make_recording_dict({'a':4,'d':42},set([]),set(['not_here']))
        sage: sorted(D.items())
        [('a', 4), ('d', 42)]
        sage: D.got
        {'not_here'}
    """
    ans = RecordingDict(D)
    ans.set = st
    ans.got = gt
    return ans

class NestedName:
    """
    Class used to construct fully qualified names based on indentation level.

    EXAMPLES::

        sage: from sage.doctest.util import NestedName
        sage: qname = NestedName('sage.categories.algebras')
        sage: qname[0] = 'Algebras'; qname
        sage.categories.algebras.Algebras
        sage: qname[4] = '__contains__'; qname
        sage.categories.algebras.Algebras.__contains__
        sage: qname[4] = 'ParentMethods'
        sage: qname[8] = 'from_base_ring'; qname
        sage.categories.algebras.Algebras.ParentMethods.from_base_ring

    TESTS::

        sage: TestSuite(qname).run()
    """
    def __init__(self, base):
        """
        INPUT:

        - base -- a string: the name of the module.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname
            sage.categories.algebras
        """
        self.all = [base]

    def __setitem__(self, index, value):
        """
        Sets the value at a given indentation level.

        INPUT:

        - index -- a positive integer, the indentation level (often a multiple of 4, but not necessarily)
        - value -- a string, the name of the class or function at that indentation level.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname[1] = 'Algebras' # indirect doctest
            sage: qname
            sage.categories.algebras.Algebras
            sage: qname.all
            ['sage.categories.algebras', None, 'Algebras']
        """
        if index < 0:
            raise ValueError
        while len(self.all) <= index:
            self.all.append(None)
        self.all[index+1:] = [value]

    def __str__(self):
        """
        Return a .-separated string giving the full name.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname[1] = 'Algebras'
            sage: qname[44] = 'at_the_end_of_the_universe'
            sage: str(qname) # indirect doctest
            'sage.categories.algebras.Algebras.at_the_end_of_the_universe'
        """
        return repr(self)

    def __repr__(self):
        """
        Return a .-separated string giving the full name.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname[1] = 'Algebras'
            sage: qname[44] = 'at_the_end_of_the_universe'
            sage: print(qname) # indirect doctest
            sage.categories.algebras.Algebras.at_the_end_of_the_universe
        """
        return '.'.join(a for a in self.all if a is not None)

    def __eq__(self, other):
        """
        Comparison is just comparison of the underlying lists.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname2 = NestedName('sage.categories.algebras')
            sage: qname == qname2
            True
            sage: qname[0] = 'Algebras'
            sage: qname2[2] = 'Algebras'
            sage: repr(qname) == repr(qname2)
            True
            sage: qname == qname2
            False
        """
        if not isinstance(other, NestedName):
            return False
        return self.all == other.all

    def __ne__(self, other):
        """
        Test for non-equality.

        EXAMPLES::

            sage: from sage.doctest.util import NestedName
            sage: qname = NestedName('sage.categories.algebras')
            sage: qname2 = NestedName('sage.categories.algebras')
            sage: qname != qname2
            False
            sage: qname[0] = 'Algebras'
            sage: qname2[2] = 'Algebras'
            sage: repr(qname) == repr(qname2)
            True
            sage: qname != qname2
            True
        """
        return not (self == other)
