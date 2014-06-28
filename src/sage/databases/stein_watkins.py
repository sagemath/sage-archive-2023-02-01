r"""
The Stein-Watkins table of elliptic curves

Sage gives access to the Stein-Watkins table of elliptic curves, via an
optional package that you must install. This is a huge database of elliptic
curves. You can install the database (a 2.6GB package) with the command

::

    sage -i database_stein_watkins

You can also automatically download a small version, which takes much less
time, using the command

::

    sage -i database_stein_watkins_mini

This database covers a wide range of conductors, but unlike the
:mod:`Cremona database <sage.databases.cremona>`, this database need not list
all curves of a given conductor. It lists the curves whose coefficients are not
"too large" (see [SteinWatkins]_).


-  The command ``SteinWatkinsAllData(n)`` returns an iterator over the curves
   in the `n`-th Stein-Watkins table, which contains elliptic curves of
   conductor between `n10^5` and `(n+1)10^5`. Here `n` can be between 0 and
   999, inclusive.

-  The command ``SteinWatkinsPrimeData(n)`` returns an iterator over the curves
   in the `n^{th}` Stein-Watkins prime table, which contains prime conductor
   elliptic curves of conductor between `n10^6` and `(n+1)10^6`. Here `n`
   varies between 0 and 99, inclusive.


EXAMPLES: We obtain the first table of elliptic curves.

::

    sage: d = SteinWatkinsAllData(0)
    sage: d
    Stein-Watkins Database a.0 Iterator

We type ``d.next()`` to get each isogeny class of
curves from ``d``::

    sage: C = d.next()                                   # optional - database_stein_watkins
    sage: C                                              # optional - database_stein_watkins
    Stein-Watkins isogeny class of conductor 11
    sage: d.next()                                       # optional - database_stein_watkins
    Stein-Watkins isogeny class of conductor 14
    sage: d.next()                                       # optional - database_stein_watkins
    Stein-Watkins isogeny class of conductor 15

An isogeny class has a number of attributes that give data about
the isogeny class, such as the rank, equations of curves,
conductor, leading coefficient of `L`-function, etc.

::

    sage: C.data                                         # optional - database_stein_watkins
    ['11', '[11]', '0', '0.253842', '25', '+*1']
    sage: C.curves                                       # optional - database_stein_watkins
    [[[0, -1, 1, 0, 0], '(1)', '1', '5'],
     [[0, -1, 1, -10, -20], '(5)', '1', '5'],
     [[0, -1, 1, -7820, -263580], '(1)', '1', '1']]
    sage: C.conductor                                    # optional - database_stein_watkins
    11
    sage: C.leading_coefficient                          # optional - database_stein_watkins
    '0.253842'
    sage: C.modular_degree                               # optional - database_stein_watkins
    '+*1'
    sage: C.rank                                         # optional - database_stein_watkins
    0
    sage: C.isogeny_number                               # optional - database_stein_watkins
    '25'

If we were to continue typing ``d.next()`` we would
iterate over all curves in the Stein-Watkins database up to
conductor `10^5`. We could also type ``for C in d:
...``

To access the data file starting at `10^5` do the
following::

    sage: d = SteinWatkinsAllData(1)
    sage: C = d.next()                                  # optional - database_stein_watkins
    sage: C                                             # optional - database_stein_watkins
    Stein-Watkins isogeny class of conductor 100002
    sage: C.curves                                      # optional - database_stein_watkins
    [[[1, 1, 0, 112, 0], '(8,1,2,1)', 'X', '2'],
     [[1, 1, 0, -448, -560], '[4,2,1,2]', 'X', '2']]

Next we access the prime-conductor data::

    sage: d = SteinWatkinsPrimeData(0)
    sage: C = d.next()                                 # optional - database_stein_watkins
    sage: C                                            # optional - database_stein_watkins
    Stein-Watkins isogeny class of conductor 11

Each call ``d.next()`` gives another elliptic curve of
prime conductor::

    sage: C = d.next()                                 # optional - database_stein_watkins
    sage: C                                            # optional - database_stein_watkins
    Stein-Watkins isogeny class of conductor 17
    sage: C.curves                                     # optional - database_stein_watkins
    [[[1, -1, 1, -1, 0], '[1]', '1', '4'],
     [[1, -1, 1, -6, -4], '[2]', '1', '2x'],
     [[1, -1, 1, -1, -14], '(4)', '1', '4'],
     [[1, -1, 1, -91, -310], '[1]', '1', '2']]
    sage: C = d.next()                                 # optional - database_stein_watkins
    sage: C                                            # optional - database_stein_watkins
    Stein-Watkins isogeny class of conductor 19

REFERENCE:

.. [SteinWatkins] William Stein and Mark Watkins, *A database of elliptic
   curves---first report*. In *Algorithmic number theory (ANTS V), Sydney,
   2002*, Lecture Notes in Computer Science 2369, Springer, 2002, p267--275.
   http://modular.math.washington.edu/papers/stein-watkins/

"""

#*****************************************************************************
#
#       Sage: Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import bz2, os

from sage.misc.misc import SAGE_SHARE

class SteinWatkinsIsogenyClass:
    def __init__(self, conductor):
        self.conductor = conductor

    def __repr__(self):
        return "Stein-Watkins isogeny class of conductor %s"%self.conductor

    def __len__(self):
        try:
            return len(self.curves)
        except AttributeError:
            return 0

    def __iter__(self):
        try:
            for E in self.curves:
                yield E
        except AttributeError:
            return


def _lines(s):
    while True:
        i = s.find("\n")
        if i == -1:
            yield ""
            return
        line = s[:i]
        s = s[i+1:]
        yield line


class SteinWatkinsAllData:
    """
    Class for iterating through one of the Stein-Watkins database files
    for all conductors.
    """
    def __init__(self, num):
        num = int(num)
        self.num = num
        if num < 0:
            raise RuntimeError("num (=%s) must be a nonnegative integer"%num)
        name = str(num)
        name = '0'*(3-len(name)) + name
        self._file = os.path.join(SAGE_SHARE, 'stein_watkins', 'a.%s.bz2'%name)
        self._iter = self.__iter__()

    def __repr__(self):
        """
        EXAMPLES::

            sage: d = SteinWatkinsAllData(1)
            sage: d
            Stein-Watkins Database a.1 Iterator
        """
        return "Stein-Watkins Database a.%s Iterator"%self.num

    def __iter__(self):
        """
        EXAMPLES::

            sage: d = SteinWatkinsAllData(0)
            sage: d = d[10:20]                         # optional - database_stein_watkins; long time
            sage: for C in d:                          # optional - database_stein_watkins; long time
            ....:     print C
            Stein-Watkins isogeny class of conductor 11
            Stein-Watkins isogeny class of conductor 14
            Stein-Watkins isogeny class of conductor 15
            Stein-Watkins isogeny class of conductor 17
            Stein-Watkins isogeny class of conductor 19
            Stein-Watkins isogeny class of conductor 20
        """
        try:
            file = bz2.BZ2File(self._file, 'r')
        except IOError:
            raise IOError("The Stein-Watkins data file %s must be installed."%self._file)
        C = None
        for L in file:
            if len(L) == 0:
                continue
            if L[0] != '[': # new curve
                if C is not None:
                    yield C
                x = L.split()
                N = int(x[0])
                C = SteinWatkinsIsogenyClass(N)
                C.rank = int(x[2])
                C.leading_coefficient = x[3]
                C.isogeny_number = x[4]
                C.modular_degree = x[5]
                C.curves = []
                C.data = x
            else:
                w = L.split()
                C.curves.append([eval(w[0]), w[1], w[2], w[3]])
        yield C

    def next(self):
        return self._iter.next()

    def __getitem__(self, N):
        """
        Return the curves of conductor N in this table. (Very slow!)
        Return all data about curves between the given levels in this
        database file.

        EXAMPLES::

            sage: d = SteinWatkinsAllData(0)
            sage: d[15:18]                             # optional - database_stein_watkins; long time
            [Stein-Watkins isogeny class of conductor 15, Stein-Watkins isogeny
             class of conductor 17]
        """
        X = []
        if isinstance(N, slice):
            min_level, max_level, step = N.indices(len(list(self)))
            for C in self:
                M = C.conductor
                if M >= min_level and M <= max_level:
                    X.append(C)
                elif M > max_level:
                    return X
        else:
            for C in self:
                M = C.conductor
                if M == N:
                    X.append(C)
                elif M > N:
                    return X
            return X

    def iter_levels(self):
        """
        Iterate through the curve classes, but grouped into lists by
        level.

        EXAMPLE::

            sage: d = SteinWatkinsAllData(1)
            sage: E = d.iter_levels()
            sage: E.next()                             # optional - database_stein_watkins
            [Stein-Watkins isogeny class of conductor 100002]
            sage: E.next()                             # optional - database_stein_watkins
            [Stein-Watkins isogeny class of conductor 100005,
            Stein-Watkins isogeny class of conductor 100005]
            sage: E.next()                             # optional - database_stein_watkins
            [Stein-Watkins isogeny class of conductor 100007]
        """
        iter = self.__iter__()
        C = []
        N = 0
        while True:
            try:
                E = iter.next()
            except StopIteration:
                if C != []:
                    yield C
                raise StopIteration
            if E.conductor != N:
                if C != []:
                    yield C
                C = [E]
                N = E.conductor
            else:
                C.append(E)
        yield C


class SteinWatkinsPrimeData(SteinWatkinsAllData):
    def __init__(self, num):
        num = int(num)
        self.num = num
        if num < 0:
            raise RuntimeError("num (=%s) must be a nonnegative integer"%num)
        name = str(num)
        name = '0'*(2-len(name)) + name
        self._file = os.path.join(SAGE_SHARE,'stein_watkins', 'p.%s.bz2'%name)
        self._iter = self.__iter__()

    def __repr__(self):
        """
        EXAMPLES::

            sage: d = SteinWatkinsPrimeData(1)
            sage: d
            Stein-Watkins Prime Conductor Database p.1 Iterator
        """
        return "Stein-Watkins Prime Conductor Database p.%s Iterator"%self.num


def ecdb_num_curves(max_level=200000):
    r"""
    Return a list whose `N`-th entry, for ``0 <= N <= max_level``, is the
    number of elliptic curves of conductor `N` in the database.

    EXAMPLES::

        sage: sage.databases.stein_watkins.ecdb_num_curves(100) # optional - database_stein_watkins
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 6, 8, 0, 4, 0, 3, 4, 6, 0, 0,
         6, 0, 5, 4, 0, 0, 8, 0, 4, 4, 4, 3, 4, 4, 5, 4, 4, 0, 6, 1, 2, 8, 2, 0,
         6, 4, 8, 2, 2, 1, 6, 4, 6, 7, 3, 0, 0, 1, 4, 6, 4, 2, 12, 1, 0, 2, 4, 0,
         6, 2, 0, 12, 1, 6, 4, 1, 8, 0, 2, 1, 6, 2, 0, 0, 1, 3, 16, 4, 3, 0, 2,
         0, 8, 0, 6, 11, 4]
    """
    i = 0
    N = 1
    d = SteinWatkinsAllData(i)
    v = [int(0) for _ in xrange(max_level+1)]
    while True:
        try:
            C = d.next()
        except StopIteration:
            i += 1
            d = SteinWatkinsAllData(i)
            continue
        N = C.conductor
        if N > max_level:
            break
        v[N] += len(C.curves)
    return v


