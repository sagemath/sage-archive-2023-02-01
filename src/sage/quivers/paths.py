#*****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import MonoidElement
from sage.rings.integer_ring import ZZ

class QuiverPath(MonoidElement):
    """
    Class for paths in a quiver.

    A path is given by two vertices, start and end, and a finite (possibly empty)
    list of edges e_1, e_2, ..., e_n such that the initial vertex of e_1 is start,
    the final vertex of e_i is the initial vertex of e_{i + 1}, and the final
    vertex of e_n is end.  In the case where no edges are specified we must have
    start = end and the path is called the trivial path at the given vertex.  This
    class can also represent an invalid path.

    INPUT:

    - ``path`` -- tuple or iterable (default: empty list), if ``path`` is a tuple then it is
      assumed to be of the form (vertex, vertex) or (start, end, label).  In the
      first case the trivial path at the given vertex is created.  In the second
      case a path consisting of just the given edge is created.  If ``path`` is not a
      tuple then it is assumed to be an iterable variable giving the edges of a
      path, where each edge is in one of the two forms above.  If these edges do
      not compose to form a valid path then an invalid path is returned.
    - ``check`` -- bool (default: True). If it is False, no sanity check will be performed
      on the given iterable.
    - ``parent`` -- the free small category associated with a quiver

    OUTPUT:

    - QuiverPath

    NOTE:

    Do *not* use this constructor directly! Instead, pass the input to the
    free small category that shall be the parent of this path.

    EXAMPLES::

    Specify a path by giving a list of edges::

        sage: Q = DiGraph({1:{2:['a','d'], 3:['e']}, 2:{3:['b']}, 3:{1:['f'], 4:['c']}})
        sage: F = Q.path_semigroup()
        sage: p = F([(1, 2, 'a'), (2, 2), (2, 3, 'b')])
        sage: p
        a*b

    Paths are not *unique*, but different representations of "the same" path
    yield *equal* paths::

        sage: q = F([(1, 1), (1, 2, 'a'), (2, 3, 'b'), (3, 3)])
        sage: p is q
        False
        sage: p == q
        True

    There is an "invalid path", that is obtained by concatenating arrows whose
    endpoints do not match.  It can be detected by converting a path to a
    Boolean. Valid paths are True, an invalid path is False::

        sage: inv1 = F([(1, 1), (2, 2)])
        sage: inv2 = F([(1, 2, 'a'), (1, 2, 'a')])
        sage: inv1 is inv2
        False
        sage: inv1 == inv2
        True
        sage: bool(p)
        True
        sage: bool(inv1)
        False

    The ``*`` operator is concatenation of paths. If the two paths do not
    compose the result is the invalid path::

        sage: bool(p*q)
        False
        sage: p*F((3, 4, 'c'))
        a*b*c
        sage: F([(2,3,'b'), (3,1,'f')])*p
        b*f*a*b

    The length of a path is the number of edges in that path.  The invalid path and
    trivial paths are therefore length 0::

        sage: len(p)
        2
        sage: triv = F((1, 1))
        sage: len(triv)
        0
        sage: len(inv1)
        0

    List index and slice notation can be used to access the edges in a path.
    QuiverPaths can also be iterated over.  Trivial paths and the invalid path have
    no elements::

        sage: for x in p: print x
        (1, 2, 'a')
        (2, 3, 'b')
        sage: triv[:]
        []
        sage: inv1[0]
        Traceback (most recent call last):
        ...
        IndexError: tuple index out of range

    There are methods giving the initial and terminal vertex of a path.  These
    return None when called on the invalid path::

        sage: inv1.initial_vertex()
        sage: inv1.terminal_vertex()
        sage: p.initial_vertex()
        1
        sage: p.terminal_vertex()
        3
    """
    def __init__(self, parent, path, check=True):
        """
        Creates a path object.  Type QuiverPath? for more information.

        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 1), (1, 1)])

        Note that QuiverPath should not be called directly, because
        the elements of the free small category associated with a quiver
        use a sub-class of QuiverPath. Nonetheless, just for test, we
        show that it *is* possible to create a path in a deprecated way::

            sage: p == QuiverPath((1, 1), parent=Q)
            True
            sage: Q([(1, 1), (1, 2, 'a'), (2, 2), (2, 3, 'b'), (3, 3)])._path
            ((1, 2, 'a'), (2, 3, 'b'))
        """
        MonoidElement.__init__(self, parent=parent)
        if path is None: # invalid path
            self._path = None
            return
        # Normalise the given path, unless it is asserted that the input is
        # fine
        if not check:
            self._path = tuple(path)
            return
        if path == 1:
            # We do not guarantee that there is only one vertex.
            # However, this element certainly exists.
            v = parent.quiver().vertices()[0]
            self._path = ((v,v),)
            return
        E = parent.quiver().edges()
        if isinstance(path, QuiverPath):
            if path.parent() is parent:
                self._path = path._path
                return
            new_path = list(path._path)
        # A tuple is assumed to be an edge, anything else is assumed to be a
        # list of edges
        elif path is None:
            new_path = []
        elif isinstance(path, tuple):
            new_path = [path]
        else:
            new_path = list(path)

        # Check that each edge in the path is valid
        good = True
        for x in new_path:
            if (len(x) < 2 or x[0] not in ZZ or x[1] not in ZZ
                           or len(x) == 2 and x[0] != x[1]
                           or len(x) == 3 and x not in E
                           or len(x) > 3):
                good = False
                break
        if not good:
            new_path = []
        # Delete trivial edges, and clear the path if not valid
        i = 0
        while i + 1 < len(new_path):
            if new_path[i][1] != new_path[i + 1][0]:
                new_path = []
            elif len(new_path[i])!=3:
                del new_path[i]
            else:
                i += 1
        if len(new_path) > 1 and len(new_path[-1])!=3:
            del new_path[-1]
        self._path = tuple(new_path)

    def _repr_(self):
        """
        Default representation of a path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: Q([(1, 2, 'a'), (1, 2, 'a')]) # indirect doctest
            invalid path
            sage: Q([(1, 2, 'a'), (2, 3, 'b')]) # indirect doctest
            a*b
            sage: Q((1, 1)) # indirect doctest
            e_1
        """
        if not self._path:
            return 'invalid path'
        elif len(self._path[0])!=3:
            return 'e_{0}'.format(self._path[0][0])
        else:
            return '*'.join([e[2] for e in self._path])

    def __len__(self):
        """
        Returns the length of the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: len(Q([(1, 2, 'a'), (1, 2, 'a')]))
            0
            sage: len(Q([(1, 2, 'a'), (2, 3, 'b')]))
            2
            sage: len(Q((1, 1)))
            0
            sage: len(Q((1, 2, 'a')))
            1
        """
        if (not self._path) or (len(self._path[0])==2):
            return 0
        else:
            return len(self._path)

    def deg(self):
        """
        Returns the degree of the path, which is the same as its length.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: Q([(1, 2, 'a'), (1, 2, 'a')]).deg()
            0
            sage: Q([(1, 2, 'a'), (2, 3, 'b')]).deg()
            2
            sage: Q((1, 1)).deg()
            0
            sage: Q((1, 2, 'a')).deg()
            1
        """
        return len(self)

    def __nonzero__(self):
        """
        Implements boolean values for the object.

        NOTE:

        The boolean value is True if and only if the path
        is not the invalid path

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: a = Q((1, 2, 'a'))
            sage: b = Q((2, 3, 'b'))
            sage: bool(a*b)
            True
            sage: bool(b*a)
            False
        """
        return bool(self._path)

    def __cmp__(self, other):
        """
        Comparison for QuiverPaths.

        If other is not of type QuiverPath then the comparison is made using `id`.
        If other is a QuiverPath then equality is tested using `is`.  If the
        QuiverPaths are unequal then one of the following data (listed in order of
        preferance) is unequal and used for comparison:

        - Length of the path
        - String representation of the path
        - Vertexes of the path (from initial to final)

        .. NOTES::

            This code is used by CombinatorialFreeModule to order the monomials
            when printing elements of QuiverAlgebras.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c'], 4:['d']}}).path_semigroup()
            sage: a = Q((1, 2, 'a'))
            sage: b = Q((1, 2, 'b'))
            sage: c = Q((2, 3, 'c'))
            sage: c2 = Q((2, 4, 'd'))
            sage: a < a*b
            True
            sage: a*b < a
            False
            sage: a < b
            True
            sage: b < a
            False
            sage: a*c < a*c2
            True
            sage: a*c2 < a*c
            False
            sage: a < a
            False
        """
        # Since QuiverPath inherits from element, it is guaranteed that
        # both arguments are elements of the same free small category

        # Compare lengths if different
        c = cmp(len(other),len(self))
        if c:
            return c

        # Compare internal tuple
        return cmp(self._path, other._path)

    def __getitem__(self, *args):
        """
        Implements index notation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: p
            a*b*c
            sage: p[0]
            (1, 2, 'a')
            sage: p[-1]
            (3, 4, 'c')
            sage: p[1:]
            ((2, 3, 'b'), (3, 4, 'c'))
        """

        if self._path and self._path[0][0] == self._path[0][1]:
            return list().__getitem__(*args)
        else:
            return self._path.__getitem__(*args)

    def __iter__(self):
        """
        Implements iteration over the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: for e in p: print e
            (1, 2, 'a')
            (2, 3, 'b')
            (3, 4, 'c')
        """

        # Return an iterator over an empty tuple for trivial paths, otherwise
        # return an iterator for _path as a list
        if not len(self):
            return list().__iter__()
        else:
            return list(self._path).__iter__()

    def _mul_(self, other):
        """
        Composes two paths.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}, 4:{5:['d']}}).path_semigroup()
            sage: x = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y = Q([(3, 4, 'c'), (4, 5, 'd')])
            sage: y*x
            invalid path
            sage: x*y
            a*b*c*d
            sage: x*Q((3, 4, 'c'))
            a*b*c
            sage: x*Q([(3, 4, 'c'), (4, 5, 'd')])
            a*b*c*d
            sage: x*6
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Free small category of Quiver on 5 vertices' and 'Integer Ring'
        """
        # By Sage's coercion model, both paths belong to the same quiver
        Q = self.parent()
        # Handle the invalid path
        if not self._path:
            return self

        if isinstance(other, QuiverPath):
            if not other._path:
                return other
            return Q(list(self._path) + list(other._path))
        if isinstance(other, tuple):
            return Q(list(self._path) + [other])
        if isinstance(other, list):
            return Q(list(self._path) + other)
        raise TypeError("QuiverPath cannot be multiplied with {0}".format(type(other)))

    def __mod__(self, other):
        """
        Returns self with other deleted from the beginning.

        If other is not the beginning of self the result is the invalid path.  Deleting
        the trivial path at vertex v from a path that begins at v does nothing.
        Deleting it from a path that does not begin at v returns the invalid path.
        Deleting the invalid path returns the invalid path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: a = Q((1, 2, 'a'))
            sage: b = Q((2, 3, 'b'))
            sage: e1 = Q((1, 1))
            sage: e2 = Q((2, 2))
            sage: p % a
            b
            sage: p % b
            invalid path
            sage: p % e1
            a*b
            sage: p % e2
            invalid path

        """
        # Handle invalid path
        if not self._path:
            return self
        Q = self.parent()
        # Convert other to a QuiverPath
        oth = Q(other)
        # Is oth invalid or equal to self?
        if not oth._path:
            return oth
        if self._path==oth._path:
            v = self._path[-1][1]
            return Q(((v, v),), check=False)

        # Handle trivial paths
        if oth._path[0][0] == oth._path[0][1]:
            if self._path[0][0] == oth._path[0][0]:
                return self
            else:
                return Q(None)

        # If other is the beginning, return the rest
        if self._path[:len(oth._path)] == oth._path:
            return Q(list(self._path[len(oth._path):]), check=False)
        else:
            return Q(None)

    def initial_vertex(self):
        """
        Returns the initial vertex of the path.

        The invalid path does not have an initial vertex, so None is returned.

        OUTPUT:

        - integer or None

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.initial_vertex()
            1
            sage: print Q([(1, 1), (2, 2)]).initial_vertex()
            None
        """

        if self._path:
            return self._path[0][0]
        else:
            return None

    def terminal_vertex(self):
        """
        Returns the terminal vertex of the path.

        The invalid path does not have a terminal vertex, so None is returned.

        OUTPUT:

        - integer or None

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.terminal_vertex()
            3
            sage: print Q([(1, 1), (2, 2)]).terminal_vertex()
            None
        """

        if self._path:
            return self._path[-1][1]
        else:
            return None

    def reverse(self):
        """
        Returns the path along the same edges in the opposite quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: p.reverse()
            b*a
            sage: p.reverse().parent() is Q.reverse()
            True

        """
        Q = self.parent().reverse()
        # Handle the invalid path and trivial paths
        if not self._path:
            return Q(None)
        if self._path[0][0] == self._path[0][1]:
            return Q((self._path[0][0],self._path[0][0]))

        # Reverse all the edges in the path, then reverse the path
        new_path = [(e[1], e[0], e[2]) for e in self._path]
        return Q(reversed(new_path))
