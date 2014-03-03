r"""
Some useful functions for the matroid class.

For direct access to the methods :meth:`newlabel`, :meth:`setprint` and
:meth:`get_nonisomorphic_matroids`, type::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.

AUTHORS:

- Stefan van Zwam (2011-06-24): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import Matrix
from sage.rings.all import ZZ, QQ, FiniteField, GF
from sage.graphs.all import BipartiteGraph
from pprint import pformat
from sage.structure.all import SageObject


def setprint(X):
    """
    Print nested data structures nicely.

    Python's data structures ``set`` and ``frozenset`` do not print nicely.
    This function can be used as replacement for ``print`` to overcome this.
    For direct access to ``setprint``, run::

        sage: from sage.matroids.advanced import *

    .. NOTE::

        This function will be redundant when Sage moves to Python 3, since the
        default ``print`` will suffice then.

    INPUT:

    - ``X`` -- Any Python object

    OUTPUT:

    ``None``. However, the function prints a nice representation of ``X``.

    EXAMPLES:

    Output looks much better::

        sage: from sage.matroids.advanced import setprint
        sage: L = [{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {4, 1, 3}]
        sage: print(L)
        [set([1, 2, 3]), set([1, 2, 4]), set([2, 3, 4]), set([1, 3, 4])]
        sage: setprint(L)
        [{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {1, 3, 4}]

    Note that for iterables, the effect can be undesirable::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.Fano().delete('efg')
        sage: M.bases()
        Iterator over a system of subsets
        sage: setprint(M.bases())
        [{'a', 'b', 'c'}, {'a', 'c', 'd'}, {'a', 'b', 'd'}]

    An exception was made for subclasses of SageObject::

        sage: from sage.matroids.advanced import setprint
        sage: G = graphs.PetersenGraph()
        sage: list(G)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: setprint(G)
        Petersen graph: Graph on 10 vertices
    """
    print setprint_s(X, toplevel=True)


def setprint_s(X, toplevel=False):
    """
    Create the string for use by ``setprint()``.

    INPUT:

    - ``X`` -- any Python object
    - ``toplevel`` -- (default: ``False``) indicates whether this is a
      recursion or not.

    OUTPUT:

    A string representation of the object, with nice notation for sets and
    frozensets.

    EXAMPLES::

        sage: from sage.matroids.utilities import setprint_s
        sage: L = [{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {4, 1, 3}]
        sage: setprint_s(L)
        '[{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {1, 3, 4}]'

    The ``toplevel`` argument only affects strings, to mimic ``print``'s
    behavior::

        sage: X = 'abcd'
        sage: setprint_s(X)
        "'abcd'"
        sage: setprint_s(X, toplevel=True)
        'abcd'
    """
    if isinstance(X, frozenset) or isinstance(X, set):
        return '{' + ', '.join([setprint_s(x) for x in sorted(X)]) + '}'
    elif isinstance(X, dict):
        return '{' + ', '.join([setprint_s(key) + ': ' + setprint_s(val) for key, val in sorted(X.iteritems())]) + '}'
    elif isinstance(X, str):
        if toplevel:
            return X
        else:
            return "'" + X + "'"
    elif hasattr(X, '__iter__') and not isinstance(X, SageObject):
        return '[' + ', '.join([setprint_s(x) for x in sorted(X)]) + ']'
    else:
        return repr(X)


def newlabel(groundset):
    r"""
    Create a new element label different from the labels in ``groundset``.

    INPUT:

    - ``groundset`` -- A set of objects.

    OUTPUT:

    A string not in the set ``groundset``.

    For direct access to ``newlabel``, run::

        sage: from sage.matroids.advanced import *

    ALGORITHM:

    #. Create a set of all one-character alphanumeric strings.
    #. Remove the string representation of each existing element from this
       list.
    #. If the list is nonempty, return any element.
    #. Otherwise, return the concatenation of the strings of each existing
       element, preceded by 'e'.

    EXAMPLES::

        sage: from sage.matroids.advanced import newlabel
        sage: S = set(['a', 42, 'b'])
        sage: newlabel(S) in S
        False

        sage: S = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
        sage: t = newlabel(S)
        sage: len(t)
        63
        sage: t[0]
        'e'

    """
    char_list = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
    char_list.difference_update([str(e) for e in groundset])
    try:
        s = char_list.pop()
    except KeyError:
        s = 'e' + ''.join([str(e) for e in groundset])
    return s


def sanitize_contractions_deletions(matroid, contractions, deletions):
    r"""
    Return a fixed version of sets ``contractions`` and ``deletions``.

    INPUT:

    - ``matroid`` -- a :class:`Matroid <sage.matroids.matroid.Matroid>`
      instance.
    - ``contractions`` -- a subset of the groundset.
    - ``deletions`` -- a subset of the groundset.

    OUTPUT:

    An independent set ``C`` and a coindependent set ``D`` such that

        ``matroid / contractions \ deletions == matroid / C \ D``

    Raise an error if either is not a subset of the groundset of ``matroid``
    or if they are not disjoint.

    This function is used by the
    :meth:`Matroid.minor() <sage.matroids.matroid.Matroid.minor>` method.

    EXAMPLES::

        sage: from sage.matroids.utilities import setprint
        sage: from sage.matroids.utilities import sanitize_contractions_deletions
        sage: M = matroids.named_matroids.Fano()
        sage: setprint(sanitize_contractions_deletions(M, 'abc', 'defg'))
        [{'a', 'b', 'c'}, {'d', 'e', 'f', 'g'}]
        sage: setprint(sanitize_contractions_deletions(M, 'defg', 'abc'))
        [{'d', 'e', 'g'}, {'a', 'b', 'c', 'f'}]
        sage: setprint(sanitize_contractions_deletions(M, [1, 2, 3], 'efg'))
        Traceback (most recent call last):
        ...
        ValueError: input contractions is not a subset of the groundset.
        sage: setprint(sanitize_contractions_deletions(M, 'efg', [1, 2, 3]))
        Traceback (most recent call last):
        ...
        ValueError: input deletions is not a subset of the groundset.
        sage: setprint(sanitize_contractions_deletions(M, 'ade', 'efg'))
        Traceback (most recent call last):
        ...
        ValueError: contraction and deletion sets are not disjoint.

    """
    if contractions is None:
        contractions = frozenset([])
    contractions = frozenset(contractions)
    if not matroid.groundset().issuperset(contractions):
        raise ValueError("input contractions is not a subset of the groundset.")

    if deletions is None:
        deletions = frozenset([])
    deletions = frozenset(deletions)
    if not matroid.groundset().issuperset(deletions):
        raise ValueError("input deletions is not a subset of the groundset.")

    if not contractions.isdisjoint(deletions):
        raise ValueError("contraction and deletion sets are not disjoint.")

    conset = matroid._max_independent(contractions)
    delset = matroid._max_coindependent(deletions)

    return conset.union(deletions.difference(delset)), delset.union(contractions.difference(conset))


def make_regular_matroid_from_matroid(matroid):
    r"""
    Attempt to construct a regular representation of a matroid.

    INPUT:

    - ``matroid`` -- a matroid.

    OUTPUT:

    Return a `(0, 1, -1)`-matrix over the integers such that, if the input is
    a regular matroid, then the output is a totally unimodular matrix
    representing that matroid.

    EXAMPLES::

        sage: from sage.matroids.utilities import make_regular_matroid_from_matroid
        sage: make_regular_matroid_from_matroid(
        ....:               matroids.CompleteGraphic(6)).is_isomorphic(
        ....:                                     matroids.CompleteGraphic(6))
        True
    """
    import sage.matroids.linear_matroid
    M = matroid
    if isinstance(M, sage.matroids.linear_matroid.RegularMatroid):
        return M
    rk = M.full_rank()
    # First create a reduced 0-1 matrix
    B = list(M.basis())
    NB = list(M.groundset().difference(B))
    dB = {}
    i = 0
    for e in B:
        dB[e] = i
        i += 1
    dNB = {}
    i = 0
    for e in NB:
        dNB[e] = i
        i += 1
    A = Matrix(ZZ, len(B), len(NB), 0)
    G = BipartiteGraph(A.transpose())  # Sage's BipartiteGraph uses the column set as first color class. This is an edgeless graph.
    for e in NB:
        C = M.circuit(B + [e])
        for f in C.difference([e]):
            A[dB[f], dNB[e]] = 1
    # Change some entries from -1 to 1
    entries = BipartiteGraph(A.transpose()).edges(labels=False)
    while len(entries) > 0:
        L = [G.shortest_path(u, v) for u, v in entries]
        mindex, minval = min(enumerate(L), key=lambda x: len(x[1]))

        # if minval = 0, there is an edge not spanned by the current subgraph. Its entry is free to be scaled any way.
        if len(minval) > 0:
            # Check the subdeterminant
            S = frozenset(L[mindex])
            rows = []
            cols = []
            for i in S:
                if i < rk:
                    rows.append(i)
                else:
                    cols.append(i - rk)
            if A[rows, cols].det() != 0:
                A[entries[mindex][0], entries[mindex][1] - rk] = -1
        G.add_edge(entries[mindex][0], entries[mindex][1])
        entries.pop(mindex)
    return sage.matroids.linear_matroid.RegularMatroid(groundset=B + NB, reduced_matrix=A)


def get_nonisomorphic_matroids(MSet):
    """
    Return non-isomorphic members of the matroids in set ``MSet``.

    For direct access to ``get_nonisomorphic_matroids``, run::

        sage: from sage.matroids.advanced import *

    INPUT:

    - ``MSet`` -- an iterable whose members are matroids.

    OUTPUT:

    A list containing one representative of each isomorphism class of
    members of ``MSet``.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: L = matroids.Uniform(3, 5).extensions()
        sage: len(list(L))
        32
        sage: len(get_nonisomorphic_matroids(L))
        5
    """
    OutSet = []
    for M in MSet:
        seen = False
        for N in OutSet:
            if N.is_isomorphic(M):
                seen = True
                break
        if not seen:
            OutSet.append(M)
    return OutSet
