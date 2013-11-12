r"""
Six Vertex Model
"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

class SixVertexConfiguration(ClonableArray):
    """
    A configuration in the six vertex model.
    """
    def check(self):
        """
        Check if ``self`` is a valid 6 vertex configuration.

        EXAMPLES::

            sage: M = SixVertexModel(3, boundary_conditions='ice')
            sage: M[0].check()
        """
        if self not in self.parent():
            raise ValueError("invalid configuration")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = SixVertexModel(3, boundary_conditions='ice')
            sage: M[0]
                ^    ^    ^
                |    |    |
            --> # <- # <- # <--
                |    ^    ^
                V    |    |  
            --> # -> # <- # <--
                |    |    ^
                V    V    |
            --> # -> # -> # <--
                |    |    |
                V    V    V
        """
        # List are in the order of URDL
        ascii = [[r'  V  ', ' -', r'  ^  ', '- '], # LR
                 [r'  |  ', ' <', r'  ^  ', '- '], # LU
                 [r'  V  ', ' <', r'  |  ', '- '], # LD
                 [r'  |  ', ' <', r'  |  ', '> '], # UD
                 [r'  |  ', ' -', r'  ^  ', '> '], # UR
                 [r'  V  ', ' -', r'  |  ', '> ']] # RD
        ret = '  '
        # Do the top line
        for entry in self[0]:
            if entry == 1 or entry == 3 or entry == 4:
                ret += '  ^  '
            else:
                ret += '  |  '

        # Do the meat of the ascii art
        for row in self:
            ret += '\n  '
            # Do the top row
            for entry in row:
                ret += ascii[entry][0]
            ret += '\n'

            # Do the left-most entry
            if row[0] == 0 or row[0] == 1 or row[0] == 2:
                ret += '<-'
            else:
                ret += '--'

            # Do the middle row
            for entry in row:
                ret += ascii[entry][3] + '#' + ascii[entry][1]

            # Do the right-most entry
            if row[-1] == 0 or row[-1] == 4 or row[-1] == 5:
                ret += '->'
            else:
                ret += '--'

            # Do the bottom row
            ret += '\n  '
            for entry in row:
                ret += ascii[entry][2]

        # Do the bottom line
        ret += '\n  '
        for entry in self[-1]:
            if entry == 2 or entry == 3 or entry == 5:
                ret += '  V  '
            else:
                ret += '  |  '

        return ret

    def to_signed_matrix(self):
        """
        Return the signed matrix of ``self``.

        The signed matrix corresponding to a six vertex configuration is
        given by `0` if there is a cross flow, a `1` if the outward arrows
        are vertical and `-1` if the outward arrows are horizonal.

        EXAMPLES::

            sage: M = SixVertexModel(3, boundary_conditions='ice')
            sage: map(lambda x: x.to_signed_matrix(), M)
            [
            [1 0 0]  [1 0 0]  [ 0  1  0]  [0 1 0]  [0 1 0]  [0 0 1]  [0 0 1]
            [0 1 0]  [0 0 1]  [ 1 -1  1]  [1 0 0]  [0 0 1]  [1 0 0]  [0 1 0]
            [0 0 1], [0 1 0], [ 0  1  0], [0 0 1], [1 0 0], [0 1 0], [1 0 0]
            ]
        """
        from sage.matrix.constructor import matrix
        # verts = ['LR', 'LU', 'LD', 'UD', 'UR', 'RD']
        def matrix_sign(x):
            if x == 0:
                return -1
            if x == 3:
                return 1
            return 0
        return matrix([map(matrix_sign, row) for row in self])

    def plot(self, color='sign'):
        """
        Return a plot of ``self``.

        INPUT:

        - ``color`` -- can be any of the following:

          * ``4`` - use 4 colors: black, red, blue, and green with each
            corresponding to up, right, down, and left respectively
          * ``2`` - use 2 colors: red for horizontal, blue for vertical arrows
          * ``'sign'`` - use red for right and down arrows, blue for left
            and up arrows
          * a list of 4 colors for each direction
          * a function which takes a direction and a boolean corresponding
            to the sign

        EXAMPLES::

            sage: M = SixVertexModel(2, boundary_conditions='ice')
            sage: print M[0].plot().description()
            Arrow from (-1.0,0.0) to (0.0,0.0)
            Arrow from (-1.0,1.0) to (0.0,1.0)
            Arrow from (0.0,0.0) to (0.0,-1.0)
            Arrow from (0.0,0.0) to (0.0,1.0)
            Arrow from (0.0,1.0) to (1.0,1.0)
            Arrow from (0.0,2.0) to (0.0,1.0)
            Arrow from (1.0,0.0) to (0.0,0.0)
            Arrow from (1.0,0.0) to (1.0,-1.0)
            Arrow from (1.0,1.0) to (1.0,0.0)
            Arrow from (1.0,1.0) to (1.0,2.0)
            Arrow from (2.0,0.0) to (1.0,0.0)
            Arrow from (2.0,1.0) to (1.0,1.0)
        """
        from sage.plot.graphics import Graphics
        from sage.plot.circle import circle
        from sage.plot.arrow import arrow

        if color == 4:
            color_list = ['black', 'red', 'blue', 'green']
            cfunc = lambda d,pm: color_list[d]
        elif color == 2:
            cfunc = lambda d,pm: 'red' if d % 2 == 0 else 'blue'
        elif color == 1 or color is None:
            cfunc = lambda d,pm: 'black'
        elif color == 'sign':
            cfunc = lambda d,pm: 'red' if pm else 'blue' # RD are True
        elif isinstance(color, (list, tuple)):
            cfunc = lambda d,pm: color[d]
        else:
            cfunc = color

        G = Graphics()
        for j,row in enumerate(self):
            for i,entry in enumerate(row):
                if entry == 0: # LR
                    G += arrow((i,j+1), (i,j), color=cfunc(2, True))
                    G += arrow((i,j), (i+1,j), color=cfunc(1, True))
                    if j == 0:
                        G += arrow((i,j-1), (i,j), color=cfunc(0, False))
                    if i == 0:
                        G += arrow((i,j), (i-1,j), color=cfunc(3, False))
                elif entry == 1: # LU
                    G += arrow((i,j+1), (i,j), color=cfunc(2, True))
                    G += arrow((i+1,j), (i,j), color=cfunc(3, False))
                    if j == 0:
                        G += arrow((i,j), (i,j-1), color=cfunc(2, True))
                    if i == 0:
                        G += arrow((i,j), (i-1,j), color=cfunc(3, False))
                elif entry == 2: # LD
                    G += arrow((i,j), (i,j+1), color=cfunc(0, False))
                    G += arrow((i+1,j), (i,j), color=cfunc(3, False))
                    if j == 0:
                        G += arrow((i,j-1), (i,j), color=cfunc(0, False))
                    if i == 0:
                        G += arrow((i,j), (i-1,j), color=cfunc(3, False))
                elif entry == 3: # UD
                    G += arrow((i,j), (i,j+1), color=cfunc(0, False))
                    G += arrow((i+1,j), (i,j), color=cfunc(3, False))
                    if j == 0:
                        G += arrow((i,j), (i,j-1), color=cfunc(2, True))
                    if i == 0:
                        G += arrow((i-1,j), (i,j), color=cfunc(1, True))
                elif entry == 4: # UR
                    G += arrow((i,j), (i,j+1), color=cfunc(0, False))
                    G += arrow((i,j), (i+1,j), color=cfunc(1, True))
                    if j == 0:
                        G += arrow((i,j-1), (i,j), color=cfunc(0, False))
                    if i == 0:
                        G += arrow((i-1,j), (i,j), color=cfunc(1, True))
                elif entry == 5: # RD
                    G += arrow((i,j+1), (i,j), color=cfunc(2, True))
                    G += arrow((i,j), (i+1,j), color=cfunc(1, True))
                    if j == 0:
                        G += arrow((i,j), (i,j-1), color=cfunc(2, True))
                    if i == 0:
                        G += arrow((i-1,j), (i,j), color=cfunc(1, True))
        return G

class SixVertexModel(Parent, UniqueRepresentation):
    """
    The six vertex model.

    We model a configuration by indicating which configuration by the
    following six configurations which are determined by the two outgoing
    arrows in the **U**p, **R**ight, **D**own, **L**eft directions:

    1 - LR
    2 - LU
    3 - LD
    4 - UD
    5 - UR
    6 - RD

    INPUT:

    - ``n`` -- the number of rows
    - ``m`` -- (optional) the number of columns, if not specified, then
      the number of columns is the number of rows
    - ``boundary_conditions`` -- (optional) a quadruple of tuples whose
      entries are either:

      * ``True`` for an inward arrow,
      * ``False`` for an outward arrow, or
      * ``None`` for no boundary condition

      There are also the following special cases:

      * ``'inward'`` - all arrows are inward arrows
      * ``'outward'`` - all arrows are outward arrows
      * ``'ice'`` - the top and bottom boundary conditions are outward and the
        left and right boundary conditions are inward

    EXAMPLES:

    Here are the six types of vertices that can be created::

        sage: M = SixVertexModel(1, 1)
        sage: list(M)
        [
            |          ^          |          ^          ^          |
            V          |          V          |          |          V
        <-- # -->  <-- # <--  <-- # <--  --> # <--  --> # -->  --> # -->
            ^          ^          |          |          ^          |
            |    ,     |    ,     V    ,     V    ,     |    ,     V
        ]

    When using the square ice model, it is known that the number of
    configurations is equal to the number of alternating sign matrices::

        sage: M = SixVertexModel(1, boundary_conditions='ice')
        sage: len(M)
        1
        sage: M = SixVertexModel(4, boundary_conditions='ice')
        sage: len(M)
        42
        sage: all(len(SixVertexModel(n, boundary_conditions='ice'))
        ....:     == AlternatingSignMatrices(n).cardinality() for n in range(1, 7))
        True

    An example with a specified non-standard boundary condition and
    non-rectangular shape::

        sage: M = SixVertexModel(2, 1, [[None],[True,True],[None],[None,None]])
        sage: list(M)
        [
            ^          ^          |          ^
            |          |          V          |
        <-- # <--  <-- # <--  <-- # <--  --> # <--
            ^          ^          |          |
            |          |          V          V
        <-- # <--  --> # <--  <-- # <--  <-- # <--
            ^          |          |          |
            |    ,     V    ,     V    ,     V
        ]
    """
    @staticmethod
    def __classcall_private__(cls, n, m=None, boundary_conditions=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: M1 = SixVertexModel(1, boundary_conditions=[[False],[True],[False],[True]])
            sage: M2 = SixVertexModel(1, 1, ((False,),(True,),(False,),(True,)))
            sage: M1 is M2
            True
        """
        if m is None:
            m = n
        if boundary_conditions is None:
            boundary_conditions = ((None,)*m, (None,)*n)*2
        elif boundary_conditions == 'inward':
            boundary_conditions = ((True,)*m, (True,)*n)*2
        elif boundary_conditions == 'outward':
            boundary_conditions = ((False,)*m, (False,)*n)*2
        elif boundary_conditions == 'ice':
            boundary_conditions = ((False,)*m, (True,)*n)*2
        else:
            boundary_conditions = tuple(tuple(x) for x in boundary_conditions)
        return super(SixVertexModel, cls).__classcall__(cls, n, m, boundary_conditions)

    def __init__(self, n, m, boundary_conditions):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = SixVertexModel(2, boundary_conditions='ice')
            sage: TestSuite(M).run()
        """
        self._nrows = n
        self._ncols = m
        self._bdry_cond = boundary_conditions # Ordered URDL
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SixVertexModel(2, boundary_conditions='ice')
            The six vertex model on a 2 by 2 grid
        """
        return "The six vertex model on a {} by {} grid".format(self._nrows, self._ncols)

    def boundary_conditions(self):
        """
        Return the boundary conditions of ``self``.

        EXAMPLES::

            sage: M = SixVertexModel(2, boundary_conditions='ice')
            sage: M.boundary_conditions()
            ((False, False), (True, True), (False, False), (True, True))
        """
        return self._bdry_cond

    def _repr_option(self, key):
        """
        Metadata about the ``_repr_()`` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: M = SixVertexModel(2, boundary_conditions='ice')
            sage: M._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return True
        return Parent._repr_option(self, key)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: M = SixVertexModel(2, boundary_conditions='ice')
            sage: M([[3,1],[5,3]])
                ^    ^
                |    |
            --> # <- # <--
                |    ^
                V    |
            --> # -> # <--
                |    |
                V    V
        """
        if isinstance(x, SixVertexConfiguration):
            if x.parent() is not self:
                return self.element_class(self, tuple(x))
            return x

        verts = ['LR', 'LU', 'LD', 'UD', 'UR', 'RD']
        elt = []
        for row in x:
            elt.append([])
            for entry in row:
                if entry in verts:
                    elt[-1].append(verts.index(entry))
                elif entry in range(6):
                    elt[-1].append(entry)
                else:
                    raise ValueError("invalid entry")
            elt[-1] = tuple(elt[-1])
        return self.element_class(self, tuple(elt))

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: M = SixVertexModel(2, boundary_conditions='ice')
            sage: list(M)
            [
                ^    ^          ^    ^
                |    |          |    |
            --> # <- # <--  --> # -> # <--
                |    ^          ^    |
                V    |          |    V
            --> # -> # <--  --> # <- # <--
                |    |          |    |
                V    V    ,     V    V
            ]
        """
        # Boundary conditions ordered URDL
        # The top row boundary condition of True is a downward arrow
        # The left condition of True is a right arrow
        # verts = ['LR', 'LU', 'LD', 'UD', 'UR', 'RD']
        next_top = [False, False, True, True, False, True]
        next_left = [True, False, False, False, True, True]
        check_top = [True, False, True, False, False, True]
        check_left = [False, False, False, True, True, True]

        bdry = [self._bdry_cond[0]]
        lbd = list(self._bdry_cond[3]) + [None] # Dummy
        left = [ [lbd[0]] ]
        cur = [[-1]]
        n = self._nrows
        m = self._ncols
        # [[3, 1], [5, 3]]
        # [[4, 3], [3, 2]]

        while len(cur) > 0:
            # If we're at the last row and all our final boundry condition is statisfied
            if len(cur) > n and all(x is not self._bdry_cond[2][i]
                                    for i,x in enumerate(bdry[-1])):
                cur.pop()
                bdry.pop()
                left.pop()
                yield self.element_class(self, tuple(tuple(x) for x in cur))

            # Find the next row
            row = cur[-1]
            l = left[-1]
            i = len(cur) - 1
            while len(row) > 0:
                row[-1] += 1
                # Check to see if we have more vertices
                if row[-1] > 5:
                    row.pop()
                    l.pop()
                    continue
                # Check to see if we can add the vertex
                if (check_left[row[-1]] is l[-1] or l[-1] is None) \
                        and (check_top[row[-1]] is bdry[-1][len(row)-1]
                             or bdry[-1][len(row)-1] is None):
                    if len(row) != m:
                        l.append(next_left[row[-1]])
                        row.append(-1)
                    # Check the right bdry condition since we are at the rightmost entry
                    elif next_left[row[-1]] is not self._bdry_cond[1][i]:
                        bdry.append(map(lambda x: next_top[x], row))
                        cur.append([-1])
                        left.append([lbd[i+1]])
                        break

            # If we've killed this row, backup
            if len(cur[-1]) == 0:
                cur.pop()
                bdry.pop()
                left.pop()

    Element = SixVertexConfiguration

