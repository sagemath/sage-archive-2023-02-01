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
        """
        pass # We're not checking for now

    def _repr_(self):
        """
        Return a string representation of ``self``.
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

    def to_digraph(self):
        """
        Return a digraph version of ``self``.
        """

    def plot(self):
        """
        Return a plot of ``self``.
        """

class SixVertexModel(Parent, UniqueRepresentation):
    """
    The six vertex model.

    We model a configuration by indicating which configuration by the
    following six configurations which are determined by the two outgoing
    arrows:

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
    """
    @staticmethod
    def __classcall_private__(cls, n, m=None, boundary_conditions=None):
        """
        Normalize input to ensure a unique representation.
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
        """
        self._nrows = n
        self._ncols = m
        self._bdry_cond = boundary_conditions # Ordered URDL
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "The six vertex model on an {} by {} grid".format(self._nrows, self._ncols)

    def boundary_conditions(self):
        """
        Return the boundary conditions of ``self``.
        """
        return self._bdry_cond

    def _repr_option(self, key):
        """
        Metadata about the ``_repr_()`` output.

        See :meth:`sage.structure.parent._repr_option` for details.
        """
        if key == 'element_ascii_art':
            return True
        return Parent._repr_option(self, key)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.
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
                if x in verts:
                    elt[-1].append(verts.index(x))
                elif x in range(6):
                    elt[-1].append(x)
                else:
                    raise ValueError("invalid entry")
            elt[-1] = tuple(elt[-1])
        return self.element_class(self, tuple(elt))

    def __iter__(self):
        """
        Iterate through ``self``.
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
                if check_left[row[-1]] is l[-1] \
                        and check_top[row[-1]] is bdry[-1][len(row)-1]:
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

