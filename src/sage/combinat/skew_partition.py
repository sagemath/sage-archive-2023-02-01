r"""
Skew Partitions

A skew partition `skp` of size `n` is a pair of
partitions `[p_1, p_2]` where `p_1` is a
partition of the integer `n_1`, `p_2` is a
partition of the integer `n_2`, `p_2` is an inner
partition of `p_1`, and `n = n_1 - n_2`. We say
that `p_1` and `p_2` are respectively the *inner*
and *outer* partitions of `skp`.

A skew partition can be depicted by a diagram made of rows of
cells, in the same way as a partition. Only the cells of the outer
partition `p_1` which are not in the inner partition
`p_2` appear in the picture. For example, this is the
diagram of the skew partition [[5,4,3,1],[3,3,1]].

::

    sage: print SkewPartition([[5,4,3,1],[3,3,1]]).diagram()
       **
       *
     **
    *

A skew partition can be *connected*, which can easily be described
in graphic terms: for each pair of consecutive rows, there are at
least two cells (one in each row) which have a common edge. This is
the diagram of the connected skew partition [[5,4,3,1],[3,1]]::

    sage: print SkewPartition([[5,4,3,1],[3,1]]).diagram()
       **
     ***
    ***
    *
    sage: SkewPartition([[5,4,3,1],[3,1]]).is_connected()
    True

The first example of a skew partition is not a connected one.

Applying a reflection with respect to the main diagonal yields the
diagram of the *conjugate skew partition*, here
[[4,3,3,2,1],[3,3,2]]::

    sage: SkewPartition([[5,4,3,1],[3,3,1]]).conjugate()
    [[4, 3, 3, 2, 1], [3, 2, 2]]
    sage: print SkewPartition([[5,4,3,1],[3,3,1]]).conjugate().diagram()
       *
      *
      *
    **
    *

The *outer corners* of a skew partition are the corners of its
outer partition. The *inner corners* are the internal corners of
the outer partition when the inner partition is taken off. Shown
below are the coordinates of the inner and outer corners.

::

    sage: SkewPartition([[5,4,3,1],[3,3,1]]).outer_corners()
    [(0, 4), (1, 3), (2, 2), (3, 0)]
    sage: SkewPartition([[5,4,3,1],[3,3,1]]).inner_corners()
    [(0, 3), (2, 1), (3, 0)]

EXAMPLES:

There are 9 skew partitions of size 3, with no empty row nor empty
column::

    sage: SkewPartitions(3).cardinality()
    9
    sage: SkewPartitions(3).list()
    [[[3], []],
     [[2, 1], []],
     [[3, 1], [1]],
     [[2, 2], [1]],
     [[3, 2], [2]],
     [[1, 1, 1], []],
     [[2, 2, 1], [1, 1]],
     [[2, 1, 1], [1]],
     [[3, 2, 1], [2, 1]]]

There are 4 connected skew partitions of size 3::

    sage: SkewPartitions(3, overlap=1).cardinality()
    4
    sage: SkewPartitions(3, overlap=1).list()
    [[[3], []], [[2, 1], []], [[2, 2], [1]], [[1, 1, 1], []]]

This is the conjugate of the skew partition [[4,3,1],[2]]

::

    sage: SkewPartition([[4,3,1],[2]]).conjugate()
    [[3, 2, 2, 1], [1, 1]]

Geometrically, we just applied a reflection with respect to the main
diagonal on the diagram of the partition. Of course, this operation
is an involution::

    sage: SkewPartition([[4,3,1],[2]]).conjugate().conjugate()
    [[4, 3, 1], [2]]

The jacobi_trudy() method computes the Jacobi-Trudi matrix. See
Macdonald I.-G., (1995), "Symmetric Functions and Hall
Polynomials", Oxford Science Publication for a definition and
discussion.

::

    sage: SkewPartition([[4,3,1],[2]]).jacobi_trudi()
    [h[2]  h[]    0]
    [h[5] h[3]  h[]]
    [h[6] h[4] h[1]]

This example shows how to compute the corners of a skew partition.

::

    sage: SkewPartition([[4,3,1],[2]]).inner_corners()
    [(0, 2), (1, 0)]
    sage: SkewPartition([[4,3,1],[2]]).outer_corners()
    [(0, 3), (1, 2), (2, 0)]
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

import sage.combinat.composition
import sage.combinat.partition
from sage.rings.all import QQ
from sage.sets.set import Set
from sage.graphs.digraph import DiGraph
from combinat import CombinatorialClass, CombinatorialObject
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.infinity import PlusInfinity
infinity = PlusInfinity()

def SkewPartition(skp):
    """
    Returns the skew partition object corresponding to skp.

    EXAMPLES::

        sage: skp = SkewPartition([[3,2,1],[2,1]]); skp
        [[3, 2, 1], [2, 1]]
        sage: skp.inner()
        [2, 1]
        sage: skp.outer()
        [3, 2, 1]
    """
    if skp not in SkewPartitions():
        raise ValueError, "invalid skew partition: %s"%skp
    else:
        return SkewPartition_class(skp)

class SkewPartition_class(CombinatorialObject):
    def __init__(self, skp):
        """
        TESTS::

            sage: skp = SkewPartition([[3,2,1],[2,1]])
            sage: skp == loads(dumps(skp))
            True
        """
        outer = sage.combinat.partition.Partition(skp[0])
        inner = sage.combinat.partition.Partition(skp[1])
        CombinatorialObject.__init__(self, [outer, inner])

    def ferrers_diagram(self):
        """
        Returns the Ferrers diagram of self.

        EXAMPLES::

            sage: print SkewPartition([[5,4,3,1],[3,3,1]]).ferrers_diagram()
               **
               *
             **
            *
            sage: print SkewPartition([[5,4,3,1],[3,1]]).diagram()
               **
             ***
            ***
            *
            sage: Partitions.global_options(diagram_str='#', convention="French")
            sage: print SkewPartition([[5,4,3,1],[3,1]]).diagram()
            #
            ###
             ###
               ##
            sage: Partitions.global_options.reset()
        """
        # TODO: Should be SkewPartitions.global_options
        char, convention = sage.combinat.partition.Partitions.global_options('diagram_str','convention')

        if convention == "English":
            L = range(len(self[0]))
        else:
            L = reversed(range(len(self[0])))
        s = ""
        for i in L:
            if len(self[1]) > i:
                s += " "*self[1][i]
                s += char*(self[0][i]-self[1][i])
            else:
                s += char*self[0][i]
            s += "\n"
        return s[:-1]

    diagram = ferrers_diagram

    def _ascii_art_(self):
        """
        TESTS::

            sage: ascii_art(SkewPartitions(3).list())
            [                        *   *   *    * ]
            [      **   **   *    *  *   *  *    *  ]
            [ ***, * , *  , **, ** , *, * , * , *   ]
            sage: Partitions.global_options(diagram_str='#', convention="French")
            sage: ascii_art(SkewPartitions(3).list())
            [                        #  #   #   #   ]
            [      #   #    ##  ##   #   #  #    #  ]
            [ ###, ##,  ##,  #,   #, #,  #,  #,   # ]
            sage: Partitions.global_options.reset()
        """
        from sage.misc.ascii_art import AsciiArt
        return AsciiArt(self.diagram().splitlines())


    def inner(self):
        """
        Returns the inner partition of the skew partition.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).inner()
            [1, 1]
        """
        return self[1]

    def outer(self):
        """
        Returns the outer partition of the skew partition.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).outer()
            [3, 2, 1]
        """
        return self[0]

    def column_lengths(self):
        """
        Returns the column lengths of the skew partition.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).column_lengths()
            [1, 2, 1]
            sage: SkewPartition([[5,2,2,2],[2,1]]).column_lengths()
            [2, 3, 1, 1, 1]
        """
        return self.conjugate().row_lengths()

    def row_lengths(self):
        """
        Returns the row lengths of the skew partition.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).row_lengths()
            [2, 1, 1]
        """
        skp = self
        o = skp[0]
        i = skp[1]+[0]*(len(skp[0])-len(skp[1]))
        return [x[0]-x[1] for x in zip(o,i)]


    def size(self):
        """
        Returns the size of the skew partition.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[1,1]]).size()
            4
        """
        return sum(self.row_lengths())

    def is_connected(self):
        """
        Returns True if self is a connected skew partition.

        A skew partition is said to be *connected* if for each pair of
        consecutive rows, there are at least two cells (one in each row)
        which have a common edge.

        EXAMPLES::

            sage: SkewPartition([[5,4,3,1],[3,3,1]]).is_connected()
            False
            sage: SkewPartition([[5,4,3,1],[3,1]]).is_connected()
            True
        """
        return self.is_overlap(1)


    def overlap(self):
        """
        Returns the overlap of self.

        The overlap of two consecutive rows in a skew partition is the
        number of pairs of cells (one in each row) that share a common
        edge. This number can be positive, zero, or negative.

        The overlap of a skew partition is the minimum of the overlap of
        the consecutive rows, or infinity in the case of at most one row.
        If the overlap is positive, then the skew partition is called
        *connected*.

        EXAMPLES::

            sage: SkewPartition([[],[]]).overlap()
            +Infinity
            sage: SkewPartition([[1],[]]).overlap()
            +Infinity
            sage: SkewPartition([[10],[]]).overlap()
            +Infinity
            sage: SkewPartition([[10],[2]]).overlap()
            +Infinity
            sage: SkewPartition([[10,1],[2]]).overlap()
            -1
            sage: SkewPartition([[10,10],[1]]).overlap()
            9
        """
        p,q = self
        if len(p) <= 1:
            return infinity
        if q == []:
            return min(p)
        q = [q[0]] + list(q)
        return min(row_lengths_aux([p,q]))


    def is_overlap(self, n):
        r"""
        Returns True if n = self.overlap()

        .. seealso::

           :meth:`overlap`

        EXAMPLES::

            sage: SkewPartition([[5,4,3,1],[3,1]]).is_overlap(1)
            True
        """
        return n <= self.overlap()

    def conjugate(self):
        """
        Returns the conjugate of the skew partition skp.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2]]).conjugate()
            [[3, 2, 1], [1, 1]]
        """
        return SkewPartition(map(lambda x: x.conjugate(), self))


    def outer_corners(self):
        """
        Returns a list of the outer corners of the skew partition skp.

        EXAMPLES::

            sage: SkewPartition([[4, 3, 1], [2]]).outer_corners()
            [(0, 3), (1, 2), (2, 0)]
        """
        return self.outer().corners()


    def inner_corners(self):
        """
        Returns a list of the inner corners of the skew partition skp.

        EXAMPLES::

            sage: SkewPartition([[4, 3, 1], [2]]).inner_corners()
            [(0, 2), (1, 0)]
            sage: SkewPartition([[4, 3, 1], []]).inner_corners()
            [(0, 0)]
        """
        inner = self.inner()
        outer = self.outer()
        if inner == []:
            if outer == []:
                return []
            else:
                return [(0,0)]
        icorners = [(0, inner[0])]
        nn = len(inner)
        for i in range(1,nn):
            if inner[i] != inner[i-1]:
                icorners += [ (i, inner[i]) ]

        icorners += [(nn, 0)]
        return icorners

    def cells(self):
        """
        Returns the coordinates of the cells of self. Coordinates are given
        as (row-index, column-index) and are 0 based.

        EXAMPLES::

            sage: SkewPartition([[4, 3, 1], [2]]).cells()
            [(0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (2, 0)]
            sage: SkewPartition([[4, 3, 1], []]).cells()
            [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (2, 0)]
            sage: SkewPartition([[2], []]).cells()
            [(0, 0), (0, 1)]
        """
        outer = self.outer()
        inner = self.inner()[:]
        inner += [0]*(len(outer)-len(inner))
        res = []
        for i in range(len(outer)):
            for j in range(inner[i], outer[i]):
                res.append( (i,j) )
        return res

    def to_list(self):
        """
        Returns self as a list of lists.

        EXAMPLES::

            sage: s = SkewPartition([[4,3,1],[2]])
            sage: s.to_list()
            [[4, 3, 1], [2]]
            sage: type(s.to_list())
            <type 'list'>
        """
        return map(list, list(self))


    def to_dag(self):
        """
        Returns a directed acyclic graph corresponding to the skew
        partition.

        EXAMPLES::

            sage: dag = SkewPartition([[3, 2, 1], [1, 1]]).to_dag()
            sage: dag.edges()
            [('0,1', '0,2', None), ('0,1', '1,1', None)]
            sage: dag.vertices()
            ['0,1', '0,2', '1,1', '2,0']
        """
        i = 0

        #Make the skew tableau from the shape
        skew = [[1]*row_length for row_length in self.outer()]
        inner = self.inner()
        for i in range(len(inner)):
            for j in range(inner[i]):
                skew[i][j] = None

        G = DiGraph()
        for row in range(len(skew)):
            for column in range(len(skew[row])):
                if skew[row][column] is not None:
                    string = "%d,%d" % (row, column)
                    G.add_vertex(string)
                    #Check to see if there is a node to the right
                    if column != len(skew[row]) - 1:
                        newstring = "%d,%d" % (row, column+1)
                        G.add_edge(string, newstring)

                    #Check to see if there is anything below
                    if row != len(skew) - 1:
                        if len(skew[row+1]) > column:
                            if skew[row+1][column] is not None:
                                newstring = "%d,%d" % (row+1, column)
                                G.add_edge(string, newstring)
        return G

    def quotient(self, k):
        """
        The quotient map extended to skew partitions.

        EXAMPLES::

            sage: SkewPartition([[3, 3, 2, 1], [2, 1]]).quotient(2)
            [[[3], []], [[], []]]
        """
        ## k-th element is the skew partition built using the k-th partition of the
        ## k-quotient of the outer and the inner partition.
        ## This bijection is only defined if the inner and the outer partition
        ## have the same core
        if self.inner().core(k) == self.outer().core(k):
            rqinner = self.inner().quotient(k)
            rqouter = self.outer().quotient(k)
            return [ SkewPartition_class([rqouter[i],rqinner[i]]) for i in range(k) ]
        else:
            raise ValueError, "quotient map is only defined for skew partitions with inner and outer partitions having the same core"

    def rows_intersection_set(self):
        """
        Returns the set of cells in the lines of lambda which intersect the
        skew partition.

        EXAMPLES::

            sage: skp = SkewPartition([[3,2,1],[2,1]])
            sage: cells = Set([ (0,0), (0, 1), (0,2), (1, 0), (1, 1), (2, 0)])
            sage: skp.rows_intersection_set() == cells
            True
        """
        res = []
        outer = self.outer()
        inner = self.inner()
        inner += [0] * int(len(outer)-len(inner))

        for i in range(len(outer)):
            for j in range(outer[i]):
                if outer[i] != inner[i]:
                    res.append((i,j))
        return Set(res)

    def columns_intersection_set(self):
        """
        Returns the set of cells in the lines of lambda which intersect the
        skew partition.

        EXAMPLES::

            sage: skp = SkewPartition([[3,2,1],[2,1]])
            sage: cells = Set([ (0,0), (0, 1), (0,2), (1, 0), (1, 1), (2, 0)])
            sage: skp.columns_intersection_set() == cells
            True
        """
        res = [ (x[1], x[0]) for x in self.conjugate().rows_intersection_set()]
        return Set(res)


    def pieri_macdonald_coeffs(self):
        """
        Computation of the coefficients which appear in the Pieri formula
        for Macdonald polynomials given in his book ( Chapter 6.6 formula
        6.24(ii) )

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2,1]]).pieri_macdonald_coeffs()
            1
            sage: SkewPartition([[3,2,1],[2,2]]).pieri_macdonald_coeffs()
            (q^2*t^3 - q^2*t - t^2 + 1)/(q^2*t^3 - q*t^2 - q*t + 1)
            sage: SkewPartition([[3,2,1],[2,2,1]]).pieri_macdonald_coeffs()
            (q^6*t^8 - q^6*t^6 - q^4*t^7 - q^5*t^5 + q^4*t^5 - q^3*t^6 + q^5*t^3 + 2*q^3*t^4 + q*t^5 - q^3*t^2 + q^2*t^3 - q*t^3 - q^2*t - t^2 + 1)/(q^6*t^8 - q^5*t^7 - q^5*t^6 - q^4*t^6 + q^3*t^5 + 2*q^3*t^4 + q^3*t^3 - q^2*t^2 - q*t^2 - q*t + 1)
            sage: SkewPartition([[3,3,2,2],[3,2,2,1]]).pieri_macdonald_coeffs()
            (q^5*t^6 - q^5*t^5 + q^4*t^6 - q^4*t^5 - q^4*t^3 + q^4*t^2 - q^3*t^3 - q^2*t^4 + q^3*t^2 + q^2*t^3 - q*t^4 + q*t^3 + q*t - q + t - 1)/(q^5*t^6 - q^4*t^5 - q^3*t^4 - q^3*t^3 + q^2*t^3 + q^2*t^2 + q*t - 1)
        """

        set_prod = self.rows_intersection_set() - self.columns_intersection_set()
        res = 1
        for s in set_prod:
            res *= self.inner().arms_legs_coeff(s[0],s[1])
            res /= self.outer().arms_legs_coeff(s[0],s[1])
        return res

    def k_conjugate(self, k):
        """
        Returns the k-conjugate of the skew partition.

        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2,1]]).k_conjugate(3)
            [[2, 1, 1, 1, 1], [2, 1]]
            sage: SkewPartition([[3,2,1],[2,1]]).k_conjugate(4)
            [[2, 2, 1, 1], [2, 1]]
            sage: SkewPartition([[3,2,1],[2,1]]).k_conjugate(5)
            [[3, 2, 1], [2, 1]]
        """
        return SkewPartition([ self.outer().k_conjugate(k), self.inner().k_conjugate(k) ])

    def jacobi_trudi(self):
        """
        EXAMPLES::

            sage: SkewPartition([[3,2,1],[2,1]]).jacobi_trudi()
            [h[1]    0    0]
            [h[3] h[1]    0]
            [h[5] h[3] h[1]]
            sage: SkewPartition([[4,3,2],[2,1]]).jacobi_trudi()
            [h[2]  h[]    0]
            [h[4] h[2]  h[]]
            [h[6] h[4] h[2]]
        """
        p = self.outer()
        q = self.inner()
        from sage.combinat.sf.sf import SymmetricFunctions
        if len(p) == 0 and len(q) == 0:
            return MatrixSpace(SymmetricFunctions(QQ).homogeneous(), 0)(0)
        nn = len(p)
        h = SymmetricFunctions(QQ).homogeneous()
        H = MatrixSpace(h, nn)

        q  = q + [0]*int(nn-len(q))
        m = []
        for i in range(1,nn+1):
            row = []
            for j in range(1,nn+1):
                v = p[j-1]-q[i-1]-j+i
                if v < 0:
                    row.append(h(0))
                elif v == 0:
                    row.append(h([]))
                else:
                    row.append(h([v]))
            m.append(row)
        return H(m)


def from_row_and_column_length(rowL, colL):
    """
    Construct a partition from its row lengths and column lengths.

    INPUT:

     - ``rowL`` -- a composition or a list of positive integers
     - ``colL`` -- a composition or a list of positive integers

    OUTPUT:

     - If it exists the unique skew-partitions with row lengths ``rowL``
       and column lengths ``colL``.
     - Raise a ``ValueError`` if ``rowL`` and ``colL`` are not compatible.

    EXAMPLES::

       sage: from sage.combinat.skew_partition import from_row_and_column_length
       sage: print from_row_and_column_length([3,1,2,2],[2,3,1,1,1]).diagram()
         ***
        *
       **
       **
       sage: from_row_and_column_length([],[])
       [[], []]
       sage: from_row_and_column_length([1],[1])
       [[1], []]
       sage: from_row_and_column_length([2,1],[2,1])
       [[2, 1], []]
       sage: from_row_and_column_length([1,2],[1,2])
       [[2, 2], [1]]
       sage: from_row_and_column_length([1,2],[1,3])
       Traceback (most recent call last):
       ...
       ValueError: Sum mismatch : [1, 2] and [1, 3]
       sage: from_row_and_column_length([3,2,1,2],[2,3,1,1,1])
       Traceback (most recent call last):
       ...
       ValueError: Incompatible row and column length : [3, 2, 1, 2] and [2, 3, 1, 1, 1]

    .. warning::

       If some rows and columns have length zero, there is no way to retrieve
       unambiguously the skew partition. We therefore raise a ValueError. For
       examples here are two skew partitions with the same row and column
       lengths::

           sage: skp1 = SkewPartition([[2,2],[2,2]])
           sage: skp2 = SkewPartition([[2,1],[2,1]])
           sage: skp1.row_lengths(), skp1.column_lengths()
           ([0, 0], [0, 0])
           sage: skp2.row_lengths(), skp2.column_lengths()
           ([0, 0], [0, 0])
           sage: from_row_and_column_length([0,0], [0,0])
           Traceback (most recent call last):
           ...
           ValueError: row and column length must be positive

    TESTS::

       sage: all(from_row_and_column_length(p.row_lengths(), p.column_lengths()) == p
       ...         for i in range(8) for p in SkewPartitions(i))
       True
    """
    if sum(rowL) != sum(colL):
        raise ValueError, "Sum mismatch : %s and %s"%(rowL, colL)
    if not all(i>0 for i in rowL) or not all(i>0 for i in colL):
        raise ValueError, "row and column length must be positive"
    if rowL == []:
        return SkewPartition([[],[]])
    colL_new = colL[:]
    resIn  = []
    resOut = []
    inPOld = len(colL)
    for row in rowL:
        inP = len(colL_new) - row
        if inP < 0 or inP > inPOld:
            raise ValueError, "Incompatible row and column length : %s and %s"%(rowL, colL)
        inPOld = inP
        resIn.append(inP)
        resOut.append(len(colL_new))
        for iCol in range(inP, len(colL_new)):
            colL_new[iCol] -= 1;
            if colL_new[iCol] < 0:
                raise ValueError, "Incompatible row and column length : %s and %s"%(rowL, colL)
        while colL_new != [] and colL_new[-1] == 0:
            colL_new.pop()
    return SkewPartition([resOut, filter(lambda x:x, resIn)])


def row_lengths_aux(skp):
    """
    EXAMPLES::

        sage: from sage.combinat.skew_partition import row_lengths_aux
        sage: row_lengths_aux([[5,4,3,1],[3,3,1]])
        [2, 1, 2]
        sage: row_lengths_aux([[5,4,3,1],[3,1]])
        [2, 3]
    """
    if skp[0] == []:
        return []
    else:
        return map(lambda x: x[0] - x[1], zip(skp[0], skp[1]))

def SkewPartitions(n=None, row_lengths=None, overlap=0):
    """
    Returns the combinatorial class of skew partitions.

    EXAMPLES::

        sage: SkewPartitions(4)
        Skew partitions of 4
        sage: SkewPartitions(4).cardinality()
        28
        sage: SkewPartitions(row_lengths=[2,1,2])
        Skew partitions with row lengths [2, 1, 2]
        sage: SkewPartitions(4, overlap=2)
        Skew partitions of 4 with overlap of 2
        sage: SkewPartitions(4, overlap=2).list()
        [[[4], []], [[2, 2], []]]
    """
    number_of_arguments = 0
    for arg in ['n', 'row_lengths']:
        if locals()[arg] is not None:
            number_of_arguments += 1

    if number_of_arguments > 1:
        raise ValueError, "you can only specify one of n or row_lengths"

    if number_of_arguments == 0:
        return SkewPartitions_all()

    if n is not None:
        return SkewPartitions_n(n, overlap)
    else:
        return SkewPartitions_rowlengths(row_lengths, overlap)

class SkewPartitions_all(CombinatorialClass):
    def __init__(self):
        """
        TESTS::

            sage: S = SkewPartitions()
            sage: S == loads(dumps(S))
            True
        """
        pass

    Element = SkewPartition_class

    def __contains__(self, x):
        """
        TESTS::

            sage: [[], []] in SkewPartitions()
            True
            sage: [[], [1]] in SkewPartitions()
            False
            sage: [[], [-1]] in SkewPartitions()
            False
            sage: [[], [0]] in SkewPartitions()
            False
            sage: [[3,2,1],[]] in SkewPartitions()
            True
            sage: [[3,2,1],[1]] in SkewPartitions()
            True
            sage: [[3,2,1],[2]] in SkewPartitions()
            True
            sage: [[3,2,1],[3]] in SkewPartitions()
            True
            sage: [[3,2,1],[4]] in SkewPartitions()
            False
            sage: [[3,2,1],[1,1]] in SkewPartitions()
            True
            sage: [[3,2,1],[1,2]] in SkewPartitions()
            False
            sage: [[3,2,1],[2,1]] in SkewPartitions()
            True
            sage: [[3,2,1],[2,2]] in SkewPartitions()
            True
            sage: [[3,2,1],[3,2]] in SkewPartitions()
            True
            sage: [[3,2,1],[1,1,1]] in SkewPartitions()
            True
            sage: [[7, 4, 3, 2], [8, 2, 1]] in SkewPartitions()
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions()
            True
            sage: [[4,2,1],[1,1,1,1]] in SkewPartitions()
            False
        """
        if isinstance(x, SkewPartition_class):
            return True

        try:
            if len(x) != 2:
                return False
        except TypeError:
            return False

        p = sage.combinat.partition.Partitions()
        if x[0] not in p:
            return False
        if x[1] not in p:
            return False

        if not sage.combinat.partition.Partition(x[0]).contains(x[1]):
            return False

        return True

    def __repr__(self):
        """
        TESTS::

            sage: repr(SkewPartitions())
            'Skew partitions'
        """
        return "Skew partitions"


    def _latex_(self):
        """
        LaTeX output as a young diagram.

        EXAMPLES::

            sage: P = Partition([2, 1])
            sage: print P._latex_young_diagram()
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }
        """
        if len(self._list) == 0:
            return "{\\emptyset}"

        from sage.combinat.output import tex_from_array
        return tex_from_array([ ["\\phantom{x}"]*row_size for row_size in self._list ])

    def list(self):
        """
        TESTS::

            sage: SkewPartitions().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class SkewPartitions_n(CombinatorialClass):
    """
    The combinatorial class of skew partitions with given size (and
    horizontal minimal overlap).
    """
    def __init__(self, n, overlap=0):
        """
        INPUT:

         - ``n`` -- a non-negative integer
         - ``overlap`` -- an integer

        Returns the set of the skew partitions of ``n`` with overlap
        at least ``overlap``, and no empty row.

        The iteration order is not specified yet.

        Caveat: this set is stable under conjugation only for overlap=
        0 or 1. What exactly happens for negative overlaps is not yet
        well specified, and subject to change (we may want to
        introduce vertical overlap constraints as well). ``overlap`` would
        also better be named ``min_overlap``.

        Todo: as is, this set is essentially the composition of
        ``Compositions(n)`` (which give the row lengths) and
        ``SkewPartition(n, row_lengths=...)``, and one would want to
        "inherit" list and cardinality from this composition.

        TESTS::

            sage: S = SkewPartitions(3)
            sage: S == loads(dumps(S))
            True
            sage: S = SkewPartitions(3, overlap=1)
            sage: S == loads(dumps(S))
            True
        """
        self.n = n
        if overlap == 'connected':
            self.overlap = 1
        else:
            self.overlap = overlap

    Element = SkewPartition_class

    def __contains__(self, x):
        """
        TESTS::

            sage: [[],[]] in SkewPartitions(0)
            True
            sage: [[3,2,1], []] in SkewPartitions(6)
            True
            sage: [[3,2,1], []] in SkewPartitions(7)
            False
            sage: [[3,2,1], []] in SkewPartitions(5)
            False
            sage: [[7, 4, 3, 2], [8, 2, 1]] in SkewPartitions(8)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(5)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(5, overlap=-1)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap=-1)
            True
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap=0)
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap='connected')
            False
            sage: [[7, 4, 3, 2], [5, 2, 1]] in SkewPartitions(8, overlap=-2)
            True
        """
        return x in SkewPartitions() and sum(x[0])-sum(x[1]) == self.n and self.overlap <= SkewPartition(x).overlap()


    def __repr__(self):
        """
        TESTS::

            sage: repr(SkewPartitions(3))
            'Skew partitions of 3'
            sage: repr(SkewPartitions(3, overlap=1))
            'Skew partitions of 3 with overlap of 1'
        """
        string = "Skew partitions of %s"%self.n
        if self.overlap:
            string += " with overlap of %s"%self.overlap
        return string


    def _count_slide(self, co, overlap=0):
        """
        Returns the number of skew partitions related to the composition co
        by 'sliding'. The composition co is the list of row lengths of the
        skew partition.

        Todo: this should be a method of SkewPartition

        EXAMPLES::

            sage: s = SkewPartitions(3)
            sage: s._count_slide([2,1])
            2
            sage: [ sp for sp in s if sp.row_lengths() == [2,1] ]
            [[[2, 1], []], [[3, 1], [1]]]
            sage: s = SkewPartitions(3, overlap=1)
            sage: s._count_slide([2,1], overlap=1)
            1
            sage: [ sp for sp in s if sp.row_lengths() == [2,1] ]
            [[[2, 1], []]]
        """
        nn = len(co)
        result = 1
        for i in range(nn-1):
            comb    = min(co[i], co[i+1])
            comb   += 1 - overlap
            result *= comb

        return result

    def cardinality(self):
        """
        Returns the number of skew partitions of the integer n.

        EXAMPLES::

            sage: SkewPartitions(0).cardinality()
            1
            sage: SkewPartitions(4).cardinality()
            28
            sage: SkewPartitions(5).cardinality()
            87
            sage: SkewPartitions(4, overlap=1).cardinality()
            9
            sage: SkewPartitions(5, overlap=1).cardinality()
            20
            sage: s = SkewPartitions(5, overlap=-1)
            sage: s.cardinality() == len(s.list())
            True
        """
        n = self.n
        overlap = self.overlap

        if n == 0:
            return 1

        if overlap > 0:
            gg = sage.combinat.composition.Compositions(n, min_part = max(1, overlap))
        else:
            gg = sage.combinat.composition.Compositions(n)

        sum_a = 0
        for co in gg:
            sum_a += self._count_slide(co, overlap=overlap)

        return sum_a

    def __iter__(self):
        """
        Returns a list of the skew partitions of n.

        EXAMPLES::

            sage: SkewPartitions(3).list()
            [[[3], []],
            [[2, 1], []],
            [[3, 1], [1]],
            [[2, 2], [1]],
            [[3, 2], [2]],
            [[1, 1, 1], []],
            [[2, 2, 1], [1, 1]],
            [[2, 1, 1], [1]],
            [[3, 2, 1], [2, 1]]]


            sage: SkewPartitions(3, overlap=0).list()
            [[[3], []],
            [[2, 1], []],
            [[3, 1], [1]],
            [[2, 2], [1]],
            [[3, 2], [2]],
            [[1, 1, 1], []],
            [[2, 2, 1], [1, 1]],
            [[2, 1, 1], [1]],
            [[3, 2, 1], [2, 1]]]
            sage: SkewPartitions(3, overlap=1).list()
            [[[3], []], [[2, 1], []], [[2, 2], [1]], [[1, 1, 1], []]]
            sage: SkewPartitions(3, overlap=2).list()
            [[[3], []]]
            sage: SkewPartitions(3, overlap=3).list()
            [[[3], []]]
            sage: SkewPartitions(3, overlap=4).list()
            []
        """
        n = self.n
        overlap = self.overlap

        for co in sage.combinat.composition.Compositions(n, min_part = max(1, overlap)):
            for sp in SkewPartitions(row_lengths=co, overlap=overlap):
                yield sp

######################################
# Skew Partitions (from row lengths) #
######################################
class SkewPartitions_rowlengths(CombinatorialClass):
    """
    The combinatorial class of all skew partitions with given row
    lengths.
    """
    def __init__(self, co, overlap=0):
        """
        TESTS::

            sage: S = SkewPartitions(row_lengths=[2,1])
            sage: S == loads(dumps(S))
            True
        """
        self.co = co
        if overlap == 'connected':
            self.overlap = 1
        else:
            self.overlap = overlap

    Element = SkewPartition_class

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[4,3,1],[2]] in SkewPartitions(row_lengths=[2,3,1])
            True
            sage: [[4,3,1],[2]] in SkewPartitions(row_lengths=[2,1,3])
            False
            sage: [[5,4,3,1],[3,3,1]] in SkewPartitions(row_lengths=[2,1,1,2])
            False
            sage: [[5,4,3,1],[3,3,1]] in SkewPartitions(row_lengths=[2,1,2,1])
            True
        """
        valid = x in SkewPartitions()
        if valid:
            o = x[0]
            i = x[1]+[0]*(len(x[0])-len(x[1]))
            return [x[0]-x[1] for x in zip(o,i)] == self.co


    def __repr__(self):
        """
        TESTS::

            sage: repr(SkewPartitions(row_lengths=[2,1]))
            'Skew partitions with row lengths [2, 1]'
        """
        return "Skew partitions with row lengths %s"%self.co

    def _from_row_lengths_aux(self, sskp, ck_1, ck, overlap=0):
        """
        EXAMPLES::

            sage: s = SkewPartitions(row_lengths=[2,1])
            sage: s._from_row_lengths_aux([[1], []], 1, 1, overlap=0)
            [[[1, 1], []], [[2, 1], [1]]]
            sage: s._from_row_lengths_aux([[1, 1], []], 1, 1, overlap=0)
            [[[1, 1, 1], []], [[2, 2, 1], [1, 1]]]
            sage: s._from_row_lengths_aux([[2, 1], [1]], 1, 1, overlap=0)
            [[[2, 1, 1], [1]], [[3, 2, 1], [2, 1]]]
            sage: s._from_row_lengths_aux([[1], []], 1, 2, overlap=0)
            [[[2, 2], [1]], [[3, 2], [2]]]
            sage: s._from_row_lengths_aux([[2], []], 2, 1, overlap=0)
            [[[2, 1], []], [[3, 1], [1]]]
        """
        lskp = []
        nn = min(ck_1, ck)
        mm = max(0, ck-ck_1)
        # nn should be >= 0.  In the case of the positive overlap,
        # the min_part condition insures ck>=overlap for all k

        nn -= overlap
        for i in range(nn+1):
            (skp1, skp2) = sskp
            skp2 += [0]*(len(skp1)-len(skp2))
            skp1 = map(lambda x: x + i + mm, skp1)
            skp1 += [ck]
            skp2 = map(lambda x: x + i + mm, skp2)
            skp2 = filter(lambda x: x != 0, skp2)
            lskp += [ SkewPartition([skp1, skp2]) ]

        return lskp


    def list(self):
        """
        Returns a list of all the skew partitions that have row lengths
        given by the composition self.co.

        EXAMPLES::

            sage: SkewPartitions(row_lengths=[2,2]).list()
            [[[2, 2], []], [[3, 2], [1]], [[4, 2], [2]]]
            sage: SkewPartitions(row_lengths=[2,2], overlap=1).list()
            [[[2, 2], []], [[3, 2], [1]]]
        """
        co = self.co
        overlap = self.overlap

        if co == []:
            return [ SkewPartition([[],[]]) ]

        nn = len(co)
        if nn == 1:
            return [ SkewPartition([[co[0]],[]]) ]

        result = []
        for sskp in SkewPartitions(row_lengths=co[:-1], overlap=overlap):
            result += self._from_row_lengths_aux(sskp, co[-2], co[-1], overlap)
        return result



