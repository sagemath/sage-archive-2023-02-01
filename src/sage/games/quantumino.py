# coding=utf-8
r"""
Family Games America's Quantumino solver

This module allows to solve the Quantumino puzzle made by Family Games
America (see [1] and [2]). This puzzle was left at the dinner room of the
Laboratoire de Combinatoire Informatique Mathématique in Montréal by Franco
Saliola during winter 2011.

The solution uses the dancing links code which is in Sage. Dancing links
were originally introduced by Donald Knuth in 2000 [3]. In particular,
Knuth used dancing links to solve tilings of a region by 2D pentaminos.
Here we extend the method for 3D pentaminos.

AUTHOR:

    - Sébastien Labbé, April 28th, 2011

DESCRIPTION (from [1]):

    "
    Pentamino games have been taken to a whole different level; a 3-D
    level, with this colorful creation! Using the original pentamino
    arrangements of 5 connected squares which date from 1907, players are
    encouraged to «think inside the box» as they try to fit 16 of the 17
    3-D pentamino pieces inside the playing perimeters. Remove a different
    piece each time you play for an entirely new challenge! Thousands of
    solutions to be found!
    Quantumino™ hands-on educational tool where players learn how shapes
    can be transformed or arranged into predefined shapes and spaces.
    Includes:
    1 wooden frame, 17 wooden blocks, instruction booklet.
    Age: 8+
    "

EXAMPLES:

Here are the 17 wooden blocks of the Quantumino puzzle numbered from 0 to
16 in the following 3d picture::

    sage: from sage.games.quantumino import quantumino_solver
    sage: quantumino_solver.show_pentaminos()

To solve the puzzle without using the block numbered 12::

    sage: s = quantumino_solver.get_solution(12)
    sage: s
    Quantumino solution without the following pentamino :
    3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)], Color: blue
    sage: s.show3d()      # long time (<1s)

To solve the puzzle without using the block numbered 7::

    sage: s = quantumino_solver.get_solution(7)
    sage: s
    Quantumino solution without the following pentamino :
    3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
    sage: s.show3d()      # long time (<1s)

To get all the solutions, use the iterator. Note that finding the first
solution is the most time consuming because it needs to creates the
complete data to describe the problem::

    sage: it = quantumino_solver.solutions_iterator(7)
    sage: it.next()                                     # (1s)
    Quantumino solution without the following pentamino :
    3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
    sage: it.next()                                     # (0.001s)
    Quantumino solution without the following pentamino :
    3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
    sage: it.next()                                     # (0.001s)
    Quantumino solution without the following pentamino :
    3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange

To get the solution inside other boxes::

    sage: s = quantumino_solver.get_solution(7, box=(4,4,5))        # long time (2s)
    sage: s.show3d()                                                # long time (<1s)

::

    sage: s = quantumino_solver.get_solution(7, box=(2,2,20))       # long time (1s)
    sage: s.show3d()                                                # long time (<1s)

EXAMPLES:

One can also solve similar problems by defining 3D polyominos directly. The
following is a 2d puzzle owned by Florent Hivert::

    sage: from sage.games.quantumino import Polyomino3d, solve_3d_puzzle
    sage: L = []
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,2,0),(1,3,0)]))
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,2,0)]))
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,3,0)]))
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,3,0)]))
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0)]))
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,1,0),(1,2,0)]))
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,1,0),(1,3,0)]))
    sage: L.append(Polyomino3d([(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,3,0)]))
    sage: L.append(Polyomino3d([(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,2,0)]))
    sage: L.append(Polyomino3d([(0,0,0),(0,1,0),(0,2,0),(1,0,0),(1,1,0),(1,2,0)]))

Solve the puzzle and show the solution::

    sage: it = solve_3d_puzzle(L, (8,8,1))
    sage: solution = it.next()
    sage: G = sum([piece.show3d() for piece in solution], Graphics())
    sage: G.show(aspect_ratio=1, viewer='tachyon')

REFERENCES:

    - [1] http://familygamesamerica.com/mainsite/consumers/productview.php?pro_id=274&search=quantumino
    - [2] Quantumino - How to Play, http://www.youtube.com/watch?v=jX_VKzakZi8
    - [3] Knuth, Donald (2000). "Dancing links". arXiv:cs/0011047.

"""
#*****************************************************************************
#       Copyright (C) 2011 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
from sage.combinat.all import WeylGroup
from sage.plot.plot import Graphics
from sage.modules.free_module_element import vector
from sage.plot.plot3d.platonic import cube
from sage.plot.plot3d.shapes2 import text3d
from sage.structure.sage_object import SageObject

##############################
# Orthogonal Rotation matrices
##############################
rotation_matrices = [w.matrix() for w in WeylGroup(['B',3]) if w.matrix().det() == 1]

##############################
# Class Polyomino3d
##############################
class Polyomino3d(SageObject):
    r"""
    Return the 3D polyomino defined by a set of 3d coordinates.

    INPUT:

    - ``coords`` - iterable of 3d tuple
    - ``color`` - string (optional, default: ``'gray'``), the color

    EXAMPLES::

        sage: from sage.games.quantumino import Polyomino3d
        sage: Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
        3D Polyomino: [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)], Color: blue
    """
    def __init__(self, coords, color='gray'):
        r"""
        INPUT:

        - ``coords`` - iterable of 3d tuple
        - ``color`` - string (optional, default: ``'gray'``), the color

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            3D Polyomino: [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)], Color: blue
        """
        assert all(len(a) == 3 for a in coords), "coord must be in dimension 3"
        assert isinstance(color, str)
        self._blocs = frozenset(tuple(c) for c in coords)
        self._color = color

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(0)
            sage: p
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
        """
        s = "3D Polyomino: %s, " % sorted(self._blocs)
        s += "Color: %s" % self._color
        return s

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: hash(p)
            2059134902
        """
        return hash(self._blocs)

    def __eq__(self, other):
        r"""
        Return whether self is equal to other.

        INPUT:

        - ``other`` - a polyomino

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: q = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='red')
            sage: p == q
            True
            sage: r = Polyomino3d([(0,0,0), (0,1,0), (1,1,0)], color='blue')
            sage: p == r
            False
        """
        return self._blocs == other._blocs

    def __ne__(self, other):
        r"""
        Return whether self is not equal to other.

        INPUT:

        - ``other`` - a polyomino

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: q = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='red')
            sage: p != q
            False
            sage: r = Polyomino3d([(0,0,0), (0,1,0), (1,1,0)], color='blue')
            sage: p != r
            True
        """
        return self._blocs != other._blocs

    def color(self):
        r"""
        Return the color of the polyomino.

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: p.color()
            'blue'
        """
        return self._color

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: it = iter(p)
            sage: it.next()
            (1, 1, 0)
        """
        return iter(self._blocs)

    def rotated(self):
        r"""
        Iterator over the 24 orthogonal rotated image of self.

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: L = list(p.rotated())
            sage: len(L)
            24
        """
        for m in rotation_matrices:
            L = [m * vector(p) for p in self]
            yield Polyomino3d(L, color=self._color)

    def rotated_canonical(self):
        r"""
        Iterator over the 24 orthogonal rotated image of self where the
        coordinates are all positive and minimal.

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: L = list(p.rotated_canonical())
            sage: len(L)
            24

        They might not be all different::

            sage: s = set(p.rotated_canonical())
            sage: len(s)
            12
        """
        for q in self.rotated():
            yield q.canonical()

    def canonical(self):
        r"""
        Returns the translated copy of self having minimal and positive
        coordinates

        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(0)
            sage: p
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: p.canonical()
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink

        TESTS::

            sage: p
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: p + (3,4,5)
            3D Polyomino: [(3, 4, 5), (4, 4, 5), (4, 5, 5), (4, 5, 6), (4, 6, 5)], Color: deeppink
            sage: (p + (3,4,5)).canonical()
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
        """
        minxyz,maxxyz = self.bounding_box()
        mx, my, mz = minxyz
        return self - (mx,my,mz)

    def is_higher_than(self, k):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(15)
            sage: p.is_higher_than(2)
            False
        """
        [(a,b,c), (aa,bb,cc)] = self.bounding_box()
        return cc - c > k-1

    def __sub__(self, arg):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(0)
            sage: p - (2,2,2)
            3D Polyomino: [(-2, -2, -2), (-1, -2, -2), (-1, -1, -2), (-1, -1, -1), (-1, 0, -2)], Color: deeppink

        """
        x,y,z = arg
        return Polyomino3d([(a-x, b-y, c-z) for a,b,c in self], color=self._color)

    def __add__(self, arg):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(0)
            sage: p + (2,2,2)
            3D Polyomino: [(2, 2, 2), (3, 2, 2), (3, 3, 2), (3, 3, 3), (3, 4, 2)], Color: deeppink
        """
        x,y,z = arg
        return Polyomino3d([(a+x, b+y, c+z) for a,b,c in self], color=self._color)

    def bounding_box(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(15)
            sage: p.bounding_box()
            [(0, 0, 0), (1, 2, 1)]
        """
        xx,yy,zz = zip(*self)
        minx = min(xx)
        miny = min(yy)
        minz = min(zz)
        maxx = max(xx)
        maxy = max(yy)
        maxz = max(zz)
        return [(minx,miny,minz), (maxx,maxy,maxz)]

    def translated(self, box):
        r"""
        Returns an iterator over the translated images of self inside a
        box.

        INPUT:

        - ``box`` - tuple of size three, size of the box

        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(0)
            sage: for t in p.translated(box=(5,8,2)): t
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            3D Polyomino: [(0, 1, 0), (1, 1, 0), (1, 2, 0), (1, 2, 1), (1, 3, 0)], Color: deeppink
            3D Polyomino: [(0, 2, 0), (1, 2, 0), (1, 3, 0), (1, 3, 1), (1, 4, 0)], Color: deeppink
            3D Polyomino: [(0, 3, 0), (1, 3, 0), (1, 4, 0), (1, 4, 1), (1, 5, 0)], Color: deeppink
            3D Polyomino: [(0, 4, 0), (1, 4, 0), (1, 5, 0), (1, 5, 1), (1, 6, 0)], Color: deeppink
            3D Polyomino: [(0, 5, 0), (1, 5, 0), (1, 6, 0), (1, 6, 1), (1, 7, 0)], Color: deeppink
            3D Polyomino: [(1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 1), (2, 2, 0)], Color: deeppink
            3D Polyomino: [(1, 1, 0), (2, 1, 0), (2, 2, 0), (2, 2, 1), (2, 3, 0)], Color: deeppink
            3D Polyomino: [(1, 2, 0), (2, 2, 0), (2, 3, 0), (2, 3, 1), (2, 4, 0)], Color: deeppink
            3D Polyomino: [(1, 3, 0), (2, 3, 0), (2, 4, 0), (2, 4, 1), (2, 5, 0)], Color: deeppink
            3D Polyomino: [(1, 4, 0), (2, 4, 0), (2, 5, 0), (2, 5, 1), (2, 6, 0)], Color: deeppink
            3D Polyomino: [(1, 5, 0), (2, 5, 0), (2, 6, 0), (2, 6, 1), (2, 7, 0)], Color: deeppink
            3D Polyomino: [(2, 0, 0), (3, 0, 0), (3, 1, 0), (3, 1, 1), (3, 2, 0)], Color: deeppink
            3D Polyomino: [(2, 1, 0), (3, 1, 0), (3, 2, 0), (3, 2, 1), (3, 3, 0)], Color: deeppink
            3D Polyomino: [(2, 2, 0), (3, 2, 0), (3, 3, 0), (3, 3, 1), (3, 4, 0)], Color: deeppink
            3D Polyomino: [(2, 3, 0), (3, 3, 0), (3, 4, 0), (3, 4, 1), (3, 5, 0)], Color: deeppink
            3D Polyomino: [(2, 4, 0), (3, 4, 0), (3, 5, 0), (3, 5, 1), (3, 6, 0)], Color: deeppink
            3D Polyomino: [(2, 5, 0), (3, 5, 0), (3, 6, 0), (3, 6, 1), (3, 7, 0)], Color: deeppink
            3D Polyomino: [(3, 0, 0), (4, 0, 0), (4, 1, 0), (4, 1, 1), (4, 2, 0)], Color: deeppink
            3D Polyomino: [(3, 1, 0), (4, 1, 0), (4, 2, 0), (4, 2, 1), (4, 3, 0)], Color: deeppink
            3D Polyomino: [(3, 2, 0), (4, 2, 0), (4, 3, 0), (4, 3, 1), (4, 4, 0)], Color: deeppink
            3D Polyomino: [(3, 3, 0), (4, 3, 0), (4, 4, 0), (4, 4, 1), (4, 5, 0)], Color: deeppink
            3D Polyomino: [(3, 4, 0), (4, 4, 0), (4, 5, 0), (4, 5, 1), (4, 6, 0)], Color: deeppink
            3D Polyomino: [(3, 5, 0), (4, 5, 0), (4, 6, 0), (4, 6, 1), (4, 7, 0)], Color: deeppink
        """
        box_x, box_y, box_z = box
        cano = self.canonical()
        x,y,z = zip(*cano)
        mx = max(x)
        my = max(y)
        mz = max(z)
        for i in range(box_x-mx):
            for j in range(box_y-my):
                for k in range(box_z-mz):
                    yield self + (i,j,k)

    def translated_rotated(self, box):
        r"""
        Return the translated and rotated of self that lies in the box.

        INPUT:

        - ``box`` - tuple of size three, size of the box

        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D
            sage: p = Pentamino3D(0)
            sage: L = list(p.translated_rotated(box=(5,8,2)))
            sage: len(L)
            360

        ::

            sage: p = Pentamino3D(6)
            sage: L = list(p.translated_rotated(box=(5,8,2)))
            sage: len(L)
            180
        """
        height = box[2]
        all_distinct_cano = set(self.rotated_canonical())
        for cano in all_distinct_cano:
            if cano.is_higher_than(height):
                continue
            for t in cano.translated(box=box):
                yield t

    def middle_of_neighbor_coords(self):
        r"""
        Return the list of middle of neighbor coords.

        This is use to draw cube in between two neighbor cubes.

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0),(0,0,1)])
            sage: list(p.middle_of_neighbor_coords())
            [(0.0, 0.0, 0.5)]

        ::

            sage: from sage.games.quantumino import Pentamino3D
            sage: L = sorted(Pentamino3D(0).middle_of_neighbor_coords())
            sage: for a in L: a
            (0.5, 0.0, 0.0)
            (1.0, 0.5, 0.0)
            (1.0, 1.0, 0.5)
            (1.0, 1.5, 0.0)
        """
        for (a,b,c), (p,q,r) in itertools.combinations(self, 2):
            if [0,0,1] == sorted(map(abs, (p-a, q-b, r-c))):
                yield ( (a+p)/2.0, (b+q)/2.0, (c+r)/2.0 )

    def show3d(self, size=0.75):
        r"""
        INPUT:

        - ``size`` - number (optional, default: ``0.75``), the size of the
          ``1 \times 1 \times 1`` cubes

        EXAMPLES::

            sage: from sage.games.quantumino import Polyomino3d
            sage: p = Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: p.show3d()
        """
        #bug trac #11272
        G = Graphics()
        for p in self:
            G += cube(color=self._color, size=size).translate(p)
        for m in self.middle_of_neighbor_coords():
            G += cube(color=self._color, size=size).translate(m)
        return G

#######################
# General puzzle solver
#######################
def solve_3d_puzzle(pieces, box=(5,8,2), verbose=False):
    r"""
    Solve the 3d puzzle, that is find the partition of the box into translated
    and rotated pieces where each pieces is used exactly once.

    INPUT:

    - ``pieces`` - iterable of Polyominos3d
    - ``box`` - tuple of size three (optional, default: ``(5,8,2)``),
      size of the box
    - ``verbose`` - bool

    OUTPUT:

    iterator of list of translated and rotated polyominos 3d that partition the box

    EXAMPLES::

        sage: from sage.games.quantumino import Pentamino3D, solve_3d_puzzle
        sage: L = map(Pentamino3D, range(17))
        sage: all_solutions = solve_3d_puzzle(L[:16], verbose=True)
        sage: K = all_solutions.next()
        Number of rows : 5664
        Number of distinct rows : 5664
        Row indices of the solution : [24, 696, 815, 1172, 1699, 1969, 5318, 5277, 4024, 3271, 4851, 2484, 4748, 4555, 2248, 2659]
        sage: for k in K: k
        3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
        3D Polyomino: [(0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 2, 1), (1, 2, 1)], Color: deeppink
        3D Polyomino: [(0, 2, 0), (0, 3, 0), (0, 4, 0), (1, 4, 0), (1, 4, 1)], Color: green
        3D Polyomino: [(0, 3, 1), (1, 3, 1), (2, 2, 0), (2, 2, 1), (2, 3, 1)], Color: green
        3D Polyomino: [(1, 3, 0), (2, 3, 0), (2, 4, 0), (2, 4, 1), (3, 4, 0)], Color: red
        3D Polyomino: [(1, 0, 1), (2, 0, 1), (2, 1, 0), (2, 1, 1), (3, 1, 1)], Color: red
        3D Polyomino: [(2, 0, 0), (3, 0, 0), (3, 0, 1), (4, 0, 1), (4, 1, 1)], Color: purple
        3D Polyomino: [(3, 1, 0), (3, 2, 0), (3, 2, 1), (4, 0, 0), (4, 1, 0)], Color: purple
        3D Polyomino: [(0, 4, 1), (0, 5, 0), (0, 5, 1), (0, 6, 1), (1, 5, 0)], Color: midnightblue
        3D Polyomino: [(3, 3, 0), (3, 3, 1), (4, 2, 0), (4, 2, 1), (4, 3, 1)], Color: yellow
        3D Polyomino: [(3, 4, 1), (3, 5, 1), (4, 3, 0), (4, 4, 0), (4, 4, 1)], Color: blue
        3D Polyomino: [(0, 7, 0), (0, 7, 1), (1, 7, 1), (2, 6, 1), (2, 7, 1)], Color: orange
        3D Polyomino: [(0, 6, 0), (1, 5, 1), (1, 6, 0), (1, 6, 1), (2, 5, 1)], Color: blue
        3D Polyomino: [(1, 7, 0), (2, 7, 0), (3, 6, 0), (3, 7, 0), (3, 7, 1)], Color: darkblue
        3D Polyomino: [(2, 5, 0), (2, 6, 0), (3, 5, 0), (4, 5, 0), (4, 5, 1)], Color: orange
        3D Polyomino: [(3, 6, 1), (4, 6, 0), (4, 6, 1), (4, 7, 0), (4, 7, 1)], Color: yellow
    """
    box_x, box_y, box_z = box
    number_of_pieces = len(pieces)

    # Bijection column indices <-> coordinates
    SPACE = list(itertools.product(range(box_x), range(box_y), range(box_z)))
    COORD_TO_INT = dict( (c,i+number_of_pieces) for i,c in enumerate(SPACE) )
    INT_TO_COORD = dict( (i+number_of_pieces,c) for i,c in enumerate(SPACE) )

    # Creation of the rows
    rows = []
    for i,p in enumerate(pieces):
        for q in p.translated_rotated(box):
            rows += [[i] + [COORD_TO_INT[coord] for coord in q]]
    if verbose:
        print "Number of rows : %s" % len(rows)
        print "Number of distinct rows : %s" % len(set(map(tuple,rows)))

    # Solve using dancing links
    from sage.combinat.matrices.dancing_links import dlx_solver
    x = dlx_solver(rows)

    while True:
        a = x.search()
        solution = x.get_solution()
        if verbose:
            print "Row indices of the solution : %s" % solution

        #Recover the pentos from the solution
        pentos = []
        for row_number in solution:
            row = rows[row_number]
            no = row[0]
            coords = [INT_TO_COORD[i] for i in row[1:]]
            p = Polyomino3d(coords, color=pieces[no].color())
            pentos.append(p)
        yield pentos

################################################
# Example:  The family games america: Quantumino
################################################
_list_pentaminos = []
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,2,0), (1,1,1)], color='deeppink'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (-1,0,0), (0,0,1)], color='deeppink'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,2,0), (0,0,1)], color='green'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (0,1,0), (0,2,0), (1,0,0), (1,0,1)], color='green'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,0,1), (1,-1,1)], color='red'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,0,1), (2,0,1)], color='red'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,2,0), (1,2,1)], color='orange'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (0,1,0), (0,2,0), (0,2,1)], color='orange'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1)], color='yellow'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,1,1), (0,0,1)], color='yellow'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (0,2,0), (1,1,1)], color='midnightblue'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,0,1), (1,2,0)], color='darkblue'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,1,1), (2,1,1)], color='blue'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,1,1), (1,2,1)], color='blue'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (2,1,0), (2,1,1)], color='purple'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (0,1,0), (1,1,0), (1,2,0), (1,2,1)], color='purple'))
_list_pentaminos.append(Polyomino3d([(0,0,0), (1,0,0), (1,1,0), (1,0,1), (1,-1,0)], color='gray'))

def Pentamino3D(i):
    r"""
    Return one of the seventeen 3D Pentamino pieces included in the game.

    INPUT:

    - ``i`` - integer, from 0 to 16.

    EXAMPLES::

        sage: from sage.games.quantumino import Pentamino3D
        sage: Pentamino3D(3)
        3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (1, 0, 0), (1, 0, 1)], Color: green
    """
    return _list_pentaminos[i]


##############################
# Class QuantuminoSolultion
##############################
class QuantuminoSolution(SageObject):
    r"""
    A solution of the Quantumino puzzle.

    INPUT:

    - ``pentos`` - list of 16 3D pentamino representing the solution
    - ``aside`` - the unused 3D pentamino

    EXAMPLES::

        sage: from sage.games.quantumino import Pentamino3D, QuantuminoSolution
        sage: p = Pentamino3D(0)
        sage: q = Pentamino3D(5)
        sage: r = Pentamino3D(11)
        sage: S = QuantuminoSolution([p,q], r)
        sage: S
        Quantumino solution without the following pentamino :
        3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 2, 0)], Color: darkblue

    ::

        sage: from sage.games.quantumino import quantumino_solver
        sage: quantumino_solver.get_solution(3)      # long time (1.5s)
        Quantumino solution without the following pentamino :
        3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (1, 0, 0), (1, 0, 1)], Color: green
    """
    def __init__(self, pentos, aside):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D, QuantuminoSolution
            sage: p = Pentamino3D(0)
            sage: q = Pentamino3D(5)
            sage: QuantuminoSolution([p], q)
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (2, 0, 1)], Color: red
        """
        self._pentos = pentos
        self._aside = aside

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D, QuantuminoSolution
            sage: p = Pentamino3D(0)
            sage: q = Pentamino3D(5)
            sage: QuantuminoSolution([p], q)
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (2, 0, 1)], Color: red
        """
        return "Quantumino solution without the following pentamino :\n%s" % self._aside

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import Pentamino3D, QuantuminoSolution
            sage: p = Pentamino3D(0)
            sage: q = Pentamino3D(5)
            sage: r = Pentamino3D(11)
            sage: S = QuantuminoSolution([p,q], r)
            sage: for a in S: a
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (2, 0, 1)], Color: red
        """
        return iter(self._pentos)

    def show3d(self):
        r"""
        Show the solution in 3D.

        EXAMPLES::

            sage: from sage.games.quantumino import quantumino_solver
            sage: s = quantumino_solver.get_solution(0)     # long time (1.5s)
            sage: s.show3d()      # long time (<1s)

        """
        G = Graphics()
        for p in self:
            G += p.show3d()
        aside_pento = self._aside.canonical() + (2.5,-4,0)
        G += aside_pento.show3d()
        G.show(aspect_ratio=1, frame=False)
        #return G


##############################
# Class QuantuminoSolver
##############################
class QuantuminoSolver(SageObject):
    r"""
    """
    def list_pentaminos(self):
        r"""
        Return the list of the 17 3-D pentaminos included in the game.

        OUTPUT:

        list

        EXAMPLES::

            sage: from sage.games.quantumino import quantumino_solver
            sage: L = quantumino_solver.list_pentaminos()
            sage: len(L)
            17
        """
        return _list_pentaminos

    def show_pentaminos(self, box=(5,8,2)):
        r"""
        Show the 17 3-D pentaminos included in the game and the `5 \times 8
        \times 2` box where 16 of them must fit.

        INPUT:

        - ``box`` - tuple of size three (optional, default: ``(5,8,2)``),
          size of the box

        OUTPUT:

        3D Graphic object

        EXAMPLES::

            sage: from sage.games.quantumino import quantumino_solver
            sage: quantumino_solver.show_pentaminos()    # long time (1s)
        """
        G = Graphics()
        for i,p in enumerate(self.list_pentaminos()):
            x = 3.5 * (i%4)
            y = 3.5 * (i/4)
            q = p + (x, y, 0)
            G += q.show3d()
            G += text3d(str(i), (x,y,2))
        G += cube(color='gray',opacity=0.5).scale(box).translate((17,6,0))
        a,b = G.bounding_box()
        a,b = map(vector, (a,b))
        G.frame_aspect_ratio(tuple(b-a))
        G.show(frame=False)
        #return G

    def get_solution(self, without, box=(5,8,2), verbose=False):
        r"""
        Return a solution without one of the block.

        INPUT:

        - ``without`` - integer, from 0 to 16.
        - ``box`` - tuple of size three (optional, default: ``(5,8,2)``),
          size of the box
        - ``verbose`` - bool (optional, default: ``False``)

        OUTPUT:

        QuantuminoSolution

        EXAMPLES::

            sage: from sage.games.quantumino import quantumino_solver
            sage: s = quantumino_solver.get_solution(8)
            sage: s
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 1, 0)], Color: yellow
            sage: for p in s: p
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            3D Polyomino: [(0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 2, 1), (1, 2, 1)], Color: deeppink
            3D Polyomino: [(0, 2, 0), (0, 3, 0), (0, 4, 0), (1, 4, 0), (1, 4, 1)], Color: green
            3D Polyomino: [(0, 3, 1), (1, 3, 1), (2, 2, 0), (2, 2, 1), (2, 3, 1)], Color: green
            3D Polyomino: [(1, 3, 0), (2, 3, 0), (2, 4, 0), (2, 4, 1), (3, 4, 0)], Color: red
            3D Polyomino: [(1, 0, 1), (2, 0, 0), (2, 0, 1), (2, 1, 0), (3, 0, 1)], Color: midnightblue
            3D Polyomino: [(0, 4, 1), (0, 5, 0), (0, 5, 1), (0, 6, 0), (1, 5, 0)], Color: red
            3D Polyomino: [(2, 1, 1), (3, 0, 0), (3, 1, 0), (3, 1, 1), (4, 0, 0)], Color: blue
            3D Polyomino: [(3, 2, 0), (4, 0, 1), (4, 1, 0), (4, 1, 1), (4, 2, 0)], Color: purple
            3D Polyomino: [(3, 2, 1), (3, 3, 0), (4, 2, 1), (4, 3, 0), (4, 3, 1)], Color: yellow
            3D Polyomino: [(3, 3, 1), (3, 4, 1), (4, 4, 0), (4, 4, 1), (4, 5, 0)], Color: blue
            3D Polyomino: [(0, 6, 1), (0, 7, 0), (0, 7, 1), (1, 5, 1), (1, 6, 1)], Color: purple
            3D Polyomino: [(1, 6, 0), (1, 7, 0), (1, 7, 1), (2, 7, 0), (3, 7, 0)], Color: darkblue
            3D Polyomino: [(2, 5, 0), (2, 6, 0), (3, 6, 0), (4, 6, 0), (4, 6, 1)], Color: orange
            3D Polyomino: [(2, 5, 1), (3, 5, 0), (3, 5, 1), (3, 6, 1), (4, 5, 1)], Color: gray
            3D Polyomino: [(2, 6, 1), (2, 7, 1), (3, 7, 1), (4, 7, 0), (4, 7, 1)], Color: orange

        Generalizations of the game inside different boxes::

            sage: quantumino_solver.get_solution(7, (4,4,5))       # long time (2s)
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
            sage: quantumino_solver.get_solution(7, (2,2,20))      # long time (1s)
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
            sage: quantumino_solver.get_solution(3, (2,2,20))      # long time (1s)
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (1, 0, 0), (1, 0, 1)], Color: green
        """
        return next(self.solutions_iterator(without, box, verbose))

    def solutions_iterator(self, without, box=(5,8,2), verbose=False):
        r"""
        Return an iterator over the solutions without the block.

        INPUT:

        - ``without`` - integer, from 0 to 16.
        - ``box`` - tuple of size three (optional, default: ``(5,8,2)``),
          size of the box
        - ``verbose`` - bool (optional, default: ``False``)

        OUTPUT:

        iterator of QuantuminoSolution

        EXAMPLES::

            sage: from sage.games.quantumino import quantumino_solver
            sage: it = quantumino_solver.solutions_iterator(0)
            sage: it.next()
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: it.next()
            Quantumino solution without the following pentamino :
            3D Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
        """
        if not  0 <= without < 17:
            raise ValueError, "without (=%s) must be between 0 and 16" % without
        all_pieces = self.list_pentaminos()
        pieces = all_pieces[:without] + all_pieces[without+1:]
        aside = all_pieces[without]
        for pentos in solve_3d_puzzle(pieces, box=box, verbose=verbose):
            yield QuantuminoSolution(pentos, aside)



quantumino_solver = QuantuminoSolver()


