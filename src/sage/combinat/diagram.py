r"""
Combinatorial diagrams

A combinatorial diagram is a collection of cells `(i,j)` indexed by pairs of
natural numbers. The positions are indexed by rows and columns. For example,
a Ferrer's diagram is a diagram obtained from a partition 
`\lambda = (\lambda_0, \lambda_1, \ldots, \lambda_\ell)` where the cells are
in rows `i` for `0 \leq i \leq \ell` and the cells in row `i` consist of
`(i,j)` for `0 \leq j < \lambda_i`. In English notation, the indices are read
from top left to bottom right as in a matrix.

Indexing conventions are the same as
:class:`~sage.combinat.partition.Partition`.

EXAMPLES:


Diagrams can be created by:
    - explictly passing an iterable of all the cells
    - providing a :class:`~sage.combinat.partition.Partition`,
      in which case the diagram is the corresponding (English notation) Ferrers
      diagram (see also :meth:`~sage.combination.partition.Partition.ferrers_diagram`)
    - providing a :class:`~sage.combinat.permutation.Permutation`, in which
      case the diagram is the corresponding Rothe diagram.
    - providing a list of *death rays* which are cells `(i_0,j_0)` not present in the
      diagram and have the property that there are no cells in the diagram that
      have the form `(i,j_0)` with `i > i_0` and no cells in the diagram that
      have the form `(i_0, j)` with `j > j_0`. The death rays kill all of the
      cells to right and below of them.


Passing a list of all cells::

    sage: from sage.combinat.diagram import Diagram
    sage: cells = [(0,0), (0,1), (1,0), (1,1), (4,4), (4,5), (4,6), (5,4), (7, 6)]
    sage: D = Diagram(cells); D
    [(0,0), (0,1), (1,0), (1,1), (4,4), (4,5), (4,6), (5,4), (7, 6)]
    
    
AUTHORS:

- Trevor K. Karn (2022-08-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Trevor K. Karn <karnx018 (at) umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.sets_cat import Sets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent

class Diagram(ClonableArray):
    def __init__(self, cells):
        self._cells = cells

class Diagrams(UniqueRepresentation, Parent):

    def __init__(self):

        Parent.__init__(self)

    def _element_constructor_(self, cells)

        return self.element_class(self, cells)

    Element = Diagram