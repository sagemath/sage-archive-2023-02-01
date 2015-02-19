"""
Coxeter matrices
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2015 Travis Scrimshaw <tscrim at ucdavis.edu>
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
from coxeter_type import CartanType
from sage.matrix.all import MatrixSpace
from sage.rings.all import ZZ

from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.rings.all import ZZ
from sage.combinat.root_system.coxeter_type import CoxeterType, CoxeterType_abstract
from sage.combinat.root_system.root_system import RootSystem
from sage.sets.family import Family

class CoxeterMatrix(Matrix_integer_dense, CoxeterType_abstract):
    """
    A Coxeter matrix.
    """
    def __init__(self, parent, data, coxeter_type, index_set):
        """
        Initialize ``self``.

        TESTS::

            sage: C = CoxeterMatrix(['A',1,1])
            sage: TestSuite(C).run(skip=["_test_category", "_test_change_ring"])
        """
        Matrix_integer_sparse.__init__(self, parent, data, False, False)
        self._coxeter_type = coxeter_type
        self._index_set = index_set
        self.set_immutable()

    ##########################################################################
    # Coxeter type methods

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',1,1])
            sage: C.index_set()
            (0, 1)
            sage: C = CoxeterMatrix(['E',6])
            sage: C.index_set()
            (1, 2, 3, 4, 5, 6)
        """
        return self._index_set

    def coxeter_type(self):
        """
        Return the Cartan type of ``self`` or ``self`` if unknown.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',4,1])
            sage: C.coxeter_type()
            ['A', 4, 1]

        If the Cartan type is unknown::

            sage: C = CoxeterMatrix([[2,-1,-2], [-1,2,-1], [-2,-1,2]])
            sage: C.coxeter_type()
            [ 2 -1 -2]
            [-1  2 -1]
            [-2 -1  2]
        """
        if self._coxeter_type is None:
            return self
        return self._coxeter_type

    def rank(self):
        r"""
        Return the rank of ``self``.

        EXAMPLES::

            sage: CoxeterMatrix(['C',3]).rank()
            3
            sage: CoxeterMatrix(["A2","B2","F4"]).rank()
            8
        """
        return self.ncols()

    def coxeter_matrix(self):
        r"""
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: CoxeterMatrix(['C',3]).coxeter_matrix()
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
        """
        return self

    def is_simply_laced(self):
        """
        Return if ``self`` is simply-laced.

        A Coxeter matrix is simply-laced if all non-diagonal entries are
        either 2 or 3.

        EXAMPLES::

            sage: cm = CoxeterMatrix([[2, -1, -1, -1], [-1, 2, -1, -1], [-1, -1, 2, -1], [-1, -1, -1, 2]])
            sage: cm.is_simply_laced()
            True
        """
        # We include 1 in this list to account for the diagonal
        L = [1,2,3]
        return all(x in L for row in self for x in row)

    def is_crystallographic(self):
        """
        Return if ``self`` is crystallographic.

        A Coxeter matrix is crystallographic if all non-diagonal entries
        are either 2, 4, or 6.

        EXAMPLES::

            sage: CoxeterMatrix(['F',4]).is_crystallographic()
            True
        """
        # We include 1 in this list to account for the diagonal
        L = [1,2,3,4,6]
        return all(x in L for row in self for x in row)

    def is_finite(self):
        """
        Return if ``self`` is a finite type or ``False`` if unknown.

        EXAMPLES::

            sage: M = CoxeterMatrix(['C',4])
            sage: M.is_finite()
            True
            sage: M = CoxeterMatrix(['D',4,1])
            sage: M.is_finite()
            False
            sage: M = CoxeterMatrix([[2, -4], [-3, 2]])
            sage: M.is_finite()
            False
        """
        if self._coxeter_type is not None:
            return self._coxeter_type.is_finite()
        return self.bilinear_form().is_positive_semidefinite()

    def is_affine(self):
        """
        Return if ``self`` is an affine type or ``False`` if unknown.

        EXAMPLES::

            sage: M = CoxeterMatrix(['C',4])
            sage: M.is_affine()
            False
            sage: M = CoxeterMatrix(['D',4,1])
            sage: M.is_affine()
            True
            sage: M = CoxeterMatrix([[2, -4], [-3, 2]])
            sage: M.is_affine()
            False
        """
        if self._coxeter_type is not None:
            return self._coxeter_type.is_affine()
        return self.bilinear_form().is_positive_semidefinite()

#####################################################################
## Other methods

def coxeter_matrix_as_function(t):
    """
    Returns the Coxeter matrix, as a function

    INPUT:

    - ``t`` -- a Cartan type

    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_matrix import coxeter_matrix_as_function
        sage: f = coxeter_matrix_as_function(['A',4])
        sage: matrix([[f(i,j) for j in range(1,5)] for i in range(1,5)])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
    """
    t = CartanType(t)
    m = t.coxeter_matrix()
    index_set = t.index_set()
    reverse = dict((index_set[i], i) for i in range(len(index_set)))
    return lambda i,j: m[reverse[i], reverse[j]]

def coxeter_matrix(t):
    """
    Returns the Coxeter matrix of type t.

    EXAMPLES::

        sage: coxeter_matrix(['A', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
        sage: coxeter_matrix(['B', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: coxeter_matrix(['C', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: coxeter_matrix(['D', 4])
        [1 3 2 2]
        [3 1 3 3]
        [2 3 1 2]
        [2 3 2 1]

    ::

        sage: coxeter_matrix(['E', 6])
        [1 2 3 2 2 2]
        [2 1 2 3 2 2]
        [3 2 1 3 2 2]
        [2 3 3 1 3 2]
        [2 2 2 3 1 3]
        [2 2 2 2 3 1]

    ::

        sage: coxeter_matrix(['F', 4])
        [1 3 2 2]
        [3 1 4 2]
        [2 4 1 3]
        [2 2 3 1]

    ::

        sage: coxeter_matrix(['G', 2])
        [1 6]
        [6 1]
    """
    # TODO: Deprecate this in favor of CoxeterMatrix
    return CoxeterMatrix(t)

