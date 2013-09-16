r"""
Strong and weak tableaux

There are two types of `k`-tableaux: strong `k`-tableaux and weak `k`-tableaux.
Standard weak `k`-tableaux correspond to saturated chains in the weak order,
whereas standard strong `k`-tableaux correspond to saturated chains in the strong Bruhat order.
For semistandard tableaux, the notion of weak and strong horizontal strip is necessary.
More information can be found in [LLMS2006]_ .

REFERENCES:

.. [LLMS2006] T. Lam, L. Lapointe, J. Morse, M. Shimozono,
   Affine insertion and Pieri rules for the affine Grassmannian,
   Memoirs of the AMS, 208 (2010), no. 977, :arxiv:`math.CO/0609110`

.. [LLMSSZ2013] T. Lam, L. Lapointe, J. Morse, A. Schilling, M. Shimozono, M. Zabrocki,
   `k`-Schur functions and affine Schubert calculus,
   preprint :arXiv:`1301.3569`

Authors:

- Anne Schilling and Mike Zabrocki (2013): initial version
- Avi Dalal and Nate Gallup (2013): implementation of `k`-charge
"""
#*****************************************************************************
#       Copyright (C) 2013 Anne Schilling <anne at math.ucdavis.edu>
#                          Mike Zabrocki  <zabrocki at mathstat.yorku.ca>
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
#****************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.list_clone import ClonableList
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.combinat.skew_tableau import SkewTableau, SemistandardSkewTableaux
from sage.combinat.partition import Partition, Partitions
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.combinat.core import Core
from sage.rings.all import ZZ
from sage.misc.misc import uniq

def WeakTableau(t, k, inner_shape = [], representation = "core"):
    r"""
    This is the dispatcher method for the element class of weak `k`-tableaux.

    Standard weak `k`-tableaux correspond to saturated chains in the weak order.
    There are three formulations of weak tableaux, one in terms of cores, one in terms
    of `k`-bounded partitions, and one in terms of factorizations of affine Grassmannian
    elements. For semistandard weak `k`-tableaux, all letters of the same value have to
    satisfy the conditions of a horizontal strip. In the affine Grassmannian formulation this
    means that all factors are cyclically decreasing elements. For more information, see
    for example [LLMSSZ2013]_.

    INPUT:

    - ``t`` -- a weak `k`-tableau in the specified representation:

      - for the 'core' representation ``t`` is a list of lists where each subtableaux
        should have a `k+1`-core shape; ``None`` is allowed as an entry for skew weak
        `k`-tableaux
      - for the 'bounded' representation ``t`` is a list of lists where each subtableaux
        should have a `k`-bounded shape; ``None`` is allowed as an entry for skew weak
        `k`-tableaux
      - for the 'factorized_permutation' representation ``t`` is either a list of
        cyclically decreasing Weyl group elements or a list of reduced words of cyclically
        decreasing Weyl group elements; to indicate a skew tableau in this representation,
        ``inner_shape`` should be the inner shape as a `(k+1)`-core

    - ``k`` -- positive integer

    - ``inner_shape`` -- this entry is only relevant for the 'factorized_permutation'
      representation and specifies the inner shape in case the tableau is skew
      (default: ``[]``)

    - ``representation`` -- 'core', 'bounded', or 'factorized_permutation'
      (default: 'core')

    EXAMPLES:

    Here is an example of a weak 3-tableau in core representation::

        sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
        sage: t.shape()
        [5, 2, 1]
        sage: t.weight()
        (2, 2, 2)
        sage: type(t)
        <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>

    And now we give a skew weak 3-tableau in core representation::

        sage: ts = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
        sage: ts.shape()
        ([5, 2, 1], [1, 1])
        sage: ts.weight()
        (2, 2)
        sage: type(ts)
        <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>

    Next we create the analogue of the first example in bounded representation::

        sage: tb = WeakTableau([[1,1,2],[2,3],[3]], 3, representation="bounded")
        sage: tb.shape()
        [3, 2, 1]
        sage: tb.weight()
        (2, 2, 2)
        sage: type(tb)
        <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
        sage: tb.to_core_tableau()
        [[1, 1, 2, 2, 3], [2, 3], [3]]
        sage: t == tb.to_core_tableau()
        True

    And the analogue of the skew example in bounded representation::

        sage: tbs = WeakTableau([[None, 1, 2], [None, 2], [1]], 3, representation = "bounded")
        sage: tbs.shape()
        ([3, 2, 1], [1, 1])
        sage: tbs.weight()
        (2, 2)
        sage: tbs.to_core_tableau()
        [[None, 1, 1, 2, 2], [None, 2], [1]]
        sage: ts.to_bounded_tableau() == tbs
        True

    Finally we do the same examples for the factorized permutation representation::

        sage: tf = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = "factorized_permutation")
        sage: tf.shape()
        [5, 2, 1]
        sage: tf.weight()
        (2, 2, 2)
        sage: type(tf)
        <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
        sage: tf.to_core_tableau() == t
        True

        sage: tfs = WeakTableau([[0,3],[2,1]], 3, inner_shape = [1,1], representation = 'factorized_permutation')
        sage: tfs.shape()
        ([5, 2, 1], [1, 1])
        sage: tfs.weight()
        (2, 2)
        sage: type(tfs)
        <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
        sage: tfs.to_core_tableau()
        [[None, 1, 1, 2, 2], [None, 2], [1]]

    Another way to pass from one representation to another is as follows::

        sage: ts
        [[None, 1, 1, 2, 2], [None, 2], [1]]
        sage: ts.parent()._representation
        'core'
        sage: ts.representation('bounded')
        [[None, 1, 2], [None, 2], [1]]

    To test whether a given semistandard tableau is a weak `k`-tableau in the bounded representation,
    one can ask::

        sage: t = Tableau([[1,1,2],[2,3],[3]])
        sage: t.is_k_tableau(3)
        True
        sage: t = SkewTableau([[None, 1, 2], [None, 2], [1]])
        sage: t.is_k_tableau(3)
        True
        sage: t = SkewTableau([[None, 1, 1], [None, 2], [2]])
        sage: t.is_k_tableau(3)
        False

    TESTS::

        sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = "bla")
        Traceback (most recent call last):
        ...
        NotImplementedError: The representation option needs to be 'core', 'bounded', or 'factorized_permuation'
    """
    if representation == "core":
        return WeakTableau_core(t, k)
    elif representation == "bounded":
        return WeakTableau_bounded(t, k)
    elif representation == "factorized_permutation":
        return WeakTableau_factorized_permutation(t, k, inner_shape = inner_shape)
    else:
        raise NotImplementedError, "The representation option needs to be 'core', 'bounded', or 'factorized_permuation'"

def WeakTableaux(k, shape , weight, representation = "core"):
    r"""
    This is the dispatcher method for the parent class of weak `k`-tableaux.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- shape of the weak `k`-tableaux; for the 'core' and
      'factorized_permutation' representation, the shape is inputted as a `(k+1)`-core;
      for the 'bounded' representation, the shape is inputted as a `k`-bounded partition;
      for skew tableaux, the shape is inputted as a tuple of the outer and inner shape
    - ``weight`` -- the weight of the weak `k`-tableaux as a list or tuple
    - ``representation`` -- 'core', 'bounded', or 'factorized_permutation' (default: 'core')

    EXAMPLES::

        sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1])
        sage: T.list()
        [[[1, 3, 4, 5, 6], [2, 6], [4]],
        [[1, 2, 4, 5, 6], [3, 6], [4]],
        [[1, 2, 3, 4, 6], [4, 6], [5]],
        [[1, 2, 3, 4, 5], [4, 5], [6]]]
        sage: T.cardinality()
        4

        sage: T = WeakTableaux(3, [[5,2,1], [2]], [1,1,1,1])
        sage: T.list()
        [[[None, None, 2, 3, 4], [1, 4], [2]],
        [[None, None, 1, 2, 4], [2, 4], [3]],
        [[None, None, 1, 2, 3], [2, 3], [4]]]

        sage: T = WeakTableaux(3, [3,2,1], [1,1,1,1,1,1], representation = 'bounded')
        sage: T.list()
        [[[1, 3, 5], [2, 6], [4]],
        [[1, 2, 5], [3, 6], [4]],
        [[1, 2, 3], [4, 6], [5]],
        [[1, 2, 3], [4, 5], [6]]]

        sage: T = WeakTableaux(3, [[3,2,1], [2]], [1,1,1,1], representation = 'bounded')
        sage: T.list()
        [[[None, None, 3], [1, 4], [2]],
        [[None, None, 1], [2, 4], [3]],
        [[None, None, 1], [2, 3], [4]]]

        sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1], representation = 'factorized_permutation')
        sage: T.list()
        [[s0, s3, s2, s1, s3, s0],
        [s0, s3, s2, s3, s1, s0],
        [s0, s2, s3, s2, s1, s0],
        [s2, s0, s3, s2, s1, s0]]

        sage: T = WeakTableaux(3, [[5,2,1], [2]], [1,1,1,1], representation = 'factorized_permutation')
        sage: T.list()
        [[s0, s3, s2, s3], [s0, s2, s3, s2], [s2, s0, s3, s2]]
    """
    if representation == "core":
        return WeakTableaux_core(k, shape, weight)
    elif representation == "bounded":
        return WeakTableaux_bounded(k, shape, weight)
    elif representation == "factorized_permutation":
        return WeakTableaux_factorized_permutation(k, shape, weight)
    else:
        raise NotImplementedError, "The representation option needs to be 'core', 'bounded', or 'factorized_permuation'"

#Abstract class for the elements of weak tableau
class WeakTableau_abstract(ClonableList):
    r"""
    Abstract class for the various element classes of WeakTableau.
    """
    __metaclass__ = ClasscallMetaclass

    def shape(self):
        r"""
        Return the shape of ``self``.

        When the tableau is straight, the outer shape is returned.
        When the tableau is skew, the tuple of the outer and inner shape is returned.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.shape()
            [5, 2, 1]
            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.shape()
            ([5, 2, 1], [2])

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.shape()
            [3, 2, 1]
            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.shape()
            ([3, 2, 1], [2])

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.shape()
            [5, 2, 1]
            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.shape()
            ([5, 2, 1], [2])
        """
        return self.parent().shape()

    def weight(self):
        r"""
        Return the weight of ``self``.

        The weight is a tuple whose `i`-th entry is the number of labels `i` in the
        bounded representation of ``self``.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.weight()
            (2, 2, 2)
            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.weight()
            (1, 1, 1, 1)
            sage: t = WeakTableau([[None,2,3],[3]],2)
            sage: t.weight()
            (0, 1, 1)

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.weight()
            (3, 2, 1)
            sage: t = WeakTableau([[1,1,2],[2,3],[3]], 3, representation = 'bounded')
            sage: t.weight()
            (2, 2, 2)
            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.weight()
            (1, 1, 1, 1)

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.weight()
            (3, 2, 1)
            sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = 'factorized_permutation')
            sage: t.weight()
            (2, 2, 2)
            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.weight()
            (2, 2)
        """
        return self.parent()._weight

    def size_of_shape(self):
        r"""
        Return the size of the shape of ``self``.

        In the bounded representation, the size of the shape is the number of boxes in the
        outer shape minus the number of boxes in the inner shape. For the core and
        factorized permutation representation, the size is the length of the outer shape
        minus the length of the inner shape.

        .. SEEALSO:: :meth:`sage.combinat.core.Core.length`

        EXAMPLES::

            sage: t = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            sage: t.shape()
            ([5, 2, 1], [1, 1])
            sage: t.size_of_shape()
            4
            sage: t = WeakTableau([[1,1,2],[2,3],[3]], 3, representation="bounded")
            sage: t.shape()
            [3, 2, 1]
            sage: t.size_of_shape()
            6
        """
        return self.parent().size_of_shape()

    def intermediate_shapes(self):
        r"""
        Return the intermediate shapes of ``self``.

        A (skew) tableau with letters `1,2,\ldots,\ell` can be viewed as a sequence of shapes,
        where the `i`-th shape is given by the shape of the subtableau on letters `1,2,\ldots,i`.
        The output is the list of these shapes.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.intermediate_shapes()
            [[], [2], [4, 1], [5, 2, 1]]

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.intermediate_shapes()
            [[2], [2, 1], [3, 1, 1], [4, 1, 1], [5, 2, 1]]

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.intermediate_shapes()
            [[], [3], [3, 2], [3, 2, 1]]

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.intermediate_shapes()
            [[2], [3], [3, 1], [3, 1, 1], [3, 2, 1]]

            sage: t = WeakTableau([[0],[3],[2],[3]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.intermediate_shapes()
            [[2], [2, 1], [3, 1, 1], [4, 1, 1], [5, 2, 1]]
        """
        if self.parent()._representation in ['core', 'bounded']:
            return intermediate_shapes(self)
        else:
            return intermediate_shapes(self.to_core_tableau())

    def pp(self):
        r"""
        Return a pretty print string of the tableau.

        EXAMPLES::

            sage: t = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            sage: t.pp()
            .  1  1  2  2
            .  2
            1
            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.pp()
            [s2*s0, s3*s2]
        """
        if self.parent()._representation in ['core', 'bounded']:
            print self._repr_diagram()
        else:
            print self

    def _latex_(self):
        r"""
        Return a latex method for the tableau.

        EXAMPLES::

            sage: t = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            sage: latex(t)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}c}\cline{1-5}
            \lr{}&\lr{1}&\lr{1}&\lr{2}&\lr{2}\\\cline{1-5}
            \lr{}&\lr{2}\\\cline{1-2}
            \lr{1}\\\cline{1-1}
            \end{array}$}
            }

            sage: t = WeakTableau([[0,3],[2,1]], 3, inner_shape = [1,1], representation = 'factorized_permutation')
            sage: latex(t)
            [s_{0}s_{3},s_{2}s_{1}]
        """
        def chi(x):
            if x is None:
                return ""
            if x in ZZ:
                return x
            return "%s"%x
        if self.parent()._representation in ['core', 'bounded']:
            t = [[chi(x) for x in row] for row in self._list]
            from output import tex_from_array
            return tex_from_array(t)
        else:
            return "["+"".join(self[i]._latex_()+',' for i in range(len(self)-1))+self[len(self)-1]._latex_()+"]"

    def representation(self, representation = 'core'):
        r"""
        Return the analogue of ``self`` in the specified representation.

        INPUT:

        - ``representation`` -- 'core', 'bounded', or 'factorized_permutation' (default: 'core')

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.parent()._representation
            'core'
            sage: t.representation('bounded')
            [[1, 1, 2, 4], [2, 3, 5], [3, 4], [5, 6], [6], [7]]
            sage: t.representation('factorized_permutation')
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]

            sage: tb = WeakTableau([[1, 1, 2, 4], [2, 3, 5], [3, 4], [5, 6], [6], [7]], 4, representation = 'bounded')
            sage: tb.parent()._representation
            'bounded'
            sage: tb.representation('core') == t
            True
            sage: tb.representation('factorized_permutation')
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]

            sage: tp = WeakTableau([[0],[3,1],[2,1],[0,4],[3,0],[4,2],[1,0]], 4, representation = 'factorized_permutation')
            sage: tp.parent()._representation
            'factorized_permutation'
            sage: tp.representation('core') == t
            True
            sage: tp.representation('bounded') == tb
            True
        """
        t = self
        if self.parent()._representation in ['bounded', 'factorized_permutation']:
            t = t.to_core_tableau()
        if representation == 'core':
            return t
        elif representation == 'bounded':
            return t.to_bounded_tableau()
        elif representation == 'factorized_permutation':
            return t.to_factorized_permutation_tableau()
        else:
            raise ValueError("The representation must be one of 'core', 'bounded', or 'factorized_permutation'")

#Abstract class for the parents of weak tableaux
class WeakTableaux_abstract(UniqueRepresentation, Parent):
    r"""
    Abstract class for the various parent classes of WeakTableaux.
    """
    def shape(self):
        r"""
        Return the shape of the tableaux of ``self``.

        When ``self`` is the class of straight tableaux, the outer shape is returned.
        When ``self`` is the class of skew tableaux, the tuple of the outer and inner
        shape is returned.

        Note that in the 'core' and 'factorized_permutation' representation, the shapes
        are `(k+1)`-cores.  In the 'bounded' representation, the shapes are `k`-bounded
        partitions.

        If the user wants to access the skew shape (even if the inner shape is empty),
        please use ``self._shape``.

        EXAMPLES::

            sage: T = WeakTableaux(3, [5,2,2], [2,2,2,1])
            sage: T.shape()
            [5, 2, 2]
            sage: T._shape
            ([5, 2, 2], [])
            sage: T = WeakTableaux(3, [[5,2,2], [1]], [2,1,2,1])
            sage: T.shape()
            ([5, 2, 2], [1])

            sage: T = WeakTableaux(3, [3,2,2], [2,2,2,1], representation = 'bounded')
            sage: T.shape()
            [3, 2, 2]
            sage: T._shape
            ([3, 2, 2], [])
            sage: T = WeakTableaux(3, [[3,2,2], [1]], [2,1,2,1], representation = 'bounded')
            sage: T.shape()
            ([3, 2, 2], [1])

            sage: T = WeakTableaux(3, [4,1], [2,2], representation = 'factorized_permutation')
            sage: T.shape()
            [4, 1]
            sage: T._shape
            ([4, 1], [])
            sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], representation = 'factorized_permutation')
            sage: T.shape()
            ([6, 2, 1], [2])
        """
        if self._skew:
            return (self._outer_shape, self._inner_shape)
        return self._outer_shape

    def size_of_shape(self):
        r"""
        Return the size of the shape.

        In the bounded representation, the size of the shape is the number of boxes in the
        outer shape minus the number of boxes in the inner shape. For the core and
        factorized permutation representation, the size is the length of the outer shape
        minus the length of the inner shape.

        EXAMPLES::

            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1])
            sage: T.size_of_shape()
            6
            sage: T = WeakTableaux(3, [3,2,1], [1,1,1,1,1,1], representation = 'bounded')
            sage: T.size_of_shape()
            6
            sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], 'factorized_permutation')
            sage: T.size_of_shape()
            5
        """
        if self._representation == 'bounded':
            return self._outer_shape.size() - self._inner_shape.size()
        else:
            return self._outer_shape.length() - self._inner_shape.length()

    def representation(self, representation = 'core'):
        r"""
        Return the analogue of ``self`` in the specified representation.

        INPUT:

        - ``representation`` -- 'core', 'bounded', or 'factorized_permutation' (default: 'core')

        EXAMPLES::

            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1])
            sage: T._representation
            'core'
            sage: T.representation('bounded')
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('factorized_permutation')
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)

            sage: T = WeakTableaux(3, [3,2,1], [1,1,1,1,1,1], representation = 'bounded')
            sage: T._representation
            'bounded'
            sage: T.representation('core')
            Core weak 3-Tableaux of (skew) core shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('bounded')
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('bounded') == T
            True
            sage: T.representation('factorized_permutation')
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('factorized_permutation') == T
            False

            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1], representation = 'factorized_permutation')
            sage: T._representation
            'factorized_permutation'
            sage: T.representation('core')
            Core weak 3-Tableaux of (skew) core shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('bounded')
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('factorized_permutation')
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
        """
        outer_shape = self._outer_shape
        inner_shape = self._inner_shape
        weight = self._weight
        if (self._representation in ['core', 'factorized_permutation']) and representation == 'bounded':
            outer_shape = outer_shape.to_bounded_partition()
            inner_shape = inner_shape.to_bounded_partition()
        if self._representation == 'bounded' and (representation in ['core', 'factorized_permutation']):
            outer_shape = outer_shape.to_core(self.k)
            inner_shape = inner_shape.to_core(self.k)
        return WeakTableaux(self.k, [outer_shape, inner_shape], weight, representation = representation)


#Weak Tableaux in terms of cores
class WeakTableau_core(WeakTableau_abstract):
    r"""
    A (skew) weak `k`-tableau represented in terms of `(k+1)`-cores.
    """
    @staticmethod
    def __classcall_private__(cls, t, k):
        r"""
        Implements the shortcut ``WeakTableau_core(t, k)`` to ``WeakTableaux_core(k, shape , weight)(t)``
        where ``shape`` is the shape of the tableau and ``weight`` is its weight.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_core
            sage: t = WeakTableau_core([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            False

            sage: t = WeakTableau_core([[None, None, 1, 1, 2], [1, 2], [2]],3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            True
        """
        if isinstance(t, cls):
            return t
        tab = SkewTableau(list(t))
        outer = Core(tab.outer_shape(),k+1)
        inner = Core(tab.inner_shape(),k+1)
        weight = WeakTableau_bounded.from_core_tableau(t,k).weight()
        return WeakTableaux_core(k, [outer, inner], weight)(t)

    def __init__(self, parent, t):
        r"""
        Initialization of weak `k`-tableau ``t`` in core representation.

        INPUT:

        - ``t`` -- weak tableau in core representation; the input is supposed to be a list
          of lists specifying the rows of the tableau;
          ``None`` is allowed as an entry for skew weak `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_core, WeakTableaux_core
            sage: T = WeakTableaux_core(3,[5,2,1],[2,2,2])
            sage: t = T([[1, 1, 2, 2, 3], [2, 3], [3]]); t
            [[1, 1, 2, 2, 3], [2, 3], [3]]
            sage: c = WeakTableau_core([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: T = WeakTableaux_core(3,[5,2,1],[2,2,2])
            sage: t = T([[1, 1, 2, 2, 3], [2, 3], [3]]); t
            [[1, 1, 2, 2, 3], [2, 3], [3]]
            sage: c == t
            True
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: t.parent()
            Core weak 3-Tableaux of (skew) core shape [5, 2, 1] and weight (2, 2, 2)
            sage: TestSuite(t).run()

            sage: t = WeakTableau_core([[None, None, 1, 1, 2], [1, 2], [2]],3);  t
            [[None, None, 1, 1, 2], [1, 2], [2]]
            sage: t.weight()
            (2, 2)
            sage: t.shape()
            ([5, 2, 1], [2])
            sage: TestSuite(t).run()
        """
        self.k = parent.k
        self._list = [r for r in t]
        ClonableList.__init__(self, parent, t)

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: print t._repr_diagram()
            .  .  2  3  4
            1  4
            2
        """
        t = SkewTableau(list(self))
        return t._repr_diagram()

    def shape_core(self):
        r"""
        Return the shape of ``self`` as a `(k+1)`-core.

        When the tableau is straight, the outer shape is returned as a core.  When the
        tableau is skew, the tuple of the outer and inner shape is returned as cores.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.shape_core()
            [5, 2, 1]

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.shape_core()
            ([5, 2, 1], [2])
        """
        return self.shape()

    def shape_bounded(self):
        r"""
        Return the shape of ``self`` as a `k`-bounded partition.

        When the tableau is straight, the outer shape is returned as a `k`-bounded
        partition.  When the tableau is skew, the tuple of the outer and inner shape is
        returned as `k`-bounded partitions.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.shape_bounded()
            [3, 2, 1]

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.shape_bounded()
            ([3, 2, 1], [2])
        """
        if self.parent()._skew:
            return tuple([r.to_bounded_partition() for r in self.shape_core()])
        return self.shape_core().to_bounded_partition()

    def check(self):
        r"""
        Check that ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2], [2]], 2)
            sage: t.check()
            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.check()

        TESTS::

            sage: T = WeakTableaux(2, [3,1], [1,1,1,1])
            sage: t = T([[1,2,3],[3]])
            Traceback (most recent call last):
            ...
            ValueError: The weight of the parent does not agree with the weight of the tableau!

            sage: t = WeakTableau([[1, 2, 2], [1]], 2)
            Traceback (most recent call last):
            ...
            ValueError: The tableau is not semistandard!
        """
        if not self.parent()._weight == WeakTableau_bounded.from_core_tableau(self._list,self.k).weight():
            raise ValueError("The weight of the parent does not agree with the weight of the tableau!")
        t = SkewTableau(list(self))
        if not t in SemistandardSkewTableaux():
            raise ValueError("The tableau is not semistandard!")
        outer = Core(t.outer_shape(),self.k+1)
        inner = Core(t.inner_shape(),self.k+1)
        if self.parent()._outer_shape != outer:
            raise ValueError("The outer shape of the parent does not agree with the outer shape of the tableau!")
        if self.parent()._inner_shape != inner:
            raise ValueError("The inner shape of the parent does not agree with the inner shape of the tableau!")
        self.to_bounded_tableau().check()

    def to_bounded_tableau(self):
        r"""
        Return the bounded representation of the weak `k`-tableau ``self``.

        Each restricted sutableaux of the output is a `k`-bounded partition.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: c = t.to_bounded_tableau(); c
            [[1, 1, 2], [2, 3], [3]]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.to_bounded_tableau()
            [[None, None, 3], [1, 4], [2]]
            sage: t.to_bounded_tableau().to_core_tableau() == t
            True
        """
        shapes = [ Core(p,self.k+1).to_bounded_partition() for p in self.intermediate_shapes() ]
        if self.parent()._skew:
            l = [[None]*i for i in shapes[0]]
        else:
            l = []
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return WeakTableau_bounded(l, self.k)

    def to_factorized_permutation_tableau(self):
        r"""
        Return the factorized permutation representation of the weak `k`-tableau ``self``.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: c = t.to_factorized_permutation_tableau(); c
            [s2*s0, s3*s2, s1*s0]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
            sage: c.to_core_tableau() == t
            True

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: c = t.to_factorized_permutation_tableau(); c
            [s0, s3, s2, s3]
            sage: c._inner_shape
            [2]
            sage: c.to_core_tableau() == t
            True

        TESTS::

            sage: t = WeakTableau([], 4)
            sage: c = t.to_factorized_permutation_tableau(); c
            [1]
            sage: c._inner_shape
            []
            sage: c.to_core_tableau() == t
            True
        """
        shapes = [ Core(p,self.k+1).to_grassmannian() for p in self.intermediate_shapes() ]
        perms = [ shapes[i]*(shapes[i-1].inverse()) for i in range(len(shapes)-1,0,-1)]
        return WeakTableau_factorized_permutation(perms, self.k, inner_shape = self.parent()._inner_shape)

    def residues_of_entries(self, v):
        r"""
        Return a list of residues of cells of weak `k`-tableau ``self`` labeled by ``v``.

        INPUT:

        - ``v`` -- a label of a cell in ``self``

        OUTPUT:

        - a list of residues

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.residues_of_entries(1)
            [0, 1]

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.residues_of_entries(1)
            [2, 3]
        """
        return uniq([(j - i)%(self.k+1) for i in range(len(self)) for j in range(len(self[i])) if self[i][j] == v])

    def dictionary_of_coordinates_at_residues(self, v):
        r"""
        Return a dictionary assigning to all residues of ``self`` with label ``v`` a list
        of cells with the given residue.

        INPUT:

        - ``v`` -- a label of a cell in ``self``

        OUTPUT:

        - dictionary assigning coordinates in ``self`` to residues

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.dictionary_of_coordinates_at_residues(3)
            {0: [(0, 4), (1, 1)], 2: [(2, 0)]}

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.dictionary_of_coordinates_at_residues(1)
            {2: [(0, 2)], 3: [(0, 3), (1, 0)]}

            sage: t = WeakTableau([], 3)
            sage: t.dictionary_of_coordinates_at_residues(1)
            {}
        """
        d = {}
        for r in self.residues_of_entries(v):
            d[r] = []
            for i in range(len(self)):
                for j in range(len(self[i])):
                    if self[i][j]==v and (j - i)%(self.k+1) == r:
                        d[r]+=[(i,j)]
        return d

    def list_of_standard_cells(self):
        r"""
        Return a list of lists of the coordinates of the standard cells of ``self``.

        INPUT:

        - ``self`` -- a weak `k`-tableau in core representation with partition weight

        OUTPUT:

        - a list of lists of coordinates

        .. WARNING::

            This method currently only works for straight weak tableaux with partition
            weight.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.list_of_standard_cells()
            [[(0, 1), (1, 0), (2, 0)], [(0, 0), (0, 2), (1, 1)]]
            sage: t = WeakTableau([[1, 1, 1, 2], [2, 2, 3]], 5)
            sage: t.list_of_standard_cells()
            [[(0, 2), (1, 1), (1, 2)], [(0, 1), (1, 0)], [(0, 0), (0, 3)]]
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.list_of_standard_cells()
            [[(0, 1), (1, 0), (2, 0), (0, 5), (3, 0), (4, 0), (5, 0)], [(0, 0), (0, 2), (1, 1), (2, 1), (1, 2), (3, 1)]]

        TESTS::

            sage: t = WeakTableau([],3)
            sage: t.list_of_standard_cells()
            []

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.list_of_standard_cells()
            Traceback (most recent call last):
            ...
            ValueError: This method only works for straight tableaux!

            sage: t = WeakTableau([[1,2],[2]], 3)
            sage: t.list_of_standard_cells()
            Traceback (most recent call last):
            ...
            ValueError: This method only works for weak tableaux with partition weight!
        """
        if self.parent()._skew:
            raise ValueError("This method only works for straight tableaux!")
        if self.weight() not in Partitions(sum(self.weight())):
            raise ValueError("This method only works for weak tableaux with partition weight!")
        if self._list == []:
            return []
        mu = Partition(self.weight()).conjugate()
        already_used = []
        out = []
        for i in range(self[0].count(1)):
            standard_cells = [(0,self[0].count(1) - i - 1)]
            r = self[0].count(1) - i - 1
            for v in range(1,mu[i]):
                D = self.dictionary_of_coordinates_at_residues(v+1)
                new_D = {a:b for (a,b) in D.iteritems() if all(x not in already_used for x in b)}
                r = (r - min([self.k+1 - (x-r)%(self.k+1) for x in new_D.keys()]))%(self.k+1)
                standard_cells.append(new_D[r][-1])
                already_used += new_D[r]
            out.append(standard_cells)
        return out

    def k_charge(self, algorithm = "I"):
        r"""
        Return the `k`-charge of ``self``.

        INPUT:

        - ``algorithm`` -- (default: "I") if "I", computes `k`-charge using the `I`
          algorithm, otherwise uses the `J`-algorithm

        OUTPUT:

        - a nonnegative integer

        For the definition of `k`-charge and the various algorithms to compute it see
        Section 3.3 of [LLMSSZ2013]_.

        .. SEEALSO:: :meth:`k_charge_I` and :meth:`k_charge_J`

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.k_charge()
            2
            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: t.k_charge()
            8
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.k_charge()
            12

        TESTS::

            sage: T = WeakTableaux(4, [13,9,5,3,2,1,1], [4,3,3,2,2,1,1,1])
            sage: T.cardinality()
            6
            sage: all(t.k_charge_I() == t.k_charge_J() for t in T)
            True
        """
        if algorithm == "I":
            return self.k_charge_I()
        return self.k_charge_J()

    def k_charge_I(self):
        r"""
        Return the `k`-charge of ``self`` using the `I`-algorithm.

        For the definition of `k`-charge and the `I`-algorithm see Section 3.3 of [LLMSSZ2013]_.

        OUTPUT:

        - a nonnegative integer

        .. SEEALSO:: :meth:`k_charge` and :meth:`k_charge_J`

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.k_charge_I()
            2
            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: t.k_charge_I()
            8
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.k_charge_I()
            12

        TESTS::

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.k_charge_I()
            Traceback (most recent call last):
            ...
            ValueError: k-charge is not defined for skew weak tableaux
        """
        if self.parent()._skew:
            raise ValueError("k-charge is not defined for skew weak tableaux")
        stt = self.list_of_standard_cells()
        kch = 0
        for sw in stt:
            Ii = 0
            for r in range(len(sw)-1):
                if sw[r][1] < sw[r+1][1]:
                    Ii += 1 + abs(self.parent().diag(sw[r+1],sw[r]))
                else:
                    Ii += - abs(self.parent().diag(sw[r],sw[r+1]))
                kch += Ii
        return kch

    def k_charge_J(self):
        r"""
        Return the `k`-charge of ``self`` using the `J`-algorithm.

        For the definition of `k`-charge and the `J`-algorithm see Section 3.3 of [LLMSSZ2013]_.

        OUTPUT:

        - a nonnegative integer

        .. SEEALSO:: :meth:`k_charge` and :meth:`k_charge_I`

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.k_charge_J()
            2
            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: t.k_charge_J()
            8
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.k_charge_J()
            12

        TESTS::

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.k_charge_I()
            Traceback (most recent call last):
            ...
            ValueError: k-charge is not defined for skew weak tableaux
        """
        if self.parent()._skew:
            raise ValueError("k-charge is not defined for skew weak tableaux")
        stt = self.list_of_standard_cells()
        kch = 0
        for sw in stt:
            Ji = 0
            for i in range(len(sw)-1):
                c = (self._height_of_restricted_subword(sw,i+2)+1,0)
                cdi = self.parent().circular_distance((-c[0])%(self.k+1),(sw[i][1]-sw[i][0])%(self.k+1))
                cdi1 = self.parent().circular_distance((-c[0])%(self.k+1),(sw[i+1][1]-sw[i+1][0])%(self.k+1))
                if (cdi > cdi1):
                    Ji += 1
                kch += Ji + self.parent().diag(sw[i+1],c)
        return kch

    def _height_of_restricted_subword(self, sw, r):
        r"""
        Return the row of the highest addable cell of the subtableau of ``self`` with letters `\le r`
        (excluding letters `r` in standard subwords before ``sw``).

        Restrict the weak `k`-tableau ``self`` to letters `\le r` and remove all letters
        `r` that appeared in a previous standard subword selected by
        :meth:`list_of_standard_cells`.

        INPUT:

        - ``sw`` -- one of the subwords of standard cells of ``self``
        - ``r`` -- nonnegative integer

        OUTPUT:

        - a nonnegative integer

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: s = t.list_of_standard_cells()[0]; s
            [(0, 1), (1, 0), (2, 0)]
            sage: t._height_of_restricted_subword(s,2)
            1

            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: s = t.list_of_standard_cells()[0]; s
            [(0, 0), (1, 0), (0, 1), (2, 0), (0, 3), (1, 1)]
            sage: t._height_of_restricted_subword(s,4)
            2

            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: s = t.list_of_standard_cells()[0]; s
            [(0, 1), (1, 0), (2, 0), (0, 5), (3, 0), (4, 0), (5, 0)]
            sage: t._height_of_restricted_subword(s,6)
            4
        """
        L = filter(lambda v: self[v[0]][v[1]] <= r, sw)
        return max([v[0] for v in L])

class WeakTableaux_core(WeakTableaux_abstract):
    r"""
    The class of (skew) weak `k`-tableaux in the core representation of shape ``shape``
    (as `k+1`-core) and weight ``weight``.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- the shape of the `k`-tableaux represented as a `(k+1)`-core; if the
      tableaux are skew, the shape is a tuple of the outer and inner shape (both as
      `(k+1)`-cores)
    - ``weight`` -- the weight of the `k`-tableaux

    EXAMPLES::

        sage: T = WeakTableaux(3, [4,1], [2,2])
        sage: T.list()
        [[[1, 1, 2, 2], [2]]]

        sage: T = WeakTableaux(3, [[5,2,1], [2]], [1,1,1,1])
        sage: T.list()
        [[[None, None, 2, 3, 4], [1, 4], [2]],
        [[None, None, 1, 2, 4], [2, 4], [3]],
        [[None, None, 1, 2, 3], [2, 3], [4]]]
    """

    @staticmethod
    def __classcall_private__(cls, k, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_core
            sage: T = WeakTableaux_core(3, [2,1], [1,1,1])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_core(3, [[5,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        if shape == [] or shape[0] in ZZ:
            shape = (Core(shape, k+1), Core([],k+1))
        else:
            shape = tuple([Core(r,k+1) for r in shape])
        return super(WeakTableaux_core, cls).__classcall__(cls, k, shape, tuple(weight))

    def __init__(self, k, shape, weight):
        r"""
        Initializes the parent class of (skew) weak `k`-tableaux in core representation.

        INPUT:

        - ``k`` -- positive integer
        - ``outer_shape`` -- the outer shape of the `k`-tableaux represented as a
          `(k+1)`-core
        - ``weight`` -- the weight of the `k`-tableaux
        - ``inner_shape`` -- the inner shape of the skew `k`-tableaux represented as a
          `(k+1)`-core;  for straight tableaux the inner shape does not need to be
          specified (default: [])

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_core
            sage: T = WeakTableaux_core(3, [4,1], [2,2])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_core(3, [[5,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        self.k = k
        self._skew = shape[1]!=[]
        self._outer_shape = shape[0]
        self._inner_shape = shape[1]
        self._shape = (self._outer_shape, self._inner_shape)
        self._weight = weight
        self._representation = 'core'
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_core
            sage: repr(WeakTableaux_core(3, [2,1], [1,1,1]))
            'Core weak 3-Tableaux of (skew) core shape [2, 1] and weight (1, 1, 1)'
            sage: repr(WeakTableaux_core(3, [[5,2,1], [2]], [1,1,1,1]))
            'Core weak 3-Tableaux of (skew) core shape ([5, 2, 1], [2]) and weight (1, 1, 1, 1)'
        """
        return "Core weak %s-Tableaux of (skew) core shape %s and weight %s"%(self.k, self.shape(), self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: T = WeakTableaux(3, [4,1], [2,2])
            sage: T.list()
            [[[1, 1, 2, 2], [2]]]
            sage: T = WeakTableaux(3, [5,2,2], [2,2,2,1])
            sage: T.list()
            [[[1, 1, 3, 3, 4], [2, 2], [3, 3]], [[1, 1, 2, 2, 3], [2, 3], [3, 4]]]
            sage: T = WeakTableaux(3, [[5,2,2], [1]], [2,1,2,1])
            sage: T.list()
            [[[None, 1, 3, 3, 4], [1, 2], [3, 3]],
            [[None, 1, 2, 3, 3], [1, 3], [2, 4]],
            [[None, 1, 1, 2, 3], [2, 3], [3, 4]]]
        """
        for t in WeakTableaux_bounded(self.k, [self._outer_shape.to_bounded_partition(), self._inner_shape.to_bounded_partition()], self._weight):
            yield t.to_core_tableau()

    def diag(self, c, ha):
        r"""
        Return the number of diagonals strictly between cells ``c`` and ``ha`` of the same residue as ``c``.

        INPUT:

        - ``c`` -- a cell in the lattice
        - ``ha`` -- another cell in the lattice with bigger row and smaller column than `c`

        OUTPUT:

        - a nonnegative integer

        EXAMPLES::

            sage: T = WeakTableaux(4, [5,2,2], [2,2,2,1])
            sage: T.diag((1,2),(4,0))
            0
        """
        return divmod((c[1]-c[0])-(ha[1]-ha[0])-1, self.k+1)[0]

    def circular_distance(self, cr, r):
        r"""
        Return the shortest counterclockwise distance between ``cr`` and ``r`` modulo `k+1`.

        INPUT:

        - ``cr``, ``r`` -- nonnegative integers between `0` and `k`

        OUTPUT:

        - a positive integer

        EXAMPLES::

            sage: T = WeakTableaux(10, [], [])
            sage: T.circular_distance(8, 6)
            2
            sage: T.circular_distance(8, 8)
            0
            sage: T.circular_distance(8, 9)
            10
        """
        return self.k - ((r+self.k-cr)%(self.k+1))

    Element = WeakTableau_core


#Weak tableaux in terms of `k`-bounded partitions
class WeakTableau_bounded(WeakTableau_abstract):
    r"""
    A (skew) weak `k`-tableau represented in terms of `k`-bounded partitions.
    """
    @staticmethod
    def __classcall_private__(cls, t, k):
        r"""
        Implements the shortcut ``WeakTableau_bounded(t, k)`` to ``WeakTableaux_bounded(k, shape, weight)(t)``
        where ``shape`` is the shape of the tableau and ``weight`` is its weight.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_bounded
            sage: t = WeakTableau_bounded([[1,1,2],[2,3],[3]],3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            False

            sage: t = WeakTableau_bounded([[None, None, 1], [1, 2], [2]], 3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            True
        """
        if isinstance(t, cls):
            return t
        tab = SkewTableau(list(t))
        outer = tab.outer_shape()
        inner = tab.inner_shape()
        weight = tuple(tab.weight())
        if outer.conjugate().length() > k:
            raise ValueError, "The shape of %s is not %s-bounded"%(t, k)
        return WeakTableaux_bounded(k, [outer, inner], weight)(t)

    def __init__(self, parent, t):
        r"""
        Initialization of (skew) weak `k`-tableau ``t`` in `k`-bounded representation.

        INPUT:

        - ``t`` -- weak tableau in `k`-bounded representation; the input is supposed to be
          a list of lists specifying the rows of the tableau;  ``None`` is allowed as an
          entry for skew weak `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_bounded, WeakTableaux_bounded
            sage: c = WeakTableau_bounded([[1,1,2],[2,3],[3]],3)
            sage: T = WeakTableaux_bounded(3,[3,2,1],[2,2,2])
            sage: t = T([[1,1,2],[2,3],[3]]); t
            [[1, 1, 2], [2, 3], [3]]
            sage: c == t
            True
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
            sage: t.parent()
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (2, 2, 2)
            sage: TestSuite(t).run()

            sage: t = WeakTableau_bounded([[None, None, 1], [2, 4], [3]], 3)
            sage: t.shape()
            ([3, 2, 1], [2])
            sage: t.weight()
            (1, 1, 1, 1)
            sage: TestSuite(t).run()

            sage: t = T([[1,1,3],[2,2],[3]])
            Traceback (most recent call last):
            ...
            ValueError: This is not a proper weak 3-tableau
        """
        k = parent.k
        self.k = k
        self._list = [r for r in t]
        if parent._outer_shape.conjugate().length() > k:
            raise ValueError, "%s is not a %s-bounded tableau"%(t, k)
        ClonableList.__init__(self, parent, t)

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: print t._repr_diagram()
            .  .  1
            2  4
            3
            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: print t._repr_diagram()
            1  1  1
            2  2
            3
        """
        t = SkewTableau(list(self))
        return t._repr_diagram()

    def shape_core(self):
        r"""
        Return the shape of ``self`` as `(k+1)`-core.

        When the tableau is straight, the outer shape is returned as a `(k+1)`-core.
        When the tableau is skew, the tuple of the outer and inner shape is returned as
        `(k+1)`-cores.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.shape_core()
            [5, 2, 1]

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.shape_core()
            ([5, 2, 1], [2])
        """
        if self.parent()._skew:
            return tuple([r.to_core(self.k) for r in self.shape_bounded()])
        return self.shape_bounded().to_core(self.k)

    def shape_bounded(self):
        r"""
        Return the shape of ``self`` as `k`-bounded partition.

        When the tableau is straight, the outer shape is returned as a `k`-bounded
        partition.  When the tableau is skew, the tuple of the outer and inner shape is
        returned as `k`-bounded partitions.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.shape_bounded()
            [3, 2, 1]

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.shape_bounded()
            ([3, 2, 1], [2])
        """
        return self.shape()

    def check(self):
        r"""
        Check that ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[1,1],[2]], 2, representation = 'bounded')
            sage: t.check()

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.check()

        TESTS::

            sage: t = WeakTableau([[1,1,3],[2,2],[3]], 3, representation = 'bounded')
            Traceback (most recent call last):
            ...
            ValueError: This is not a proper weak 3-tableau

            sage: T = WeakTableaux(3, [3,1], [2,1], representation = 'bounded')
            sage: t = T([[None, 1,1], [2]])
            Traceback (most recent call last):
            ...
            ValueError: The inner shape of the parent does not agree with the inner shape of the tableau!

            sage: t = WeakTableau([[1,1],[1]], 3, representation = 'bounded')
            Traceback (most recent call last):
            ...
            ValueError: The tableaux is not semistandard!
        """
        t = SkewTableau(list(self))
        if not t in SemistandardSkewTableaux():
            raise ValueError("The tableaux is not semistandard!")
        if not self.parent()._weight == tuple(t.weight()):
            raise ValueError("The weight of the parent does not agree with the weight of the tableau!")
        outer = t.outer_shape()
        inner = t.inner_shape()
        if self.parent()._outer_shape != outer:
            raise ValueError("The outer shape of the parent does not agree with the outer shape of the tableau!")
        if self.parent()._inner_shape != inner:
            raise ValueError("The inner shape of the parent does not agree with the inner shape of the tableau!")
        if not t.is_k_tableau(self.k):
            raise ValueError("This is not a proper weak %s-tableau"%(self.k))

    def _is_k_tableau(self):
        r"""
        Checks whether ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t._is_k_tableau()
            True

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t._is_k_tableau()
            True
        """
        shapes = self.intermediate_shapes()
        kshapes = [ la.k_conjugate(self.k) for la in shapes ]
        return all( kshapes[i+1].contains(kshapes[i]) for i in range(len(shapes)-1) )

    def to_core_tableau(self):
        r"""
        Return the weak `k`-tableau ``self`` where the shape of each restricted tableau is a `(k+1)`-core.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,2,4],[2,3,5],[3,4],[5,6],[6],[7]], 4, representation = 'bounded')
            sage: c = t.to_core_tableau(); c
            [[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: t = WeakTableau([], 4, representation = 'bounded')
            sage: t.to_core_tableau()
            []

            sage: from sage.combinat.k_tableau import WeakTableau_bounded
            sage: t = WeakTableau([[1,1,2],[2,3],[3]], 3, representation = 'bounded')
            sage: WeakTableau_bounded.from_core_tableau(t.to_core_tableau(),3)
            [[1, 1, 2], [2, 3], [3]]
            sage: t == WeakTableau_bounded.from_core_tableau(t.to_core_tableau(),3)
            True

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.to_core_tableau()
            [[None, None, 1, 2, 4], [2, 4], [3]]
            sage: t == WeakTableau_bounded.from_core_tableau(t.to_core_tableau(),3)
            True
        """
        shapes = [ p.to_core(self.k) for p in self.intermediate_shapes() ]
        if self.parent()._skew:
            l = [[None]*i for i in shapes[0]]
        else:
            l = []
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return WeakTableau_core(l, self.k)

    @classmethod
    def from_core_tableau(cls, t, k):
        r"""
        Construct weak `k`-bounded tableau from in `k`-core tableau.

        EXAMPLES::

            sage: from sage.combinat.k_tableau import WeakTableau_bounded
            sage: WeakTableau_bounded.from_core_tableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            [[1, 1, 2], [2, 3], [3]]

            sage: WeakTableau_bounded.from_core_tableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            [[None, None, 3], [1, 4], [2]]

            sage: WeakTableau_bounded.from_core_tableau([[None,2,3],[3]], 2)
            [[None, 2], [3]]
        """
        t = SkewTableau(list(t))
        shapes = [ Core(p, k+1).to_bounded_partition() for p in intermediate_shapes(t) ]#.to_chain() ]
        if t.inner_shape() == Partition([]):
            l = []
        else:
            l = [[None]*i for i in shapes[0]]
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return cls(l, k)

    def k_charge(self, algorithm = 'I'):
        r"""
        Return the `k`-charge of ``self``.

        INPUT:

        - ``algorithm`` -- (default: "I") if "I", computes `k`-charge using the `I`
          algorithm, otherwise uses the `J`-algorithm

        OUTPUT:

        - a nonnegative integer

        For the definition of `k`-charge and the various algorithms to compute it see Section 3.3 of [LLMSSZ2013]_.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2], [2, 3], [3]], 3, representation = 'bounded')
            sage: t.k_charge()
            2
            sage: t = WeakTableau([[1, 3, 5], [2, 6], [4]], 3, representation = 'bounded')
            sage: t.k_charge()
            8
            sage: t = WeakTableau([[1, 1, 2, 4], [2, 3, 5], [3, 4], [5, 6], [6], [7]], 4, representation = 'bounded')
            sage: t.k_charge()
            12
        """
        return self.to_core_tableau().k_charge(algorithm = algorithm)


class WeakTableaux_bounded(WeakTableaux_abstract):
    r"""
    The class of (skew) weak `k`-tableaux in the bounded representation of shape ``shape``
    (as `k`-bounded partition or tuple of `k`-bounded partitions in the skew case) and
    weight ``weight``.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- the shape of the `k`-tableaux represented as a `k`-bounded partition;
      if the tableaux are skew, the shape is a tuple of the outer and inner shape each
      represented as a `k`-bounded partition
    - ``weight`` -- the weight of the `k`-tableaux

    EXAMPLES::

        sage: T = WeakTableaux(3, [3,1], [2,2], representation = 'bounded')
        sage: T.list()
        [[[1, 1, 2], [2]]]

        sage: T = WeakTableaux(3, [[3,2,1], [2]], [1,1,1,1], representation = 'bounded')
        sage: T.list()
        [[[None, None, 3], [1, 4], [2]],
        [[None, None, 1], [2, 4], [3]],
        [[None, None, 1], [2, 3], [4]]]
    """
    @staticmethod
    def __classcall_private__(cls, k, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_bounded
            sage: T = WeakTableaux_bounded(3, [2,1], [1,1,1])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_bounded(3, [[3,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        if shape == [] or shape[0] in ZZ:
            shape = (Partition(shape), Partition([]))
        else:
            shape = tuple([Partition(r) for r in shape])
        return super(WeakTableaux_bounded, cls).__classcall__(cls, k, shape, tuple(weight))

    def __init__(self, k, shape, weight):
        r"""
        Initializes the parent class of (skew) weak `k`-tableaux in bounded representation.

        INPUT:

        - ``k`` -- positive integer
        - ``shape`` -- the shape of the `k`-tableaux represented as a `k`-bounded
          partition; if the tableaux are skew, the shape is a tuple of the outer and inner
          shape each represented as a `k`-bounded partition
        - ``weight`` -- the weight of the `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_bounded
            sage: T = WeakTableaux_bounded(3, [3,1], [2,2])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_bounded(3, [[3,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        self.k = k
        self._skew = shape[1]!=[]
        self._outer_shape = Partition(shape[0])
        self._inner_shape = Partition(shape[1])
        self._shape = (self._outer_shape, self._inner_shape)
        self._weight = tuple(weight)
        self._representation = 'bounded'
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_bounded
            sage: repr(WeakTableaux_bounded(3, [2,1], [1,1,1]))
            'Bounded weak 3-Tableaux of (skew) 3-bounded shape [2, 1] and weight (1, 1, 1)'
            sage: repr(WeakTableaux_bounded(3, [[3,2,1], [2]], [1,1,1,1]))
            'Bounded weak 3-Tableaux of (skew) 3-bounded shape ([3, 2, 1], [2]) and weight (1, 1, 1, 1)'
        """
        return "Bounded weak %s-Tableaux of (skew) %s-bounded shape %s and weight %s"%(self.k, self.k, self.shape(), self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: T = WeakTableaux(3, [3,1], [2,2], representation = 'bounded')
            sage: T.list()
            [[[1, 1, 2], [2]]]
            sage: T = WeakTableaux(3, [3,2,2], [2,2,2,1], representation = 'bounded')
            sage: T.list()
            [[[1, 1, 4], [2, 2], [3, 3]], [[1, 1, 2], [2, 3], [3, 4]]]
            sage: T = WeakTableaux(3, [[3,2,2], [1]], [2,1,2,1], representation = 'bounded')
            sage: T.list()
            [[[None, 1, 4], [1, 2], [3, 3]],
            [[None, 1, 3], [1, 3], [2, 4]],
            [[None, 1, 1], [2, 3], [3, 4]]]
        """
        for t in SemistandardSkewTableaux([self._outer_shape, self._inner_shape], self._weight):
            if t.is_k_tableau(self.k):
                yield self(t)

    Element = WeakTableau_bounded

#Weak tableaux in terms of factorized permutations
class WeakTableau_factorized_permutation(WeakTableau_abstract):
    r"""
    A weak (skew) `k`-tableau represented in terms of factorizations of affine
    permutations into cyclically decreasing elements.
    """
    @staticmethod
    def straighten_input(t, k):
        r"""
        Straightens input.

        INPUT:

        - ``t`` -- a list of reduced words or a list of elements in the Weyl group of type
          `A_k^{(1)}`
        - ``k`` -- a positive integer

        EXAMPLES::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: WeakTableau_factorized_permutation.straighten_input([[2,0],[3,2],[1,0]], 3)
            (s2*s0, s3*s2, s1*s0)
            sage: W = WeylGroup(['A',4,1])
            sage: WeakTableau_factorized_permutation.straighten_input([W.an_element(),W.an_element()], 4)
            (s0*s1*s2*s3*s4, s0*s1*s2*s3*s4)

        TESTS::

            sage: WeakTableau_factorized_permutation.straighten_input([W.an_element(),W.an_element()], 3)
            Traceback (most recent call last):
            ...
            ValueError: a matrix from Full MatrixSpace of 5 by 5 dense matrices over Rational Field cannot be converted to a matrix in Full MatrixSpace of 4 by 4 dense matrices over Rational Field!
        """
        W = WeylGroup(['A', k, 1], prefix='s')
        if len(t) > 0:
            if isinstance(t[0], list) or isinstance(t[0], tuple):
                w_tuple = tuple(W.from_reduced_word(p) for p in t)
            elif t[0] not in W:
                raise ValueError, "The input must be a list of reduced words or Weyl group elements"
            else:
                w_tuple = tuple(W(r) for r in t)
        else:
            w_tuple = tuple([W.one()])
        return w_tuple

    @staticmethod
    def __classcall_private__(cls, t, k, inner_shape = []):
        r"""
        Implements the shortcut ``WeakTableau_factorized_permutation(t, k)`` to
        ``WeakTableaux_factorized_permutation(k, shape, weight)(t)``
        where ``shape`` is the shape of the tableau as a `(k+1)`-core (or a tuple of
        `(k+1)`-cores if the tableau is skew) and ``weight`` is its weight.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: t = WeakTableau_factorized_permutation([[2,0],[3,2],[1,0]], 3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
            sage: TestSuite(t).run()

            sage: t = WeakTableau_factorized_permutation([[0,3],[2,1]], 3, inner_shape = [1,1])
            sage: t.check()
            sage: TestSuite(t).run()

            sage: t = WeakTableau_factorized_permutation([], 3); t
            [1]
            sage: t.check()
            sage: TestSuite(t).run()
        """
        if isinstance(t, cls):
            return t
        W = WeylGroup(['A', k, 1], prefix='s')
        w = cls.straighten_input(t, k)
        weight =  tuple(w[i].length() for i in range(len(w)-1,-1,-1))
        inner_shape = Core(inner_shape, k+1)
        outer_shape = (W.prod(w)*W(inner_shape.to_grassmannian())).affine_grassmannian_to_core()
        return WeakTableaux_factorized_permutation(k, [outer_shape, inner_shape], weight)(w)

    def __init__(self, parent, t):
        r"""
        Initialization of (skew) weak `k`-tableau ``t`` in factorized permutation representation.

        INPUT:

        - ``t`` -- (skew) weak tableau in factorized permutation representation; the input
          can either be a list of reduced words of cyclically decreasing elements, or a
          list of cyclically decreasing elements;  when the tableau is skew, the inner
          shape needs to be specified as a `(k+1)`-core

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation, WeakTableaux_factorized_permutation
            sage: c = WeakTableau_factorized_permutation([[2,0],[3,2],[1,0]], 3)
            sage: T = WeakTableaux_factorized_permutation(3, [5,2,1],[2,2,2])
            sage: t = T([[2,0],[3,2],[1,0]]); t
            [s2*s0, s3*s2, s1*s0]
            sage: c == t
            True
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
            sage: t.parent()
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (2, 2, 2)
            sage: TestSuite(t).run()

            sage: t = WeakTableau_factorized_permutation([[2,0],[3,2]], 3, inner_shape = [2]); t
            [s2*s0, s3*s2]
            sage: t._inner_shape
            [2]
            sage: t.weight()
            (2, 2)
            sage: t.shape()
            ([5, 2, 1], [2])
            sage: TestSuite(t).run()

            sage: t = T([[3,0],[0,3],[1,0]])
            Traceback (most recent call last):
            ...
            ValueError: The outer shape of the parent does not agree with the outer shape of the tableau!

            sage: t = WeakTableau_factorized_permutation([], 3); t
            [1]
            sage: t.parent()._outer_shape
            []
            sage: t.parent()._weight
            (0,)
        """
        self.k = parent.k
        self._inner_shape = parent._inner_shape
        ClonableList.__init__(self, parent, self.straighten_input(t, parent.k))

    def shape_core(self):
        r"""
        Return the shape of ``self`` as a `(k+1)`-core.

        When the tableau is straight, the outer shape is returned as a core.
        When the tableau is skew, the tuple of the outer and inner shape is returned as
        cores.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.shape_core()
            [5, 2, 1]

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.shape()
            ([5, 2, 1], [2])
        """
        return self.shape()

    def shape_bounded(self):
        r"""
        Return the shape of ``self`` as a `k`-bounded partition.

        When the tableau is straight, the outer shape is returned as a `k`-bounded
        partition.  When the tableau is skew, the tuple of the outer and inner shape is
        returned as `k`-bounded partitions.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.shape_bounded()
            [3, 2, 1]

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.shape_bounded()
            ([3, 2, 1], [2])
        """
        if self.parent()._skew:
            return tuple([r.to_bounded_partition() for r in self.shape_core()])
        return self.shape_core().to_bounded_partition()

    def check(self):
        r"""
        Check that ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.check()

        TESTS::

            sage: t = WeakTableau([[2,0],[3,2]], 3, representation = 'factorized_permutation')
            Traceback (most recent call last):
            ...
            ValueError: Error! this only works on type 'A' affine Grassmannian elements

            sage: T = WeakTableaux(3, [4,1], [2,1], representation = 'factorized_permutation')
            sage: t = T([[2],[1],[0]])
            Traceback (most recent call last):
            ...
            ValueError: The weight of the parent does not agree with the weight of the tableau!
        """
        weight =  tuple(self[i].length() for i in range(len(self)-1,-1,-1))
        if not self.parent()._weight == weight:
            raise ValueError("The weight of the parent does not agree with the weight of the tableau!")
        W = self[0].parent()
        outer = (W.prod(self)*W((self._inner_shape).to_grassmannian())).affine_grassmannian_to_core()
        if self.parent()._outer_shape != outer:
            raise ValueError("The outer shape of the parent does not agree with the outer shape of the tableau!")
        if not self._is_k_tableau():
            raise ValueError("This is not a proper weak %s-tableau"%(self.k))

    def _is_k_tableau(self):
        r"""
        Checks whether ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t._is_k_tableau()
            True

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t._is_k_tableau()
            True
        """
        W = self[0].parent()
        if (W.prod(self)*W(self.parent()._inner_shape.to_grassmannian())).is_affine_grassmannian():
            return all( r.is_pieri_factor() for r in self )
        return False

    def to_core_tableau(self):
        r"""
        Return the weak `k`-tableau ``self`` where the shape of each restricted tableau is a `(k+1)`-core.

        EXAMPLES::

            sage: t = WeakTableau([[0], [3,1], [2,1], [0,4], [3,0], [4,2], [1,0]], 4, representation = 'factorized_permutation'); t
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]
            sage: c = t.to_core_tableau(); c
            [[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: t = WeakTableau([[]], 4, representation = 'factorized_permutation'); t
            [1]
            sage: t.to_core_tableau()
            []

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = 'factorized_permutation')
            sage: WeakTableau_factorized_permutation.from_core_tableau(t.to_core_tableau(), 3)
            [s2*s0, s3*s2, s1*s0]
            sage: t == WeakTableau_factorized_permutation.from_core_tableau(t.to_core_tableau(), 3)
            True

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.to_core_tableau()
            [[None, None, 1, 1, 2], [1, 2], [2]]
            sage: t == WeakTableau_factorized_permutation.from_core_tableau(t.to_core_tableau(), 3)
            True
        """
        W = self[0].parent()
        factor = W(self._inner_shape.to_grassmannian())
        shapes = [factor]
        for i in range(len(self)-1,-1,-1):
            factor = self[i]*factor
            shapes += [factor.affine_grassmannian_to_core()]
        if self.parent()._skew:
            l = [[None]*i for i in self._inner_shape]
        else:
            l = []
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return WeakTableau_core(l, self.k)

    @classmethod
    def from_core_tableau(cls, t, k):
        r"""
        Construct weak factorized affine permutation tableau from a `k`-core tableau.

        EXAMPLES::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: WeakTableau_factorized_permutation.from_core_tableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            [s2*s0, s3*s2, s1*s0]
            sage: WeakTableau_factorized_permutation.from_core_tableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]
            sage: WeakTableau_factorized_permutation.from_core_tableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            [s0*s3, s2*s1]
        """
        t = SkewTableau(list(t))
        shapes = [ Core(p, k+1).to_grassmannian() for p in intermediate_shapes(t) ] #t.to_chain() ]
        perms = [ shapes[i]*(shapes[i-1].inverse()) for i in range(len(shapes)-1,0,-1)]
        return cls(perms, k, inner_shape = t.inner_shape())

    def k_charge(self, algorithm = 'I'):
        r"""
        Return the `k`-charge of ``self``.

        OUTPUT:

        - a nonnegative integer

        EXAMPLES::

            sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = 'factorized_permutation')
            sage: t.k_charge()
            2
            sage: t = WeakTableau([[0],[3],[2],[1],[3],[0]], 3, representation = 'factorized_permutation')
            sage: t.k_charge()
            8
            sage: t = WeakTableau([[0],[3,1],[2,1],[0,4],[3,0],[4,2],[1,0]], 4, representation = 'factorized_permutation')
            sage: t.k_charge()
            12
        """
        return self.to_core_tableau().k_charge(algorithm = algorithm)


class WeakTableaux_factorized_permutation(WeakTableaux_abstract):
    r"""
    The class of (skew) weak `k`-tableaux in the factorized permutation representation of shape ``shape`` (as `k+1`-core
    or tuple of `(k+1)`-cores in the skew case) and weight ``weight``.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- the shape of the `k`-tableaux represented as a `(k+1)`-core;
      in the skew case the shape is a tuple of the outer and inner shape both as `(k+1)`-cores
    - ``weight`` -- the weight of the `k`-tableaux

    EXAMPLES::

        sage: T = WeakTableaux(3, [4,1], [2,2], representation = 'factorized_permutation')
        sage: T.list()
        [[s3*s2, s1*s0]]

        sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], representation = 'factorized_permutation')
        sage: T.list()
        [[s0, s4, s3, s4*s2], [s0, s3, s4, s3*s2], [s3, s0, s4, s3*s2]]
    """
    @staticmethod
    def __classcall_private__(cls, k, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_factorized_permutation
            sage: T = WeakTableaux_factorized_permutation(3, [2,1], [1,1,1])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_factorized_permutation(4, [[6,2,1], [2]], [2,1,1,1])
            sage: TestSuite(T).run()
        """
        if shape == [] or shape[0] in ZZ:
            shape = (Core(shape, k+1), Core([],k+1))
        else:
            shape = tuple([Core(r,k+1) for r in shape])
        return super(WeakTableaux_factorized_permutation, cls).__classcall__(cls, k, shape, tuple(weight))

    def __init__(self, k, shape, weight):
        r"""
        Initializes the parent class of weak `k`-tableaux in factorized permutation representation.

        INPUT:

        - ``k`` -- positive integer
        - ``shape`` -- the shape of the `k`-tableaux represented as a `(k+1)`-core;
          in the skew case the shape is a tuple of the outer and inner shape both as
          `(k+1)`-cores
        - ``weight`` -- the weight of the `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_factorized_permutation
            sage: T = WeakTableaux_factorized_permutation(3, [4,1], [2,2])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_factorized_permutation(4, [[6,2,1], [2]], [2,1,1,1])
            sage: TestSuite(T).run()
        """
        self.k = k
        self._skew = shape[1]!=[]
        self._outer_shape = Core(shape[0], k+1)
        self._inner_shape = Core(shape[1], k+1)
        self._shape = (self._outer_shape, self._inner_shape)
        self._weight = weight
        self._representation = 'factorized_permutation'
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_factorized_permutation
            sage: repr(WeakTableaux_factorized_permutation(3, [2,1], [1,1,1]))
            'Factorized permutation (skew) weak 3-Tableaux of shape [2, 1] and weight (1, 1, 1)'
            sage: repr(WeakTableaux_factorized_permutation(4, [[6,2,1], [2]], [2,1,1,1]))
            'Factorized permutation (skew) weak 4-Tableaux of shape ([6, 2, 1], [2]) and weight (2, 1, 1, 1)'
        """
        return "Factorized permutation (skew) weak %s-Tableaux of shape %s and weight %s"%(self.k, self.shape(), self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: T = WeakTableaux(3, [4,1], [2,2], representation = 'factorized_permutation')
            sage: T.list()
            [[s3*s2, s1*s0]]
            sage: T = WeakTableaux(3, [5,2,2], [2,2,2,1], representation = 'factorized_permutation')
            sage: T.list()
            [[s0, s3*s2, s0*s3, s1*s0], [s3, s2*s0, s3*s2, s1*s0]]
            sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], representation = 'factorized_permutation')
            sage: T.list()
            [[s0, s4, s3, s4*s2], [s0, s3, s4, s3*s2], [s3, s0, s4, s3*s2]]
        """
        for t in WeakTableaux_core(self.k, self.shape(), self._weight):
            yield WeakTableau_factorized_permutation.from_core_tableau(t, self.k)

    Element = WeakTableau_factorized_permutation

def intermediate_shapes(t):
    r"""
    Return the intermediate shapes of tableau ``t``.

    A (skew) tableau with letters `1,2,\ldots,\ell` can be viewed as a sequence of shapes,
    where the `i`-th shape is given by the shape of the subtableau on letters `1,2,\ldots,i`.
    The output is the list of these shapes.

    EXAMPLES::

        sage: from sage.combinat.k_tableau import intermediate_shapes
        sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
        sage: intermediate_shapes(t)
        [[], [2], [4, 1], [5, 2, 1]]

        sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
        sage: intermediate_shapes(t)
        [[2], [2, 1], [3, 1, 1], [4, 1, 1], [5, 2, 1]]
   """
    shapes = []
    t = SkewTableau(list(t))
    for i in range(len(t.weight())+1):
        shapes += [ t.restrict(i).outer_shape()]
    return shapes
