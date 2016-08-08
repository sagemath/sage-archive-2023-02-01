r"""
`k`-regular Sequences
"""

import itertools
from sage.misc.cachefunc import cached_function, cached_method
from sage.structure.element import Element

class kRegularSequence(Element):

    def __init__(self, parent, matrices, initial=None, selection=None,
                 output_function=None, transpose=False):
        r"""
        TESTS::

            sage: from sage.combinat.k_regular_sequence import kRegularSequences
            sage: Seq2 = kRegularSequences(2, ZZ)
            sage: Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:      Matrix([[0, 1]]), Matrix([[1], [0]]),
            ....:      lambda o: o[0, 0], transpose=True)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...
        """
        super(kRegularSequence, self).__init__(parent=parent)

        def tr(M):
            try:
                return M.transpose() if transpose else M
            except AttributeError:
                return M

        self.matrices = tuple(tr(M) for M in matrices)
        self.k = len(self.matrices)
        self.d = self.matrices[0].nrows()
        if not all(M.dimensions() == (self.d, self.d) for M in self.matrices):
            raise ValueError

        if not transpose:
            self.initial = initial
            self.selection = selection
        else:
            self.initial = tr(selection)
            self.selection = tr(initial)

        if output_function is None:
            self.output_function = lambda o: o
        else:
            self.output_function = output_function


    def _repr_(self):
        # TODO
        from sage.arith.srange import xsrange
        return '{}-regular sequence '.format(self.parent().k) +\
            ', '.join(repr(self[n]) for n in xsrange(10)) + ', ...'


    def info(self):
        r"""
        EXAMPLES:

            sage: from sage.combinat.k_regular_sequence import kRegularSequences
            sage: Seq2 = kRegularSequences(2, ZZ)
            sage: Seq2.guess(lambda n: sum(n.digits(2))).info()
            matrices:
            (
            [1 0]  [ 0 -1]
            [0 1], [ 1  2]
            )
            initial:
            (0, 1)
            selection:
            (1, 0)
        """
        from sys import displayhook
        print('matrices:')
        displayhook(self.matrices)
        print('initial:')
        displayhook(self.initial)
        print('selection:')
        displayhook(self.selection)


    @cached_method
    def __getitem__(self, n):
        result = self.product_of_matrices(n)
        if self.initial is not None:
            result = self.initial * result
        if self.selection is not None:
            result = result * self.selection
        return self.output_function(result)


    @cached_method
    def product_of_matrices(self, m):
        k = self.parent().k
        if m < 0:
            raise ValueError
        if 0 <= m < k:
            return self.matrices[m]
        n = m // k
        r = m - n*k
        return self.matrices[r] * self.product_of_matrices(n)


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

class kRegularSequences(UniqueRepresentation, Parent):

    Element = kRegularSequence

    def __init__(self, k, base, category=None):
        r"""
        TESTS::

            sage: from sage.combinat.k_regular_sequence import kRegularSequences
            sage: kRegularSequences(2, ZZ)
            Set of 2-regular sequences over Integer Ring
        """
        from sage.categories.sets_cat import Sets
        self.k = k
        super(kRegularSequences, self).__init__(category=category or Sets(),
                                                base=base)


    def _repr_(self):
        return 'Set of {}-regular sequences over {}'.format(self.k, self.base())

