r"""
Fully commutative Coxeter group elements
"""
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import NormalizedClonableList
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from .coxeter_matrix import CoxeterMatrix
from .cartan_type import CartanType
from collections import deque
from sage.combinat.posets.posets import Poset
import itertools

class FullyCommutativeCoxeterElement(NormalizedClonableList):
    r"""
    A (reduced word of a) fully commutative element in the Coxeter system.
    """

    # Methods required as a subclass of NormalizedClonableList:
    def check(self):
        r"""
        Determine if the word ``self`` actually represents a fully commutative
        reduced word. Raises a ValueError if this is not the case.
        """
        if not self.parent()._is_fully_commutative(self._get_list()):
            raise ValueError('list does not represent a fully commutative word.')

    def normalize(self):
        r"""
        Normalize ``self`` to Cartier--Foata normal form.
        """
        self._require_mutable()

        out_word = []
        
        while len(self) > 0:
            for s in self.parent().index_set():
                i = self.find_left_descent(s)
                if i is not None:
                    out_word.append(s)
                    self.pop(i)

        self._set_list(out_word)

    # Operations on fully commutative elements:

    def find_left_descent(self, s):
        r"""
        Determine if ``s`` is a left descent of ``self``, and find the index of
        its occurrence in ``self``.

        INPUT:

        - ``s`` -- integer; the index of the generator (element of
          ``self.parent().index_set()``).

        OUTPUT: 
        
        The index of ``s`` in ``self`` if ``s`` is a left-descent,
        ``None`` if ``s`` is not a left descent.
        """
        m = self.parent().coxeter_matrix()
        for (i, t) in enumerate(self):
            if t == s and not any(m[x, t] > 2 for x in self[:i]):
                return i
        return None

    def has_left_descent(self, s):
        r"""
        Determine if ``s`` is a left descent of ``self``.
        """
        return self.find_left_descent(s) is not None

    def left_descents(self):
        r"""
        Obtain the set of left descents of ``self``.

        OUTPUT:

        A set of integers representing the indices of generators which are left
        descents of ``self``.
        """
        m = self.parent().coxeter_matrix()
        out = set()
        for (i, t) in enumerate(self):
            if not any(m[x, t] > 2 for x in self[:i]):
                out.add(t)
        return out

    def still_reduced_fc_after_prepending(self, s):
        r"""
        Determine if ``self`` prepended with ``s`` is still a reduced word of a
        fully commutative element in the Coxeter system.

        INPUT:

        - ``s`` -- integer; the index of the generator (element of
          ``self.parent().index_set()``).
        """
        m = self.parent().coxeter_matrix()
        if self.has_left_descent(s):
            return False

        # Find the first letter in that doesn't commute with s.
        try:
            (j, t) = next((i, x) for (i, x) in enumerate(self) if m[s, x] >= 3)
        except StopIteration:
            return True

        u = self.clone()
        u._set_list(self[j:])
        for c in range(m[s, t] - 1):
            letter = t if c % 2 == 0 else s
            i = u.find_left_descent(letter)
            if i is not None:
                u.pop(i)
            else:
                return True

        return False

    def heap(self, **kargs):
        r"""
        Create the heap of ``self``.

        The heap of a fully commutative word `w` is a partial order on the
        indices of a reduced word for `w`.

        OPTIONAL ARGUMENTS:

        - ``one_index`` -- boolean (default False); make the Poset on the
          indices {1, 2, ..., n} instead of {0, 1, ..., n-1}.

        - ``display_labeling`` -- boolean (default False); make the elements of
          the resulting Poset be tuples `(i, w_i)` where `w_i` is the letter at
          index `i` in ``self``.

        .. NOTE:

        Fully commutative elements have the property that the heap of any
        reduced word is unique up to Poset isomorphism. Of course, the
        particular one returned from this method is created from the
        Cartier--Foata normal form of the reduced word.

        OUTPUT: Poset; the partially ordered set on the indices of ``self``.

        EXAMPLES:

        Create the heap of a fully commutative element in `A_5` ::

            sage: FC = FullyCommutativeCoxeterElements(['A', 5])
            sage: FC([1, 4, 3, 5, 4, 2]).heap().cover_relations()
            [[1, 2], [1, 3], [2, 4], [3, 4], [3, 5], [0, 5]]
            sage: FC([1, 4, 3, 5, 4, 2]).heap(one_index=True).cover_relations()
            [[2, 3], [2, 4], [3, 5], [4, 5], [4, 6], [1, 6]]
        """
        m = self.parent().coxeter_matrix()

        one_index = kargs.get('one_index', False)
        display_labeling = kargs.get('display_labeling', False)
        
        # Elements are the indices {1, 2, ..., n}
        elements = list(range(1, len(self) + 1)) if one_index else list(range(len(self)))
        
        # Get the letter of w at a certain index, respecting the indexing convention
        def letter(index):
            return self[index-1] if one_index else self[index]

        # The partial order is generated by the relations i \prec j for all i < j with m(w_i, w_j) != 2
        relations = [(i, j) for [i, j] in itertools.product(elements, repeat=2) if i < j and m[letter(i), letter(j)] != 2]
        
        # Create the poset from the transitive closure of these relations
        p = Poset((elements, relations))
        if not display_labeling:
            return p
        else:
            return p.relabel(lambda i: (i, letter(i)))

    def n_value(self):
        r"""
        Calculate the n-value of ``self``.

        The *n-value* of a fully commutative element is the *width* (length
        of the longest antichain) of its heap.
        """
        return self.heap().width()
    
    def plot(self):
        r"""
        Display a semantically helpful rendering of the heap of ``self``.

        Specifically, if ``self`` is the word `w`, render the Hasse diagram of
        the heap of `w` where each vertex is labeled by `w_i` instead of `i`.
        Furthermore, vertices are arranged in level sets, where the "level" of
        an element is the maximum length of a chain ending in that element.

        OUTPUT: GraphicsObject
        """
        import sage.plot.all as plot
        from sage.modules.free_module_element import vector

        m = self.parent().coxeter_matrix()
        letters = self.parent().index_set()
        graphics = []
        
        # Returns a line between two points that respects the space of the "nodes"
        def short_line_between(a, b):
            a, b = vector(a), vector(b)
            d = (b - a).normalized()
            a2 = a + 0.1*d
            b2 = b - 0.1*d
            return plot.line([a2, b2], color='black')
        
        h = self.heap()
        levels = h.level_sets()
        letters_at_level = [set(self[i] for i in level) for level in levels]

        for (level_zero_index, members) in enumerate(levels):
            level = level_zero_index + 1
            for i in members:
                x = self[i]

                # Draw the node
                graphics.append(plot.circle((x, level), 0.1))
                graphics.append(plot.text(str(x), (x, level)))

                neighbors = {z for z in letters if m[x, z] >= 3}
                for other in neighbors:
                    highest_level = max((j + 1 for j in range(level_zero_index) if other in letters_at_level[j]), default=None)
                    if highest_level:
                        graphics.append(short_line_between((other, highest_level), (x, level)))
        
        g = sum(graphics)
        g.axes(False)
        return g


class FullyCommutativeCoxeterElements(Parent):
    r"""
    The combinatorial class of reduced words of fully commutative elements in a
    Coxeter group.

    Given a Coxeter system with Coxeter matrix `m`, a *fully commutative*
    element `w` is an element with the property that any reduced word for `w`
    can be transformed into any other reduced word for `w` by applying
    commutation relations only. Equivalently, `w` is an element such that no
    reduced word for `w` contains a long braid.

    Fully commutative elements of a Coxeter system have a unique canonical form
    called the *Cartier--Foata normal form*. All elements of this class are
    automatically normalized to Cartier--Foata form, and equality of elements is
    decied by the equality of their normal forms.

    Certain Coxeter systems are *FC-finite*, meaning they contain finitely many
    fully commutative elements. The seven FC-finite Coxeter group families are
    (non-affine) `A_n`, `B_n = C_n`, `D_n`, `E_n`, `F_n`, `H_n` and `I_2(m)`.

    INPUT:

    - ``data`` -- either data describing a Cartan type, or a CoxeterMatrix. Same
      datum as the :func:`sage.combinat.root_system.coxeter_group.CoxeterGroup`
      constructor.

    OUTPUT:

    The class of reduced words of fully commutative elements in the Coxeter
    group described by ``data``. This will belong to either the category of
    infinite enumerated sets or finite enumerated sets depending on if the group
    is FC-finite.

    EXAMPLES:

    Enumerate the reduced words of fully commutative elements in `A_3` ::

        sage: FC = FullyCommutativeReducedCoxeterWords(['A', 3])
        sage: FC.category()
        Category of finite enumerated sets
        sage: FC.list()
        [[],
         [2],
         [3],
         [1],
         [2, 3],
         [1, 2],
         [3, 2],
         [1, 3],
         [2, 1],
         [1, 2, 3],
         [3, 2, 1],
         [1, 3, 2],
         [2, 3, 1],
         [2, 3, 1, 2]]

    Count the FC elements in `B_8` ::

        sage: FC = FullyCommutativeCoxeterElements(['B', 8])
        sage: len(FC) # long time (7 seconds)
        14299

    Iterate through the fully commutative elements of length up to 3 in the non
    FC-finite group affine `A_2` ::

        sage: FC = FullyCommutativeCoxeterElements(['A', 2, 1])
        sage: FC.category()
        Category of infinite enumerated sets
        sage: list(FC.iterate_to_length(2))
        [[], [1], [2], [0], [1, 0], [0, 2], [0, 1], [1, 2], [2, 1], [2, 0]]

    Constructing an element that is not fully commutative throws an error ::

        sage: FC = FullyCommutativeCoxeterElements(['A', 3])
        sage: FC([1,2,1])
        Traceback (most recent call last):
        ...
        ValueError: list does not represent a fully commutative word.

    Elements are normalized to Cartier--Foata normal form upon construction ::

        sage: FC = FullyCommutativeCoxeterElements(['A', 3])
        sage: FC([2, 1, 3, 2])
        [2, 3, 1, 2]

    """
    def __init__(self, data):
        if isinstance(data, CoxeterMatrix):
            self._matrix = data
        else:
            try:
                t = CartanType(data)
            except (TypeError, ValueError):
                raise ValueError('Input must be a CoxeterMatrix or data describing a Cartan type.')
            self._matrix = t.coxeter_matrix()

        self._index_set = sorted(self._matrix.index_set())

        # Determine if this group is FC-finite.
        category = InfiniteEnumeratedSets()
        if self._matrix.is_finite():
            category = FiniteEnumeratedSets()
        else:
            cartan_type = self._matrix.coxeter_type().cartan_type()
            family, rank, affine = cartan_type.type(), cartan_type.rank(), cartan_type.is_affine()
            if not affine and (rank == 2 or family in {'A', 'B', 'C', 'D', 'E', 'F', 'H'}):
                category = FiniteEnumeratedSets()

        Parent.__init__(self, category=category)

    def _element_constructor_(self, lst):
        return self.element_class(self, lst)

    Element = FullyCommutativeCoxeterElement

    def coxeter_matrix(self):
        r"""
        Obtain the Coxeter matrix of the associated Coxter system.

        OUTPUT: CoxeterMatrix
        """
        return self._matrix

    def index_set(self):
        r"""
        Obtain the index set of the generators / simple reflections of the
        associated Coxeter system.

        OUTPUT: iterable of integers
        """
        return self._index_set

    def __iter__(self):
        r"""
        Enumerate the elements of this set by length, starting with the empty
        word and, if the group is FC-finite, ending with the longest fully
        commutative element.
        """
        m = self.coxeter_matrix()

        empty_word = self.element_class(self, [], check=False)
        letters = self.index_set()

        recent_words = {empty_word}
        yield empty_word

        length = 1
        while True:
            new_words = set()
            for w in recent_words:
                for s in letters:
                    if w.still_reduced_fc_after_prepending(s):
                        sw = self.element_class(self, [s] + list(w), check=False)
                        new_words.add(sw)

            if len(new_words) == 0:
                return
            for w in new_words:
                yield w
            recent_words = new_words
            length += 1

    def iterate_to_length(self, length):
        r"""
        Iterate through the elements of this class up to a maximum length.

        INPUT:

        - ``length`` -- integer; maximum length of element to generate.
        """
        assert length >= 0
        for w in self:
            if len(w) > length:
                break
            yield w

    def _is_fully_commutative(self, w):
        matrix = self.coxeter_matrix()

        def contains_long_braid(w):
            for i in range(0, len(w)-2):
                a = w[i]
                b = w[i+1]
                m = matrix[a, b]
                if m > 2 and i+m <= len(w):
                    ab_braid = (a, b) * (m // 2) + ((a,) if m % 2 == 1 else ())
                    if w[i:i+m] == ab_braid:
                        return True
            return False

        def commute_once(word, i):
            return word[:i] + (word[i+1], word[i]) + word[i+2:]

        w = tuple(w)

        if contains_long_braid(w):
            return False
        else:
            l, checked, queue = len(w), {w}, deque([w]) 
            while queue:
                word = queue.pop()
                for i in range(l-1):
                    a, b = word[i], word[i+1]
                    if matrix[a, b] == 2:
                        new_word = commute_once(word, i)
                        if new_word not in checked:
                            if contains_long_braid(new_word):
                                return False
                            else:
                                checked.add(new_word)
                                queue.appendleft(new_word)
            return True
