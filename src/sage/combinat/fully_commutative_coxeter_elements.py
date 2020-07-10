r"""
Fully commutative Coxeter group elements
"""
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import NormalizedClonableList
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from .root_system.coxeter_matrix import CoxeterMatrix
from .root_system.cartan_type import CartanType
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

        EXAMPLES:

        Observe the normal form of three equivalent FC words in `B_5` ::

            sage: FC = FullyCommutativeCoxeterElements(['B', 5])
            sage: FC([1, 4, 3, 5, 2, 4, 3])
            [1, 4, 3, 5, 2, 4, 3]
            sage: FC([4, 1, 3, 5, 2, 4, 3])
            [1, 4, 3, 5, 2, 4, 3]
            sage: FC([4, 1, 5, 3, 4, 2, 3])
            [1, 4, 3, 5, 2, 4, 3]
            sage: FC([1, 4, 3, 5, 2, 4, 3]) == FC([4, 1, 3, 5, 2, 4, 3]) == FC([4, 1, 5, 3, 4, 2, 3])
            True
        """
        self._require_mutable()

        out_word = []
        
        while len(self) > 0:
            fronts = self.descents()
            out_word.extend(sorted(fronts))
            for s in fronts:
                self.remove(s)

        self._set_list(out_word)

    # Operations on fully commutative elements:

    def find_descent(self, s, side='left'):
        r"""
        Determine if ``s`` is a left descent of ``self``, and find the index of
        its occurrence in ``self``.

        INPUT:

        - ``s`` -- integer; the index of the generator (element of
          ``self.parent().index_set()``).

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default 'left') if 'right', find ``s`` as a right
          descent.

        OUTPUT: 

        The index of ``s`` in ``self`` if ``s`` is a descent on the appropriate
        side, ``None`` if ``s`` is not a descent.
        """
        m = self.parent().coxeter_matrix()
        view = list(self) if side == 'left' else self[::-1]
        for (i, t) in enumerate(view):
            if t == s and not any(m[x, t] > 2 for x in view[:i]):
                return i
        return None

    def has_descent(self, s, side='left'):
        r"""
        Determine if ``s`` is a descent on the appropriate side of ``self``.

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default 'left') if 'right', determine if ``self``
          has ``s`` as a right descent.

        OUTPUT: boolean
        """
        return self.find_descent(s, side=side) is not None

    def descents(self, side='left'):
        r"""
        Obtain the set of left or right descents of ``self``.

        OUTPUT:

        A set of integers representing the indices of generators which are left
        descents of ``self``.

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default 'left') if 'right', find right descents.

        EXAMPLES:

        Determine the left descents of an FC element in `A_5`::

            sage: FC = FullyCommutativeCoxeterElements(['A', 5])
            sage: sorted(FC([1, 4, 3, 5, 2, 4, 3]).descents())
            [1, 4]
        """
        view = list(self) if side == 'left' else self[::-1]
        m = self.parent().coxeter_matrix()
        out = set()
        for (i, t) in enumerate(view):
            if not any(m[x, t] > 2 for x in view[:i]):
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
        if self.has_descent(s):
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
            i = u.find_descent(letter)
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
          indices `\{1, 2, \dots, n\}` instead of `\{0, 1, \dots, n-1\}`.

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
            sage: FC([1, 4, 3, 5, 2, 4]).heap().cover_relations()
            [[1, 2], [1, 3], [2, 5], [2, 4], [3, 5], [0, 4]]
            sage: FC([1, 4, 3, 5, 4, 2]).heap(one_index=True).cover_relations()
            [[2, 3], [2, 4], [3, 6], [3, 5], [4, 6], [1, 5]]
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

        m = self.parent().coxeter_matrix()
        letters = self.parent().index_set()
        graphics = []
        
        h = self.heap()
        levels = h.level_sets()
        letters_at_level = [set(self[i] for i in level) for level in levels]

        for (level_zero_index, members) in enumerate(levels):
            level = level_zero_index + 1
            for i in members:
                x = self[i]

                # Draw the node
                graphics.append(plot.circle((x, level), 0.1, fill=True, facecolor='white', edgecolor='blue', zorder=1))
                graphics.append(plot.text(str(x), (x, level), color='blue', zorder=2))

                neighbors = {z for z in letters if m[x, z] >= 3}
                for other in neighbors:
                    highest_level = max((j + 1 for j in range(level_zero_index) if other in letters_at_level[j]), default=None)
                    if highest_level:
                        graphics.append(plot.line([(other, highest_level), (x, level)], color='black', zorder=0))
        
        g = sum(graphics)
        g.axes(False)
        return g

    def coset_decomposition(self, J, side='left'):
        r"""
        Perform a coset decomposition of ``self`` with repsect to the subgroup
        generated by ``J``.

        If self is `w`, then this method returns a tuple `(w_J, w^J)` such that
        `w_J \in \langle J \rangle`, `w^J \not\in \langle J \rangle`, and `w =
        w_J w^J`.

        .. NOTE:

        If ``J`` generates a parabolic subgroup (i.e., ``J`` is a pair of
        non-commuting generators), then the decomposition returned by this
        function is unique.

        INPUT:

        - ``J`` -- collection of integers; the indices of the generators that
          generate the subgroup with repsect to which to perform the
          decomposition

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default 'left') if 'right', perform a right coset
          decomposition and return the tuple `(w^J, w_J)` where `w = w^J w_J`.

        EXAMPLES:

        Do a left and right coset decomposition in `B_6` with respect to the
        parabolic subgroup generated by {5, 6} ::

            sage: FC = FullyCommutativeCoxeterElements(['B', 6])
            sage: x = FC([1, 6, 2, 5, 4, 6, 5])
            sage: x.coset_decomposition((5, 6))
            ([6, 5, 6], [1, 2, 4, 5])
            sage: x.coset_decomposition((5, 6), side='right')
            ([1, 6, 2, 5, 4], [6, 5])
        """
        string = []  # The J-string
        remaining = self.clone()  # The remainder

        if side == 'right':
            remaining._set_list(remaining[::-1])

        while True:
            x = next((x for x in J if remaining.has_descent(x, side='left')), None)
            if x is not None:
                string.append(x)
                remaining.remove(x)
            else:
                break
        
        if side == 'right':
            remaining._set_list(remaining[::-1])
            string = string[::-1]

        string = self.parent().element_class(self.parent(), string, check=False)
        remaining.set_immutable()

        return (string, remaining) if side == 'left' else (remaining, string)

    def _star_operation_inner(self, J, direction, side):
        mst = self.parent().coxeter_matrix()[J[0], J[1]]
        # Decompose and assign the appropriate names based on the side of the decomposition.
        if side == 'left':
            (string, remaining) = self.coset_decomposition(J, side=side)
        elif side == 'right':
            (remaining, string) = self.coset_decomposition(J, side=side)

        cur_string = list(string)
        
        if direction == 'down' and 2 <= len(string) <= mst - 1:
            # Decrease the length of the J-string
            new_string = cur_string[1:] if side == 'left' else cur_string[:-1]
        elif direction == 'up' and 1 <= len(string) <= mst - 2:
            # Increase the length of the J-string
            ending_letter = cur_string[0] if side == 'left' else cur_string[-1]
            other = next(x for x in J if x != ending_letter)
            new_string = [other] + cur_string if side == 'left' else cur_string + [other]
        else:
            return None

        # concatenate in the appropriate order
        combined_data = new_string + list(remaining) if side == 'left' else list(remaining) + new_string
        
        return self.parent().element_class(self.parent(), combined_data, check=False)

    def upper_star(self, J, side='left'):
        r"""
        Perform an upper star operation on ``self`` with respect to the
        parabolic subgroup generated by ``J``.

        If ``self`` is `w` and `w` has the (unique) coset decomposition `w_J
        w^J` (where `w_J` is necessarily an alternating string in `\{s, t\}`
        where `J = \{s, t\}`), an upper star operation on `w` increases the
        length of `w_J` by adding element on the left (if such an addition would
        not result in a long braid), and a lower star operation decreases the
        length of `w_J` by removing the left-most element (if such an removal
        would not result in `w_J` being empty).

        INPUT:

        - ``J`` -- tuple of integers; a *pair* of indices of *non-commuting*
          generators, such that ``J`` generates a parabolic subgroup.

        OPTIONAL ARGUEMNTS:

        - ``side`` -- string (default 'left') if 'right', perform a right upper
          star operation.

        OUTPUT: 

        The result of the star operation if it is defined on ``self``, otherwise
        ``None``.

        EXAMPLES:

        Perform right star operations on an element in `B_6`, with repsect to
        the parabolic subgroup `\langle 5, 6 \rangle` (the long braid pair). Here both
        star operations are defined ::

            sage: FC = FullyCommutativeCoxeterElements(['B', 6])
            sage: x = FC([1, 6, 2, 5, 4, 6, 5])
            sage: x.upper_star((5, 6), side='right')
            [1, 6, 2, 5, 4, 6, 5, 6]
            sage: x.lower_star((5, 6), side='right')
            [1, 6, 2, 5, 4, 6]

        For the same element, a left upper star operation is not defined ::

            sage: FC = FullyCommutativeCoxeterElements(['B', 6])
            sage: x = FC([1, 6, 2, 5, 4, 6, 5])
            sage: x.upper_star((5, 6)) is None
            True
        """
        return self._star_operation_inner(J, 'up', side)

    def lower_star(self, J, side='left'):
        r"""
        Perform a lower star operation on ``self``; see ``.upper_star`` for
        definitions and documentation.
        """
        return self._star_operation_inner(J, 'down', side)


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

        sage: FC = FullyCommutativeCoxeterElements(['A', 3])
        sage: FC.category()
        Category of finite enumerated sets
        sage: FC.list()
        [[],
         [1],
         [2],
         [3],
         [2, 1],
         [1, 3],
         [1, 2],
         [3, 2],
         [2, 3],
         [3, 2, 1],
         [2, 1, 3],
         [1, 3, 2],
         [1, 2, 3],
         [2, 1, 3, 2]]

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
        [[], [0], [1], [2], [1, 0], [2, 0], [0, 1], [2, 1], [0, 2], [1, 2]]

    Constructing an element that is not fully commutative throws an error ::

        sage: FC = FullyCommutativeCoxeterElements(['A', 3])
        sage: FC([1,2,1])
        Traceback (most recent call last):
        ...
        ValueError: list does not represent a fully commutative word.

    Elements are normalized to Cartier--Foata normal form upon construction ::

        sage: FC = FullyCommutativeCoxeterElements(['A', 3])
        sage: FC([2, 3, 1, 2])
        [2, 1, 3, 2]

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

        # In the following, we use a dictionary's keys as a replacement for a
        # set. this is because dictinary keys are guaranteed be ordered in
        # Python 3.7+, so this is an easy way to make this iterator
        # deterministic.

        recent_words = {empty_word: True}
        yield empty_word

        length = 1
        while True:
            new_words = {}
            for w in recent_words.keys():
                for s in letters:
                    if w.still_reduced_fc_after_prepending(s):
                        sw = self.element_class(self, [s] + list(w), check=False)
                        # "Add" sw to the "set"
                        new_words[sw] = True

            if len(new_words) == 0:
                return
            for w in new_words.keys():
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
