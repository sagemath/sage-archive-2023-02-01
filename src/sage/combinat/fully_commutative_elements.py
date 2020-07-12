r"""
Fully commutative elements of Coxeter groups

An element `w` in a Coxeter system (W,S) is fully commutative (FC) if
every two reduced word of w can be related by a sequence of only
commutation relations, i.e., relations of the form `st=ts` where `s,t` are
commuting generators in `S`. See [Ste1996]_.
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

class FullyCommutativeElement(NormalizedClonableList):
    r"""
    A (reduced word of a) fully commutative (FC) element in a Coxeter system.

    An element `w` in a Coxeter system (W,S) is fully commutative (FC) if every
    two reduced word of w can be related by a sequence of only commutation
    relations, i.e., relations of the form `st=ts` where `s,t` are commuting
    generators in `S`. Equivalently, `w` is FC if and only if for every pair of
    generators `s,t \in S` for which `m(s,t)>2`, no reduced word of `w` contains
    the "braid" word `sts...` of length `m(s,t)` as a contiguous subword. See
    [Ste1996]_. We will use the braid-avoidance criterion to check if an element
    is FC.

    Every FC element has a canonical reduced word called its Cartier--Foata
    form. See [Gre2006]_. We will normalize each FC element to this form. 

    """
  
    # Methods required as a subclass of NormalizedClonableList:
    def check(self):
        r"""
        Called automatically when an element is created. Alias of
        :func:`is_fully_commutative`
        """
        return self.is_fully_commutative()

    def normalize(self):
        r"""
        Called automatically when an element is created. Alias of
        :func:`cartier_foata_form`
        """
        return self.cartier_foata_form()


    def is_fully_commutative(self):
        r"""
        Check if ``self`` is the reduced word of an FC element.

        :func:`check` is an alias of this method, and is called automatically
        when an element is created.

        EXAMPLES:

        To construct an FC element, first call the parent class
        FullyCommutativeElements. The parent class contains information about
        the Coxeter matrix of the ambient Coxeter system :: 

            sage: FC = FullyCommutativeElements(['B', 3])
            sage: FC.coxeter_matrix()
            [1 3 2]
            [3 1 4]
            [2 4 1]

        We can construct FC elements as follows ::

            sage: FC([])
            []
            sage: FC([1,2])
            [1, 2] 
            sage: FC([2,3,2])
            [2, 3, 2]
            sage: FC([3,2,3])
            [3, 2, 3]

        The output is the normalized form of ``self``, which may be a different
        reduced word of the element represented by the input ::

            sage: FC([3,1])
            [1, 3]
            sage: FC([2,3,1])
            [2, 1, 3] 
            sage: FC([1,3]) == FC([3,1])
            True

        If the input is not the reduced word of an FC element, return a
        ValueEror ::

            sage: FC([1,2,1])
            Traceback (most recent call last):
            ...
            ValueError: The input is not a reduced word of a fully commutative element. 

            sage: FC([2,3,2,3])
            Traceback (most recent call last):
            ...
            ValueError: The input is not a reduced word of a fully commutative element.
        """
        matrix = self.parent().coxeter_matrix()
        w = tuple(self)

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

        not_fc = ValueError('The input is not a reduced word of a fully commutative element.')
        if contains_long_braid(w):
            raise not_fc
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
                                raise not_fc
                            else:
                                checked.add(new_word)
                                queue.appendleft(new_word)
            return True


    def cartier_foata_form(self):
        r"""
        Return the Cartier--Foata form of ``self``.

        :func:`normalize` is an alias of this method, and is called
        automatically when an element is created.

        EXAMPLES:

        The following reduced words express the same FC elements in `B_5` ::

            sage: FC = FullyCommutativeElements(['B', 5])
            sage: FC([1, 4, 3, 5, 2, 4, 3])
            [1, 4, 3, 5, 2, 4, 3]
            sage: FC([4, 1, 3, 5, 2, 4, 3])
            [1, 4, 3, 5, 2, 4, 3]
            sage: FC([4, 3, 1, 5, 4, 2, 3])
            [1, 4, 3, 5, 2, 4, 3]

        .. NOTE::

            The Cartier--Foata form of a reduced word of an FC element w
            can be found recursively by repeatedly moving left descents of
            elements to the left and ordering the left descents from small to
            large. In the above example, the left descents of the element are 4
            and 1, therefore the Cartier--Foata form of the element is the
            concatenation of [1,4] with the Cartier--Foata form of the remaining
            part of the word. See [Gre2006]_.

        .. SEEALSO:: :func:`descents`
        """
        self._require_mutable()

        out_word = []
        
        while len(self) > 0:
            fronts = self.descents()
            out_word.extend(sorted(fronts))
            for s in fronts:
                self.remove(s)

        self._set_list(out_word)


    ########## Descents ########## 

    # The following three functions deal with descents of FC elements.
    # Finding descents for FC elements is easier than doing so for general
    # elements but also extremely useful: repeatedly searching for descents of
    # FC elements can allow us to compute the Cartier--Foata forms and certain
    # coset decompositions of FC elements. The coset decompositions will
    # in turn be useful for computing so-called star operations on FC elements
    # and for generating FC elements by increasing length.

    def find_descent(self, s, side='left'):
        r"""
        Return the index of a left descent ``s`` of ``self``.

        INPUT:

        - ``s`` -- integer representing a generator of the Coxeter system.

        OUTPUT: 

        Determine if the generator ``s`` is a left descent of ``self``; return
        the index of the leftmost occurrence of ``s`` in ``self`` if so and
        return ``None`` if not.  

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default: 'left'); if the argument is set to
          'right', the function checks if ``s`` is a right descent of ``self``
          and finds the index of the rightmost occurrence of ``s`` if so.


        EXAMPLES::

            sage: FC = FullyCommutativeElements(['B', 5])
            sage: w = FC([1, 4, 3, 5, 2, 4, 3])
            sage: w.find_descent(1)
            0
            sage: w.find_descent(1, side='right')
            <BLANKLINE>
            sage: w.find_descent(4)
            1
            sage: w.find_descent(4, side='right')
            <BLANKLINE>
            sage: w.find_descent(3)
            <BLANKLINE>

        .. NOTE:: 

            A generator `s` is a left descent of an FC element `w` if and
            only if for one (equivalently, every) reduced word of `w`, `s`
            appears to in the word and every generator to the left of the
            leftmost `s` in the word commutes with `s`. A similar result holds
            for right descents of FC elements. 

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

        OUTPUT: a boolean value

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default: 'left'); if set to 'right', determine if
          ``self`` has ``s`` as a right descent.

        EXAMPLES::

            sage: FC = FullyCommutativeElements(['B', 5])
            sage: w = FC([1, 4, 3, 5, 2, 4, 3])
            sage: w.has_descent(1)
            True
            sage: w.has_descent(1, side='right')
            False
            sage: w.has_descent(4)
            True
            sage: w.has_descent(4, side='right')
            False

        .. SEEALSO:: :func:`find_descent`
        """
        return self.find_descent(s, side=side) is not None

    def descents(self, side='left'):
        r"""
        Obtain the set of descents on the appropriate side of ``self``.

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default: 'left'); if set to 'right', find the
          right descents.

        EXAMPLES::


            sage: FC = FullyCommutativeElements(['B', 5])
            sage: w = FC([1, 4, 3, 5, 2, 4, 3])
            sage: sorted(w.descents())
            [1, 4]
            sage: w.descents(side='right')
            {3}
        
        .. NOTE::
            See the chacterization of descents for FC elements in
            :func:`find_descent`.

        """
        view = list(self) if side == 'left' else self[::-1]
        m = self.parent().coxeter_matrix()
        out = set()
        for (i, t) in enumerate(view):
            if not any(m[x, t] > 2 for x in view[:i]):
                out.add(t)
        return out

    # coset decomposition for FC elements
    def coset_decomposition(self, J, side='left'):
        r"""
        Return the coset decomposition of ``self`` with repsect to the parabolic
        subgroup generated by ``J``.

        INPUT: 

        - ``J`` -- subset of the generating set `S` of the Coxeter system.

        OUTPUT:

        The tuple of elements `(w_J, w^J)` such that `w=w_J \cdot w^J`, `w_J` is
        generated by the elements in `J`, and `w^J` has no left descent from
        `J`. This tuple is unique and satisfies the equation `l(w) = l(w_J) +
        l(w^J)` where `l` denotes Coxeter length by general theory; see
        Proposition 2.4.4 of [BB2005]_. 

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default: 'left'); if the value is set to 'right',
          then the function returns the tuple `(w'^J, w'_J)` from the coset
          decomposition `w = w'^J \cdot w'_J` of `w` with respect to `J`.  

        EXAMPLES::

            sage: FC = FullyCommutativeElements(['B', 6])
            sage: w = FC([1, 6, 2, 5, 4, 6, 5])
            sage: w.coset_decomposition({1})
            ([1], [6, 2, 5, 4, 6, 5])
            sage: w.coset_decomposition({1}, side = 'right')
            ([1, 6, 2, 5, 4, 6, 5], [])
            sage: w.coset_decomposition({5, 6})
            ([6, 5, 6], [1, 2, 4, 5])
            sage: w.coset_decomposition({5, 6}, side='right')
            ([1, 6, 2, 5, 4], [6, 5])

        .. NOTE::

            The factor `w_J` of the coset decomposition `w = w_J \cdot
            w^J` can be obtained by greedily "pulling left descents of `w` that
            are in `J` to the left"; see the proof of [BB2005]_. This greedy
            algorithm works for all elements in Coxeter group, but it becomes
            especially simple for FC elements because descents are easier to
            find for FC elements. 

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

    # First application of coset decompositions: 
    # help enumerate FC elements of a Coxeter system by induction on length: 

    def still_reduced_fc_after_prepending(self, s):
        r"""
        Determine if ``self`` prepended with ``s`` is still a reduced word of an
        FC element in the Coxeter system.

        INPUT:

        - ``s`` -- integer representing a generator of the Coxeter system.
        - ``self`` -- a reduced word of an FC element


        EXAMPLES::

            sage: FCB3 = FullyCommutativeElements(['B', 3])
            sage: FCB3.coxeter_matrix()
            [1 3 2]
            [3 1 4]
            [2 4 1]
            sage: w = FCB3([1,2])
            sage: w.still_reduced_fc_after_prepending(1)
            False
            sage: w.still_reduced_fc_after_prepending(2)
            False
            sage: w.still_reduced_fc_after_prepending(3)
            True
            sage: u = FCB3([3,1,2])            
            sage: u.still_reduced_fc_after_prepending(1)
            False
            sage: u.still_reduced_fc_after_prepending(2)
            True
            sage: u.still_reduced_fc_after_prepending(3)
            False

            sage: FCA5 = FullyCommutativeElements(['A',5])
            sage: w = FCA5([2,4,1,3,2,5])
            sage: w.still_reduced_fc_after_prepending(5)
            False


        .. NOTE:: 

            If `w` is a reduced word of an element, then the concatenation
            `sw` is still a reduced word if and only if `s` is not a left
            descent of `w` by general Coxeter group theory. So now assume `w` is
            a reduced word of an FC element and `s` is not a left descent `w`.
            In this case, Lemma 4.1 of [Ste1996]_ implies that `sw` is not a
            reduced word of an FC element if and only if some letter in `w` does
            not commute with `s` and the following conditions hold
            simultaneously for the leftmost such letter `t`:


            (1) `t` is left descent of the word `u_1` obtained by removing
            all letters to the left of the aforementioned `t` from `w`;
            (this condition is automatically true)

            (2) `s` is left descent of the word `u_2`  obtained by
            removing the leftmost `t` from `u_1`;

            (3) `t` is left descent of the word `u_3`  obtained by
            removing the leftmost `s` from `u_2`;
            ...
            (m-1) the appropriate element in `\{s, t\}` is a left descent
            of the word `u_{m-1}` obtained by removing the leftmost letter
            required to be a descent in Condition (m-2) from `u_{m-2}`.

            In the last example above, we have `s=5`, `t=4`, Condition (1)
            holds, but Condition (2) fails, therefore `5w` is still a
            reduced word of an FC element.

            Note that the conditions (1)--(m-1) are equivalent to the
            condition that the parabolic factor `u_J` from the coset
            decomposition `u_1 = u_J \cdot u^J` of `u_1` with respect to
            `J := \{s, t\}` is the element `tst...` of length `m(s,t)-1`.


        REFERENCES:

        See Lemma 4.1 of  [Ste1996]_.
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
        x, y= u.coset_decomposition({s, t})
        return len(x) != m[s,t]-1


    # Second application of coset decompositions: star operations.

    # Star operations were first defined on elements of Coxeter groups by
    # Kazhdan and Lusztig in [KL1979]_ with respect to pair of generators `s,t`
    # such that `m(s,t)=3`. Later, Lusztig generalized the definition in
    # [Lus1985]_, via coset decompositions, to allow star operations with
    # respect to any pair of generators `s,t` such that `m(s,t)\ge 3`. Given
    # such a pair, we can potentially perform four types of star operations:
    # left upper, left lower, right upper and right lower; see [Gre2006]_. 

    # Template for all four types of star operations: 

    def _star_operation_inner(self, J, direction, side):
        assert len(J) == 2, 'J needs to contain a pair of generators.' 
        s, t = J
        mst = self.parent().coxeter_matrix()[s,t]
        # Perform the coset decomposition on the specified side:
        if side == 'left':
            (string, remaining) = self.coset_decomposition(J, side=side)
        elif side == 'right':
            (remaining, string) = self.coset_decomposition(J, side=side)

        cur_string = list(string)

        # From the coset decomposition, perform the upper or lower operation:  

        if direction == 'down' and 2 <= len(string) <= mst - 1:
            # The lower star operation is defined if the parabolic factor `w_J`
            # from the coset decomposition has length at least 2; when defined, it
            # removes the outmost letter in `w_J`:
            # decrease the length of the J-string:
            new_string = cur_string[1:] if side == 'left' else cur_string[:-1]
        elif direction == 'up' and 1 <= len(string) <= mst - 2:
            # The upper star operation is defined if the parabolic factor `w_J`
            # from the coset decomposition has length at most `m(s,t)-2`; when
            # defined, it adds the appropriate letter `x` from `J` to `w_J`,
            # extending its length by 1.
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

        INPUT:

        - ``J`` -- a set of two integers representing two non-commuting
          generators of the Coxeter system. 

        OUTPUT: 

        The result of the star operation if it is defined on ``self``, ``None``
        otherwise.

        OPTIONAL ARGUEMNTS:

        - ``side`` -- string (default: 'left'); if set to 'right', perform a
          right upper star operation.


        EXAMPLES::

            sage: FC = FullyCommutativeElements(['B', 6])
            sage: w = FC([1, 6, 2, 5, 4, 6, 5])
            sage: w.coset_decomposition({5, 6})
            ([6, 5, 6], [1, 2, 4, 5])
            sage: w.upper_star({5,6})
            <BLANKLINE>
            sage: w.coset_decomposition({5, 6}, side='right')
            ([1, 6, 2, 5, 4], [6, 5])
            sage: w.upper_star((5, 6), side='right')
            [1, 6, 2, 5, 4, 6, 5, 6]
            sage: v = FC([1, 6, 2, 5, 4, 6])
            sage: v.coset_decomposition({5, 6}, side='right')
            ([1, 6, 2, 5, 4], [6])
            sage: v.upper_star({5, 6}, side='right')
            [1, 6, 2, 5, 4, 6, 5]

        .. NOTE:: 

            As the examples illustrate, the left upper star operation is
            defined if and only if in the coset decomposition `w = w_J \cdot
            {}^J w`, the parabolic part has length `1\le l(w_J) \le m(s,t)-2`.
            When this is the case, the operation returns `x\cdot w_J\cdot w^J`
            where `x` is the letter `J` different from the leftmost letter of
            `w_J`. Similar facts hold for right upper star operations. See
            [Gre2006]_.
        """
        return self._star_operation_inner(J, 'up', side)

    def lower_star(self, J, side='left'):
        r"""
        Perform a lower star operation on ``self`` with respect to the parabolic
        subgroup generated by ``J``.

        INPUT:

        - ``J`` -- a set of two integers representing two non-commuting
          generators of the Coxeter system. 

        OUTPUT: 

        The result of the star operation if it is defined on ``self``, ``None``
        otherwise.

        OPTIONAL ARGUEMNTS:

        - ``side`` -- string (default: 'left'); if set to 'right', perform a
          right lower star operation.


        EXAMPLES::

            sage: FC = FullyCommutativeElements(['B', 6])
            sage: w = FC([1, 6, 2, 5, 4, 6, 5])
            sage: w.coset_decomposition({5, 6})
            ([6, 5, 6], [1, 2, 4, 5])
            sage: w.lower_star({5,6})
            [1, 5, 2, 4, 6, 5]
            sage: w.lower_star({5,6}).coset_decomposition({5, 6})
            ([5, 6], [1, 2, 4, 5])
            sage: w.coset_decomposition({5, 6}, side='right')
            ([1, 6, 2, 5, 4], [6, 5])
            sage: w.lower_star((5, 6), side='right')
            [1, 6, 2, 5, 4, 6]
            sage: v = FC([1, 6, 2, 5, 4, 6])
            sage: v.coset_decomposition({5, 6}, side='right')
            ([1, 6, 2, 5, 4], [6])
            sage: v.lower_star({5, 6}, side='right')
            <BLANKLINE>

        .. NOTE::
        
            As the examples illustrate, the left lower star operation is
            defined if and only if in the coset decomposition `w = w_J \cdot
            {}^J w`, the parabolic part has length `2\le l(w_J) \le m(s,t)-1`.
            When this is the case, the operation removes the leftmost letter of
            `w_J` from `w`.  Similar facts hold for right upper star operations.
            See [Gre2006]_.

        """
        return self._star_operation_inner(J, 'down', side)
    

    # Star operations are important to the study of Kazhdan--Lusztig cells in
    # Coxeter groups; see, for example, [KL1979]_ and [GX2020]_. The orbits of
    # suitable FC elements under star operations are often certain
    # Kazhdan--Lusztig cells. For example, ...


    def star_orbit(self, side='left', **kargs):
        r"""
        Compute the star operation orbit of ``self``.

        OUTPUT: The set containing all elements that can be obtained from
        ``self`` via a sequence of star operation of the specified type.  

        OPTIONAL ARGUMENTS:

        - ``side`` -- string (default: 'left'); if set to 'right', the function
          compute the orbit for the specified type of right star operations. 
        - ``upper_only`` -- boolean (default: False); if passed, compute only
          the set of elements that can be obtained from ``self`` via upper star
          operations on the specified side.
        - ``lower_only`` -- boolean (default: False); if passed, compute only
          the set of elements that can be obtained from ``self`` via lower star
          operations on the specified side.


        EXAMPLES:

        Compute the left star closure of [1] in the group `I_8`. This should be
        set of `\{1,2\}`-braids of lengths 1 through 7, and should be the same
        as the left *upper* star closure ::

            sage: FC = FullyCommutativeElements(['I', 8])
            sage: sorted(FC([1]).star_orbit())
            [[1],
             [1, 2, 1],
             [1, 2, 1, 2, 1],
             [1, 2, 1, 2, 1, 2, 1],
             [2, 1],
             [2, 1, 2, 1],
             [2, 1, 2, 1, 2, 1]]
            sage: FC([1]).star_orbit() == FC([1]).star_orbit(upper_only=True)
            True
        """
        m = self.parent().coxeter_matrix()
        adjacent_pairs = [(a, b) for (a, b) in itertools.product(self.parent().index_set(), repeat=2) if a < b and m[a,b] > 2]
        
        directions = {'up', 'down'}
        if 'upper_only' in kargs and kargs['upper_only']:
            directions = {'up'}
        elif 'lower_only' in kargs and kargs['lower_only']:
            directions = {'down'}

        closure = {self}
        recent_words = {self}
        while True:
            new_words = set()
            for w in recent_words:
                for J in adjacent_pairs:
                    for d in directions:
                        n = w._star_operation_inner(J, d, side)
                        if n is not None and n not in closure:
                            new_words.add(n)
            if len(new_words) == 0:
                break
            closure.update(new_words)
            recent_words = new_words
        return closure


    ########## Heaps ########## 
    # Heaps are certain posets that are very useful for representing and
    # studying FC elements; see, for example, [Ste1996]_ and [GX2020]_.

    def heap(self, **kargs):
        r"""
        Create the heap poset of ``self``.

        The heap of an FC element `w` is a labeled poset that can be defined
        from any reduced word of `w`. Different reduced words yield isomorphic
        labeled posets, so the heap is well defined. 

        INPUT:

        - ``self`` -- list, a reduced word `w=s_0... s_{k-1}` of an FC element.

        OUTPUT: A labeled poset where the underlying set is `\{0,1,...,k-1\}`
        and where each element `i` carries `s_i` as its label.

        OPTIONAL ARGUMENTS:

        - ``one_index`` -- boolean (default: False). Setting the value to True
          will change the underlying set of the poset to `\{1, 2, \dots, n\}`.

        - ``display_labeling`` -- boolean (default: False). Setting the value to
          True will display the label `s_i` for each element `i` of the poset.

        EXAMPLES::

            sage: FC = FullyCommutativeElements(['A', 5])
            sage: FC([1, 4, 3, 5, 2, 4]).heap().cover_relations()
            [[1, 2], [1, 3], [2, 5], [2, 4], [3, 5], [0, 4]]
            sage: FC([1, 4, 3, 5, 4, 2]).heap(one_index=True).cover_relations()
            [[2, 3], [2, 4], [3, 6], [3, 5], [4, 6], [1, 5]]

        .. NOTE::

            The partial order in the heap is defined by declaring `i\prec
            j` if `i<j` and `m(s_i,s_j)\neq 2`. See [Ste1996]_.
        """
        m = self.parent().coxeter_matrix()

        one_index = kargs.get('one_index', False)
        display_labeling = kargs.get('display_labeling', False)
        # elements of the poset: 
        elements = list(range(1, len(self) + 1)) if one_index else list(range(len(self))) 
        # get the label of each poset element:
        def letter(index):
            return self[index-1] if one_index else self[index]
        # specify the partial order:
        relations = [(i, j) for [i, j] in itertools.product(elements, repeat=2) if i < j and m[letter(i), letter(j)] != 2]
        p = Poset((elements, relations))
        if not display_labeling:
            return p
        else:
            return p.relabel(lambda i: (i, letter(i)))

    def n_value(self):
        r"""
        Calculate the n-value of ``self``.

        The *n-value* of a fully commutative element is the *width* (length
        of any longest antichain) of its heap. The n-value is important as it
        coincides with Lusztig's a-value for FC elements in all Weyl and
        affine Weyl groups as well as so-called star-reducible groups; see
        [GX2020]_. 

        EXAMPLES::

            sage: FC = FullyCommutativeElements(['A', 5])
            sage: FC([1,3]).n_value()
            2
            sage: FC([1,2,3]).n_value()
            1
            sage: FC([1,3,2]).n_value()
            2
            sage: FC([1,3,2,5]).n_value()
            3

        """
        return self.heap().width()
    
    def plot(self):
        r"""
        Display the Hasse diagram of the heap of ``self``.

        The Hasse diagram is rendered in the lattice `S \times \mathbb{N}`, with
        every element `i` in the poset drawn as a point labelled by its label
        `s_i`. Every point is placed in the column for its label at a certain
        level. The levels start at 0 and the level k of an element `i` is the
        maximal number `k` such that the heap contains a chain `i_0\prec
        i_1\prec ... \prec i_k` where `i_k=i`. See [Ste1996]_ and [GX2020]_. 


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

class FullyCommutativeElements(Parent):
    r"""
    Class for the set of fully commutative (FC) elements of a Coxeter systems

    Coxeter systems with finitely many FC elements, or *FC-finite* Coxeter
    systems, are classfied by Stembridge in [Ste1996]_. They fall into seven
    families, namely the groups of types `A_n, B_n, D_n, E_n, F_n, H_n` and
    `I_2(m)`. 

    INPUT:

    - ``data`` -- Coxeter matrix or data describing the Cartan type; the
      latter should be formatted in the same way as required by the
      CoxeterGroup constructor from
      :func:`sage.combinat.root_system.coxeter_group.CoxeterGroup`.

    OUTPUT:

    The class of fully commutative elements in the Coxeter group constructed
    from ``data``. This will belong to either the category of infinite
    enumerated sets or finite enumerated sets depending on if the group is
    FC-finite.

    EXAMPLES:

    Enumerate the FC elements in `A_3` in their Cartier--Foata forms ::

        sage: FCA3 = FullyCommutativeElements(['A', 3])
        sage: FCA3.category()
        Category of finite enumerated sets
        sage: FCA3.list()
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

        sage: FCB8 = FullyCommutativeElements(['B', 8])
        sage: len(FCB8) # long time (7 seconds)
        14299

    Iterate through the FC elements of length up to 3 in the non FC-finite
    group affine `A_2` ::

        sage: FCAffineA2 = FullyCommutativeElements(['A', 2, 1])
        sage: FCAffineA2.category()
        Category of infinite enumerated sets
        sage: list(FCAffineA2.iterate_to_length(2))
        [[], [0], [1], [2], [1, 0], [2, 0], [0, 1], [2, 1], [0, 2], [1, 2]]

    Constructing an element that is not fully commutative throws an error ::

        sage: FCA3([1,2,1])
        Traceback (most recent call last):
        ...
        ValueError: The input is not a reduced word of a fully commutative element. 

    Elements are normalized to Cartier--Foata normal form upon construction ::

        sage: FCA3([2, 3, 1, 2])
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

    Element = FullyCommutativeElement

    def coxeter_matrix(self):
        r"""
        Obtain the Coxeter matrix of the associated Coxter system.

        OUTPUT: CoxeterMatrix
        """
        return self._matrix

    def index_set(self):
        r"""
        Obtain the set of the generators / simple reflections of the
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
        # set. Dictinary keys are guaranteed to be ordered in Python 3.7+, so
        # this is an easy way to make this iterator deterministic.

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

