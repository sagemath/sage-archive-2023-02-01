from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import NormalizedClonableList, ClonableArray
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets 
from ..words.words import FiniteWords
from ..words.word import FiniteWord_list
from .coxeter_matrix import CoxeterMatrix
from .cartan_type import CartanType
from collections import deque
from sage.rings.all import Infinity

def FullyCommutativeReducedCoxeterWords(data):
    if isinstance(data, CoxeterMatrix):
        return FullyCommutativeReducedCoxeterWords_class(data)
    else:
        try:
            t = CartanType(data)
        except (TypeError, ValueError):
            raise ValueError('Input must be a CoxeterMatrix or data for a Cartan type.')
        return FullyCommutativeReducedCoxeterWords_class(t.coxeter_matrix())

class FullyCommutativeReducedCoxeterWord(NormalizedClonableList):
    def check(self):
        if not self.parent()._is_fully_commutative(self._get_list()):
            raise ValueError('list does not represent a fully commutative word.')

    def normalize(self):
        self._require_mutable()
        normalized = self.parent()._cartier_foata(self._get_list())
        self._set_list(normalized)

    # Provide public methods on the elements that wrap the private element
    # methods in the parent.

    def find_left_descent(self, s):
        return self.parent()._find_left_descent(s, list(self))

    def has_left_descent(self, s):
        return self.find_left_descent(s) is not None

    def left_descents(self):
        return self.parent()._left_descents(list(self))

    def still_reduced_fc_after_prepending(self, s):
        return self.parent()._still_reduced_fc_after_prepending(s, list(self))


class FullyCommutativeReducedCoxeterWords_class(Parent):
    def __init__(self, coxeter_matrix):
        self._matrix = coxeter_matrix

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

    Element = FullyCommutativeReducedCoxeterWord

    def coxeter_matrix(self):
        return self._matrix

    def __iter__(self):
        m = self.coxeter_matrix()

        empty_word = self.element_class(self, [], check=False)
        letters = self.coxeter_matrix().index_set()

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

    def _find_left_descent(self, s, w):
        m = self.coxeter_matrix()
        for (i, t) in enumerate(w):
            if t == s and not any(m[x, t] > 2 for x in w[:i]):
                return i
        return None

    def _left_descents(self, w):
        m = self.coxeter_matrix()
        out = set()
        for (i, t) in enumerate(w):
            if not any(m[x, t] > 2 for x in w[:i]):
                out.add(t)
        return out

    def _cartier_foata(self, w):
        cur_word = list(w)
        out_word = []
        
        while len(cur_word) > 0:
            for s in self.coxeter_matrix().index_set():
                i = self._find_left_descent(s, cur_word)
                if i is not None:
                    out_word.append(s)
                    # cur_word = cur_word[:i] + cur_word[i+1:]
                    cur_word.pop(i)
            
        return out_word

    def _still_reduced_fc_after_prepending(self, s, w):
        m = self.coxeter_matrix()
        if self._find_left_descent(s, w) is not None:
            return False

        # Find the first letter in that doesn't commute with s.
        try:
            (j, t) = next((i, x) for (i, x) in enumerate(w) if m[s, x] >= 3)
        except StopIteration:
            return True

        u = w[j:]
        for c in range(m[s, t] - 1):
            letter = t if c % 2 == 0 else s
            i = self._find_left_descent(letter, u)
            if i is not None:
                # Remove letter
                # u = u[:i] + u[i+1:]
                # u = u.with_index_removed(i)
                u.pop(i)
            else:
                return True

        return False
