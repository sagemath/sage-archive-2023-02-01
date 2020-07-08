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

        out_word = []
        
        while len(self) > 0:
            for s in self.parent().index_set():
                i = self.find_left_descent(s)
                if i is not None:
                    out_word.append(s)
                    self.pop(i)

        self._set_list(out_word)

    # Provide public methods on the elements that wrap the private element
    # methods in the parent.

    def find_left_descent(self, s):
        m = self.parent().coxeter_matrix()
        for (i, t) in enumerate(self):
            if t == s and not any(m[x, t] > 2 for x in self[:i]):
                return i
        return None

    def has_left_descent(self, s):
        return self.find_left_descent(s) is not None

    def left_descents(self):
        m = self.parent().coxeter_matrix()
        out = set()
        for (i, t) in enumerate(self):
            if not any(m[x, t] > 2 for x in self[:i]):
                out.add(t)
        return out

    def still_reduced_fc_after_prepending(self, s):
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


class FullyCommutativeReducedCoxeterWords_class(Parent):
    def __init__(self, coxeter_matrix):
        self._matrix = coxeter_matrix
        self._index_set = sorted(coxeter_matrix.index_set())

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

    def index_set(self):
        return self._index_set

    def __iter__(self):
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
