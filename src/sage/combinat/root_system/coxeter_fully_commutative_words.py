from ..words.words import FiniteWords
from ..words.word import FiniteWord_list
from .coxeter_matrix import CoxeterMatrix
from .cartan_type import CartanType

from sage.misc.lazy_attribute import lazy_attribute

def FullyCommutativeReducedCoxeterWords(data):
    if isinstance(data, CoxeterMatrix):
        return FullyCommutativeReducedCoxeterWords_class(data)
    else:
        try:
            t = CartanType(data)
        except (TypeError, ValueError):
            raise ValueError('Input must be a CoxeterMatrix or data for a Cartan type.')
        return FullyCommutativeReducedCoxeterWords_class(t.coxeter_matrix())

class FullyCommutativeReducedCoxeterWords_class(FiniteWords):
    def __init__(self, coxeter_matrix):
        self._matrix = coxeter_matrix
        FiniteWords.__init__(self, sorted(self._matrix.index_set()))

    @lazy_attribute
    def _element_classes(self):
        return {'list': FullyCommutativeReducedCoxeterWord_list}

    def coxeter_matrix(self):
        return self._matrix

    def __iter__(self):
        m = self.coxeter_matrix()

        empty_word = self._element_classes['list'](self, [])
        letters = self.coxeter_matrix().index_set()

        recent_words = {empty_word}
        yield empty_word

        length = 1
        while True:
            new_words = set()
            for w in recent_words:
                for s in letters:
                    if w.still_reduced_fc_after_prepending(s):
                        sw = self._element_classes['list'](self, [s] + list(w)).cartier_foata()
                        new_words.add(sw)

            if len(new_words) == 0:
                return
            for w in new_words:
                yield w
            recent_words = new_words
            length += 1

    # Operations on or between letters and words are implemented on the parent,
    # all taking letters as integers and words as lists.
    # This is for performance reasons; methods on the elements are provided
    # that delegate to these functions on the parent.

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
            for s in self.alphabet():
                i = self._find_left_descent(s, cur_word)
                if i is not None:
                    out_word.append(s)
                    # cur_word = cur_word[:i] + cur_word[i+1:]
                    cur_word.pop(i)
            
        return self._element_classes['list'](self, out_word)

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

class FullyCommutativeReducedCoxeterWord_list(FiniteWord_list):
    def _repr_(self):
        return 'FC ' + super()._repr_()

    # Provide public methods on the elements that wrap the private element
    # methods in the parent.

    def find_left_descent(self, s):
        return self.parent()._find_left_descent(s, list(self))

    def has_left_descent(self, s):
        return self.find_left_descent(s) is not None

    def left_descents(self):
        return self.parent()._left_descents(list(self))

    def cartier_foata(self):
        return self.parent()._cartier_foata(list(self))

    def still_reduced_fc_after_prepending(self, s):
        return self.parent()._still_reduced_fc_after_prepending(s, list(self))
