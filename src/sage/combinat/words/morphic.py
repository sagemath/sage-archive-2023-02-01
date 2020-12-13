from sage.combinat.words.word_datatypes import WordDatatype
from sage.rings.all import Infinity
from sage.modules.free_module_element import vector
import itertools


class WordDatatype_morphic(WordDatatype):
    r"""
    Datatype for a word defined by a callable.
    """
    def __init__(self, parent, morphism, letter, coding=None):
        r"""
        INPUT:

        - ``parent`` - a parent
        - ``morphism`` - a word morphism
        - ``letter`` - a starting letter
        - ``coding`` - dict (default: ``None``), if ``None`` 
          the identity map is used for the coding

        EXAMPLES::

            sage: f = lambda n : 'x' if n % 2 == 0 else 'y'
            sage: w = Word(f, length=9, caching=False); w
            word: xyxyxyxyx
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_callable'>
            sage: w.length()
            9

        ::

            sage: w = Word(f, caching=False); w
            word: xyxyxyxyxyxyxyxyxyxyxyxyxyxyxyxyxyxyxyxy...
            sage: type(w)
            <class 'sage.combinat.words.word.InfiniteWord_callable'>
            sage: w.length() is None
            False
            sage: w.length()
            +Infinity

        TESTS::

            sage: from sage.combinat.words.word_infinite_datatypes import WordDatatype_callable
            sage: WordDatatype_callable(Words(),lambda n:n%3)
            <sage.combinat.words.word_infinite_datatypes.WordDatatype_callable object at ...>
            sage: WordDatatype_callable(Words([0,1,2]),lambda n:n%3)
            <sage.combinat.words.word_infinite_datatypes.WordDatatype_callable object at ...>
        """
        self._len = Infinity
        self._parent = parent
        # for hashing
        self._hash = None
        
        self._morphism = morphism
        self._letter = letter
        self._alphabet = self._morphism.domain().alphabet()
        if coding == None:
            self._coding = {a:a for a in self._alphabet}
        else:
            self._coding = coding
            

    def representation(self, n):
        """
        EXAMPLES::
        
            sage: from sage.combinat.words.morphic import WordDatatype_morphic
            sage: m = WordMorphism('a->ab,b->a')
            sage: W = m.domain()
            sage: w = WordDatatype_morphic(W,m,'a')
            sage: w.representation(5)
            [1, 0, 0, 0]
        """
        letters_to_int =  {a:i for (i,a) in enumerate(self._alphabet)}
        position = letters_to_int[self._letter]
        M = self._morphism.incidence_matrix()
        vMk = vector([1]*len(self._alphabet))
        length_of_images = [vMk]
        while vMk[position] <= n:
            vMk = vMk*M
            length_of_images.append(vMk)
        length_of_images.pop()
        k = len(length_of_images)  
        letter_k = self._letter
        n_k = n
        path = []
        while k > 0:
            m_letter_k = self._morphism(letter_k)
            S = 0
            j = 0
            while S <= n_k:
                a = m_letter_k[j]
                i = letters_to_int[a]
                pile_length = length_of_images[k-1][i]
                S += pile_length
                j += 1
            path.append(j-1)
            n_k -= S - pile_length
            letter_k = a
            k -= 1
        return path
    def __getitem__(self, key):
        letter = self._letter
        for a in self.representation(key):
            letter = (self._morphism(letter))[a]
        if key == 0: 
            return self._coding[letter]
        return self._coding[letter]
