# -*- coding: utf-8 -*-
r"""
Morphic words

This modules implements morphic words (letter-to-letter coding of fixed
point of a morphism).

AUTHORS:

- Jana Lepsova (January 2021): initial version

EXAMPLES:

Creation of a morphism::

    sage: m = WordMorphism('a->abc,b->baba,c->ca')

"""

from sage.combinat.words.word_infinite_datatypes import WordDatatype_callable
from sage.rings.all import Infinity
from sage.modules.free_module_element import vector

class WordDatatype_morphic(WordDatatype_callable):
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

            sage: m = WordMorphism('a->ab,b->a')
            sage: w = m.fixed_point('a')
            sage: w
            word: abaababaabaababaababaabaababaabaababaaba...
            sage: w[555:1000]
            word: abaababaabaababaababaabaababaabaababaaba...
        
        ::

            sage: m = WordMorphism('a->abc,b->baba,c->ca')
            sage: m.fixed_point('a')
            word: abcbabacababaabcbabaabccaabcbabaabcbabaa...
            sage: w = m.fixed_point('a')
            sage: w[7]
            'c'
            sage: w[2:7]
            word: cbaba
            sage: w[500:503]
            word: caa

            
        TESTS::

            sage: from sage.combinat.words.word_infinite_datatypes import WordDatatype_callable
            sage: WordDatatype_callable(Words(),lambda n:n%3)
            <sage.combinat.words.word_infinite_datatypes.WordDatatype_callable object at ...>
            sage: WordDatatype_callable(Words([0,1,2]),lambda n:n%3)
            <sage.combinat.words.word_infinite_datatypes.WordDatatype_callable object at ...>
        """
        self._len = Infinity
        self._parent = parent
        # self._func = callable
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
        r"""
        Return the representation of the integer n in the numeration system
        associated to the morphism.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT:

        list
        
        EXAMPLES::
        
            sage: m = WordMorphism('a->ab,b->a')
            sage: w = m.fixed_point('a')
            sage: w.representation(5)
            [1, 0, 0, 0]

        TESTS:

        Accessing this method from an instance of the current class (no using
        the inherited word classes)::
            
            sage: m = WordMorphism('a->ab,b->a')
            sage: W = m.domain()
            sage: from sage.combinat.words.morphic import WordDatatype_morphic
            sage: w = WordDatatype_morphic(W, m, 'a')
            sage: type(w)
            <class 'sage.combinat.words.morphic.WordDatatype_morphic'>           
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

    def _func(self, key):
        """
        Returns a letter of a fixed point of a morphism on position ``key``.

        INPUT:

        - ``self`` - a fixed point of a morphism
        - ``key`` - an integer, the position

        OUTPUT:

        - a letter

        EXAMPLES::
        
            sage: m = WordMorphism("a->ab,b->a")
            sage: w = m.fixed_point("a")
            sage: w[0]
            'a'
            sage: w[5]
            'a'
            sage: w[10000]
            'a'

        TESTS:

        Accessing this method from an instance of the current class (no using
        the inherited word classes)::
            
            sage: m = WordMorphism('a->ab,b->a')
            sage: W = m.domain()
            sage: from sage.combinat.words.morphic import WordDatatype_morphic
            sage: w = WordDatatype_morphic(W, m, 'a')
            sage: w._func(5)
            'a'
       
        """
        letter = self._letter
        for a in self.representation(key):
            letter = (self._morphism(letter))[a]
        if key == 0: 
            return self._coding[letter]
        return self._coding[letter]

    def __iter__(self):
        r"""
        Returns an iterator of the letters of the fixed point of ``self``
        starting with ``letter``.

        If w is the iterated word, then this iterator: outputs the elements
        of morphism[ w[i] ], appends morphism[ w[i+1] ] to w, increments i.

        INPUT:

        - ``self`` - an endomorphism, must be prolongable on
           letter

        - ``letter`` - a letter in the domain of ``self``

        OUTPUT:

        - iterator of the fixed point

        EXAMPLES::

            sage: m = WordMorphism("a->ab,b->a")
            sage: w = m.fixed_point("a")
            sage: it = iter(w)
            sage: [next(it) for _ in range(10)]
            ['a', 'b', 'a', 'a', 'b', 'a', 'b', 'a', 'a', 'b']

        Works with erasing morphisms::

            sage: m = WordMorphism('a->abc,b->,c->')
            sage: w = m.fixed_point("a")
            sage: list(w)
            ['a', 'b', 'c']

        The morphism must be prolongable on the letter or the iterator will
        be empty::

            sage: list(m.fixed_point("b"))
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on b

        The morphism must be an endomorphism::

            sage: m = WordMorphism('a->ac,b->aac')
            sage: w = m.fixed_point('a')
            Traceback (most recent call last):
            ...
            TypeError: self (=a->ac, b->aac) is not an endomorphism

        We check that :trac:`8595` is fixed::

            sage: s = WordMorphism({('a', 1):[('a', 1), ('a', 2)], ('a', 2):[('a', 1)]})
            sage: w = s.fixed_point(('a', 1))
            sage: it = iter(w)
            sage: next(it)
            ('a', 1)

        This shows that ticket :trac:`13668` has been resolved::

            sage: s = WordMorphism({1:[1,2],2:[2,3],3:[4],4:[5],5:[6],6:[7],7:[8],8:[9],9:[10],10:[1]})
            sage: (s^7).fixed_points()
            [word: 1223234234523456234567234567823456789234...,
             word: 2,3,4,5,6,7,8,9,10,1,1,2,1,2,2,3,1,2,2,3,2,3,4,1,2,2,3,2,3,4,2,3,4,5,1,2,2,3,2,3,...]
            sage: (s^7).reversal().fixed_points()
            []
        """
        from itertools import chain
        w = iter(self._morphism.image(self._letter))
        while True:
            try:
                for a in self._morphism.image(next(w)):
                    yield self._coding[a]
                else:
                    next_w = next(w)
                    w = chain([next_w], w, self._morphism.image(next_w))
            except StopIteration:
                return


