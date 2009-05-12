# coding=utf-8
"""
Word generators
"""
#*****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>,
#                          Sébastien Labbé <slabqc@gmail.com>,
#                          Arnaud Bergeron <abergeron@gmail.com>,
#                          Amy Glen <amy.glen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import cycle, count
from random import randint
from sage.structure.sage_object import SageObject
from sage.rings.arith import gcd
from sage.combinat.words.word import FiniteWord_over_OrderedAlphabet, InfiniteWord_over_OrderedAlphabet, is_Word, is_FiniteWord
from sage.combinat.words.words import Words
from sage.combinat.words.alphabet import OrderedAlphabet
from sage.combinat.words.utils import *
from sage.combinat.words.morphism import WordMorphism

def _build_tab(sym, tab, W):
    r"""
    Internal function building a coding table for the phi_inv_tab
    function.

    TESTS::

        sage: from sage.combinat.words.word_generators import _build_tab
        sage: _build_tab(1, [], Words([1, 2]))
        [1]
        sage: _build_tab(1, [1], Words([1, 2]))
        [1, 2]
        sage: _build_tab(2, [1], Words([1, 2]))
        [2, 2]
        sage: _build_tab(2, [1, 2], Words([1, 2]))
        [2, 2, 1]
        sage: _build_tab(1, [2, 2], Words([1, 2]))
        [1, 1, 2]
    """
    res = [sym]
    if len(tab) == 0:
        return res
    if sym == 1:
        res += tab
        res[1] = (res[1] % W.size_of_alphabet()) + 1
        return res
    w = W([sym]).delta_inv(W, tab[0])
    w = w[1:]
    res.append((w[-1] % W.size_of_alphabet()) + 1)
    for i in xrange(1, len(tab)):
        w = w.delta_inv(W, tab[i])
        res.append((w[-1] % W.size_of_alphabet()) + 1)
    return res

class WordGenerator(object):
    r"""
    A class consisting of constructors for several famous words.

    TESTS::

        sage: type(loads(dumps(words)))
        <class 'sage.combinat.words.word_generators.WordGenerator'>
    """

    def ThueMorseWord(self, alphabet=(0, 1)):
        r"""
        Returns the Thue-Morse word over the given two-letter alphabet.

        There are several ways to define the Thue-Morse word `t`.
        We use the following definition: `t[n]` is the sum modulo 2
        of the digits in the binary expansion of `n`.

        INPUT:


        -  ``alphabet`` - any container of length two that is
           suitable to build an instance of OrderedAlphabet (list, tuple, str,
           ...)


        EXAMPLES::

            sage: t = words.ThueMorseWord(); t
            Thue-Morse word on the alphabet [0, 1]
            sage: print t[:50]
            word: 01101001100101101001011001101001100101100110100101

        ::

            sage: t = words.ThueMorseWord('ab'); t
            Thue-Morse word on the alphabet ['a', 'b']
            sage: print t[:50]
            word: abbabaabbaababbabaababbaabbabaabbaababbaabbabaabab

        ::

            sage: t = words.ThueMorseWord(['L1', 'L2']); t
            Thue-Morse word on the alphabet ['L1', 'L2']
            sage: words.ThueMorseWord(['L1','L2'])[:8]
            word: L1,L2,L2,L1,L2,L1,L1,L2

        TESTS::

            sage: words.ThueMorseWord("abc")
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
        """
        if len(set(alphabet)) != 2:
            raise TypeError, "alphabet does not contain two distinct elements"
        from functools import partial
        f = partial(self._ThueMorseWord_nth_digit, alphabet=alphabet)
        w = InfiniteWord_over_OrderedAlphabet(Words(alphabet), f, format= 'function')
        w.rename("Thue-Morse word on the alphabet %s" % w.alphabet().string_rep())
        return w

    def _ThueMorseWord_nth_digit(self, n, alphabet=(0,1)):
        r"""
        Returns the digit in position n of the Thue-Morse word (0 or 1).

        The `n`-th digit of the Thue-Morse word can be defined as
        the number of bits in the 2-complement representation of the
        position modulo 2 which is what this function uses. The running
        time is `O(\log n)` where `n` is the position
        desired.

        INPUT:


        -  ``n`` - integer, the position


        OUTPUT:


        -  ``0 or 1`` - the digit at the position


        TESTS::

            sage: from sage.combinat.words.word_generators import WordGenerator
            sage: WordGenerator()._ThueMorseWord_nth_digit(0)
            0
            sage: WordGenerator()._ThueMorseWord_nth_digit(3)
            0
            sage: WordGenerator()._ThueMorseWord_nth_digit(32)
            1
        """
        for tn in count():
            if n == 0:
                return alphabet[tn & 1]
            n &= n - 1

    def ChristoffelWord(self, p, q, alphabet=(0,1)):
        r"""
        Returns the lower Christoffel word of slope `p/q`, where
        `p` and `q` are relatively prime non-negative
        integers, over the given two-letter alphabet.

        The *Christoffel word of slope `p/q`* is obtained from the
        Cayley graph of `\ZZ/(p+q)\ZZ` with generator
        `q` as follows. If `u \rightarrow v` is an edge in
        the Cayley graph, then `v = u + p \mod{p+q}`. Label the
        edge `u \rightarrow v` by alphabet[1] if `u < v`
        and alphabet[0] otherwise. The Christoffel word is the word
        obtained by reading the edge labels along the cycle beginning from
        0.

        INPUT:


        -  ``alphabet`` - any container of length two that is
           suitable to build an instance of OrderedAlphabet (list, tuple, str,
           ...)


        EXAMPLES::

            sage: w = words.ChristoffelWord(1,0); w
            Lower Christoffel word of slope 1/0 over the alphabet [0, 1]
            sage: print w
            word: 1

        ::

            sage: w = words.ChristoffelWord(0,1); w
            Lower Christoffel word of slope 0/1 over the alphabet [0, 1]
            sage: print w
            word: 0

        ::

            sage: w = words.ChristoffelWord(1,1); w
            Lower Christoffel word of slope 1/1 over the alphabet [0, 1]
            sage: print w
            word: 01

        ::

            sage: w = words.ChristoffelWord(4,7); w
            Lower Christoffel word of slope 4/7 over the alphabet [0, 1]
            sage: print w
            word: 00100100101

        TESTS::

            sage: words.ChristoffelWord(51,43,"abc")
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
        """
        if len(set(alphabet)) != 2:
            raise TypeError, "alphabet does not contain two distinct elements"
        return ChristoffelWord_Lower(p, q, alphabet=alphabet)

    LowerChristoffelWord = ChristoffelWord

    def UpperChristoffelWord(self, p, q, alphabet=(0,1)):
        r"""
        Returns the upper Christoffel word of slope `p/q`, where
        `p` and `q` are relatively prime non-negative
        integers, over the given alphabet.

        The *upper Christoffel word of slope `p/q`* is equal to the
        reversal of the lower Christoffel word of slope `p/q`.
        Equivalently, if `xuy` is the lower Christoffel word of
        slope `p/q`, where `x` and `y` are letters,
        then `yux` is the upper Christoffel word of slope
        `p/q` (because `u` is a palindrome).

        INPUT:


        -  ``alphabet`` - any container of length two that is
           suitable to build an instance of OrderedAlphabet (list, tuple, str,
           ...)


        EXAMPLES::

            sage: w = words.UpperChristoffelWord(1,0); w
            Upper Christoffel word of slope 1/0 over the alphabet [0, 1]
            sage: print w
            word: 1

        ::

            sage: w = words.UpperChristoffelWord(0,1); w
            Upper Christoffel word of slope 0/1 over the alphabet [0, 1]
            sage: print w
            word: 0

        ::

            sage: w = words.UpperChristoffelWord(1,1); w
            Upper Christoffel word of slope 1/1 over the alphabet [0, 1]
            sage: print w
            word: 10

        ::

            sage: w = words.UpperChristoffelWord(4,7); w
            Upper Christoffel word of slope 4/7 over the alphabet [0, 1]
            sage: print w
            word: 10100100100

        TESTS::

            sage: words.UpperChristoffelWord(51,43,"abc")
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
        """
        if len(set(alphabet)) != 2:
            raise TypeError, "alphabet does not contain two distinct elements"
        w = words.LowerChristoffelWord(p, q, alphabet=alphabet).reversal()
        w.rename("Upper Christoffel word of slope %s/%s over the alphabet %s" % (p, q, w.parent().alphabet().string_rep()))
        return w

    def FibonacciWord(self, alphabet=(0, 1), construction_method="recursive"):
        r"""
        Returns the Fibonacci word on the given two-letter alphabet.

        INPUT:


        -  ``alphabet`` - any container of length two that is
           suitable to build an instance of OrderedAlphabet (list, tuple, str,
           ...)

        -  ``construction_method`` - can be any of the
           following: "recursive", "fixed point" (see below for
           definitions).


        Recursive construction: the Fibonacci word is the limit of the
        following sequence of words: `S_0 = 0`,
        `S_1 = 01`, `S_n = S_{n-1} S_{n-2}` for
        `n \geq 2`.

        Fixed point construction: the Fibonacci word is the fixed point of
        the morphism: `0 \mapsto 01` and `1 \mapsto 0`.
        Hence, it can be constructed by the following read-write process:

        beginning at the first letter of `01`, if the next letter
        is `0`, append `01` to the word; if the next letter
        is `1`, append `1` to the word; move to the next
        letter of the word.

        EXAMPLES::

            sage: w = words.FibonacciWord(construction_method="recursive"); w
            Fibonacci word over [0, 1], defined recursively
            sage: print w[:99]
            word: 010010100100101001010010010100100101001010010010100101001001010010010100101001001010010010100101001

        ::

            sage: v = words.FibonacciWord(construction_method="recursive", alphabet='ab'); v
            Fibonacci word over ['a', 'b'], defined recursively
            sage: print v[:99]
            word: abaababaabaababaababaabaababaabaababaababaabaababaababaabaababaabaababaababaabaababaabaababaababaab

        ::

            sage: u = words.FibonacciWord(construction_method="fixed_point"); u
            Fibonacci word over [0, 1], defined as the fixed point of a morphism
            sage: print u[:99]
            word: 010010100100101001010010010100100101001010010010100101001001010010010100101001001010010010100101001

        ::

            sage: x = words.FibonacciWord(construction_method="fixed_point", alphabet=[4, 1]); x
            Fibonacci word over [1, 4], defined as the fixed point of a morphism
            sage: print x[:99]
            word: 414414144144141441414414414144144141441414414414144141441441414414414144141441441414414414144141441

        TESTS::

            sage: from math import floor, sqrt
            sage: golden_ratio = (1 + sqrt(5))/2.0
            sage: a = golden_ratio / (1  + 2*golden_ratio)
            sage: wn = lambda n : int(floor(a*(n+2)) - floor(a*(n+1)))
            sage: f = Words([0,1])(wn); f
            Infinite word over [0, 1]
            sage: f[:10000] == w[:10000]
            True
            sage: f[:10000] == u[:10000] #long time
            True
            sage: words.FibonacciWord("abc")
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
        """
        if len(set(alphabet)) != 2:
            raise TypeError, "alphabet does not contain two distinct elements"
        if construction_method == "recursive":
            w = InfiniteWord_over_OrderedAlphabet(Words(alphabet), \
                    self._FibonacciWord_RecursiveConstructionIterator(alphabet), \
                    format='iterator')
            w.rename("Fibonacci word over %s, defined recursively" % w.alphabet().string_rep())
            return w
        elif construction_method == "fixed_point":
            d = {alphabet[1]:[alphabet[0]],alphabet[0]:[alphabet[0],alphabet[1]]}
            w = self.FixedPointOfMorphism(d, alphabet[0])
            w.rename("Fibonacci word over %s, defined as the fixed point of a morphism" % w.alphabet().string_rep())
            return w
        else:
            raise NotImplementedError

    def _FibonacciWord_RecursiveConstructionIterator(self,alphabet=(0,1)):
        r"""
        Iterates over the symbols of the Fibonacci word, as defined by the
        following recursive construction: the Fibonacci word is the limit
        of the sequence `S_0 = 0`, `S_1 = 01`,
        `S_n = S_{n-1}
                S_{n-2}` for `n \geq 2`.

        TESTS::

            sage: from sage.combinat.words.word_generators import WordGenerator
            sage: from itertools import islice
            sage: it = WordGenerator()._FibonacciWord_RecursiveConstructionIterator()
            sage: list(islice(it,13))
            [0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1]
        """
        Fib0 = [0]
        Fib1 = [0,1]
        n = 0
        while True:
            it = iter(Fib1[n:])
            for i in it:
                n += 1
                yield alphabet[i]
            else:
                Fib1, Fib0 = Fib1 + Fib0, Fib1

    def FixedPointOfMorphism(self, morphism, first_letter):
        r"""
        Returns the fixed point of the morphism beginning with
        first_letter.

        A *fixed point* of a morphism `\varphi` is a word
        `w` such that `\varphi(w) = w`.

        INPUT:


        -  ``morphism`` - endomorphism prolongable on
           first_letter. It must be something that WordMorphism's constructor
           understands (dict, str, ...).

        -  ``first_letter`` - the first letter of the fixed
           point


        OUTPUT:


        -  ``word`` - the fixed point of the morphism beginning
           with first_letter


        EXAMPLES::

            sage: mu = {0:[0,1], 1:[1,0]}
            sage: tm = words.FixedPointOfMorphism(mu,0); tm
            Fixed point beginning with 0 of the morphism WordMorphism: 0->01, 1->10
            sage: TM = words.ThueMorseWord()
            sage: tm[:1000] == TM[:1000]
            True

        ::

            sage: mu = {0:[0,1], 1:[0]}
            sage: f = words.FixedPointOfMorphism(mu,0); f
            Fixed point beginning with 0 of the morphism WordMorphism: 0->01, 1->0
            sage: F = words.FibonacciWord(); F
            Fibonacci word over [0, 1], defined recursively
            sage: f[:1000] == F[:1000]
            True

        ::

            sage: fp = words.FixedPointOfMorphism('a->abc,b->,c->','a'); fp
            Fixed point beginning with 'a' of the morphism WordMorphism: a->abc, b->, c->
            sage: fp[:10]
            word: abc
        """
        return WordMorphism(morphism).fixed_point(letter=first_letter)

    def CodingOfRotationWord(self, alpha, beta, x=0, alphabet=(0,1)):
        r"""
        Returns the infinite word obtained from the coding of rotation of
        parameters `(\alpha,\beta, x)` over the given two-letter
        alphabet.

        The *coding of rotation* corresponding to the parameters
        `(\alpha,\beta, x)` is the symbolic sequence
        `u = (u_n)_{n\geq 0}` defined over the binary alphabet
        `\{0, 1\}` by `u_n = 1` if
        `x+n\alpha\in[0, \beta[` and `u_n = 0` otherwise.
        See [1].

        EXAMPLES::

            sage: alpha = 0.45
            sage: beta = 0.48
            sage: w = words.CodingOfRotationWord(0.45, 0.48); w
            Coding of rotation with parameters (0.450000000000000, 0.480000000000000, 0)
            sage: w[:100]
            Finite word of length 100 over [0, 1]
            sage: w[:30]
            word: 110101010100101010101101010101

        ::

            sage: u = words.CodingOfRotationWord(0.45, 0.48, alphabet='xy')
            sage: u[:100]
            Finite word of length 100 over ['x', 'y']
            sage: u[:30]
            word: yyxyxyxyxyxxyxyxyxyxyyxyxyxyxy

        TESTS::

            sage: words.CodingOfRotationWord(0.51,0.43,alphabet=[1,0,2])
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements

        REFERENCES:

        - [1] B. Adamczewski, J. Cassaigne, On the transcendence of
          real numbers with a regular expansion, J. Number Theory 103
          (2003) 27-37.
        """
        if len(set(alphabet)) != 2:
            raise TypeError, "alphabet does not contain two distinct elements"
        from functools import partial
        f = partial(self._CodingOfRotationWord_function,alpha=alpha,beta=beta,x=x,alphabet=alphabet)
        w = InfiniteWord_over_OrderedAlphabet(Words(alphabet), f, format='function')
        w.rename("Coding of rotation with parameters (%s, %s, %s)"%(alpha,beta,x))
        return w

    def _CodingOfRotationWord_function(self, n, alpha, beta, x=0, alphabet=(0,1)):
        r"""
        Internal function that returns the symbol in position `n`
        of the coding of rotation word corresponding to the parameters
        `\alpha`, `\beta`, and `x`.

        TESTS::

            sage: alpha, beta = 0.45, 0.48
            sage: words._CodingOfRotationWord_function(3, alpha, beta)
            1
            sage: words._CodingOfRotationWord_function(10, alpha, beta)
            0
            sage: words._CodingOfRotationWord_function(17, alpha, beta)
            0
        """
        hauteur = x + n * alpha
        fracH = hauteur.frac()
        if fracH < 0:
            fracH += 1
        if 0 <= fracH < beta:
            return alphabet[1]
        else:
            return alphabet[0]

    def CharacteristicSturmianWord(self, cf, alphabet=(0, 1), repeat=True):
        r"""
        Returns the characteristic Sturmian word over the given two-letter
        alphabet with slope given by the continued fraction [0, cf[1],
        cf[2], ...]. Here cf should be an iterator.

        The *characteristic Sturmian word of slope*
        `\alpha = [0, d[1] + 1, d[2], d[3], \ldots]` is the limit
        of the sequence: `s_0 = 1`, `s_1 = 0`, ...,
        `s_{n+1} = s_n^{d_n} s_{n-1}` for `n > 0`.

        Equivalently, the `n`-th term of the characteristic
        Sturmian word of slope `\alpha` is
        `\lfloor\alpha(n+1)\rfloor - \lfloor\alpha n\rfloor + \lfloor\alpha\rfloor`.

        INPUT:


        -  ``alphabet`` - any container of length two that is
           suitable to build an instance of OrderedAlphabet (list, tuple, str,
           ...)

        -  ``cf`` - an iterator outputting integers (thought of
           as a continued fraction)


        EXAMPLES::

            sage: def cf():
            ...     yield 0
            ...     while True: yield 1
            ...
            sage: F = words.CharacteristicSturmianWord(cf()); F
            Characteristic Sturmian word over [0, 1], defined recursively
            sage: Fib = words.FibonacciWord(); Fib
            Fibonacci word over [0, 1], defined recursively
            sage: F[:10000] == Fib[:10000]
            True

        ::

            sage: def cf():
            ...     yield 0
            ...     while True: yield 1
            ...
            sage: G = words.CharacteristicSturmianWord(cf(),'rs'); G
            Characteristic Sturmian word over ['r', 's'], defined recursively
            sage: print G[:50]
            word: rsrrsrsrrsrrsrsrrsrsrrsrrsrsrrsrrsrsrrsrsrrsrrsrsr

        TESTS::

            sage: words.CharacteristicSturmianWord([1,1,1],'xyz')
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
        """
        if len(set(alphabet)) != 2:
            raise TypeError, "alphabet does not contain two distinct elements"
        if not repeat:
            d = iter(cf)
        else:
            d = cycle(cf)
        w = InfiniteWord_over_OrderedAlphabet(Words(alphabet), \
                self._CharacteristicSturmianWord_LetterIterator(d,alphabet), \
                format='iterator')
        w.rename("Characteristic Sturmian word over %s, defined recursively" \
                % w.alphabet().string_rep())
        return w

    def _CharacteristicSturmianWord_LetterIterator(self, d, alphabet=(0,1)):
        r"""
        Internal function iterating over the symbols of the characteristic
        Sturmian word of slope `[0, d[1] + 1, d[2], d[3], \ldots]`.
        This word is the limit of the sequence:
        `s_0 = 1, s_1 = 0, \ldots, s_{n+1} = s_n^{d_n} s_{n-1}`
        for `n > 0`.

        TESTS::

            sage: d = iter([1,1,1,1,1,1])
            sage: list(words._CharacteristicSturmianWord_LetterIterator(d))
            [0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1]
        """
        s0 = [0]
        s1 = [1]
        s1, s0 = s1*(d.next()-1) + s0, s1
        n = 0
        while True:
            for i in s1[n:]:
                n += 1
                yield alphabet[i]
            else:
                s1, s0 = s1*d.next() + s0, s1

    def StandardEpisturmianWord(self, directive_word):
        r"""
        Returns the standard episturmian word (or epistandard word)
        directed by directive_word. Over a 2-letter alphabet, this
        function gives characteristic Sturmian words.

        An infinite word `w` over a finite alphabet `A` is
        said to be *standard episturmian* (or *epistandard*) iff there
        exists an infinite word `x_1x_2x_3\cdots` over
        `A` (called the *directive  word* of `w`) such that
        `w` is the limit as `n` goes to infinity of
        `Pal(x_1\cdots x_n)`, where `Pal` is the iterated
        palindromic closure function.

        Note that an infinite word is *episturmian* if it has the same set
        of factors as some epistandard word.

        See for instance [1], [2], and [3].

        INPUT:


        -  ``directive word`` - an infinite word or a period of
           a periodic infinite word


        EXAMPLES::

            sage: Fibonacci = words.StandardEpisturmianWord(Word('ab')); Fibonacci
            Standard episturmian word over ['a', 'b']
            sage: Fibonacci[:25]
            word: abaababaabaababaababaabaa
            sage: Tribonacci = words.StandardEpisturmianWord(Word('abc')); Tribonacci
            Standard episturmian word over ['a', 'b', 'c']
            sage: Tribonacci[:25]
            word: abacabaabacababacabaabaca
            sage: S = words.StandardEpisturmianWord(Word('aabcabada')); S
            Standard episturmian word over ['a', 'b', 'c', 'd']
            sage: print S[:75]
            word: aabaacaabaaabaacaabaabaacaabaaabaacaabaaabaacaabaabaacaabaaabaacaabaadaabaa
            sage: S = words.StandardEpisturmianWord(Fibonacci); S
            Standard episturmian word over ['a', 'b']
            sage: S[:25]
            word: abaabaababaabaabaababaaba
            sage: S = words.StandardEpisturmianWord(Tribonacci); S
            Standard episturmian word over ['a', 'b', 'c']
            sage: S[:25]
            word: abaabacabaabaabacabaababa
            sage: words.StandardEpisturmianWord(123)
            Traceback (most recent call last):
            ...
            TypeError: directive_word is not a word, so it cannot be used to build an episturmian word
            sage: words.StandardEpisturmianWord(Words('ab'))
            Traceback (most recent call last):
            ...
            TypeError: directive_word is not a word, so it cannot be used to build an episturmian word

        REFERENCES:

        - [1] X. Droubay, J. Justin, G. Pirillo, Episturmian words and
          some constructions of de Luca and Rauzy, Theoret. Comput.
          Sci. 255 (2001) 539-553.

        - [2] J. Justin, G. Pirillo, Episturmian words and episturmian
          morphisms, Theoret. Comput. Sci. 276 (2002) 281-313.

        - [3] A. Glen, J. Justin, Episturmian words: a survey,
          Preprint, 2007, arXiv:0801.1655.
        """
        if not is_Word(directive_word):
           raise TypeError, "directive_word is not a word, so it cannot be used to build an episturmian word"
        epistandard = directive_word.parent()(\
                self._StandardEpisturmianWord_LetterIterator(directive_word), \
                format='iterator')
        epistandard.rename("Standard episturmian word over %s" % \
                epistandard.alphabet().string_rep())
        return epistandard

    def _StandardEpisturmianWord_LetterIterator(self, directive_word):
        r"""
        Internal iterating over the symbols of the standard episturmian
        word defined by the (directive) word directive_word.

        An infinite word `w` over a finite alphabet `A` is
        standard episturmian (or epistandard) iff there exists an infinite
        word `x_1x_2x_3\ldots` over `A` (called the
        directive word of `w`) such that `w` is the limit
        as `n` goes to infinity of
        `Pal(x_1x_2\cdots x_n)`, where `Pal` is the
        iterated palindromic closure function.

        INPUT:


        -  ``directive word`` - an infinite word or a finite
           word. If directive_word is finite, then it is repeated to give an
           infinite word.


        TESTS::

            sage: import itertools
            sage: it = words._StandardEpisturmianWord_LetterIterator(Word('ab'))
            sage: list(itertools.islice(it, 13))
            ['a', 'b', 'a', 'a', 'b', 'a', 'b', 'a', 'a', 'b', 'a', 'a', 'b']
        """
        if is_FiniteWord(directive_word):
           d = cycle(directive_word)
        else:
           d = iter(directive_word)
        W = directive_word.parent()
        w = W(d.next())
        n = 0
        while True:
              for x in w[n:]:
                  n += 1
                  yield x
              else:
                  w = W(w*W(d.next())).palindromic_closure()

    def MinimalSmoothPrefix(self, n):
        r"""
        This function finds and returns the minimal smooth prefix of length
        n.

        See [1] for a definition.

        INPUT:


        -  ``n`` - the desired length of the prefix


        OUTPUT:


        -  ``word`` - the prefix


        Be patient, this function can take a really long time if asked for
        a large prefix.

        EXAMPLES::

            sage: words.MinimalSmoothPrefix(10)
            word: 1212212112

        REFERENCES: [1] S. Brlek, G. Melançon, G. Paquin, Properties of the
        extremal infinite smooth words, Discrete Math. Theor. Comput. Sci.
        9 (2007) 33-49.
        """
        tab = []
        W = Words([1, 2])
        suff1 = W([1, 2, 2]).phi_inv()
        suff2 = W([2, 2]).phi_inv()
        w = [1]
        tab = _build_tab(1, tab, W)
        for k in xrange(1, n):
            if suff1._phi_inv_tab(tab) < suff2._phi_inv_tab(tab):
                w.append(1)
                tab = _build_tab(1, tab, W)
            else:
                w.append(2)
                tab = _build_tab(2, tab, W)
        return W(w)

    def RandomWord(self, n, m=2, alphabet=None):
        """
        Returns a random word of length n over the given m-letter alphabet
        (default alphabet is 0,1,...,m-1).

        INPUT:


        -  ``n`` - integer, the length of the word

        -  ``m`` - integer (default 2), the size of the output
           alphabet

        -  ``alphabet`` - any container of length m that is
           suitable to build an instance of OrderedAlphabet (list, tuple, str,
           ...)


        EXAMPLES::

            sage: words.RandomWord(10)         # random results
            word: 0110100101
            sage: words.RandomWord(10, 4)      # random results
            word: 0322313320
            sage: words.RandomWord(100, 7)
            Finite word of length 100 over [0, 1, 2, 3, 4, 5, 6]
            sage: words.RandomWord(100, 7, range(-3,4))
            Finite word of length 100 over [-3, -2, -1, 0, 1, 2, 3]
            sage: words.RandomWord(100, 5, "abcde")
            Finite word of length 100 over ['a', 'b', 'c', 'd', 'e']
            sage: words.RandomWord(17, 5, "abcde")     # random results
            word: dcacbbecbddebaadd

        TESTS::

            sage: words.RandomWord(2,3,"abcd")
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain 3 distinct elements
        """
        if alphabet is None:
            alphabet = range(m)
        if len(set(alphabet)) != m:
            raise TypeError, "alphabet does not contain %s distinct elements" % m
        return FiniteWord_over_OrderedAlphabet(Words(alphabet), [alphabet[randint(0,m-1)] for i in xrange(n)])

words = WordGenerator()

class ChristoffelWord_Lower(FiniteWord_over_OrderedAlphabet):
    def __init__(self, p, q, alphabet=(0,1)):
        r"""
        Returns the lower Christoffel word of slope `p/q`, where
        `p` and `q` are relatively prime non-negative
        integers, over the given two-letter alphabet.

        The *Christoffel word of slope `p/q`* is obtained from the
        Cayley graph of `\ZZ/(p+q)\ZZ` with generator
        `q` as follows. If `u \rightarrow v` is an edge in
        the Cayley graph, then `v = u + p \mod{p+q}`. Label the
        edge `u \rightarrow v` by alphabet[1] if `u < v`
        and alphabet[0] otherwise. The Christoffel word is the word
        obtained by reading the edge labels along the cycle beginning from
        0.

        EXAMPLES::

            sage: from sage.combinat.words.word_generators import ChristoffelWord_Lower
            sage: w = ChristoffelWord_Lower(1,0); w
            Lower Christoffel word of slope 1/0 over the alphabet [0, 1]
            sage: print w
            word: 1

        ::

            sage: w = ChristoffelWord_Lower(0,1,'xy'); w
            Lower Christoffel word of slope 0/1 over the alphabet ['x', 'y']
            sage: print w
            word: x

        ::

            sage: w = ChristoffelWord_Lower(1,1); w
            Lower Christoffel word of slope 1/1 over the alphabet [0, 1]
            sage: print w
            word: 01

        ::

            sage: w = ChristoffelWord_Lower(4,7); w
            Lower Christoffel word of slope 4/7 over the alphabet [0, 1]
            sage: print w
            word: 00100100101

        ::

            sage: w = ChristoffelWord_Lower(4,7,alphabet='ab'); w
            Lower Christoffel word of slope 4/7 over the alphabet ['a', 'b']
            sage: print w
            word: aabaabaabab

        TESTS::

            sage: ChristoffelWord_Lower(4,8)
            Traceback (most recent call last):
            ...
            ValueError: 4 and 8 are not relatively prime
            sage: ChristoffelWord_Lower(17, 39, 'xyz')
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
            sage: w = ChristoffelWord_Lower(4,7)
            sage: w2 = loads(dumps(w))
            sage: w == w2
            True
            sage: type(w2)
            <class 'sage.combinat.words.word_generators.ChristoffelWord_Lower'>
            sage: _ = w2.standard_factorization() # hackish test for self.__p and self.__q
        """
        if len(set(alphabet)) != 2:
            raise TypeError, "alphabet does not contain two distinct elements"
        # Compute gcd of p, q; raise TypeError if not 1.
        if gcd(p,q) != 1:
            raise ValueError, "%s and %s are not relatively prime" % (p, q)
        # Compute the Christoffel word
        w = []
        u = 0
        if (p, q) == (0, 1):
            w = [alphabet[0]]
        else:
            for i in range(p + q):
                v = (u+p) % (p+q)
                new_letter = alphabet[0] if u < v else alphabet[1]
                w.append(new_letter)
                u = v
        super(ChristoffelWord_Lower, self).__init__(Words(alphabet), w)
        self.__p = p
        self.__q = q
        self.rename("Lower Christoffel word of slope %s/%s over the alphabet %s" % (self.__p, self.__q, self.alphabet().string_rep()))

    def markoff_number(self):
        r"""
        Returns the Markoff number associated to the Christoffel word
        self.

        The *Markoff number* of a Christoffel word `w` is
        `1/3 trace M(w)`, where `M(w)` is the
        `2\times 2` matrix obtained by applying the morphism: 0 -
        matrix(2,[2,1,1,1]) 1 - matrix(2,[5,2,2,1])

        EXAMPLES::

            sage: from sage.combinat.words.word_generators import ChristoffelWord_Lower
            sage: w0 = ChristoffelWord_Lower(4,7)
            sage: w1, w2 = w0.standard_factorization()
            sage: (m0,m1,m2) = (w.markoff_number() for w in (w0,w1,w2))
            sage: (m0,m1,m2)
            (294685, 13, 7561)
            sage: m0**2 + m1**2 + m2**2 == 3*m0*m1*m2
            True
        """
        from sage.matrix.constructor import matrix
        eta = {0:matrix(2,[2,1,1,1]), 1:matrix(2,[5,2,2,1])}
        M = matrix(2,[1,0,0,1])
        for a in self._word_content:
            M *= eta[a]
        return M.trace()/3

    def standard_factorization(self):
        r"""
        Returns the standard factorization of the Christoffel word self.

        The *standard factorization* of a Christoffel word `w` is
        the unique factorization of `w` into two Christoffel
        words.

        EXAMPLES::

            sage: from sage.combinat.words.word_generators import ChristoffelWord_Lower
            sage: w = ChristoffelWord_Lower(5,9)
            sage: print w
            word: 00100100100101
            sage: w1, w2 = w.standard_factorization()
            sage: print w1
            word: 001
            sage: print w2
            word: 00100100101

        ::

            sage: w = ChristoffelWord_Lower(51,37)
            sage: w1, w2 = w.standard_factorization()
            sage: w1
            Lower Christoffel word of slope 11/8 over the alphabet [0, 1]
            sage: w2
            Lower Christoffel word of slope 40/29 over the alphabet [0, 1]
            sage: w1 * w2 == w
            True
        """
        p, q = self.__p, self.__q
        index = 0
        u = 0
        for i in range(p + q):
            v = (u+p) % (p+q)
            if v == 1:
                index = i
                break
            u = v
        w1, w2 = self[:index+1], self[index+1:]
        return Factorization([ChristoffelWord_Lower(w1.count(1),w1.count(0)),
                ChristoffelWord_Lower(w2.count(1),w2.count(0))])
