# coding=utf-8
#*****************************************************************************
#       Copyright (C) 2008 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from utils import isint, is_iterable
from itertools import imap, ifilterfalse
from sage.structure.sage_object import SageObject
from sage.rings.infinity import infinity
from sage.matrix.constructor import Matrix
from sage.rings.integer_ring import IntegerRing
from sage.misc.latex import latex
from sage.combinat.words.word import is_FiniteWord, FiniteWord_over_OrderedAlphabet
from sage.combinat.words.words import is_Words, Words
from sage.combinat.words.alphabet import OrderedAlphabet

class WordMorphism(SageObject):
    r"""
    TESTS:
        sage: wm = WordMorphism('a->ab,b->ba')
        sage: wm == loads(dumps(wm))
        True
    """
    def __init__(self, data, codomain=None):
        r"""
        Construction of the morphism.

        EXAMPLES:

        1. If data is a str
            sage: print WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: print WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: print WordMorphism('a->abc,b->bca,c->cab')
            WordMorphism: a->abc, b->bca, c->cab
            sage: print WordMorphism('a->abdsf,b->hahdad,c->asdhasd')
            WordMorphism: a->abdsf, b->hahdad, c->asdhasd
            sage: print WordMorphism('(->(),)->)(')
            WordMorphism: (->(), )->)(
            sage: print WordMorphism('a->53k,b->y5?,$->49i')
            WordMorphism: $->49i, a->53k, b->y5?

        An erasing morphism:
            sage: print WordMorphism('a->ab,b->')
            WordMorphism: a->ab, b->

        Use the arrows ('->') correctly:
            sage: WordMorphism('a->ab,b-')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '-')
            sage: WordMorphism('a->ab,b')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '')
            sage: WordMorphism('a->ab,a-]asdfa')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '-]')

        Each letter must be defined only once:
            sage: WordMorphism('a->ab,a->ba')
            Traceback (most recent call last):
            ...
            ValueError: The image of 'a' is defined twice.

        2. From a dictionary.
            sage: print WordMorphism({"a":"ab","b":"ba"})
            WordMorphism: a->ab, b->ba
            sage: print WordMorphism({2:[4,5,6],3:[1,2,3]})
            WordMorphism: 2->456, 3->123
            sage: print WordMorphism({'a':['a',6,'a'],6:[6,6,6,'a']})
            WordMorphism: 6->666a, a->a6a

        The image of a letter can be a set, but the order is not preserved:
            sage: print WordMorphism({2:[4,5,6],3:set([4,1,8])}) #random results
            WordMorphism: 2->456, 3->814

        If the image of a letter is not iterable, it is considered as a letter:
            sage: print WordMorphism({0:1, 1:0})
            WordMorphism: 0->1, 1->0
            sage: print WordMorphism({0:123, 1:789})
            WordMorphism: 0->123, 1->789
            sage: print WordMorphism({2:[4,5,6], 3:123})
            WordMorphism: 2->456, 3->123

        3. From a WordMorphism
            sage: print WordMorphism(WordMorphism('a->ab,b->ba'))
            WordMorphism: a->ab, b->ba

        TESTS:
            sage: print WordMorphism(',,,a->ab,,,b->ba,,')
            WordMorphism: a->ab, b->ba
        """
        if isinstance(data, WordMorphism):
            self._domain = data._domain
            self._codomain = data._codomain
            self._morph = data._morph
        else:
            if isinstance(data, str):
                data = self._build_dict(data)
            elif not isinstance(data, dict):
                raise NotImplementedError

            if codomain is None:
                codomain = self._build_codomain(data)

            if not is_Words(codomain):
                raise TypeError, "the codomain must be a Words domain"
            self._codomain = codomain

            self._morph = {}

            dom_alph = list()
            cod = self._codomain
            for (key,val) in data.iteritems():
                dom_alph.append(key)
                self._morph[key] = cod(val)

            dom_alph.sort()
            self._domain = Words(dom_alph)

    def _build_dict(self, s):
        r"""
        Parse the string input to WordMorphism and build the dictionary
        it represents.

        TESTS:
            sage: wm = WordMorphism('a->ab,b->ba')
            sage: wm._build_dict('a->ab,b->ba') == {'a': 'ab', 'b': 'ba'}
            True
            sage: wm._build_dict('a->ab,a->ba')
            Traceback (most recent call last):
            ...
            ValueError: The image of 'a' is defined twice.
            sage: wm._build_dict('a->ab,b>ba')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '>b')
        """
        tmp_dict = {}
        for fleche in s.split(','):
            if len(fleche) == 0:
                continue

            if len(fleche) < 3 or fleche[1:3] != '->':
                raise ValueError, "The second and third characters must be '->' (not '%s')"%fleche[1:3]

            lettre = fleche[0]
            image  = fleche[3:]

            if lettre in tmp_dict:
                raise ValueError, "The image of %r is defined twice." %lettre

            tmp_dict[lettre] = image
        return tmp_dict

    def _build_codomain(self, data):
        r"""
        Returns a Words domain containing all the letter in the keys of
        data (which must be a dictionnary).

        TESTS:
            1. If the image of all the letters are iterable :
            sage: wm = WordMorphism('a->ab,b->ba')
            sage: wm._build_codomain({'a': 'ab', 'b': 'ba'})
            Words over Ordered Alphabet ['a', 'b']
            sage: wm._build_codomain({'a': 'dcb', 'b': 'a'})
            Words over Ordered Alphabet ['a', 'b', 'c', 'd']
            sage: wm._build_codomain({2:[4,5,6],3:[1,2,3]})
            Words over Ordered Alphabet [1, 2, 3, 4, 5, 6]
            sage: wm._build_codomain({2:[4,5,6],3:set([4,1,8])})
            Words over Ordered Alphabet [1, 4, 5, 6, 8]

            2. If the image of a letter is not iterable, it is considered as
            a letter:
            sage: wm._build_codomain({2:[4,5,6],3:123})
            Words over Ordered Alphabet [4, 5, 6, 123]
            sage: wm._build_codomain({0:1, 1:0, 2:2})
            Words over Ordered Alphabet [0, 1, 2]
        """
        codom_alphabet = set()
        for key,val in data.iteritems():
            if is_iterable(val):
                codom_alphabet.update(val)
            else:
                codom_alphabet.add(val)
        return Words(sorted(codom_alphabet))

    def __eq__(self, other):
        r"""
        Returns True if self is equal to other.

        EXAMPLES:
            sage: n = WordMorphism('a->a,b->aa,c->aaa')
            sage: n**3 == n**1
            True
            sage: WordMorphism('b->ba,a->ab') == WordMorphism('a->ab,b->ba')
            True
            sage: WordMorphism('b->ba,a->ab') == WordMorphism({"a":"ab","b":"ba"})
            True
            sage: m = WordMorphism({0:[1,2,3],1:[4,5,6]}); print m
            WordMorphism: 0->123, 1->456
            sage: o = WordMorphism('0->123,1->456'); print o
            WordMorphism: 0->123, 1->456
            sage: m == o
            False
        """
        return self._morph == other._morph

    def __repr__(self):
        r"""
        Returns the morphism in str (for display).

        EXAMPLES:
            sage: WordMorphism('a->ab,b->ba')
            Morphism from Words over Ordered Alphabet ['a', 'b'] to Words over Ordered Alphabet ['a', 'b']
            sage: d = {0:[0,1],1:[1,0]}
            sage: WordMorphism(d)
            Morphism from Words over Ordered Alphabet [0, 1] to Words over Ordered Alphabet [0, 1]
        """
        return "Morphism from %s to %s" %(self.domain(),self.codomain())

    def __str__(self):
        r"""
        Returns the morphism in str (for display).

        EXAMPLES:
            sage: print WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: print WordMorphism('b->ba,a->ab')
            WordMorphism: a->ab, b->ba
            sage: d = {0:[0,1],1:[1,0]}
            sage: print WordMorphism(d)
            WordMorphism: 0->01, 1->10
        """
        l = [str(lettre) + '->' + image.string_rep() for lettre,image in self._morph.iteritems()]

        return "WordMorphism: %s" % ', '.join(sorted(l))

    def __call__(self, w, order=1):
        r"""
        Returns the image of w under self to the given order.

        INPUT:
            w -- finite word in the domain of self, must be of length one if
                 order is infinity
            order -- integer or plus infinity (default: 1)

        OUTPUT:
            word -- order-th iterated image under self of w

        EXAMPLES:
        The image of a word under a morphism:

        1. The image of a finite word under a morphism:
            sage: tm = WordMorphism ('a->ab,b->ba')
            sage: tm('a')
            word: ab
            sage: tm('aabababb')
            word: ababbaabbaabbaba

        2. The iterated image of a word:
            sage: tm('a', 2)
            word: abba
            sage: tm('aba', 3)
            word: abbabaabbaababbaabbabaab

        3. The infinitely iterated image of a letter:
            sage: tm('a', oo)
            Fixed point beginning with 'a' of the morphism WordMorphism: a->ab, b->ba

        4. The image of an infinite word:
            sage: t = words.ThueMorseWord()
            sage: n = WordMorphism({0:[0, 1], 1:[1, 0]})
            sage: n(t)
            Infinite word over [0, 1]
            sage: n(t, 3)
            Infinite word over [0, 1]
            sage: n(t)[:39]
            word: 011010011001011010010110011010011001011
            sage: n(t)[:1000] == t[:1000]
            True

            sage: w = words.FibonacciWord()
            sage: m = WordMorphism({0:'a', 1:'b'})
            sage: m(w)
            Infinite word over ['a', 'b']
            sage: f = words.FibonacciWord('ab')
            sage: f[:1000] == m(w)[:1000]
            True

            sage: w = words.FibonacciWord("ab")
            sage: m = WordMorphism('a->01,b->101')
            sage: m(w)
            Infinite word over ['0', '1']
            sage: m(w)[:39]
            word: 011010101101011010101101010110101101010

        The word must be in the domain of self:
            sage: tm('0021')
            Traceback (most recent call last):
            ...
            TypeError: cannot apply to incompatible set

        The order must be a positive integer or plus infinity:
            sage: tm('a', -1)
            Traceback (most recent call last):
            ...
            TypeError: order (-1) must be a positive integer or plus infinity
            sage: tm('a', 6.7)
            Traceback (most recent call last):
            ...
            TypeError: order (6.70000000000000) must be a positive integer or plus infinity

        Infinitely iterated image of a word is defined only for those of length one:
            sage: tm('aba',oo)
            Traceback (most recent call last):
            ...
            TypeError: For infinite powers, the length of the word must be 1 (not 3)

        self must be prolongable on the given letter for infinitely iterated image:
            sage: m = WordMorphism('a->ba,b->ab')
            sage: m('a', oo)
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on a

        TESTS:
            sage: for i in range(8):
            ...     tm('a', i)
            ...
            word: a
            word: ab
            word: abba
            word: abbabaab
            word: abbabaabbaababba
            word: abbabaabbaababbabaababbaabbabaab
            Finite word of length 64 over ['a', 'b']
            Finite word of length 128 over ['a', 'b']
            sage: m = WordMorphism('a->,b->')
            sage: m('')
            word:
        """
        try:
            w = self._domain(w)
        except:
            raise TypeError, "cannot apply to incompatible set"

        if order is infinity:
            if len(w) != 1:
                raise TypeError, "For infinite powers, the length of the word must be 1 (not %s)"%len(w)
            return self.fixed_point(letter=w[0])

        if not isint(order) or order < 0 :
            raise TypeError, "order (%s) must be a positive integer or plus infinity" % order
        elif order == 0:
            return w
        elif order == 1:
            if is_FiniteWord(w):
                if len(w) == 1:
                    return self._morph[w[0]]
                else:
                    s = sum(Matrix(self)*Matrix(w.parikh_vector()).transpose())
                    if isint(s):
                        image_length = s
                    else:
                        image_length = s[0]
                    return self.codomain()((x for y in w for x in self(y)),
                                           part=slice(image_length))
            else:
                return self.codomain()(x for y in w for x in self(y))
        elif order > 1:
            return self(self(w, order-1))

    def __mul__(self, other):
        r"""
        Returns the morphism self*other

        EXAMPLES:
            sage: m = WordMorphism('a->ab,b->ba')
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: print fibo*m
            WordMorphism: a->aba, b->aab
            sage: print fibo*fibo
            WordMorphism: a->aba, b->ab
            sage: print m*fibo
            WordMorphism: a->abba, b->ab
            sage: m = WordMorphism('a->ab,b->ba')
            sage: n = WordMorphism('a->a,b->aa,c->aaa')
            sage: p1 = n*m; p2 = m*n
            sage: print p1
            WordMorphism: a->aaa, b->aaa
            sage: p1.domain()
            Words over Ordered Alphabet ['a', 'b']
            sage: p1.codomain()
            Words over Ordered Alphabet ['a']
            sage: print p2
            WordMorphism: a->ab, b->abab, c->ababab
            sage: p2.domain()
            Words over Ordered Alphabet ['a', 'b', 'c']
            sage: p2.codomain()
            Words over Ordered Alphabet ['a', 'b']
            sage: m = WordMorphism('a->b,b->c,c->a')
            sage: WordMorphism('')*m
            Traceback (most recent call last):
            ...
            TypeError: cannot apply to incompatible set
            sage: print m*WordMorphism('')
            WordMorphism:
        """
        #TODO : Est-ce que c'est le comportement que l'on veut pour le produit
        #par le morphisme vide? Voir lignes ci-haut.
        return WordMorphism(dict((key, self(w)) for (key, w) in other._morph.iteritems()))

    def __pow__(self, exp):
        r"""
        Returns the power of self with exponent = exp.

        EXAMPLES:
            sage: m = WordMorphism('a->ab,b->ba')
            sage: m^1.5
            Traceback (most recent call last):
              ...
            ValueError: exponent (1.50000000000000) must be an integer
            sage: print m^2
            WordMorphism: a->abba, b->baab
            sage: print m^3
            WordMorphism: a->abbabaab, b->baababba
            sage: print m^-2
            Traceback (most recent call last):
              ...
            ValueError: exponent (-2) must be strictly positive
            sage: n = WordMorphism('a->ba,b->abc')
            sage: print n^1
            WordMorphism: a->ba, b->abc
            sage: n^2
            Traceback (most recent call last):
              ...
            TypeError: cannot apply to incompatible set
        """
        #If exp is not an integer
        if not isint(exp):
            raise ValueError, "exponent (%s) must be an integer" %exp

        #If exp is negative
        elif exp <= 0:
            raise ValueError, "exponent (%s) must be strictly positive" %exp

        #Base of induction
        elif exp == 1:
            return self

        else:
            nexp = int(exp / 2)
            over = exp % 2
            res = (self * self) ** nexp
            if over == 1:
                res *= self
            return res

    def _matrix_(self, R=None):
        r"""
        Returns the incidence matrix of the morphism over the specified ring.

        EXAMPLES:
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: tm = WordMorphism('a->ab,b->ba')
            sage: Mfibo = (fibo)._matrix_(); Mfibo
            [1 1]
            [1 0]
            sage: Mtm = (tm)._matrix_(); Mtm
            [1 1]
            [1 1]
            sage: Mtm * Mfibo == (tm*fibo)._matrix_()
            True
            sage: Mfibo * Mtm == (fibo*tm)._matrix_()
            True
            sage: Mfibo.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: p = Mfibo.charpoly(); p
            x^2 - x - 1
            sage: p.roots(ring=RR, multiplicities=False)
            [-0.618033988749895, 1.61803398874989]
        """
        if R is None:
            return self.incidence_matrix()
        else:
            return self.incidence_matrix().change_ring(R)

    def incidence_matrix(self):
        r"""
        Returns the incidence matrix of the morphism. The order of the rows and
        column are given by the order defined on the alphabet of the domain
        and the codomain.

        The matrix returned is over the integers.  If a different ring
        is desired, use either the change_ring function or the matrix
        function.

        EXAMPLES:
            sage: m = WordMorphism('a->abc,b->a,c->c')
            sage: m.incidence_matrix()
            [1 1 0]
            [1 0 0]
            [1 0 1]
            sage: m = WordMorphism('a->abc,b->a,c->c,d->abbccccabca,e->abc')
            sage: m.incidence_matrix()
            [1 1 0 3 1]
            [1 0 0 3 1]
            [1 0 1 5 1]
        """
        L = []
        for b in self.domain().alphabet():
            w = self(b)
            L.append(w.parikh_vector())

        M = Matrix(IntegerRing(), L).transpose()
        return M


    def domain(self):
        r"""
        Returns domain of self.

        EXAMPLES:
            sage: WordMorphism('a->ab,b->a').domain()
            Words over Ordered Alphabet ['a', 'b']
            sage: WordMorphism('b->ba,a->ab').domain()
            Words over Ordered Alphabet ['a', 'b']
            sage: WordMorphism('6->ab,y->5,0->asd').domain()
            Words over Ordered Alphabet ['0', '6', 'y']
        """
        return self._domain

    def codomain(self):
        r"""
        Returns the domain of the images.

        EXAMPLES :
            sage: WordMorphism('a->ab,b->a').codomain()
            Words over Ordered Alphabet ['a', 'b']
            sage: WordMorphism('6->ab,y->5,0->asd').codomain()
            Words over Ordered Alphabet ['5', 'a', 'b', 'd', 's']
        """
        return self._codomain

    def is_endomorphism(self):
        r"""
        Returns True if the codomain is a subset of the domain.

        EXAMPLES :
            sage: WordMorphism('a->ab,b->a').is_endomorphism()
            True
            sage: WordMorphism('6->ab,y->5,0->asd').is_endomorphism()
            False
            sage: WordMorphism('a->a,b->aa,c->aaa').is_endomorphism()
            True
            sage: Wabc = Words('abc')
            sage: m = WordMorphism('a->a,b->aa,c->aaa',codomain = Wabc)
            sage: m.is_endomorphism()
            True
        """
        return self.codomain() <= self.domain()

    def images(self):
        r"""
        Returns the list of all the images of the letters of the alphabet
        under self.

        EXAMPLES :
            sage: WordMorphism('a->ab,b->a').images()
            [word: ab, word: a]
            sage: WordMorphism('6->ab,y->5,0->asd').images()
            [word: 5, word: asd, word: ab]
        """
        return self._morph.values()

    def reversal(self):
        r"""
        Returns the reversal of self.

        EXAMPLES :
            sage: print WordMorphism('6->ab,y->5,0->asd').reversal()
            WordMorphism: 0->dsa, 6->ba, y->5
            sage: print WordMorphism('a->ab,b->a').reversal()
            WordMorphism: a->ba, b->a
        """
        return WordMorphism(dict((key, w.reversal()) for (key, w) in self._morph.iteritems()))

    def is_empty(self):
        r"""
        Returns True if the cardinality of the domain is zero and False otherwise.

        EXAMPLES:
            sage: WordMorphism('').is_empty()
            True
            sage: WordMorphism('a->a').is_empty()
            False
        """
        return len(self._morph) == 0

    def is_erasing(self):
        r"""
        Returns True if self is an erasing morphism, i.e.
        the image of a letter is the empty word.

        EXAMPLES :
            sage: WordMorphism('a->ab,b->a').is_erasing()
            False
            sage: WordMorphism('6->ab,y->5,0->asd').is_erasing()
            False
            sage: WordMorphism('6->ab,y->5,0->asd,7->').is_erasing()
            True
            sage: WordMorphism('').is_erasing()
            False
        """
        for image in self.images():
            if len(image) == 0:
                return True
        return False

    def is_identity(self):
        r"""
        Returns True if self is the identity morphism.

        EXAMPLES:
            sage: m = WordMorphism('a->a,b->b,c->c,d->e')
            sage: m.is_identity()
            False
            sage: WordMorphism('a->a,b->b,c->c').is_identity()
            True
            sage: WordMorphism('a->a,b->b,c->cb').is_identity()
            False
            sage: m = WordMorphism('a->b,b->c,c->a')
            sage: (m^2).is_identity()
            False
            sage: (m^3).is_identity()
            True
            sage: (m^4).is_identity()
            False
            sage: WordMorphism('').is_identity()
            True
            sage: WordMorphism({0:[0],1:[1]}).is_identity()
            True
        """
        if self.domain() != self.codomain():
            return False

        for letter in self.domain().alphabet():
            if self(letter) != self.codomain()(letter):
                return False
        return True

    def is_involution(self):
        r"""
        Returns True if self is an involution, i.e. its square
        is the identity.

        EXAMPLES:
            sage: WordMorphism('a->b,b->a').is_involution()
            True
            sage: WordMorphism('').is_involution()
            True
            sage: WordMorphism({0:[1],1:[0]}).is_involution()
            True
        """
        return (self*self).is_identity()

    def _check_primitive(self):
        r"""
        Returns True if all the letters of the domain appear in all the images
        of letters of the domain.

        INPUT:
            self -- the codomain must be an instance of Words

        EXAMPLES:
            sage: m = WordMorphism('a->ab,b->ba')
            sage: m._check_primitive()
            True
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo._check_primitive()
            False
            sage: print WordMorphism({2:[4,5,6],3:[4,1,8]})
            WordMorphism: 2->456, 3->418
            sage: WordMorphism({2:[4,5,6],3:[4,1,8]})._check_primitive()
            False
    	"""
        if not is_Words(self.codomain()):
            raise TypeError, "codomain of self(=%s) must be an instance of Words"%self

        dom_alphabet = self.domain().alphabet()

        for image in self.images():
            if not dom_alphabet <= OrderedAlphabet(set(image)):
                return False
        else:
            return True

    def is_primitive(self):
        r"""
        Returns True if self is primitive.

        INPUT:
            self -- an endomorphism

        -In English:
        A morphism $\varphi$ is \emph{primitive} if there exists an positive
        integer $k$ such that for all $\alpha\in\Sigma$, $\varphi^k(\alpha)$
        contains all the letters of $\Sigma$.

        -En français:
        Un morphisme $\varphi$ est \emph{primitif} s'il existe un nombre naturel
        $k$ tel que pour tout $\alpha\in\Sigma, \varphi^k(\alpha)$ contient
        toutes les lettres de $\Sigma$.

        EXAMPLES:
            sage: tm = WordMorphism('a->ab,b->ba')
            sage: tm.is_primitive()
            True
            sage: fibo = WordMorphism('a->ab,b->a');
            sage: fibo.is_primitive()
            True
            sage: m = WordMorphism('a->bb,b->aa')
            sage: m.is_primitive()
            False
            sage: f = WordMorphism({0:[1],1:[0]})
            sage: f.is_primitive()
            False
            sage: m = WordMorphism('a->bb,b->aac')
            sage: m.is_primitive()
            Traceback (most recent call last):
            ...
            TypeError: self (=WordMorphism: a->bb, b->aac) is not a endomorphism
            sage: m = WordMorphism('a->,b->',codomain=Words('ab'))
            sage: m.is_primitive()
            False
            sage: m = WordMorphism('a->,b->')
            sage: m.is_primitive()
            False
        """
        if not self.is_endomorphism():
            raise TypeError, "self (=%s) is not a endomorphism"%self

        dom_alphabet = self.domain().alphabet()
        dim = self.domain().size_of_alphabet()

        power = self
        for k in range(dim):
            if power._check_primitive():
                return True
            power *= self
        else:
            return False

    def is_prolongable(self, letter):
        r"""
        Returns True if self is prolongable on letter.

        A morphism $\varphi$ is prolongable on a letter $a$ if $a$ is a prefix
        of $\varphi(a)$.

        INPUT:
            self -- the codomain must be an instance of Words
            letter -- in the domain alphabet

        OUTPUT:
            boolean -- if self is prolongable on letter.

        EXAMPLES:
            sage: WordMorphism('a->ab,b->a').is_prolongable(letter='a')
            True
            sage: WordMorphism('a->ab,b->a').is_prolongable(letter='b')
            False
            sage: WordMorphism('a->ba,b->ab').is_prolongable(letter='b')
            False
            sage: (WordMorphism('a->ba,b->ab')^2).is_prolongable(letter='b')
            True
            sage: WordMorphism('a->ba,b->').is_prolongable(letter='b')
            False
            sage: WordMorphism('a->bb,b->aac').is_prolongable(letter='a')
            False

        TESTS:
            sage: WordMorphism('a->ab,b->b,c->ba').is_prolongable(letter='d')
            Traceback (most recent call last):
            ...
            TypeError: letter (=d) is not in the domain alphabet (=Ordered Alphabet ['a', 'b', 'c'])
            sage: n0, n1 = matrix(2,[1,1,1,0]), matrix(2,[2,1,1,0])
            sage: n = {'a':n0, 'b':n1}
            sage: WordMorphism(n).is_prolongable(letter='a') #todo: not implemented
            Traceback (most recent call last):
            ...
            TypeError: codomain of self must be an instance of Words
        """
        if not is_Words(self.codomain()):
            raise TypeError, "codomain of self must be an instance of Words"


        if letter not in self.domain().alphabet():
            raise TypeError, "letter (=%s) is not in the domain alphabet (=%s)"\
                                %(letter, self.domain().alphabet())
        image = self(letter)
        return not image.is_empty() and letter == image[0]

    def letter_iterator(self, letter):
        r"""
        Returns an iterator of the letters of the fixed point of self starting
        with letter.

        If w is the word, then this iterator:
            outputs the elements of morphism[ w[i] ],
            appends morphism[ w[i+1] ] to w,
            increments i.

        INPUT:
            self -- an endomorphism, must be prolongable on letter
            letter -- in the domain of self

        OUTPUT:
            iterator -- iterator of the fixed point

        EXAMPLES:
            sage: m = WordMorphism('a->abc,b->,c->')
            sage: list(m.letter_iterator('b'))
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on b
            sage: list(m.letter_iterator('a'))
            ['a', 'b', 'c']
            sage: m = WordMorphism('a->aa,b->aac')
            sage: list(m.letter_iterator('a'))
            Traceback (most recent call last):
            ...
            TypeError: self (=WordMorphism: a->aa, b->aac) is not a endomorphism
        """
        if not self.is_endomorphism():
            raise TypeError, "self (=%s) is not a endomorphism"%self

        if not self.is_prolongable(letter=letter):
            raise TypeError, "self must be prolongable on %s"%letter

        w = self(letter)
        i = 0
        while w:
            for a in self(w[i]):
                yield a
            else:
                i += 1
                if i >= len(w):
                    raise StopIteration
                else:
                    w *= self(w[i])

        ## This appears to be slower than the above, and I'm not sure why.
        #while w:
        #    for a in self[w.pop(0)]:
        #        yield a
        #    w += self[w[0]]

    def fixed_point(self, letter):
        r"""
        Returns the fixed point of self beginning by the given letter.

        A fixed point of morphism $\varphi$ is a word $w$ such that
        $\varphi(w) = w$.

        INPUT:
            self -- an endomorphism, must be prolongable on letter
            letter -- in the domain of self, the first letter of
                      the fixed point.

        OUTPUT:
            word -- the fixed point of self beginning with letter.

        EXAMPLES:

        1. Infinite fixed point:
            sage: WordMorphism('a->ab,b->ba').fixed_point(letter='a')
            Fixed point beginning with 'a' of the morphism WordMorphism: a->ab, b->ba
            sage: WordMorphism('a->ab,b->a').fixed_point(letter='a')
            Fixed point beginning with 'a' of the morphism WordMorphism: a->ab, b->a
            sage: WordMorphism('a->ab,b->b,c->ba').fixed_point(letter='a')
            Fixed point beginning with 'a' of the morphism WordMorphism: a->ab, b->b, c->ba

        2. Infinite fixed point of an erasing morphism:
            sage: WordMorphism('a->ab,b->,c->ba').fixed_point(letter='a')
            Fixed point beginning with 'a' of the morphism WordMorphism: a->ab, b->, c->ba

        3. Finite fixed point:
            sage: WordMorphism('a->ab,b->b,c->ba').fixed_point(letter='b')
            word: b

        4. Finite fixed point of an erasing morphism:
            sage: m = WordMorphism('a->abc,b->,c->')
            sage: fp = m.fixed_point('a'); fp
            Fixed point beginning with 'a' of the morphism WordMorphism: a->abc, b->, c->
            sage: fp[:10]
            word: abc
            sage: m = WordMorphism('a->ba,b->')
            sage: m('ba')
            word: ba
            sage: m.fixed_point('a') #todo: not implemented
            word: ba

        5. Fixed point of a power of a morphism:
            sage: m = WordMorphism('a->ba,b->ab')
            sage: (m^2).fixed_point(letter='a')
            Fixed point beginning with 'a' of the morphism WordMorphism: a->abba, b->baab

        TESTS:
            sage: WordMorphism('a->ab,b->,c->ba').fixed_point(letter='b')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on b
            sage: WordMorphism('a->ab,b->,c->ba').fixed_point(letter='c')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on c
            sage: WordMorphism('a->ab,b->,c->ba').fixed_point(letter='d')
            Traceback (most recent call last):
            ...
            TypeError: letter (=d) is not in the domain alphabet (=Ordered Alphabet ['a', 'b', 'c'])
            sage: WordMorphism('a->aa,b->aac').fixed_point(letter='a')
            Traceback (most recent call last):
            ...
            TypeError: self (=WordMorphism: a->aa, b->aac) is not a endomorphism
        """
        if not self.is_endomorphism():
            raise TypeError, "self (=%s) is not a endomorphism"%self

        if not self.is_prolongable(letter=letter):
            raise TypeError, "self must be prolongable on %s"%letter

        image = self(letter)

        if len(image) == 1:
            return image

        # Construct the word.
        w = self.codomain()(self.letter_iterator(letter), format='iterator')
        w.rename("Fixed point beginning with %r of the morphism %s"%(letter,self))
        return w

    def list_fixed_points(self):
        r"""
        Returns the list of all fixed points of self.

        EXAMPLES
            sage: WordMorphism('a->ab,b->ba').list_fixed_points() #not implemented
            [Fixed point beginning with 'a' of the morphism WordMorphism: a->ab, b->ba,
            Fixed point beginning with 'b' of the morphism WordMorphism: a->ab, b->ba]
        """
        raise NotImplementedError

    def conjugate(self, pos):
        r"""
        Returns the morphism where each image of self is conjugate
        of parameter pos.

        Retourne le morphisme dont toutes les images ont etes conjuguees de pos.

        EXAMPLES :
            sage: m = WordMorphism('a->abcde')
            sage: m.conjugate(0) == m
            True
            sage: print m.conjugate(1)
            WordMorphism: a->bcdea
            sage: print m.conjugate(3)
            WordMorphism: a->deabc
            sage: print WordMorphism('').conjugate(4)
            WordMorphism:
            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: print m.conjugate(2)
            WordMorphism: a->cdeab, b->zxy
        """
        return WordMorphism(dict((key, w.conjugate(pos)) for (key, w) in self._morph.iteritems()))

    def has_left_conjugate(self):
        r"""
        Returns true if all the non empty images of self begins with the same letter.

        Retourne vrai si toutes les images non vides commencent par la meme lettre.

        EXAMPLES:
            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.has_left_conjugate()
            False
            sage: WordMorphism('b->xyz').has_left_conjugate()
            True
            sage: WordMorphism('').has_left_conjugate()
            True
            sage: WordMorphism('a->,b->xyz').has_left_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb').has_left_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb,c->').has_left_conjugate()
            True
        """
        I = ifilterfalse(FiniteWord_over_OrderedAlphabet.is_empty, self.images())

        try:
            letter = I.next()[0]
        except StopIteration:
            return True

        #Compare the first letter of all the non empty images
        for image in I:
            if image[0] != letter:
                return False

        return True

    def has_right_conjugate(self):
        r"""
        Returns true if all the non empty images of self ends with the same letter.

        Retourne vrai si toutes les images non vide terminent par la meme lettre.

        EXAMPLES:
            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.has_right_conjugate()
            False
            sage: WordMorphism('b->xyz').has_right_conjugate()
            True
            sage: WordMorphism('').has_right_conjugate()
            True
            sage: WordMorphism('a->,b->xyz').has_right_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb').has_right_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb,c->').has_right_conjugate()
            True
        """
        return self.reversal().has_left_conjugate()

    def list_of_conjugates(self):
        r"""
        Retourne une liste des morphismes conjugues du receveur
        obtenus en conjuguant les prefixes (et suffixes) communs
        de chacune des images.

        Returns the list of all the conjugate morphisms of self obtain
        by taking the conjugate of the common prefixes and suffixes of all
        the images.

        DEFINITION (from [1]):
        Recall from Lothaire [2] (Section 2.3.4) that $\varphi$ is \emph{right conjugate} of $\varphi'$, noted $\varphi\triangleleft\varphi'$, if there exists $u \in \Sigma^*$ such that $$\varphi(\alpha)u = u\varphi'(\alpha), \quad \textrm{for all } \alpha \in \Sigma$$, or equivalently that $\varphi(x)u = u\varphi'(x)$, for all words $x \in \Sigma^*$. Clearly, this relation is not symmetric so that we say that two morphisms $\varphi$ and $\varphi'$ are \emph{conjugate}, noted $\varphi\bowtie\varphi'$, if $\varphi\triangleleft\varphi'$ or $\varphi'\triangleleft\varphi$. It is easy to see that conjugacy of morphisms is an equivalence relation.

        REFERENCES:
            [1] A. Blondin-Mass\'e, S. Brlek, S. Labb\'e, Palindromic lacunas
            of the Thue-Morse word, 2008, to appear.
            [2] M. Lothaire, Algebraic Combinatorics on words, Cambridge
            University Press, 2002.

        EXAMPLES:
            sage: m = WordMorphism('a->abbab,b->abb')
            sage: map(str, m.list_of_conjugates())
            ['WordMorphism: a->babba, b->bab',
            'WordMorphism: a->abbab, b->abb',
            'WordMorphism: a->bbaba, b->bba',
            'WordMorphism: a->babab, b->bab',
            'WordMorphism: a->ababb, b->abb',
            'WordMorphism: a->babba, b->bba',
            'WordMorphism: a->abbab, b->bab']
            sage: m = WordMorphism('a->aaa,b->aa')
            sage: map(str, m.list_of_conjugates())
            ['WordMorphism: a->aaa, b->aa']
            sage: WordMorphism('').list_of_conjugates()
            [Morphism from Words over Ordered Alphabet [] to Words over Ordered Alphabet []]
            sage: m = WordMorphism('a->aba,b->aba')
            sage: map(str, m.list_of_conjugates())
            ['WordMorphism: a->baa, b->baa',
            'WordMorphism: a->aab, b->aab',
            'WordMorphism: a->aba, b->aba']
            sage: m = WordMorphism('a->abb,b->abbab,c->')
            sage: map(str, m.list_of_conjugates())
            ['WordMorphism: a->bab, b->babba, c->',
            'WordMorphism: a->abb, b->abbab, c->',
            'WordMorphism: a->bba, b->bbaba, c->',
            'WordMorphism: a->bab, b->babab, c->',
            'WordMorphism: a->abb, b->ababb, c->',
            'WordMorphism: a->bba, b->babba, c->',
            'WordMorphism: a->bab, b->abbab, c->']
        """
        if self.is_empty():
            return [self]

        #Construire la liste c des morphismes conjugues
        c = []
        m = self
        c.append(m)
        while(m.has_left_conjugate()):
            m = m.conjugate(1)
            if m == self:
                break
            c.append(m)
        m = self
        while(m.has_right_conjugate()):
            m = m.conjugate(-1)
            if m == self:
                break
            c.insert(0, m)

        #Construire la liste d des morphismes distincts
        d = []
        for m in c:
            if m not in d:
                d.append(m)
        return d

    def is_in_classP(self, f=None):
        r"""
        Returns True if self is in class $P$ (or $f$-$P$).

        DEFINITION:
            "[Let $A$ be an alphabet]
            We say that a primitive substitution $S$ is in the \emph{class P}
            if there exists a palindrome $p$ and for each $b\in A$ a palindrome
            $q_b$ such that $S(b)=pq_b$ for all $b\in A$." [1]

            Let $\Sigma$ be an alphabet and $f$ be an involution on $\Sigma$.
            "We say that a morphism $\varphi$ is in {\em class $f$-P if
            there exists an $f$-palindrome $p\in\hbox{$f$-\rm Pal}(\Sigma^*)$
            and for each $\alpha \in \Sigma$ there exists an $f$-palindrome
            $q_\alpha\in\hbox{$f$-\rm Pal}(\Sigma^*)$ such that
            $\varphi(\alpha)=pq_\alpha$." [2]

        INPUT:
            f -- involution (default: None) on the alphabet of self. It must
                 be something that WordMorphism's contructor understand
                 (dict, str, ...).

        OUTPUT:
            boolean -- True if f is None and if self is in class P;
                       True otherwise and if self is in class $f$-P;
                       False otherwise.

        REFERENCES:
            [1] Hof, A., O. Knill et B. Simon, Singular continuous spectrum for
            palindromic Schr\"odinger operators, {\em Commun. Math. Phys.}
            {\bf 174} (1995) 149--159.
            [2] Labbé, Sébastien. Propriétés combinatoires des $f$-palindromes,
            Mémoire de maîtrise en Mathématiques, Montréal, UQAM, 2008, 109 pages.

        EXAMPLES:
            sage: WordMorphism('a->bbaba,b->bba').is_in_classP()
            True
            sage: tm = WordMorphism('a->ab,b->ba')
            sage: tm.is_in_classP()
            False
            sage: f='a->b,b->a'
            sage: tm.is_in_classP(f=f)
            True
            sage: (tm^2).is_in_classP()
            True
            sage: (tm^2).is_in_classP(f=f)
            False
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo.is_in_classP()
            True
            sage: fibo.is_in_classP(f=f)
            False
            sage: (fibo^2).is_in_classP()
            False
            sage: f='a->b,b->a,c->c'
            sage: WordMorphism('a->acbcc,b->acbab,c->acbba').is_in_classP (f)
            True
        """
        if self.is_empty():
            return True

        #Compute the longest common prefix of all the images of letters
        images = self.images()
        lcp = images[0]
        for image in images:
            lcp = lcp.longest_common_prefix(image)

        #Find a common palindrome prefix
        for i in range(len(lcp)+1):
            if lcp[:i].is_palindrome(f=f):

                #If all the suffixes are palindromes,
                for image in images:
                    if not image[i:].is_palindrome(f=f):
                        break
                else:
                    return True

        return False

    def has_conjugate_in_classP(self, f=None):
        r"""
        Returns True if self has a conjugate in class $f$-$P$.

        DEFINITION:
            Let $\Sigma$ be an alphabet and $f$ be an involution on $\Sigma$.
            "We say that a morphism $\varphi$ is in {\em class $f$-P} if
            there exists an $f$-palindrome $p\in\hbox{$f$-\rm Pal}(\Sigma^*)$
            and for each $\alpha \in \Sigma$ there exists an $f$-palindrome
            $q_\alpha\in\hbox{$f$-\rm Pal}(\Sigma^*)$ such that
            $\varphi(\alpha)=pq_\alpha$." [2]

            We say that a morphism $\varphi'$ est de {\em classe $f-P'$}
            if there exists a morphism $\varphi$ conjugate of $\varphi'$
            such that $\varphi$ is in class $f$-P. [2]

        INPUT:
            f -- involution (default: None) on the alphabet of self. It must
                 be something that WordMorphism's contructor understand
                 (dict, str, ...).

        OUTPUT:
            boolean -- True if f is None and if self has a conjugate in class P;
                       True otherwise and if self has a conjugate in class $f$-P;
                       False otherwise.

        REFERENCES:
            [1] Hof, A., O. Knill et B. Simon, Singular continuous spectrum for
            palindromic Schr\"odinger operators, {\em Commun. Math. Phys.}
            {\bf 174} (1995) 149--159.
            [2] Labbé, Sébastien. Propriétés combinatoires des $f$-palindromes,
            Mémoire de maîtrise en Mathématiques, Montréal, UQAM, 2008, 109 pages.

        EXAMPLES:
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo.has_conjugate_in_classP()
            True
            sage: (fibo^2).is_in_classP()
            False
            sage: (fibo^2).has_conjugate_in_classP()
            True
        """
        for k in self.list_of_conjugates():
            if k.is_in_classP(f=f):
                return True
        return False

#    def _latex_(self):
#        r"""
#        Returns latex representation of this morphism.
#
#        EXAMPLES:
#            sage: tm = WordMorphism('a->ab,b->ba'); print tm._latex_()
#            \begin{array}{cccl}
#            \varphi: & [\text{a}, \text{b}]^* \cup [\text{a}, \text{b}]^\omega & \rightarrow & [\text{a}, \text{b}]^* \cup [\text{a}, \text{b}]^\omega \\
#            & \text{a} & \mapsto & \text{ab} \\
#            & \text{b} & \mapsto & \text{ba}
#            \end{array}
#            sage: fibo = WordMorphism('a->ab,b->a'); print (fibo^3)._latex_()
#            \begin{array}{cccl}
#            \varphi: & [\text{a}, \text{b}]^* \cup [\text{a}, \text{b}]^\omega & \rightarrow & [\text{a}, \text{b}]^* \cup [\text{a}, \text{b}]^\omega \\
#            & \text{a} & \mapsto & \text{abaab} \\
#            & \text{b} & \mapsto & \text{aba}
#            \end{array}
#            sage: print (WordMorphism(''))._latex_()
#            \begin{array}{cccl}
#            \varphi: & []^* \cup []^\omega & \rightarrow & []^* \cup []^\omega
#            \end{array}
#            sage: print (WordMorphism('b->baa,a->ab'))._latex_()
#            \begin{array}{cccl}
#            \varphi: & [\text{a}, \text{b}]^* \cup [\text{a}, \text{b}]^\omega & \rightarrow & [\text{a}, \text{b}]^* \cup [\text{a}, \text{b}]^\omega \\
#            & \text{a} & \mapsto & \text{ab} \\
#            & \text{b} & \mapsto & \text{baa}
#            \end{array}
#        """
#        dom = self.domain()
#        cod = self.codomain()
#
#        lines = ["\\varphi: & %s & \\rightarrow & %s" %(latex(dom), latex(cod))]
#        lines += ['& %s & \mapsto & %s' %(latex(b), latex(self(b).string_rep()))  for b in dom.alphabet()]
#        s = " \\\\\n".join(lines)
#
#        return "\\begin{array}{cccl}\n" + s + "\n\\end{array}"


