r"""
Factor-enumerable words

Words having an algorithm that enumerates the factor of length n which
includes finite words and some family of infinite words. This file gathers
methods (e.g. ``rauzy_graph``) that depends only on the existence of such
an algorithm.

AUTHORS:

- Sebastien Labbe (February 24, 2010) : initial version
- Julien Leroy (March 2010): reduced_rauzy_graph

EXAMPLES:

Enumeration of factors::

    sage: w = Word([4,5,6])^7
    sage: it = w.factor_iterator(4)
    sage: it.next()
    word: 6456
    sage: it.next()
    word: 5645
    sage: it.next()
    word: 4564
    sage: it.next()
    Traceback (most recent call last):
    ...
    StopIteration

The set of factors::

    sage: sorted(w.factor_set(3))
    [word: 456, word: 564, word: 645]
    sage: sorted(w.factor_set(4))
    [word: 4564, word: 5645, word: 6456]
    sage: w.factor_set().cardinality()
    61

Rauzy graphs::

    sage: f = words.FibonacciWord()[:30]
    sage: f.rauzy_graph(4)
    Looped digraph on 5 vertices
    sage: f.reduced_rauzy_graph(4)
    Looped multi-digraph on 2 vertices

Left-special and bispecial factors::

    sage: f.number_of_left_special_factors(7)
    1
    sage: f.bispecial_factors()
    [word: , word: 0, word: 010, word: 010010, word: 01001010010]
"""
#*****************************************************************************
#       Copyright (C) 2010 Sebastien Labbe <slabqc@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.combinat.words.abstract_word import Word_class
from itertools import product, chain, tee
from sage.sets.set import Set
from sage.rings.all import Infinity

class Word_nfactor_enumerable(Word_class):
    def factor_iterator(self, n=None):
        r"""
        Returns an iterator over the factor of lenght `n`.

        INPUT:

        -  ``n`` - an integer, or ``None``.

        OUTPUT:

            If ``n`` is an integer, returns an iterator over all distinct
            factors of length ``n``. If ``n`` is ``None``, returns an iterator
            generating all distinct factors.

        .. NOTE:

            This function must be implemented in a more specific class
            (e.g. finite word class).

        EXAMPLES::

            sage: w = Word(range(10))
            sage: sorted(w.factor_iterator(7))
            [word: 0123456, word: 1234567, word: 2345678, word: 3456789]
        """
        msg = 'This function must be defined in a inherited class'
        raise NotImplementedError, msg

    def number_of_factors(self,n):
        r"""
        Return the number of distinct factors of self.

        INPUT:

        -  ``n`` - an integer

        OUTPUT:

            integer

        EXAMPLES::

            sage: w = Word(range(10))
            sage: w.number_of_factors(6)
            5
        """
        return len(list(self.factor_iterator(n)))

    def factor_set(self, n=None):
        r"""
        Returns the set of factors (of length n) of self.

        INPUT:

        - ``n`` - an integer or ``None`` (default: None).

        OUTPUT:

            If ``n`` is an integer, returns the set of all distinct
            factors of length ``n``. If ``n`` is ``None``, returns the set
            of all distinct factors.

        EXAMPLES::

            sage: w = Word('121')
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: 1, word: 12, word: 121, word: 2, word: 21]

        ::

            sage: w = Word('1213121')
            sage: for i in range(w.length()): sorted(w.factor_set(i))
            [word: ]
            [word: 1, word: 2, word: 3]
            [word: 12, word: 13, word: 21, word: 31]
            [word: 121, word: 131, word: 213, word: 312]
            [word: 1213, word: 1312, word: 2131, word: 3121]
            [word: 12131, word: 13121, word: 21312]
            [word: 121312, word: 213121]

        ::

            sage: w = Word([1,2,1,2,3])
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: 1, word: 12, word: 121, word: 1212, word: 12123, word: 123, word: 2, word: 21, word: 212, word: 2123, word: 23, word: 3]

        TESTS::

            sage: w = Word("xx")
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: x, word: xx]

        ::

            sage: Set(Word().factor_set())
            {word: }
        """
        return Set(set(self.factor_iterator(n)))

    def topological_entropy(self, n):
        r"""
        Return the topological entropy for the factors of length n.

        The topological entropy of a sequence `u` is defined as the
        exponential growth rate of the complexity of `u` as the length
        increases: `H_{top}(u)=\lim_{n\to\infty}\frac{\log_d(p_u(n))}{n}`
        where `d` denotes the cardinality of the alphabet and `p_u(n)` is
        the complexity function, i.e. the number of factors of length `n`
        in the sequence `u` [1].

        INPUT:

        - ``self`` - a word defined over a finite alphabet
        -  ``n`` - positive integer

        OUTPUT:

        real number (a symbolic expression)

        EXAMPLES::

            sage: W = Words([0, 1])
            sage: w = W([0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1])
            sage: t = w.topological_entropy(3); t
            1/3*log(7)/log(2)
            sage: n(t)
            0.935784974019201

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: topo = w.topological_entropy
            sage: for i in range(0, 41, 5): print i, n(topo(i), digits=5)
            0 1.0000
            5 0.71699
            10 0.48074
            15 0.36396
            20 0.28774
            25 0.23628
            30 0.20075
            35 0.17270
            40 0.14827

        If no alphabet is specified, an error is raised::

            sage: w = Word(range(20))
            sage: w.topological_entropy(3)
            Traceback (most recent call last):
            ...
            TypeError: The word must be defined over a finite alphabet

        The following is ok::

            sage: W = Words(range(20))
            sage: w = W(range(20))
            sage: w.topological_entropy(3)
            1/3*log(18)/log(20)

        REFERENCES:

           [1] N. Pytheas Fogg, Substitutions in Dynamics, Arithmetics,
           and Combinatorics, Lecture Notes in Mathematics 1794, Springer
           Verlag. V. Berthe, S. Ferenczi, C. Mauduit and A. Siegel, Eds.
           (2002).
        """
        d = self.parent().size_of_alphabet()
        if d is Infinity:
            raise TypeError, "The word must be defined over a finite alphabet"
        if n == 0:
            return 1
        pn = self.number_of_factors(n)
        from sage.functions.all import log
        return log(pn, base=d)/n

    @cached_method
    def rauzy_graph(self, n):
        r"""
        Returns the Rauzy graph of the factors of length n of self.

        The vertices are the factors of length `n` and there is an edge from
        `u` to `v` if `ua = bv` is a factor of length `n+1` for some letters
        `a` and `b`.

        INPUT:

        - ``n`` - integer

        EXAMPLES::

            sage: w = Word(range(10)); w
            word: 0123456789
            sage: g = w.rauzy_graph(3); g
            Looped digraph on 8 vertices
            sage: WordOptions(identifier='')
            sage: g.vertices()
            [012, 123, 234, 345, 456, 567, 678, 789]
            sage: g.edges()
            [(012, 123, 3),
             (123, 234, 4),
             (234, 345, 5),
             (345, 456, 6),
             (456, 567, 7),
             (567, 678, 8),
             (678, 789, 9)]
            sage: WordOptions(identifier='word: ')

        ::

            sage: f = words.FibonacciWord()[:100]
            sage: f.rauzy_graph(8)
            Looped digraph on 9 vertices

        ::

            sage: w = Word('1111111')
            sage: g = w.rauzy_graph(3)
            sage: g.edges()
            [(word: 111, word: 111, word: 1)]

        ::

            sage: w = Word('111')
            sage: for i in range(5) : w.rauzy_graph(i)
            Looped multi-digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 0 vertices

        Multi-edges are allowed for the empty word::

            sage: W = Words('abcde')
            sage: w = W('abc')
            sage: w.rauzy_graph(0)
            Looped multi-digraph on 1 vertex
            sage: _.edges()
            [(word: , word: , word: a),
             (word: , word: , word: b),
             (word: , word: , word: c)]
        """
        from sage.graphs.digraph import DiGraph
        multiedges = True if n == 0 else False
        g = DiGraph(loops=True, multiedges=multiedges)
        if n == self.length():
            g.add_vertex(self)
        else:
            for w in self.factor_iterator(n+1):
                u = w[:-1]
                v = w[1:]
                a = w[-1:]
                g.add_edge(u,v,a)
        return g

    def reduced_rauzy_graph(self, n):
        r"""
        Returns the reduced Rauzy graph of order `n` of self.

        INPUT:

        - ``n`` - non negative integer. Every vertex of a reduced
          Rauzy graph of order `n` is a factor of length `n` of self.

        OUTPUT:

        Looped multi-digraph

        DEFINITION:

        For infinite periodic words (resp. for finite words of type `u^i
        u[0:j]`), the reduced Rauzy graph of order `n` (resp. for `n`
        smaller or equal to `(i-1)|u|+j`) is the directed graph whose
        unique vertex is the prefix `p` of length `n` of self and which has
        an only edge which is a loop on `p` labelled by `w[n+1:|w|] p`
        where `w` is the unique return word to `p`.

        In other cases, it is the directed graph defined as followed.  Let
        `G_n` be the Rauzy graph of order `n` of self. The vertices are the
        vertices of `G_n` that are either special or not prolongable to the
        right or to the left. For each couple (`u`, `v`) of such vertices
        and each directed path in `G_n` from `u` to `v` that contains no
        other vertices that are special, there is an edge from `u` to `v`
        in the reduced Rauzy graph of order `n` whose label is the label of
        the path in `G_n`.

        .. NOTE::

            In the case of infinite recurrent non periodic words, this
            definition correspond to the following one that can be found in
            [1] and [2]  where a simple path is a path that begins with a
            special factor, ends with a special factor and contains no
            other vertices that are special:

            The reduced Rauzy graph of factors of length `n` is obtained
            from `G_n` by replacing each simple path `P=v_1 v_2 ...
            v_{\ell}` with an edge `v_1 v_{\ell}` whose label is the
            concatenation of the labels of the edges of `P`.

        EXAMPLES::

            sage: w = Word(range(10)); w
            word: 0123456789
            sage: g = w.reduced_rauzy_graph(3); g
            Looped multi-digraph on 2 vertices
            sage: g.vertices()
            [word: 012, word: 789]
            sage: g.edges()
            [(word: 012, word: 789, word: 3456789)]

        For the Fibonacci word::

            sage: f = words.FibonacciWord()[:100]
            sage: g = f.reduced_rauzy_graph(8);g
            Looped multi-digraph on 2 vertices
            sage: g.vertices()
            [word: 01001010, word: 01010010]
            sage: g.edges()
            [(word: 01001010, word: 01010010, word: 010), (word: 01010010, word: 01001010, word: 01010), (word: 01010010, word: 01001010, word: 10)]

        For periodic words::

            sage: from itertools import cycle
            sage: w = Word(cycle('abcd'))[:100]
            sage: g = w.reduced_rauzy_graph(3)
            sage: g.edges()
            [(word: abc, word: abc, word: dabc)]

        ::

            sage: w = Word('111')
            sage: for i in range(5) : w.reduced_rauzy_graph(i)
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped multi-digraph on 1 vertex
            Looped multi-digraph on 0 vertices

        For ultimately periodic words::

            sage: sigma = WordMorphism('a->abcd,b->cd,c->cd,d->cd')
            sage: w = sigma.fixed_point('a')[:100]; w
            word: abcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd...
            sage: g = w.reduced_rauzy_graph(5)
            sage: g.vertices()
            [word: abcdc, word: cdcdc]
            sage: g.edges()
            [(word: abcdc, word: cdcdc, word: dc), (word: cdcdc, word: cdcdc, word: dc)]

        AUTHOR:

        Julien Leroy (March 2010): initial version

        REFERENCES:

        - [1] M. Bucci et al.  A. De Luca, A. Glen, L. Q. Zamboni, A
          connection between palindromic and factor complexity using
          return words," Advances in Applied Mathematics 42 (2009) 60-74.

        - [2] L'ubomira Balkova, Edita Pelantova, and Wolfgang Steiner.
          Sequences with constant number of return words. Monatsh. Math,
          155 (2008) 251-263.
        """
        from sage.graphs.all import DiGraph
        from copy import copy
        g = copy(self.rauzy_graph(n))
        # Otherwise it changes the rauzy_graph function.
        l = [v for v in g if g.in_degree(v)==1 and g.out_degree(v)==1]
        if g.num_verts() !=0 and len(l)==g.num_verts():
            # In this case, the Rauzy graph is simply a cycle.
            g = DiGraph()
            g.allow_loops(True)
            g.add_vertex(self[:n])
            g.add_edge(self[:n],self[:n],self[n:n+len(l)])
        else:
            g.allow_loops(True)
            g.allow_multiple_edges(True)
            for v in l:
                [i] = g.neighbors_in(v)
                [o] = g.neighbors_out(v)
                g.add_edge(i,o,g.edge_label(i,v)[0]*g.edge_label(v,o)[0])
                g.delete_vertex(v)
        return g

    def left_special_factors_iterator(self, n=None):
        r"""
        Returns an iterator over the left special factors (of length n).

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           an iterator over all left special factors.

        EXAMPLES::

            sage: alpha, beta, x = 0.54, 0.294, 0.1415
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: sorted(w.left_special_factors_iterator(3))
            [word: 000, word: 010]
            sage: sorted(w.left_special_factors_iterator(4))
            [word: 0000, word: 0101]
            sage: sorted(w.left_special_factors_iterator(5))
            [word: 00000, word: 01010]
        """
        if n is None:
            for i in range(self.length()):
                for w in self.left_special_factors_iterator(i):
                    yield w
        else:
            g = self.rauzy_graph(n)
            in_d = g.in_degree
            for v in g:
                if in_d(v) > 1:
                    yield v

    def left_special_factors(self, n=None):
        r"""
        Returns the left special factors (of length n).

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           all left special factors.

        OUTPUT:

        A list of words.

        .. WARNING::

            This may not halt for infinite words having an infinite number
            of such factors. Use the iterator version of this function
            instead.

        EXAMPLES::

            sage: alpha, beta, x = 0.54, 0.294, 0.1415
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: for i in range(5): print i, sorted(w.right_special_factors(i))
            0 [word: ]
            1 [word: 0]
            2 [word: 00, word: 10]
            3 [word: 000, word: 010]
            4 [word: 0000, word: 1010]
        """
        return list(self.left_special_factors_iterator(n))

    def right_special_factors_iterator(self, n=None):
        r"""
        Returns an iterator over the right special factors (of length n).

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           an iterator over all right special factors.

        EXAMPLES::

            sage: alpha, beta, x = 0.61, 0.54, 0.3
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: sorted(w.right_special_factors_iterator(3))
            [word: 010, word: 101]
            sage: sorted(w.right_special_factors_iterator(4))
            [word: 0101, word: 1010]
            sage: sorted(w.right_special_factors_iterator(5))
            [word: 00101, word: 11010]
        """
        if n is None:
            for i in range(self.length()):
                for w in self.right_special_factors_iterator(i):
                    yield w
        else:
            g = self.rauzy_graph(n)
            out_d = g.out_degree
            for v in g:
                if out_d(v) > 1:
                    yield v

    def right_special_factors(self, n=None):
        r"""
        Returns the right special factors (of length n).

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           all right special factors.

        OUTPUT:

        A list of words.

        .. WARNING::

            This may not halt for infinite words having an infinite number
            of such factors. Use the iterator version of this function
            instead.

        EXAMPLES::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(5): print i, sorted(w.right_special_factors(i))
            0 [word: ]
            1 [word: 0, word: 1]
            2 [word: 01, word: 10]
            3 [word: 001, word: 010, word: 101, word: 110]
            4 [word: 0110, word: 1001]
        """
        return list(self.right_special_factors_iterator(n))

    def bispecial_factors_iterator(self, n=None):
        r"""
        Returns an iterator over the bispecial factors (of length n).

        A factor `u` of a word `w` is *bispecial* if it is right special
        and left special.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           an iterator over all bispecial factors.

        EXAMPLES::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(10):
            ...     for u in sorted(w.bispecial_factors_iterator(i)):
            ...         print i,u
            0 word:
            1 word: 0
            1 word: 1
            2 word: 01
            2 word: 10
            3 word: 010
            3 word: 101
            4 word: 0110
            4 word: 1001
            6 word: 011001
            6 word: 100110
            8 word: 10010110

        ::

            sage: for u in sorted(w.bispecial_factors_iterator(), key=lambda u:(len(u),u)): print u
            word:
            word: 0
            word: 1
            word: 01
            word: 10
            word: 010
            word: 101
            word: 0110
            word: 1001
            word: 011001
            word: 100110
            word: 10010110
        """
        if n is None:
            for i in range(self.length()):
                for w in self.bispecial_factors_iterator(i):
                    yield w
        else:
            g = self.rauzy_graph(n)
            in_d = g.in_degree
            out_d = g.out_degree
            for v in g:
                if out_d(v) > 1 and in_d(v) > 1:
                    yield v

    def bispecial_factors(self, n=None):
        r"""
        Returns the bispecial factors (of length n).

        A factor `u` of a word `w` is *bispecial* if it is right special
        and left special.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           all bispecial factors.

        OUTPUT:

        A list of words.

        .. WARNING::

            This may not halt for infinite words having an infinite number
            of such factors. Use the iterator version of this function
            instead.

        EXAMPLES::

            sage: w = words.FibonacciWord()[:30]
            sage: w.bispecial_factors()
            [word: , word: 0, word: 010, word: 010010, word: 01001010010]

        ::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(10): print i, sorted(w.bispecial_factors(i))
            0 [word: ]
            1 [word: 0, word: 1]
            2 [word: 01, word: 10]
            3 [word: 010, word: 101]
            4 [word: 0110, word: 1001]
            5 []
            6 [word: 011001, word: 100110]
            7 []
            8 [word: 10010110]
            9 []
        """
        return list(self.bispecial_factors_iterator(n))

    def number_of_left_special_factors(self, n):
        r"""
        Returns the number of left special factors of length n.

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer

        OUTPUT:

        Non negative integer

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.number_of_left_special_factors(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_left_special_factors(i) for i in range(10)]
            [1, 2, 2, 4, 2, 4, 4, 2, 2, 4]
        """
        L = self.rauzy_graph(n).in_degree()
        return sum(1 for i in L if i>1)

    def number_of_right_special_factors(self, n):
        r"""
        Returns the number of right special factors of length n.

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer

        OUTPUT:

        Non negative integer

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.number_of_right_special_factors(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_right_special_factors(i) for i in range(10)]
            [1, 2, 2, 4, 2, 4, 4, 2, 2, 4]
        """
        L = self.rauzy_graph(n).out_degree()
        return sum(1 for i in L if i>1)

