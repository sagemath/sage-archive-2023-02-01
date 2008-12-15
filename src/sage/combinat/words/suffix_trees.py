# Suffix Tries and Suffix Trees
#*****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import chain
from sage.structure.sage_object import SageObject
from sage.graphs.graph import Graph, DiGraph
from sage.graphs.graph_fast import spring_layout_fast
from sage.misc.flatten import flatten
from sage.sets.set import Set
from sage.combinat.words.word import Words, FiniteWord_over_OrderedAlphabet
from sage.combinat.words.word_content import BuildWordContent
from sage.combinat.words.utils import reverse_map, isint

################################################################################
# Suffix Tries
################################################################################

# new unique object for the end of string symbol.
end_of_string = object()

class SuffixTrie(SageObject):
    def __init__(self, word):
        r"""
        Construct the suffix trie of the word w.

        The suffix trie of a finite word w is a data structure representing
        the factors of w. It is a tree whose edges are labelled with
        letters of w, and whose leafs correspond to suffixes of w.

        This is a straightforward implementation of Algorithm 1 from [1].
        It constructs the suffix trie of w[:i] from that of w[:i-1].

        A suffix trie is modelled as a deterministic finite-state automaton
        together with the suffix_link map. The set of states corresponds to
        factors of the word (below we write x' for the state corresponding
        to x); these are always 0, 1, .... The state 0 is the initial
        state, and it corresponds to the empty word.  For the purposes of
        the algorithm, there is also an auxiliary state -1. The transition
        function t is defined as
                t(-1,a) = 0 for all letters a; and
                t(x',a) = y' for all x',y' \in Q such that y = xa,
        and the suffix link function is defined as
                suffix_link(0) = -1;
                suffix_link(x') = y', if x = ay for some letter a.

        REFERENCES:
            [1] E. Ukkonen, "On-line construction of suffix trees",
            Algorithmica, 1995, volume 14, number 3, pages 249--260.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w); t
            Suffix Trie of the word: cacao

            sage: e = Words("ab")(format="empty")
            sage: t = SuffixTrie(e); t
            Suffix Trie of the word:
            sage: t.process_letter("a"); t
            Suffix Trie of the word: a
            sage: t.process_letter("b"); t
            Suffix Trie of the word: ab
            sage: t.process_letter("a"); t
            Suffix Trie of the word: aba

        TESTS:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: s = SuffixTrie(w)
            sage: loads(dumps(s))
            Suffix Trie of the word: cacao
        """
        # Create the suffix trie for the empty word.
        self._active_state = 0
        self._transition_function = {}
        self._suffix_link = [-1]
        self._alphabet = word.alphabet()

        # Process each letter, in order.
        W = Words(word.alphabet())
        w = W()
        for letter in word:
            self._process_letter(W([letter]))

    def _process_letter(self,letter):
        r"""
        Process a letter. That is, modify the current suffix trie producing
        the suffix trie for self.word() + letter.

        NOTE:
            letter must occur within the alphabet of the word.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: t = SuffixTrie(Word("ababba"))
            sage: t._process_letter(Words("ab")("b")); t
            Suffix Trie of the word: ababbab
        """
        r = self._active_state
        # While r is not the auxiliary vertex, or
        # there is not transition from r along letter, ...
        while r != -1 and \
                not self._transition_function.has_key((r,letter)):
            # adjoin a new state s
            s = len(self._suffix_link)
            self._suffix_link.append(None)
            # create a transition from r to s along letter
            self._transition_function[(r,letter)] = s
            if r != self._active_state:
                # update the suffix link
                self._suffix_link[old_s] = s
            old_s = s
            r = self._suffix_link[r]
        # update the suffix link for the last visited state
        if r == -1:
            self._suffix_link[old_s] = 0
        else:
            self._suffix_link[old_s] = self._transition_function[(r,letter)]
        # update the active state
        self._active_state = \
                self._transition_function[(self._active_state, letter)]

    def process_letter(self,letter):
        r"""
        Modify self to produce the suffix trie for self.word() + letter.

        NOTE:
            letter must occur within the alphabet of the word.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("ab")("ababba")
            sage: t = SuffixTrie(w); t
            Suffix Trie of the word: ababba
            sage: t.process_letter("a"); t
            Suffix Trie of the word: ababbaa

        TESTS:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w); t
            Suffix Trie of the word: cacao
            sage: t.process_letter("d")
            Traceback (most recent call last):
                ...
            IndexError: letter not in alphabet: 'd'
        """
        # Make certain that letter is a word containing one letter.
        letter = Words(self._alphabet)([letter])
        self._process_letter(letter)

    #####
    # The following are not necessary for constructing the suffix trie (just
    # the __init__ and process_letter are needed). They just add additional
    # functionality to the class.
    #####

    def _repr_(self):
        """
        TESTS:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: SuffixTrie(Word("abcba"))._repr_()
            'Suffix Trie of the word: abcba'
        """
        return 'Suffix Trie of the %s' % self.word()

    def node_to_word(self, state=0):
        r"""
        Returns the word obtained by reading the edge labels from 0 to
        state.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("abc")("abcba")
            sage: t = SuffixTrie(w)
            sage: t.node_to_word(10)
            word: abcba
            sage: t.node_to_word(7)
            word: abcb
        """
        if state == 0:
            return Words(self._alphabet)()
        # We first invert the transition function
        tf_inv = reverse_map(self._transition_function)
        # Starting from the active state,
        # read labels along the unique path to the root.
        (u,letter) = tf_inv[state]
        w = letter
        s = u
        while s != 0:
            (u,letter) = tf_inv[s]
            w = letter * w
            s = u
        return w

    def word(self):
        r"""
        Returns the word whose suffix tree this is.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("abc")("abcba")
            sage: t = SuffixTrie(w)
            sage: t.word()
            word: abcba
            sage: t.word() == w
            True
        """
        return self.node_to_word(self._active_state)

    def __eq__(self,other):
        r"""
        If self and other have the same transition function, the same
        suffix link, and the same word, then they are equal.

        TEST:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: SuffixTrie(Word("cacao")) == SuffixTrie(Word("ababc"))
            False
            sage: W = Words("cao")
            sage: s = SuffixTrie(W("cacao"))
            sage: t = SuffixTrie(W())
            sage: t.process_letter("c")
            sage: t.process_letter("a")
            sage: t.process_letter("c")
            sage: t.process_letter("a")
            sage: t.process_letter("o")
            sage: t == s
            True
        """
        if not isinstance(other,SuffixTrie):
            return False
        return self._transition_function == other._transition_function \
            and self._suffix_link == other._suffix_link \
            and self.word() == other.word()

    def transition_function(self, node, word):
        r"""
        Returns the state reached by beginning at node and following the
        arrows in the transition graph labelled by the letters of word.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words([0,1])([0,1,0,1,1])
            sage: t = SuffixTrie(w)
            sage: [t.transition_function(u,letter) == v \
                    for ((u,letter),v) in t._transition_function.iteritems()] \
                    == [True] * len(t._transition_function)
            True
        """
        if node == -1:
            return self.transition_function(0, word[1:])
        if len(word) == 0:
            return 0
        if len(word) == 1:
            return self._transition_function[(node,word)]
        else:
            return self.transition_function( \
                    self._transition_function[(node,word[0:1])], word[1:])

    def states(self):
        r"""
        Returns the states of the automaton defined by the suffix trie.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words([0,1])([0,1,1])
            sage: t = SuffixTrie(w)
            sage: t.states()
            [0, 1, 2, 3, 4]

            sage: u = Words("aco")("cacao")
            sage: s = SuffixTrie(u)
            sage: s.states()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        """
        return range(len(self._transition_function))

    def suffix_link(self, state):
        r"""
        Evaluates the suffix link map of the suffix trie on state.
        Note that the suffix link map is not defined on -1.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: map(t.suffix_link, range(13))
            [-1, 0, 3, 0, 5, 1, 7, 2, 9, 10, 11, 12, 0]
            sage: t.suffix_link(0)
            -1


        TESTS:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: t.suffix_link([1])
            Traceback (most recent call last):
                ...
            TypeError: [1] is not an integer
            sage: t.suffix_link(-1)
            Traceback (most recent call last):
                ...
            TypeError: suffix link is not defined for -1
            sage: t.suffix_link(17)
            Traceback (most recent call last):
                ...
            TypeError: 17 is not a state
        """
        if not isint(state):
            raise TypeError, "%s is not an integer" % state
        if state == -1:
            raise TypeError, "suffix link is not defined for -1"
        if state not in range(len(self._suffix_link)):
            raise TypeError, "%s is not a state" % state
        return self._suffix_link[state]

    def active_state(self):
        r"""
        Returns the active state of the suffix trie. This is the state
        corresponding to the word as a suffix of itself.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: t.active_state()
            8

            sage: u = Words([0,1])([0,1,1,0,1,0,0,1])
            sage: s = SuffixTrie(u)
            sage: s.active_state()
            22
        """
        return self._active_state

    def final_states(self):
        r"""
        Returns the set of final states of the suffix trie. These are the
        states corresponding to the suffixes of self.word(). They are
        obtained be repeatedly following the suffix link from the active
        state until we reach 0.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: t.final_states() == Set([8, 9, 10, 11, 12, 0])
            True
        """
        s = self._active_state
        F = [s]
        while s != 0:
            s = self._suffix_link[s]
            F.append(s)
        return Set(F)

    def has_suffix(self,word):
        r"""
        Return True if and only if word is a suffix of self.word().

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: [t.has_suffix(w[i:]) for i in range(len(w)+1)]
            [True, True, True, True, True, True]
            sage: [t.has_suffix(w[:i]) for i in range(len(w)+1)]
            [True, False, False, False, False, True]
        """
        # Find the state corresponding to word, and
        # check to see if s is a final state.
        s = self.transition_function(0, word)
        q = self._active_state
        if q == s:
            return True
        else:
            while q != 0:
                q = self._suffix_link[q]
                if q == s:
                    return True
        return False

    def to_digraph(self):
        r"""
        Returns a DiGraph object of the transition graph of the suffix trie.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cac")
            sage: t = SuffixTrie(w)
            sage: d = t.to_digraph(); d
            Digraph on 6 vertices
            sage: d.adjacency_matrix()
            [0 1 0 1 0 0]
            [0 0 1 0 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
            [0 0 0 0 0 0]
            [0 0 0 0 0 0]
        """
        dag = {}
        for ((u,letter),v) in self._transition_function.iteritems():
            dag.setdefault(u, {})[v] = letter
        return DiGraph(dag)

    def plot(self, layout='tree', tree_root=0, tree_orientation='up',
            vertex_colors=None, edge_labels=True, *args, **kwds):
        r"""
        Returns a Graphics object corresponding to the transition graph of
        the suffix trie.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: SuffixTrie(Word("cacao")).plot()

        TESTS:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: type(SuffixTrie(Word("cacao")).plot())
            <class 'sage.plot.plot.Graphics'>
        """
        tree = self.to_digraph()
        for (u,v,label) in tree.edge_iterator():
            tree.set_edge_label(u, v, label.string_rep())
        if vertex_colors is None:
            suffix_nodes = self.final_states()
            non_suffix_nodes = list(set(self.states()) - set(suffix_nodes))
            vertex_colors = {'#fec7b8':suffix_nodes,'#ffffff':non_suffix_nodes}
        return tree.plot(layout=layout, tree_root=tree_root,
                tree_orientation=tree_orientation,
                vertex_colors=vertex_colors, edge_labels=edge_labels,
                *args, **kwds)

    def show(self, *args, **kwds):
        r"""
        Displays the output of self.plot().

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cac")
            sage: t = SuffixTrie(w)
            sage: t.show()
        """
        self.plot(*args, **kwds).show()
        return

################################################################################
# Suffix Trees
################################################################################

class ImplicitSuffixTree(SageObject):
    def __init__(self, word):
        r"""
        Construct the implicit suffix tree of a word w.

        The suffix tree of a word w is a compactification of the suffix
        trie for w. The compactification removes all nodes that have
        exactly one incoming edge and exactly one outgoing edge. It
        consists of two components: a tree and a word. Thus, instead of
        labelling the edges by factors of w, we can labelled them by
        indices of the occurrence of the factors in w.

        The following is a straightforward implementation of Ukkonen's
        on-line algorithm for constructing the implicit suffix tree [1].
        It constructs the suffix tree for w[:i] from that of w[:i-1].

        GENERAL IDEA. The suffix tree of w[:i+1] can be obtained from that
        of w[:i] by visiting each node corresponding to a suffix of w[:i]
        and modifying the tree by applying one of two rules (either append
        a new node to the tree, or split an edge into two). The "active
        state" is the node where the algorithm begins and the "suffix link"
        carries us to the next node that needs to be dealt with.

        TREE. The tree is modelled as an automaton, which is stored as a
        dictionary of dictionaries: it is keyed by the nodes of the tree,
        and the corresponding dictionary is keyed by pairs (i,j) of
        integers representing the word w[i-1:j]. This makes it faster to
        look up a particular transition beginning at a specific node.

        STATES/NODES. The states will always be -1, 0, 1, ..., n. The state
        -1 is special and is only used for the purposes of the algorithm.
        All transitions map -1 to 0, so this information is not explicitly
        stored in the transition function.

        EXPLICIT/IMPLICIT NODES. By definition, some of the nodes will not
        be states, but merely locations along an edge; these are called
        implicit nodes. A node r (implicit or explicit) is referenced as a
        pair (s,(k,p)) where s is an ancestor of r and w[k-1:p] is the word
        read by transitioning from s to r in the suffix trie. A reference
        pair is canonical if s is the closest ancestor of r.

        SUFFIX LINK. The algorithm makes use of a map from (some) nodes to
        other nodes, called the suffix link. This is stored as a
        dictionary.

        ACTIVE STATE. We store as ._active_state the active state of the
        tree, the state where the algorithm will begin when processing the
        next letter.

        RUNNING TIME. The running time and storage space of the algorithm
        is linear in the length of the word w (whereas for a suffix tree it
        is quadratic).

        REFERENCES:
            [1] E. Ukkonen, "On-line construction of suffix trees",
            Algorithmica, 1995, volume 14, number 3, pages 249--260.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("cacao")
            sage: t = ImplicitSuffixTree(w); t
            Implicit Suffix Tree of the word: cacao
            sage: ababb = Words([0,1])([0,1,0,1,1])
            sage: s = ImplicitSuffixTree(ababb); s
            Implicit Suffix Tree of the word: 01011

        TESTS:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("cao")("cacao")
            sage: s = ImplicitSuffixTree(w)
            sage: loads(dumps(s))
            Implicit Suffix Tree of the word: cacao
        """
        # For constructing the suffix tree.
        self._transition_function = {0:{}}
        self._suffix_link = {0:-1}
        self._active_state = (0,(1,1))
        self._word_content = []
        for letter in word._word_content:
            self._word_content.append(letter)
            self._process_letter(letter)
        # _word is not needed for constructing the suffix tree,
        # but it is useful for the other methods.
        self._word = word

    def _process_letter(self, letter):
        r"""
        This is the main part of Ukkonen's algorithm. This corresponds to
        the algorithm "update" in [1].

        REFERENCES:
            [1] E. Ukkonen, "On-line construction of suffix trees",
            Algorithmica, 1995, volume 14, number 3, pages 249--260.

        TESTS:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("caca")
            sage: t = ImplicitSuffixTree(w); t
            Implicit Suffix Tree of the word: caca
            sage: new_letter = w.alphabet().rank("o")
            sage: t._word_content.append(new_letter)
            sage: t._process_letter(new_letter)
            sage: t
            Implicit Suffix Tree of the word: cacao

            sage: W = Words([0,1])
            sage: s = ImplicitSuffixTree(W([0,1,0,1])); s
            Implicit Suffix Tree of the word: 0101
            sage: new_letter = W.alphabet().rank(1)
            sage: s._word_content.append(new_letter)
            sage: s._process_letter(new_letter)
            sage: s
            Implicit Suffix Tree of the word: 01011
        """
        (s,(k,i)) = self._active_state
        old_r = 0
        (end_state, r) = self._test_and_split(s,(k,i-1),letter)
        while end_state == False:
            # adjoin a new state rr and create a transition from r to rr
            rr = len(self._transition_function)
            self._transition_function[rr] = {}
            self._transition_function[r][(i,None)] = rr
            # update the suffix link, if necessary
            if old_r != 0:
                self._suffix_link[old_r] = r
            old_r = r
            # follow the suffix link to the next state
            (s, k) = self._canonize(self._suffix_link[s], (k,i-1))
            (end_state, r) = self._test_and_split(s, (k,i-1), letter)
        # update the suffix link, if necessary
        if old_r != 0:
            self._suffix_link[old_r] = s
        # set the active state
        (s,k) = self._canonize(s,(k,i))
        self._active_state = (s, (k, i+1))
        return

    def _test_and_split(self, s, (k, p), letter):
        r"""
        Helper function for _process_letter. Tests to see whether an edge
        needs to be split. Returns (True, state), where state is the next
        state to process (either a newly created state or the original s).

        TESTS:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("caca")
            sage: t = ImplicitSuffixTree(w)
            sage: t._word_content.append(w.alphabet().rank("o"))
            sage: t._test_and_split(0, (4,5), w.alphabet().rank("o"))
            (False, 3)
        """
        if k <= p:
            # find the transition from s that begins with k-th letter
            ((kk,pp), ss) = self._find_transition(s, self._word_content[k-1])
            if letter == self._word_content[kk + p - k]:
                return (True, s)
            else:
                # replace transition above by transitions
                del self._transition_function[s][(kk,pp)]
                r = len(self._transition_function)
                self._transition_function[r] = {}
                self._transition_function[s][(kk, kk+p-k)] = r
                self._transition_function[r][(kk+p-k+1, pp)] = ss
                return (False, r)
        else:
            transition = self._find_transition(s, letter)
            if transition is None:
                return (False, s)
            else:
                return (True, s)

    def _canonize(self, s, (k, p)):
        r"""
        Given an implicit or explicit reference pair for a node, returns
        the canonical reference pair.

        Recall that a node r is referenced as (s, (k,p)), where s is an
        ancestor or r and w[k-1:p] is the word obtained by reading the edge
        labels along the path from s to r. A reference pair is canonical if
        s is the closest ancestor of r.

        TESTS:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: t = ImplicitSuffixTree(Word("cacao"))
            sage: t._canonize(0,(3,5))
            (3, 5)
            sage: t._canonize(0,(2,5))
            (5, 3)
        """
        if p < k:
            return (s,k)
        else:
            ((kk,pp), ss) = self._find_transition(s, self._word_content[k-1])
            while pp is not None and pp - kk <= p - k:
                k = k + pp - kk + 1
                s = ss
                if k <= p:
                    ((kk,pp), ss) = self._find_transition(s, self._word_content[k-1])
            return (s, k)

    def _find_transition(self, state, letter):
        r"""
        Returns the transition from state that begins with letter. Returns
        None if no such transition exists.

        The transitions are stored as a dictionary of dictionaries: keyed
        by the nodes, with the corresponding dictionary keyed by pairs
        (i,j) of integers representing the word w[i-1:j].

            ._transition_function = {..., node: {(i,j): target_node, ...} }

        TESTS:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: t = ImplicitSuffixTree(Word("cacao"))
            sage: t._find_transition(-1, 1)
            ((0, 0), 0)
            sage: t._find_transition(0, 0)
            ((2, 2), 5)
            sage: t._find_transition(0, 1)
            ((1, 2), 3)
            sage: t._find_transition(5, 1)
            ((3, None), 2)
            sage: t._find_transition(5, 0)

            sage: t = ImplicitSuffixTree(Word([0,1,0,1,1]))
            sage: t._find_transition(3, 1)
            ((5, None), 4)
        """
        if state == -1:
            return ((0, 0), 0)
        else:
            if self._transition_function.has_key(state):
                for ((k,p),s) in self._transition_function[state].iteritems():
                    if self._word_content[k-1] == letter:
                        return ((k,p), s)
            return None

    #####
    # The following are not necessary for constructing the implicit suffix
    # tree; they add additional functionality to the class.
    #####

    #####
    # Visualization
    #####

    def _repr_(self):
        r"""
        TESTS:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: ImplicitSuffixTree(Word("abcba"))._repr_()
            'Implicit Suffix Tree of the word: abcba'
        """
        return 'Implicit Suffix Tree of the %s' % self.word()

    def word(self):
        r"""
        Returns the word whose implicit suffix tree this is.

        TEST:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: ImplicitSuffixTree(Word([0,1,0,1,0])).word() == Word([0,1,0,1,0])
            True
        """
        return self._word.parent()(map(self._word.alphabet().unrank, self._word_content))

    def transition_function_dictionary(self):
        r"""
        Returns the transition function as a dictionary of dictionaries.
        The format is consistent with the input format for DiGraph.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words("aco")
            sage: t = ImplicitSuffixTree(W("cac"))
            sage: t.transition_function_dictionary()
            {0: {1: (0, None), 2: (1, None)}}

            sage: W = Words([0,1])
            sage: t = ImplicitSuffixTree(W([0,1,0]))
            sage: t.transition_function_dictionary()
            {0: {1: (0, None), 2: (1, None)}}
        """
        d = {}
        for (u,v,(i,j)) in self.edge_iterator():
            d.setdefault(u, {})[v] = (i,j)
        return d

    def to_digraph(self, word_labels=False):
        r"""
        Returns a DiGraph object of the transition graph of the suffix tree.

            word_labels -- if False, labels the edges by pairs (i, j);
                           if True, labels the edges by word[i:j].

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.to_digraph()
            Digraph on 8 vertices
        """
        if self._word_content == []:
            d = {0:{}}
            return DiGraph(d)
        d = self.transition_function_dictionary()
        for u in d:
            for (v,(i,j)) in d[u].iteritems():
                if word_labels:
                    d[u][v] = self._word[i:j]
                elif j == None:
                    d[u][v] = (i,len(self._word_content))
        return DiGraph(d)

    def plot(self, word_labels=False, layout='tree', tree_root=0,
            tree_orientation='up', vertex_colors=None, edge_labels=True,
            *args, **kwds):
        r"""
        Returns a Graphics object corresponding to the transition graph of
        the suffix tree.

        INPUT:
            word_labels -- if False, labels the edges by pairs (i, j);
                           if True, labels the edges by word[i:j].

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: ImplicitSuffixTree(Word('cacao')).plot(word_labels=True)
            sage: ImplicitSuffixTree(Word('cacao')).plot(word_labels=False)

        TESTS:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: type(ImplicitSuffixTree(Word('cacao')).plot(word_labels=True))
            <class 'sage.plot.plot.Graphics'>
            sage: type(ImplicitSuffixTree(Word('cacao')).plot(word_labels=False))
            <class 'sage.plot.plot.Graphics'>
        """
        tree = self.to_digraph(word_labels=word_labels)
        if word_labels:
            for (u,v,label) in tree.edge_iterator():
                tree.set_edge_label(u, v, label.string_rep())
        if vertex_colors is None:
            veretex_colors = {'#fec7b8':tree.vertices()}
        return tree.plot(layout=layout, tree_root=tree_root,
                tree_orientation=tree_orientation,
                vertex_colors=vertex_colors, edge_labels=edge_labels,
                *args, **kwds)

    def show(self, word_labels=None, *args, **kwds):
        r"""
        Displays the output of self.plot().

        INPUT:
            word_labels -- if False, labels the edges by pairs (i, j);
                           if True, labels the edges by word[i:j].

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("cao")("cacao")
            sage: t = ImplicitSuffixTree(w)
            sage: t.show(word_labels=True)
            sage: t.show(word_labels=False)
        """
        self.plot(word_labels=word_labels, *args, **kwds).show()
        return


    #####
    # Various methods
    #####

    def __eq__(self,other):
        r"""
        If self and other have the same transition function and the
        same word, then they are equal.

        TEST:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words([0,1,2])([0,1,0,1,2])
            sage: u = Words([0,1,2])(iter([0,1,0,1,2]))[:5]
            sage: ImplicitSuffixTree(w) == ImplicitSuffixTree(u)
            True
        """
        if not isinstance(other,ImplicitSuffixTree):
            return False
        return self._transition_function == other._transition_function \
            and self._word_content == other._word_content

    def transition_function(self,word,node=0):
        r"""
        Returns the node obtained by starting from node and following the
        edges labelled by the letters of word. Returns ("explicit",
        end_node) if we end at end_node, or ("implicit", (edge, d)) if we
        end d spots along an edge.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.transition_function(W([0,1,0]))
            ('implicit', (3, 1), 1)
            sage: t.transition_function(W([0,1,2]))
            ('explicit', 4)
            sage: t.transition_function(W([0,1,2]), 5)
            ('explicit', 2)
            sage: t.transition_function(W([0,1]), 5)
            ('implicit', (5, 2), 2)
        """
        if len(word) == 0:
            return "explicit", node
        ((k,p),s) = self._find_transition(node, word._word_content[0])
        if p is None:
            # test that word is a prefix of self._word_content[k-1:]
            if word == self._word[k-1:(k-1)+len(word)]:
                if len(word) == len(self._word_content) - k + 1:
                    return "explicit", s
                else:
                    edge = (node,s)
                    return "implicit", edge, len(word)
        else:
            # find longest common prefix
            m = min(p-k+1,len(word))
            i = 0
            while i < m and self._word[k-1+i] == word[i]:
                i += 1
            if i == p-k+1:
                return self.transition_function(word[p-k+1:],s)
            else:
                edge = (node,s)
                return "implicit", edge, i
            return "explicit", node

    def states(self):
        r"""
        Returns the states (explicit nodes) of the suffix tree.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.states()
            [0, 1, 2, 3, 4, 5, 6, 7]
        """
        return range(len(self._transition_function))

    def suffix_link(self, state):
        r"""
        Evaluates the suffix link map of the implicit suffix tree on state.
        Note that the suffix link is not defined for all states.

        The suffix link of a state x' that corresponds to the suffix x is
        defined to be -1 is x' is the root (0) and y' otherwise, where y'
        is the state corresponding to the suffix x[1:].

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.suffix_link(3)
            5
            sage: t.suffix_link(5)
            0
            sage: t.suffix_link(0)
            -1
            sage: t.suffix_link(-1)
            Traceback (most recent call last):
                ...
            TypeError: there is no suffix link from -1
        """
        if self._suffix_link.has_key(state):
            return self._suffix_link[state]
        else:
            raise TypeError, "there is no suffix link from %s" % state

    def active_state(self):
        r"""
        Returns the active state of the suffix tree.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.active_state()
            (0, (6, 6))
        """
        return self._active_state

    def process_letter(self, letter):
        r"""
        Modifies the current implicit suffix tree producing the implicit
        suffix tree for self.word() + letter.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("cacao")
            sage: t = ImplicitSuffixTree(w[:-1]); t
            Implicit Suffix Tree of the word: caca
            sage: t.process_letter(w[-1]); t
            Implicit Suffix Tree of the word: cacao

            sage: W = Words([0,1])
            sage: s = ImplicitSuffixTree(W([0,1,0,1])); s
            Implicit Suffix Tree of the word: 0101
            sage: s.process_letter(W([1])[0]); s
            Implicit Suffix Tree of the word: 01011
        """
        letter2int = self._word.alphabet().rank
        self._word_content.append(letter2int(letter))
        self._process_letter(self._word_content[-1])

    def to_explicit_suffix_tree(self):
        r"""
        Converts self to an explicit suffix tree. It is obtained by
        processing an end of string letter as if it were a regular
        letter, except that no new leaf nodes are created (thus, the only
        thing that happens is that some implicit nodes become explicit).

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("cacao")
            sage: t = ImplicitSuffixTree(w[:-1])
            sage: t.to_explicit_suffix_tree()

            sage: W = Words([0,1])
            sage: s = ImplicitSuffixTree(W([0,1,0,1, 1]))
            sage: s.to_explicit_suffix_tree()
        """
        # test whether end_of_string appears in the word
        alphabet = self._word.alphabet()
        if end_of_string in alphabet:
            raise TypeError, "word contains the end of string symbol"
        # append the end of string symbol to the word and process the new letter
        self._word_content.append(end_of_string)
        (s,(k,i)) = self._active_state
        old_r = 0
        (end_state, r) = self._test_and_split(s,(k,i-1), end_of_string)
        while end_state == False:
            (s, k) = self._canonize(self._suffix_link[s], (k,i-1))
            (end_state, r) = self._test_and_split(s, (k,i-1), end_of_string)
        self._word_content = self._word_content[:-1]
        return

    def edge_iterator(self):
        r"""
        Returns an iterator over the edges of the suffix tree. The
        edge from u to v labelled by (i,j) is returned as the tuple
        (u,v,(i,j)).

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: sorted( ImplicitSuffixTree(Word("aaaaa")).edge_iterator() )
            [(0, 1, (0, None))]
            sage: sorted( ImplicitSuffixTree(Word([0,1,0,1])).edge_iterator() )
            [(0, 1, (0, None)), (0, 2, (1, None))]
            sage: sorted( ImplicitSuffixTree(Word()).edge_iterator() )
            []
        """
        queue = [0]
        while queue:
            v=queue.pop()
            for ((i,j),u) in self._transition_function[v].iteritems():
                yield (v,u,(i-1,j))
                queue.append(u)

    def number_of_factors(self,n=None):
        r"""
        Count the number of distinct factors of self.word().

        INPUT:
            n -- an integer, or None.

        OUTPUT:
            If n is an integer, returns the number of distinct factors
            of length n. If n is None, returns the total number of
            distinct factors.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: t = ImplicitSuffixTree(Word([1,2,1,3,1,2,1]))
            sage: t.number_of_factors()
            22
            sage: t.number_of_factors(1)
            3
            sage: t.number_of_factors(9)
            0
            sage: t.number_of_factors(0)
            1

            sage: t = ImplicitSuffixTree(Word("cacao"))
            sage: t.number_of_factors()
            13
            sage: map(t.number_of_factors, range(10))
            [1, 3, 3, 3, 2, 1, 0, 0, 0, 0]

            sage: t = ImplicitSuffixTree(Word("c"*1000))
            sage: t.number_of_factors()
            1001
            sage: t.number_of_factors(17)
            1
            sage: t.number_of_factors(0)
            1

            sage: ImplicitSuffixTree(Word()).number_of_factors()
            1

            sage: blueberry = ImplicitSuffixTree(Word("blueberry"))
            sage: blueberry.number_of_factors()
            43
            sage: map(blueberry.number_of_factors, range(10))
            [1, 6, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        if n is None:
            length_word = len(self.word())
            num_factors = 1 # empty word
            for (u,v,(i,j)) in self.edge_iterator():
                if j == None:
                    num_factors += length_word - i
                else:
                    num_factors += j - i
        elif isint(n):
            length_word = len(self.word())
            num_factors = 0
            queue = [(0, 0)]
            while queue:
                (v,l) = queue.pop()
                if l == n:
                    num_factors += 1
                if l < n:
                    if self._transition_function[v] != {}:
                        for ((i,j),u) in self._transition_function[v].iteritems():
                            if j == None:
                                j = len(self.word())
                            if j - i >= n - l:
                                num_factors += 1
                            else:
                                queue.append((u,l+j-i+1))
        else:
            raise TypeError, "not an integer or None: %s" %s
        return num_factors

    def factor_iterator(self,n=None):
        r"""
        Generate distinct factors of self.

        INPUT:
            n -- an integer, or None.

        OUTPUT:
            If n is an integer, returns an iterator over all distinct
            factors of length n. If n is None, returns an iterator
            generating all distinct factors.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: sorted( ImplicitSuffixTree(Word("cacao")).factor_iterator() )
            [word: , word: a, word: ac, word: aca, word: acao, word: ao, word: c, word: ca, word: cac, word: caca, word: cacao, word: cao, word: o]
            sage: sorted( ImplicitSuffixTree(Word("cacao")).factor_iterator(1) )
            [word: a, word: c, word: o]
            sage: sorted( ImplicitSuffixTree(Word("cacao")).factor_iterator(2) )
            [word: ac, word: ao, word: ca]
            sage: sorted( ImplicitSuffixTree(Word([0,0,0])).factor_iterator() )
            [word: , word: 0, word: 00, word: 000]
            sage: sorted( ImplicitSuffixTree(Word([0,0,0])).factor_iterator(2) )
            [word: 00]
            sage: sorted( ImplicitSuffixTree(Word([0,0,0])).factor_iterator(0) )
            [word: ]
            sage: sorted( ImplicitSuffixTree(Word()).factor_iterator() )
            [word: ]
            sage: sorted( ImplicitSuffixTree(Word()).factor_iterator(2) )
            []
        """
        # Every factor is a prefix of a suffix, so we do a depth
        # first search of the implicit suffix tree of the word.
        if n is None:
            queue = [(0, self._word.parent()())]
            while queue:
                (v,w) = queue.pop()
                yield w
                if self._transition_function[v] != {}:
                    for ((i,j),u) in self._transition_function[v].iteritems():
                        if j == None:
                            j = len(self.word())
                        for k in range(i,j):
                            yield w * self.word()[i-1:k]
                        queue.append((u,w*self.word()[i-1:j]))
        elif isint(n):
            queue = [(0, self._word.parent()())]
            while queue:
                (v,w) = queue.pop()
                length_w = len(w)
                if length_w == n:
                    yield w
                if length_w < n:
                    if self._transition_function[v] != {}:
                        for ((i,j),u) in self._transition_function[v].iteritems():
                            if j == None:
                                j = len(self.word())
                            if j - i >= n - length_w:
                                yield w*self.word()[i-1:i-1+n-length_w]
                            else:
                                queue.append((u,w*self.word()[i-1:j]))
        else:
            raise TypeError, "not an integer or None: %s" %s

    #####
    # Miscellaneous methods
    #####

    def uncompactify(self):
        r"""
        Returns the tree obtained from self by splitting edges so that they
        are labelled by exactly one letter. The resulting tree is
        isomorphic to the suffix trie.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree, SuffixTrie
            sage: abbab = Words("ab")("abbab")
            sage: s = SuffixTrie(abbab)
            sage: t = ImplicitSuffixTree(abbab)
            sage: t.uncompactify().is_isomorphic(s.to_digraph())
            True
        """
        tree = self.to_digraph(word_labels=True)
        newtree = DiGraph()
        newtree.add_vertices(range(tree.order()))
        new_node = tree.order() + 1
        for (u,v,label) in tree.edge_iterator():
            if len(label) == 1:
                newtree.add_edge(u,v)
            else:
                newtree.add_edge(u,new_node,label[0]);
                for w in label[1:-1]:
                    newtree.add_edge(new_node,new_node+1,w)
                    new_node += 1
                newtree.add_edge(new_node,v,label[-1])
                new_node += 1
        return newtree

    def trie_type_dict(self):
        r"""
        Returns a dictionary in a format compatible with that of the suffix
        trie transition function.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree, SuffixTrie
            sage: W = Words("ab")
            sage: t = ImplicitSuffixTree(W("aba"))
            sage: t.trie_type_dict() == dict([[(0, W("a")), 4], [(0, W("b")), 3], [(3, W("a")), 2], [(4, W("b")), 5], [(5, W("a")), 1]])
            True
        """
        d = {}
        new_node = len(self._transition_function)
        for (u, dd) in self._transition_function.iteritems():
            for (sl, v) in dd.iteritems():
                w = self._word[sl[0]-1:sl[1]]
                if len(w) == 1:
                    d[u,w] = v
                else:
                    d[u,w[0:1]] = new_node
                    for i in range(1,len(w)-1):
                        d[new_node, w[i:i+1]] = new_node + 1
                        new_node += 1
                    d[new_node,w[-1:]] = v
                    new_node += 1
        return d

################################################################################
# For testing
################################################################################

def NaiveSuffixTree(word):
    r"""
    Given a word, constructs the suffix tree of the word using a naive
    algorithm. Useful for testing.

    INPUT:
        word -- any word

    EXAMPLES:
        sage: from sage.combinat.words.suffix_trees import NaiveSuffixTree
        sage: W = Words('abco')
        sage: s = NaiveSuffixTree(W("abacabaabaca")); s
        Suffix Tree of the word: abacabaabaca

        sage: s = NaiveSuffixTree(W("cacao")); s
        Suffix Tree of the word: cacao
        sage: s.edges()
        [(0, 3, word: ca), (0, 5, word: a), (0, 7, word: o), (3, 1, word: cao), (3, 4, word: o), (5, 2, word: cao), (5, 6, word: o)]

        sage: W = Words([0,1])
        sage: s = NaiveSuffixTree(W([0,1,0,1,1])); s
        Suffix Tree of the word: 01011
        sage: s.edges()
        [(0, 3, word: 01), (0, 5, word: 1), (3, 1, word: 011), (3, 4, word: 1), (5, 2, word: 011), (5, 6, word: 1)]
    """
    # The algorithm begins from the suffix tree of the empty word (one node).
    ST = NaiveSuffixTreeClass({0:[]})
    ST.set_vertex(0,{'position':0,'suffix':True})

    # test whether the end of string symbol appears in the word
    alphabet = word.alphabet()
    if end_of_string in alphabet:
        raise TypeError, "word contains the end of string symbol"
    # append end_of_string to the word and process the new letter
    W = Words(chain(alphabet, [end_of_string]))
    word *= W([end_of_string])

    # run the algorithm
    for i in range(len(word)):
        ST.naive_process_suffix(word[i:])

    # prune the tree (remove end_of_string).
    for (u,v,l) in ST.edges():
        if l[-1] == end_of_string:
            if l == W([end_of_string]):
                ST.delete_edge(u,v)
                ST.delete_vertex(v)
            else:
                ST.set_edge_label(u,v,l[:-1])
    return ST

class NaiveSuffixTreeClass(DiGraph):

    def word_to_node(self,node,word):
        r"""
        Returns the node obtained by starting from node and following
        the edges labelled by the letters of word. Returns ("node",
        node2, word2) if we end at node2, or ("implicit", edge, word2)
        if we end in the middle of an edge (called an implicit node),
        where word2 is the suffix of word that does not match any more
        edge labels.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import NaiveSuffixTree
            sage: W = Words("abc")
            sage: NST = NaiveSuffixTree(W("abacabaabaca"))
            sage: NST.word_to_node(0, W("ca"))
            ('node', 17, word: )
            sage: NST.word_to_node(0, W("cab"))
            ('implicit', ((17, 5, word: baabaca), 1), word: )
            sage: NST.word_to_node(17, W("b"))
            ('implicit', ((17, 5, word: baabaca), 1), word: )
            sage: NST.word_to_node(15, W("baabacabca"))
            ('node', 4, word: bca)
        """
        if len(word) == 0:
            return "node", node, word
        for v in self.successor_iterator(node):
            label = self.edge_label(node,v)
            if label[0] != word[0]:
                continue
            # find longest common prefix
            m = min(len(label),len(word))
            i = 0
            while i < m and label[i] == word[i]:
                i += 1
            lcp = label[:i]
            if len(lcp) == len(label):
                return self.word_to_node(v,word[len(label):])
            else:
                edge = ((node,v,label),i)
                return "implicit", edge, word[len(lcp):]
        return "node", node, word

    def naive_process_suffix(self,suffix):
        r"""
        Helper function for constructing the suffix tree.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import NaiveSuffixTreeClass
            sage: W = Words([0,1,2])
            sage: NST = NaiveSuffixTreeClass({0:[]})
            sage: NST.set_vertex(0,{'position':0,'suffix':True})
            sage: NST.naive_process_suffix(W([0,1,1,2,0])); NST
            Suffix Tree of the word: 01120
        """
        if self.order() == 1:
            self.add_vertex(1)
            self.add_edge(0,1)
            self.set_edge_label(0,1,suffix)
        else:
            typ, data, unmatched_suffix = self.word_to_node(0,suffix)
            if unmatched_suffix == "":
                return
            if typ == "node":
                node = data
            elif typ == "implicit":
                ((u,v,label),i) = data
                w = self.order()
                self.add_vertex(w)
                self.add_edge(u,w)
                self.set_edge_label(u,w,label[:i])
                self.add_edge(w,v)
                self.set_edge_label(w,v,label[i:])
                self.delete_edge(u,v)
                node = w
            w = self.order()
            self.add_vertex(w)
            self.add_edge(node,w)
            self.set_edge_label(node,w,unmatched_suffix)
        return

    def plot(self, layout='tree', tree_root=0, tree_orientation='up',
            vertex_colors=None, edge_labels=True, *args, **kwds):
        r"""
        Returns a Graphics object associated to self.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import NaiveSuffixTree
            sage: NaiveSuffixTree(Word("abacabaabaca")).plot()

        TESTS:
            sage: from sage.combinat.words.suffix_trees import NaiveSuffixTree
            sage: type(NaiveSuffixTree(Word("abacabaabaca")).plot())
            <class 'sage.plot.plot.Graphics'>
        """
        if vertex_colors is None:
            veretex_colors = {'#fec7b8':self.vertices()}
        return super(NaiveSuffixTreeClass,self).plot(layout=layout,
                tree_root=tree_root, tree_orientation=tree_orientation,
                vertex_colors=vertex_colors, edge_labels=edge_labels,
                *args, **kwds)

    def node_to_word(self,node):
        r"""
        Returns the word obtained by reading the edge labels from the root to
        node.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import NaiveSuffixTree
            sage: NaiveSuffixTree(Word("abacabaabaca")).node_to_word(17)
            word: ca
            sage: NaiveSuffixTree(Word("cacao")).node_to_word(0)
            word:
            sage: NaiveSuffixTree(Word([0,1,1,0,1,0,0,1])).node_to_word(10)
            word: 001
        """
        s = Words("")()
        while node != 0:
            parent = self.predecessors(node)[0]
            s = self.edge_label(parent,node) * s
            node = parent
        return s

    def word(self):
        r"""
        Returns the word whose suffix tree this is.

        EXAMPLES:
            sage: from sage.combinat.words.suffix_trees import NaiveSuffixTree
            sage: NaiveSuffixTree(Word("abcba")).word()
            word: abcba
        """
        return self.node_to_word(1)

    def _repr_(self):
        """
        TESTS:
            sage: from sage.combinat.words.suffix_trees import NaiveSuffixTree
            sage: NaiveSuffixTree(Word("abcba"))._repr_()
            'Suffix Tree of the word: abcba'
        """
        return 'Suffix Tree of the %s' % self.word()

###########################################################################
# The following would be very useful to have, but we need to first settle
# on a definition of suffix tree---there are slight differences that would
# cause the following to fail.
###########################################################################
#
#    def is_valid(self):
#        r"""
#        This verifies that the construction is correct.
#        For each suffix we compute the index of the
#        first occurrence of the suffix in self.word()
#        and compare this with the position label of
#        the corresponding node.
#        """
#        w = self.word()
#        W = Words(set(w))(w)
#        # Tests that the node labels are set correctly for each suffix:
#        #   position is the first occurrence of the suffix in the word;
#        #   suffix is True
#        for i in range(len(w)):
#            typ, data, unmatched_suffix = self.word_to_node(0,w[i:])
#            if typ == "implicit":
#                print "implicit node for suffix in position %s" % i
#            if type == "node":
#                v = data
#                d = self.get_vertex(v)
#                if d['position'] != W.first_pos_in(W[i:]) or not d['suffix']:
#                    raise ValueError, "position of node corresponding to suffix
#%s is wrong" % i
#
#        # Tests that the node labels are set correctly (for each node):
#        # let u = node_to_word(node); then
#        #   position is the first occurrence u in word;
#        #   suffix is True if u is a suffix, False otherwise.
#        for v in self.vertices():
#            d = self.get_vertex(v)
#            u = self.node_to_word(v)
#            U = Words(set(u))(u)
#            if d['position'] != U.first_pos_in(W):
#                raise ValueError, "position of node %s is wrong" % v
#            if d['suffix'] != U.is_suffix_of(W):
#                raise ValueError, "suffix bit of node %s is wrong" % v
#        return True
