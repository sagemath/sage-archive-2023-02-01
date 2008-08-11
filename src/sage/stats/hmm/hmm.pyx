r"""
Hidden Markov Models

AUTHOR: William Stein
"""

#############################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version at your option.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

import math

from sage.matrix.all import is_Matrix, matrix
from sage.rings.all  import RDF
from sage.misc.randstate import random

include "../../ext/stdsage.pxi"

include "misc.pxi"

cdef class HiddenMarkovModel:
    """
    Abstract base class for all Hidden Markov Models.

    INPUT:
        A -- matrix or list
        B -- matrix or list
        pi -- list of floats

    EXAMPLES:
    One shouldn't directly call this constructor since this class is an abstract
    base class.
        sage: sage.stats.hmm.hmm.HiddenMarkovModel([[0.4,0.6],[1,0]], [[1,0],[0.5,0.5]], [0.5,0.5])
        <sage.stats.hmm.hmm.HiddenMarkovModel object at ...>
    """
    def __init__(self, A, B, pi=None):
        """
        INPUT:
            A -- matrix or list
            B -- matrix or list
            pi -- list of floats

        EXAMPLES:
            sage: hmm.DiscreteHiddenMarkovModel([[1]], [[0.0,1.0]], [1])
            Discrete Hidden Markov Model with 1 States and 2 Emissions
            Transition matrix:
            [1.0]
            Emission matrix:
            [0.0 1.0]
            Initial probabilities: [1.0]
        """
        cdef Py_ssize_t n

        # Convert A to a matrix
        if not is_Matrix(A):
            n = len(A)
            A = matrix(RDF, n, len(A[0]) if n > 0 else 0, A)

        # Convert B to a matrix
        if not is_Matrix(B):
            n = len(B)
            B = matrix(RDF, n, len(B[0]) if n > 0 else 0, B)

        # Some consistency checks
        if not A.is_square():
            print A.parent()
            raise ValueError, "A must be square"
        if A.nrows() != B.nrows():
            raise ValueError, "number of rows of A and B must be the same"

        # Make sure A and B are over RDF.
        if A.base_ring() != RDF:
            A = A.change_ring(RDF)
        if B.base_ring() != RDF:
            B = B.change_ring(RDF)

        # Make sure the initial probabilities are all floats.
        if pi is None:
            if A.nrows() == 0:
                self.pi = []
            else:
                self.pi = [1.0/A.nrows()]*A.nrows()
        else:
            self.pi = [float(x) for x in pi]
            if len(self.pi) != A.nrows():
                raise ValueError, "length of pi must equal number of rows of A"

        # Record the now validated matrices A and B as attributes.
        # They get used later as attributes in the constructors for
        # derived classes.
        self.A = A
        self.B = B

        # Check that all entries of A are nonnegative and raise a
        # ValueError otherwise. Note that we don't check B since it
        # has entries that are potentially negative in the continuous
        # case.  But GHMM prints clear warnings when the emission
        # probabilities are negative, i.e., it does not silently give
        # wrong results like it does for negative transition
        # probabilities.
        cdef Py_ssize_t i, j
        for i from 0 <= i < self.A._nrows:
            for j from 0 <= j < self.A._ncols:
                if self.A.get_unsafe_double(i,j) < 0:
                    raise ValueError, "each transition probability must be nonnegative"



cdef class DiscreteHiddenMarkovModel(HiddenMarkovModel):
    """
    Create a discrete hidden Markov model.

    hmm.DiscreteHiddenMarkovModel(A, B, pi=None, emission_symbols=None, name=None, normalize=True)
n
    INPUTS:
        A  -- square matrix of doubles; the state change probabilities
        B  -- matrix of doubles; emission probabilities
        pi -- list of floats; probabilities for each initial state
        emission_state -- list of B.ncols() symbols (just used for printing)
        name -- (optional) name of the model
        normalize -- (optional; default=True) whether or not to normalize
                     the model so the probabilities add to 1

    EXAMPLES:
    We create a discrete HMM with 2 internal states on an alphabet of size 2.
        sage: hmm.DiscreteHiddenMarkovModel([[0.2,0.8],[0.5,0.5]], [[1,0],[0,1]], [0,1])
        Discrete Hidden Markov Model with 2 States and 2 Emissions
        Transition matrix:
        [0.2 0.8]
        [0.5 0.5]
        Emission matrix:
        [1.0 0.0]
        [0.0 1.0]
        Initial probabilities: [0.0, 1.0]

    The transition probabilities must be nonnegative:
        sage: hmm.DiscreteHiddenMarkovModel([[-0.2,0.8],[0.5,0.5]], [[1,0],[0,1]], [0,1])
        Traceback (most recent call last):
        ...
        ValueError: each transition probability must be nonnegative

    The transition probabilities are by default automatically normalized:
        sage: a = hmm.DiscreteHiddenMarkovModel([[0.2,0.3],[0.5,0.5]], [[1,0],[0,1]], [0,1])
        sage: a.transition_matrix()
        [0.4 0.6]
        [0.5 0.5]
    """
    def __init__(self, A, B, pi=None, emission_symbols=None, name=None, normalize=True):
        """
        EXAMPLES:
            sage: hmm.DiscreteHiddenMarkovModel([[1]], [[0.0,1.0]], [1])
            Discrete Hidden Markov Model with 1 States and 2 Emissions
            Transition matrix:
            [1.0]
            Emission matrix:
            [0.0 1.0]
            Initial probabilities: [1.0]
        """
        # memory has not all been setup yet.
        self.initialized = False

        # This sets self.A, self.B and pi after doing appropriate coercions, etc.
        HiddenMarkovModel.__init__(self, A, B, pi)

        self.set_emission_symbols(emission_symbols)

        self.m = <ghmm_dmodel*> safe_malloc(sizeof(ghmm_dmodel))

        self.m.label = to_int_array(range(len(self._emission_symbols)))

        # Set all pointers to NULL
        self.m.s = NULL; self.m.name = NULL; self.m.silent = NULL
        self.m.tied_to = NULL; self.m.order = NULL; self.m.background_id = NULL
        self.m.bp = NULL; self.m.topo_order = NULL; self.m.pow_lookup = NULL;
        self.m.label_alphabet = NULL; self.m.alphabet = NULL

        # Set number of states and number of outputs
        self.m.N = self.A.nrows()
        self.m.M = len(self._emission_symbols)
        # Set the model type to discrete
        self.m.model_type = GHMM_kDiscreteHMM

        # Set that no a prior model probabilities are set.
        self.m.prior = -1

        # Assign model identifier if specified
        if name is not None:
            name = str(name)
            self.m.name = <char*> safe_malloc(len(name))
            strcpy(self.m.name, name)
        else:
            self.m.name = NULL

        # Allocate and initialize states
        cdef ghmm_dstate* states = <ghmm_dstate*> safe_malloc(sizeof(ghmm_dstate) * self.m.N)
        cdef ghmm_dstate* state

        silent_states = []
        tmp_order     = []

        cdef Py_ssize_t i, j, k

        for i from 0 <= i < self.m.N:
            v = self.B[i]

            # Get a reference to the i-th state for convenience of the notation below.
            state = &(states[i])
            state.desc = NULL

            # Compute state order
            if self.m.M > 1:
                order = math.log( len(v), self.m.M ) - 1
            else:
                order = len(v) - 1

            # Check for valid number of emission parameters
            order = int(order)
            if self.m.M**(order+1) == len(v):
                tmp_order.append(order)
            else:
                raise ValueError, "number of columns (= %s) of B must be a power of the number of emission symbols (= %s)"%(
                    self.B.ncols(), len(emission_symbols))

            state.b = to_double_array(v)
            state.pi = self.pi[i]

            silent_states.append( 1 if sum(v) == 0 else 0)

            # Set "out" probabilities, i.e., the probabilities to
            # transition to another hidden state from this state.
            v = self.A[i]
            k = self.m.N
            state.out_states = k
            state.out_id = <int*> safe_malloc(sizeof(int)*k)
            state.out_a  = <double*> safe_malloc(sizeof(double)*k)
            for j from 0 <= j < k:
                state.out_id[j] = j
                state.out_a[j]  = v[j]

            # Set "in" probabilities
            v = self.A.column(i)
            state.in_states = k
            state.in_id = <int*> safe_malloc(sizeof(int)*k)
            state.in_a  = <double*> safe_malloc(sizeof(double)*k)
            for j from 0 <= j < k:
                state.in_id[j] = j
                state.in_a[j]  = v[j]

            state.fix = 0

        self.m.s = states
        self.initialized = True
        if normalize:
            self.normalize()

    def __cmp__(self, other):
        """
        Compare two Discrete HMM's.

        INPUT:
            self, other -- objects; if other is not a Discrete HMM compare types.
        OUTPUT:
            -1,0,1

        The transition matrices are compared, then the emission
        parameters, and the initial probabilities.  The names are not
        compared, so two GHMM's with the same parameters but different
        names compare to be equal.

        EXAMPLES:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[1,1]], [1,2], normalize=False)
            sage: m.__cmp__(m)
            0

        Note that the name doesn't matter:
            sage: n = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[1,1]], [1,2], normalize=False)
            sage: m.__cmp__(n)
            0

        Normalizing fixes the initial probabilities, hence m and n are no longer equal.
            sage: n.normalize()
            sage: m.__cmp__(n)
            1
        """
        if not isinstance(other, DiscreteHiddenMarkovModel):
            return cmp(type(self), type(other))

        if self is other: return 0  # easy special case

        cdef DiscreteHiddenMarkovModel o = other
        if self.m.N < o.m.N:
            return -1
        elif self.m.N > o.m.N:
            return 1
        cdef Py_ssize_t i, j

        # This code is similar to code in chmm.pyx, but with several small differences.

        # Compare transition matrices
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.m.s[i].out_states:
                if self.m.s[i].out_a[j] < o.m.s[i].out_a[j]:
                    return -1
                elif self.m.s[i].out_a[j] > o.m.s[i].out_a[j]:
                    return 1

        # Compare emission matrices
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.B._ncols:
                if self.m.s[i].b[j] < o.m.s[i].b[j]:
                    return -1
                elif self.m.s[i].b[j] > o.m.s[i].b[j]:
                    return 1

        # Compare initial state probabilities
        for 0 <= i < self.m.N:
            if self.m.s[i].pi < o.m.s[i].pi:
                return -1
            elif self.m.s[i].pi > o.m.s[i].pi:
                return 1

        # Compare emission symbols
        return cmp(self._emission_symbols, o._emission_symbols)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[1,1]], [0,1], ['a','b'], name='test model')
            sage: f,g = m.__reduce__()
            sage: f(*g) == m
            True
        """
        return unpickle_discrete_hmm_v0, (self.transition_matrix(), self.emission_matrix(),
                      self.initial_probabilities(), self._emission_symbols, self.name())

    def __dealloc__(self):
        """
        Deallocate memory allocated by the HMM.

        EXAMPLES:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.2,0.8],[0.5,0.5]], [[1,0],[0,1]], [0,1])  # indirect doctest
            sage: del m
        """
        if self.initialized:
            ghmm_dmodel_free(&self.m)

    def fix_emissions(self, Py_ssize_t i, bint fixed=True):
        """
        Sets the i-th emission parameters to be either fixed or not
        fixed.  If it is fixed, then running the Baum-Welch algorithm
        will not change the emission parameters for the i-th state.

        INPUT:
            i -- nonnegative integer < self.m.N
            fixed -- bool

        EXAMPLES:
        First without calling fix_emissions:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],[[.5,.5],[.5,.5]])
            sage: m.baum_welch([0,0,0,1,1,1])
            sage: m.emission_matrix()
            [              1.0               0.0]
            [3.92881039079e-05    0.999960711896]

        We call fix_emissions on the first state and notice that the first
        row of the emission matrix does not change:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],[[.5,.5],[.5,.5]])
            sage: m.fix_emissions(0)
            sage: m.baum_welch([0,0,0,1,1,1])
            sage: m.emission_matrix()
            [              0.5               0.5]
            [0.000542712675606    0.999457287324]

        We call fix_emissions on the second state and notice that the second
        row of the emission matrix does not change:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],[[.5,.5],[.5,.5]])
            sage: m.fix_emissions(1)
            sage: m.baum_welch([0,0,0,1,1,1])
            sage: m.emission_matrix()
            [   0.999999904763 9.52366620142e-08]
            [              0.5               0.5]

        TESTS:
        Make sure that out of range indices are handled correctly with an IndexError.
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]],[[.5,.5],[.5,.5]])
            sage: m.fix_emissions(2)
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: m.fix_emissions(-1)
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i < 0 or i >= self.m.N:
            raise IndexError, "index out of range"
        self.m.s[i].fix = fixed

    def __repr__(self):
        """
        Return string representation of this HMM.

        OUTPUT:
            string

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[0.9,0.1],[0.1,0.9]], [0.5,0.5], [3/4, 'abc'])
            sage: a.__repr__()
            "Discrete Hidden Markov Model with 2 States and 2 Emissions\nTransition matrix:\n[0.1 0.9]\n[0.1 0.9]\nEmission matrix:\n[0.9 0.1]\n[0.1 0.9]\nInitial probabilities: [0.5, 0.5]\nEmission symbols: [3/4, 'abc']"
        """
        s = "Discrete Hidden Markov Model%s with %s States and %s Emissions"%(
            ' ' + self.m.name if self.m.name else '',
            self.m.N, self.m.M)
        s += '\nTransition matrix:\n%s'%self.transition_matrix()
        s += '\nEmission matrix:\n%s'%self.emission_matrix()
        s += '\nInitial probabilities: %s'%self.initial_probabilities()
        if self._emission_symbols_dict:
            s += '\nEmission symbols: %s'%self._emission_symbols
        return s

    def name(self):
        """
        Return the name of this model.

        OUTPUT:
            string or None

        EXAMPLES:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[1,1]], [1,2], name='test model')
            sage: m.name()
            'test model'

        If the model is not explicitly named then this function returns None:
            sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[1,1]], [1,2])
            sage: m.name() is None
            True
        """
        if self.m.name:
            s = str(self.m.name)
            return s
        else:
            return None

    def initial_probabilities(self):
        """
        Return the list of initial state probabilities.

        OUTPUT:
            list of floats

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.9,0.1],[0.9,0.1]], [[0.5,0.5,0,0],[0,0,.5,.5]], [0.5,0.5], [1,-1,3,-3])
            sage: a.initial_probabilities()
            [0.5, 0.5]
        """
        cdef Py_ssize_t i
        return [self.m.s[i].pi for i in range(self.m.N)]

    def transition_matrix(self, list_only=True):
        """
        Return the hidden state transition matrix.

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.9,0.1],[0.9,0.1]], [[0.5,0.5,0,0],[0,0,.5,.5]], [0.5,0.5], [1,-1,3,-3])
            sage: a.transition_matrix()
            [0.9 0.1]
            [0.9 0.1]
        """
        cdef Py_ssize_t i, j
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.m.s[i].out_states:
                self.A.set_unsafe_double(i,j,self.m.s[i].out_a[j])
        return self.A

    def emission_matrix(self, list_only=True):
        """
        Return the emission probability matrix.

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.9,0.1],[0.9,0.1]], [[0.5,0.5,0,0],[0,0,.5,.5]], [0.5,0.5], [1,-1,3,-3])
            sage: a.emission_matrix()
            [0.5 0.5 0.0 0.0]
            [0.0 0.0 0.5 0.5]
        """
        cdef Py_ssize_t i, j
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.B._ncols:
                self.B.set_unsafe_double(i,j,self.m.s[i].b[j])
        return self.B

    def normalize(self):
        """
        Normalize the transition and emission probabilities, if applicable.

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,1],[1.2,0.9]], [[1,0.5],[0.5,1]], [0.1,1.2])
            sage: a.normalize()
            sage: a
            Discrete Hidden Markov Model with 2 States and 2 Emissions
            Transition matrix:
            [0.333333333333 0.666666666667]
            [0.571428571429 0.428571428571]
            Emission matrix:
            [0.666666666667 0.333333333333]
            [0.333333333333 0.666666666667]
            Initial probabilities: [0.076923076923076927, 0.92307692307692302]
        """
        ghmm_dmodel_normalize(self.m)

    def sample(self, long length, number=None):
        """
        Return number samples from this HMM of given length.

        INPUT:
            length -- positive integer
            number -- (default: None) if given, compute list of this many sample sequences

        OUTPUT:
            if number is not given, return a single TimeSeries.
            if number is given, return a list of TimeSeries.

        EXAMPLES:
            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[1,0],[0,1]], [0,1])
            sage: print a.sample(10, 3)
            [[1, 0, 1, 1, 0, 1, 1, 0, 1, 0], [1, 1, 1, 1, 1, 1, 1, 1, 1, 0], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
            sage: a.sample(15)
            [1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1]
            sage: list(a.sample(1000)).count(0)
            95

        If the emission symbols are set
            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.1,0.9]], [[1,0],[0,1]], [0,1], ['up', 'down'])
            sage: a.sample(10)
            ['down', 'up', 'down', 'down', 'up', 'down', 'down', 'up', 'down', 'up']
        """
        seed = random()
        if number is None:
            number = 1
            single = True
        else:
            single = False
        cdef ghmm_dseq *d = ghmm_dmodel_generate_sequences(self.m, seed, length, number, length)
        cdef Py_ssize_t i, j
        v = [[d.seq[j][i] for i in range(length)] for j in range(number)]
        if self._emission_symbols_dict:
            w = self._emission_symbols
            v = [[w[i] for i in z] for z in v]
        if single:
            return v[0]
        return v

    def emission_symbols(self):
        """
        Return a copy of the list of emission symbols of self.

        Use set_emission_symbols to set the list of emission symbols.

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.1,0.9]], [[1,0],[0,1]], [0,1], ['up', -3/179])
            sage: a.emission_symbols()
            ['up', -3/179]
        """
        return list(self._emission_symbols)

    def set_emission_symbols(self, emission_symbols):
        """
        Set the list of emission symbols.

        EXAMPLES:
            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.1,0.9]], [[1,0],[0,1]], [0,1], ['up', 'down'])
            sage: a.set_emission_symbols([3,5])
            sage: a.emission_symbols()
            [3, 5]
            sage: a.sample(10)
            [5, 3, 5, 5, 3, 5, 5, 3, 5, 3]
            sage: a.set_emission_symbols([pi,5/9+e])
            sage: a.sample(10)
            [e + 5/9, e + 5/9, e + 5/9, e + 5/9, e + 5/9, e + 5/9, pi, pi, e + 5/9, pi]
        """
        if emission_symbols is None:
            self._emission_symbols = range(self.B.ncols())
            self._emission_symbols_dict = None
        else:
            self._emission_symbols = list(emission_symbols)
            if self._emission_symbols != range(self.B.ncols()):
                self._emission_symbols_dict = dict([(x,i) for i, x in enumerate(emission_symbols)])


    ####################################################################
    # HMM Problem 1 -- Computing likelihood: Given the parameter set
    # lambda of an HMM model and an observation sequence O, determine
    # the likelihood P(O | lambda).
    ####################################################################
    def log_likelihood(self, seq):
        r"""
        HMM Problem 1: Likelihood. Return $\log ( P[seq | model] )$,
        the log of the probability of seeing the given sequence given
        this model, using the forward algorithm and assuming
        independance of the sequence seq.

        INPUT:
            seq -- a list; sequence of observed emissions of the HMM

        OUTPUT:
            float -- the log of the probability of seeing this sequence
                     given the model

        WARNING: By convention we return -inf for 0 probability events.

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[1,0],[0,1]], [0,1])
            sage: a.log_likelihood([1,1])
            -0.10536051565782635
            sage: a.log_likelihood([1,0])
            -2.3025850929940459

        Notice that any sequence starting with 0 can't occur, since
        the above model always starts in a state that produces 1 with
        probability 1.  By convention log(probability 0) is -inf.
            sage: a.log_likelihood([0,0])
            -inf

        Here's a special case where each sequence is equally probable.
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.5,0.5]], [[1,0],[0,1]], [0.5,0.5])
            sage: a.log_likelihood([0,0])
            -1.3862943611198906
            sage: log(0.25)
            -1.38629436111989
            sage: a.log_likelihood([0,1])
            -1.3862943611198906
            sage: a.log_likelihood([1,0])
            -1.3862943611198906
            sage: a.log_likelihood([1,1])
            -1.3862943611198906
        """
        if self._emission_symbols_dict:
            seq = [self._emission_symbols_dict[z] for z in seq]
        cdef double log_p
        cdef int* O = to_int_array(seq)
        cdef int ret = ghmm_dmodel_logp(self.m, O, len(seq), &log_p)
        sage_free(O)
        if ret == -1:
            # forward returned -1: sequence can't be built
            return -float('Inf')
        return log_p

    ####################################################################
    # HMM Problem 2 -- Decoding: Given the complete parameter set that
    # defines the model and an observation sequence seq, determine the
    # best hidden sequence Q.  Use the Viterbi algorithm.
    ####################################################################
    def viterbi(self, seq):
        """
        HMM Problem 2: Decoding.  Determine a hidden sequence of
        states that is most likely to produce the given sequence seq
        of obserations.

        INPUT:
            seq -- sequence of emitted symbols

        OUTPUT:
            list -- a most probable sequence of hidden states, i.e., the
                    Viterbi path.
            float -- log of the probability that the sequence of hidden
                     states actually produced the given sequence seq.

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[0.9,0.1],[0.1,0.9]], [0.5,0.5])
            sage: a.viterbi([1,0,0,1,0,0,1,1])
            ([1, 0, 0, 1, 1, 0, 1, 1], -11.062453224772216)

        We predict the state sequence when the emissions are 3/4 and 'abc'.
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[0.9,0.1],[0.1,0.9]], [0.5,0.5], [3/4, 'abc'])

        Note that state 0 is common below, despite the model trying hard to
        switch to state 1:
            sage: a.viterbi([3/4, 'abc', 'abc'] + [3/4]*10)
            ([0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], -25.299405845367794)
        """
        if len(seq) == 0:
            return [], 0.0
        if self._emission_symbols_dict:
            seq = [self._emission_symbols_dict[z] for z in seq]
        cdef int* path
        cdef int* O = to_int_array(seq)
        cdef int pathlen
        cdef double log_p

        path = ghmm_dmodel_viterbi(self.m, O, len(seq), &pathlen, &log_p)
        sage_free(O)
        if not path:
            raise RuntimeError, "error computing viterbi path"
        p = [path[i] for i in range(pathlen)]
        sage_free(path)

        return p, log_p

    ####################################################################
    # HMM Problem 3 -- Learning: Given an observation sequence O and
    # the set of states in the HMM, improve the HMM to increase the
    # probability of observing O.
    ####################################################################
    def baum_welch(self, training_seqs, nsteps=None, log_likelihood_cutoff=None):
        """
        HMM Problem 3: Learning.  Given an observation sequence O and
        the set of states in the HMM, improve the HMM using the
        Baum-Welch algorithm to increase the probability of observing O.

        INPUT:
            training_seqs -- a list of lists of emission symbols (or a single list)
            nsteps -- integer or None (default: None) maximum number
                      of Baum-Welch steps to take
            log_likehood_cutoff -- positive float or None (default:
                      None); the minimal improvement in likelihood
                      with respect to the last iteration required to
                      continue. Relative value to log likelihood

        OUTPUT:
            changes the model in places, or raises a RuntimError
            exception on error

        EXAMPLES:
        We make a model that has two states and is equally likely to output
        0 or 1 in either state and transitions back and forth with equal
        probability.
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.5,0.5]], [[0.5,0.5],[0.5,0.5]], [0.5,0.5])

        We give the model some training data this much more likely to
        be 1 than 0.
            sage: a.baum_welch([[1,1,1,1,0,1], [1,0,1,1,1,1]])

        Now the model's emission matrix changes since it is much
        more likely to emit 1 than 0.
            sage: a
            Discrete Hidden Markov Model with 2 States and 2 Emissions
            Transition matrix:
            [0.5 0.5]
            [0.5 0.5]
            Emission matrix:
            [0.166666666667 0.833333333333]
            [0.166666666667 0.833333333333]
            Initial probabilities: [0.5, 0.5]

        Note that 1/6 = 1.666...:
            sage: 1.0/6
            0.166666666666667

        We run Baum-Welch on a single sequence:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.5,0.5]], [[0.5,0.5],[0.5,0.5]], [0.5,0.5])
            sage: a.baum_welch([1,0,1]*10)
            sage: a
            Discrete Hidden Markov Model with 2 States and 2 Emissions
            Transition matrix:
            [0.5 0.5]
            [0.5 0.5]
            Emission matrix:
            [0.333333333333 0.666666666667]
            [0.333333333333 0.666666666667]
            Initial probabilities: [0.5, 0.5]

        TESTS:
        We test training with non-default string symbols:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.5,0.5]], [[0.5,0.5],[0.5,0.5]], [0.5,0.5], ['up','down'])
            sage: a.baum_welch([['up','up'], ['down','up']])
            sage: a
            Discrete Hidden Markov Model with 2 States and 2 Emissions
            Transition matrix:
            [0.5 0.5]
            [0.5 0.5]
            Emission matrix:
            [0.75 0.25]
            [0.75 0.25]
            Initial probabilities: [0.5, 0.5]
            Emission symbols: ['up', 'down']

        NOTE: Training for models including silent states is not yet supported.

        REFERENCES:
            Rabiner, L.R.: "`A Tutorial on Hidden Markov Models and Selected
            Applications in Speech Recognition"', Proceedings of the IEEE,
            77, no 2, 1989, pp 257--285.
        """
        if len(training_seqs) > 0 and not isinstance(training_seqs[0], (list, tuple)):
            training_seqs = [training_seqs]

        if self._emission_symbols_dict:
            seqs = [[self._emission_symbols_dict[z] for z in x] for x in training_seqs]
        else:
            seqs = training_seqs

        cdef ghmm_dseq* d = malloc_ghmm_dseq(seqs)

        if ghmm_dmodel_baum_welch(self.m, d):
            raise RuntimeError, "error running Baum-Welch algorithm"

        ghmm_dseq_free(&d)


##################################################################################
# Helper Functions
##################################################################################

cdef ghmm_dseq* malloc_ghmm_dseq(seqs) except NULL:
    """
    Allocate a discrete sequence of samples.

    INPUT:
        seqs -- a list of sequences

    OUTPUT:
        C pointer to ghmm_dseq
    """
    cdef ghmm_dseq* d = ghmm_dseq_calloc(len(seqs))
    if d == NULL:
        raise MemoryError
    cdef int i, j, m, n
    m = len(seqs)
    d.seq_number = m
    d.capacity = m
    d.total_w = m
    for i from 0 <= i < m:
        v = seqs[i]
        n = len(v)
        d.seq[i] = <int*> safe_malloc(sizeof(int) * n)
        for j from 0 <= j < n:
            d.seq[i][j] = v[j]
        d.seq_len[i] = n
        d.seq_id[i] = i
        d.seq_w[i] = 1
    d.flags = 0
    return d


def unpickle_discrete_hmm_v0(A, B, pi, emission_symbols,name):
    """
    TESTS:
        sage: m = hmm.DiscreteHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [[0.0,1.0],[0.5,0.5]], [1,0], name='test model')
        sage: loads(dumps(m)) == m
        True
        sage: loads(dumps(m)).name()
        'test model'
        sage: sage.stats.hmm.hmm.unpickle_discrete_hmm_v0(m.transition_matrix(), m.emission_matrix(), m.initial_probabilities(), ['a','b'], m.name())
        Discrete Hidden Markov Model test model with 2 States and 2 Emissions
        Transition matrix:
        [0.4 0.6]
        [0.1 0.9]
        Emission matrix:
        [0.0 1.0]
        [0.5 0.5]
        Initial probabilities: [1.0, 0.0]
        Emission symbols: ['a', 'b']
    """
    return DiscreteHiddenMarkovModel(A,B,pi,emission_symbols,name)
