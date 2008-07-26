r"""
Hidden Markov Models

AUTHOR: William Stein

EXAMPLES:

"""

import math

from sage.matrix.all import is_Matrix, matrix
from sage.rings.all  import RDF
from sage.misc.randstate import random

from sage.matrix.matrix_real_double_dense cimport Matrix_real_double_dense

include "../../ext/stdsage.pxi"

cdef extern from "ghmm/ghmm.h":
    cdef int GHMM_kSilentStates
    cdef int GHMM_kHigherOrderEmissions
    cdef int GHMM_kDiscreteHMM

cdef extern from "ghmm/model.h":
    ctypedef struct ghmm_dstate:
        # Initial probability
        double pi
        # Output probability
        double *b

        # ID's of the following states
        int *out_id
        # ID's of the previous states
        int *in_id

        # Transition probabilities to successor states.
        double *out_a
        # Transition probabilities from predecessor states.
        double *in_a

        # Number of successor states
        int out_states
        # Number of precursor states
        int in_states
        # if fix == 1 then output probabilities b stay fixed during the training
        int fix
        # Contains a description of the state (null terminated utf-8)
        char * desc
        # x coordinate position for graph representation plotting
        int xPosition
        # y coordinate position for graph representation plotting
        int yPosition

    ctypedef struct ghmm_dbackground:
        pass

    ctypedef struct ghmm_alphabet:
        pass

    # Discrete HMM's.
    ctypedef struct ghmm_dmodel:
        # Number of states
        int N
        # Number of outputs
        int M
        # Vector of states
        ghmm_dstate *s
        # Contains bit flags for varios model extensions such as
        # kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
        int model_type
        # The a priori probability for the model.
        # A value of -1 indicates that no prior is defined.
        # Note: this is not to be confused with priors on emission distributions.
        double prior
        # An arbitrary name for the model (null terminated utf-8)
        char * name

        # Flag variables for each state indicating whether it is emitting or not.
        # Note: silent != NULL iff (model_type & kSilentStates) == 1
        int *silent

        # Int variable for the maximum level of higher order emissions.
        int maxorder

        # saves the history of emissions as int,
        # the nth-last emission is (emission_history * |alphabet|^n+1) % |alphabet|
        int emission_history

        # Flag variables for each state indicating whether the states emissions
        # are tied to another state. Groups of tied states are represented
        # by their tie group leader (the lowest numbered member of the group).
        # tied_to[s] == kUntied  : s is not a tied state
        # tied_to[s] == s        : s is a tie group leader
        # tied_to[t] == s        : t is tied to state s (t>s)
        # Note: tied_to != NULL iff (model_type & kTiedEmissions) != 0  */
        int *tied_to

        # Note: State store order information of the emissions.
        # Classical HMMS have emission order 0; that is, the emission probability
        # is conditioned only on the state emitting the symbol.
        # For higher order emissions, the emission are conditioned on the state s
        # as well as the previous emission_order observed symbols.
        # The emissions are stored in the state's usual double* b.
        # Note: order != NULL iff (model_type & kHigherOrderEmissions) != 0  */
        int * order

        # ghmm_dbackground is a pointer to a
        # ghmm_dbackground structure, which holds (essentially) an
        # array of background distributions (which are just vectors of floating
        # point numbers like state.b).
        # For each state the array background_id indicates which of the background
        # distributions to use in parameter estimation. A value of kNoBackgroundDistribution
        # indicates that none should be used.
        # Note: background_id != NULL iff (model_type & kHasBackgroundDistributions) != 0
        int *background_id
        ghmm_dbackground *bp

        # Topological ordering of silent states
        #  Condition: topo_order != NULL iff (model_type & kSilentStates) != 0
        int *topo_order
        int topo_order_length

        # pow_lookup is a array of precomputed powers
        # It contains in the i-th entry M (alphabet size) to the power of i
        # The last entry is maxorder+1.
        int *pow_lookup

        # Store for each state a class label. Limits the possibly state sequence
        # Note: label != NULL iff (model_type & kLabeledStates) != 0
        int* label
        ghmm_alphabet* label_alphabet
        ghmm_alphabet* alphabet


    # Discrete sequences
    ctypedef struct ghmm_dseq:
        # sequence array. sequence[i] [j] = j-th symbol of i-th seq
        int **seq
        # matrix of state ids, can be used to save the viterbi path during sequence generation.
        # ATTENTION: is NOT allocated by ghmm_dseq_calloc
        int **states
        # array of sequence length
        int *seq_len
        # array of state path lengths
        int *states_len
        # array of sequence IDs
        double *seq_id
        # positive sequence weights.  default is 1 = no weight
        double *seq_w
        # total number of sequences
        long seq_number
        # reserved space for sequences is always >= seq_number
        long capacity
        # sum of sequence weights
        double total_w
        # matrix of state labels corresponding to seq
        int **state_labels
        # number of labels for each sequence
        int *state_labels_len
        # internal flags
        unsigned int flags

    ghmm_dseq *ghmm_dmodel_label_generate_sequences(
        ghmm_dmodel * mo, int seed, int global_len, long seq_number, int Tmax)
    ghmm_dseq *ghmm_dmodel_generate_sequences(ghmm_dmodel* mo, int seed, int global_len,
                                          long seq_number, int Tmax)


    int ghmm_dmodel_normalize(ghmm_dmodel *m)
    int ghmm_dmodel_free(ghmm_dmodel **m)
    int ghmm_dmodel_logp(ghmm_dmodel *m, int *O, int len, double *log_p)
    int *ghmm_dmodel_viterbi (ghmm_dmodel *m, int *O, int len, int *pathlen, double *log_p)
    int ghmm_dmodel_baum_welch (ghmm_dmodel *m, ghmm_dseq *sq)
    int ghmm_dmodel_baum_welch_nstep (ghmm_dmodel * m, ghmm_dseq *sq, int max_step,
                                      double likelihood_delta)


cdef class HiddenMarkovModel:
    pass


cdef class DiscreteHiddenMarkovModel:
    cdef ghmm_dmodel* m
    cdef bint initialized
    cdef object _emission_symbols, _emission_symbols_dict
    cdef Matrix_real_double_dense A, B

    def __init__(self, A, B, pi, emission_symbols=None, name=None):
        """
        INPUTS:
            A  -- square matrix of doubles; the state change probabilities
            B  -- matrix of doubles; emission probabilities
            pi -- list of doubles; probabilities for each initial state
            emission_state -- list of B.ncols() symbols (just used for printing)
            name -- (optional) name of the model

        EXAMPLES:
        We create a discrete HMM with 2 internal states on an alphabet of
        size 2.

            sage: a = hmm.DiscreteHiddenMarkovModel([[0.2,0.8],[0.5,0.5]], [[1,0],[0,1]], [0,1])

        """
        if not is_Matrix(A):
            A = matrix(RDF, len(A), len(A[0]), A)
        if not is_Matrix(B):
            B = matrix(RDF, len(B), len(B[0]), B)
        if not A.is_square():
            raise ValueError, "A must be square"
        if A.nrows() != B.nrows():
            raise ValueError, "number of rows of A and B must be the same"
        if A.base_ring() != RDF:
            A = A.change_ring(RDF)
        if B.base_ring() != RDF:
            B = B.change_ring(RDF)
        if len(pi) != A.nrows():
            raise ValueError, "length of pi must equal number of rows of A"

        cdef Py_ssize_t i, j, k

        self.A = A
        self.B = B
        self.set_emission_symbols(emission_symbols)

        self.m = <ghmm_dmodel*> safe_malloc(sizeof(ghmm_dmodel))

        self.m.label = to_int_array(range(len(self._emission_symbols)))

        # Set all pointers to NULL
        self.m.s = NULL; self.m.name = NULL; self.m.silent = NULL
        self.m.tied_to = NULL; self.m.order = NULL; self.m.background_id = NULL
        self.m.bp = NULL; self.m.topo_order = NULL; self.m.pow_lookup = NULL;
        self.m.label_alphabet = NULL; self.m.alphabet = NULL

        # Set number of states and number of outputs
        self.m.N = A.nrows()
        self.m.M = len(self._emission_symbols)
        # Set the model type to discrete
        self.m.model_type = GHMM_kDiscreteHMM

        # Set that no a prior model probabilities are set.
        self.m.prior = -1
        # Assign model identifier if specified
        if name is not None:
            name = str(name)
            self.m.name = name
        else:
            self.m.name = NULL

        # Allocate and initialize states
        cdef ghmm_dstate* states = <ghmm_dstate*> safe_malloc(sizeof(ghmm_dstate) * self.m.N)
        cdef ghmm_dstate* state

        silent_states = []
        tmp_order     = []

        for i in range(self.m.N):
            v = B[i]

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
                    B.ncols(), len(emission_symbols))

            state.b = to_double_array(v)
            state.pi = pi[i]

            silent_states.append( 1 if sum(v) == 0 else 0)

            # Set out probabilities, i.e., the probabilities that each
            # symbol will be emitted from this state.
            v = A[i]
            nz = v.nonzero_positions()
            k = len(nz)
            state.out_states = k
            state.out_id = <int*> safe_malloc(sizeof(int)*k)
            state.out_a  = <double*> safe_malloc(sizeof(double)*k)
            for j in range(k):
                state.out_id[j] = nz[j]
                state.out_a[j]  = v[nz[j]]

            # Set "in" probabilities
            v = A.column(i)
            nz = v.nonzero_positions()
            k = len(nz)
            state.in_states = k
            state.in_id = <int*> safe_malloc(sizeof(int)*k)
            state.in_a  = <double*> safe_malloc(sizeof(double)*k)
            for j in range(k):
                state.in_id[j] = nz[j]
                state.in_a[j]  = v[nz[j]]

            state.fix = 0

        self.m.s = states
        self.initialized = True; return

        if sum(silent_states) > 0:
            self.m.model_type |= GHMM_kSilentStates
            self.m.silent = to_int_array(silent_states)

        self.m.maxorder = max(tmp_order)
        if self.m.maxorder > 0:
            self.m.model_type |= GHMM_kHigherOrderEmissions
            self.m.order = to_int_array(tmp_order)

        # Initialize lookup table for powers of the alphabet size,
        # which speeds up models with higher order states.
        powLookUp = [1] * (self.m.maxorder+2)
        for i in range(1,len(powLookUp)):
            powLookUp[i] = powLookUp[i-1] * self.m.M
        self.m.pow_lookup = to_int_array(powLookUp)

        self.initialized = True

    def __dealloc__(self):
        if self.initialized:
            ghmm_dmodel_free(&self.m)

    def __repr__(self):
        """
        Return string representation of this HMM.

        OUTPUT:
            string

        EXAMPLES:
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[0.9,0.1],[0.1,0.9]], [0.5,0.5], [3/4, 'abc'])
            sage: a.__repr__()
            "Discrete Hidden Markov Model (2 states, 2 outputs)\n\nInitial probabilities: [0.5, 0.5]\nTransition matrix:\n[0.1 0.9]\n[0.1 0.9]\nEmission matrix:\n[0.9 0.1]\n[0.1 0.9]\nEmission symbols: [3/4, 'abc']"
        """
        s = "Discrete Hidden Markov Model%s (%s states, %s outputs)\n"%(
            ' ' + self.m.name if self.m.name else '',
            self.m.N, self.m.M)
        s += '\nInitial probabilities: %s'%self.initial_probabilities()
        s += '\nTransition matrix:\n%s'%self.transition_matrix()
        s += '\nEmission matrix:\n%s'%self.emission_matrix()
        if self._emission_symbols_dict:
            s += '\nEmission symbols: %s'%self._emission_symbols
        return s

    def initial_probabilities(self):
        """
        Return the list of initial state probabilities.

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
            Discrete Hidden Markov Model (2 states, 2 outputs)
            Initial probabilities: [0.076923076923076927, 0.92307692307692302]
            Transition matrix:
            [0.333333333333 0.666666666667]
            [0.571428571429 0.428571428571]
            Emission matrix:
            [0.666666666667 0.333333333333]
            [0.333333333333 0.666666666667]
        """
        ghmm_dmodel_normalize(self.m)

    def sample_single(self, long length):
        """
        Return a single sample computed using this Hidden Markov Model.

        EXAMPLE:
            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[1,0],[0,1]], [0,1])
            sage: a.sample_single(20)
            [1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: a.sample_single(1000).count(0)
            113

        If the emission symbols are set
            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.1,0.9]], [[1,0],[0,1]], [0,1], ['up', 'down'])
            sage: print a.sample_single(10)
            ['down', 'up', 'down', 'down', 'up', 'down', 'down', 'up', 'down', 'up']

        """
        seed = random()
        cdef ghmm_dseq *d = ghmm_dmodel_generate_sequences(self.m, seed, length, 1, length)
        cdef Py_ssize_t i
        v = [d.seq[0][i] for i in range(length)]
        if self._emission_symbols_dict:
            w = self._emission_symbols
            return [w[i] for i in v]
        else:
            return v

    def sample(self, long length, long number):
        """
        Return number samples from this HMM of given length.

        EXAMPLES:
            sage: set_random_seed(0)
            sage: a = hmm.DiscreteHiddenMarkovModel([[0.1,0.9],[0.1,0.9]], [[1,0],[0,1]], [0,1])
            sage: print a.sample(10, 3)
            [[1, 0, 1, 1, 0, 1, 1, 0, 1, 0], [1, 1, 1, 1, 1, 1, 1, 1, 1, 0], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        """
        seed = random()
        cdef ghmm_dseq *d = ghmm_dmodel_generate_sequences(self.m, seed, length, number, length)
        cdef Py_ssize_t i, j
        v = [[d.seq[j][i] for i in range(length)] for j in range(number)]
        if self._emission_symbols_dict:
            w = self._emission_symbols
            return [[w[i] for i in z] for z in v]
        else:
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
            sage: a.sample_single(10)
            [5, 3, 5, 5, 3, 5, 5, 3, 5, 3]
            sage: a.set_emission_symbols([pi,5/9+e])
            sage: a.sample_single(10)
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
        HMM Problem 1: Return $\log ( P[seq | model] )$, the log of
        the probability of seeing the given sequence given this model,
        using the forward algorithm and assuming independance of the
        sequence seq.

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
        HMM Problem 2: Determine a hidden sequence of states that is
        most likely to produce the given sequence seq of obserations.

        INPUT:
            seq -- sequence of emitted symbols

        OUTPUT:
            list -- a most probable sequence of hidden states, i.e., the
                    Viterbi path.
            float -- log of the probability that the sequence of hidden
                     states actually produced the given sequence seq.
                     [[TODO: I do not understand precisely what this means.]]

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
        if self._emission_symbols_dict:
            seq = [self._emission_symbols_dict[z] for z in seq]
        cdef int* path
        cdef int* O = to_int_array(seq)
        cdef int pathlen
        cdef double log_p

        path = ghmm_dmodel_viterbi(self.m, O, len(seq), &pathlen, &log_p)
        sage_free(O)
        p = [path[i] for i in range(pathlen)]
        sage_free(path)

        return p, log_p

    ####################################################################
    # HMM Problem 3 -- Learning: Given an observation sequence O and
    # the set of states in the HMM, improve the HMM to increase the
    # probability of observing O.
    ####################################################################
    def baum_welch(self, training_seqs, nsteps=None, log_likelihood_cutoff=None):
        pass
## /** Baum-Welch-Algorithm for parameter reestimation (training) in
##     a discrete (discrete output functions) HMM. Scaled version
##     for multiple sequences, alpha and beta matrices are allocated with
## 	ighmm_cmatrix_stat_alloc
##     New parameters set directly in hmm (no storage of previous values!).
##     For reference see:
##     Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
##                 Applications in Speech Recognition"', Proceedings of the IEEE,
## 	77, no 2, 1989, pp 257--285
##   @return            0/-1 success/error
##   @param mo          initial model
##   @param sq          training sequences
##   */
## /** Just like reestimate_baum_welch, but you can limit
##     the maximum number of steps
##   @return            0/-1 success/error
##   @param mo          initial model
##   @param sq          training sequences
##   @param max_step    maximal number of Baum-Welch steps
##   @param likelihood_delta minimal improvement in likelihood required for carrying on. Relative value
##   to log likelihood
##   */







##################################################################################
# Helper Functions
##################################################################################

cdef double* to_double_array(v) except NULL:
    cdef double x
    cdef double* w = <double*> safe_malloc(sizeof(double)*len(v))
    cdef Py_ssize_t i = 0
    for x in v:
        w[i] = x
        i += 1
    return w

cdef int* to_int_array(v) except NULL:
    cdef int x
    cdef int* w = <int*> safe_malloc(sizeof(int)*len(v))
    cdef Py_ssize_t i = 0
    for x in v:
        w[i] = x
        i += 1
    return w

cdef void* safe_malloc(int bytes) except NULL:
    """
    malloc the given bytes of memory and check that the malloc
    succeeds -- if not raise a MemoryError.
    """
    cdef void* t = sage_malloc(bytes)
    if not t:
        raise MemoryError, "error allocating memory for Hidden Markov Model"
    return t

