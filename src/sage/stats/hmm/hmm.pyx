r"""
Hidden Markov Models

AUTHOR: William Stein

EXAMPLES:

"""

import math

from sage.matrix.all import is_Matrix, matrix
from sage.rings.all  import RDF

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


    int ghmm_dmodel_normalize(ghmm_dmodel *m)
    int ghmm_dmodel_free(ghmm_dmodel **m)


cdef class HiddenMarkovModel:
    pass


cdef class DiscreteHiddenMarkovModel:
    cdef ghmm_dmodel* m
    cdef bint initialized
    cdef object emission_symbols
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

        if emission_symbols is None:
            self.emission_symbols = range(B.ncols())
        else:
            self.emission_symbols = list(emission_symbols)

        self.m = <ghmm_dmodel*> sage_malloc(sizeof(ghmm_dmodel))
        if not self.m: raise MemoryError

        # Set all pointers to NULL
        self.m.s = NULL; self.m.name = NULL; self.m.silent = NULL
        self.m.tied_to = NULL; self.m.order = NULL; self.m.background_id = NULL
        self.m.bp = NULL; self.m.topo_order = NULL; self.m.pow_lookup = NULL;
        self.m.label = NULL; self.m.label_alphabet = NULL; self.m.alphabet = NULL

        # Set number of states and number of outputs
        self.m.N = A.nrows()
        self.m.M = len(self.emission_symbols)
        # Set the model type to discrete
        self.m.model_type = GHMM_kDiscreteHMM

        self.A = A
        self.B = B

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


    def __repr__(self):
        s = "Discrete Hidden Markov Model%s (%s states, %s outputs)\n"%(
            ' ' + self.m.name if self.m.name else '',
            self.m.N, self.m.M)
        s += 'Initial probabilities: %s\n'%self.initial_probabilities()
        s += 'Transition matrix:\n%s\n'%self.transition_matrix()
        s += 'Emission matrix:\n%s\n'%self.emission_matrix()
        return s

    def initial_probabilities(self):
        cdef Py_ssize_t i
        return [self.m.s[i].pi for i in range(self.m.N)]

    def transition_matrix(self, list_only=True):
        cdef Py_ssize_t i, j
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.m.s[i].out_states:
                self.A.set_unsafe_double(i,j,self.m.s[i].out_a[j])
        return self.A

    def emission_matrix(self, list_only=True):
        cdef Py_ssize_t i, j
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.B._ncols:
                self.B.set_unsafe_double(i,j,self.m.s[i].b[j])
        return self.B

    def normalize(self):
        """
        Normalize the transition and emission probabilities, if applicable.
        """
        ghmm_dmodel_normalize(self.m)

    def __dealloc__(self):
        if self.initialized:
            ghmm_dmodel_free(&self.m)




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

