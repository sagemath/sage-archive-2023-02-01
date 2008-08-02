

include "../../ext/stdsage.pxi"

include "misc.pxi"

cdef class ContinuousHiddenMarkovModel(HiddenMarkovModel):
    def __init__(self, A, B, pi, name=None):
        self.initialized = False
        HiddenMarkovModel.__init__(self, A, B, pi)
        self.m = <ghmm_cmodel*> safe_malloc(sizeof(ghmm_cmodel))
        # Set number of states
        self.m.N = self.A.nrows()

cdef class GaussianHiddenMarkovModel(ContinuousHiddenMarkovModel):
    """
    Create a Gaussian Hidden Markov Model.  The probability
    distribution associated with each state is a Gaussian
    distribution.

    GaussianHiddenMarkovModel(A, B, pi, name)

    INPUT:
        A  -- matrix; the transition matrix (n x n)
        B  -- list of n pairs (mu, sigma) that define the
              Gaussian distributions associated to each state
        pi -- list of floats that sums to 1.0; these are
              the initial probabilities of each hidden state
        name -- (default: None);

    EXAMPLES:
    Define the transition matrix:
        sage: A = [[0,1,0],[0.5,0,0.5],[0.3,0.3,0.4]]

    Parameters of the normal emission distributions in pairs of (mu, sigma):
        sage: B = [(0,1), (-1,0.5), (1,0.2)]

    The initial probabilities per state:
        sage: pi = [1,0,0]

    Create the continuous Gaussian hidden Markov model:
        sage: m = hmm.GaussianHiddenMarkovModel(A, B, pi); m
    """
    def __init__(self, A, B, pi, name=None):
        ContinuousHiddenMarkovModel.__init__(self, A, B, pi)

        # Set number of outputs.  This is 1 here because each
        # output is a single Gaussian distribution.
        self.m.M = 1

        # Set the model type to continuous
        self.m.model_type = GHMM_kContinuousHMM

        # Assign model identifier if specified
        if name is not None:
            name = str(name)
            self.m.name = name
        else:
            self.m.name = NULL

        # 1 transition matrix
        self.m.cos   =  1
        # Set that no a prior model probabilities are set.
        self.m.prior = -1
        # Dimension is 1
        self.m.dim   =  1

        # Allocate and initialize states
        cdef ghmm_cstate* states = <ghmm_cstate*> safe_malloc(sizeof(ghmm_cstate) * self.m.N)
        cdef ghmm_cstate* state
        cdef ghmm_c_emission* e

        for i in range(self.m.N):
            # Parameters of normal distribution
            mu, sigma   = self.B[i]
            # Get a reference to the i-th state for convenience of the notation below.
            state = &(states[i])
            state.M     = 1
            state.pi    = pi[i]
            state.desc  = NULL
            state.out_states = 0
            state.in_states = 0
            e = <ghmm_c_emission*> safe_malloc(sizeof(ghmm_c_emission))
            e.type      = 0  # normal
            e.dimension = 1
            e.mean.val  = mu
            e.variance.val = sigma
            # fixing of emissions is deactivated by default
            e.fixed     = 0
            e.sigmacd   = NULL
            e.sigmainv  = NULL
            state.e     = e
            state.c     = to_double_array([1.0])
            state.in_a  = ighmm_cmatrix_alloc(1, self.m.N)
            state.out_a = ighmm_cmatrix_alloc(1, self.m.N)

        # Set states
        self.m.s = states

        self.m.class_change = NULL

        self.initialized = True

    def __dealloc__(self):
        if self.initialized:
            ghmm_cmodel_free(&self.m)

    def __repr__(self):
        """
        Return string representation of this Continuous HMM.

        OUTPUT:
            string

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [1,0,0])
            sage: a.__repr__()
            "Discrete Hidden Markov Model (2 states, 2 outputs)\nInitial probabilities: [0.5, 0.5]\nTransition matrix:\n[0.1 0.9]\n[0.1 0.9]\nEmission matrix:\n[0.9 0.1]\n[0.1 0.9]\nEmission symbols: [3/4, 'abc']"
        """
        s = "Gaussian Hidden Markov Model%s (%s states, %s outputs)"%(
            ' ' + self.m.name if self.m.name else '',
            self.m.N, self.m.M)
        s += '\nInitial probabilities: %s'%self.initial_probabilities()
        s += '\nTransition matrix:\n%s'%self.transition_matrix()
        s += '\nEmission parameters:\n%s'%self.emission_parameters()
        return s

    def initial_probabilities(self):
        """
        Return the list of initial state probabilities.

        OUTPUT:
            list of floats

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [0.4,0.3,0.3])
            sage: m.initial_probabilities()
            [0.4, 0.3, 0.3]
        """
        cdef Py_ssize_t i
        return [self.m.s[i].pi for i in range(self.m.N)]

    def transition_matrix(self, list_only=True):
        """
        Return the hidden state transition matrix.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [1,0,0])
            sage: m.transition_matrix()
            [0.9 0.1]
            [0.9 0.1]
        """
        cdef Py_ssize_t i, j
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.m.s[i].out_states:
                self.A.set_unsafe_double(i,j,self.m.s[i].out_a[0][j])
        return self.A

    def emission_parameters(self):
        """
        Return the emission probability matrix.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [0.1,0.4,0.5])
            sage: m.emission_parameters()
            [(0.0, 1.0), (-1.0, 0.5), (1.0, 0.20000...)]
        """
        cdef Py_ssize_t i, j
        v = []
        for i from 0 <= i < self.m.N:
            v.append((self.m.s[i].e.mean.val, self.m.s[i].e.variance.val))
        return v
