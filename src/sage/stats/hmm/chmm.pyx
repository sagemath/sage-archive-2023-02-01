

include "../../ext/stdsage.pxi"

include "misc.pxi"

cdef class ContinuousHiddenMarkovModel(HiddenMarkovModel):
    def __init__(self, A, B, pi, name=None):
        self.initialized = False
        HiddenMarkovModel.__init__(self, A, B, pi)
        self.m = <ghmm_cmodel*> safe_malloc(sizeof(ghmm_cmodel))
        # Set number of states
        self.m.N = self.A.nrows()

cdef class NormalHiddenMarkovModel(ContinuousHiddenMarkovModel):
    """
    EXAMPLES:
    The transition matrix:
        sage: A = [[0,1,0],[0.5,0,0.5],[0.3,0.3,0.4]]

    Parameters of the normal emission distributions in pairs of (mu, sigma):
        sage: B = [(0,1), (-1,0.5), (1,0.2)]

    The initial probabilities per state:
        sage: pi = [1,0,0]

    Create the continuous HMM:
        sage: m = hmm.NormalHiddenMarkovModel(A, B, pi); m
    """
    def __init__(self, A, B, pi, name=None):
        ContinuousHiddenMarkovModel.__init__(self, A, B, pi)

        # Set number of outputs
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
            mu, sigma = self.B[i]
            # Get a reference to the i-th state for convenience of the notation below.
            state = &(states[i])
            state.desc = NULL
            state.M = 1
            e = <ghmm_c_emission*> safe_malloc(sizeof(ghmm_c_emission))
            e.type = 0  # normal
            e.dimension = 1
            e.mean.val = mu
            e.variance.val = sigma
            # fixing of emissions is deactivated by default
            e.fixed = 0
            e.sigmacd = NULL
            e.sigmainv = NULL
            state.e = e
            state.c = to_double_array([1.0])

        # Set states
        self.m.s = states

        self.initialized = True

    def __dealloc__(self):
        return
        if self.initialized:
            ghmm_cmodel_free(&self.m)

    def __repr__(self):
        """
        Return string representation of this Continuous HMM.

        OUTPUT:
            string

        EXAMPLES:
            sage: m = hmm.NormalHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [1,0,0])
            sage: a.__repr__()
            "Discrete Hidden Markov Model (2 states, 2 outputs)\nInitial probabilities: [0.5, 0.5]\nTransition matrix:\n[0.1 0.9]\n[0.1 0.9]\nEmission matrix:\n[0.9 0.1]\n[0.1 0.9]\nEmission symbols: [3/4, 'abc']"
        """
        s = "Normal Hidden Markov Model%s (%s states, %s outputs)"%(
            ' ' + self.m.name if self.m.name else '',
            self.m.N, self.m.M)
        #s += '\nInitial probabilities: %s'%self.initial_probabilities()
        #s += '\nTransition matrix:\n%s'%self.transition_matrix()
        #s += '\nEmission matrix:\n%s'%self.emission_matrix()
        return s

