"""
Continuous Hidden Markov Models

AUTHORS:
    -- William Stein (2008)
    -- The authors of GHMM http://ghmm.sourceforge.net/
"""

#############################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version at your option.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

include "../../ext/stdsage.pxi"

include "misc.pxi"

cdef class ContinuousHiddenMarkovModel(HiddenMarkovModel):
    """
    Abstract base class for continuous hidden Markov models.


    INPUT:
        A -- square matrix (or list of lists)
        B -- list or matrix with numbers of rows equal to number of rows of A;
             each row defines the emissions
        pi -- list
        name -- (default: None); a string

    EXAMPLES:
        sage: sage.stats.hmm.chmm.ContinuousHiddenMarkovModel([[1.0]], [(-0.0,10.0)], [1], "model")
        <sage.stats.hmm.chmm.ContinuousHiddenMarkovModel object at ...>
    """
    def __init__(self, A, B, pi, name):
        """
        Constructor for continuous Hidden markov model abstract base class.

        EXAMPLES:
        This class is an abstract base class, so shouldn't ever be constructed by users.
            sage: sage.stats.hmm.chmm.ContinuousHiddenMarkovModel([[1.0]], [(0.0,1.0)], [1], None)
            <sage.stats.hmm.chmm.ContinuousHiddenMarkovModel object at ...>
        """
        self.initialized = False
        HiddenMarkovModel.__init__(self, A, B, pi)
        self.m = <ghmm_cmodel*> safe_malloc(sizeof(ghmm_cmodel))
        # Set number of states
        self.m.N = self.A.nrows()
        # Assign model identifier (the name) if specified
        if name is not None:
            name = str(name)
            self.m.name = <char*> safe_malloc(len(name))
            strcpy(self.m.name, name)
        else:
            self.m.name = NULL


    def name(self):
        """
        Return the name of this model.

        OUTPUT:
            string or None

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [(0.0,1.0),(1,1)], [1,0], 'Test Model')
            sage: m.name()
            'Test Model'

        If the model is not explicitly named then this function returns None:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [(0.0,1.0),(1,1)], [1,0])
            sage: m.name()
            sage: m.name() is None
            True
        """
        if self.m.name:
            s = str(self.m.name)
            return s
        else:
            return None


cdef class GaussianHiddenMarkovModel(ContinuousHiddenMarkovModel):
    r"""
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
        name -- (default: None); a string

    EXAMPLES:
    Define the transition matrix:
        sage: A = [[0,1,0],[0.5,0,0.5],[0.3,0.3,0.4]]

    Parameters of the normal emission distributions in pairs of (mu, sigma):
        sage: B = [(0,1), (-1,0.5), (1,0.2)]

    The initial probabilities per state:
        sage: pi = [1,0,0]

    We create the continuous Gaussian hidden Markov model defined by $A,B,\pi$:
        sage: hmm.GaussianHiddenMarkovModel(A, B, pi)
        Gaussian Hidden Markov Model with 3 States
        Transition matrix:
        [0.0 1.0 0.0]
        [0.5 0.0 0.5]
        [0.3 0.3 0.4]
        Emission parameters:
        [(0.0, 1.0), (-1.0, 0.5), (1.0, 0.20000000000000001)]
        Initial probabilities: [1.0, 0.0, 0.0]
    """
    def __init__(self, A, B, pi, name=None):
        """
        EXAMPLES:
        We make a very simple model:
            sage: hmm.GaussianHiddenMarkovModel([[1]], [(0,1)], [1], 'simple')
            Gaussian Hidden Markov Model simple with 1 States
            Transition matrix:
            [1.0]
            Emission parameters:
            [(0.0, 1.0)]
            Initial probabilities: [1.0]

        We test a degenerate case:
            sage: hmm.GaussianHiddenMarkovModel([], [], [], 'simple')
            Gaussian Hidden Markov Model simple with 0 States
            Transition matrix:
            []
            Emission parameters:
            []
            Initial probabilities: []
        """
        ContinuousHiddenMarkovModel.__init__(self, A, B, pi, name=name)

        # Set number of outputs.  This is 1 here because each
        # output is a single Gaussian distribution.
        self.m.M = 1

        # Set the model type to continuous
        self.m.model_type = GHMM_kContinuousHMM

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
        cdef Py_ssize_t i, j, k

        for i in range(self.m.N):
            # Parameters of normal distribution
            mu, sigma   = self.B[i]
            # Get a reference to the i-th state for convenience of the notation below.
            state = &(states[i])
            state.M     = 1
            state.pi    = pi[i]
            state.desc  = NULL
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

            #########################################################
            # Initialize state transition data.
            # NOTE: This code is similar to a block of code in hmm.pyx.

            # Set "out" probabilities, i.e., the probabilities to
            # transition to another hidden state from this state.
            v = self.A[i]
            k = self.m.N
            state.out_states = k
            state.out_id = <int*> safe_malloc(sizeof(int)*k)
            state.out_a  = ighmm_cmatrix_alloc(1, k)
            for j in range(k):
                state.out_id[j] = j
                state.out_a[0][j]  = v[j]

            # Set "in" probabilities
            v = self.A.column(i)
            state.in_states = k
            state.in_id = <int*> safe_malloc(sizeof(int)*k)
            state.in_a  = ighmm_cmatrix_alloc(1, k)
            for j in range(k):
                state.in_id[j] = j
                state.in_a[0][j]  = v[j]

            #########################################################


        # Set states
        self.m.s = states

        self.m.class_change = NULL

        self.initialized = True

    def __dealloc__(self):
        """
        Dealloc the memory used by this Gaussian HMM, but only if the
        HMM was successfully initialized.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[1.0]], [(0.0,1.0)], [1])   # implicit doctest
            sage: del m
        """
        if self.initialized:
            ghmm_cmodel_free(&self.m)

    def __copy__(self):
        """
        Return a copy of this Gaussian HMM.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [(0.0,1.0),(1,1)], [1,0], 'NAME')
            sage: copy(m)
            Gaussian Hidden Markov Model NAME with 2 States
            Transition matrix:
            [0.4 0.6]
            [0.1 0.9]
            Emission parameters:
            [(0.0, 1.0), (1.0, 1.0)]
            Initial probabilities: [1.0, 0.0]
        """

        return GaussianHiddenMarkovModel(self.transition_matrix(), self.emission_parameters(),
                                         self.initial_probabilities(), self.name())

    def __cmp__(self, other):
        """
        Compare two Gaussian HMM's.

        INPUT:
            self, other -- objects; if other is not a Gaussian HMM compare types.
        OUTPUT:
            -1,0,1

        The transition matrices are compared, then the emission
        parameters, and the initial probabilities.  The names are not
        compared, so two GHMM's with the same parameters but different
        names compare to be equal.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [(0.0,1.0),(1,1)], [1,2], "Test 1")
            sage: m.__cmp__(m)
            0

        Note that the name doesn't matter:
            sage: n = hmm.GaussianHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [(0.0,1.0),(1,1)], [1,2], "Test 2")
            sage: m.__cmp__(n)
            0

        Normalizing fixes the initial probabilities, hence m and n are no longer equal.
            sage: n.normalize()
            sage: m.__cmp__(n)
            1
        """
        if not isinstance(other, GaussianHiddenMarkovModel):
            return cmp(type(self), type(other))

        if self is other: return 0  # easy special case

        cdef GaussianHiddenMarkovModel o = other
        if self.m.N < o.m.N:
            return -1
        elif self.m.N > o.m.N:
            return 1
        cdef Py_ssize_t i, j

        # The code below is somewhat long and tedious because I want it to be
        # very fast.  All it does is explicitly loop through the transition
        # matrix, emission parameters and initial state probabilities checking
        # that they agree, and if not returning -1 or 1.
        # Compare transition matrices
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.m.s[i].out_states:
                if self.m.s[i].out_a[0][j] < o.m.s[i].out_a[0][j]:
                    return -1
                elif self.m.s[i].out_a[0][j] > o.m.s[i].out_a[0][j]:
                    return 1

        # Compare emissions parameters
        for i from 0 <= i < self.m.N:
            if self.m.s[i].e.mean.val < o.m.s[i].e.mean.val:
                return -1
            elif self.m.s[i].e.mean.val > o.m.s[i].e.mean.val:
                return 1
            if self.m.s[i].e.variance.val < o.m.s[i].e.variance.val:
                return -1
            elif self.m.s[i].e.variance.val > o.m.s[i].e.variance.val:
                return 1

        # Compare initial state probabilities
        for 0 <= i < self.m.N:
            if self.m.s[i].pi < o.m.s[i].pi:
                return -1
            elif self.m.s[i].pi > o.m.s[i].pi:
                return 1

        return 0

    def __repr__(self):
        """
        Return string representation of this Continuous HMM.

        OUTPUT:
            string

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [1,0,0])
            sage: m.__repr__()
            'Gaussian Hidden Markov Model with 3 States\nTransition matrix:\n[0.0 1.0 0.0]\n[0.5 0.0 0.5]\n[0.3 0.3 0.4]\nEmission parameters:\n[(0.0, 1.0), (-1.0, 0.5), (1.0, 0.20000000000000001)]\nInitial probabilities: [1.0, 0.0, 0.0]'
        """
        s = "Gaussian Hidden Markov Model%s with %s States"%(
            ' ' + self.m.name if self.m.name else '',
            self.m.N)
        s += '\nTransition matrix:\n%s'%self.transition_matrix()
        s += '\nEmission parameters:\n%s'%self.emission_parameters()
        s += '\nInitial probabilities: %s'%self.initial_probabilities()
        return s

    def initial_probabilities(self):
        """
        Return the list of initial state probabilities.

        OUTPUT:
            list of floats

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [0.4,0.3,0.3])
            sage: m.initial_probabilities()
            [0.40000000000000002, 0.29999999999999999, 0.29999999999999999]
        """
        cdef Py_ssize_t i
        return [self.m.s[i].pi for i in range(self.m.N)]

    def transition_matrix(self):
        """
        Return the hidden state transition matrix.

        OUTPUT:
            matrix whose rows give the transition probabilities of the
            Hidden Markov Model states.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.0,1.0,0],[0.5,0.0,0.5],[0.3,0.3,0.4]], [(0.0,1.0), (-1.0,0.5), (1.0,0.2)], [1,0,0])
            sage: m.transition_matrix()
            [0.0 1.0 0.0]
            [0.5 0.0 0.5]
            [0.3 0.3 0.4]
        """
        cdef Py_ssize_t i, j
        # Update the state of the "immutable" matrix A, then return a reference to it.
        for i from 0 <= i < self.m.N:
            for j from 0 <= j < self.m.s[i].out_states:
                self.A.set_unsafe_double(i,j,self.m.s[i].out_a[0][j])
        return self.A

    def emission_parameters(self):
        """
        Return the emission parameters list.

        OUTPUT:
            list of tuples (mu, sigma) that define Gaussian distributions associated to each state.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[0.4,0.6],[0.1,0.9]], [(1.5,2),(-1,3)], [1,0], 'NAME')
            sage: m.emission_parameters()
            [(1.5, 2.0), (-1.0, 3.0)]
        """
        cdef Py_ssize_t i
        return [(self.m.s[i].e.mean.val, self.m.s[i].e.variance.val) for i in range(self.m.N)]

    def normalize(self):
        """
        Normalize the transition and emission probabilities, if
        applicable.  This changes self in place.

        EXAMPLES:
            sage: m = hmm.GaussianHiddenMarkovModel([[1.0,1.2],[0,0.1]], [(0.0,1.0),(1,1)], [1,2])
            sage: m.normalize()
            sage: m
            Gaussian Hidden Markov Model with 2 States
            Transition matrix:
            [0.454545454545 0.545454545455]
            [           0.0            1.0]
            Emission parameters:
            [(0.0, 1.0), (1.0, 1.0)]
            Initial probabilities: [0.33333333333333331, 0.66666666666666663]
        """
        if ghmm_cmodel_normalize(self.m):
            raise RuntimeError, "error normalizing model (note that model may still have been partly changed)"
