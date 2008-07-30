
from hmm cimport HiddenMarkovModel

cdef extern from "ghmm/ghmm.h":
    cdef int GHMM_kContinuousHMM

cdef extern from "ghmm/smodel.h":
    #############################################################
    # Continuous emission density types
    #############################################################
    cdef enum ghmm_density_t:
        normal,
        normal_right,   # right tail
        normal_approx,
        normal_left,    # left tail
        uniform,
        binormal,
        multinormal,
        density_number

    #############################################################
    # Continuous emission
    #############################################################
    # Mean and variance types
    cdef union mean_t:
        double val
        double *vec

    cdef union variance_t:
        double val
        double *mat

    cdef struct ghmm_c_emission:
        # Flag for density function for each component of the mixture
        #   0: normal density
        #   1: truncated normal (right side) density
        #   2: approximated normal density
        #   3: truncated normal (left side)
        #   4: uniform distribution
        #   6: multivariate normal
        int type
        # dimension > 1 for multivariate normals
        int dimension
        # mean for output functions (pointer to mean vector for multivariate)
        mean_t mean
        variance_t variance
        # pointer to inverse of covariance matrix if multivariate normal;
        # else NULL
        double *sigmainv
        # determinant of covariance matrix for multivariate normal
        double det
        # Cholesky decomposition of covariance matrix A,
        #   if A = G*G' sigmacd only holds G
        double *sigmacd
        # minimum of uniform distribution or left boundary for
        # right-tail gaussians
        double min
        # maximum of uniform distribution or right boundary for
        # left-tail gaussians
        double max
        # if fixed != 0 the parameters of the density are fixed
        int fixed

    #############################################################
    # Continous Emission States
    #############################################################
    ctypedef struct ghmm_cstate:
        # Number of output densities per state
        int M
        # initial prob.
        double pi
        # IDs of successor states
        int *out_id
        # IDs of predecessor states
        int *in_id
        # transition probs to successor states. It is a
        # matrix in case of mult. transition matrices (COS > 1)
        double **out_a
        # transition probs from predecessor states. It is a
        # matrix in case of mult. transition matrices (COS > 1)
        double **in_a
        # number of  successor states
        int out_states
        # number of  predecessor states
        int in_states
        # weight vector for output function components
        double *c
        # flag for fixation of parameter. If fix = 1 do not change parameters of
        # output functions, and if fix = 0 do normal training. Default is 0.
        int fix
        # vector of ghmm_c_emission
        # (type and parameters of output function components)
        ghmm_c_emission *e
        # contains a description of the state (null terminated utf-8)
        char *desc
        # x coordinate position for graph representation plotting
        int xPosition
        # y coordinate position for graph representation plotting
        int yPosition

    #############################################################
    # Continous Hidden Markov Model
    #############################################################
    ctypedef struct ghmm_cmodel

    ctypedef struct ghmm_cmodel_class_change_context:
        # Names of class change module/function (for python callback)
        char *python_module
        char *python_function
        # index of current sequence
        int k
        # pointer to class function
        int (*get_class) (ghmm_cmodel*, double*, int, int)
        # space for any data necessary for class switch, USER is RESPONSIBLE
        void *user_data

    ctypedef struct ghmm_cmodel:
        # Number of states
        int N
        # Maximun number of components in the states
        int M
        # Number of dimensions of the emission components.
        # NOTE: All emissions must have the same number of dimensions.
        int dim
        # The ghmm_cmodel includes two continuous models
        #   cos=1: one transition matrix
        #   cos>1: (integer) extension for models with several matrices.
        int cos
        # Prior for a priori probability of the model.
        # -1 means no prior specified (all models have equal prob. a priori.)
        double prior
        # Contains a arbitrary name for the model (null terminated utf-8)
        char * name
        # Contains bit flags for various model extensions such as
        # kSilentStates (see ghmm.h for a complete list).
        int model_type
        # All states of the model. Transition probabilities are part of the states.
        ghmm_cstate *s
        # pointer to a ghmm_cmodel_class_change_context struct
        # necessary for multiple transition classes
        ghmm_cmodel_class_change_context *class_change

    int ghmm_cmodel_class_change_alloc(ghmm_cmodel * smo)

    #############################################################
    # Continous Sequence of states
    #
    # Sequence structure for double sequences.  Contains an array of
    # sequences and corresponding data like sequence labels, sequence
    # weights, etc. Sequences may have different length.
    # multi-dimension sequences are linearized
    #############################################################
    ctypedef struct ghmm_cseq:
        # sequence array. sequence[i][j] = j-th symbol of i-th seq.
        # sequence[i][D * j] = first dimension of j-th observation of i-th sequence
        double **seq
        # array of sequence length
        int *seq_len
        # array of sequence IDs
        double *seq_id
        # positive! sequence weights.  default is 1 = no weight
        double *seq_w
        # total number of sequences
        long seq_number
        # reserved space for sequences is always >= seq_number
        long capacity
        # sum of sequence weights
        double total_w
        # total number of dimensions
        int dim
        # flags (internal)
        unsigned int flags


    #############################################################
    # Continous Baum-Welch algorithm context
    #############################################################
    ctypedef struct ghmm_cmodel_baum_welch_context:
        # pointer of continuous model
        ghmm_cmodel *smo
        # sequence pointer
        ghmm_cseq *sqd
        # calculated log likelihood
        double *logp
        # leave reestimation loop if diff. between successive logp values
        # is smaller than eps
        double eps
        # max. no of iterations
        int max_iter


    int ghmm_cmodel_free (ghmm_cmodel ** smo)


cdef class ContinuousHiddenMarkovModel(HiddenMarkovModel):
    cdef ghmm_cmodel* m
    cdef bint initialized

cdef class NormalHiddenMarkovModel(ContinuousHiddenMarkovModel):
    pass


