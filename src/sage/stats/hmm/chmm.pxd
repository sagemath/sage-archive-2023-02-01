#############################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version at your option.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

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
        # matrix in case of mult. transition matrices (cos > 1)
        double **out_a
        # transition probs from predecessor states. It is a
        # matrix in case of mult. transition matrices (cos > 1)
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


    double **ighmm_cmatrix_alloc (int nrows, int ncols)

    #############################################################
    # Continous Hidden Markov Model
    #############################################################
    ctypedef struct ghmm_cmodel
    ctypedef struct ghmm_cmodel_class_change_context

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
        # array of sequence lengths
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
        # obsolete
        int* seq_label

    int ghmm_cseq_free (ghmm_cseq ** sq)

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

    # Normalizes initial and transition probabilities and mixture weights.
    # 0 on success / -1 on error
    int ghmm_cmodel_normalize(ghmm_cmodel *smo)


    # Produces sequences computed using a given model. All memory that
    # is needed for the sequences is allocated inside the function. It
    # is possible to define the length of the sequences global
    # (global_len > 0) or it can be set inside the function, when a
    # final state in the model is reached (a state with no output). If
    # the model has no final state, the sequences will have length
    # MAX_SEQ_LEN.
    #
    #  INPUT:
    #      smo       -- model
    #      seed      -- initial parameter for the random number generator;
    #                   if seed == 0, don't reinitialize random number generator
    #      global_en -- length of sequence.  If 0, sequence ended automatically
    #                   when the final state is reached
    #      seq_number-- number of sequences
    #      Tmax      -- maximal sequence length (or -1)
    #
    #  OUTPUT:
    #      pointer to an array of sequences
    #

    ghmm_cseq *ghmm_cmodel_generate_sequences (ghmm_cmodel * smo, int seed,
                                             int global_len, long seq_number,
                                             int Tmax)

    # Computes sum over all sequence of seq_w * log( P ( O|lambda )).
    # If a sequence can't be generated by smo error cost of seq_w *
    # PENALTY_LOGP are imposed.
    #   OUTPUT: number of evaluated sequences; -1: error
    int ghmm_cmodel_likelihood (ghmm_cmodel * smo, ghmm_cseq * sqd, double *log_p)


    # Viterbi algorithm: calculation of the viterbi path (best possible
    # state sequence for a given sequence for a given model (smo)). Also
    # calculates log(p) according to this path, the matrices in the local_store
    # struct are allocated using stat_matrix_d_alloc.
    int *ghmm_cmodel_viterbi (ghmm_cmodel * smo, double *o, int n, double *log_p)

cdef extern from "ghmm/sreestimate.h":
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

    # Baum-Welch Algorithm for continuous HMMs
    # Training of model parameter with multiple double sequences (incl. scaling).
    # New parameters set directly in hmm (no storage of previous values!). Matrices
    # are allocated with stat_matrix_d_alloc.
    int ghmm_cmodel_baum_welch (ghmm_cmodel_baum_welch_context * cs)


cdef class ContinuousHiddenMarkovModel(HiddenMarkovModel):
    cdef ghmm_cmodel* m
    cdef bint initialized

cdef class GaussianHiddenMarkovModel(ContinuousHiddenMarkovModel):
    pass


