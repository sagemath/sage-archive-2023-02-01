#############################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version at your option.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

from sage.matrix.matrix_real_double_dense cimport Matrix_real_double_dense

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
    int ghmm_dseq_free (ghmm_dseq ** sq)
    ghmm_dseq *ghmm_dseq_calloc (long seq_number)

    int ghmm_dmodel_normalize(ghmm_dmodel *m)
    int ghmm_dmodel_free(ghmm_dmodel **m)

    # Problem 1: Likelihood
    int ghmm_dmodel_logp(ghmm_dmodel *m, int *O, int len, double *log_p)

    # Problem 2: Decoding
    int *ghmm_dmodel_viterbi (ghmm_dmodel *m, int *O, int len, int *pathlen, double *log_p)

    # Problem 3: Learning
    int ghmm_dmodel_baum_welch (ghmm_dmodel *m, ghmm_dseq *sq)


cdef class HiddenMarkovModel:
    cdef Matrix_real_double_dense A, B
    cdef list pi

cdef class DiscreteHiddenMarkovModel(HiddenMarkovModel):
    cdef ghmm_dmodel* m
    cdef bint initialized
    cdef object _emission_symbols, _emission_symbols_dict
