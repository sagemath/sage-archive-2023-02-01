
include '../ext/cdefs.pxi'

cdef class BinaryCode:
    cdef unsigned int *basis
#    cdef unsigned int *columns
    cdef unsigned int *words
    cdef int ncols
    cdef int nrows
    cdef int radix
    cdef unsigned int nwords
    cdef int is_one(self, unsigned int, int)
    cdef int is_automorphism(self, int *, unsigned int *)

cdef class OrbitPartition:
    cdef unsigned int nwords
    cdef int ncols
    cdef unsigned int *wd_parent
    cdef unsigned int *wd_rank
    cdef unsigned int *wd_min_cell_rep
    cdef unsigned int *wd_size
    cdef int *col_parent
    cdef int *col_rank
    cdef int *col_min_cell_rep
    cdef int *col_size

    cdef unsigned int wd_find(self, unsigned int)
    cdef void wd_union(self, unsigned int, unsigned int)
    cdef int col_find(self, int)
    cdef void col_union(self, int, int)
    cdef int merge_perm(self, int *, unsigned int *)

cdef class PartitionStack:
    cdef unsigned int *wd_ents
    cdef int *wd_lvls
    cdef int *col_ents
    cdef int *col_lvls
    cdef unsigned int nwords
    cdef int ncols
    cdef int radix
    cdef unsigned int *col_degs   #
    cdef int *col_counts          #
    cdef int *col_output          #
    cdef int *wd_degs             #
    cdef unsigned int *wd_counts  # These are just for scratch space...
    cdef unsigned int *wd_output  #

    cdef int is_discrete(self, int)
    cdef int num_cells(self, int)
    cdef int sat_225(self, int)
    cdef unsigned int min_cell_reps(self, int)
    cdef unsigned int fixed_cols(self, unsigned int, int)
    cdef unsigned int first_smallest_nontrivial(self, int)
    cdef void col_percolate(self, int, int)
    cdef void wd_percolate(self, unsigned int, unsigned int)
    cdef int split_column(self, int, int)
    cdef unsigned int col_degree(self, BinaryCode, int, unsigned int, int)
    cdef int wd_degree(self, BinaryCode, unsigned int, int, int)
    cdef int sort_cols(self, int, int)
    cdef unsigned int sort_wds(self, unsigned int, int)

################################################################################
################################################################################
################################################################################

    cdef unsigned int refine(self, int, unsigned int *, int *, BinaryCode)
    cdef void get_permutation(self, PartitionStack, int *, int *)
    cdef int cmp(self, PartitionStack, BinaryCode)

cdef class BinaryCodeClassifier:
    cdef int *ham_wts


















