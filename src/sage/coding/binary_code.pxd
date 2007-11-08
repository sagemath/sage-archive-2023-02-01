
include '../ext/cdefs.pxi'

cdef class BinaryCode:
    cdef unsigned int *columns
    cdef int *ham_wts
    cdef int ncols
    cdef int nrows
    cdef int radix
    cdef unsigned int nwords
    cdef int is_one(self, unsigned int, int)
    cdef int is_automorphism(self, int *, unsigned int *)

cdef class OrbitPartition:
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
    cdef int merge_perm(self, int *, unsigned int *, int, int)

cdef class PartitionStack:
    cdef int *wd_ents
    cdef int *wd_lvls
    cdef int *col_ents
    cdef int *col_lvls
    cdef int nwords
    cdef int ncols
    cdef int is_discrete(self, int)
    cdef int num_cells(self, int)
    cdef int sat_225(self, int)
    cdef int is_min_cell_rep(self, int, int, int)
    cdef int is_fixed(self, int, int, int)
    cdef void col_percolate(self, int start, int end)
    cdef void wd_percolate(self, int start, int end)
    cdef int split_vertex(self, int, int, int)
    cdef int col_degree(self, BinaryCode, int, int, int)
    cdef int wd_degree(self, BinaryCode, int, int, int)
    cdef int sort_cols(self, int, int *, int)
    cdef int sort_wds(self, int, int *, int)
    cdef int refine(self, int, int *, int *, BinaryCode)
    cdef void get_permutation(self, PartitionStack, int *, int *)
    cdef int cmp(self, PartitionStack, BinaryCode)