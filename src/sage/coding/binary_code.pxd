
cdef class BinaryCodeGraph:

    cdef int *columns
    cdef int *ham_wts
    cdef int ncols
    cdef int nrows
    cdef int radix
    cdef int ptn_mask
    cdef int nwords
    cdef int has_edge_bip(self, int, int)
    cdef int has_edge(self, int, int)

cdef class PartitionStack:
    cdef int *wd_ents
    cdef int *col_ents
    cdef int *wd_lvls
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
    cdef int col_degree(self, BinaryCodeGraph, int, int, int)
    cdef int wd_degree(self, BinaryCodeGraph, int, int, int)
    cdef int sort_cols(self, int, int *, int)
    cdef int sort_wds(self, int, int *, int)
    cdef int refine(self, int, int *, int *, BinaryCodeGraph)
    cdef void get_permutation(self, PartitionStack, int *, int *)
    cdef int cmp(self, PartitionStack, BinaryCodeGraph)