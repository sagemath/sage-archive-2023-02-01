cdef extern from *:
    ctypedef int uint_fast8_t "uint_fast8_t"
    ctypedef int uint_fast8 "uint_fast8"
    ctypedef int uint_fast32_t "uint_fast32_t"
    ctypedef int subset_t "uint_least32_t"

cdef extern from "rankwidth/rw.h":
    int init_rw_dec(uint_fast8_t n)
    void destroy_rw()
    void calculate_level(uint_fast8_t subset_size)
    uint_fast8_t get_rw()
    subset_t *cslots
    subset_t *slots
    subset_t *adjacency_matrix

cdef extern from "rankwidth/rw.c":
    uint_fast8_t subset_size
    uint_fast8_t num_vertices


cdef void set_am(int i, int j, int val)
cdef void print_rank_dec(subset_t s, int l)



