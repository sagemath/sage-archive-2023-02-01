cdef extern from "cliquer/graph.h":
    struct graph_t:
       pass

cdef extern from "cliquer/cliquer.h":
    struct clique_options:
       pass


cdef extern from "cliquer/reorder.h":
     cdef int *reorder_by_greedy_coloring(graph_t *g, bool weighted)
     cdef int *reorder_by_degree(graph_t *g, bool weighted)

cdef extern from "cliquer/cliquer.h":
     bool clique_print_time(intlevel, int i, int n, int max, double cputime, double realtime, clique_options *opts)


cdef extern from "cliquer/cl.h":
     cdef int sage_clique_max(graph_t *g, int ** list)
     cdef int sage_all_clique_max(graph_t *g, int ** list)
     cdef int sage_clique_number(graph_t *g)
     cdef void parse_input(char *str,graph_t *g)

cdef extern from "cliquer/graph.h":
     cdef graph_t * graph_new(int n)
     cdef void graph_print(graph_t *g)

