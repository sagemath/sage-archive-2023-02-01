#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from libcpp.vector cimport vector

cdef extern from "boost/graph/adjacency_list.hpp" namespace "boost":
    cdef cppclass vecS:
        pass
    cdef cppclass listS:
        pass
    cdef cppclass slistS:
        pass
    cdef cppclass setS:
        pass
    cdef cppclass multisetS:
        pass
    cdef cppclass hash_setS:
        pass
    cdef cppclass undirectedS:
        pass
    cdef cppclass directedS:
        pass
    cdef cppclass bidirectionalS:
        pass

cdef extern from "boost_interface.cpp":
    ctypedef unsigned long v_index
    ctypedef unsigned long e_index

    cdef cppclass result_ec:
        v_index ec
        vector[v_index] edges

    cdef cppclass result_cc:
        float average_clustering_coefficient
        vector[float] clust_of_v

    cdef cppclass BoostGraph[OutEdgeListS, VertexListS, DirectedS, EdgeListS]:
        BoostGraph()
        void add_vertex()
        v_index num_verts()
        void add_edge(v_index u, v_index v)
        e_index num_edges()
        result_ec edge_connectivity()
        double clustering_coeff(v_index v)
        result_cc clustering_coeff_all()
        vector[v_index] dominator_tree(v_index v)

ctypedef BoostGraph[vecS, vecS, undirectedS,    vecS] BoostVecGraph
ctypedef BoostGraph[vecS, vecS, bidirectionalS, vecS] BoostVecDiGraph

ctypedef BoostGraph[setS, vecS, undirectedS,    vecS] BoostSetGraph

ctypedef fused BoostVecGenGraph:
    BoostVecGraph
    BoostVecDiGraph

ctypedef fused BoostGenGraph:
    BoostVecGraph
    BoostVecDiGraph
    BoostGraph[vecS, vecS, directedS,      vecS]
    BoostSetGraph
    BoostGraph[setS, vecS, directedS,      vecS]
    BoostGraph[setS, vecS, bidirectionalS, vecS]
