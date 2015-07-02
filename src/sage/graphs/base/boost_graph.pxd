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
    cdef cppclass result_ec:
        int ec
        vector[int] edges
    cdef cppclass BoostGraph[OutEdgeListS, VertexListS, DirectedS, EdgeListS]:
        BoostGraph()
        void add_vertex()
        int num_verts()
        void add_edge(int u, int v)
        int num_edges()
        result_ec edge_connectivity()

ctypedef BoostGraph[vecS, vecS, undirectedS,    vecS] BoostVecGraph
ctypedef BoostGraph[vecS, vecS, bidirectionalS, vecS] BoostVecDiGraph

ctypedef fused BoostVecGenGraph:
    BoostVecGraph
    BoostVecDiGraph

ctypedef fused BoostGenGraph:
    BoostVecGraph
    BoostVecDiGraph
    BoostGraph[vecS, vecS, directedS,      vecS]
    BoostGraph[setS, vecS, undirectedS,    vecS]
    BoostGraph[setS, vecS, directedS,      vecS]
    BoostGraph[setS, vecS, bidirectionalS, vecS]
