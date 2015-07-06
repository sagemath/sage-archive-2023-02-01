#ifndef BOOSTGRAPH
#define BOOSTGRAPH
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/edge_connectivity.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/clustering_coefficient.hpp>

#include <iostream>

using namespace std;
using namespace boost;

// This struct is the output of the edge connectivity Boost algorithm.
typedef struct {
    int ec; // The edge connectivity
    vector<int> edges; // The edges in a minimum cut, stored as a list of
                       // nodes. For instance, if the minimum cut is
                       // {(1,2),(3,4)}, the output vector will be (1,2,3,4).
} result_ec;


// This struct is the output of the clustering coefficient Boost algorithm.
typedef struct {
    double average_clustering_coefficient; // The average clustering coefficient
    vector<double> clust_of_v;             // The clustering coefficient of each node.
} result_cc;

template <class OutEdgeListS, // How neighbors are stored
          class VertexListS,  // How vertices are stored
          class DirectedS,    // The kind of graph (undirectedS, directedS, or bidirectionalS)
          class EdgeListS>    // How the list of edges is stored
class BoostGraph
/*
 * This generic class wraps a Boost graph, in order to make it Cython-friendly.
 *
 * In particular, it allows to "keep together" the Boost graph and the vector
 * *vertices: these two variables are generic, and Cython is not able to deal
 * with them properly, since it does not support generic classes.
 *
 * Vertices are numbers from 0 to n-1, where n is the total number of vertices:
 * this class takes care of the relation between number i and the corresponding
 * Boost vertex descriptor (which might be any object, depending on the value of
 * VertexListS). In particular, (*vertices)[i] contains the Boost vertex
 * corresponding to number i, while to transform a Boost vertex v into a number
 * we use vertex properties, and the syntax is (*graph)[v].
*/
{
    typedef typename boost::adjacency_list<OutEdgeListS, VertexListS, DirectedS,
    int, no_property, no_property, EdgeListS> adjacency_list;
    typedef typename boost::graph_traits<adjacency_list>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<adjacency_list>::edge_descriptor edge_descriptor;
    typedef typename std::vector<edge_descriptor> edge_container;

public:
    adjacency_list *graph;
    vector<vertex_descriptor> *vertices;

    BoostGraph() {
        graph = new adjacency_list();
        vertices = new vector<vertex_descriptor>();
    }

    ~BoostGraph() {
        delete graph;
        delete vertices;
    }

    int num_verts() {
        return num_vertices(*graph);
    }

    int num_edges() {
        return boost::num_edges(*graph);
    }

    void add_vertex() {
        (*vertices).push_back(boost::add_vertex((*vertices).size(), *graph));
    }

    void add_edge(int u, int v) {
        boost::add_edge((*vertices)[u], (*vertices)[v], *graph);
    }

    result_ec edge_connectivity() {
        result_ec to_return;
        edge_container disconnecting_set;
        back_insert_iterator<edge_container> inserter(disconnecting_set);
        to_return.ec = boost::edge_connectivity(*graph, inserter);

        for (unsigned i = 0; i < disconnecting_set.size(); i++) {
            edge_descriptor edge = disconnecting_set[i];
            to_return.edges.push_back((*graph)[boost::source(edge, *graph)]);
            to_return.edges.push_back((*graph)[boost::target(edge, *graph)]);
        }
        return to_return;
    }

    double clustering_coeff(int v) {
        return clustering_coefficient(*graph, (*vertices)[v]);
    }

    result_cc clustering_coeff_all() {
        result_cc to_return;

        typedef typename exterior_vertex_property<adjacency_list, double>::container_type ClusteringContainer;
        typedef typename exterior_vertex_property<adjacency_list, double>::map_type ClusteringMap;

        ClusteringContainer coefs(num_vertices(*graph));
        ClusteringMap cm(coefs, *graph);

        to_return.average_clustering_coefficient = all_clustering_coefficients(*graph, cm);
        for (unsigned i = 0; i < num_verts(); i++) {
            to_return.clust_of_v.push_back(cm[(*graph)[i]]);
        }

        return to_return;
    }
};



#endif // BOOSTGRAPH
