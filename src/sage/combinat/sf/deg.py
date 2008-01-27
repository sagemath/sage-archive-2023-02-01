"""
Dual equivalence graphs
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.combinatorial_algebra import CombinatorialAlgebra
from sage.combinat.word import Words
from sage.combinat.tableau import StandardTableaux, Tableau
from sage.combinat.permutation import Permutation
from sage.combinat.partition import Partition, Partitions
from sage.graphs.graph import Graph
from sage.rings.all import QQ

###################################
#Quasisymmetric function utilities#
###################################
class GesselQuasisymmetricFunctions(CombinatorialAlgebra):
    def __init__(self, R):
        self._combinatorial_class = Words()
        self._name = "Quasisymmetric functions"
        self._one = []
        self._prefix = "Q"
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        raise NotImplementedError

Q = GesselQuasisymmetricFunctions(QQ)

def is_positive(q):
    """
    Takes in a combinatorial algebra element q, and returns True
    if all the coefficients of q are positive.
    """
    return all(map(lambda x: x > 0, q.monomial_coefficients().values()))

def q_to_s(qsum, n):
    """
    Given a Gessel quasisymmetric function, return it as a Schur symmetric
    function if possible.
    """
    s = SFASchur(QQ)
    r = qsum
    sam = 0
    for mu in Partitions(n):
        s_mu = s(mu)
        new_r = r - s_to_q(mu)
        while is_positive(new_r):
            sam += s_mu
            r = new_r
            new_r = r - s_to_q(mu)
        if r == 0:
            break

    if r == 0:
        return sam
    else:
        raise ValueError, "cannot convert q to a Schur function"
def s_to_q(part):
    """
    Convert a Schur symmetric function to a Gessel quasisymmetric
    function.
    """
    res = 0
    for st in StandardTableaux(part):
        #Construct the permutation by reading the entries of
        #the standard tableau
        perm = Permutation( st.to_word_by_reading_order() )

        #Get the idescents of the permutation
        idescents = perm.idescents()

        #Add one to each of them since Python's indexing starts
        #at 0 and not one.
        idescents = map(lambda i: i+1, idescents)

        #Add the corresponding quasisymmetric function
        res += Q( idescents )
    return res


#############
#Sami Graphs#
#############
class SamiGraph(Graph):
    """
    A subclass of graph that provides extra functionality
    relating to dual equivalence graphs.
    """
    def __init__(self, *args, **kwargs):
        kwargs['multiedges'] = True
        Graph.__init__(self, *args, **kwargs)

    def i_edges(self, i):
        """
        Return a list of all edges labeled by i.
        """
        return self.edges(labels=i)

    def i_subgraph(self, i):
        """
        Returns the subgraph of self whose edges are labeled
        by i.
        """
        return SamiGraph(self.subgraph(edges=self.i_edges(i)))

    def show_connected_components(self, *args, **kwargs):
        """
        Show each of the connected components of self seperately.
        """
        for cc in self.connected_components():
            self.subgraph(vertices=cc).plot(*args, **kwargs).show()

    def subgraph(self, *args, **kwargs):
        """
        Provides a wrapper around Graph.subgraph which returns a SamiGraph
        instead of a Graph.
        """
        return SamiGraph( Graph.subgraph(self, *args, **kwargs) )

    def satisfies_axiom_1(self):
        """
        Check whether or not self satisfies axiom 1.
        """
        for w in self.vertices():
            signature = w.idescents_signature()
            #i is the index
            for i in range(1, len(w)-1):
                i_neighbors = self.neighbors(w, labels=i+1)
                if signature[i-1] == - signature[i] and (len(i_neighbors)!=1):
                    return False
        return True

    def satisfies_axiom_2(self):
        """
        Check whether or not self satisfies axiom 1.
        """
        for (w,x,i) in self.edges():
            sigma_w = w.idescents_signature()
            sigma_x = x.idescents_signature()

            for j in [i-1, i-2]:
                if sigma_w[j] != -sigma_x[j]:
                    return False
            for h in range(0,len(w)):
                if not ( h < i-2 or h > i+1 ):
                    continue
                if sigma_w[h-1] != sigma_x[h-1]:
                    return False
        return True

    def satisfies_axiom_3(self):
        """
        Check whether or not self satisfies axiom 3.
        """
        for (w,x,i) in self.edges():
            sigma_w = w.idescents_signature()
            sigma_x = x.idescents_signature()
            if sigma_w[i-3] == -sigma_x[i-3] and sigma_w[i-3] != -sigma_w[i-2]:
                return False
            if sigma_w[i] == -sigma_x[i] and sigma_w[i] != -sigma_w[i-1]:
                return False
        return True


    def satisfies_axiom_4(self):
        """
        Check whether or not self satisfies axiom 4.
        """
        w = self.vertices()[0]
        for i in range(4, len(w)):
            #Check the first part
            subgraph = self.i_subgraph([i,i-1,i-2])
            for cc in subgraph.connected_components():
                try:
                    a = q_to_s( subgraph.subgraph(vertices = cc).to_quasisymmetric_function(), len(w) )
                except ValueError:
                    return False

                if len(a) != 1:
                    return False

            #Check the second part
            if i < 5:
                continue
            subgraph = self.i_subgraph(range(2,i+1))
            for cc in subgraph.connected_components():
                try:
                    a = q_to_s( subgraph.subgraph(vertices = cc).to_quasisymmetric_function(), len(w) )
                except ValueError:
                    return False

                if len(a) != 1:
                    return False

        return True

    def satisfies_axiom_4_weak(self):
        """
        Check whether or not self satisfies the weak version of axiom 4.
        """
        w = self.vertices()[0]
        for i in range(3, len(w)-1):
            subgraph = self.i_subgraph([i,i-1])
            for cc in subgraph.connected_components():
               #Make sure that self's generating function on cc is Schur positive
                try:
                    a = q_to_s( subgraph.subgraph(vertices = cc).to_quasisymmetric_function(), len(w) )
                except ValueError:
                    return False

        for i in range(4, len(w)-1):
            subgraph = self.i_subgraph([i,i-1,i-2])
            for cc in subgraph.connected_components():
                #Make sure that self's generating function on cc is Schur positive
                try:
                    a = q_to_s( subgraph.subgraph(vertices = cc).to_quasisymmetric_function(), len(w) )
                except ValueError:
                    return False

        return True


    def satisfies_axiom_5(self):
        w = self.vertices()[0]
        for i in range(2, len(w)):
            for j in range(2, i-3) + range(i+3, len(w)):
                for x in self.vertices():
                    xi_neighbors = self.neighbors(x, labels=i)
                    xj_neighbors = self.neighbors(x, labels=j)
                    for w in xj_neighbors:
                        for y in xi_neighbors:
                            if len(Set(self.neighbors(w,labels=j)).intersection(Set(self.neighbors(y,labels=i)))) == 0:
                                return False
        return True

    def neighbors(self, v, labels=None):
        """
        Returns the neighbors of the vertex v in self.  If labels is not None, then
        it returns the neighbors of v which are connected by a label in labels.
        """
        if labels is None:
            return Graph.neighbors(self, v)
        else:
            if not isinstance(labels, list):
                labels = [labels]
            return [ e[1] for e in self.edges_incident(v) if e[2] in labels ]


    def to_quasisymmetric_function(self):
        """
        Returns the generating function of self.
        """
        res = 0
        for w in self.vertices():
            res += Q( map(lambda i: i+1, w.idescents()) )
        return res

    def print_axioms(self):
        print self.satisfies_axiom_1()
        print self.satisfies_axiom_2()
        print self.satisfies_axiom_3()
        print self.satisfies_axiom_4()
        print self.satisfies_axiom_4_weak()
        print self.satisfies_axiom_5()

    def is_dual_equivalence_graph(self):
        """
        Returns True if and only if self is a dual equivalence graph.
        """
        return self.satisfies_axiom_1() and \
                self.satisfies_axiom_2() and \
                self.satisfies_axiom_3() and \
                self.satisfies_axiom_4() and \
                self.satisfies_axiom_5()
    is_deg = is_dual_equivalence_graph

    def is_dgraph(self):
        """
        Returns True if and only if self is a D dgraph.
        """
        return self.satisfies_axiom_1() and \
                self.satisfies_axiom_2() and \
                self.satisfies_axiom_3() and \
                self.satisfies_axiom_4_weak() and \
                self.satisfies_axiom_5()

    is_dg = is_dgraph

    def has_edge(self, u, v, label=None):
        """
        A wrapper around Graph.has_edge since it has a bug which
        ignores the label.
        """
        return self._nxg.has_edge(u,v,label)

    def apply_transformation_1(self):
        """
        Applies transformation 1 to self.
        """
        w = self.vertices()[0]
        for i in range(3, len(w)-1):
            subgraph = self.i_subgraph(i)

        raise NotImplementedError

    def apply_transformation_2(self):
        """
        Applies transformation 2 to self.
        """
        raise NotImplementedError

    def apply_transformation_3(self):
        """
        """
        raise NotImplementedError

#####################
#Adjacency functions#
#####################
def imove(word, content, i, k):
    """
    Return either an i-switch or i-shift of word based on k and
    the content.

    If no move applies, then return word itself.
    """
    state, pos_im1, pos_i, pos_ip1 = Permutation(word)._icondition(i)
    if state is None:
        return word


    positions = [pos_im1, pos_i, pos_ip1]
    positions.sort()
    l = positions[0]
    r = positions[-1]

    if content[r]-content[l] > k:
        return word.iswitch(i)
    else:
        return word.ishift(i)

def macdonald_adjacency(shape):
    """
    Given a partition shape, return the adjacency function
    for the D-graph of H_shape.
    """
    k = shape[0]
    content = []
    start = 1
    for row_length in reversed(shape):
        content += range(start, start+row_length)
        start += k

    def f(perm):
        edges = []
        for i in range(2, len(perm)):
            im = imove(perm, content, i, k)
            if im is not perm:
                edges.append((im,i))
        return edges
    return f

def iswitches(perm):
    """
    The adjacency function for i-switches.
    """
    edges = []
    for i in range(2, len(perm)):
        i_switch = perm.iswitch(i)
        if i_switch is not perm:
            edges.append((i_switch, i))
    return edges

def ishifts(perm):
    """
    The adjacency function for i-shifts.
    """
    edges = []
    for i in range(2, len(perm)):
        i_shift = perm.ishift(i)
        if i_shift is not perm:
            edges.append((i_shift, i))
    return edges


##########################
#Graph creation functions#
##########################
def create_graph(words, f):
    """
    Creates a graph whose vertices are words and the edges intersecting
    vertex w are f(w).
    """
    g = SamiGraph()

    #Add the vertices and edges
    for w in words:
        try:
            hash(w)
        except TypeError:
            w = tuple(w)

        g.add_vertex(w)
        for (w2, color) in f(w):
            try:
                hash(w2)
            except TypeError:
                w2 = tuple(w2)

            if not g.has_edge(w, w2, label=color):
                g.add_edge(w, w2, color)
    return g


def permutation_graph(perm, f):
    """
    Given a permutation perm, return the connected component containing
    perm in the graph defined by the adjacency function f.
    """
    perm = Permutation(perm)
    g = SamiGraph()
    g.add_vertex(perm)
    todo = [ perm ]
    while todo != []:
        y = todo[0]
        del todo[0]

        for z,i, in f(y):
            if not g.has_vertex(z):
                g.add_vertex(z)
                todo.append(z)
            if not g.has_edge(y,z,i):
                g.add_edge(y,z,i)

    return g

def standard_deg(shape):
    perms = [ Permutation(st.to_word_by_reading_order()) for st in StandardTableaux(shape) ]
    return create_graph(perms, iswitches)

def macdonald_dgraph(shape):
    n = sum(shape)
    return create_graph(Permutations(n), macdonald_adjacency(shape))

#########################
#Miscellaneous Functions#
#########################
def tableau_signature(t):
    """
    Returnt the idescent signature of the word obtained from
    the reading order of t.
    """
    return Tableau(t).to_permutation_by_reading_order().idescents_signature()

########
#Checks#
########
def check_standard_tableaux(n):
    for p in Partitions(n):
        sd = standard_deg(p)
        assert sd.is_deg() is True
        assert sd.satisfies_axiom_4_weak() is True


def check_macdonald_polynomials(n):
    for mu in Partitions(n):
        f = macdonald_adjacency(mu)
        dgraph = create_graph(Permutations(n), f)
        assert dgraph.satisfies_axiom_1() is True
        assert dgraph.satisfies_axiom_2() is True
        assert dgraph.satisfies_axiom_3() is True
        assert dgraph.satisfies_axiom_5() is True

        assert dgraph.satisfies_axiom_4_weak() is True
        #assert dgraph.satisfies_axiom_4() is False


        #Verify that the major index is constant
        #on the connected components of the dgraph
        #of H_\mu
        for cc in dgraph.connected_components():
            value = cc[0].to_tableau_by_shape(mu).major_index()
            assert all([ value == w.to_tableau_by_shape(mu).major_index() for w in cc ])


        #Verify that the inversion number is constant
        #on the connected components of the dgraph of
        #H_\mu
        for cc in dgraph.connected_components():
            value = cc[0].to_tableau_by_shape(mu).inversion_number()
            print ([ w.to_tableau_by_shape(mu).inversion_number() for w in cc ])
