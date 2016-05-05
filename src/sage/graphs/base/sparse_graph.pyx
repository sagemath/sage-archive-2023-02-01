r"""
Fast sparse graphs

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Usage Introduction
------------------

::

    sage: from sage.graphs.base.sparse_graph import SparseGraph

Sparse graphs are initialized as follows::

    sage: S = SparseGraph(nverts = 10, expected_degree = 3, extra_vertices = 10)

This example initializes a sparse graph with room for twenty vertices, the first
ten of which are in the graph. In general, the first ``nverts`` are "active."
For example, see that 9 is already in the graph::

    sage: S.verts()
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    sage: S.add_vertex(9)
    9
    sage: S.verts()
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

But 10 is not, until we add it::

    sage: S.add_vertex(10)
    10
    sage: S.verts()
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

You can begin working with unlabeled arcs right away as follows::

    sage: S.add_arc(0,1)
    sage: S.add_arc(1,2)
    sage: S.add_arc(1,0)
    sage: S.has_arc(7,3)
    False
    sage: S.has_arc(0,1)
    True
    sage: S.in_neighbors(1)
    [0]
    sage: S.out_neighbors(1)
    [0, 2]
    sage: S.del_all_arcs(0,1)
    sage: S.all_arcs(0,1)
    []
    sage: S.all_arcs(1,2)
    [0]
    sage: S.del_vertex(7)
    sage: S.all_arcs(7,3)
    Traceback (most recent call last):
    ...
    LookupError: Vertex (7) is not a vertex of the graph.

Sparse graphs support multiple edges and labeled edges, but requires that the
labels be positive integers (the case label = 0 is treated as no label).

::

    sage: S.add_arc_label(0,1,-1)
    Traceback (most recent call last):
    ...
    ValueError: Label (-1) must be a nonnegative integer.
    sage: S.add_arc(0,1)
    sage: S.arc_label(0,1)
    0

Note that ``arc_label`` only returns the first edge label found in the specified
place, and this can be in any order (if you want all arc labels, use
``all_arcs``)::

    sage: S.add_arc_label(0,1,1)
    sage: S.arc_label(0,1)
    1
    sage: S.all_arcs(0,1)
    [0, 1]

Zero specifies only that there is no labeled arc::

    sage: S.arc_label(1,2)
    0

So do not be fooled::

    sage: S.all_arcs(1,2)
    [0]
    sage: S.add_arc(1,2)
    sage: S.arc_label(1,2)
    0

Instead, if you work with unlabeled edges, be sure to use the right functions::

    sage: T = SparseGraph(nverts = 3, expected_degree = 2)
    sage: T.add_arc(0,1)
    sage: T.add_arc(1,2)
    sage: T.add_arc(2,0)
    sage: T.has_arc(0,1)
    True

Sparse graphs are by their nature directed. As of this writing, you need to do
operations in pairs to treat the undirected case (or use a backend or a Sage
graph)::

    sage: T.has_arc(1,0)
    False

Multiple unlabeled edges are also possible::

    sage: for _ in range(10): S.add_arc(5,4)
    sage: S.all_arcs(5,4)
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

The curious developer is encouraged to check out the ``unsafe`` functions,
which do not check input but which run in pure C.

Underlying Data Structure
-------------------------

The class ``SparseGraph`` contains the following variables which are inherited
from ``CGraph`` (for explanation, refer to the documentation there)::

        cdef int num_verts
        cdef int num_arcs
        cdef int *in_degrees
        cdef int *out_degrees
        cdef bitset_t active_vertices

It also contains the following variables::

        cdef int hash_length
        cdef int hash_mask
        cdef SparseGraphBTNode **vertices

For each vertex ``u``, a hash table of length ``hash_length`` is instantiated.
An arc ``(u, v)`` is stored at ``u * hash_length + hash(v)`` of the array
``vertices``, where ``hash`` should be thought of as an arbitrary but fixed hash
function which takes values in ``0 <= hash < hash_length``. Each address may
represent different arcs, say ``(u, v1)`` and ``(u, v2)`` where
``hash(v1) == hash(v2)``. Thus, a binary tree structure is used at this step to
speed access to individual arcs, whose nodes (each of which represents a pair
``(u,v)``) are instances of the following type::

    cdef struct SparseGraphBTNode:
        int vertex
        int number
        SparseGraphLLNode *labels
        SparseGraphBTNode *left
        SparseGraphBTNode *right

Which range of the ``vertices`` array the root of the tree is in determines
``u``, and ``vertex`` stores ``v``. The integer ``number`` stores only the
number of unlabeled arcs from ``u`` to ``v``.

Currently, labels are stored in a simple linked list, whose nodes are instances
of the following type::

    cdef struct SparseGraphLLNode:
        int label
        int number
        SparseGraphLLNode *next

The int ``label`` must be a positive integer, since 0 indicates no label, and
negative numbers indicate errors. The int ``number`` is the number of arcs with
the given label.

TODO: Optimally, edge labels would also be represented by a binary tree, which
would help performance in graphs with many overlapping edges. Also, a more
efficient binary tree structure could be used, although in practice the trees
involved will usually have very small order, unless the degree of vertices
becomes significantly larger than the ``expected_degree`` given, because this is
the size of each hash table. Indeed, the expected size of the binary trees is
`\frac{\text{actual degree}}{\text{expected degree}}`. Ryan Dingman, e.g., is
working on a general-purpose Cython-based red black tree, which would be optimal
for both of these uses.
"""

#*****************************************************************************
#       Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from libc.string cimport memset
include 'sage/data_structures/bitset.pxi'

cdef enum:
    BT_REORDERING_CONSTANT = 145533211
    # Since the binary tree will often see vertices coming in already sorted,
    # we don't use the normal ordering on integers, instead multiplying by a
    # randomly chosen number and (after reducing mod the size of integers)
    # comparing the result. This isn't necessarily the most efficient way to do
    # things, but it may just be on binary trees that are never bigger than two
    # or three nodes.

cdef inline int compare(int a, int b):
    # Here we rely on the fact that C performs arithmetic on unsigned
    # ints modulo 2^wordsize.
    cdef unsigned int aa = a, bb = b # signed ints lead to badness like a>b>c>a...
    if aa*BT_REORDERING_CONSTANT > bb*BT_REORDERING_CONSTANT:
        return 1
    elif aa*BT_REORDERING_CONSTANT < bb*BT_REORDERING_CONSTANT:
        return -1
    return 0

cdef class SparseGraph(CGraph):
    """
    Compiled sparse graphs.

    ::

        sage: from sage.graphs.base.sparse_graph import SparseGraph

    Sparse graphs are initialized as follows::

        sage: S = SparseGraph(nverts = 10, expected_degree = 3, extra_vertices = 10)

    INPUT:

     - ``nverts`` - non-negative integer, the number of vertices.
     - ``expected_degree`` - non-negative integer (default: 16), expected upper
        bound on degree of vertices.
     - ``extra_vertices`` - non-negative integer (default: 0), how many extra
        vertices to allocate.
     - ``verts`` - optional list of vertices to add
     - ``arcs`` - optional list of arcs to add

    The first ``nverts`` are created as vertices of the graph, and the next
    ``extra_vertices`` can be freely added without reallocation. See top level
    documentation for more details. The input ``verts`` and ``arcs`` are mainly
    for use in pickling.

    """

    def __cinit__(self, int nverts, int expected_degree = 16, int extra_vertices = 10, verts=None, arcs=None):
        """
        Allocation and initialization happen in one place.

        Memory usage is roughly

        O(  (nverts + extra_vertices)*expected_degree + num_arcs  ).

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts = 10, expected_degree = 3, extra_vertices = 10)

        TESTS::

            sage: Graph(-1)
            Traceback (most recent call last):
            ...
            ValueError: The number of vertices cannot be strictly negative!
        """
        cdef int i = 1
        if nverts < 0:
            raise ValueError("The number of vertices cannot be strictly negative!")
        if nverts == 0 and extra_vertices == 0:
            raise RuntimeError('Sparse graphs must allocate space for vertices!')
        self.num_verts = nverts
        nverts += extra_vertices
        self.num_arcs = 0
        while i < expected_degree:
            i = i << 1
        self.hash_length = i
        self.hash_mask = i - 1

        # Allocating memory
        self.vertices = <SparseGraphBTNode **> \
          sig_malloc(nverts * self.hash_length * sizeof(SparseGraphBTNode *))
        self.in_degrees = <int *> sig_malloc(nverts * sizeof(int))
        self.out_degrees = <int *> sig_malloc(nverts * sizeof(int))

        # Checking the memory was actually allocated
        if not self.vertices or not self.in_degrees or not self.out_degrees:
            if self.vertices: sig_free(self.vertices)
            if self.in_degrees: sig_free(self.in_degrees)
            if self.out_degrees: sig_free(self.out_degrees)
            raise RuntimeError("Failure allocating memory.")

        # Initializing variables:
        #
        # self.vertices[i] = 0
        memset(self.vertices, <int> NULL, nverts * self.hash_length * sizeof(SparseGraphBTNode *))

        # self.in_degrees[i] = 0
        memset(self.in_degrees, 0, nverts * sizeof(int))

        # self.out_degrees[i] = 0
        memset(self.out_degrees, 0, nverts * sizeof(int))

        bitset_init(self.active_vertices, self.num_verts + extra_vertices)
        bitset_set_first_n(self.active_vertices, self.num_verts)

        if verts is not None:
            self.add_vertices(verts)

        if arcs is not None:
            for u,v,l in arcs:
                self.add_arc_label(u,v,l)

    def __dealloc__(self):
        """
        New and dealloc are both tested at class level.
        """
        cdef SparseGraphBTNode **temp
        cdef SparseGraphLLNode *label_temp
        cdef int i

        # Freeing the list of arcs attached to each vertex
        for i from 0 <= i < self.active_vertices.size * self.hash_length:
            temp = &(self.vertices[i])

            # While temp[0]=self.vertices[i] is not NULL, find a leaf in the
            # tree rooted at temp[0] and free it. Then go back to temp[0] and do
            # it again. When self.vertices[i] is NULL, go for self.vertices[i+1]
            while temp[0] != NULL:
                if temp[0].left != NULL:
                    temp = &(temp[0].left)
                elif temp[0].right != NULL:
                    temp = &(temp[0].right)
                else:
                    label_temp = temp[0].labels
                    while label_temp != NULL:
                        temp[0].labels = label_temp.next
                        sig_free(label_temp)
                        label_temp = temp[0].labels
                    sig_free(temp[0])
                    temp[0] = NULL
                    temp = &(self.vertices[i])

        sig_free(self.vertices)
        sig_free(self.in_degrees)
        sig_free(self.out_degrees)
        bitset_free(self.active_vertices)

    cpdef realloc(self, int total):
        """
        Reallocate the number of vertices to use, without actually adding any.

        INPUT:

         - ``total`` - integer, the total size to make the array

        Returns -1 and fails if reallocation would destroy any active vertices.

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.current_allocation()
            8
            sage: S.add_vertex(6)
            6
            sage: S.current_allocation()
            8
            sage: S.add_vertex(10)
            10
            sage: S.current_allocation()
            16
            sage: S.add_vertex(40)
            Traceback (most recent call last):
            ...
            RuntimeError: Requested vertex is past twice the allocated range: use realloc.
            sage: S.realloc(50)
            sage: S.add_vertex(40)
            40
            sage: S.current_allocation()
            50
            sage: S.realloc(30)
            -1
            sage: S.current_allocation()
            50
            sage: S.del_vertex(40)
            sage: S.realloc(30)
            sage: S.current_allocation()
            30

        """
        if total == 0:
            raise RuntimeError('Sparse graphs must allocate space for vertices!')
        cdef bitset_t bits
        if total < self.active_vertices.size:
            bitset_init(bits, self.active_vertices.size)
            bitset_set_first_n(bits, total)
            if not bitset_issubset(self.active_vertices, bits):
                bitset_free(bits)
                return -1
            bitset_free(bits)

        self.vertices = <SparseGraphBTNode **> sig_realloc(self.vertices, total * self.hash_length * sizeof(SparseGraphBTNode *))
        self.in_degrees = <int *> sig_realloc(self.in_degrees, total * sizeof(int))
        self.out_degrees = <int *> sig_realloc(self.out_degrees, total * sizeof(int))

        cdef int new_vertices = total - self.active_vertices.size

        # Initializing the entries corresponding to new vertices if any
        if new_vertices>0:

            # self.vertices
            memset(self.vertices+self.active_vertices.size *  self.hash_length,
                   <int> NULL,
                   new_vertices * self.hash_length * sizeof(SparseGraphBTNode *))

            # self.int_degrees
            memset(self.in_degrees+self.active_vertices.size, 0, new_vertices * sizeof(int))

            # self.out_degrees
            memset(self.out_degrees+self.active_vertices.size, 0, new_vertices * sizeof(int))

        # self.active_vertices
        bitset_realloc(self.active_vertices, total)

    ###################################
    # Unlabeled arc functions
    ###################################

    cdef int add_arc_unsafe(self, int u, int v) except -1:
        """
        Adds arc (u, v) to the graph with no label.

        INPUT:
            u, v -- non-negative integers
        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared
        cdef SparseGraphBTNode **ins_pt = &(self.vertices[i])
        while ins_pt[0] != NULL:
            compared = compare(ins_pt[0].vertex, v)
            if compared > 0:
                ins_pt = &(ins_pt[0].left)
            elif compared < 0:
                ins_pt = &(ins_pt[0].right)
            else:
                ins_pt[0].number += 1
                break
        if ins_pt[0] == NULL:
            ins_pt[0] = <SparseGraphBTNode *> sig_malloc(sizeof(SparseGraphBTNode))
            if not ins_pt[0]:
                raise RuntimeError("Failure allocating memory.")
            ins_pt[0].vertex = v
            ins_pt[0].number = 1
            ins_pt[0].left = NULL
            ins_pt[0].right = NULL
            ins_pt[0].labels = NULL
        self.in_degrees[v] += 1
        self.out_degrees[u] += 1
        self.num_arcs += 1

    cpdef add_arc(self, int u, int v):
        """
        Adds arc ``(u, v)`` to the graph with no label.

        INPUT:

         - ``u, v`` -- non-negative integers, must be in self

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(4,7)
            Traceback (most recent call last):
            ...
            LookupError: Vertex (7) is not a vertex of the graph.
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        self.check_vertex(u)
        self.check_vertex(v)
        self.add_arc_unsafe(u,v)

    cdef int has_arc_unsafe(self, int u, int v) except -1:
        """
        Checks whether arc (u, v) is in the graph.

        INPUT:
            u, v -- non-negative integers, must be in self

        OUTPUT:
            0 -- False
            1 -- True

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef SparseGraphBTNode *temp = self.vertices[i]
        while temp != NULL:
            if temp.vertex == v:
                return 1
            if compare(temp.vertex, v) > 0:
                temp = temp.left
            else: # note compare < 0
                temp = temp.right
        return 0

    cpdef bint has_arc(self, int u, int v) except -1:
        """
        Checks whether arc ``(u, v)`` is in the graph.

        INPUT:
         - ``u, v`` - integers

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1)
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        if u < 0 or u >= self.active_vertices.size or not bitset_in(self.active_vertices, u):
            return False
        if v < 0 or v >= self.active_vertices.size or not bitset_in(self.active_vertices, v):
            return False
        return self.has_arc_unsafe(u,v)

    cdef int del_arc_unsafe(self, int u, int v) except -1:
        """
        Deletes *all* arcs from u to v.

        INPUT:
            u, v -- non-negative integers, must be in self

        OUTPUT:
            0 -- No error.
            1 -- No arc to delete.

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared, left_len, right_len
        cdef SparseGraphBTNode *temp
        cdef SparseGraphBTNode **left_child
        cdef SparseGraphBTNode **right_child
        cdef SparseGraphBTNode **parent = &self.vertices[i]
        cdef SparseGraphLLNode *labels

        # Assigning to parent the SparseGraphBTNode corresponding to arc (u,v)
        while parent[0] != NULL:
            compared = compare(parent[0].vertex, v)
            if compared > 0:
                parent = &(parent[0].left)
            elif compared < 0:
                parent = &(parent[0].right)
            else:# if parent[0].vertex == v:
                break

        # If not found, there is no arc to delete !
        if parent[0] == NULL:
            return 1

        # now parent[0] points to the BT node corresponding to (u,v)
        labels = parent[0].labels
        i = parent[0].number
        self.in_degrees[v] -= i
        self.out_degrees[u] -= i
        self.num_arcs -= i

        # Freeing the labels
        while labels != NULL:
            i = labels.number
            parent[0].labels = parent[0].labels.next
            sig_free(labels)
            labels = parent[0].labels
            self.in_degrees[v] -= i
            self.out_degrees[u] -= i
            self.num_arcs -= i

        # Now, if the SparseGraphBTNode element is to be removed, it has to be
        # replaced in the binary tree by one of its children.

        # If there is no left child
        if parent[0].left == NULL:
            temp = parent[0]
            parent[0] = parent[0].right
            sig_free(temp)
            return 0

        # If there is no right child
        elif parent[0].right == NULL:
            temp = parent[0]
            parent[0] = parent[0].left
            sig_free(temp)
            return 0

        # Both children
        else:
            left_len = 0
            right_len = 0
            left_child = &(parent[0].left)
            right_child = &(parent[0].right)

            # left_len is equal to the maximum length of a path LR...R. The
            # last element of this path is the value of left_child

            while left_child[0].right != NULL:
                left_len += 1
                left_child = &(left_child[0].right)
            # right_len is equal to the maximum length of a path RL...L. The
            # last element of this path is the value of right_child

            while right_child[0].left != NULL:
                right_len += 1
                right_child = &(right_child[0].left)

            # According to the respective lengths, replace parent by the left or
            # right child and place the other child at its expected place.
            if left_len > right_len:
                left_child[0].right = parent[0].right
                temp = parent[0]
                parent[0] = left_child[0]
                left_child[0] = left_child[0].left
                parent[0].left = temp.left
                sig_free(temp)
                return 0
            else:
                right_child[0].left = parent[0].left
                temp = parent[0]
                parent[0] = right_child[0]
                right_child[0] = right_child[0].right
                parent[0].right = temp.right
                sig_free(temp)
                return 0

    cpdef del_all_arcs(self, int u, int v):
        """
        Deletes all arcs from ``u`` to ``v``.

        INPUT:
         - ``u, v`` - integers

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1,0)
            sage: G.add_arc_label(0,1,1)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,3)
            sage: G.del_all_arcs(0,1)
            sage: G.has_arc(0,1)
            False
            sage: G.arc_label(0,1)
            0
            sage: G.del_all_arcs(0,1)

        """
        self.check_vertex(u)
        self.check_vertex(v)
        self.del_arc_unsafe(u,v)

    ###################################
    # Neighbor functions
    ###################################

    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except -2:
        """
        Gives all v such that (u, v) is an arc of the graph.

        INPUT:

        - ``u`` -- non-negative integer, must be in self neighbors -- must be a
            pointer to an (allocated) integer array size -- the length of the
            array

        OUTPUT:

            nonnegative integer -- the number of v such that (u, v) is an arc -1
            -- indicates that the array has been filled with neighbors, but
            there were more

        """
        cdef int i, num_nbrs = 0, current_nbr = 0
        if self.out_degrees[u] == 0:
            return 0

        cdef SparseGraphBTNode ** pointers[1]
        cdef list l = []
        cdef int n_neighbors = self.out_neighbors_BTNode_unsafe(u, pointers)
        if size >= n_neighbors:
            for i in range(n_neighbors):
                neighbors[i] = pointers[0][i].vertex
        else:
            for i in range(size):
                neighbors[i] = pointers[0][i].vertex
            n_neighbors = -1

        sig_free(pointers[0])
        return n_neighbors

    cdef int out_neighbors_BTNode_unsafe(self, int u, SparseGraphBTNode *** p_pointers):
        """
        Lists the out-neighbors of a vertex as BTNodes

        Technically, this function transforms a binary tree into a list. The
        information it returns is a list of pointers toward a
        ``SparseGraphBTNode``, thus a ``SparseGraphBTNode **``.

        INPUT:

        - ``u`` -- the vertex to consider

        - ``p_pointers`` -- a pointer toward a ``SparseGraphBTNode **``, i.e. a
          ``SparseGraphBTNode ***``. When the function terminates,
          ``p_pointers[0]`` points toward a filled ``SparseGraphBTNode **``. It
          returns the length of this array.

        .. NOTE::

            Don't forget to free ``p_pointers[0]``  !
        """
        cdef int i, num_nbrs = 0, current_nbr = 0
        cdef int degree = self.out_degrees[u]
        if degree == 0:
            p_pointers[0] = NULL
            return 0
        cdef SparseGraphBTNode **pointers = <SparseGraphBTNode **> sig_malloc(degree * sizeof(SparseGraphBTNode *))
        p_pointers[0] = pointers
        if pointers == NULL:
            raise RuntimeError("Failure allocating memory.")
        for i from u * self.hash_length <= i < (u+1) * self.hash_length:
            if self.vertices[i] == NULL:
                continue
            pointers[num_nbrs] = self.vertices[i]
            num_nbrs += 1

            # While all the neighbors have not been added to the list, explore
            # element pointers[current_nbr] and append its children to the end
            # of pointers if necessary, the increment current_nbr.
            while current_nbr < num_nbrs:
                if pointers[current_nbr].left != NULL:
                    pointers[num_nbrs] = pointers[current_nbr].left
                    num_nbrs += 1
                if pointers[current_nbr].right != NULL:
                    pointers[num_nbrs] = pointers[current_nbr].right
                    num_nbrs += 1
                current_nbr += 1
        return num_nbrs

    cpdef list out_neighbors(self, int u):
        """
        Gives all ``v`` such that ``(u, v)`` is an arc of the graph.

        INPUT:
         - ``u`` - integer

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(1,2)
            sage: G.add_arc(1,3)
            sage: G.out_neighbors(0)
            [1]
            sage: G.out_neighbors(1)
            [2, 3]

        """
        cdef int i, num_nbrs
        self.check_vertex(u)
        if self.out_degrees[u] == 0:
            return []
        cdef int size = self.out_degrees[u]
        cdef int *neighbors = <int *> sig_malloc(size * sizeof(int))
        if not neighbors:
            raise RuntimeError("Failure allocating memory.")
        num_nbrs = self.out_neighbors_unsafe(u, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sig_free(neighbors)
        return output

    cpdef int out_degree(self, int u):
        """
        Returns the out-degree of ``v``

        INPUT:

         - ``u`` - integer

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(1,2)
            sage: G.add_arc(1,3)
            sage: G.out_degree(0)
            1
            sage: G.out_degree(1)
            2
        """
        return self.out_degrees[u]

    cdef list out_arcs_unsafe(self, int u, bint labels):
        r"""
        Builds the list of arcs leaving a vertex.

        Note that the source of each edge is *NOT* returned.

        INPUT:

        - ``u`` -- the vertex to consider

        - ``labels`` -- whether to return the labels alors with the outneighbor.
          If set to ``True``, the function returns a list of pairs
          ``(destination, label)`` for each arc leaving `u`. If set to
          ``False``, it returns a list of outneighbors (with multiplicity if
          several edges link two vertices).
        """
        cdef SparseGraphBTNode ** pointers[1]
        cdef SparseGraphBTNode * node
        cdef int neighbors = self.out_neighbors_BTNode_unsafe(u, pointers)
        cdef SparseGraphLLNode *label
        cdef int i,j
        cdef list l = []
        if labels:
            for i in range(neighbors):
                node = pointers[0][i]
                for j in range(node.number):
                    l.append((node.vertex, 0))
                label = node.labels
                while label != NULL:
                    for k in range(label.number):
                        l.append((node.vertex, label.label))
                    label = label.next
        else:
            for i in range(neighbors):
                node = pointers[0][i]
                for j in range(node.number):
                    l.append(node.vertex)
                label = node.labels
                while label != NULL:
                    for k in range(label.number):
                        l.append(node.vertex)
                    label = label.next

        if pointers[0] != NULL:
            sig_free(pointers[0])

        return l

    cdef int in_neighbors_unsafe(self, int v, int *neighbors, int size) except -2:
        """
        Gives all u such that (u, v) is an arc of the graph.

        INPUT:
            v -- non-negative integer, must be in self
            neighbors -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            nonnegative integer -- the number of u such that (u, v) is an arc
            -1 -- indicates that the array has been filled with neighbors, but
        there were more

        NOTE: Due to the implementation of SparseGraph, this method is much more
        expensive than out_neighbors_unsafe.

        """
        cdef int i, num_nbrs = 0
        if self.in_degrees[v] == 0:
            return 0
        for i from 0 <= i < self.active_vertices.size:
            if not bitset_in(self.active_vertices, i): continue
            if self.has_arc_unsafe(i, v):
                if num_nbrs == size:
                    return -1
                neighbors[num_nbrs] = i
                num_nbrs += 1
        return num_nbrs

    cpdef list in_neighbors(self, int v):
        """
        Gives all ``u`` such that ``(u, v)`` is an arc of the graph.

        INPUT:
         - ``v`` - integer

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(3,1)
            sage: G.add_arc(1,3)
            sage: G.in_neighbors(1)
            [0, 3]
            sage: G.in_neighbors(3)
            [1]

        NOTE: Due to the implementation of SparseGraph, this method is much more
        expensive than neighbors_unsafe.
        """
        cdef int i, num_nbrs
        self.check_vertex(v)
        if self.in_degrees[v] == 0:
            return []
        cdef int size = self.in_degrees[v]
        cdef int *neighbors = <int *> sig_malloc(size * sizeof(int))
        if not neighbors:
            raise RuntimeError("Failure allocating memory.")
        num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sig_free(neighbors)
        return output

    cpdef int in_degree(self, int u):
        """
        Returns the in-degree of ``v``

        INPUT:
         - ``u`` - integer

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(1,2)
            sage: G.add_arc(1,3)
            sage: G.in_degree(0)
            0
            sage: G.in_degree(1)
            1
        """
        return self.in_degrees[u]


    ###################################
    # Labeled arc functions
    ###################################

    cdef int add_arc_label_unsafe(self, int u, int v, int l) except -1:
        """
        Adds arc (u, v) to the graph with label l.

        INPUT:
            u, v -- non-negative integers
            l -- a positive integer label, or zero for no label

        OUTPUT:
            0 -- No error.

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared
        cdef SparseGraphBTNode **ins_pt = &(self.vertices[i])
        cdef SparseGraphLLNode *label_ptr
        while ins_pt[0] != NULL:
            compared = compare(ins_pt[0].vertex, v)
            if compared > 0:
                ins_pt = &(ins_pt[0].left)
            elif compared < 0:
                ins_pt = &(ins_pt[0].right)
            else:
                break
        if ins_pt[0] == NULL:
            ins_pt[0] = <SparseGraphBTNode *> sig_malloc(sizeof(SparseGraphBTNode))
            if not ins_pt[0]:
                raise RuntimeError("Failure allocating memory.")
            ins_pt[0].number = 0
            ins_pt[0].vertex = v
            ins_pt[0].left = NULL
            ins_pt[0].right = NULL
            ins_pt[0].labels = NULL
        if l:
            label_ptr = ins_pt[0].labels
            while label_ptr != NULL and label_ptr.label != l:
                label_ptr = label_ptr.next
            if label_ptr == NULL:
                label_ptr = <SparseGraphLLNode *> sig_malloc(sizeof(SparseGraphLLNode))
                if not label_ptr:
                    sig_free(ins_pt[0])
                    raise RuntimeError("Failure allocating memory.")
                label_ptr.label = l
                label_ptr.number = 1
                label_ptr.next = ins_pt[0].labels
                ins_pt[0].labels = label_ptr
            else:
                label_ptr.number += 1
        else:
            ins_pt[0].number += 1
        self.in_degrees[v] += 1
        self.out_degrees[u] += 1
        self.num_arcs += 1

    def add_arc_label(self, int u, int v, int l=0):
        """
        Adds arc ``(u, v)`` to the graph with label ``l``.

        INPUT:
         - ``u, v`` - non-negative integers, must be in self
         - ``l`` - a positive integer label, or zero for no label

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1)
            sage: G.add_arc_label(4,7)
            Traceback (most recent call last):
            ...
            LookupError: Vertex (7) is not a vertex of the graph.
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True
            sage: G.add_arc_label(1,2,2)
            sage: G.arc_label(1,2)
            2

        """
        self.check_vertex(u)
        self.check_vertex(v)
        if l < 0:
            raise ValueError("Label ({0}) must be a nonnegative integer.".format(l))
        self.add_arc_label_unsafe(u,v,l)

    cdef int arc_label_unsafe(self, int u, int v):
        """
        Retrieves the first label found associated with (u, v) (a positive
        integer).

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        OUTPUT:
            positive integer -- indicates that there is a label on (u, v).
            0 -- either the arc (u, v) is unlabeled, or there is no arc at all.

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared
        cdef SparseGraphBTNode *temp = self.vertices[i]
        while temp != NULL:
            compared = compare(temp.vertex, v)
            if compared > 0:
                temp = temp.left
            elif compared < 0:
                temp = temp.right
            else:
                break
        if temp == NULL or temp.labels == NULL:
            return 0
        return temp.labels.label

    cpdef int arc_label(self, int u, int v):
        """
        Retrieves the first label found associated with ``(u, v)``.

        INPUT:
         - ``u, v`` - non-negative integers, must be in self

        OUTPUT:
         - positive integer - indicates that there is a label on ``(u, v)``.
         - 0 - either the arc ``(u, v)`` is unlabeled, or there is no arc at all.

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(3,4,7)
            sage: G.arc_label(3,4)
            7

        NOTES:

        To this function, an unlabeled arc is indistinguishable from a non-arc::

            sage: G.add_arc_label(1,0)
            sage: G.arc_label(1,0)
            0
            sage: G.arc_label(1,1)
            0

        This function only returns the *first* label it finds from ``u`` to ``v``::

            sage: G.add_arc_label(1,2,1)
            sage: G.add_arc_label(1,2,2)
            sage: G.arc_label(1,2)
            2

        """
        self.check_vertex(u)
        self.check_vertex(v)
        return self.arc_label_unsafe(u,v)

    cdef int all_arcs_unsafe(self, int u, int v, int *arc_labels, int size):
        """
        Gives the labels of all arcs (u, v).

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            arc_labels -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            integer -- the number of arcs (u, v)
            -1 -- indicates that the array has been filled with labels, but
        there were more

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask), j
        cdef int compared, num_arcs
        cdef SparseGraphBTNode *temp = self.vertices[i]
        cdef SparseGraphLLNode *label
        while temp != NULL:
            compared = compare(temp.vertex, v)
            if compared > 0:
                temp = temp.left
            elif compared < 0:
                temp = temp.right
            else: # temp.vertex == v:
                break
        if temp == NULL:
            return 0
        j = 0
        num_arcs = temp.number
        while j < num_arcs and j < size:
            arc_labels[j] = 0
            j += 1
        label = temp.labels
        while label != NULL:
            num_arcs += label.number
            while j < num_arcs and j < size:
                arc_labels[j] = label.label
                j += 1
            label = label.next
        if j == size and label != NULL:
            return -1
        return num_arcs

    cpdef list all_arcs(self, int u, int v):
        """
        Gives the labels of all arcs ``(u, v)``. An unlabeled arc is interpreted as
        having label 0.

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(1,2,1)
            sage: G.add_arc_label(1,2,2)
            sage: G.add_arc_label(1,2,2)
            sage: G.add_arc_label(1,2,2)
            sage: G.add_arc_label(1,2,3)
            sage: G.add_arc_label(1,2,3)
            sage: G.add_arc_label(1,2,4)
            sage: G.all_arcs(1,2)
            [4, 3, 3, 2, 2, 2, 1]

        """
        cdef int size, num_arcs, i
        cdef int *arc_labels
        cdef list output
        self.check_vertex(u)
        self.check_vertex(v)
        if self.in_degrees[v] < self.out_degrees[u]:
            size = self.in_degrees[v]
        else:
            size = self.out_degrees[u]
        arc_labels = <int *> sig_malloc(size * sizeof(int))
        if not arc_labels:
            raise RuntimeError("Failure allocating memory.")
        num_arcs = self.all_arcs_unsafe(u, v, arc_labels, size)
        if num_arcs == -1:
            sig_free(arc_labels)
            raise RuntimeError("There was an error: there seem to be more arcs than self.in_degrees or self.out_degrees indicate.")
        output = [arc_labels[i] for i from 0 <= i < num_arcs]
        sig_free(arc_labels)
        return output

    cdef int del_arc_label_unsafe(self, int u, int v, int l):
        """
        Delete an arc (u, v) with label l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label, or zero for no label

        OUTPUT:
            0 -- No error.
            1 -- No arc with label l.

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared
        cdef SparseGraphBTNode **parent = &self.vertices[i]
        cdef SparseGraphLLNode **labels
        cdef SparseGraphLLNode *label
        while parent[0] != NULL:
            compared = compare(parent[0].vertex, v)
            if compared > 0:
                parent = &(parent[0].left)
            elif compared < 0:
                parent = &(parent[0].right)
            else: # if parent[0].vertex == v:
                break
        if parent[0] == NULL:
            return 1 # indicate an error
        if l == 0:
            if parent[0].number > 1: parent[0].number -= 1
            elif parent[0].number == 1:
                if parent[0].labels == NULL:
                    self.del_arc_unsafe(u, v)
                    return 0
                else: parent[0].number -= 1
            else: return 1 # indicate an error
        else:
            labels = &(parent[0].labels)
            while labels[0] != NULL and labels[0].label != l:
                labels = &(labels[0].next)
            if labels[0] == NULL:
                return 1
            label = labels[0]
            if label.number > 1:
                label.number -= 1
            else:
                labels[0] = labels[0].next
                sig_free(label)
                if labels == &(parent[0].labels) and labels[0] == NULL and parent[0].number == 0:
                    # here we need to delete an "empty" binary tree node
                    self.del_arc_unsafe(u, v)
        self.in_degrees[v] -= 1
        self.out_degrees[u] -= 1
        self.num_arcs -= 1

    cpdef del_arc_label(self, int u, int v, int l):
        """
        Delete an arc ``(u, v)`` with label ``l``.

        INPUT:
         - ``u, v`` - non-negative integers, must be in self
         - ``l`` - a positive integer label, or zero for no label

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1,0)
            sage: G.add_arc_label(0,1,1)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,3)
            sage: G.del_arc_label(0,1,2)
            sage: G.all_arcs(0,1)
            [0, 3, 2, 1]
            sage: G.del_arc_label(0,1,0)
            sage: G.all_arcs(0,1)
            [3, 2, 1]

        """
        self.check_vertex(u)
        self.check_vertex(v)
        if l < 0:
            raise ValueError("Label ({0}) must be a nonnegative integer.".format(l))
        self.del_arc_label_unsafe(u,v,l)

    cdef int has_arc_label_unsafe(self, int u, int v, int l):
        """
        Indicates whether there is an arc (u, v) with label l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label, or zero for no label

        OUTPUT:
            0 -- False
            1 -- True

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared
        cdef SparseGraphBTNode *temp = self.vertices[i]
        cdef SparseGraphLLNode *label
        while temp != NULL:
            compared = compare(temp.vertex, v)
            if compared > 0:
                temp = temp.left
            elif compared < 0:
                temp = temp.right
            else:# if temp.vertex == v:
                break
        if temp == NULL:
            return 0
        if l == 0 and temp.number > 0:
            return 1
        label = temp.labels
        while label != NULL:
            if label.label == l:
                return 1
            label = label.next
        return 0

    cpdef bint has_arc_label(self, int u, int v, int l):
        """
        Indicates whether there is an arc ``(u, v)`` with label ``l``.

        INPUT:
         - ``u, v`` -- non-negative integers, must be in self
         - ``l`` -- a positive integer label, or zero for no label

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1,0)
            sage: G.add_arc_label(0,1,1)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,2)
            sage: G.has_arc_label(0,1,1)
            True
            sage: G.has_arc_label(0,1,2)
            True
            sage: G.has_arc_label(0,1,3)
            False

        """
        self.check_vertex(u)
        self.check_vertex(v)
        if l < 0:
            raise ValueError("Label ({0}) must be a nonnegative integer.".format(l))
        return self.has_arc_label_unsafe(u,v,l) == 1

##############################
# Further tests. Unit tests for methods, functions, classes defined with cdef.
##############################

def _test_adjacency_sequence_out():
    """
    Randomly test the method ``SparseGraph.adjacency_sequence_out()``. No output
    indicates that no errors were found.

    TESTS::

        sage: from sage.graphs.base.sparse_graph import _test_adjacency_sequence_out
        sage: _test_adjacency_sequence_out()  # long time
    """
    from sage.graphs.digraph import DiGraph
    from sage.graphs.graph_generators import GraphGenerators
    from sage.misc.prandom import randint, random
    low = 0
    high = 1000
    randg = DiGraph(GraphGenerators().RandomGNP(randint(low, high), random()))
    n = randg.order()
    # set all labels to 0
    E = [(u, v, 0) for u, v in randg.edges(labels=False)]
    cdef SparseGraph g = SparseGraph(n,
                                     verts=randg.vertices(),
                                     arcs=E)
    assert g.num_verts == randg.order(), (
        "Graph order mismatch: %s vs. %s" % (g.num_verts, randg.order()))
    assert g.num_arcs == randg.size(), (
        "Graph size mismatch: %s vs. %s" % (g.num_arcs, randg.size()))
    M = randg.adjacency_matrix()
    cdef int *V = <int *>sig_malloc(n * sizeof(int))
    cdef int i = 0
    for v in randg.vertex_iterator():
        V[i] = v
        i += 1
    cdef int *seq = <int *> sig_malloc(n * sizeof(int))
    for 0 <= i < randint(50, 101):
        u = randint(low, n - 1)
        g.adjacency_sequence_out(n, V, u, seq)
        A = [seq[k] for k in range(n)]
        try:
            assert A == list(M[u])
        except AssertionError:
            sig_free(V)
            sig_free(seq)
            raise AssertionError("Graph adjacency mismatch")
    sig_free(seq)
    sig_free(V)

###########################################
# Sparse Graph Backend
###########################################

cdef class SparseGraphBackend(CGraphBackend):
    """
    Backend for Sage graphs using SparseGraphs.

    ::

        sage: from sage.graphs.base.sparse_graph import SparseGraphBackend

    This class is only intended for use by the Sage Graph and DiGraph class.
    If you are interested in using a SparseGraph, you probably want to do
    something like the following example, which creates a Sage Graph instance
    which wraps a SparseGraph object::

        sage: G = Graph(30, implementation="c_graph", sparse=True)
        sage: G.add_edges([(0,1), (0,3), (4,5), (9, 23)])
        sage: G.edges(labels=False)
        [(0, 1), (0, 3), (4, 5), (9, 23)]

    Note that Sage graphs using the backend are more flexible than SparseGraphs
    themselves. This is because SparseGraphs (by design) do not deal with Python
    objects::

        sage: G.add_vertex((0,1,2))
        sage: G.vertices()
        [0,
        ...
         29,
         (0, 1, 2)]
        sage: from sage.graphs.base.sparse_graph import SparseGraph
        sage: SG = SparseGraph(30)
        sage: SG.add_vertex((0,1,2))
        Traceback (most recent call last):
        ...
        TypeError: an integer is required

    """

    def __init__(self, int n, directed=True):
        """
        Initialize a sparse graph with n vertices.

        EXAMPLE:

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        """
        self._cg = SparseGraph(n)
        self._cg_rev = SparseGraph(n) if directed else self._cg
        self._directed = directed
        self.vertex_labels = {}
        self.vertex_ints = {}
        self.edge_labels = {}
        self.edge_labels_max = 1
        self.edge_labels_available_ids = []

    cdef inline int new_edge_label(self, object l):
        """
        Returns a new unique int representing the arbitrary label l.
        """
        if l is None:
            return 0

        cdef int l_int
        if self.edge_labels_available_ids:
            l_int = self.edge_labels_available_ids.pop(-1)
        else:
            l_int = self.edge_labels_max
            self.edge_labels_max += 1

        self.edge_labels[l_int] = l
        return l_int

    def add_edge(self, object u, object v, object l, bint directed):
        """
        Adds the edge ``(u,v)`` to self.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label
         - ``directed`` - if False, also add ``(v,u)``

        EXAMPLE::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        TESTS::

            sage: D = DiGraph(implementation='c_graph', sparse=True)
            sage: D.add_edge(0,1,2)
            sage: D.add_edge(0,1,3)
            sage: D.edges()
            [(0, 1, 3)]

        """
        if u is None: u = self.add_vertex(None)
        if v is None: v = self.add_vertex(None)

        cdef int u_int = self.check_labelled_vertex(u, self._directed)
        cdef int v_int = self.check_labelled_vertex(v, self._directed)

        cdef int l_int
        if l is None:
            l_int = 0
        else:
            l_int = self.new_edge_label(l)

        if (not self.loops(None)) and u_int == v_int:
            return
        if not self.multiple_edges(None):
            if self._cg.has_arc_label(u_int, v_int, l_int):
                return
            else:
                self._cg.del_all_arcs(u_int, v_int)
                if not directed:
                    self._cg.del_all_arcs(v_int, u_int)
        if directed:
            self._cg.add_arc_label(u_int, v_int, l_int)
            self._cg_rev.add_arc_label(v_int, u_int, l_int)
        elif u_int == v_int:
            self._cg.add_arc_label(u_int, v_int, l_int)
        else:
            self._cg.add_arc_label(u_int, v_int, l_int)
            self._cg.add_arc_label(v_int, u_int, l_int)

    def add_edges(self, object edges, bint directed):
        """
        Add edges from a list.

        INPUT:

         - ``edges`` - the edges to be added - can either be of the form
           ``(u,v)`` or ``(u,v,l)``
         - ``directed`` - if False, add ``(v,u)`` as well as ``(u,v)``

        EXAMPLE::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]

        """
        cdef object u,v,l,e
        for e in edges:
            try:
                u,v,l = e
            except Exception:
                u,v = e
                l = None
            self.add_edge(u,v,l,directed)

    def del_edge(self, object u, object v, object l, bint directed):
        """
        Delete edge ``(u,v,l)``.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label
         - ``directed`` - if False, also delete ``(v,u,l)``

        EXAMPLE::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]
            sage: D.del_edge(0,1,None,True)
            sage: list(D.iterator_out_edges(range(9), True))
            [(1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]

        TESTS::

            sage: G = Graph(implementation='c_graph', sparse=True)
            sage: G.add_edge(0,1,2)
            sage: G.delete_edge(0,1)
            sage: G.edges()
            []

            sage: G = Graph(multiedges=True, implementation='c_graph', sparse=True)
            sage: G.add_edge(0,1,2)
            sage: G.add_edge(0,1,None)
            sage: G.delete_edge(0,1)
            sage: G.edges()
            [(0, 1, 2)]

        Do we remove loops correctly? (:trac:`12135`)::

            sage: g=Graph({0:[0,0,0]}, implementation='c_graph', sparse=True)
            sage: g.edges(labels=False)
            [(0, 0), (0, 0), (0, 0)]
            sage: g.delete_edge(0,0); g.edges(labels=False)
            [(0, 0), (0, 0)]
        """
        if not ( self.has_vertex(u) and self.has_vertex(v) ):
            return
        cdef int u_int = self.check_labelled_vertex(u, self._directed)
        cdef int v_int = self.check_labelled_vertex(v, self._directed)

        if l is None:
            if self._cg.has_arc_label(u_int, v_int, 0):
                l_int = 0
            else:
                l_int = self._cg.arc_label(u_int, v_int)
        else:
            for l_int in self.edge_labels:
                if self.edge_labels[l_int] == l and self._cg.has_arc_label(u_int, v_int, l_int):
                    break
            else:
                return

        if directed:
            self._cg.del_arc_label(u_int, v_int, l_int)
            self._cg_rev.del_arc_label(v_int, u_int, l_int)
            if l_int:
                self.edge_labels.pop(l_int)
                self.edge_labels_available_ids.append(l_int)
        else:
            self._cg.del_arc_label(u_int, v_int, l_int)
            if v_int != u_int: self._cg.del_arc_label(v_int, u_int, l_int)
            if l_int:
                self.edge_labels.pop(l_int)
                self.edge_labels_available_ids.append(l_int)

    def get_edge_label(self, object u, object v):
        """
        Returns the edge label for ``(u,v)``.

        INPUT:

         - ``u,v`` - the vertices of the edge

        EXAMPLE::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1,1), (2,3,2), (4,5,3), (5,6,2)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, 1), (2, 3, 2), (4, 5, 3), (5, 6, 2)]
            sage: D.get_edge_label(3,2)
            2

        """
        cdef int l_int
        if not self.has_vertex(u):
            raise LookupError("({0}) is not a vertex of the graph.".format(repr(u)))
        if not self.has_vertex(v):
            raise LookupError("({0}) is not a vertex of the graph.".format(repr(v)))
        cdef int u_int = self.get_vertex(u)
        cdef int v_int = self.get_vertex(v)
        if not (<SparseGraph>self._cg).has_arc_unsafe(u_int, v_int):
            raise LookupError("({0}, {1}) is not an edge of the graph.".format(repr(u),repr(v)))
        if self.multiple_edges(None):
            return [self.edge_labels[l_int] if l_int != 0 else None
                     for l_int in self._cg.all_arcs(u_int, v_int)]
        l_int = self._cg.arc_label(u_int, v_int)
        return self.edge_labels[l_int] if l_int else None

    def has_edge(self, object u, object v, object l):
        """
        Returns whether this graph has edge ``(u,v)`` with label ``l``. If ``l``
        is ``None``, return whether this graph has an edge ``(u,v)`` with any
        label.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label, or ``None``

        EXAMPLE::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: D.has_edge(0,1,None)
            True

        """
        if not ( self.has_vertex(u) and self.has_vertex(v) ):
            return False
        cdef int u_int = self.get_vertex(u)
        cdef int v_int = self.get_vertex(v)
        if l is None:
            return self._cg.has_arc(u_int, v_int)
        for l_int in self._cg.all_arcs(u_int, v_int):
            if l_int and self.edge_labels[l_int] == l:
                return True
        return False

    def iterator_edges(self, object vertices, bint labels):
        """
        Iterate over the edges incident to a sequence of vertices. Edges are
        assumed to be undirected.

        INPUT:

        - ``vertices`` - a list of vertex labels
        - ``labels`` - boolean, whether to return labels as well

        EXAMPLE::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,3,False)
            sage: list(G.iterator_edges(range(9), False))
            [(1, 2)]
            sage: list(G.iterator_edges(range(9), True))
            [(1, 2, 3)]

        TEST::

            sage: g = graphs.PetersenGraph()
            sage: g.edges_incident([0,1,2])
            [(0, 1, None),
             (0, 4, None),
             (0, 5, None),
             (1, 2, None),
             (1, 6, None),
             (2, 3, None),
             (2, 7, None)]
        """
        cdef object u, v, l
        cdef int u_int, v_int, l_int
        cdef FrozenBitset b_vertices

        # ALL edges
        if not isinstance(vertices, list):
            if labels:
                for v in self.iterator_verts():
                    v_int = self.get_vertex(v)
                    for u_int, l_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, True):
                        if u_int >= v_int:
                            u = self.vertex_label(u_int)
                            l = self.edge_labels[l_int] if l_int else None
                            yield (v, u, l) if v<=u else (u, v, l)

            else:
                for v in self.iterator_verts():
                    v_int = self.get_vertex(v)
                    for u_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, False):
                        if u_int >= v_int:
                            u = self.vertex_label(u_int)
                            yield (v, u) if v <= u else (u, v)

        # One vertex
        elif len(vertices) == 1:
            v = vertices[0]
            v_int = self.get_vertex(v)

            if labels:
                for u_int, l_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, True):
                    u = self.vertex_label(u_int)
                    l = self.edge_labels[l_int] if l_int else None
                    yield (v, u, l) if v<=u else (u, v, l)
            else:
                for u_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, False):
                    u = self.vertex_label(u_int)
                    yield (v, u) if v <= u else (u, v)

        # Several vertices (nonempty list)
        elif vertices:
            b_vertices = FrozenBitset([self.get_vertex(v) for v in vertices])
            if labels:
                for v in vertices:
                    v_int = self.get_vertex(v)

                    for u_int, l_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, True):
                        if u_int >= v_int or u_int not in b_vertices:
                            u = self.vertex_label(u_int)
                            l = self.edge_labels[l_int] if l_int else None
                            yield (v, u, l) if v<=u else (u, v, l)
            else:
                for v in vertices:
                    v_int = self.get_vertex(v)
                    for u_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, False):
                        if u_int >= v_int or u_int not in b_vertices:
                            u = self.vertex_label(u_int)
                            yield (v, u) if v <= u else (u, v)

    def iterator_in_edges(self, object vertices, bint labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` - a list of vertex labels
        - ``labels`` - boolean, whether to return labels as well

        EXAMPLE::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,3,True)
            sage: list(G.iterator_in_edges([1], False))
            []
            sage: list(G.iterator_in_edges([2], False))
            [(1, 2)]
            sage: list(G.iterator_in_edges([2], True))
            [(1, 2, 3)]

        """
        cdef object u, v, L, l
        vertices = [self.get_vertex(v) for v in vertices if self.has_vertex(v)]
        cdef int u_int, v_int, l_int
        if self.multiple_edges(None):
            if labels:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int, l_int in (<SparseGraph> self._cg_rev).out_arcs_unsafe(v_int, True):
                        u = self.vertex_label(u_int)
                        l = self.edge_labels[l_int] if l_int else None
                        yield (u, v, l)
            else:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int in (<SparseGraph> self._cg_rev).out_arcs_unsafe(v_int, False):
                        u = self.vertex_label(u_int)
                        yield (u, v)
        else:
            if labels:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int in self._cg_rev.out_neighbors(v_int):
                        l_int = self._cg.arc_label(u_int, v_int)
                        yield (self.vertex_label(u_int),
                               v,
                               None if l_int == 0 else self.edge_labels[l_int])
            else:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int in self._cg_rev.out_neighbors(v_int):
                        yield (self.vertex_label(u_int),
                               v)

    def iterator_out_edges(self, object vertices, bint labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:
         - ``vertices`` - a list of vertex labels
         - ``labels`` - boolean, whether to return labels as well

        EXAMPLE::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,3,True)
            sage: list(G.iterator_out_edges([2], False))
            []
            sage: list(G.iterator_out_edges([1], False))
            [(1, 2)]
            sage: list(G.iterator_out_edges([1], True))
            [(1, 2, 3)]

        """
        cdef object u, v, L, l
        vertices = [self.get_vertex(v) for v in vertices if self.has_vertex(v)]
        cdef int u_int, v_int, l_int
        if self.multiple_edges(None):
            if labels:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int, l_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, True):
                        u = self.vertex_label(u_int)
                        l = self.edge_labels[l_int] if l_int else None
                        yield (v, u, l)
            else:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int in (<SparseGraph> self._cg).out_arcs_unsafe(v_int, False):
                        u = self.vertex_label(u_int)
                        yield (v, u)
        else:
            if labels:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int in self._cg.out_neighbors(v_int):
                        l_int = self._cg.arc_label(v_int, u_int)
                        yield (v,
                               self.vertex_label(u_int),
                               None if l_int == 0 else self.edge_labels[l_int])
            else:
                for v_int in vertices:
                    v = self.vertex_label(v_int)
                    for u_int in self._cg.out_neighbors(v_int):
                        yield (v,
                               self.vertex_label(u_int))

    def multiple_edges(self, new):
        """
        Get/set whether or not ``self`` allows multiple edges.

        INPUT:

         - ``new`` - boolean (to set) or ``None`` (to get)

        EXAMPLES::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.multiple_edges(True)
            sage: G.multiple_edges(None)
            True
            sage: G.multiple_edges(False)
            sage: G.multiple_edges(None)
            False
            sage: G.add_edge(0,1,0,True)
            sage: G.add_edge(0,1,0,True)
            sage: list(G.iterator_edges(range(9), True))
            [(0, 1, 0)]

        """
        if new is None:
            return self._multiple_edges
        self._multiple_edges = bool(new)

    def set_edge_label(self, object u, object v, object l, bint directed):
        """
        Label the edge ``(u,v)`` by ``l``.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label
         - ``directed`` - if False, also set ``(v,u)`` with label ``l``

        EXAMPLE::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,None,True)
            sage: G.set_edge_label(1,2,'a',True)
            sage: list(G.iterator_edges(range(9), True))
            [(1, 2, 'a')]

        Note that it fails silently if there is no edge there::

            sage: G.set_edge_label(2,1,'b',True)
            sage: list(G.iterator_edges(range(9), True))
            [(1, 2, 'a')]

        """
        if not self.has_edge(u, v, None):
            return
        if self.multiple_edges(None):
            if len(self.get_edge_label(u, v)) > 1:
                raise RuntimeError("Cannot set edge label, since there are multiple edges from %s to %s."%(u,v))
        # now we know there is exactly one edge from u to v
        cdef int l_int, ll_int
        if l is None:
            l_int = 0
        else:
            l_int = self.new_edge_label(l)
        cdef int u_int = self.get_vertex(u)
        cdef int v_int = self.get_vertex(v)
        if not (<SparseGraph>self._cg).has_arc_unsafe(u_int, v_int):
            return
        ll_int = (<SparseGraph>self._cg).arc_label_unsafe(u_int, v_int)
        if ll_int:
            self.edge_labels.pop(ll_int)
            self.edge_labels_available_ids.append(ll_int)
        if directed:
            self._cg.del_arc_label(u_int, v_int, ll_int)
            self._cg_rev.del_arc_label(v_int, u_int, ll_int)
            self._cg.add_arc_label(u_int, v_int, l_int)
            self._cg_rev.add_arc_label(v_int, u_int, l_int)
        elif u_int == v_int:
            self._cg.del_arc_label(u_int, v_int, ll_int)
            self._cg.add_arc_label(u_int, v_int, l_int)
        else:
            self._cg.del_arc_label(u_int, v_int, ll_int)
            self._cg.del_arc_label(v_int, u_int, ll_int)
            self._cg.add_arc_label(u_int, v_int, l_int)
            self._cg.add_arc_label(v_int, u_int, l_int)
