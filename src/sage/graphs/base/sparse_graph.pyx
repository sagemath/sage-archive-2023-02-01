
#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

from c_graph import CGraphBackend

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
    Implements compiled sparse graphs via an array of hash tables.

    Creating a new sparse graph instance:
        G = SparseGraph(int nverts, int expected_degree = 16)

    INPUT:
        nverts -- non-negative integer, number of vertices.
        expected_degree -- non-negative integer, expected loose bound on
    degree of vertices. This need not be very accurate, but it does affect
    the footprint. You may also get a speed improvement from this.

    NOTES:
        SparseGraph does not distinguish whether it is directed or not. In fact,
    the datastructure itself is directed. An edge is simply an arc in both
    directions.

        If you are planning on using multiple arcs, it is a good idea to
    give all of the arcs labels. When you add a labeled arc (u, v, l) when
    there is already an unlabeled arc (u,v), it simply adds the label l to the
    arc without making a new one. Multiple unlabeled arcs are not allowed.

    """

    def __new__(self, int nverts, int expected_degree = 16):
        cdef int i = 1
        self.num_verts = nverts
        self.num_arcs = 0
        while i < expected_degree:
            i = i << 1
        self.hash_length = i
        self.hash_mask = i - 1
        self.vertices = <SparseGraphBTNode **> \
          sage_malloc(nverts * self.hash_length * sizeof(SparseGraphBTNode *))
        self.in_degrees = <int *> sage_malloc(nverts * sizeof(int))
        self.out_degrees = <int *> sage_malloc(nverts * sizeof(int))
        if not self.vertices or not self.in_degrees or not self.out_degrees:
            if self.vertices: sage_free(self.vertices)
            if self.in_degrees: sage_free(self.in_degrees)
            if self.out_degrees: sage_free(self.out_degrees)
            raise RuntimeError("Failure allocating memory.")
        for i from 0 <= i < nverts * self.hash_length:
            self.vertices[i] = NULL
        for i from 0 <= i < nverts:
            self.in_degrees[i] = 0
            self.out_degrees[i] = 0

    def __dealloc__(self):
        cdef SparseGraphBTNode **temp
        cdef SparseGraphLLNode *label_temp
        cdef int i
        for i from 0 <= i < self.num_verts * self.hash_length:
            temp = &(self.vertices[i])
            while temp[0] != NULL:
                if temp[0].left != NULL:
                    temp = &(temp[0].left)
                elif temp[0].right != NULL:
                    temp = &(temp[0].right)
                else:
                    label_temp = temp[0].labels
                    while label_temp != NULL:
                        temp[0].labels = label_temp.next
                        sage_free(label_temp)
                        label_temp = temp[0].labels
                    sage_free(temp[0])
                    temp[0] = NULL
                    temp = &(self.vertices[i])
        sage_free(self.vertices)
        sage_free(self.in_degrees)
        sage_free(self.out_degrees)

    cdef int add_arc_unsafe(self, int u, int v) except? -1:
        """
        Adds arc (u, v) to the graph with no label.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        OUTPUT:
            0 -- No error.
            1 -- Arc already exists (nothing was done).

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
                return 1
        ins_pt[0] = <SparseGraphBTNode *> sage_malloc(sizeof(SparseGraphBTNode))
        if not ins_pt[0]:
            raise RuntimeError("Failure allocating memory.")
        ins_pt[0].vertex = v
        ins_pt[0].left = NULL
        ins_pt[0].right = NULL
        ins_pt[0].labels = NULL
        self.in_degrees[v] += 1
        self.out_degrees[u] += 1
        self.num_arcs += 1

    def add_arc(self, int u, int v):
        """
        Adds arc (u, v) to the graph with no label.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(4,7)
            Traceback (most recent call last):
            ...
            RuntimeError: Second vertex (7) is not a vertex of the graph.
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        self.add_arc_unsafe(u,v)

    cdef int has_arc_unsafe(self, int u, int v):
        """
        Checks whether arc (u, v) is in the graph.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

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

    def has_arc(self, int u, int v):
        """
        Checks whether arc (u, v) is in the graph.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1)
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        if u < 0 or u >= self.num_verts or v < 0 or v >= self.num_verts:
            return False
        return self.has_arc_unsafe(u,v) == 1

    cdef int del_arc_unsafe(self, int u, int v):
        """
        Deletes *all* arcs from u to v.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        OUTPUT:
            0 -- No error.
            1 -- No arc to delete.

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared, left_len, right_len
        cdef SparseGraphBTNode *temp, **left_child, **right_child
        cdef SparseGraphBTNode **parent = &self.vertices[i]
        cdef SparseGraphLLNode *labels
        while parent[0] != NULL:
            compared = compare(parent[0].vertex, v)
            if compared > 0:
                parent = &(parent[0].left)
            elif compared < 0:
                parent = &(parent[0].right)
            else:# if parent[0].vertex == v:
                break
        if parent[0] == NULL:
            return 1 # indicate there is no arc to delete
        labels = parent[0].labels
        if labels == NULL:
            self.in_degrees[v] -= 1
            self.out_degrees[u] -= 1
            self.num_arcs -= 1
        while labels != NULL:
            parent[0].labels = parent[0].labels.next
            sage_free(labels)
            labels = parent[0].labels
            self.in_degrees[v] -= 1
            self.out_degrees[u] -= 1
            self.num_arcs -= 1
        if parent[0].left == NULL:
            temp = parent[0]
            parent[0] = parent[0].right
            sage_free(temp)
            return 0
        elif parent[0].right == NULL:
            temp = parent[0]
            parent[0] = parent[0].left
            sage_free(temp)
            return 0
        else:
            left_len = 0
            right_len = 0
            left_child = &(parent[0].left)
            while left_child[0].right != NULL:
                left_len += 1
                left_child = &(left_child[0].right)
            right_child = &(parent[0].right)
            while right_child[0].left != NULL:
                right_len += 1
                right_child = &(right_child[0].left)
            if left_len > right_len:
                left_child[0].right = parent[0].right
                temp = parent[0]
                parent[0] = left_child[0]
                left_child[0] = left_child[0].left
                parent[0].left = temp.left
                sage_free(temp)
                return 0
            else:
                right_child[0].left = parent[0].left
                temp = parent[0]
                parent[0] = right_child[0]
                right_child[0] = right_child[0].right
                parent[0].right = temp.right
                sage_free(temp)
                return 0

    def del_all_arcs(self, int u, int v):
        """
        Deletes all arcs from u to v.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        EXAMPLE:
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
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        self.del_arc_unsafe(u,v)

    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2:
        """
        Gives all v such that (u, v) is an arc of the graph.

        INPUT:
            u -- integer from 0, ..., n-1, where n is the number of vertices
            neighbors -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            nonnegative integer -- the number of v such that (u, v) is an arc
            -1 -- indicates that the array has been filled with neighbors, but
        there were more

        """
        cdef int i, num_nbrs = 0, current_nbr = 0
        if self.out_degrees[u] == 0:
            return 0
        cdef SparseGraphBTNode **pointers = <SparseGraphBTNode **> \
          sage_malloc(size * sizeof(SparseGraphBTNode *))
        if not pointers:
            raise RuntimeError("Failure allocating memory.")
        for i from u * self.hash_length <= i < (u+1) * self.hash_length:
            if self.vertices[i] == NULL:
                continue
            if num_nbrs == size:
                sage_free(pointers)
                return -1
            pointers[num_nbrs] = self.vertices[i]
            neighbors[num_nbrs] = self.vertices[i].vertex
            num_nbrs += 1
            while current_nbr < num_nbrs:
                if pointers[current_nbr].left != NULL:
                    if num_nbrs == size:
                        sage_free(pointers)
                        return -1
                    pointers[num_nbrs] = pointers[current_nbr].left
                    neighbors[num_nbrs] = pointers[current_nbr].left.vertex
                    num_nbrs += 1
                if pointers[current_nbr].right != NULL:
                    if num_nbrs == size:
                        sage_free(pointers)
                        return -1
                    pointers[num_nbrs] = pointers[current_nbr].right
                    neighbors[num_nbrs] = pointers[current_nbr].right.vertex
                    num_nbrs += 1
                current_nbr += 1
        sage_free(pointers)
        return num_nbrs

    def out_neighbors(self, int u):
        """
        Gives all v such that (u, v) is an arc of the graph.

        INPUT:
            u -- integer from 0, ..., n-1, where n is the number of vertices

        EXAMPLES:
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
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("Vertex (%d) is not a vertex of the graph."%u)
        if self.out_degrees[u] == 0:
            return []
        cdef int size = self.out_degrees[u]
        cdef int *neighbors = <int *> sage_malloc(size * sizeof(int))
        if not neighbors:
            raise RuntimeError("Failure allocating memory.")
        num_nbrs = self.out_neighbors_unsafe(u, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sage_free(neighbors)
        return output

    cdef int in_neighbors_unsafe(self, int v, int *neighbors, int size):
        """
        Gives all u such that (u, v) is an arc of the graph.

        INPUT:
            v -- integer from 0, ..., n-1, where n is the number of vertices
            neighbors -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            nonnegative integer -- the number of u such that (u, v) is an arc
            -1 -- indicates that the array has been filled with neighbors, but
        there were more

        """
        cdef int i, num_nbrs = 0
        if self.in_degrees[v] == 0:
            return 0
        for i from 0 <= i < self.num_verts:
            if self.has_arc_unsafe(i, v):
                if num_nbrs == size:
                    return -1
                neighbors[num_nbrs] = i
                num_nbrs += 1
        return num_nbrs

    def in_neighbors(self, int v):
        """
        Gives all u such that (u, v) is an arc of the graph.

        INPUT:
            v -- integer from 0, ..., n-1, where n is the number of vertices

        EXAMPLES:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(3,1)
            sage: G.add_arc(1,3)
            sage: G.in_neighbors(1)
            [0, 3]
            sage: G.in_neighbors(3)
            [1]

        """
        cdef int i, num_nbrs
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Vertex (%d) is not a vertex of the graph."%v)
        if self.in_degrees[v] == 0:
            return []
        cdef int size = self.in_degrees[v]
        cdef int *neighbors = <int *> sage_malloc(size * sizeof(int))
        if not neighbors:
            raise RuntimeError("Failure allocating memory.")
        num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sage_free(neighbors)
        return output

    cdef int add_arc_label_unsafe(self, int u, int v, int l) except? -1:
        """
        Adds arc (u, v) to the graph with label l. To add an arc without a
        label, use l = 0, or add_arc_unsafe or add_arc. Note that if there is an
        unlabeled arc (u,v), this function simply labels that arc with l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label, or zero for no label

        OUTPUT:
            0 -- No error.
            1 -- No label specified, and arc already exists (nothing was done).

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
            ins_pt[0] = <SparseGraphBTNode *> sage_malloc(sizeof(SparseGraphBTNode))
            if not ins_pt[0]:
                raise RuntimeError("Failure allocating memory.")
            if l:
                label_ptr = <SparseGraphLLNode *> sage_malloc(sizeof(SparseGraphLLNode))
                if not label_ptr:
                    sage_free(ins_pt[0])
                    raise RuntimeError("Failure allocating memory.")
            ins_pt[0].vertex = v
            ins_pt[0].left = NULL
            ins_pt[0].right = NULL
            ins_pt[0].labels = NULL
            self.in_degrees[v] += 1
            self.out_degrees[u] += 1
            self.num_arcs += 1
        elif not l:
            return 1
        elif ins_pt[0].labels != NULL:
            label_ptr = <SparseGraphLLNode *> sage_malloc(sizeof(SparseGraphLLNode))
            if not label_ptr:
                raise RuntimeError("Failure allocating memory.")
            self.in_degrees[v] += 1
            self.out_degrees[u] += 1
            self.num_arcs += 1
        else:
            label_ptr = <SparseGraphLLNode *> sage_malloc(sizeof(SparseGraphLLNode))
            if not label_ptr:
                raise RuntimeError("Failure allocating memory.")
        if l:
            label_ptr.next = ins_pt[0].labels
            ins_pt[0].labels = label_ptr
            ins_pt[0].labels.label = l

    def add_arc_label(self, int u, int v, int l=0):
        """
        Adds arc (u, v) to the graph with label l. To add an arc without a
        label (the default), use l = 0, or add_arc. Note that if there is an
        unlabeled arc (u,v), this function simply labels that arc with l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label, or zero for no label

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1)
            sage: G.add_arc_label(4,7)
            Traceback (most recent call last):
            ...
            RuntimeError: Second vertex (7) is not a vertex of the graph.
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True
            sage: G.add_arc_label(1,2,2)
            sage: G.arc_label(1,2)
            2

        """
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        if l < 0:
            raise RuntimeError("Label (%d) must be a nonnegative integer."%l)
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

    def arc_label(self, int u, int v):
        """
        Retrieves the first label found associated with (u, v) (a positive
        integer).

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        OUTPUT:
            positive integer -- indicates that there is a label on (u, v).
            0 -- either the arc (u, v) is unlabeled, or there is no arc at all.

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(3,4,7)
            sage: G.arc_label(3,4)
            7

        NOTES:
        To this function, an unlabeled arc is indistinguishable from a non-arc:
            sage: G.add_arc_label(1,0)
            sage: G.arc_label(1,0)
            0
            sage: G.arc_label(1,1)
            0

        This function only returns the *first* label it finds from u to v:
            sage: G.add_arc_label(1,2,1)
            sage: G.add_arc_label(1,2,2)
            sage: G.arc_label(1,2)
            2

        """
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        return self.arc_label_unsafe(u,v)

    cdef int all_arcs_unsafe(self, int u, int v, int *arc_labels, int size):
        """
        Gives the labels of all arcs (u, v).

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            arc_labels -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            positive integer -- the number of labels on the arc (u, v)
            0 -- if there is an unlabeled arc (u, v)
            -1 -- if there is no arc (u, v) at all
            -2 -- indicates that the array has been filled with labels, but
        there were more

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared, num_arcs = 0
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
            return -1
        label = temp.labels
        while label != NULL and num_arcs < size:
            arc_labels[num_arcs] = label.label
            label = label.next
            num_arcs += 1
        if num_arcs == size and label != NULL:
            return -2
        return num_arcs

    def all_arcs(self, int u, int v):
        """
        Gives the labels of all arcs (u, v). An unlabeled arc is interpreted as
        having label None.

        EXAMPLE:
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
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        if self.in_degrees[v] < self.out_degrees[u]:
            size = self.in_degrees[v]
        else:
            size = self.out_degrees[u]
        arc_labels = <int *> sage_malloc(size * sizeof(int))
        if not arc_labels:
            raise RuntimeError("Failure allocating memory.")
        num_arcs = self.all_arcs_unsafe(u, v, arc_labels, size)
        if num_arcs == -2:
            sage_free(arc_labels)
            raise RuntimeError("There was an error: there seem to be more arcs than self.in_degrees or self.out_degrees indicate.")
        if num_arcs == -1:
            sage_free(arc_labels)
            return []
        if num_arcs == 0:
            sage_free(arc_labels)
            return [None]
        output = [arc_labels[i] for i from 0 <= i < num_arcs]
        sage_free(arc_labels)
        return output

    cdef int del_arc_label_unsafe(self, int u, int v, int l):
        """
        Delete an arc (u, v) with label l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label

        OUTPUT:
            0 -- No error.
            1 -- No arc with label l.

        """
        cdef int i = (u * self.hash_length) + (v & self.hash_mask)
        cdef int compared
        cdef SparseGraphBTNode **parent = &self.vertices[i]
        cdef SparseGraphLLNode **labels, *label
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
        labels = &(parent[0].labels)
        while labels[0] != NULL and labels[0].label != l:
            labels = &(labels[0].next)
        if labels[0] == NULL:
            return 1
        label = labels[0]
        labels[0] = labels[0].next
        sage_free(label)
        if labels == &(parent[0].labels) and labels[0] == NULL:
            self.del_arc_unsafe(u, v)
        else:
            self.in_degrees[v] -= 1
            self.out_degrees[u] -= 1
            self.num_arcs -= 1

    def del_arc_label(self, int u, int v, int l):
        """
        Delete an arc (u, v) with label l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1,0)
            sage: G.add_arc_label(0,1,1)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,3)
            sage: G.del_arc_label(0,1,2)
            sage: G.all_arcs(0,1)
            [3, 2, 1]

        """
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        if l < 1:
            raise RuntimeError("Label (%d) must be a positive integer."%l)
        self.del_arc_label_unsafe(u,v,l)

    cdef int has_arc_label_unsafe(self, int u, int v, int l):
        """
        Indicates whether there is an arc (u, v) with label l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label

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
        label = temp.labels
        while label != NULL:
            if label.label == l:
                return 1
            label = label.next
        return 0

    def has_arc_label(self, int u, int v, int l):
        """
        Indicates whether there is an arc (u, v) with label l.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices
            l -- a positive integer label

        EXAMPLE:
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
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        if l < 1:
            raise RuntimeError("Label (%d) must be a positive integer."%l)
        return self.has_arc_label_unsafe(u,v,l) == 1

    cdef int del_vertex_unsafe(self, int v):
        """
        Deletes the vertex v, along with all edges incident to it. The vertices
        larger than v, i.e. {v+1, ..., n-1} are renamed to {v, ..., n-2}.

        INPUT:
            v -- integer from 0, ..., n-1, where n is the number of vertices

        """
        cdef int size = 0, num_nbrs, i, j, w, *neighbors, compared
        cdef SparseGraphBTNode *temp
        cdef SparseGraphLLNode *label_temp
        for w from v <= w < self.num_verts:
            if self.in_degrees[w] > size:
                size = self.in_degrees[w]
            if self.out_degrees[w] > size:
                size = self.out_degrees[w]
        neighbors = <int *> sage_malloc(size * sizeof(int))
        if not neighbors:
            raise RuntimeError("Failure allocating memory.")
        # delete each arc incident with v
        num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
        for i from 0 <= i < num_nbrs:
            self.del_arc_unsafe(neighbors[i], v)
        num_nbrs = self.out_neighbors_unsafe(v, neighbors, size)
        for i from 0 <= i < num_nbrs:
            self.del_arc_unsafe(v, neighbors[i])
        # readjust the arcs greater than v
        for w from v+1 <= w < self.num_verts:
            # move the incoming arcs from w to w-1
            num_nbrs = self.in_neighbors_unsafe(w, neighbors, size)
            self.in_degrees[w-1] = self.in_degrees[w] - num_nbrs # adjust for arc additions
            for i from 0 <= i < num_nbrs:
                j = (neighbors[i] * self.hash_length) + (w & self.hash_mask)
                temp = self.vertices[j]
                while True:
                    compared = compare(temp.vertex, w)
                    if compared > 0:
                        temp = temp.left
                    elif compared < 0:
                        temp = temp.right
                    else:
                        break
                label_temp = temp.labels
                temp.labels = NULL
                self.del_arc_unsafe(neighbors[i], w)
                self.add_arc_unsafe(neighbors[i], w-1)
                j = (neighbors[i] * self.hash_length) + ((w-1) & self.hash_mask)
                temp = self.vertices[j]
                while True:
                    compared = compare(temp.vertex, w-1)
                    if compared > 0:
                        temp = temp.left
                    elif compared < 0:
                        temp = temp.right
                    else:
                        break
                temp.labels = label_temp
            # move the outgoing arcs from w to w-1
            num_nbrs = self.out_neighbors_unsafe(w, neighbors, size)
            self.out_degrees[w-1] = self.out_degrees[w] - num_nbrs # adjust for arc additions
            for i from 0 <= i < num_nbrs:
                j = (w * self.hash_length) + (neighbors[i] & self.hash_mask)
                temp = self.vertices[j]
                while True:
                    compared = compare(temp.vertex, neighbors[i])
                    if compared > 0:
                        temp = temp.left
                    elif compared < 0:
                        temp = temp.right
                    else:
                        break
                label_temp = temp.labels
                temp.labels = NULL
                self.del_arc_unsafe(w, neighbors[i])
                self.add_arc_unsafe(w-1, neighbors[i])
                j = ((w-1) * self.hash_length) + (neighbors[i] & self.hash_mask)
                temp = self.vertices[j]
                while True:
                    compared = compare(temp.vertex, neighbors[i])
                    if compared > 0:
                        temp = temp.left
                    elif compared < 0:
                        temp = temp.right
                    else:
                        break
                temp.labels = label_temp
        self.num_verts -= 1
        self.vertices = <SparseGraphBTNode **> \
          sage_realloc(self.vertices, self.num_verts * self.hash_length * sizeof(SparseGraphBTNode *))
        self.in_degrees = <int *> sage_realloc(self.in_degrees, self.num_verts * sizeof(int))
        self.out_degrees = <int *> sage_realloc(self.out_degrees, self.num_verts * sizeof(int))
        sage_free(neighbors)

    def del_vertex(self, v):
        """
        Deletes the vertex v, along with all edges incident to it. The vertices
        larger than v, i.e. {v+1, ..., n-1} are renamed to {v, ..., n-2}.

        EXAMPLES:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3)
            sage: G.add_arc(0,1)
            sage: G.add_arc(0,2)
            sage: G.add_arc(1,2)
            sage: G.add_arc(2,0)
            sage: G.del_vertex(2)
            sage: for i in range(2):
            ...    for j in range(2):
            ...        if G.has_arc(i,j):
            ...            print i,j
            0 1
            sage: G = SparseGraph(3)
            sage: G.add_arc(0,1)
            sage: G.add_arc(0,2)
            sage: G.add_arc(1,2)
            sage: G.add_arc(2,0)
            sage: G.del_vertex(1)
            sage: for i in range(2):
            ...    for j in range(2):
            ...        if G.has_arc(i,j):
            ...            print i,j
            0 1
            1 0

        """
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Vertex (%d) is not a vertex of the graph."%v)
        self.del_vertex_unsafe(v)

    cdef int add_vertices_unsafe(self, int k):
        """
        Adds k vertices to the graph. If the vertices are {0, ..., n-1}, then the
        vertices {n, n+1, ..., n+k-1} are added.

        INPUT:
            k -- positive integer

        """
        cdef int i, *in_degrees, *out_degrees
        cdef SparseGraphBTNode **vertices
        vertices = <SparseGraphBTNode **> \
          sage_realloc(self.vertices, (self.num_verts+k) * self.hash_length * sizeof(SparseGraphBTNode *))
        in_degrees = <int *> sage_realloc(self.in_degrees, (self.num_verts+k) * sizeof(int))
        out_degrees = <int *> sage_realloc(self.out_degrees, (self.num_verts+k) * sizeof(int))
        if vertices: self.vertices = vertices
        if in_degrees: self.in_degrees = in_degrees
        if out_degrees: self.out_degrees = out_degrees
        if not vertices or not in_degrees or not out_degrees:
            raise RuntimeError("Error reallocating memory.")
        for i from self.num_verts <= i < self.num_verts + k:
            self.in_degrees[i] = 0
            self.out_degrees[i] = 0
        for i from self.num_verts * self.hash_length <= i < (self.num_verts + k) * self.hash_length:
            self.vertices[i] = NULL
        self.num_verts += k

    def add_vertices(self, k):
        """
        Adds k vertices to the graph. If the vertices are {0, ..., n-1}, then the
        vertices {n, n+1, ..., n+k-1} are added.

        INPUT:
            k -- positive integer

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3)
            sage: G.add_vertices(6)
            sage: G.add_arc(2,5)
            sage: G.add_arc(1,3)
            sage: G.has_arc(1,3)
            True
            sage: G.has_arc(2,3)
            False

        """
        if k < 0:
            raise RuntimeError("Input (%d) must be a positive integer."%k)
        if k == 0:
            return
        self.add_vertices_unsafe(k)


def random_stress():
    """
    Randomly search for mistakes in the code.

    EXAMPLE:
    No output indicates that no errors were found.
        sage: from sage.graphs.base.sparse_graph import random_stress
        sage: for _ in xrange(20):
        ...    random_stress()

    """
    cdef int i, j, k, l, m, n
    cdef SparseGraph Gnew
    num_verts = 10
    Gnew = SparseGraph(num_verts)
    from random import randint
    from sage.graphs.graph import DiGraph
    from sage.misc.misc import uniq
    Gold = DiGraph(multiedges=True, loops=True, implementation='networkx')
    Gold.add_vertices(xrange(num_verts))
    for n from 0 <= n < 100:
        i = randint(0,num_verts-1)
        j = randint(0,num_verts-1)
        l = randint(1,num_verts-1)
        k = randint(0,num_verts-1)
        if k > 7:
#            print 'G.add_arc_label(%d,%d,%d);'%( i, j, l ) + ' Gold.add_edge(%d,%d,%d)'%( i, j, l )
            Gold.add_edge(i,j,l)
            Gnew.add_arc_label_unsafe(i,j,l)
        elif k > 5:
            m = randint(1,7)
#            print 'G.add_vertices(%d); '%m + ' Gold.add_vertices(range(num_verts, num_verts+%d));'%m + ' num_verts += %d'%m
            Gold.add_vertices(range(num_verts, num_verts+m))
            Gnew.add_vertices(m)
            num_verts += m
        elif k > 3:
            m = randint(0,num_verts-1)
#            print 'G.del_vertex(%d); '%m + ' Gold.delete_vertex(%d); '%m + ' Gold.relabel(range(%d) + [\'unused\'] + range(%d,num_verts-1)); num_verts -= 1'%(m,m)
            Gold.delete_vertex(m)
            Gold.relabel(range(m) + ['unused'] + range(m,num_verts-1))
            Gnew.del_vertex(m)
            num_verts -= 1
        elif k > 1:
#            print 'G.del_all_arcs(%d,%d);'%( i, j ) + ' Gold.delete_edges([(u,v,ll) for u,v,ll in Gold.edges() if u==%d and v==%d])'%(i,j)
            Gold.delete_edges([(u,v,ll) for u,v,ll in Gold.edges() if u==i and v==j])
            Gnew.del_arc_unsafe(i,j)
        else:
#            print 'G.del_arc_label(%d,%d,%d);'%( i, j, l ) + ' Gold.delete_edge(%d,%d,%d)'%( i, j, l )
            Gold.delete_edge(i,j,l)
            Gnew.del_arc_label_unsafe(i,j,l)
    if Gnew.num_arcs != Gold.size():
        #print Gnew.num_arcs, Gold.size()
        raise RuntimeError( "NO:size" )
    for i from 0 <= i < num_verts:
        if Gnew.out_degrees[i] != Gold.out_degree(i):
            raise RuntimeError( "NO:out degree" )
        if Gnew.in_degrees[i] != Gold.in_degree(i):
            raise RuntimeError( "NO:in degree" )
        if sorted(Gnew.out_neighbors(i)) != uniq([v for u,v,l in Gold.outgoing_edge_iterator(i)]):
            raise RuntimeError( "NO:out neighbors" )
        if sorted(Gnew.in_neighbors(i)) != uniq([u for u,v,l in Gold.incoming_edge_iterator(i)]):
            raise RuntimeError( "NO:in neighbors" )
        for j from 0 <= j < num_verts:
            l = Gnew.arc_label_unsafe(i,j)
            if l != 0:
                if not Gold.has_edge(i,j,l):
                    raise RuntimeError( "NO:has_edge" )
            else:
                if Gold.has_edge(i,j):
                    raise RuntimeError( "NO:has_edge" )
            list1 = Gnew.all_arcs(i,j)
            list2 = [l for (u,v,l) in Gold.edges() if u==i and v==j]
            if sorted(list1) != sorted(list2):
                raise RuntimeError("NO:edges")
            for l from 1 <= l < num_verts:
                if Gold.has_edge(i,j,l) != Gnew.has_arc_label(i,j,l):
                    raise RuntimeError("NO:edges")
    Gnew = SparseGraph(num_verts)
    Gold = DiGraph(loops=True, implementation='networkx')
    Gold.add_vertices(xrange(num_verts))
    for n from 0 <= n < 100:
        i = randint(0,num_verts-1)
        j = randint(0,num_verts-1)
        k = randint(0,num_verts-1)
        if k != 0:
            Gold.add_edge(i,j)
            Gnew.add_arc_unsafe(i,j)
        else:
            Gold.delete_edge(i,j)
            Gnew.del_arc_unsafe(i,j)
    if Gnew.num_arcs != Gold.size():
        raise RuntimeError( "NO" )
    for i from 0 <= i < num_verts:
        if Gnew.out_degrees[i] != Gold.out_degree(i):
            raise RuntimeError( "NO" )
        if Gnew.in_degrees[i] != Gold.in_degree(i):
            raise RuntimeError( "NO" )
        if sorted(Gnew.out_neighbors(i)) != uniq([v for u,v,_ in Gold.outgoing_edge_iterator(i)]):
            raise RuntimeError( "NO" )
        if sorted(Gnew.in_neighbors(i)) != uniq([u for u,v,_ in Gold.incoming_edge_iterator(i)]):
            raise RuntimeError( "NO" )
        for j from 0 <= j < num_verts:
            if Gnew.has_arc_unsafe(i,j) != Gold.has_edge(i,j):
                raise RuntimeError( "NO" )


class SparseGraphBackend(CGraphBackend):

    def __init__(self, n):
        """
        Initialize a dense graph with n vertices.

        EXAMPLE:
            sage: import sage.graphs.base.sparse_graph
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None), (1, 0, None)]

        """
        self._cg = SparseGraph(n)

    def add_edge(self, u, v, l, directed):
        """
        Adds the edge u,v to self.

        EXAMPLE:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None), (1, 0, None)]

        """
        if not self.has_vertex(u) or not self.has_vertex(v):
            n = self.num_verts()
            if (u == n and v == n+1) or (v == n and u == n+1):
                self._cg.add_vertices(2)
            elif u == n or v == n:
                self._cg.add_vertices(1)
        if l is None:
            if directed:
                self._cg.add_arc(u, v)
            else:
                self._cg.add_arc(u, v)
                self._cg.add_arc(v, u)
        else:
            if directed:
                self._cg.add_arc_label(u, v, l)
            else:
                self._cg.add_arc_label(u, v, l)
                self._cg.add_arc_label(v, u, l)

    def add_edges(self, edges, directed):
        """
        Add edges from a list.

        EXAMPLE:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None),
             (1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]

        """
        for e in edges:
            try:
                u,v,l = e
            except:
                u,v = e
                l = None
            self.add_edge(u,v,l,directed)

    def add_vertex(self, name):
        """
        Add a labelled vertex to self.

        INPUT:
            name: vertex label

        DOCTEST:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_vertex(10)
        """
        if not self.has_vertex(name):
            self._cg.add_vertices(1)

    def add_vertices(self, vertices):
        """
        Add labelled vertices to self.

        INPUT:
            vertices: iterator of vertex labels

        DOCTEST:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(1)
            sage: D.add_vertices([1,2,3])
        """
        i = 0
        vertices = list(vertices)
        while i < len(vertices):
            if self.has_vertex(vertices[i]):
                vertices.pop(i)
            else:
                i += 1
        self._cg.add_vertices(len(list(vertices)))

    def del_edge(self, u, v, l, directed):
        """
        Delete edge u,v.

        EXAMPLE:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None),
             (1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]
            sage: D.del_edge(0,1,None,True)
            sage: list(D.iterator_edges(range(9), True, True))
            [(1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]

        """
        if v is None:
            try:
                u1,v,l = u[:3]
                u = u1
            except:
                u, v = u[:2]
                l = None
        if l is None:
            if directed:
                self._cg.del_all_arcs(u, v)
            else:
                self._cg.del_all_arcs(u, v)
                self._cg.del_all_arcs(v, u)
        else:
            if directed:
                self._cg.del_arc_label(u, v, l)
            else:
                self._cg.del_arc_label(u, v, l)
                self._cg.del_arc_label(v, u, l)

    def del_vertex(self, v):
        """
        Delete a labelled vertex in self.

        INPUT:
            v: vertex label

        DOCTEST:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.del_vertex(0)
        """
        self._cg.del_vertex(v)

    def del_vertices(self, vertices):
        """
        Delete labelled vertices in self.

        INPUT:
            vertices: iterator of vertex labels

        DOCTEST:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.del_vertices([1,2,3])
        """
        verts = sorted(vertices)
        i = 0
        while i < len(verts):
            self.del_vertex(verts[i]-i)
            i += 1

    def get_edge_label(self, u, v):
        """
        Returns the edge label for u,v.

        EXAMPLE:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1,1), (2,3,2), (4,5,3), (5,6,2)], False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, 1), (1, 0, 1), (2, 3, 2), (3, 2, 2), (4, 5, 3), (5, 4, 3), (5, 6, 2), (6, 5, 2)]
            sage: D.del_edge(0,1,None,True)
            sage: list(D.iterator_edges(range(9), True, True))
            [(1, 0, 1), (2, 3, 2), (3, 2, 2), (4, 5, 3), (5, 4, 3), (5, 6, 2), (6, 5, 2)]
            sage: D.get_edge_label(1,0)
            1

        """
        if self.multiple_edges(None):
            return self._cg.all_arcs(u, v)
        if not self.has_edge(u, v, None):
            raise RuntimeError("%s, %s not an edge of the graph."%(u,v))
        l = self._cg.arc_label(u, v)
        if l == 0: return None
        return l

    def has_edge(self, u, v, l):
        """
        Returns whether this graph has edge u,v.

        EXAMPLE:
            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: D.has_edge(0,1,None)
            True

        """
        if l is None:
            return self._cg.has_arc(u, v)
        else:
            return self._cg.has_arc_label(u, v, l)

    def iterator_edges(self, vertices, labels, not_directed):
        """
        Iterate over the edges incident to a sequence of vertices.

        INPUT:
            vertices:     a list of vertex labels
            labels:       boolean
            not_directed: boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        DOCTEST:
            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.iterator_edges([],True,True)
            <listiterator object at ...>
        """
        if not_directed and labels:
            return iter([(v,u,l) for v in vertices for u in self._cg.out_neighbors(v) for l in self._cg.all_arcs(v, u)])
        elif not_directed:
            return iter([(v,u) for v in vertices for u in self._cg.out_neighbors(v) for l in self._cg.all_arcs(v, u)])
        elif labels:
            return iter([(v,u,l) for v in vertices for u in self._cg.out_neighbors(v) if u >= v or u not in vertices for l in self._cg.all_arcs(v, u)])
        else:
            return iter([(v,u) for v in vertices for u in self._cg.out_neighbors(v) if u >= v or u not in vertices for l in self._cg.all_arcs(v, u)])

    def iterator_in_edges(self, vertices, labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:
            vertices:     a list of vertex labels
            labels:       boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        DOCTEST:
            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.iterator_in_edges([],True)
            <listiterator object at ...>
        """
        if self.multiple_edges(None):
            if labels:
                return iter([(u,v,l) for v in vertices for u in self._cg.in_neighbors(v) for l in self.get_edge_label(u, v)])
            else:
                return iter([(u,v) for v in vertices for u in self._cg.in_neighbors(v) for l in self.get_edge_label(u, v)])
        else:
            if labels:
                return iter([(u,v,self.get_edge_label(u, v)) for v in vertices for u in self._cg.in_neighbors(v)])
            else:
                return iter([(u,v) for v in vertices for u in self._cg.in_neighbors(v)])

    def iterator_out_edges(self, vertices, labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:
            vertices:     a list of vertex labels
            labels:       boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        DOCTEST:
            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.iterator_out_edges([],True)
            <listiterator object at ...>
        """
        if self.multiple_edges(None):
            if labels:
                return iter([(v,u,l) for v in vertices for u in self._cg.out_neighbors(v) for l in self.get_edge_label(v, u)])
            else:
                return iter([(v,u) for v in vertices for u in self._cg.out_neighbors(v) for l in self.get_edge_label(v, u)])
        else:
            if labels:
                return iter([(v,u,self.get_edge_label(v, u)) for v in vertices for u in self._cg.out_neighbors(v)])
            else:
                return iter([(v,u) for v in vertices for u in self._cg.out_neighbors(v)])

    def multiple_edges(self, new):
        """
        Get/set whether or not self allows multiple edges.

        INPUT:
            new: boolean or None

        DOCTEST:
            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.multiple_edges(True)
            sage: G.multiple_edges(None)
            True
        """
        if new is None:
            return self._multiple_edges
        if new:
            self._multiple_edges = True
        else:
            self._multiple_edges = False

    def set_edge_label(self, u, v, l, directed):
        """
        Label the edge (u,v) by l.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.set_edge_label(1,2,'a',True)
        """
        if not self.has_edge(u, v, None):
            return
        if self.multiple_edges(None):
            if len(self.get_edge_label(u, v)) > 1:
                raise RuntimeError("Cannot set edge label, since there are multiple edges from %s to %s."%(u,v))
        # now we know there is exactly one edge from u to v
        if directed:
            ll = self.get_edge_label(u,v)
            if ll is not None:
                self._cg.del_arc_label(u, v, ll)
            self._cg.add_arc_label(u, v, l)
        else:
            ll = self.get_edge_label(u,v)
            if ll is not None:
                self._cg.del_arc_label(u, v, ll)
                self._cg.del_arc_label(v, u, ll)
            self._cg.add_arc_label(u, v, l)
            self._cg.add_arc_label(v, u, l)



