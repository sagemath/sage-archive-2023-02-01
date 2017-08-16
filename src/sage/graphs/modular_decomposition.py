# -*- coding: utf-8 -*-
r"""
Modular Decomposition

This module implements the function for computing the modular decomposition
of undirected graphs.


#*****************************************************************************
#       Copyright (C) 2017 Lokesh Jain <lokeshj1703@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""

from sage.graphs.graph import Graph
from collections import deque

PRIME = 0
SERIES = 1
PARALLEL = 2
NORMAL = 3
FOREST = -1
LEFT_SPLIT = 1
RIGHT_SPLIT = 2
BOTH_SPLIT = 3
NO_SPLIT = 0
LEFT_OF_SOURCE = -1
RIGHT_OF_SOURCE = 1
SOURCE = 0

class NodeInfo:
    """
    Node class stores information about the node type, node split and index 
    of the node in the parent tree. 

    Node type can be PRIME, SERIES, PARALLEL, NORMAL or FOREST. Node split can 
    be NO_SPLIT, LEFT_SPLIT, RIGHT_SPLIT or BOTH_SPLIT. A node is split in the 
    refinement phase and the split used is propagated to the ancestors.

    - ``node_type`` -- Specifies the type of node

        * ``"PARALLEL"`` -- indicates the node is a parallel module

        * ``"SERIES"`` -- indicates the node is a series module

        * ``"PRIME"`` -- indicates the node is a prime module

        * ``"FOREST"`` -- indicates a forest containing trees

        * ``"NORMAL"`` -- indicates the node is normal containing a vertex

    - ``node_split`` -- Specifies the type of splits which have occurred in 
                        the node and its descendants

        * ``"LEFT_SPLIT"`` -- indicates a left split has occurred

        * ``"RIGHT_SPLIT"`` -- indicates a right split has occurred

        * ``"BOTH_SPLIT"`` -- indicates both left and right split have occurred

        * ``"NO_SPLIT"`` -- indicates no split has occurred

    - ``index_in_root`` -- specifies the index of the node in the forest 
                           obtained after promotion phase

    - ``comp_num`` -- specifies the number given to nodes in a (co)component 
                      before refinement

    - ``is_separated`` -- specifies whether a split has occurred with the node 
                          as the root

    """

    def __init__(self, node_type):
        self.node_type = node_type
        self.node_split = NO_SPLIT
        self.index_in_root = -1
        self.comp_num = -1
        self.is_separated = False

    def set_node_split(self, node_split):
        """
        Add node_split to the node split of self. 

        LEFT_SPLIT and RIGHT_SPLIT can exist together in self as BOTH_SPLIT.

        INPUT:

        - ``node_split`` -- node_split to be added to self

        """
        if self.node_split == NO_SPLIT:
            self.node_split = node_split
        elif ((self.node_split == LEFT_SPLIT and
                       node_split == RIGHT_SPLIT) or
                  (self.node_split == RIGHT_SPLIT and
                           node_split == LEFT_SPLIT)):
            self.node_split = BOTH_SPLIT

    def has_left_split(self):
        """
        Return true if self has LEFT_SPLIT

        """
        return self.node_split == LEFT_SPLIT or self.node_split == BOTH_SPLIT

    def has_right_split(self):
        """
        Return true if self has RIGHT_SPLIT

        """
        return self.node_split == RIGHT_SPLIT or self.node_split == BOTH_SPLIT

    def __str__(self):
        if self.node_type == SERIES:
            return "SERIES"
        elif self.node_type == PARALLEL:
            return "PARALLEL"
        elif self.node_type == PRIME:
            return "PRIME"
        elif self.node_type == FOREST:
            return "FOREST"
        else:
            return "NORMAL"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return self.node_type == other.node_type


def modular_decomposition(graph):
    """
    Compute the modular decomposition tree for the input graph.

    The tree structure is represented in form of nested lists. A tree node is 
    a list with two elements. The first element is object of class NodeInfo 
    and second element is a list which contains other tree nodes.

    INPUT:

    - ``graph`` -- The graph for which modular decompostion
                   tree needs to be computed

    OUTPUT:

    A nested list representing the modular decomposition tree computed
    for the graph

    EXAMPLES:

    The Icosahedral graph is Prime::

        sage: from sage.graphs.modular_decomposition import \
              modular_decomposition, test_modular_decomposition, print_md_tree
        sage: print_md_tree(modular_decomposition(graphs.IcosahedralGraph()))
        PRIME
              8
              5
              1
              11
              7
              0
              6
              9
              2
              4
              10
              3

    The Octahedral graph is not Prime::

        sage: print_md_tree(modular_decomposition(graphs.OctahedralGraph()))
        SERIES
              PARALLEL
                2
                3
              PARALLEL
                1
                4
              PARALLEL
                0
                5

    Tetrahedral Graph is Series::

        sage: print_md_tree(modular_decomposition(graphs.TetrahedralGraph()))
        SERIES
              3
              2
              1
              0

    Modular Decomposition tree containing combination of parallel and series modules::

        sage: d = {2:[4,3,5], 1:[4,3,5], 5:[3,2,1,4], 3:[1,2,5], 4:[1,2,5]}
        sage: g = Graph(d)
        sage: print_md_tree(modular_decomposition(g))
        SERIES
              5
              PARALLEL
                3
                4
              PARALLEL
                1
                2

    TESTS:

    Bad Input::

        sage: g = DiGraph()
        sage: modular_decomposition(g)
        Traceback (most recent call last):
        ...
        ValueError: Graph must be undirected

    Empty Graph is Prime::

        sage: g = Graph()
        sage: modular_decomposition(g)
        [PRIME, []]

    Graph from Marc Tedder implementation of modular decomposition::

        sage: d = {1:[5,4,3,24,6,7,8,9,2,10,11,12,13,14,16,17], 2:[1], \
                    3:[24,9,1], 4:[5,24,9,1], 5:[4,24,9,1], 6:[7,8,9,1], \
                    7:[6,8,9,1], 8:[6,7,9,1], 9:[6,7,8,5,4,3,1], 10:[1], \
                    11:[12,1], 12:[11,1], 13:[14,16,17,1], 14:[13,17,1], \
                    16:[13,17,1], 17:[13,14,16,18,1], 18:[17], 24:[5,4,3,1]}
        sage: g = Graph(d)
        sage: test_modular_decomposition(modular_decomposition(g), g)
        True

    Graph from wikipedia link :wikipedia:`Modular_decomposition`::

        sage: d2 = {1:[2,3,4], 2:[1,4,5,6,7], 3:[1,4,5,6,7], 4:[1,2,3,5,6,7], \
                    5:[2,3,4,6,7], 6:[2,3,4,5,8,9,10,11], \
                    7:[2,3,4,5,8,9,10,11], 8:[6,7,9,10,11], 9:[6,7,8,10,11], \
                    10:[6,7,8,9], 11:[6,7,8,9]}
        sage: g = Graph(d)
        sage: test_modular_decomposition(modular_decomposition(g), g)
        True

    Tetrahedral Graph is Series::

        sage: print_md_tree(modular_decomposition(graphs.TetrahedralGraph()))
        SERIES
              3
              2
              1
              0

    Modular Decomposition tree containing combination of parallel and series modules::

        sage: d = {2:[4,3,5], 1:[4,3,5], 5:[3,2,1,4], 3:[1,2,5], 4:[1,2,5]}
        sage: g = Graph(d)
        sage: print_md_tree(modular_decomposition(g))
        SERIES
              5
              PARALLEL
                3
                4
              PARALLEL
                1
                2

    """
    if graph.is_directed():
        raise ValueError("Graph must be undirected")

    if graph.order() == 0:  #Empty Graph
        return create_prime_node()

    if graph.order() == 1:  # Single vertex graph
        root = create_normal_node(next(graph.vertex_iterator()))
        return root

    if not graph.is_connected():

        # Parallel case:- The tree contains the MD trees of its connected
        # components as subtrees
        components = graph.connected_components()
        root = create_parallel_node()
        for component in components:
            root[1].append(modular_decomposition(graph.subgraph(component)))
        return root
    elif graph.complement().is_connected():     #Prime Graph
        root = create_prime_node()
    else:
        root = create_series_node()     #Series Graph

    bfs_generator = graph.breadth_first_search(next(graph.vertex_iterator()),
                                               report_distance=True)

    prev_level_distance = -1  # used as a demarker for different levels in bfs
    prev_level_list = []  # stores the vertices in previous level

    # dictionary stores the distance of vertices from the SOURCE
    vertex_dist = {}

    # dictionary stores the position of vertices w.r.t SOURCE
    vertex_status = {}
    vertex_status[next(graph.vertex_iterator())] = SOURCE

    # Different distances from the source vertex are considered
    # as different levels in the algorithm
    for (vertex, distance) in bfs_generator:
        vertex_dist[vertex] = distance

        # Mark the neighbours of source as LEFT_OF_SOURCE as they appear to
        # left of source in the forest, other vertices are marked as
        # RIGHT_OF_SOURCE
        if distance == 1:
            vertex_status[vertex] = LEFT_OF_SOURCE
        elif distance != 0:
            vertex_status[vertex] = RIGHT_OF_SOURCE

        if distance != prev_level_distance:  # On start of new level in BFS
            if prev_level_list:
                # MD Tree is computed for each level and added to the forest
                root[1].append(modular_decomposition(
                                   graph.subgraph(prev_level_list))
                              )
            prev_level_list = []
            prev_level_distance = distance
        prev_level_list.append(vertex)

    # The last level is left out in the above loop
    root[1].append(modular_decomposition(graph.subgraph(prev_level_list)))

    # The MD tree for the neighbours of source marked as LEFT_OF_SOURCE
    # are placed left of Source in the forest. root[1][1] is the source
    # and root[1][0] is the MD tree for the neighbours therefore the
    # the first two elements in the list are replaced
    root[1][0], root[1][1] = root[1][1], root[1][0]

    root[0].node_type = FOREST
    clear_node_split_info(root)
    number_cocomponents(root, vertex_status)
    number_components(root, vertex_status)
    refine(graph, root, vertex_dist, vertex_status)
    promote_left(root)
    promote_right(root)
    promote_child(root)
    assembly(graph, root, vertex_status, vertex_dist)

    if root[0].node_type == FOREST:
        return root[1][0]
    else:
        return root

def number_components(root, vertex_status):
    """
    Function to number the components to the right of SOURCE vertex in the
    forest input to the assembly phase

    INPUT:

    - ``root`` -- the forest which contains the components and cocomponents
    - ``vertex_status`` -- dictionary which stores the position of vertex
                           w.r.t SOURCE

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, number_components
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        create_normal_node(3), create_normal_node(1), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], [NodeInfo(PARALLEL), \
                        [create_normal_node(6), create_normal_node(7)]]]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: RIGHT_OF_SOURCE, \
                               5: RIGHT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: number_components(forest, vertex_status)
        sage: forest[1][-1][1][0][0].comp_num
        2
        sage: forest[1][-1][1][1][0].comp_num
        3

    TESTS:

        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        create_normal_node(3), create_normal_node(1), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], [NodeInfo(PARALLEL), \
                        [create_normal_node(6), create_normal_node(7)]]]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: RIGHT_OF_SOURCE, \
                               5: RIGHT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: number_components(forest, vertex_status)
        sage: forest[1][-1][1][0][0].comp_num == 2 and \
              forest[1][-1][1][1][0].comp_num == 3
        True
        sage: forest[1][-2][1][0][0].comp_num == 1 and \
              forest[1][-2][1][1][0].comp_num == 1
        True

    """
    comp_num = 0
    flag = False

    if not root:    #root is empty
        return ValueError("Input forest {} is empty".format(root))

    for tree in root[1]:

        #flag set to True after source vertex is encountered
        if tree[0].node_type == NORMAL and \
                        vertex_status[tree[1][0]] == SOURCE:
            flag = True
            continue

        if not flag:  # Cocomponents are skipped
            continue

        comp_num += recursively_number_cocomponents(tree, comp_num, PARALLEL)

def number_cocomponents(root, vertex_status):
    """
    Function to number the cocomponents to the left of SOURCE vertex in the
    forest input to the assembly phase

    INPUT:

    - ``root`` -- the forest which contains the cocomponents and components
    - ``vertex_status`` -- dictionary which stores the position of vertex
                           w.r.t SOURCE

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, number_cocomponents
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], [NodeInfo(PARALLEL), \
                        [create_normal_node(6), create_normal_node(7)]], \
                        create_normal_node(3), create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: LEFT_OF_SOURCE, \
                               7: LEFT_OF_SOURCE}
        sage: number_cocomponents(forest, vertex_status)
        sage: forest[1][1][1][0][0].comp_num
        1
        sage: forest[1][1][1][1][0].comp_num
        2

    TESTS:

        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], [NodeInfo(PARALLEL), \
                        [create_normal_node(6), create_normal_node(7)]], \
                        create_normal_node(3), create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: LEFT_OF_SOURCE, \
                               7: LEFT_OF_SOURCE}
        sage: number_cocomponents(forest, vertex_status)
        sage: forest[1][1][1][0][0].comp_num == 1 and \
              forest[1][1][1][1][0].comp_num == 2
        True
        sage: forest[1][2][1][0][0].comp_num == 3 and \
              forest[1][2][1][1][0].comp_num == 3
        True

    """
    cocomp_num = 0
    for tree in root[1]:
        # Only cocomponents are numbered
        if tree[0].node_type == NORMAL and \
                        vertex_status[tree[1][0]] == SOURCE:
            break
        cocomp_num += recursively_number_cocomponents(tree, cocomp_num,
                                                      SERIES)


def recursively_number_cocomponents(tree, cocomp_num, by_type):
    """
    Recursively number the nodes in the (co)components. 

    If the tree node_type is same as by_type then cocomp_num is incremented 
    before assigning to the subtree else entire tree is numbered by cocomp_num

    INPUT:

    - ``tree`` -- the forest which contains the cocomponents and components
    - ``cocomp_num`` -- input number to be used as reference for numbering
                        the (co)components
    - ``by_type`` -- type which determines how numbering is done

    OUTPUT:

    The value incremented to cocomp_num

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, recursively_number_cocomponents
        sage: tree = [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]]
        sage: recursively_number_cocomponents(tree, 1, SERIES)
        2
        sage: tree[0].comp_num 
        1
        sage: tree[1][0][0].comp_num
        1
        sage: tree[1][1][0].comp_num
        2

    TESTS:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, recursively_number_cocomponents
        sage: tree = [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]]
        sage: recursively_number_cocomponents(tree, 1, SERIES)
        2
        sage: tree[0].comp_num == 1 and tree[1][0][0].comp_num == 1 and tree[1][1][0].comp_num == 2
        True

    """
    
    # inner function
    def number_subtree(tree, number):
        """
        set the ``comp_num`` for all the nodes in the subtree to ``number``

        INPUT:

        - ``tree`` -- tree to be numbered
        - ``number`` -- number assigned to the tree

        """
        tree[0].comp_num = number
        if tree[0].node_type != NORMAL:
            for subtree in tree[1]:
                number_subtree(subtree, number)

    orig_cocomp_num = cocomp_num

    if tree[0].node_type==by_type:
        # if node_type is same as tree's node_type then cocomp_num is
        # incremented before assigning to each subtree
        tree[0].comp_num = cocomp_num
        for subtree in tree[1]:
            number_subtree(subtree, cocomp_num)
            cocomp_num += 1
    else:
        # entire tree is numbered by cocomp_num
        number_subtree(tree, cocomp_num)
        cocomp_num+=1
    return cocomp_num - orig_cocomp_num

def assembly(graph, root, vertex_status, vertex_dist):
    """
    Assemble the forest obtained after the promotion phase into a modular 
    decomposition tree.

    INPUT:

    - ``graph`` -- graph whose MD tree is to be computed
    - ``root`` -- Forest which would be assembled into a MD tree
    - ``vertex_status`` -- Dictionary which stores the position of
                           vertex with respect to the source

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, number_cocomponents, number_components, \
              assembly
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: forest[1][0][0].comp_num = 1
        sage: forest[1][1][0].comp_num = 1
        sage: forest[1][1][1][0][0].comp_num = 1
        sage: forest[1][1][1][1][0].comp_num = 1
        sage: number_components(forest, vertex_status)
        sage: assembly(g, forest, vertex_status, vertex_dist)
        sage: forest[1]
        [[PRIME, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [NORMAL, [3]], [PARALLEL, [[NORMAL, [6]], [NORMAL, [7]]]], [NORMAL, [1]]]]]

        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: number_cocomponents(forest, vertex_status)
        sage: assembly(g, forest, vertex_status, vertex_dist)
        sage: forest[1]
        [[PRIME, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]], [NORMAL, [3]]]], [PARALLEL, [[NORMAL, [6]], [NORMAL, [7]]]], [NORMAL, [1]]]]]

    TESTS:

        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: number_cocomponents(forest, vertex_status)
        sage: number_components(forest, vertex_status)
        sage: assembly(g, forest, vertex_status, vertex_dist)
        sage: forest[1]
        [[PRIME, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]], [NORMAL, [3]]]], [PARALLEL, [[NORMAL, [6]], [NORMAL, [7]]]], [NORMAL, [1]]]]]

    """

    # Maps index to the mu (co)component computed for the tree at the index
    mu = {}

    # Stores index of tree containing the source vertex in the forest
    source_index = -1

    # Maps index to list of vertices in the tree at the index in the forest
    vertices_in_component = {}

    # comp_num of parent should be equal to comp_num of its first child
    update_comp_num(root)

    for index, component in enumerate(root[1]):

        if component[0].node_type == NORMAL and \
                        vertex_status[component[1][0]] == SOURCE:
            source_index = root[1].index(component)

        vertices_in_component[index] = get_vertices(component)
        component[0].index_in_root = index

    # compute mu values for (co)components
    for index, component in enumerate(root[1]):
        if index < source_index:
            mu[index] = compute_mu_for_co_component(graph, index,
                                                    source_index, root,
                                                    vertices_in_component)
        elif index > source_index:
            mu[index] = compute_mu_for_component(graph, index,
                                                 source_index, root,
                                                 vertices_in_component)

    mu[source_index] = root[1][source_index]

    # stores the leftmost cocomponent included in the module containing
    # source_index

    left = root[1][source_index]

    # stores the rightmost component included in the module containing
    # source_index
    right = root[1][source_index]

    while len(root[1]) != 1:
        # source_index is changed everytime a new module is formed therefore
        # updated left or right are changed every time module is formed.

        # First series module is attempted
        result, source_index = check_series(root, left, right,
                                              source_index, mu)
        if result:
            left = root[1][source_index][1][0]
            continue

        # If series module cant be formed, parallel is tried
        result, source_index = check_parallel(graph, root, left, right,
                                                source_index, mu, vertex_dist,
                                                vertices_in_component)
        if result:
            right = root[1][source_index][1][-1]
            continue

        # Finally a prime module is formed if both
        # series and parallel can not be created
        result, source_index = check_prime(graph, root, left, right,
                                             source_index, mu, vertex_dist,
                                             vertices_in_component)
        if result:
            if root[1][source_index][1][0][0].index_in_root != -1:
                left = root[1][source_index][1][0]
            if root[1][source_index][1][-1][0].index_in_root != -1:
                right = root[1][source_index][1][-1]


def update_comp_num(root):
    """
    Set the comp_num of the parent to the comp_num of its first child

    INPUT:

    - ``root`` -- root of the tree whose nodes comp_num needs to be updated

    """
    if root[0].node_type != NORMAL:
        root[0].comp_num = root[1][0][0].comp_num
        for child in root[1]:
            update_comp_num(child)


def check_prime(graph, root, left, right,
                source_index, mu, vertex_dist,
                vertices_in_component):
    """
    Assemble the forest to create a prime module.

    INPUT:

    - ``root`` - forest which needs to be assembled
    - ``left`` - The leftmost fragment of the last module
    - ``right`` - The rightmost fragment of the last module
    - ``source_index`` - index of the tree containing the source vertex
    - ``mu`` - dictionary which maps the (co)components with their mu values.

    OUTPUT:

    ``[module_formed, source_index]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``source_index`` is the index of the
    new module which contains the source vertex

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, number_cocomponents, number_components, \
              check_prime, get_vertices, compute_mu_for_co_component, \
              compute_mu_for_component
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: source_index = 2
        sage: vertices_in_component = {}
        sage: mu = {}
        sage: left = right = forest[1][2]
        sage: for index, component in enumerate(forest[1]):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component[0].index_in_root = index
        sage: for index, component in enumerate(forest[1]):
        ....:     if index < source_index:
        ....:         mu[index] = compute_mu_for_co_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        ....:     elif index > source_index:
        ....:         mu[index] = compute_mu_for_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        sage: forest[1][0][0].comp_num = 1
        sage: forest[1][1][0].comp_num = 1
        sage: forest[1][1][1][0][0].comp_num = 1
        sage: forest[1][1][1][1][0].comp_num = 1
        sage: number_components(forest, vertex_status)
        sage: check_prime(g, forest, left, right,
        ....:              source_index, mu, vertex_dist,
        ....:              vertices_in_component)
        [True, 0]
        sage: forest[1]
        [[PRIME, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [NORMAL, [3]], [PARALLEL, [[NORMAL, [6]], [NORMAL, [7]]]], [NORMAL, [1]]]]]

    """
    # stores the index of rightmost component included in the prime module
    new_right_index = source_index + 1 if source_index + 1 < len(root[1]) \
                                       else source_index

    # stores the index of leftmost component included in the prime module
    new_left_index = source_index - 1 if source_index - 1 >= 0 \
                                      else source_index

    # stores the indices of the cocomponents included in the prime module
    # the cocomponents are extracted one by one for adding more components
    left_queue = deque()

    # stores the indices of the components included in the prime module
    # the components are extracted one by one for adding more cocomponents
    right_queue = deque()

    if new_left_index != source_index:
        left_queue.append(new_left_index)
    if new_right_index != source_index:
        right_queue.append(new_right_index)

    while left_queue or right_queue:

        if left_queue:

            # cocomponent indices extracted from the queue
            left_index = left_queue.popleft()

            # more components added based on the below condition
            while new_right_index < len(root[1]) - 1 and \
                            root[1][new_right_index][0].index_in_root < \
                            mu[left_index][0].index_in_root:
                new_right_index += 1
                right_queue.append(new_right_index)

            # cocomponent added while cocomponent at left_index
            # has cocomponent to its left with same comp_num
            while has_left_cocomponent_fragment(root, left_index):
                if left_index >= 1:
                    left_index -= 1
                    if new_left_index > left_index:
                        left_queue.append(left_index)
                    new_left_index = min(left_index, new_left_index)

        if right_queue:

            # component indices extracted from the queue
            right_index = right_queue.popleft()

            # more cocomponents added based on the below condition
            while new_left_index > 0 and \
                            root[1][new_left_index][0].index_in_root > \
                            mu[right_index][0].index_in_root:
                new_left_index -= 1
                left_queue.append(new_left_index)

            # component is added while component at right_index
            # has component to its right with same comp_num
            # or has a connected component with vertices at different
            # level from the source vertex
            while has_right_component_fragment(root, right_index) or \
                    has_right_layer_neighbor(graph, root,
                                             right_index, vertex_dist,
                                             vertices_in_component):

                if has_right_layer_neighbor(graph, root,
                                            right_index, vertex_dist,
                                            vertices_in_component):
                    new_left_index = 0
                    new_right_index = len(root[1]) - 1
                    break

                if right_index + 1 < len(root[1]):
                    right_index += 1
                    if new_right_index < right_index:
                        right_queue.append(right_index)
                    new_right_index = max(right_index, new_right_index)

    node = create_prime_node()

    # vertices or modules are added in the prime_module
    for temp in range(new_left_index, new_right_index + 1):
        node[1].append(root[1][temp])

    # list elements included in the prime module
    # are removed from the forest
    root[1][new_left_index:new_right_index + 1] = []

    #insert the newly created prime module in the forest
    root[1].insert(new_left_index, node)

    return [True, new_left_index]


def check_parallel(graph, root, left, right,
                   source_index, mu, vertex_dist,
                   vertices_in_component):
    """
    Assemble the forest to create a parallel module.

    INPUT:

    - ``root`` -- forest which needs to be assembled
    - ``left`` -- The leftmost fragment of the last module
    - ``right`` -- The rightmost fragment of the last module
    - ``source_index`` -- index of the tree containing the source vertex
    - ``mu`` -- dictionary which maps the (co)components with their mu values.

    OUTPUT:

    ``[module_formed, source_index]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``source_index`` is the index of the
    new module which contains the source vertex

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, number_cocomponents, number_components, \
              check_parallel, get_vertices, compute_mu_for_co_component, \
              compute_mu_for_component
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(4, 1)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7), create_normal_node(1)]]]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 2}
        sage: source_index = 2
        sage: vertices_in_component = {}
        sage: mu = {}
        sage: left = right = forest[1][2]
        sage: for index, component in enumerate(forest[1]):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component[0].index_in_root = index
        sage: for index, component in enumerate(forest[1]):
        ....:     if index < source_index:
        ....:         mu[index] = compute_mu_for_co_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        ....:     elif index > source_index:
        ....:         mu[index] = compute_mu_for_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        sage: number_components(forest, vertex_status)
        sage: check_parallel(g, forest, left, right,
        ....:              source_index, mu, vertex_dist,
        ....:              vertices_in_component)
        [True, 2]
        sage: forest[1]
        [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [PARALLEL, [[NORMAL, [3]], [NORMAL, [6]], [NORMAL, [7]], [NORMAL, [1]]]]]

    """

    # stores the index of rightmost component included in the parallel module
    new_right_index = source_index

    while new_right_index + 1 < len(root[1]):

        # component at new_right_index + 1 is added only if it doesn't have
        # a component to its right with same comp_num
        if has_right_component_fragment(root, new_right_index + 1):
            break

        # component at new_right_index + 1 is added only if it doesn't have a
        # connected component to its right with vertices at different level
        # from its vertices
        if has_right_layer_neighbor(graph, root, new_right_index + 1,
                                    vertex_dist, vertices_in_component):
            break

        # stores the index in root of new component to be added in the
        # parallel module
        i = root[1][new_right_index + 1][0].index_in_root

        # condition for adding more components in the parallel module
        if mu[i][0].index_in_root >= left[0].index_in_root:
            new_right_index += 1
        else:
            break

    # if new_right_index > source_index then only parallel
    # module can be formed
    if source_index != new_right_index:
        node = create_parallel_node()
        temp = source_index
        for temp in range(source_index, new_right_index + 1):

            # if module X to be included in the new parallel module Y is also
            # parallel then children of X and not X are included in Y
            if root[1][temp][0].node_type == PARALLEL:
                for tree in root[1][temp][1]:
                    node[1].append(tree)
                    tree[0].index_in_root = root[1][temp][0].index_in_root
            else:
                node[1].append(root[1][temp])

        # list elements included in the parallel module are removed from the
        # forest
        root[1][source_index:new_right_index + 1] = []

        # insert the newly created parallel module into the forest
        root[1].insert(source_index, node)

        return [True, source_index]

    # no parallel module was formed
    return [False, source_index]


def check_series(root, left, right, source_index, mu):
    """
    Assemble the forest to create a series module.

    - ``root`` -- forest which needs to be assembled
    - ``left`` -- The leftmost fragment of the last module
    - ``right`` -- The rightmost fragment of the last module
    - ``source_index`` -- index of the tree containing the source vertex
    - ``mu`` -- dictionary which maps the (co)components with their mu values.

    OUTPUT:

    ``[module_formed, source_index]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``source_index`` is the index of the
    new module which contains the source vertex

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, number_cocomponents, number_components, \
              check_series, get_vertices, compute_mu_for_co_component, \
              compute_mu_for_component
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: source_index = 2
        sage: vertices_in_component = {}
        sage: mu = {}
        sage: left = right = forest[1][2]
        sage: for index, component in enumerate(forest[1]):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component[0].index_in_root = index
        sage: for index, component in enumerate(forest[1]):
        ....:     if index < source_index:
        ....:         mu[index] = compute_mu_for_co_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        ....:     elif index > source_index:
        ....:         mu[index] = compute_mu_for_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        sage: number_cocomponents(forest, vertex_status)
        sage: number_components(forest, vertex_status)
        sage: check_series(forest, left, right,
        ....:              source_index, mu)
        [True, 1]
        sage: forest[1]
        [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]], [NORMAL, [3]]]], [PARALLEL, [[NORMAL, [6]], [NORMAL, [7]]]], [NORMAL, [1]]]

    """

    # stores the index of leftmost component included in the parallel module
    new_left_index = source_index

    while new_left_index > 0:

        # cocomponent at new_left_index - 1 is added only if it doesn't have
        # a cocomponent to its left with same comp_num
        if has_left_cocomponent_fragment(root, new_left_index - 1):
            break

        # stores the index in root of new cocomponent to be added in the
        # series module
        i = root[1][new_left_index - 1][0].index_in_root

        # condition for adding more cocomponents in the series module
        if mu[i][0].index_in_root <= right[0].index_in_root:
            new_left_index -= 1
        else:
            break

    # if new_left_index < source_index then only series module can be formed
    if source_index != new_left_index:
        node = create_series_node()
        for temp in range(new_left_index, source_index + 1):

            if root[1][temp][0].node_type == SERIES:
                # if module X to be included in the new series module Y is
                # also series then children of X and not X are included in Y
                for tree in root[1][temp][1]:
                    tree[0].index_in_root = root[1][temp][0].index_in_root
                    node[1].append(tree)
            else:
                node[1].append(root[1][temp])

        # list elements included in the series module
        # are removed from the forest
        root[1][new_left_index:source_index + 1] = []

        # insert the newly created series module into the forest
        root[1].insert(new_left_index, node)

        return [True, new_left_index]

    # no series module could be formed
    return [False, new_left_index]


def has_left_cocomponent_fragment(root, cocomp_index):
    """
    Return True if cocomponent at cocomp_index has a cocomponent to its left 
    with same comp_num

    INPUT:

    - ``root`` -- The forest to which cocomponent belongs
    - ``cocomp_index`` -- Index at which cocomponent is present in root

    OUTPUT:

    ``True`` if cocomponent at  cocomp_index has a cocomponent
    to its left with same comp_num else ``False``

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, create_normal_node, \
              has_left_cocomponent_fragment
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: forest[1][0][0].comp_num = 1
        sage: forest[1][1][0].comp_num = 1
        sage: forest[1][1][1][0][0].comp_num = 1
        sage: forest[1][1][1][1][0].comp_num = 1
        sage: has_left_cocomponent_fragment(forest, 1)
        True
        sage: has_left_cocomponent_fragment(forest, 0)
        False

    """
    for index in range(cocomp_index):
        if root[1][index][0].comp_num == root[1][cocomp_index][0].comp_num:
            return True
    return False


def has_right_component_fragment(root, comp_index):
    """
    Return True if component at comp_index has a component to its right with 
    same comp_num

    INPUT:

    - ``root`` -- The forest to which component belongs
    - ``comp_index`` -- Index at which component is present in root

    OUTPUT:

    ``True`` if component at  comp_index has a component
    to its right with same comp_num else ``False``

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, create_normal_node, \
              has_right_component_fragment
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: forest[1][3][0].comp_num = 1
        sage: forest[1][4][0].comp_num = 1
        sage: has_right_component_fragment(forest, 3)
        True
        sage: has_right_component_fragment(forest, 4)
        False

    """
    for index in range(comp_index + 1, len(root[1])):
        if root[1][index][0].comp_num == root[1][comp_index][0].comp_num:
            return True
    return False


def has_right_layer_neighbor(graph, root, comp_index,
                             vertex_dist, vertices_in_component):
    """
    Return True if component at comp_index has a connected component to its 
    right with vertices at different level from the source vertex

    INPUT:

    - ``root`` -- The forest to which component belongs
    - ``comp_index`` -- Index at which component is present in root

    OUTPUT:

    ``True`` if component at comp_index has a right layer neighbor
    else ``False``

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, create_normal_node, \
              get_vertices, has_right_layer_neighbor
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest[1]):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component[0].index_in_root = index
        sage: has_right_layer_neighbor(g, forest, 3, vertex_dist, 
        ....:                          vertices_in_component)
        True

    """
    for index in range(comp_index + 1, len(root[1])):

        # check vertex in component at index has different level from vertex
        # in component at comp_index and are connected to each other
        if ((vertex_dist[get_vertex_in(root[1][index])] >
                 vertex_dist[get_vertex_in(root[1][comp_index])]
             ) and
            (is_component_connected(graph, root[1][index][0].index_in_root,
                                    root[1][comp_index][0].index_in_root,
                                    vertices_in_component)
            )):
            return True

    return False


def get_vertex_in(tree):
    """
    Return the first vertex encountered in the depth-first traversal of the 
    tree

    INPUT:

    - ``tree`` -- input modular decomposition tree

    OUTPUT:

    Return the first vertex encountered in recursion

    """
    while tree[0].node_type != NORMAL:
        tree = tree[1][0]
    return tree[1][0]

def compute_mu_for_co_component(graph, component_index, source_index,
                                root, vertices_in_component):
    """
    Compute the mu value for co-component

    INPUT:

    - ``graph`` -- Graph whose MD tree needs to be computed
    - ``component_index`` -- index of the co-component
    - ``source_index`` -- index of the source in the forest
    - ``root`` -- the forest which needs to be assembled into a MD tree
    - ``vertices_in_component`` -- Dictionary which maps index i to list of
                                  vertices in the tree at index i in the forest

    OUTPUT:

    The mu value (component in the forest) for the co-component

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, create_normal_node, \
              get_vertices, compute_mu_for_co_component
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest[1]):
        ....:     vertices_in_component[index] = get_vertices(component)
        sage: compute_mu_for_co_component(g, 0, 2, forest, 
        ....:                             vertices_in_component)
        [NORMAL, [1]]
        sage: compute_mu_for_co_component(g, 1, 2, forest, 
        ....:                             vertices_in_component)
        [NORMAL, [3]]
        
    """

    for index in range(len(root[1]) - 1, source_index, -1):
        if is_component_connected(graph, component_index,
                                  index, vertices_in_component):
            return root[1][index]

    # return the default value
    return root[1][source_index]


def compute_mu_for_component(graph, component_index, source_index,
                             root, vertices_in_component):
    """
    Compute the mu value for component

    INPUT:

    - ``graph`` -- Graph whose MD tree needs to be computed
    - ``component_index`` -- index of the component
    - ``source_index`` -- index of the source in the forest
    - ``root`` -- the forest which needs to be assembled into a MD tree
    - ``vertices_in_component`` -- Dictionary which maps index i to list of
                                   vertices in the tree at the index i in the
                                   forest

    OUTPUT:

    The mu value (co-component in the forest) for the component

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, create_normal_node, \
              get_vertices, compute_mu_for_component
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest[1]):
        ....:     vertices_in_component[index] = get_vertices(component)
        sage: compute_mu_for_component(g, 3, 2, forest, 
        ....:                          vertices_in_component)
        [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]]
        sage: compute_mu_for_component(g, 4, 2, forest, 
        ....:                          vertices_in_component)
        [NORMAL, [2]]

    """

    # default mu value for a component
    mu_for_component = root[1][0]

    for index in range(0, source_index):
        if mu_for_component == root[1][index] and \
                is_component_connected(graph, component_index,
                                       index, vertices_in_component):
            mu_for_component = root[1][index + 1]

    # return the default value
    return mu_for_component


def is_component_connected(graph, index1, index2, vertices_in_component):
    """
    Return True if two (co)components are connected else False

    INPUT:

    - ``graph`` -- Graph whose MD tree needs to be computed
    - ``index1`` -- index of the first (co)component
    - ``index2`` -- index of the second (co)component
    - ``vertices_in_component`` -- Dictionary which maps index i to list of
                                   vertices in the tree at the index i in the
                                   forest

    OUTPUT:

    ``True`` if the (co)components are connected else ``False``

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, create_normal_node, \
              get_vertices, is_component_connected
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest[1]):
        ....:     vertices_in_component[index] = get_vertices(component)
        sage: is_component_connected(g, 0, 1, vertices_in_component)
        False
        sage: is_component_connected(g, 0, 3, vertices_in_component)
        True

    """
    vertices = vertices_in_component[index1]
    index2_vertices_set = set(vertices_in_component[index2])

    for vertex in vertices:
        neighbors = graph.neighbors(vertex)
        if not index2_vertices_set.isdisjoint(neighbors):
            return True
    return False


def get_vertices(component):
    """
    Compute the list of vertices in the (co)component

    INPUT:

    - ``component`` -- (co)component whose vertices need to be returned as a
                       list

    OUTPUT:

    list of vertices in the (co)component
    """
    vertices = []

    # inner recursive function to recurse over the elements in the 
    # ``component``
    def recurse_component(root, vertices):
        if root[0].node_type == NORMAL:
            vertices.append(root[1][0])
            return
        for tree in root[1]:
            recurse_component(tree, vertices)

    recurse_component(component, vertices)
    return vertices

def promote_left(root):
    """
    Perform the promotion phase on the forest root. 

    If child and parent both are marked by LEFT_SPLIT then child is removed 
    and placed just before the parent

    INPUT:

    - ``root`` -- The forest which needs to be promoted

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, maximal_subtrees_with_leaves_in_x, promote_left
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(4, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: x = {u for u in g.neighbor_iterator(2) 
        ....:            if vertex_dist[u] != vertex_dist[2]}
        sage: maximal_subtrees_with_leaves_in_x(forest, 2, x, vertex_status, 
        ....:                                   False, 0)
        sage: promote_left(forest)
        sage: forest
        [FOREST, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [NORMAL, [3]], [PARALLEL, [[NORMAL, [6]]]], [PARALLEL, [[NORMAL, [7]]]], [PARALLEL, []], [NORMAL, [1]]]]

    """
    q = deque()

    # q has [parent, child] elements as parent needs to be modified
    for tree in root[1]:
        q.append([root, tree])

    while q:

        parent, child = q.popleft()

        if child[0].node_type == NORMAL:
            continue

        # stores the elements to be removed from the child
        to_remove = []

        # stores the index of child in parent list
        index = parent[1].index(child)

        for tree in child[1]:

            # if tree and child both have LEFT_SPLIT then tree from
            # child is inserted just before child in the parent
            if tree[0].has_left_split() and child[0].has_left_split():
                parent[1].insert(index, tree)
                index += 1
                to_remove.append(tree)
                q.append([parent, tree])
            else:
                q.append([child, tree])

        for tree in to_remove:
            child[1].remove(tree)


def promote_right(root):
    """
    Perform the promotion phase on the forest root. 

    If child and parent both are marked by RIGHT_SPLIT then child is removed 
    and placed just after the parent

    INPUT:

    - ``root`` -- The forest which needs to be promoted

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, refine, promote_right
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(4, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: refine(g, forest, vertex_dist, vertex_status)
        sage: promote_right(forest)
        sage: forest
        [FOREST, [[NORMAL, [2]], [SERIES, [[SERIES, [[NORMAL, [4]]]], [SERIES, [[NORMAL, [5]]]]]], [NORMAL, [3]], [PARALLEL, []], [PARALLEL, [[NORMAL, [7]]]], [PARALLEL, [[NORMAL, [6]]]], [NORMAL, [1]]]]

    """
    q = deque()

    # q has [parent, child] elements as parent needs to be modified
    for tree in root[1]:
        q.append([root, tree])

    while q:

        parent, child = q.popleft()

        if child[0].node_type == NORMAL:
            continue

        # stores the elements to be removed from the child
        to_remove = []

        # stores the index of child in parent list
        index = parent[1].index(child)

        for tree in child[1]:

            # if tree and child both have RIGHT_SPLIT then tree from
            # child is inserted just after child in the parent
            if tree[0].has_right_split() and child[0].has_right_split():
                parent[1].insert(index + 1, tree)
                to_remove.append(tree)
                q.append([parent, tree])
            else:
                q.append([child, tree])

        for tree in to_remove:
            child[1].remove(tree)


def promote_child(root):
    """
    Perform the promotion phase on the forest `root`. 

    If marked parent has no children it is removed, if it has one child then 
    it is replaced by its child

    INPUT:

    - ``root`` -- The forest which needs to be promoted

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, refine, promote_right, promote_child
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(4, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: refine(g, forest, vertex_dist, vertex_status)
        sage: promote_right(forest)
        sage: promote_child(forest)
        sage: forest
        [FOREST, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [NORMAL, [3]], [NORMAL, [7]], [NORMAL, [6]], [NORMAL, [1]]]]

    """
    q = deque()

    # q has [parent, child] elements as parent needs to be modified
    for tree in root[1]:
        q.append([root, tree])

    while q:

        parent, child = q.popleft()

        if child[0].node_type == NORMAL:
            continue

        # if child node itself has only one child
        if len(child[1]) == 1 and (child[0].node_split != NO_SPLIT or
                                           child[0].node_type == FOREST):
            # replace child node by its own child

            tree = child[1][0]
            index = parent[1].index(child)
            parent[1].insert(index, tree)
            parent[1].remove(child)
            q.append([parent, tree])
        # if child node has no children
        elif ((not child[1]) and child[0].node_split != NO_SPLIT):
            # remove the child node
            parent[1].remove(child)
        else:
            for tree in child[1]:
                q.append([child, tree])


def clear_node_split_info(root):
    """
    Set the node_split of nodes to NO_SPLIT

    INPUT:

    - ``root`` -- The forest which needs to be cleared of split information

    """

    root[0].node_split = NO_SPLIT

    if root[0].node_type != NORMAL:
        for subroot in root[1]:
            clear_node_split_info(subroot)


def refine(graph, root, vertex_dist, vertex_status):
    """
    Refine the forest based on the active edges

    INPUT:

    - ``graph`` -- graph whose MD tree needs to be computed
    - ``root`` -- the forest which needs to be assembled into a MD tree
    - ``vertex_dist`` -- dictionary mapping the vertex with distance from the
                         source
    - ``vertex_status`` -- dictionary mapping the vertex to the position
                           w.r.t source

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, refine
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: refine(g, forest, vertex_dist, vertex_status)
        sage: forest
        [FOREST, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [NORMAL, [3]], [PARALLEL, [[PARALLEL, [[NORMAL, [6]]]], [PARALLEL, [[NORMAL, [7]]]]]], [NORMAL, [1]]]]

    """
    x_used = []

    # active edges of each vertex in the graph is used to refine the forest
    for v in graph.vertices():
        if v in vertex_status and vertex_status[v] == SOURCE:
            continue

        # set of vertices connected through active edges to v
        x = {u for u in graph.neighbor_iterator(v) 
                        if vertex_dist[u] != vertex_dist[v]}

        if x not in x_used:
            x_used.append(x)
            maximal_subtrees_with_leaves_in_x(root, v, x,
                                              vertex_status, False, 0)

    get_child_splits(root)


def get_child_splits(root):
    """
    Add the node_split of children to the parent node

    INPUT:

    - ``root`` -- input modular decomposition tree

    """
    if root[0].node_type != NORMAL:
        for tree in root[1]:
            get_child_splits(tree)
            root[0].set_node_split(tree[0].node_split)


def maximal_subtrees_with_leaves_in_x(root, v, x, vertex_status,
                                      tree_left_of_source, level):
    """
    Refine the forest based on the active edges(x) of vertex v

    INPUT:

    - ``root`` -- the forest which needs to be assembled into a MD tree
    - ``v`` -- the vertex used to refine
    - ``x`` -- set of vertices connected to v and at different distance
               from source compared to v
    - ``vertex_status`` -- dictionary mapping the vertex to the position
                           w.r.t source
    - ``tree_left_of_source`` -- flag indicating whether tree is
    - ``level`` -- indicates the recursion level, 0 for root

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import NodeInfo, \
              FOREST, SERIES, PARALLEL, LEFT_OF_SOURCE, SOURCE, RIGHT_OF_SOURCE, \
              create_normal_node, maximal_subtrees_with_leaves_in_x
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = [NodeInfo(FOREST), [create_normal_node(2), \
                        [NodeInfo(SERIES), [create_normal_node(4), \
                        create_normal_node(5)]], create_normal_node(3), \
                        [NodeInfo(PARALLEL), [create_normal_node(6), \
                        create_normal_node(7)]], create_normal_node(1)]]
        sage: vertex_status = {2: LEFT_OF_SOURCE, 3: SOURCE, \
                               1: RIGHT_OF_SOURCE, 4: LEFT_OF_SOURCE, \
                               5: LEFT_OF_SOURCE, 6: RIGHT_OF_SOURCE, \
                               7: RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: x = {u for u in g.neighbor_iterator(2) 
        ....:            if vertex_dist[u] != vertex_dist[2]}
        sage: maximal_subtrees_with_leaves_in_x(forest, 2, x, vertex_status, 
        ....:                                   False, 0)
        sage: forest
        [FOREST, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [NORMAL, [3]], [PARALLEL, [[NORMAL, [6]], [NORMAL, [7]]]], [NORMAL, [1]]]]
        sage: x = {u for u in g.neighbor_iterator(1) 
        ....:            if vertex_dist[u] != vertex_dist[1]}
        sage: maximal_subtrees_with_leaves_in_x(forest, 1, x, vertex_status, 
        ....:                                   False, 0)
        sage: forest
        [FOREST, [[NORMAL, [2]], [SERIES, [[NORMAL, [4]], [NORMAL, [5]]]], [NORMAL, [3]], [PARALLEL, [[PARALLEL, [[NORMAL, [6]]]], [PARALLEL, [[NORMAL, [7]]]]]], [NORMAL, [1]]]]


    """

    def update_node_info(node, node_type, node_split, comp_num, subtree_list):
        """
        Set the various fields for a tree node and update its subtrees
        
        - ``node`` -- node whose fields need to be updated
        - ``node_type`` -- node_type to be set
        - ``node_split`` -- node_split to be set
        - ``comp_num`` -- comp_num to be set
        - ``subtree_list`` -- list containing the subtrees

        """
        node[0].node_type = node_type
        node[0].node_split = node_split
        node[0].comp_num = comp_num
        node[1] = subtree_list

    return_split = NO_SPLIT     #initialize split to NO_SPLIT

    # all trees in a forest are refined using x
    if root[0].node_type == FOREST:

        # indicates whether tree is left of source, True if left of source
        left_flag = True

        for tree in root[1]:
            if tree[0].node_type == NORMAL and tree[1][0] in vertex_status \
                    and vertex_status[tree[1][0]] == SOURCE:
                left_flag = False
            subtree_result = maximal_subtrees_with_leaves_in_x(tree, v, x,
                                                               vertex_status,
                                                               left_flag,
                                                               level)
            if subtree_result:
                # Mark the ancestors
                root[0].set_node_split(subtree_result[1])

    # handles the prime, series and parallel cases
    elif root[0].node_type != NORMAL:

        flag = True  # indicates the entire root is contained in x
        split_flag = False  # indicates a split is required
        Ta = []  # contains subtrees with leaves in x
        Tb = []  # contains subtrees with leaves not in x

        for subtree in root[1]:

            # refines the children of root
            subtree_result = maximal_subtrees_with_leaves_in_x(subtree, v, x,
                                                               vertex_status,
                                                               tree_left_of_source,
                                                               level + 1)

            if subtree_result:
                flag = flag and subtree_result[0]

                # add the node split of children to root
                root[0].set_node_split(subtree_result[1])

                if subtree_result[0]:
                    Ta.append(subtree)
                    split_flag = True
                else:
                    Tb.append(subtree)

        if root[0].node_type == PRIME:
            # mark all the children of prime nodes
            for prime_subtree in root[1]:
                prime_subtree[0].set_node_split(root[0].node_split)

        if flag:
            # return if all subtrees are in x, no split required
            return [True, root[0].node_split]
        elif split_flag:  # split required`

            split = LEFT_SPLIT

            # if v is right of source and tree is also right of source then
            # RIGHT_SPLIT
            if vertex_status[v] == RIGHT_OF_SOURCE and not tree_left_of_source:
                split = RIGHT_SPLIT

            # add the split to root node_split
            root[0].set_node_split(split)

            if root[0].node_type == PRIME:
                # mark all the children of prime nodes
                for subtree in root[1]:
                    subtree[0].set_node_split(split)
                return [False, split]

            if root[0].is_separated:
                # if root has already been split then further split not
                # required
                return [flag, root[0].node_split]

            node_type = root[0].node_type
            root[0].is_separated = True

            # root[1] would now contain Ta and Tb
            root[1] = []

            # add two nodes for Ta and Tb
            a = create_parallel_node()
            update_node_info(a, node_type, root[0].node_split, 
                             Ta[0][0].comp_num, Ta)
            b = create_parallel_node()
            update_node_info(b, node_type, root[0].node_split, 
                             Tb[0][0].comp_num, Tb)
            root[1].append(a)
            root[1].append(b)

        return_split = root[0].node_split
        return [flag, return_split]
    # root is a vertex and is contained in x
    elif root[1][0] in x:
        return [True, root[0].node_split]
    # root is a vertex and is not contained in x
    else:
        return [False, root[0].node_split]


def create_prime_node():
    """
    Return a prime node with no children

    """

    return [NodeInfo(PRIME), []]


def create_parallel_node():
    """
    Return a parallel node with no children

    """
    return [NodeInfo(PARALLEL), []]


def create_series_node():
    """
    Return a series node with no children

    """
    return [NodeInfo(SERIES), []]


def create_normal_node(vertex):
    """
    Return a normal node with no children

    """
    return [NodeInfo(NORMAL), [vertex]]

def print_md_tree(root):
    """
    Print the modular decomposition tree
    
    - ``root`` -- root of the modular decomposition tree

    """

    def recursive_print_md_tree(root, level):
        """
        Print the modular decomposition tree at root
        
        - ``root`` -- root of the modular decomposition tree 
        - ``level`` -- indicates the depth of root in the original modular 
                       decomposition tree

        """
        if root[0].node_type != NORMAL:
            print("{}{}".format(level,str(root[0])))
            for tree in root[1]:
                recursive_print_md_tree(tree, level + " ")
        else:
            print("{}{}".format(level,str(root[1][0])))

    recursive_print_md_tree(root, "")

#=============================================================================

    # Below functions are implemented to test the modular decomposition tree

#=============================================================================

#Function implemented for testing
def test_modular_decomposition(tree, graph):
    """
    This function tests the input modular decomposition tree using recursion.

    INPUT:

    - ``tree`` -- Modular decomposition tree to be tested
    - ``graph`` -- Graph whose modular decomposition tree needs to be tested

    OUTPUT:

    ``True`` if input ``tree`` is a modular decomposition else ``False``

    EXAMPLES:

        sage: from sage.graphs.modular_decomposition import \
              modular_decomposition, test_modular_decomposition
        sage: g = graphs.HexahedralGraph()
        sage: test_modular_decomposition(modular_decomposition(g), g)
        True

    """
    if tree[0].node_type != NORMAL:
        for module in tree[1]:
            if not test_module(module, graph):
                # test whether modules pass the defining
                # characteristics of modules
                return False
            if not test_modular_decomposition(module,
                                              graph.subgraph(
                                                  get_vertices(module))):
                # recursively test the modular decompostion subtrees
                return False

        if not test_maximal_modules(tree, graph):
            # test whether the mdoules are maximal in nature
            return False

    return True

#Function implemented for testing
def test_maximal_modules(tree, graph):
    """
    This function tests maximal nature of modules in a modular decomposition
    tree. 

    Suppose the module M = [M1, M2, ..., Mn] is the input modular
    decomposition tree. Algorithm forms pairs like (M1, M2), (M1, M3),
    ...(M1, Mn); (M2, M3), (M2, M4), ...(M2, Mn); ... and so on and tries to
    form a module using the pair. If the module formed has same type as M and
    is of type SERIES or PARALLEL then the formed module is not considered
    maximal. Otherwise it is considered maximal and M is not a modular
    decomposition tree.

    INPUT:

    - ``tree`` -- Modular decomposition tree whose modules are tested for
                  maximal nature
    - ``graph`` -- Graph whose modular decomposition tree is tested

    OUTPUT:

    ``True`` if all modules at first level in the modular ddecomposition tree
    are maximal in nature

    """
    if tree[0].node_type != NORMAL:
        for index, module in enumerate(tree[1]):
            for other_index in range(index + 1, len(tree[1])):

                # compute the module formed using modules at index and
                # other_index
                module_formed = form_module(index, other_index, tree, graph)

                if module_formed[0]:
                    # Module formed and the parent of the formed module
                    # should not both be of type SERIES or PARALLEL
                    if ((get_node_type(graph.subgraph(module_formed[1])) ==
                             tree[0].node_type
                         ) and
                        (tree[0].node_type == PARALLEL or
                             tree[0].node_type == SERIES
                        )):
                        continue
                    return False
    return True


#Function implemented for testing
def get_node_type(graph):
    """
    Return the module type of the root of modular decomposition tree for the
    input graph

    INPUT:

    - ``graph`` -- Input sage graph

    OUTPUT:

    ``PRIME`` if graph is PRIME, ``PARALLEL`` if graph is PARALLEL and
    ``SERIES`` if graph is of type SERIES

    """
    if not graph.is_connected():
        return PARALLEL
    elif graph.complement().is_connected():
        return PRIME
    return SERIES


#Function implemented for testing
def form_module(index, other_index, tree, graph):
    """
    This function forms a module out of the modules in the module pair. 

    Let modules input be M1 and M2. Let V be the set of vertices in these
    modules. Suppose x is a neighbor of subset of the vertices in V but not
    all the vertices and x does not belong to V. Then the set of modules also
    include the module which contains x. This process is repeated until a
    module is formed and the formed module if subset of V is returned.

    INPUT:

    - ``index`` -- First module in the module pair
    - ``other_index`` -- Second module in the module pair
    - ``tree`` -- Modular decomposition tree which contains the modules in
                  the module pair
    - ``graph`` -- Graph whose modular decomposition tree is created

    OUTPUT:

    ``[module_formed, vertices]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``vertices`` is a list of verices
    included in the formed module

    """
    vertices = set(get_vertices(tree[1][index]) +
                   get_vertices(tree[1][other_index]))

    # stores all neighbors which are common for all vertices in V
    common_neighbors = set()

    # stores all neighbors of vertices in V which are outside V
    all_neighbors = set()

    while True:
        # remove vertices from all_neighbors and common_neighbors
        all_neighbors = all_neighbors - set(vertices)
        common_neighbors = common_neighbors - set(vertices)

        for v in vertices:
            # stores the neighbors of v which are outside the set of vertices
            neighbor_list = set(graph.neighbors(v))
            neighbor_list = neighbor_list - vertices

            # update all_neighbors and common_neighbors using the
            # neighbor_list
            all_neighbors = all_neighbors | neighbor_list
            common_neighbors = common_neighbors & neighbor_list

        if all_neighbors == common_neighbors:  # indicates a module is formed

            # module formed covers the entire graph
            if len(vertices) == graph.order():
                return [False, vertices]

            return [True, vertices]

        # add modules containing uncommon neighbors into the formed module
        for v in (all_neighbors - common_neighbors):
            for index in range(len(tree[1])):
                if v in get_vertices(tree[1][index]):
                    vertices = vertices | set(get_vertices(tree[1][index]))
                    break

#Function implemented for testing
def test_module(module, graph):
    """
    Test whether input module is actually a module

    INPUT:

    - ``module`` -- Module which needs to be tested
    - ``graph`` -- Input sage graph which contains the module

    OUTPUT:

    ``True`` if input module is a module by definition else ``False``

    """

    # A single vertex is a module
    if module[0].node_type == NORMAL:
        return True

    #vertices contained in module
    vertices_in_module = get_vertices(module)

    #vertices outside module
    vertices_outside = list(set(graph.vertices()) - set(vertices_in_module))

    # Nested module with only one child
    if module[0].node_type != NORMAL and len(module[1]) == 1:
        return False

    # If children of SERIES module are all SERIES modules
    if module[0].node_type == SERIES:
        if children_node_type(module, SERIES):
            return False

    # If children of PARALLEL module are all PARALLEL modules
    if module[0].node_type == PARALLEL:
        if children_node_type(module, PARALLEL):
            return False

    # check the module by definition. Vertices in a module should all either
    # be connected or disconnected to any vertex outside module
    for v in vertices_outside:
        if not either_connected_or_not_connected(v, vertices_in_module,
                                                 graph):
            return False
    return True


#Function implemented for testing
def children_node_type(module, node_type):
    """
    Test whether node_type of children of a node is same as input node_type

    INPUT:

    - ``module`` -- module which is tested
    - ``node_type`` -- input node_type

    OUTPUT:

    ``True`` if node_type of children of module is same as input node_type
    else ``False``

    """
    for tree in module[1]:
        if tree[0].node_type != node_type:
            return False
    return True


#Function implemented for testing
def either_connected_or_not_connected(v, vertices_in_module, graph):
    """
    Test whether v is connected or disconnected to all vertices in the module

    INPUT:

    - ``v`` -- vertex tested
    - ``vertices_in_module`` -- list containing vertices in the module
    - ``graph`` -- graph to which the vertices belong

    OUTPUT:

    ``True`` if v is either connected or disconnected to all the vertices in
    the module else ``False``

    """

    # marks whether vertex v is connected to first vertex in the module
    connected = graph.has_edge(vertices_in_module[0], v)

    # if connected is True then all vertices in module should be connected to
    # v else disconnected
    for u in vertices_in_module:
        if (graph.has_edge(u,v) != connected):
            return False
    return True
