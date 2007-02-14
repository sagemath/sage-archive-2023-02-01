r"""
A module for dealing with lists of graphs.

AUTHOR:
    -- Robert L. Miller (2007-02-10): initial version
    -- Emily A. Kirkman (2007-02-13): added show functions (to_graphics_array and show_graphs)
"""

#*****************************************************************************
#           Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

def from_whatever(data):
    if isinstance(data, file):
        if data.name[data.name.rindex('.'):] == '.g6':
            return from_graph6(data)
        elif data.name[data.name.rindex('.'):] == '.s6':
            return from_sparse6(data)
        else: # convert to list of lines, do each separately
            L = data.readlines()
            return from_whatever(L)
    if isinstance(data, list):
        l = []
        for d in data:
            if isinstance(d, str):
                if len(d) == d.find('\n') + 1 or d.find('\n') == -1:
                    l.append(Graph(d[:d.rfind('\n')]))
                else:
                    l.append(from_whatever(d))
            else:
                l.append(from_whatever(d))
        return l
    if isinstance(data, str):
        data = data.split('\n')
        return from_whatever(data)

def from_graph6(data):
    """
    Returns a list of SAGE Graphs, given a list of graph6 data.

    INPUT:
        data -- can be a string, a list of strings, or a file stream.
    """
    from sage.graphs.graph import Graph
    if isinstance(data,str):
        data = data.split('\n')
        l = []
        for d in data:
            l.append(Graph(d, format = 'graph6'))
        return l
    elif isinstance(data,list):
        l = []
        for d in data:
            if isinstance(d, str):
                if len(d) == d.find('\n') + 1 or d.find('\n') == -1:
                    l.append(Graph(d[:d.rfind('\n')], format='graph6'))
                else:
                    l.append(from_graph6(d))
            else:
                l.append(from_graph6(d))
        return l
    elif isinstance(data,file):
        strlist = data.readlines()
        l = []
        for s in strlist:
            l.append(Graph(s[:s.rfind('\n')], format='graph6'))
        return l

def from_sparse6(data):
    """
    Returns a list of SAGE Graphs, given a list of sparse6 data.

    INPUT:
        data -- can be a string, a list of strings, or a file stream.
    """
    from sage.graphs.graph import Graph
    if isinstance(data,str):
        data = data.split('\n')
        l = []
        for d in data:
            if not d == '':
                l.append(Graph(d, format = 'sparse6'))
        return l
    elif isinstance(data,list):
        l = []
        for d in data:
            if isinstance(d, str):
                if len(d) == d.find('\n') + 1 or d.find('\n') == -1:
                    l.append(Graph(d[:d.rfind('\n')], format='sparse6'))
                else:
                    l.append(from_sparse6(d))
            else:
                l.append(from_sparse6(d))
        return l
    elif isinstance(data,file):
        strlist = data.readlines()
        l = []
        for s in strlist:
            l.append(Graph(s[:s.rfind('\n')], format='sparse6'))
        return l

def to_graph6(list, file = None):
    """
    Converts a list of SAGE graphs to a single string of graph6 graphs. If file is
    specified, then the string will be written quietly to the file.

    INPUT:
    list -- a Python list of SAGE Graphs
    file -- (optional) a file stream to write to (must be in 'w' mode)
    """
    l = ''
    for G in list:
        l.append(G.graph6_string() + '\n')
    if file is None:
        return l
    else:
        file.write(l)
        file.flush()

def to_sparse6(list, file = None):
    """
    Converts a list of SAGE graphs to a single string of sparse6 graphs. If file is
    specified, then the string will be written quietly to the file.

    INPUT:
    list -- a Python list of SAGE Graphs
    file -- (optional) a file stream to write to (must be in 'w' mode)
    """
    l = ''
    for G in list:
        l.append(G.sparse6_string() + '\n')
    if file is None:
        return l
    else:
        file.write(l)
        file.flush()

def to_graphics_arrays(list):
    """
    Returns a list of SAGE graphics arrays containing the graphs in list.
    The maximum number of graphs per array is 20 (5 rows of 4).  Use this
    function if there are too many graphs for the show_graphs function.
    The graphics arrays will contain 20 graphs each except potentially the
    last graphics array in the list.

    INPUT:
        list -- a list of SAGE graphs

    GRAPH PLOTTING:
    Uses the circular layout for graphs.  This allows for a nicer display
    in a small area and takes much less time to compute than the spring-
    layout algorithm for many graphs.

    EXAMPLES:
        sage: glist = []
        sage: for i in range(999):
        ... glist.append(graphs.RandomGNP(6,.45))
        ...
        sage: garray = graphs_list.to_graphics_arrays(glist)

        # Display the first graphics array in the list.
        sage.: garray[0].show()

        # Display the last graphics array in the list.
        sage.: garray[len(garray)-1].show()
    """
    from sage.plot.plot import graphics_array
    from sage.graphs import graph
    plist = []
    g_arrays = []
    for i in range (len(list)):
        if ( isinstance( list[i], graph.Graph ) ):
            pos = list[i].__get_pos__()
            if ( pos is None ):
                plist.append(list[i].plot(layout='circular', node_size=50, vertex_labels=False, graph_border=True))
            else: plist.append(list[i].plot(pos=pos, node_size=50, vertex_labels=False, graph_border=True))
        else:  raise TypeError, 'Param list must be a list of SAGE graphs.'

    num_arrays = len(plist)/20
    if ( len(plist)%20 > 0 ): num_arrays += 1
    rows = 5
    cols = 4

    for i in range (num_arrays-1):
        glist = []
        for j in range (rows*cols):
            glist.append(plist[ i*rows*cols + j ])
        ga = graphics_array(glist, rows, cols)
        ga.__set_figsize__([8,10])
        g_arrays.append(ga)

    last = len(plist)%20
    index = (num_arrays-1)*rows*cols
    last_rows = last/cols
    if ( last%cols > 0 ): last_rows += 1

    glist = []
    for i in range (last):
        glist.append(plist[ i + index])
    ga = graphics_array(glist, last_rows, cols)
    ga.__set_figsize__([8, 2*last_rows])
    g_arrays.append(ga)

    return g_arrays

def show_graphs(list):
    """
    Shows a maximum of 20 graphs from list in a sage graphics array.

    Note:  Raises ValueError if too many graphs are in the list.
    If the user wants to display more than 20 graphs, use the function
    to_graphics_arrays and show each graphics array individually.

    INPUT:
        list -- a list of SAGE graphs

    GRAPH PLOTTING:
    Uses the circular layout for graphs.  This allows for a nicer display
    in a small area and takes much less time to compute than the spring-
    layout algorithm for many graphs.

    EXAMPLES:
        # Create a list of graphs:
        sage: glist = []
        sage: glist.append(graphs.CompleteGraph(6))
        sage: glist.append(graphs.CompleteBipartiteGraph(4,5))
        sage: glist.append(graphs.BarbellGraph(7,4))
        sage: glist.append(graphs.CycleGraph(15))
        sage: glist.append(graphs.DiamondGraph())
        sage: glist.append(graphs.HouseGraph())
        sage: glist.append(graphs.HouseXGraph())
        sage: glist.append(graphs.KrackhardtKiteGraph())
        sage: glist.append(graphs.LadderGraph(5))
        sage: glist.append(graphs.LollipopGraph(5,6))
        sage: glist.append(graphs.PathGraph(15))
        sage: glist.append(graphs.PetersenGraph())
        sage: glist.append(graphs.StarGraph(17))
        sage: glist.append(graphs.WheelGraph(9))

        # Check that length is <= 20:
        sage: len(glist)
        14

        # Show the graphs in a graphics array:
        sage.: graphs_list.show_graphs(glist)

        # But more than 20 graphs will raise an exception:
        sage: glist21 = []
        sage: for i in range(21):
        ... glist21.append(graphs.RandomGNP(6,.45))
        ...
        sage: graphs_list.show_graphs(glist21)
        Traceback (most recent call last):
        ...
        ValueError: List is too long to display in a graphics array.  Try using the to_graphics_arrays function.

        # In this case, use to_graphics_arrays instead:
        sage: garrays = graphs_list.to_graphics_arrays(glist21)
        sage.: garrays[0].show()
        sage.: garrays[1].show()
    """
    if ( len(list) > 20 ):
        raise ValueError, 'List is too long to display in a graphics array.  Try using the to_graphics_arrays function.'

    from sage.plot.plot import graphics_array
    from sage.graphs import graph

    plist = []
    for i in range (len(list)):
        if ( isinstance( list[i], graph.Graph ) ):
            pos = list[i].__get_pos__()
            if ( pos is None ):
                plist.append(list[i].plot(layout='circular', node_size=50, vertex_labels=False, graph_border=True))
            else: plist.append(list[i].plot(pos=pos, node_size=50, vertex_labels=False, graph_border=True))
        else:  raise TypeError, 'Param list must be a list of SAGE graphs.'

    rows = len(list)/4
    if ( len(list)%4 > 0 ): rows += 1

    ga = graphics_array(plist, rows, 4)
    ga.__set_figsize__([8, 2*rows])

    ga.show()
    return

