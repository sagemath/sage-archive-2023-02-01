r"""
A module for dealing with lists of graphs.

AUTHOR:
    -- Robert L. Miller (2007-02-10): initial version
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
    use this if too many graphs for show function
    max graphs per array = 20 (5 rows of 4)
    appends graphics arrays in list
    then enter list and index and .show()
    """
    from sage.plot.plot import graphics_array
    from sage.graphs import graph
    plist = []
    g_arrays = []
    for i in range (len(list)):
        if ( isinstance( list[i], graph.Graph ) ):
            plist.append(list[i].plot(node_size=50, vertex_labels=False, graph_border=True))
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
        g_arrays.append(ga)

    last = len(plist)%20
    index = (num_arrays-1)*rows*cols
    last_rows = last/cols
    if ( last%cols > 0 ): last_rows += 1

    glist = []
    for i in range (last):
        glist.append(plist[ i + index])
    ga = graphics_array(glist, last_rows, cols)
    g_arrays.append(ga)

    return g_arrays

def show_graphs(list):
    if ( len(list) > 20 ):
        raise ValueError, 'List is too long to display in a graphics array.  Try using the to_graphics_array function.'
    return

