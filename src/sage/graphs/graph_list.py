r"""
A module for dealing with lists of graphs.

AUTHOR:
    -- Robert L. Miller (2007-02-10): initial version
"""

#*****************************************************************************
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
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

