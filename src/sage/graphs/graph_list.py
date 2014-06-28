r"""
Lists of graphs

AUTHORS:

- Robert L. Miller (2007-02-10): initial version

- Emily A. Kirkman (2007-02-13): added show functions
  (to_graphics_array and show_graphs)
"""

#*****************************************************************************
#           Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

def from_whatever(data):
    """
    Returns a list of Sage Graphs, given a list of whatever kind of
    data.

    INPUT:


    -  ``data`` - can be a string, a list of strings, or a
       file stream, or whatever.


    EXAMPLE::

        sage: l = ['N@@?N@UGAGG?gGlKCMO',':P_`cBaC_ACd`C_@BC`ABDHaEH_@BF_@CHIK_@BCEHKL_BIKM_BFGHI']
        sage: graphs_list.from_whatever(l)
        [Graph on 15 vertices, Looped multi-graph on 17 vertices]
    """
    from sage.graphs.graph import Graph
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
                nn = d.rfind('\n')
                if nn == -1:
                    sparse = bool(d[0] == ':')
                    l.append(Graph(d, sparse=sparse))
                elif len(d) == nn + 1:
                    sparse = bool(d[0] == ':')
                    l.append(Graph(d[:nn], sparse=sparse))
                else:
                    l.append(from_whatever(d))
            else:
                l.append(from_whatever(d))
        return l
    if isinstance(data, str):
        data = data.split('\n')
        l = []
        for d in data:
            if not d == '':
                sparse = bool(d[0] == ':')
                l.append(Graph(d, sparse=sparse))
        return l

def from_graph6(data):
    """
    Returns a list of Sage Graphs, given a list of graph6 data.

    INPUT:


    -  ``data`` - can be a string, a list of strings, or a
       file stream.


    EXAMPLE::

        sage: l = ['N@@?N@UGAGG?gGlKCMO','XsGGWOW?CC?C@HQKHqOjYKC_uHWGX?P?~TqIKA`OA@SAOEcEA??']
        sage: graphs_list.from_graph6(l)
        [Graph on 15 vertices, Graph on 25 vertices]
    """
    from sage.graphs.graph import Graph
    if isinstance(data,str):
        data = data.split('\n')
        l = []
        for d in data:
            if not d == '':
                l.append(Graph(d, format = 'graph6'))
        return l
    elif isinstance(data,list):
        l = []
        for d in data:
            if isinstance(d, str):
                nn = d.rfind('\n')
                if nn == -1:
                    l.append(Graph(d,format='graph6'))
                elif len(d) == nn + 1:
                    l.append(Graph(d[:nn], format='graph6'))
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
    Returns a list of Sage Graphs, given a list of sparse6 data.

    INPUT:


    -  ``data`` - can be a string, a list of strings, or a
       file stream.


    EXAMPLE::

        sage: l = [':P_`cBaC_ACd`C_@BC`ABDHaEH_@BF_@CHIK_@BCEHKL_BIKM_BFGHI',':f`??KO?B_OOSCGE_?OWONDBO?GOJBDO?_SSJdApcOIG`?og_UKEbg?_SKFq@[CCBA`p?oYMFp@gw]Qaa@xEMHDb@hMCBCbQ@ECHEcAKKQKFPOwo[PIDQ{KIHEcQPOkVKEW_WMNKqPWwcRKOOWSKIGCqhWt??___WMJFCahWzEBa`xOu[MpPPKqYNoOOOKHHDBPs|??__gWMKEcAHKgTLErqA?A@a@G{kVLErs?GDBA@XCs\\NggWSOJIDbHh@?A@aF']
        sage: graphs_list.from_sparse6(l)
        [Looped multi-graph on 17 vertices, Looped multi-graph on 39 vertices]
    """
    from sage.graphs.graph import Graph
    if isinstance(data,str):
        data = data.split('\n')
        l = []
        for d in data:
            if not d == '':
                l.append(Graph(d, format = 'sparse6', sparse=True))
        return l
    elif isinstance(data,list):
        l = []
        for d in data:
            if isinstance(d, str):
                nn = d.rfind('\n')
                if nn == -1:
                    l.append(Graph(d, format='sparse6', sparse=True))
                elif len(d) == nn + 1:
                    l.append(Graph(d[:nn], format='sparse6', sparse=True))
                else:
                    l.append(from_sparse6(d))
            else:
                l.append(from_sparse6(d))
        return l
    elif isinstance(data,file):
        strlist = data.readlines()
        l = []
        for s in strlist:
            l.append(Graph(s[:s.rfind('\n')], format='sparse6', sparse=True))
        return l

def to_graph6(list, file = None, output_list=False):
    r"""
    Converts a list of Sage graphs to a single string of graph6 graphs.
    If file is specified, then the string will be written quietly to
    the file. If output_list is True, then a list of strings will be
    returned, one string per graph.

    INPUT:


    -  ``list`` - a Python list of Sage Graphs

    -  ``file`` - (optional) a file stream to write to
       (must be in 'w' mode)

    -  ``output_list`` - False - output is a string True -
       output is a list of strings (ignored if file gets specified)


    EXAMPLE::

        sage: l = [graphs.DodecahedralGraph(), graphs.PetersenGraph()]
        sage: graphs_list.to_graph6(l)
        'ShCHGD@?K?_@?@?C_GGG@??cG?G?GK_?C\nIheA@GUAo\n'
    """
    l = ''
    for G in list:
        l += G.graph6_string() + '\n'
    if file is None:
        if output_list:
            a = l.split('\n')
            a = a[:len(a)-1]
            return a
        else:
            return l
    else:
        file.write(l)
        file.flush()

def to_sparse6(list, file = None, output_list=False):
    r"""
    Converts a list of Sage graphs to a single string of sparse6
    graphs. If file is specified, then the string will be written
    quietly to the file. If output_list is True, then a list of
    strings will be returned, one string per graph.

    INPUT:


    -  ``list`` - a Python list of Sage Graphs

    -  ``file`` - (optional) a file stream to write to
       (must be in 'w' mode)

    -  ``output_list`` - False - output is a string True -
       output is a list of strings (ignored if file gets specified)


    EXAMPLE::

        sage: l = [graphs.DodecahedralGraph(), graphs.PetersenGraph()]
        sage: graphs_list.to_sparse6(l)
        ':S_`abcaDe`Fg_HijhKfLdMkNcOjP_BQ\n:I`ES@obGkqegW~\n'
    """
    l = ''
    for G in list:
        l += G.sparse6_string() + '\n'
    if file is None:
        if output_list:
            a = l.split('\n')
            a = a[:len(a)-1]
            return a
        else:
            return l
    else:
        file.write(l)
        file.flush()

def to_graphics_arrays(list, **kwds):
    """
    Returns a list of Sage graphics arrays containing the graphs in
    list. The maximum number of graphs per array is 20 (5 rows of 4).
    Use this function if there are too many graphs for the show_graphs
    function. The graphics arrays will contain 20 graphs each except
    potentially the last graphics array in the list.

    INPUT:


    -  ``list`` - a list of Sage graphs


    GRAPH PLOTTING: Defaults to circular layout for graphs. This allows
    for a nicer display in a small area and takes much less time to
    compute than the spring- layout algorithm for many graphs.

    EXAMPLES::

        sage: glist = []
        sage: for i in range(999):
        ...    glist.append(graphs.RandomGNP(6,.45))
        ...
        sage: garray = graphs_list.to_graphics_arrays(glist)

    Display the first graphics array in the list.

    ::

        sage: garray[0].show()

    Display the last graphics array in the list.

    ::

        sage: garray[len(garray)-1].show()

    See the .plot() or .show() documentation for an individual graph
    for options, all of which are available from to_graphics_arrays

    ::

        sage: glist = []
        sage: for _ in range(10):
        ...       glist.append(graphs.RandomLobster(41, .3, .4))
        sage: w = graphs_list.to_graphics_arrays(glist, layout='spring', vertex_size=20)
        sage: len(w)
        1
        sage: w[0]
    """
    from sage.plot.plot import graphics_array
    from sage.graphs import graph
    plist = []
    g_arrays = []
    for i in range (len(list)):
        if ( isinstance( list[i], graph.GenericGraph ) ):
            pos = list[i].get_pos()
            if ( pos is None ):
                if 'layout' not in kwds:
                    kwds['layout'] = 'circular'
                if 'vertex_size' not in kwds:
                    kwds['vertex_size'] = 50
                if 'vertex_labels' not in kwds:
                    kwds['vertex_labels'] = False
                kwds['graph_border'] = True
                plist.append(list[i].plot(**kwds))
            else: plist.append(list[i].plot(pos=pos, vertex_size=50, vertex_labels=False, graph_border=True))
        else:  raise TypeError('Param list must be a list of Sage (di)graphs.')

    num_arrays = len(plist) // 20
    if ( len(plist)%20 > 0 ): num_arrays += 1
    rows = 5
    cols = 4

    for i in range (num_arrays-1):
        glist = []
        for j in range (rows*cols):
            glist.append(plist[ i*rows*cols + j ])
        ga = graphics_array(glist, rows, cols)
        ga._set_figsize_([8,10])
        g_arrays.append(ga)

    last = len(plist)%20
    if ( last == 0 and len(plist) != 0 ): last = 20
    index = (num_arrays-1)*rows*cols
    last_rows = last/cols
    if ( last%cols > 0 ): last_rows += 1

    glist = []
    for i in range (last):
        glist.append(plist[ i + index])
    ga = graphics_array(glist, last_rows, cols)
    ga._set_figsize_([8, 2*last_rows])
    g_arrays.append(ga)

    return g_arrays

def show_graphs(list, **kwds):
    """
    Shows a maximum of 20 graphs from list in a sage graphics array. If
    more than 20 graphs are given in the list argument, then it will
    display one graphics array after another with each containing at
    most 20 graphs.

    Note that if to save the image output from the notebook, you must
    save each graphics array individually. (There will be a small space
    between graphics arrays).

    INPUT:


    -  ``list`` - a list of Sage graphs


    GRAPH PLOTTING: Defaults to circular layout for graphs. This allows
    for a nicer display in a small area and takes much less time to
    compute than the spring- layout algorithm for many graphs.

    EXAMPLES: Create a list of graphs::

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

    Check that length is = 20::

        sage: len(glist)
        14

    Show the graphs in a graphics array::

        sage: graphs_list.show_graphs(glist)

    Here's an example where more than one graphics array is used::

        sage: gq = GraphQuery(display_cols=['graph6'],num_vertices=5)
        sage: g = gq.get_graphs_list()
        sage: len(g)
        34
        sage: graphs_list.show_graphs(g)

    See the .plot() or .show() documentation for an individual graph
    for options, all of which are available from to_graphics_arrays

    ::

        sage: glist = []
        sage: for _ in range(10):
        ...       glist.append(graphs.RandomLobster(41, .3, .4))
        sage: graphs_list.show_graphs(glist, layout='spring', vertex_size=20)
    """
    ga_list = to_graphics_arrays(list, **kwds)

    for i in range (len(ga_list)):
        (ga_list[i]).show()

    return


