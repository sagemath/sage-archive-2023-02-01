# -*- coding: utf-8 -*-
r"""
Computing Lovasz theta-function of graphs

AUTHORS:

- Dima Pasechnik (2015-06-30): Initial version

REFERENCE:

.. [Lovasz1979] László Lovász,
  "On the Shannon capacity of a graph",
  IEEE Trans. Inf. Th. 25(1979), 1-7.

Methods
-------
"""

def lovasz_theta(graph):
    """
    Returns value of Lovasz theta-function of graph,
    see :wikipedia:`Lovász_number`. For a graph :math:`G` it is denoted
    by :math:`\\theta(G)`, and the latter can be computed in polynomial time.
    Mathematically, the important property of it is 
    :math:`\\alpha(G)\leq\\theta(G)\leq\chi(\overline{G})`, for :math:`\\alpha(G)` 
    and :math:`\chi(\overline{G})` being, respectively, the maximum size of an independent set of :math:`G`
    and the chromatic number of the complement :math:`\overline{G}` of :math:`G`. 

    Implemented for undirected graphs only. Use
    to_undirected to convert a digraph to an undirected graph.
    
    Currently, the only solver to compute this value is the
    one from the optional/experimental package csdp,
    which you might need to install with 'sage -i csdp'

    EXAMPLES::

          sage: C=graphs.PetersenGraph()
          sage: lovasz_theta(C)                             # optional csdp
          4.0
          sage: lovasz_theta(graphs.CycleGraph(5))          # optional csdp
          2.236068

    TEST::

        sage: g = Graph()
        sage: lovasz_theta(g)
        0
    """
    n = graph.order()
    if n == 0:
        return 0

    from networkx import write_edgelist
    from sage.misc.temporary_file import tmp_filename
    import os, subprocess
    from sage.env import SAGE_LOCAL
    from sage.misc.package import is_package_installed

    if not is_package_installed('csdp'):
        raise NotImplementedError("Package csdp is required. Please install it with 'sage -i csdp'.")
 
    g = graph.relabel(inplace=False, perm=range(1,n+1)).networkx_graph()
    tf_name = tmp_filename()
    tf = open(tf_name, 'wb')
    tf.write(str(n)+'\n'+str(g.number_of_edges())+'\n')
    write_edgelist(g, tf, data=False)
    tf.close()
    lines = subprocess.check_output([os.path.join(SAGE_LOCAL, 'bin', 'theta'), tf_name])
    return float(lines.split()[-1])
