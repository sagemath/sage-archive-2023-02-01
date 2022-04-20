r"""
Graph editor widget

The ``phitigra`` package should be installed on your Sage installation.

AUTHORS:

- Radoslav Kirov (2009) -- initial version of the editor for use with sagenb
- Jean-Florent Raymond (2022) -- replacement with the ``phitigra`` package
"""
# ****************************************************************************
#      Copyright (C) 2022 Jean-Florent Raymond
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.features import PythonModule
lazy_import('phitigra', 'GraphEditor',
            feature=PythonModule('phitigra', spkg='phitigra'))


def graph_editor(graph=None, **display_options):
    """
    Return a graph editor widget.

    INPUT:

    - ``graph`` - a :class:`Graph` instance (default: ``None``); the
      graph to edit.

    - ``display_options`` - options for the widget.

    The graph editor widget can be displayed with Jupyter or JupyterLab.
    It is provided by the ``phitigra`` optional package, see
    https://github.com/jfraymond/phitigra for details about the
    possible options (changing the width/height of the canvas, the
    default size and color of vertices, etc.).

    EXAMPLES::

        sage: e = graph_editor()
        sage: e.show()                      # not tested

    Opening an existing graph::

        sage: G = graphs.PetersenGraph()
        sage: e = graph_editor(G)
        sage: e.show()                      # not tested

    Retrieving a copy of the drawn graph::

        sage: G = graphs.RandomGNP(10, 0.5)
        sage: e = graph_editor(G)
        sage: H = e.get_graph()
        sage: H == G and not H is G
        True

    Using different display options::

        sage: e = graph_editor(graphs.PetersenGraph(), width=300, height=300, default_radius=12, default_vertex_color='orange', default_edge_color='#666', show_vertex_labels=False)
        sage: e.show()                      # not tested
    """

    from .graph import Graph

    if graph is None:
        graph = Graph(0)

    return GraphEditor(graph, **display_options)
