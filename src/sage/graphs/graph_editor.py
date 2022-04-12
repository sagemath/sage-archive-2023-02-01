r"""
Graph editor
"""
# ****************************************************************************
#      Copyright (C) 2021 Jean-Florent Raymond
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


def graph_editor(graph=None, **editor_options):
    """
    Open a graph editor in the Jupyter notebook.

    INPUT:

    - ``graph`` - a :class:`Graph` instance (default: `None'); the
      graph to edit.

    - ``layout_options`` - options for the editor. 

    EXAMPLES::

        sage: g = graphs.CompleteGraph(3)
        sage: graph_editor(g)                       # not tested
        sage: graph_editor(graphs.HouseGraph())     # not tested
        sage: h = graphs.StarGraph(6)
    """
    try:
        from phitigra import GraphEditor
    except ModuleNotFoundError:
        raise ModuleNotFoundError("The graph editor is provided by the optional package phitigra, which currently not installed.")

    from .graph import Graph

    if graph is None:
        graph = Graph(0)

    return GraphEditor(graph, **editor_options)
