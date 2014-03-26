r"""
This file contains some utility functions for the interface with dot2tex
"""
#*****************************************************************************
#      Copyright (C) 2010   Nicolas M. Thiery <nicolas.thiery at u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import re
from sage.misc.latex import latex

def have_dot2tex():
    """
    Returns whether ``dot2tex`` >= 2.8.7 and graphviz are installed
    and functional

    EXAMPLES::

        sage: sage.graphs.dot2tex_utils.have_dot2tex() # optional - dot2tex graphviz
        True
        sage: sage.graphs.dot2tex_utils.have_dot2tex() in [True, False]
        True
    """
    try:
        import dot2tex
        # Test for this required feature from dot2tex 2.8.7
        return dot2tex.dot2tex("graph {}", format = "positions") == {}
    except Exception:
        return False
    return True

def assert_have_dot2tex():
    """
    Tests whether ``dot2tex`` >= 2.8.7 and graphviz are installed and
    functional, and raises an error otherwise

    EXAMPLES::

        sage: sage.graphs.dot2tex_utils.assert_have_dot2tex() # optional - dot2tex graphviz
    """
    check_error_string = """
An error occurs while testing the dot2tex installation.

Please see :meth:`sage.graphs.generic_graph.GenericGraph.layout_graphviz`
and check the installation of graphviz and the dot2tex spkg.

For support, please contact <sage-combinat-devel at googlegroups.com>.
"""
    missing_error_string = """
dot2tex not available.

Please see :meth:`sage.graphs.generic_graph.GenericGraph.layout_graphviz`
for installation instructions.
"""
    try:
        import dot2tex
        if dot2tex.dot2tex("graph {}", format = "positions") != {}:
            raise RuntimeError(check_error_string)
    except ImportError:
        raise RuntimeError(missing_error_string)

def quoted_latex(x):
    """
    Strips the latex representation of ``x`` to make it suitable for a
    ``dot2tex`` string.

    EXAMPLES::

        sage: sage.graphs.dot2tex_utils.quoted_latex(matrix([[1,1],[0,1],[0,0]]))
        '\\left(\\begin{array}{rr}1 & 1 \\\\0 & 1 \\\\0 & 0\\end{array}\\right)'
    """
    return re.sub("\"|\r|(%[^\n]*)?\n","", latex(x))

def quoted_str(x):
    """
    Strips the string representation of ``x`` to make it suitable for
    a ``dot2tex`` string, and especially a node label (``dot2tex``
    gets confused by newlines, and braces)

    EXAMPLES::

        sage: sage.graphs.dot2tex_utils.quoted_str(matrix([[1,1],[0,1],[0,0]]))
        '[1 1]\\n\\\n[0 1]\\n\\\n[0 0]'
        sage: print sage.graphs.dot2tex_utils.quoted_str(matrix([[1,1],[0,1],[0,0]]))
        [1 1]\n\
        [0 1]\n\
        [0 0]
    """
    return re.sub("\n",r"\\n\\"+"\n", re.sub("\"|\r|}|{","", str(x)))

def key(x):
    r"""
    Strips the string representation of ``x`` of quotes and newlines
    to get a key suitable for naming vertices in a dot2tex string.

    EXAMPLES::

        sage: sage.graphs.dot2tex_utils.key(matrix([[1,1],[0,1],[0,0]]))
        '110100'
        sage: sage.graphs.dot2tex_utils.key("blah{bleh}\nblih{")
        'blahblehblih'
    """
    return re.sub("[\\\'\"\[\]() \t\r\n{}]","", str(x))

def key_with_hash(x):
    """
    Same as :function:`key`, except that the hash of the object is
    prepended to better ensure uniqueness. This requires the object to
    be hashable, but this is anyway necessary to use it as a vertex of
    a graph.

    EXAMPLES::

        sage: sage.graphs.dot2tex_utils.key_with_hash(3)
        '3_3'
        sage: sage.graphs.dot2tex_utils.key_with_hash((1,2,3))
        '1,2,3_...'

    """
    return key(x)+"_"+str(hash(x))
