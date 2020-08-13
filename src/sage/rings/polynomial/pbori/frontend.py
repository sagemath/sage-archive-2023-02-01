# Import basic functionality
r"""
This module defines an initial ring, and patches the declare_ring to use
a given context.


>>> x(0)
x(0)
>>> x(0)*x(0)
x(0)
>>> x(0) + x(0)
0
>>> x(9999)
x(9999)
>>> x(9999)*x(9999)
x(9999)
>>> x(9999) + x(9999)
0

>>> from sage.rings.polynomial.pbori.brial.frontend import *
>>> context = dict(globals())
>>> polybori_start(context) # doctest: +ELLIPSIS
ipbori...
>>> r = context['declare_ring']('abc')
>>> context['a']
a
>>> r.variable(0)
a
"""

from __future__ import print_function

from . import *
from .blocks import declare_ring as orig_declare_ring
from os import environ as env, path as os_path


def block_scheme_names(blocks):
    """Helper for Singular interface."""

    context = dict()
    from .blocks import declare_block_scheme
    declare_block_scheme(blocks, context)

    return list(context.keys())

ipbname = 'ipbori'


def polybori_copyright():
    print("""Copyright (c) 2007-2011 by The PolyBoRi Team.
  Michael Brickenstein (MFO)  brickenstein@mfo.de
  Alexander Dreyer (ITWM) alexander.dreyer@itwm.fraunhofer.de

The PolyBoRi Team is a joint project of
  Mathematisches Forschungsinstitut Oberwolfach (MFO), Germany
  Department of Mathematics, University of Kaiserslautern, Germany, and
  Fraunhofer Institute for Industrial Mathematics (ITWM), Kaiserslautern, Germany.

PolyBoRi incorporates the following works:
  The CU Decision Diagrams Package  Release 2.4.1 (CUDD) by Fabio Somenzi,
    Copyright (c) 1995-2004, Regents of the University of Colorado. All Rights Reserved.
  The M4RI Library - http://m4ri.sagemath.org
    Copyright (C) 2007-2010, Martin Albrecht, Gregory Bard, and The M4RI Team""")


def polybori_license():
    print("""ipbori and the PolyBoRi framework are licensed under the terms of
the GNU General Public License (GPL) version 2 or later.
See http://www.gnu.org/licenses/ for details.""")


def polybori_start(global_context):
    def declare_ring(blocks, context=None):
        if context is None:
            context = global_context

        return orig_declare_ring(blocks, context)
    declare_ring.__doc__ = orig_declare_ring.__doc__
    global_context["declare_ring"] = declare_ring

    print(ipbname + """ -- The interactive command line tool of PolyBoRi %s
Type "polybori_copyright()" or "polybori_license()" for more information.
""" % global_context.get("polybori_version", ''))

# Here come the defaults
r = Ring(10000)
x = VariableFactory(r)
