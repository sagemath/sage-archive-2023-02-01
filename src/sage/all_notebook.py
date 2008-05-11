"""nodoctest"""

from __future__ import with_statement


from sage.all import *
preparser(on=True)

from sage.calculus.predefined import x


sage_mode = 'notebook'

from sage.misc.latex import Latex, pretty_print_default, typeset, JSMath
latex = Latex(density=130)
latex_debug = Latex(debug=True, density=130)
slide = Latex(slide=True, density=256)
slide_debug = Latex(slide=True, debug=True, density=256)
# we need a global instance of this in order to get %jsmath to work...
jsmath = JSMath()

from sage.misc.python import python

from sage.misc.html import html

from sage.server.support import help






