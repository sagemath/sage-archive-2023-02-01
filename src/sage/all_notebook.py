from sage.all import *
from sagenb.notebook.all import *

preparser(on=True)

from sage.calculus.predefined import x


sage_mode = 'notebook'

from sage.misc.latex import Latex, pretty_print_default, MathJax
latex = Latex(density=130)
latex_debug = Latex(debug=True, density=130)
slide = Latex(slide=True, density=256)
slide_debug = Latex(slide=True, debug=True, density=256)
pdflatex = Latex(density=130, pdflatex=True)
pdflatex_debug = Latex(density=130, pdflatex=True, debug=True)

from sage.misc.python import python

from sage.misc.html import html

from sage.server.support import help

from sagenb.misc.support import automatic_names

sage.misc.session.init()




