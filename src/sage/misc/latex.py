"""
Latex printing support

In order to support latex formating, an object should define a special
method _latex_(self) that returns a string.
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

from misc import tmp_dir

def list_function(x):
    return "[" + ", \n ".join([latex(v) for v in x]) + "]"

def tuple_function(x):
    return "(" + ", \n ".join([latex(v) for v in x]) + ")"

# One can add to the latex_table in order to install latexing
# functionality for other types.  (Suggested by Robert Kerns of UCSD.)

latex_table = {list: list_function, tuple:tuple_function}


def latex(x):
    """
    Output x formated for inclusion in a LaTeX document.

    The output should compile correctly if inserted into any latex
    document in math mode, assuming the amsmath package is included.
    No special macros should be required.
    """
    try:

        return x._latex_()

    except AttributeError:
        for k, f in latex_table.iteritems():
            if isinstance(x, k):
                return f(x)

        if x is None:
            return "\\mbox{\\rm None}"

        return str(x)


def view(objects, title='SAGE', zoom=4, expert=True, debug=False, \
         sep='$$ $$', tiny=False, center=False):
    """
    Compute a latex representation of each object in objects, compile, and display
    using xdvi.  (Requires latex and xdvi be installed.)

    INPUT:
        objects -- list (or object)
        size -- latex size of document ('small', 'tiny')

    OUTPUT:
        Pops up xdvi with the objects displayed.
    """
    process = True
    if hasattr(objects, '_latex_'):
        objects = [objects]

    if hasattr(objects, '__doc__') and hasattr(objects, 'func_name'):
        process = False
        title = "\\begin{verbatim}%s\\end{verbatim}"%objects.func_name
        objects = [objects.__doc__]

    if not isinstance(objects, list):
        objects = [objects]

    if expert:
        expert='-expert'
    else:
        expert=''

    if tiny:
        size='tiny'
    else:
        size='small'

    if center:
        center0 = '\\begin{center}'
        center1 = '\\end{center}'
    else:
        center0 =''
        center1 = ''

    if not debug:
        direct = '1>/dev/null 2>/dev/null'
    else:
        direct = ''
    s = '\\documentclass{article}\n\\usepackage{amsmath}\\usepackage{amssymb}\\usepackage{amsfonts}\n\\input{macros}\n\n\n\\pagestyle{empty}\n\\usepackage{fullpage}\n\\begin{document}\n\\begin{center}{\\Large\\bf %s}\\end{center}\n\\thispagestyle{empty}\n %s\\%s '%(title, center0, size)

    #s += "(If something is missing it may be on the next page or there may be errors in the latex.  Use view with {\\tt debug=True}.)\\vfill"
    s += '\\vfill'
    if process:
        for i in range(len(objects)):
            x = objects[i]
            s += '\\thispagestyle{empty}\\pagestyle{empty} $$%s$$'%latex(x)
            if i < len(objects)-1:
                s += '\n\n%s\n\n'%sep
    else:
        s += "\n\n".join([str(x) for x in objects])

    s += '\\vfill %s\\vfill\\end{document}'%center1
    if debug:
        print s
    SAGE_ROOT = os.environ['SAGE_ROOT']
    tmp = tmp_dir('sage_viewer')

    open('%s/sage.tex'%tmp,'w').write(s)
    os.system('cp %s/devel/doc/commontex/macros.tex %s'%(SAGE_ROOT, tmp))
    O = open('%s/go'%tmp,'w')
    #O.write('export TEXINPUTS=%s/doc/commontex:.\n'%SAGE_ROOT)
    O.write('latex \\\\nonstopmode \\\\input{sage.tex}; xdvi -noscan -offsets 0.3 -paper 100000x100000 -s %s sage.dvi ; rm sage.* macros.* go ; cd .. ; rmdir %s'%(zoom,tmp))
    O.close()
    os.system('cd %s; chmod +x go; ./go %s&'%(tmp,direct))
    #return os.popen('cd %s; chmod +x go; ./go %s & '%(tmp,direct), 'r').read()


def coeff_repr(c):
    try:
        return c._latex_coeff_repr()
    except AttributeError:
        pass
    if isinstance(c, (int, long, float)):
        return str(c)
    s = latex(c)
    if s.find("+") != -1 or s.find("-") != -1:
        return "(%s)"%s
    return s

def repr_lincomb(symbols, coeffs):
    """
    Compute a latex representation of a linear combination of some
    formal symbols.

    INPUT:
        symbols -- list of symbols
        coeffs -- list of coefficients of the symbols

    OUTPUT:
        str -- a string

    EXAMPLES:
        sage: t = PolynomialRing(Q, 't').0
        sage: from sage.misc.latex import repr_lincomb
        sage: repr_lincomb(['a', 's', ''], [-t, t - 2, t^12 + 2])
        '-ta + (t - 2)s + (t^{12} + 2)'
    """
    s = ""
    first = True
    i = 0

    all_atomic = True
    for c in coeffs:
        b = latex(symbols[i])
        if c != 0:
            coeff = coeff_repr(c)
            if not first:
                coeff = " + %s"%coeff
            else:
                coeff = "%s"%coeff
            s += "%s%s"%(coeff, b)
            first = False
        i += 1
    if first:
        s = "0"
    s = s.replace("+ -","- ")
    return s


