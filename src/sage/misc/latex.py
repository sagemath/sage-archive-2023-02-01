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

__doc_exclude = ['_latex_file_', 'list_function', 'tuple_function', \
                 'bool_function', 'str_function', 'tmp_dir']

import os

import os.path

from misc import tmp_dir

def list_function(x):
    K = [latex(v) for v in x]
    if len(K) > 0 and max([len(r) for r in K]) > 60:
        sep = ', \n\\\\'
    else:
        sep = ', \n '
    return "\\left[" + sep.join([latex(v) for v in x]) + "\\right]"

def tuple_function(x):
    return "\\left(" + ", \n ".join([latex(v) for v in x]) + "\\right)"

def bool_function(x):
    if x:
        return "\\text{True}"
    else:
        return "\\text{False}"

def str_function(x):
    return "\\text{%s}"%x

    # this messes up too many things.

    #if not '#' in x:
    #    delim = '#'
    #elif not '@' in x:
    #    delim = '@'
    #elif not '~' in x:
    #    delim = '~'
    #return "\\verb%s%s%s"%(delim, x, delim)
    #return "\\begin{verbatim}%s\\end{verbatim}"%x

# One can add to the latex_table in order to install latexing
# functionality for other types.  (Suggested by Robert Kerns of Enthought.)

latex_table = {list: list_function, tuple:tuple_function, bool:bool_function,
               str: str_function}


def latex(x):
    """
    Output x formated for inclusion in a LaTeX document.

    The output should compile correctly if inserted into any latex
    document in math mode, assuming the amsmath package is included.
    No special macros should be required.
    """
    try:

        return x._latex_()

    except (AttributeError, TypeError):

        for k, f in latex_table.iteritems():
            if isinstance(x, k):
                return f(x)

        if x is None:
            return "\\mbox{\\rm None}"

        return str_function(str(x))


def _latex_file_(objects, title='SAGE', expert=True, debug=False, \
                 sep='$$ $$', tiny=False, center=False, math_left='$$',
                 math_right='$$',
                 extra_preamble='', brk=0):
    """
    Compute a latex file that defines a representation of each object in
    objects.

    INPUT:
        objects -- list (or object)
        size -- latex size of document ('small', 'tiny')
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

    s = '\\documentclass{article}\n\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\n\\input{macros}\n\n\n\\pagestyle{empty}\n%s\n\\begin{document}\n\\begin{center}{\\Large\\bf %s}\\end{center}\n\\thispagestyle{empty}\n %s\\%s '%(extra_preamble, title, center0, size)

    #s += "(If something is missing it may be on the next page or there may be errors in the latex.  Use view with {\\tt debug=True}.)\\vfill"
    s += '\\vfill'
    if process:
        for i in range(len(objects)):
            x = objects[i]
            L = latex(x)
            if not '\\begin{verbatim}' in L:
                s += '\\thispagestyle{empty}\\pagestyle{empty}\n\n %s %s %s'%(math_left, latex(x), math_right)
            else:
                s += '\\thispagestyle{empty}\\pagestyle{empty}\n\n %s'%latex(x)
            if i < len(objects)-1:
                s += '\n\n%s\n\n'%sep
    else:
        s += "\n\n".join([str(x) for x in objects])

    s += '\n\n\\vfill %s\\vfill\\end{document}'%center1
    if debug:
        print s

    # Finally break input so there is whitespace every brk characters, assuming brk > 0
    if brk > 0:
        # add a space to any block of brk characters or more.
        i = 0
        j = 0
        while i < len(s):
            if s[i] in ['\n', '\t', ' ']:
                j = i
            else:
                if i - j > brk:
                    s = s[:i] + ' ' + s[i:]
                    j = i
            i += 1

    return s


def view(objects, title='SAGE', zoom=4, expert=True, debug=False, \
         sep='$$ $$', tiny=False,  center=False):
    """
    Compute a latex representation of each object in objects, compile, and display
    using xdvi.  (Requires latex and xdvi be installed.)

    INPUT:
        objects -- list (or object)
        size -- latex size of document ('small', 'tiny')

    OUTPUT:
        Pops up xdvi with the objects displayed.
    """
    s = _latex_file_(objects, title=title, expert=expert,
                     debug=debug, sep=sep, tiny=tiny, center=center)

    SAGE_ROOT = os.environ['SAGE_ROOT']
    tmp = tmp_dir('sage_viewer')
    open('%s/sage.tex'%tmp,'w').write(s)
    os.system('ln -sf %s/devel/doc/commontex/macros.tex %s'%(SAGE_ROOT, tmp))
    O = open('%s/go'%tmp,'w')
    #O.write('export TEXINPUTS=%s/doc/commontex:.\n'%SAGE_ROOT)
    O.write('latex \\\\nonstopmode \\\\input{sage.tex}; xdvi -noscan -offsets 0.3 -paper 100000x100000 -s %s sage.dvi ; rm sage.* macros.* go ; cd .. ; rmdir %s'%(zoom,tmp))
    O.close()
    if not debug:
        direct = '1>/dev/null 2>/dev/null'
    else:
        direct = ''
    os.system('cd %s; chmod +x go; ./go %s&'%(tmp,direct))
    #return os.popen('cd %s; chmod +x go; ./go %s & '%(tmp,direct), 'r').read()

def png(x, filename, density=150, debug=False, brk=0):
    """
    Create a png image representation of x and save to the given
    filename.
    """
    s = _latex_file_([x], math_left='$\\displaystyle', math_right='$', title='',
                     debug=debug, tiny=False, extra_preamble='\\textheight=2\\textheight',
                     brk=brk)
    abs_path_to_png = os.path.abspath(filename)

    SAGE_ROOT = os.environ['SAGE_ROOT']
    tmp = tmp_dir('sage_viewer')
    open('%s/sage.tex'%tmp,'w').write(s)
    os.system('ln -sf %s/devel/doc/commontex/macros.tex %s'%(SAGE_ROOT, tmp))
    O = open('%s/go'%tmp,'w')
    go = 'latex \\\\nonstopmode \\\\input{sage.tex}; dvips -l =1 -f < sage.dvi > sage.ps ; convert -density %sx%s -trim sage.ps "%s";'%(density, density, abs_path_to_png)
    go += ' rm sage.* macros.* go ; cd .. ; rmdir %s'%tmp
    if debug:
        print go
    O.write(go)
    O.close()
    if not debug:
        direct = '1>/dev/null 2>/dev/null'
    else:
        direct = ''
    os.system('cd %s; chmod +x go; ./go %s&'%(tmp,direct))
    return s

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
        '-t\\text{a} + (t - 2)\\text{s} + (t^{12} + 2)\\text{}'
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


