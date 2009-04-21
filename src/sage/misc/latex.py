"""
Latex printing support

In order to support latex formating, an object should define a
special method _latex_(self) that returns a string.
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

EMBEDDED_MODE = False

LATEX_HEADER='\\documentclass{article}\\usepackage{fullpage}\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\\usepackage{graphicx}\usepackage{pstricks}\pagestyle{empty}\n'

#SLIDE_HEADER='\\documentclass[landscape]{slides}\\usepackage{fullpage}\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\\usepackage{graphicx}\usepackage{pstricks}\pagestyle{empty}\n'
SLIDE_HEADER='\\documentclass[a0,8pt]{beamer}\\usepackage{fullpage}\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\\usepackage{graphicx}\usepackage{pstricks}\pagestyle{empty}\n\\textwidth=1.1\\textwidth\\textheight=2\\textheight'


import os, shutil

import os.path

import random

from misc import tmp_dir
import sage_eval
from sage.misc.misc import SAGE_DOC

_have_dvipng = None
def have_dvipng():
    """
    Return True if this computer has the program dvipng.

    The first time it is run, this function caches its result in the
    variable ``_have_dvipng``, and any subsequence time, it just
    checks the value of the variable.

    EXAMPLES::

        sage: from sage.misc.latex import have_dvipng
        sage: sage.misc.latex._have_dvipng is None
        True
        sage: have_dvipng() # random
        True
        sage: sage.misc.latex._have_dvipng is None
        False
        sage: sage.misc.latex._have_dvipng == have_dvipng()
        True
    """
    global _have_dvipng
    if _have_dvipng is None:
        _have_dvipng = not bool(os.system('which dvipng >/dev/null'))
    return _have_dvipng

def list_function(x):
    r"""
    Returns the LaTeX code for a list x.

    INPUT: x - a list

    EXAMPLES::

        sage: from sage.misc.latex import list_function
        sage: list_function([1,2,3])
        '\\left[1, \n 2, \n 3\\right]'
        sage: latex([1,2,3])  # indirect doctest
        \left[1,
        2,
        3\right]
        sage: latex([Matrix(ZZ,3,range(9)), Matrix(ZZ,3,range(9))]) # indirect doctest
        \left[\left(\begin{array}{rrr}
        0 & 1 & 2 \\
        3 & 4 & 5 \\
        6 & 7 & 8
        \end{array}\right),
        \\\left(\begin{array}{rrr}
        0 & 1 & 2 \\
        3 & 4 & 5 \\
        6 & 7 & 8
        \end{array}\right)\right]
    """
    K = [latex(v) for v in x]
    if len(K) > 0 and sum([len(r) for r in K]) > 80:
        if EMBEDDED_MODE:
            s = '\\begin{array}{l}'
            K[0] = '[' + K[0]
            K[-1] = K[-1] + ']'
            s += ',\\\\\n'.join(K)
            s += '\\end{array}'
            return s
        else:
            sep = ', \n\\\\'
    else:
        sep = ', \n '
    return "\\left[" + sep.join(K) + "\\right]"

def tuple_function(x):
    r"""
    Returns the LaTeX code for a tuple x.

    INPUT: x - a tuple

    EXAMPLES::

        sage: from sage.misc.latex import tuple_function
        sage: tuple_function((1,2,3))
        '\\left(1, \n 2, \n 3\\right)'
    """
    return "\\left(" + ", \n ".join([latex(v) for v in x]) + "\\right)"

def bool_function(x):
    r"""
    Returns the LaTeX code for a tuple x.

    INPUT: x - boolean

    EXAMPLES::

        sage: from sage.misc.latex import bool_function
        sage: bool_function(2==3)
        '\\mbox{\\rm False}'
        sage: bool_function(3==(2+1))
        '\\mbox{\\rm True}'
    """
    if x:
        s = "\\mbox{\\rm True}"
    else:
        s = "\\mbox{\\rm False}"
    if EMBEDDED_MODE:
        return s[5:]
    return s

def str_function(x):
    r"""
    Returns the LaTeX code for a string x.

    INPUT: x - a string

    EXAMPLES::

        sage: from sage.misc.latex import str_function
        sage: str_function('hello world')
        '\\text{hello world}'
    """
    #if EMBEDDED_MODE:
    return '\\text{%s}'%(x.replace('_','\\_'))
    #return "\\mbox{\\rm %s}"%x'

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
               str: str_function, int:str, long:str, float:str}

class LatexExpr(str):
    def __init__(self, x):
        str.__init__(self, x)

    def __repr__(self):
        return str(self)

    def _latex_(self):
        return str(self)

from sage.structure.sage_object import SageObject

class _Latex_prefs_object(SageObject):
    """
    An object that holds LaTeX global preferences.
    """
    def __init__(self, bb=False, delimiters=["(", ")"]):
        """
        Define an object that holds LaTeX global preferences.
        """
        self._option = {}
        self._option["blackboard_bold"] = bb
        self._option["matrix_delimiters"] = list(delimiters)
        self._option["vector_delimiters"] = list(delimiters)

_Latex_prefs = _Latex_prefs_object()

##############################################################
# The Latex class is used to make slides and latex output in
# the SAGE Notebook
#########################################

class Latex:
    r"""nodetex
    Enter, e.g.,

    ::

                %latex
                The equation $y^2 = x^3 + x$ defines an elliptic curve.
                We have $2006 = Sage{factor(2006)}$.


    in an input cell to get a typeset version (care of slitex). Use
    ``%latex_debug`` to get debugging output.

    Use ``latex(...)`` to typeset a Sage object.

    Use ``%slide`` instead to typeset slides.

    .. warning::

       You must have dvipng (or dvips and convert) installed
       on your operating system, or this command won't work.

    """
    def __init__(self, debug=False, slide=False, density=150):
        self.__debug = debug
        self.__slide = slide
        self.__density = density

    def __call__(self, x):
        if hasattr(x, '_latex_'):
            return LatexExpr(x._latex_())

        for k, f in latex_table.iteritems():
            if isinstance(x, k):
                return LatexExpr(f(x))

        if x is None:
            return LatexExpr("\\mbox{\\rm None}")

        return LatexExpr(str_function(str(x)))

    def _latex_preparse(self, s, locals):
        """
        Replace instances of '\sage{x}' in s with the LaTeX version of
        x in the running session.

        EXAMPLES:
            sage: s = 2
            sage: sage.misc.latex.Latex()._latex_preparse('\sage{s}', locals())
            '2'
        """
        i0 = -1
        while True:
            i = s.find('\\sage{')
            if i == -1 or i == i0:
                return s
            i0 = i
            t = s[i+6:]
            j = t.find('}')
            if j == -1:
                return s

            var = t[:j]
            try:
                k = latex(sage_eval.sage_eval(var, locals))
            except Exception, msg:
                print msg
                k = '\\mbox{\\rm [%s undefined]}'%var
            s = s[:i] + k + t[j+1:]

    def eval(self, x, globals, strip=False, filename=None, debug=None,
             density=None, locals={}):
        """
        INPUT:
            globals -- a globals dictionary


        -  ``x`` - string to evaluate.

        -  ``strip`` - ignored

        -  ``filename`` - output filename

        -  ``debug`` - whether to print verbose debugging
           output

        -  ``density`` - how big output image is.

        -  ``locals`` - extra local variables used when
           evaluating Sage.. code in x.

        .. warning::

           You must have dvipng (or dvips and convert) installed on
           your operating system, or this command won't work.

        """
        from sage.misc.latex_macros import sage_latex_macros
        MACROS="\n" + "\n".join(sage_latex_macros) + "\n"

        if density is None:
            density = self.__density
        if filename is None:
            filename = 'sage%s'%random.randint(1,100) # to defeat browser caches
        else:
            filename = os.path.splitext(filename)[0]  # get rid of extension
        base = tmp_dir()
        orig_base, filename = os.path.split(os.path.abspath(filename))
        if len(filename.split()) > 1:
            raise ValueError, "filename must contain no spaces"
        if debug is None:
            debug = self.__debug
        x = self._latex_preparse(x, locals)
        O = open('%s/%s.tex'%(base,filename),'w')
        if self.__slide:
            O.write(SLIDE_HEADER)
            O.write(MACROS)
            O.write('\\begin{document}\n\n')
        else:
            O.write(LATEX_HEADER)
            O.write(MACROS)
            O.write('\\begin{document}\n')

        O.write(x)
        if self.__slide:
            O.write('\n\n\\end{document}')
        else:
            O.write('\n\n\\end{document}\n')

        O.close()
        if not debug:
            redirect=' 2>/dev/null 1>/dev/null '
        else:
            redirect=''
        lt = 'cd "%s"&& sage-native-execute latex \\\\nonstopmode \\\\input{%s.tex} %s'%(base, filename, redirect)
        if have_dvipng():
            dvipng = 'sage-native-execute dvipng -q -T bbox -D %s %s.dvi -o %s.png'%(density, filename, filename)
            cmd = ' && '.join([lt, dvipng])

        else:
            dvips = 'sage-native-execute dvips %s.dvi %s'%(filename, redirect)
            convert = 'sage-native-execute convert -density %sx%s -trim %s.ps %s.png %s '%\
                      (density,density, filename, filename, redirect)
            cmd = ' && '.join([lt, dvips, convert])
        if debug:
            print cmd
        e = os.system(cmd + ' ' + redirect)
        if e:
            print "An error occured."
            try:
                print open(filename + '.log').read()
            except IOError:
                pass
            return 'Error latexing slide.'
        shutil.copy(base + '/' + filename + '.png', orig_base + '/'+filename + '.png')
        shutil.rmtree(base)
        return ''

    def blackboard_bold(self, t = None):
        """
        Controls whether Sage uses blackboard bold or ordinary bold
        face for typesetting ZZ, RR, etc.

        INPUT:

        - ``t`` -- boolean or None

        OUTPUT: if t is None, return the current setting (True or False).

        If t == True, use blackboard bold (\\mathbb); otherwise use
        boldface (\\mathbf).

        EXAMPLES::

            sage: latex.blackboard_bold()
            False
            sage: latex.blackboard_bold(True)
            sage: latex.blackboard_bold()
            True
            sage: latex.blackboard_bold(False)
        """
        if t is None:
            return _Latex_prefs._option["blackboard_bold"]
        from latex_macros import sage_latex_macros, sage_jsmath_macros, sage_configurable_latex_macros, convert_latex_macro_to_jsmath
        global sage_latex_macros
        global sage_jsmath_macros
        old = _Latex_prefs._option["blackboard_bold"]
        _Latex_prefs._option["blackboard_bold"] = bool(t)
        if bool(old) != bool(t):
            if old:
                old_macro = "\\newcommand{\\Bold}[1]{\\mathbb{#1}}"
            else:
                old_macro = "\\newcommand{\\Bold}[1]{\\mathbf{#1}}"
            if bool(t):
                macro = "\\newcommand{\\Bold}[1]{\\mathbb{#1}}"
            else:
                macro = "\\newcommand{\\Bold}[1]{\\mathbf{#1}}"
            sage_latex_macros.remove(old_macro)
            sage_configurable_latex_macros.remove(old_macro)
            sage_latex_macros.append(macro)
            sage_configurable_latex_macros.append(macro)
            sage_jsmath_macros.remove(convert_latex_macro_to_jsmath(old_macro))
            sage_jsmath_macros.append(convert_latex_macro_to_jsmath(macro))

    def matrix_delimiters(self, left=None, right=None):
        r"""
        Change the left and right delimiters for the LaTeX representation
        of matrices

        INPUT:

        - ``left``, ``right`` - strings or None

        If both ``left`` and ``right`` are ``None``, then return the
        current delimiters.  Otherwise, set the left and/or right
        delimiters, whichever are specified.

        Good choices for ``left`` and ``right`` are any delimiters which
        LaTeX understands and knows how to resize; some examples are:

        - parentheses: '(', ')'
        - brackets: '[', ']'
        - braces: '\\{', '\\}'
        - vertical lines: '|'
        - angle brackets: '\\langle', '\\rangle'

        .. note::

           Putting aside aesthetics, you may combine these in any way
           imaginable; for example, you could set ``left`` to be a
           right-hand bracket ']' and ``right`` to be a right-hand
           brace '\\}', and it will be typeset correctly.

        EXAMPLES::

            sage: a = matrix(1, 1, [17])
            sage: latex(a)
            \left(\begin{array}{r}
            17
            \end{array}\right)
            sage: latex.matrix_delimiters("[", "]")
            sage: latex(a)
            \left[\begin{array}{r}
            17
            \end{array}\right]
            sage: latex.matrix_delimiters(left="\\{")
            sage: latex(a)
            \left\{\begin{array}{r}
            17
            \end{array}\right]
            sage: latex.matrix_delimiters()
            ['\\{', ']']

        Restore defaults::

            sage: latex.matrix_delimiters("(", ")")
        """
        if left is None and right is None:
            return _Latex_prefs._option['matrix_delimiters']
        else:
            if left is not None:
                _Latex_prefs._option['matrix_delimiters'][0] = left
            if right is not None:
                _Latex_prefs._option['matrix_delimiters'][1] = right

    def vector_delimiters(self, left=None, right=None):
        r"""
        Change the left and right delimiters for the LaTeX representation
        of vectors

        INPUT:

        - ``left``, ``right`` - strings or None

        If both ``left`` and ``right`` are ``None``, then return the
        current delimiters.  Otherwise, set the left and/or right
        delimiters, whichever are specified.

        Good choices for ``left`` and ``right`` are any delimiters which
        LaTeX understands and knows how to resize; some examples are:

        - parentheses: '(', ')'
        - brackets: '[', ']'
        - braces: '\\{', '\\}'
        - vertical lines: '|'
        - angle brackets: '\\langle', '\\rangle'

        .. note::

           Putting aside aesthetics, you may combine these in any way
           imaginable; for example, you could set ``left`` to be a
           right-hand bracket ']' and ``right`` to be a right-hand
           brace '\\}', and it will be typeset correctly.

        EXAMPLES::

            sage: a = vector(QQ, [1,2,3])
            sage: latex(a)
            \left(1,2,3\right)
            sage: latex.vector_delimiters("[", "]")
            sage: latex(a)
            \left[1,2,3\right]
            sage: latex.vector_delimiters(right="\\}")
            sage: latex(a)
            \left[1,2,3\right\}
            sage: latex.vector_delimiters()
            ['[', '\\}']

        Restore defaults::

            sage: latex.vector_delimiters("(", ")")
        """
        if left is None and right is None:
            return _Latex_prefs._option['vector_delimiters']
        else:
            if left is not None:
                _Latex_prefs._option['vector_delimiters'][0] = left
            if right is not None:
                _Latex_prefs._option['vector_delimiters'][1] = right

# Note: latex used to be a separate function, which by default was
# only loaded in command-line mode: in the notebook, all_notebook.py
# defined (and still defines) latex by 'latex = Latex(density=130)'.
# Meanwhile, the __call__ method for Latex used to call the latex
# function.  This has been changed around so that the contents of the
# old latex function are now in Latex.__call__; thus the following
# assignment.

latex = Latex()
#########################################

def _latex_file_(objects, title='SAGE', debug=False, \
                 sep='', tiny=False, math_left='\\[',
                 math_right='\\]',
                 extra_preamble=''):
    r"""nodetex
    Produce a string to be used as a LaTeX file, containing a
    representation of each object in objects.

    INPUT:


    -  ``objects`` - list (or object)

    -  ``title`` - string (default: 'Sage'): title for the document

    -  ``math_left`` - string (default: '\\['), left delimiter for math mode

    -  ``math_right`` - string (default: '\\]'), right delimiter for math mode

    -  ``debug`` - bool (default: False): print verbose output

    -  ``sep`` - string (default: ''): separator between math objects

    -  ``tiny`` - bool (default: False): use 'tiny' font.

    -  ``extra_preamble`` - string (default: ''): extra LaTeX commands,
       inserted before "\\begin{document}"

    This creates a string intended to be a LaTeX file containing the
    LaTeX representations of objects. It contains the following:

      - a header (with documentclass and usepackage commands)

      - ``extra_preamble``

      - the title (centered)

      - a size specification if ``tiny`` is True

      - LaTeX representation of the first element of ``objects``,
        surrounded by ``math_left`` and ``math_right``

    Then if ``objects`` contains more than one element, for each
    remaining element:

      - the string ``sep``: you can use this, for example, to add
        vertical space between objects with ``sep='\\vspace{15mm}'``,
        or to add a horizontal line between objects with
        ``sep='\\hrule'``, or to insert a page break between objects
        with ``sep='\\newpage'``.

      - the LaTeX representation of the element

    The string ends with '\\end{document}'.

    EXAMPLES::

        sage: from sage.misc.latex import _latex_file_
        sage: _latex_file_(3, title="The number three")
        '\\documentclass{article}\\usepackage{fullpage}\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\\usepackage{graphicx}\\usepackage{pstricks}\\pagestyle{empty}\n\n\n\\newcommand{\\ZZ}{\\Bold{Z}}\n\\newcommand{\\RR}{\\Bold{R}}\n\\newcommand{\\CC}{\\Bold{C}}\n\\newcommand{\\QQ}{\\Bold{Q}}\n\\newcommand{\\QQbar}{\\overline{\\QQ}}\n\\newcommand{\\GF}[1]{\\Bold{F}_{#1}}\n\\newcommand{\\Zp}[1]{\\ZZ_{#1}}\n\\newcommand{\\Qp}[1]{\\QQ_{#1}}\n\\newcommand{\\Zmod}[1]{\\ZZ/#1\\ZZ}\n\\newcommand{\\CDF}{\\text{Complex Double Field}}\n\\newcommand{\\CIF}{\\Bold{C}}\n\\newcommand{\\CLF}{\\Bold{C}}\n\\newcommand{\\RDF}{\\Bold{R}}\n\\newcommand{\\RIF}{\\I \\R}\n\\newcommand{\\RLF}{\\Bold{R}}\n\\newcommand{\\RQDF}{\\Bold{R}}\n\\newcommand{\\CFF}{\\Bold{CFF}}\n\\newcommand{\\Bold}[1]{\\mathbf{#1}}\n\n\\begin{document}\n\\begin{center}{\\Large\\bf The number three}\\end{center}\n\\vspace{40mm}\\[3\\]\n\\end{document}'
        sage: _latex_file_([7, 8, 9], title="Why was six afraid of seven?", sep='\\vfill\\hrule\\vfill')
        '\\documentclass{article}\\usepackage{fullpage}\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\\usepackage{graphicx}\\usepackage{pstricks}\\pagestyle{empty}\n\n\n\\newcommand{\\ZZ}{\\Bold{Z}}\n\\newcommand{\\RR}{\\Bold{R}}\n\\newcommand{\\CC}{\\Bold{C}}\n\\newcommand{\\QQ}{\\Bold{Q}}\n\\newcommand{\\QQbar}{\\overline{\\QQ}}\n\\newcommand{\\GF}[1]{\\Bold{F}_{#1}}\n\\newcommand{\\Zp}[1]{\\ZZ_{#1}}\n\\newcommand{\\Qp}[1]{\\QQ_{#1}}\n\\newcommand{\\Zmod}[1]{\\ZZ/#1\\ZZ}\n\\newcommand{\\CDF}{\\text{Complex Double Field}}\n\\newcommand{\\CIF}{\\Bold{C}}\n\\newcommand{\\CLF}{\\Bold{C}}\n\\newcommand{\\RDF}{\\Bold{R}}\n\\newcommand{\\RIF}{\\I \\R}\n\\newcommand{\\RLF}{\\Bold{R}}\n\\newcommand{\\RQDF}{\\Bold{R}}\n\\newcommand{\\CFF}{\\Bold{CFF}}\n\\newcommand{\\Bold}[1]{\\mathbf{#1}}\n\n\\begin{document}\n\\begin{center}{\\Large\\bf Why was six afraid of seven?}\\end{center}\n\\vspace{40mm}\\[7\\]\n\n\\vfill\\hrule\\vfill\n\n\\[8\\]\n\n\\vfill\\hrule\\vfill\n\n\\[9\\]\n\\end{document}'
    """
    from sage.misc.latex_macros import sage_latex_macros
    MACROS="\n" + "\n".join(sage_latex_macros) + "\n"

    process = True
    if hasattr(objects, '_latex_'):
        objects = [objects]

    if hasattr(objects, '__doc__') and hasattr(objects, 'func_name'):
        process = False
        title = "\\begin{verbatim}%s\\end{verbatim}"%objects.func_name
        objects = [objects.__doc__]

    if not isinstance(objects, list):
        objects = [objects]

    if tiny:
        size='\\tiny\n'
    else:
        size=''

    s = LATEX_HEADER + '\n' + MACROS
    s += '%s\n\\begin{document}\n\\begin{center}{\\Large\\bf %s}\\end{center}\n%s'%(
        extra_preamble, title, size)

    #s += "(If something is missing it may be on the next page or there may be errors in the latex.  Use view with {\\tt debug=True}.)\\vfill"
    s += '\\vspace{40mm}'
    if process:
        for i in range(len(objects)):
            x = objects[i]
            L = latex(x)
            if not '\\begin{verbatim}' in L:
                s += '%s%s%s'%(math_left, latex(x), math_right)
            else:
                s += '%s'%latex(x)
            if i < len(objects)-1:
                s += '\n\n%s\n\n'%sep
    else:
        s += "\n\n".join([str(x) for x in objects])

    s += '\n\\end{document}'
    if debug:
        print s

    return s

class JSMathExpr:
    '''
    An arbitrary JSMath expression that can be nicely concatenated.
    '''
    def __init__(self, y):
        self.__y = y

    def __repr__(self):
        return str(self.__y)

    def __add__(self, y):
        return JSMathExpr(self.__y + y)

    def __radd__(self, y):
        return JSMathExpr(y + self.__y)

class JSMath:
    '''
    A simple object for rendering LaTeX input using JSMath.

    '''

    def __call__(self, x):
        return self.eval(x)

    def eval(self, x, globals=None, locals=None, mode='display'):
        try:
            # try to get a latex representation of the object
            x = x._latex_()
        except AttributeError:
            # otherwise just get the string representation
            x = str(x)


        # in JSMath:
        # inline math: <span class="math">...</span>
        # displaymath: <div class="math">...</div>
        from sage.misc.latex_macros import sage_configurable_latex_macros
        if 'display' == mode:
            return JSMathExpr('<html><div class="math">'
                              + ''.join(sage_configurable_latex_macros)
                              + '%s</div></html>'%x)
        elif 'inline' == mode:
            return JSMathExpr('<html><span class="math">'
                              + ''.join(sage_configurable_latex_macros)
                              + '%s</span></html>'%x)
        else:
            # what happened here?
            raise ValueError, "mode must be either 'display' or 'inline'"

def jsmath(x, mode='display'):
    r'''
    Attempt to nicely render an arbitrary SAGE object wih jsmath typesetting.
    Tries to call ._latex_() on x. If that fails, it will render a string
    representation of x.

    .. warning::

        2009-04: This function is deprecated; use ``html`` instead:
        replace ``jsmath('MATH', mode='display')`` with ``html('$$MATH$$')``,
        and replace ``jsmath('MATH', mode='inline')`` with ``html('$MATH$')``.

    INPUT:
        x -- the object to render
        mode -- 'display' for displaymath or 'inline' for inline math

    OUTPUT:
        A string of html that contains the LaTeX represntation of x. In the
        notebook this gets embedded into the cell.

    EXAMPLES::

        sage: from sage.misc.latex import jsmath
        sage: f = maxima('1/(x^2+1)')
        sage: g = f.integrate()
        sage: jsmath(f)
        ... DeprecationWarning: The jsmath function is deprecated.  Use html('$math$') for inline mode or html('$$math$$') for display mode.
        # -*- coding: utf-8 -*-
        <html><font color='black'><div class="math">{{1}\over{x^2+1}}</div></font></html>
        sage: jsmath(g, 'inline')
        <html><font color='black'><span class="math">\tan^{-1} x</span></font></html>
        sage: jsmath('\int' + latex(f) + '\ dx=' + latex(g))
        <html><font color='black'><div class="math">\int{{1}\over{x^2+1}}\ dx=\tan^{-1} x</div></font></html>

    AUTHORS:

    - William Stein (2006-10): general layout (2006-10)

    - Bobby Moretti (2006-10): improvements, comments, documentation
    '''
    from sage.misc.misc import deprecation
    from sage.misc.html import html
    deprecation("The jsmath function is deprecated.  Use html('$math$') for inline mode or html('$$math$$') for display mode.")
    if mode == 'display':
        delimiter = '$$'
    elif mode == 'inline':
        delimiter = '$'
    else:
        raise ValueError, "mode must be either 'display' or 'inline'"
    try:
        # try to get a latex representation of the object
        x = x._latex_()
    except AttributeError:
        # otherwise just get the string representation
        x = str(x)
    return html(delimiter + x + delimiter)

def typeset(x):
    return JSMath().eval(x, mode='inline')

def view(objects, title='SAGE', debug=False, sep='', tiny=False,  **kwds):
    r"""nodetex
    Compute a latex representation of each object in objects, compile,
    and display typeset. If used from the command line, this requires
    that latex be installed.

    INPUT:


    -  ``objects`` - list (or object)

    -  ``title`` - string (default: 'Sage'): title for the
       document

    -  ``debug`` - bool (default: False): print verbose
       output

    -  ``sep`` - string (default: ''): separator between
       math objects

    -  ``tiny`` - bool (default: False): use tiny font.


    OUTPUT: Display typeset objects.

    This function behaves differently depending on whether in notebook
    mode or not.

    If not in notebook mode, this opens up a window displaying a dvi
    (or pdf) file, displaying the following: the title string is
    printed, centered, at the top. Beneath that, each object in objects
    is typeset on its own line, with the string sep inserted between
    these lines.

    The value of ``sep`` is inserted between each element of the list
    ``objects``; you can, for example, add vertical space between
    objects with ``sep='\\vspace{15mm}'``, while ``sep='\\hrule'``
    adds a horizontal line between objects, and ``sep='\\newpage'``
    inserts a page break between objects.

    If in notebook mode, this uses jmath to display the output in the
    notebook. Only the first argument, objects, is relevant; the others
    are ignored. If objects is a list, the result is typeset as a
    Python list, e.g. [12, -3.431] - each object in the list is not
    printed on its own line.

    EXAMPLES::

        sage: sage.misc.latex.EMBEDDED_MODE = True
        sage: view(3)
        <html><span class="math">\newcommand{\Bold}[1]{\mathbf{#1}}3</span></html>
        sage: sage.misc.latex.EMBEDDED_MODE = False
    """
    if EMBEDDED_MODE:
        print typeset(objects)
        return

    if isinstance(objects, LatexExpr):
        s = str(objects)
    else:
        s = _latex_file_(objects, title=title, debug=debug, sep=sep, tiny=tiny)
    from sage.misc.viewer import dvi_viewer
    viewer = dvi_viewer()
    tmp = tmp_dir('sage_viewer')
    open('%s/sage.tex'%tmp,'w').write(s)
    os.system('ln -sf %s/common/macros.tex %s'%(SAGE_DOC, tmp))
    O = open('%s/go'%tmp,'w')
    #O.write('export TEXINPUTS=%s/doc/commontex:.\n'%SAGE_ROOT)
    # O.write('latex \\\\nonstopmode \\\\input{sage.tex}; xdvi -noscan -offsets 0.3 -paper 100000x100000 -s %s sage.dvi ; rm sage.* macros.* go ; cd .. ; rmdir %s'%(zoom,tmp))

    # Added sleep 1 to allow viewer to open the file before it gets removed
    # Yi Qiang 2008-05-09
    O.write('latex \\\\nonstopmode \\\\input{sage.tex}; %s sage.dvi ; sleep 1 rm sage.* macros.* go ; cd .. ; rmdir %s' % (viewer, tmp))
    O.close()
    if not debug:
        direct = '1>/dev/null 2>/dev/null'
    else:
        direct = ''
    os.system('cd %s; chmod +x go; ./go %s&'%(tmp,direct))
    #return os.popen('cd %s; chmod +x go; ./go %s & '%(tmp,direct), 'r').read()


def png(x, filename, density=150, debug=False,
        do_in_background=True, tiny=False):
    """
    Create a png image representation of x and save to the given
    filename.

    INPUT:


    -  ``x`` - object to be displayed

    -  ``filename`` - file in which to save the image

    -  ``density`` - integer (default: 150)

    -  ``debug`` - bool (default: False): print verbose
       output

    -  ``do_in_background`` - bool (default: True): create the
       file in the background

    -  ``tiny`` - bool (default: False): use 'tiny' font

    """
    import sage.plot.all
    if sage.plot.all.is_Graphics(x):
        x.save(filename)
        return
    s = _latex_file_([x], math_left='$\\displaystyle', math_right='$', title='',
                     debug=debug, tiny=tiny,
                     extra_preamble='\\textheight=2\\textheight')
    abs_path_to_png = os.path.abspath(filename)

    tmp = tmp_dir('sage_viewer')
    open('%s/sage.tex'%tmp,'w').write(s)
    os.system('ln -sf %s/common/macros.tex %s'%(SAGE_DOC, tmp))
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
    if do_in_background:
        background = '&'
    else:
        background = ''
    os.system('cd %s; chmod +x go; ./go %s%s'%(tmp,direct,background))
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


    -  ``symbols`` - list of symbols

    -  ``coeffs`` - list of coefficients of the symbols


    OUTPUT:


    -  ``str`` - a string


    EXAMPLES::

        sage: t = PolynomialRing(QQ, 't').0
        sage: from sage.misc.latex import repr_lincomb
        sage: repr_lincomb(['a', 's', ''], [-t, t - 2, t^12 + 2])
        '-t\\text{a} + \\left(t - 2\\right)\\text{s} + \\left(t^{12} + 2\\right)\\text{}'
    """
    s = ""
    first = True
    i = 0

    all_atomic = True
    for c in coeffs:
        b = latex(symbols[i])
        if c != 0:
            if c == 1:
                s += b
            else:
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


def print_or_typeset(object):
    r"""
    'view' or 'print' the object depending on the situation.

    In particular, if in notebook mode with the typeset box checked,
    view the object. Otherwise, print the object.

    INPUT: object: anything

    EXAMPLES::

        sage: sage.misc.latex.print_or_typeset(3)
        3
        sage: sage.misc.latex.EMBEDDED_MODE=True
        sage: sage.misc.latex.print_or_typeset(3)
        3
        sage: TEMP = sys.displayhook
        sage: sys.displayhook = sage.misc.latex.pretty_print
        sage: sage.misc.latex.print_or_typeset(3)
        <html><span class="math">\newcommand{\Bold}[1]{\mathbf{#1}}3</span></html>
        sage: sage.misc.latex.EMBEDDED_MODE=False
        sage: sys.displayhook = TEMP
    """
    import sys
    if EMBEDDED_MODE and sys.displayhook == pretty_print:
        view(object)
    else:
        print(object)


def pretty_print (object):
    """
    Try to pretty print the object in an intelligent way. For many
    things, this will convert the object to latex inside of html and
    rely on a latex-aware front end (like jsMath) to render the text
    """
    if object is None:
        return
    import __builtin__
    __builtin__._=object

    from sage.plot.plot import Graphics
    from sage.plot.plot3d.base import Graphics3d
    if isinstance(object, (Graphics, Graphics3d)):
        print repr(object)
        return
    else:
        try:
            print typeset(object)
        except:
            import sys
            sys.__displayhook__(object)


def pretty_print_default(enable=True):
    """
    Enable or disable default pretty printing. Pretty printing means
    rendering things so that jsMath or some other latex-aware front end
    can render real math.
    """
    import sys
    if enable:
        sys.displayhook = pretty_print
    else:
        sys.displayhook = sys.__displayhook__



common_varnames = ['alpha',
                   'beta',
                   'gamma',
                   'Gamma',
                   'delta',
                   'Delta',
                   'epsilon',
                   'zeta',
                   'eta',
                   'theta',
                   'Theta',
                   'iota',
                   'kappa',
                   'lambda',
                   'Lambda',
                   'mu',
                   'nu',
                   'xi',
                   'Xi',
                   'pi',
                   'Pi',
                   'rho',
                   'sigma',
                   'Sigma',
                   'tau',
                   'upsilon',
                   'phi',
                   'Phi',
                   'varphi',
                   'chi',
                   'psi',
                   'Psi',
                   'omega',
                   'Omega']

def latex_varify(a):
    if a in common_varnames:
        return "\\" + a
    elif len(a) == 1:
        return a
    else:
        return '\\mbox{%s}'%a

def latex_variable_name(x):
    r"""
    Return latex version of a variable name.

    Here are some guiding principles for usage of this function:

    1. If the variable is a single letter, that is the latex version.

    2. If the variable name is suffixed by a number, we put the number
       in the subscript.

    3. If the variable name contains an '_' we start the subscript at
       the underscore. Note that #3 trumps rule #2.

    4. If a component of the variable is a greek letter, escape it
       properly.

    5. Recurse nicely with subscripts.

    Refer to the examples section for how these rules might play out in
    practice.

    EXAMPLES::

        sage: from sage.misc.latex import latex_variable_name
        sage: latex_variable_name('a')
        'a'
        sage: latex_variable_name('abc')
        '\\mbox{abc}'
        sage: latex_variable_name('sigma')
        '\\sigma'
        sage: latex_variable_name('sigma_k')
        '\\sigma_{k}'
        sage: latex_variable_name('sigma389')
        '\\sigma_{389}'
        sage: latex_variable_name('beta_00')
        '\\beta_{00}'
        sage: latex_variable_name('Omega84')
        '\\Omega_{84}'
        sage: latex_variable_name('sigma_alpha')
        '\\sigma_{\\alpha}'
        sage: latex_variable_name('nothing1')
        '\\mbox{nothing}_{1}'
        sage: latex_variable_name('nothing_abc')
        '\\mbox{nothing}_{\\mbox{abc}}'
        sage: latex_variable_name('alpha_beta_gamma12')
        '\\alpha_{\\beta_{\\gamma_{12}}}'

    AUTHORS:

    - Joel B. Mohler: drastic rewrite and many doc-tests
    """
    underscore = x.find("_")
    if underscore == -1:
        import re
        # * The "\d|[.,]" means "decimal digit" or period or comma
        # * The "+" means "1 or more"
        # * The "$" means "at the end of the line"
        m = re.search('(\d|[.,])+$',x)
        if m is None:
            prefix = x
            suffix = None
        else:
            prefix = x[:m.start()]
            suffix = x[m.start():]
    else:
        prefix = x[:underscore]
        suffix = x[underscore+1:]
    if suffix and len(suffix) > 0:
        # handle the suffix specially because it very well might be numeric
        # I use strip to avoid using regex's -- It makes it a bit faster (and the code is more comprehensible to non-regex'ed people)
        if suffix.strip("1234567890")!="":
            suffix = latex_variable_name(suffix) # recurse to deal with recursive subscripts
        return '%s_{%s}'%(latex_varify(prefix), suffix)
    else:
        return latex_varify(prefix)
