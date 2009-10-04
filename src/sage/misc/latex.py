"""
LaTeX printing support

In order to support latex formatting, an object should define a
special method ``_latex_(self)`` that returns a string.
"""

#*****************************************************************************
#
#   Sage: System for Algebra and Geometry Experimentation
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

EMBEDDED_MODE = False

COMMON_HEADER='\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\\usepackage{graphicx}\\usepackage{pstricks}\\pagestyle{empty}\\usepackage[utf8]{inputenc}\n'

LATEX_HEADER='\\documentclass{article}' + COMMON_HEADER + '\\oddsidemargin 0.0in\n\\evensidemargin 0.0in\n\\textwidth 6.45in\n\\topmargin 0.0in\n\\headheight 0.0in\n\\headsep 0.0in\n\\textheight 9.0in\n'

SLIDE_HEADER='\\documentclass[a0,8pt]{beamer}' + COMMON_HEADER + '\\textwidth=1.1\\textwidth\\textheight=2\\textheight\n'

#SLIDE_HEADER='\\documentclass[landscape]{slides}\\usepackage{fullpage}\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage{amsfonts}\\usepackage{graphicx}\usepackage{pstricks}\pagestyle{empty}\n'

import os, shutil

import os.path

import random

from misc import tmp_dir, graphics_filename
import sage_eval
from sage.misc.misc import SAGE_DOC

_have_latex = None
def have_latex():
    """
    Return True if this computer has the program latex.

    The first time it is run, this function caches its result in the
    variable ``_have_latex``, and any subsequent time, it just
    checks the value of the variable.

    If this computer doesn't have latex installed, you may obtain it
    from http://ctan.org/.

    EXAMPLES::

        sage: from sage.misc.latex import have_latex
        sage: have_latex() # random
        True
        sage: sage.misc.latex._have_latex is None
        False
        sage: sage.misc.latex._have_latex == have_latex()
        True
    """
    global _have_latex
    if _have_latex is None:
        _have_latex = not bool(os.system('which latex >/dev/null'))
    return _have_latex

_have_pdflatex = None
def have_pdflatex():
    """
    Return True if this computer has the program pdflatex.

    The first time it is run, this function caches its result in the
    variable ``_have_pdflatex``, and any subsequent time, it just
    checks the value of the variable.

    If this computer doesn't have pdflatex installed, you may obtain it
    from http://ctan.org/.

    EXAMPLES::

        sage: from sage.misc.latex import have_pdflatex
        sage: have_pdflatex() # random
        True
        sage: sage.misc.latex._have_pdflatex is None
        False
        sage: sage.misc.latex._have_pdflatex == have_pdflatex()
        True
    """
    global _have_pdflatex
    if _have_pdflatex is None:
        _have_pdflatex = not bool(os.system('which pdflatex >/dev/null'))
    return _have_pdflatex

_have_dvipng = None
def have_dvipng():
    """
    Return True if this computer has the program dvipng.

    The first time it is run, this function caches its result in the
    variable ``_have_dvipng``, and any subsequent time, it just
    checks the value of the variable.

    If this computer doesn't have dvipng installed, you may obtain it
    from http://sourceforge.net/projects/dvipng/

    EXAMPLES::

        sage: from sage.misc.latex import have_dvipng
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

_have_convert = None
def have_convert():
    """
    Return True if this computer has the program convert.

    The first time it is run, this function caches its result in the
    variable ``_have_convert``, and any subsequent time, it just
    checks the value of the variable.

    If this computer doesn't have convert installed, you may obtain it
    (along with the rest of the ImageMagick suite) from
    http://www.imagemagick.org

    EXAMPLES::

        sage: from sage.misc.latex import have_convert
        sage: have_convert() # random
        True
        sage: sage.misc.latex._have_convert is None
        False
        sage: sage.misc.latex._have_convert == have_convert()
        True
    """
    global _have_convert
    if _have_convert is None:
        _have_convert = not bool(os.system('which convert >/dev/null'))
    return _have_convert

def list_function(x):
    r"""
    Returns the LaTeX code for a list ``x``.

    INPUT: ``x`` - a list

    EXAMPLES::

        sage: from sage.misc.latex import list_function
        sage: list_function([1,2,3])
        '\\left[1, 2, 3\\right]'
        sage: latex([1,2,3])  # indirect doctest
        \left[1, 2, 3\right]
        sage: latex([Matrix(ZZ,3,range(9)), Matrix(ZZ,3,range(9))]) # indirect doctest
        \left[\left(\begin{array}{rrr}
        0 & 1 & 2 \\
        3 & 4 & 5 \\
        6 & 7 & 8
        \end{array}\right), \left(\begin{array}{rrr}
        0 & 1 & 2 \\
        3 & 4 & 5 \\
        6 & 7 & 8
        \end{array}\right)\right]
    """
    return "\\left[" + ", ".join([latex(v) for v in x]) + "\\right]"

def tuple_function(x):
    r"""
    Returns the LaTeX code for a tuple ``x``.

    INPUT: ``x`` - a tuple

    EXAMPLES::

        sage: from sage.misc.latex import tuple_function
        sage: tuple_function((1,2,3))
        '\\left(1, 2, 3\\right)'
    """
    return "\\left(" + ", ".join([latex(v) for v in x]) + "\\right)"

def bool_function(x):
    r"""
    Returns the LaTeX code for a boolean ``x``.

    INPUT: ``x`` - boolean

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
    Returns the LaTeX code for a string ``x``.

    INPUT: ``x`` - a string

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

def dict_function(x):
    r"""
    Returns the LaTeX code for a dictionary ``x``.

    INPUT: ``x`` - a dictionary

    EXAMPLES::

        sage: from sage.misc.latex import dict_function
        sage: x,y,z = var('x,y,z')
        sage: dict_function({x/2: y^2})
        '\\left\\{\\frac{1}{2} \\, x:\\: y^{2}\\right\\}'
        sage: d = {(1,2,x^2): [sin(z^2), y/2]}
        sage: latex(d)
        \left\{\left(1, 2, x^{2}\right):\: \left[\sin\left(z^{2}\right), \frac{1}{2} \, y\right]\right\}
    """
    return "\\left\\{" + ", ".join([latex(key) + ":\\: " + latex(value) for key,value in x.iteritems()]) + "\\right\\}"

# One can add to the latex_table in order to install latexing
# functionality for other types.  (Suggested by Robert Kerns of Enthought.)

latex_table = {list: list_function, tuple:tuple_function, bool:bool_function,
               str: str_function, int:str, long:str, float:str, dict: dict_function}

class LatexExpr(str):
    """
    LaTeX expression.  This is a string.

    EXAMPLES::

        sage: from sage.misc.latex import LatexExpr
        sage: LatexExpr(3)
        3
        sage: LatexExpr(False)
        False
        sage: LatexExpr(pi)
        pi
        sage: LatexExpr(3.14)
        3.14000000000000
        sage: LatexExpr("abc")
        abc
    """
    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.misc.latex import LatexExpr
            sage: LatexExpr("abc").__repr__()
            'abc'
        """
        return str(self)

    def _latex_(self):
        """
        EXAMPLES::

            sage: from sage.misc.latex import LatexExpr
            sage: LatexExpr("abc")._latex_()
            'abc'
        """
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
        self._option["macros"] = ""
        self._option["preamble"] = ""
        self._option["pdflatex"] = False
        self._option["jsmath_avoid"] = []

_Latex_prefs = _Latex_prefs_object()

##############################################################
# The Latex class is used to make slides and LaTeX output in
# the Sage Notebook
#########################################

def latex_extra_preamble():
    r"""
    Return the string containing the user-configured preamble,
    ``sage_latex_macros``, and any user-configured macros.  This is
    used in the :meth:`~Latex.eval` method for the :class:`Latex`
    class, and in :func:`_latex_file_`; it follows either
    ``LATEX_HEADER`` or ``SLIDE_HEADER`` (defined at the top of this
    file) which is a string containing the documentclass and standard
    usepackage commands.

    EXAMPLES::

        sage: from sage.misc.latex import latex_extra_preamble
        sage: latex_extra_preamble()
        '\n\\newcommand{\\ZZ}{\\Bold{Z}}\n\\newcommand{\\RR}{\\Bold{R}}\n\\newcommand{\\CC}{\\Bold{C}}\n\\newcommand{\\QQ}{\\Bold{Q}}\n\\newcommand{\\QQbar}{\\overline{\\QQ}}\n\\newcommand{\\GF}[1]{\\Bold{F}_{#1}}\n\\newcommand{\\Zp}[1]{\\ZZ_{#1}}\n\\newcommand{\\Qp}[1]{\\QQ_{#1}}\n\\newcommand{\\Zmod}[1]{\\ZZ/#1\\ZZ}\n\\newcommand{\\CDF}{\\text{Complex Double Field}}\n\\newcommand{\\CIF}{\\Bold{C}}\n\\newcommand{\\CLF}{\\Bold{C}}\n\\newcommand{\\RDF}{\\Bold{R}}\n\\newcommand{\\RIF}{\\Bold{I} \\Bold{R}}\n\\newcommand{\\RLF}{\\Bold{R}}\n\\newcommand{\\CFF}{\\Bold{CFF}}\n\\newcommand{\\Bold}[1]{\\mathbf{#1}}\n'
    """
    from sage.misc.latex_macros import sage_latex_macros
    return (_Latex_prefs._option['preamble'] + "\n"
                + "\n".join(sage_latex_macros) + "\n"
                + _Latex_prefs._option['macros'])

def _run_latex_(filename, debug=False, density=150,
                       pdflatex=None, png=False, do_in_background=False):
    """
    This runs LaTeX on the TeX file "filename.tex".  It produces files
    "filename.dvi" (or "filename.pdf"` if ``pdflatex`` is ``True``)
    and if ``png`` is True, "filename.png".  If ``png`` is True and
    dvipng can't convert the dvi file to png (because of postscript
    specials or other issues), then dvips is called, and the PS file
    is converted to a png file.

    INPUT:

    -  ``filename`` - string: file to process, including full path

    -  ``debug`` - bool (optional, default False): whether to print
       verbose debugging output

    -  ``density`` - integer (optional, default 150): how big output
       image is.

    -  ``pdflatex`` - bool (optional, default False): whether to use
       pdflatex.

    -  ``png`` - bool (optional, default False): whether to produce a
       png file.

    -  ``do_in_background`` - bool (optional, default False): whether
       to run in the background.

    OUTPUT: string, which could be a string starting with 'Error' (if
    there was a problem), or it could be 'pdf' or 'dvi'.  If
    ``pdflatex`` is False, then a dvi file is created, but if there
    appear to be problems with it (because of PS special commands, for
    example), then a pdf file is created instead.  The function
    returns 'dvi' or 'pdf' to indicate which type of file is created.
    (Detecting problems requires that dvipng be installed; if it is
    not, then the dvi file is not checked for problems and 'dvi' is
    returned.)  If ``pdflatex`` is True and there are no errors, then
    'pdf' is returned.

    .. warning::

       If ``png`` is True, then when using latex (the default), you
       must have 'dvipng' (or 'dvips' and 'convert') installed on your
       operating system, or this command won't work.  When using
       pdflatex, you must have 'convert' installed.

    EXAMPLES::

        sage: from sage.misc.latex import _run_latex_, _latex_file_
        sage: file = os.path.join(SAGE_TMP, "temp.tex")
        sage: O = open(file, 'w')
        sage: O.write(_latex_file_([ZZ[x], RR])); O.close()
        sage: _run_latex_(file) # random - depends on whether latex is installed
        'dvi'
    """
    if pdflatex is None:
        pdflatex = _Latex_prefs._option["pdflatex"]
    if not pdflatex and not have_latex():
        print "Error: LaTeX does not seem to be installed.  Download it from"
        print "ctan.org and try again."
        return "Error"
    if pdflatex and not have_pdflatex():
        print "Error: PDFLaTeX does not seem to be installed.  Download it from"
        print "ctan.org and try again."
        return "Error"
    # if png output + latex, check to see if dvipng or convert is installed.
    if png:
        if not pdflatex and not (have_dvipng() or have_convert()):
            print ""
            print "Error: neither dvipng nor convert (from the ImageMagick suite)"
            print "appear to be installed. Displaying LaTeX or PDFLaTeX output"
            print "requires at least one of these programs, so please install"
            print "and try again."
            print ""
            print "Go to http://sourceforge.net/projects/dvipng/ and"
            print "http://www.imagemagick.org to download these programs."
            return "Error"
    # if png output + pdflatex, check to see if convert is installed.
        elif pdflatex and not have_convert():
            print ""
            print "Error: convert (from the ImageMagick suite) does not"
            print "appear to be installed. Displaying PDFLaTeX output"
            print "requires this program, so please install and try again."
            print ""
            print "Go to http://www.imagemagick.org to download it."
            return "Error"
    # check_validity: check to see if the dvi file is okay by trying
    # to convert to a png file.  if this fails, return_suffix will be
    # set to "pdf".  return_suffix is the return value for this
    # function.
    #
    # thus if not png output, check validity of dvi output if dvipng
    # or convert is installed.
    else:
        check_validity = have_dvipng()
    # set up filenames, other strings:
    base, filename = os.path.split(filename)
    filename = os.path.splitext(filename)[0]  # get rid of extension
    if len(filename.split()) > 1:
        raise ValueError, "filename must contain no spaces"
    if not debug:
        redirect=' 2>/dev/null 1>/dev/null '
    else:
        redirect=''
    if do_in_background:
        background = ' &'
    else:
        background = ''
    if pdflatex:
        command = "pdflatex"
        # 'suffix' is used in the string 'convert' ...
        suffix = "pdf"
        return_suffix = "pdf"
    else:
        command = "latex"
        suffix = "ps"
        return_suffix = "dvi"
    # Define the commands to be used:
    lt = 'cd "%s"&& sage-native-execute %s \\\\nonstopmode \\\\input{%s.tex} %s'%(base, command, filename, redirect)
    # dvipng is run with the 'picky' option: this means that if
    # there are warnings, no png file is created.
    dvipng = 'cd "%s"&& sage-native-execute dvipng --picky -q -T tight -D %s %s.dvi -o %s.png'%(base, density, filename, filename)
    dvips = 'sage-native-execute dvips %s.dvi %s'%(filename, redirect)
    ps2pdf = 'sage-native-execute ps2pdf %s.ps %s'%(filename, redirect)
    # We seem to need a larger size when using convert compared to
    # when using dvipng:
    density = int(1.4 * density / 1.3)
    convert = 'sage-native-execute convert -density %sx%s -trim %s.%s %s.png %s '%\
        (density,density, filename, suffix, filename, redirect)

    e = 1 # it is possible to get through the following commands
          # without running a program, so in that case we force error
    if pdflatex:
        if png:
            cmd = ' && '.join([lt, convert])
        else:
            cmd = lt
        if debug:
            print cmd
        e = os.system(cmd + ' ' + redirect + background)
    else:  # latex, not pdflatex
        if (png or check_validity):
            if have_dvipng():
                cmd = ' && '.join([lt, dvipng])
                if debug:
                    print cmd
                e = os.system(cmd + ' ' + redirect)
                dvipng_error = not os.path.exists(base + '/' + filename + '.png')
                # If there is no png file, then either the latex
                # process failed or dvipng failed.  Assume that dvipng
                # failed, and try running dvips and convert.  (If the
                # latex process failed, then dvips and convert will
                # fail also, so we'll still catch the error.)
                if dvipng_error:
                    if png:
                        if have_convert():
                            cmd = ' && '.join(['cd "%s"'%(base,), dvips, convert])
                            if debug:
                                print "'dvipng' failed; trying 'convert' instead..."
                                print cmd
                            e = os.system(cmd + ' ' + redirect + background)
                        else:
                            print "Error: 'dvipng' failed and 'convert' is not installed."
                            return "Error: dvipng failed."
                    else:  # not png, i.e., check_validity
                        return_suffix = "pdf"
                        cmd = ' && '.join(['cd "%s"'%(base,), dvips, ps2pdf])
                        if debug:
                            print "bad dvi file; running dvips and ps2pdf instead..."
                            print cmd
                        e = os.system(cmd)
                        if e:  # error running dvips and/or ps2pdf
                            command = "pdflatex"
                            lt = 'cd "%s"&& sage-native-execute %s \\\\nonstopmode \\\\input{%s.tex} %s'%(base, command, filename, redirect)
                            if debug:
                                print "error running dvips and ps2pdf; trying pdflatex instead..."
                                print cmd
                            e = os.system(cmd + background)
            else:  # don't have dvipng, so must have convert.  run latex, dvips, convert.
                cmd = ' && '.join([lt, dvips, convert])
                if debug:
                    print cmd
                e = os.system(cmd + ' ' + redirect + background)
    if e:
        print "An error occurred."
        try:
            print open(base + '/' + filename + '.log').read()
        except IOError:
            pass
        return "Error latexing slide."
    return return_suffix

class Latex:
    r"""nodetex
    Enter, e.g.,

    ::

                %latex
                The equation $y^2 = x^3 + x$ defines an elliptic curve.
                We have $2006 = \sage{factor(2006)}$.


    in an input cell in the notebook to get a typeset version. Use
    ``%latex_debug`` to get debugging output.

    Use ``latex(...)`` to typeset a Sage object.

    Use ``%slide`` instead to typeset slides.

    .. warning::

       You must have dvipng (or dvips and convert) installed
       on your operating system, or this command won't work.

    """
    def __init__(self, debug=False, slide=False, density=150, pdflatex=None):
        self.__debug = debug
        self.__slide = slide
        self.__pdflatex = pdflatex
        self.__density = density

    def __call__(self, x):
        r"""
        Return a :class:`LatexExpr` built out of the argument ``x``.

        INPUT:

        - ``x`` - a Sage object

        OUTPUT: a LatexExpr built from ``x``

        EXAMPLES::

            sage: latex(Integer(3))  # indirect doctest
            3
            sage: latex(1==0)
            \mbox{\rm False}
            sage: print latex([x,2])
            \left[x, 2\right]
        """
        if hasattr(x, '_latex_'):
            return LatexExpr(x._latex_())

        try:
            f = latex_table[type(x)]
            return LatexExpr(f(x))

        except KeyError:
            if x is None:
                return LatexExpr("\\mbox{\\mathrm{None}}")

            return LatexExpr(str_function(str(x)))

    def _relation_symbols(self):
        """
        Returns a dictionary whose keys are attributes of the
        :mod:`operator` module and whose values are the corresponding
        LaTeX expressions.

        EXAMPLES::

            sage: import operator
            sage: latex._relation_symbols()[operator.ge]
            ' \\geq '
        """
        import operator
        return {operator.lt:' < ', operator.le:' \\leq ',
                operator.eq:' = ', operator.ne:' \\neq ',
                operator.ge:' \\geq ', operator.gt:' > '}


    def _latex_preparse(self, s, locals):
        r"""
        Replace instances of '\sage{x}' in s with the LaTeX version of
        x in the running session.

        EXAMPLES::

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
             density=None, pdflatex=None, locals={}):
        """
        INPUT:

        -  ``globals`` -- a globals dictionary

        -  ``x`` - string to evaluate.

        -  ``strip`` - ignored

        -  ``filename`` - output filename

        -  ``debug`` - whether to print verbose debugging
           output

        -  ``density`` - how big output image is.

        -  ``pdflatex`` - whether to use pdflatex.

        -  ``locals`` - extra local variables used when
           evaluating Sage code in ``x``.

        .. warning::

           When using latex (the default), you must have 'dvipng' (or
           'dvips' and 'convert') installed on your operating system,
           or this command won't work.  When using pdflatex, you must
           have 'convert' installed.
        """
        MACROS = latex_extra_preamble()

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
        O = open(os.path.join(base, filename + ".tex"), 'w')
        if self.__slide:
            O.write(SLIDE_HEADER)
            O.write(MACROS)
            O.write('\\begin{document}\n\n')
        else:
            O.write(LATEX_HEADER)
            O.write(MACROS)
            O.write('\\begin{document}\n')

        O.write(x.encode('utf-8'))
        if self.__slide:
            O.write('\n\n\\end{document}')
        else:
            O.write('\n\n\\end{document}\n')

        O.close()
        if not debug:
            redirect=' 2>/dev/null 1>/dev/null '
        else:
            redirect=''
        if pdflatex is None:
            if self.__pdflatex is None:
                pdflatex = _Latex_prefs._option["pdflatex"]
            else:
                pdflatex = bool(self.__pdflatex)
        e = _run_latex_(os.path.join(base, filename + ".tex"), debug=debug,
                               density=density, pdflatex=pdflatex, png=True)
        if e.find("Error") == -1:
            shutil.copy(os.path.join(base, filename + ".png"),
                        os.path.join(orig_base, filename + ".png"))
            shutil.rmtree(base)
            return ''
        else:
            return

    def blackboard_bold(self, t = None):
        r"""nodetex
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
        r"""nodetex
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
        r"""nodetex
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

    def extra_macros(self, macros=None):
        r"""nodetex
        String containing extra LaTeX macros to use with %latex,
        %html, and %jsmath.

        INPUT: ``macros`` - string

        If ``macros`` is None, return the current string.  Otherwise,
        set it to ``macros``.  If you want to *append* to the string
        of macros instead of replacing it, use :meth:`latex.add_macro <Latex.add_macro>`.

        EXAMPLES::

            sage: latex.extra_macros("\\newcommand{\\foo}{bar}")
            sage: latex.extra_macros()
            '\\newcommand{\\foo}{bar}'
            sage: latex.extra_macros("")
            sage: latex.extra_macros()
            ''
        """
        if macros is None:
            return _Latex_prefs._option['macros']
        else:
            _Latex_prefs._option['macros'] = macros

    def add_macro(self, macro):
        r"""nodetex
        Append to the string of extra LaTeX macros, for use with
        %latex, %html, and %jsmath.

        INPUT: ``macro`` - string

        EXAMPLES::

            sage: latex.extra_macros()
            ''
            sage: latex.add_macro("\\newcommand{\\foo}{bar}")
            sage: latex.extra_macros()
            '\\newcommand{\\foo}{bar}'
            sage: latex.extra_macros("")  # restore to default
        """
        current = latex.extra_macros()
        if current.find(macro) == -1:
            _Latex_prefs._option['macros'] += macro

    def extra_preamble(self, s=None):
        r"""nodetex
        String containing extra preamble to be used with %latex.
        Anything in this string won't be processed by %jsmath.

        INPUT: ``s`` - string or ``None``

        If ``s`` is None, return the current preamble.  Otherwise, set
        it to ``s``.  If you want to *append* to the current extra
        preamble instead of replacing it, use
        :meth:`latex.add_to_preamble <Latex.add_to_preamble>`.

        EXAMPLES::

            sage: latex.extra_preamble("\\DeclareMathOperator{\\Ext}{Ext}")
            sage: latex.extra_preamble()
            '\\DeclareMathOperator{\\Ext}{Ext}'
            sage: latex.extra_preamble("")
            sage: latex.extra_preamble()
            ''
        """
        if s is None:
            return _Latex_prefs._option['preamble']
        else:
            _Latex_prefs._option['preamble'] = s

    def add_to_preamble(self, s):
        r"""nodetex
        Append to the string of extra LaTeX macros, for use with
        %latex.  Anything in this string won't be processed by
        %jsmath.

        EXAMPLES::

            sage: latex.extra_preamble()
            ''
            sage: latex.add_to_preamble("\\DeclareMathOperator{\\Ext}{Ext}")

        At this point, a notebook cell containing

        ::

          %latex
          $\Ext_A^*(\GF{2}, \GF{2}) \Rightarrow \pi_*^s*(S^0)$

        will be typeset correctly.

        ::

            sage: latex.add_to_preamble("\\usepackage{xypic}")
            sage: latex.extra_preamble()
            '\\DeclareMathOperator{\\Ext}{Ext}\\usepackage{xypic}'

        Now one can put various xypic diagrams into a %latex cell, such as

        ::

          %latex
          \[ \xymatrix{ \circ \ar `r[d]^{a} `[rr]^{b} `/4pt[rr]^{c} `[rrr]^{d}
          `_dl[drrr]^{e} [drrr]^{f} & \circ & \circ & \circ \\ \circ & \circ &
          \circ & \circ } \]

        Reset the preamble to its default, the empty string::

            sage: latex.extra_preamble('')
            sage: latex.extra_preamble()
            ''
        """
        current = latex.extra_preamble()
        if current.find(s) == -1:
            _Latex_prefs._option['preamble'] += s

    def jsmath_avoid_list(self, L=None):
        r"""nodetex
        List of strings which signal that jsMath should not
        be used when 'view'ing.

        INPUT: ``L`` - list or ``None``

        If ``L`` is ``None``, then return the current list.
        Otherwise, set it to ``L``.  If you want to *append* to the
        current list instead of replacing it, use
        :meth:`latex.add_to_jsmath_avoid_list <Latex.add_to_jsmath_avoid_list>`.

        EXAMPLES::

            sage: latex.jsmath_avoid_list(["\\mathsf", "pspicture"])
            sage: latex.jsmath_avoid_list()  # display current setting
            ['\\mathsf', 'pspicture']
            sage: latex.jsmath_avoid_list([])  # reset to default
            sage: latex.jsmath_avoid_list()
            []
        """
        if L is None:
            return _Latex_prefs._option['jsmath_avoid']
        else:
            _Latex_prefs._option['jsmath_avoid'] = L

    def add_to_jsmath_avoid_list(self, s):
        r"""nodetex
        Add to the list of strings which signal that jsMath should not
        be used when 'view'ing.

        INPUT: ``s`` - string -- add ``s`` to the list of 'jsMath avoid' strings

        If you want to replace the current list instead of adding to
        it, use :meth:`latex.jsmath_avoid_list <Latex.jsmath_avoid_list>`.

        EXAMPLES::

            sage: latex.add_to_jsmath_avoid_list("\\mathsf")
            sage: latex.jsmath_avoid_list()  # display current setting
            ['\\mathsf']
            sage: latex.add_to_jsmath_avoid_list("tkz-graph")
            sage: latex.jsmath_avoid_list()  # display current setting
            ['\\mathsf', 'tkz-graph']
            sage: latex.jsmath_avoid_list([])  # reset to default
            sage: latex.jsmath_avoid_list()
            []
        """
        current = latex.jsmath_avoid_list()
        if s not in current:
            _Latex_prefs._option['jsmath_avoid'].append(s)

    def pdflatex(self, t = None):
        """
        Controls whether Sage uses PDFLaTeX or LaTeX when typesetting
        with :func:`view`, in ``%latex`` cells, etc.

        INPUT:

        - ``t`` -- boolean or None

        OUTPUT: if t is None, return the current setting (True or False).

        If t == True, use PDFLaTeX; otherwise use LaTeX.

        EXAMPLES::

            sage: latex.pdflatex()
            False
            sage: latex.pdflatex(True)
            sage: latex.pdflatex()
            True
            sage: latex.pdflatex(False)
        """
        if t is None:
            return _Latex_prefs._option["pdflatex"]
        _Latex_prefs._option["pdflatex"] = bool(t)

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
        '\\documentclass{article}...\\begin{document}\n\\begin{center}{\\Large\\bf The number three}\\end{center}\n\\vspace{40mm}\\[3\\]\n\\end{document}'
        sage: _latex_file_([7, 8, 9], title="Why was six afraid of seven?", sep='\\vfill\\hrule\\vfill')
        '\\documentclass{article}...\\begin{document}\n\\begin{center}{\\Large\\bf Why was six afraid of seven?}\\end{center}\n\\vspace{40mm}\\[7\\]\n\n\\vfill\\hrule\\vfill\n\n\\[8\\]\n\n\\vfill\\hrule\\vfill\n\n\\[9\\]\n\\end{document}'
    """
    MACROS = latex_extra_preamble()

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
    """
    An arbitrary JSMath expression that can be nicely concatenated.

    EXAMPLES::

        sage: from sage.misc.latex import JSMathExpr
        sage: JSMathExpr("a^{2}") + JSMathExpr("x^{-1}")
        a^{2}x^{-1}
    """
    def __init__(self, y):
        """
        Initialize a JSMath expression.

        INPUT:

        - ``y`` - a string

        Note that no error checking is done on the type of ``y``.

        EXAMPLES::

            sage: from sage.misc.latex import JSMathExpr
            sage: js = JSMathExpr(3); js  # indirect doctest
            3
        """
        self.__y = y

    def __repr__(self):
        """
        Print representation.

        EXAMPLES::

            sage: from sage.misc.latex import JSMathExpr
            sage: js = JSMathExpr('3')
            sage: js.__repr__()
            '3'
        """
        return str(self.__y)

    def __add__(self, y):
        """
        'Add' JSMathExpr ``self`` to ``y``.  This concatenates them
        (assuming that they're strings).

        EXAMPLES::

            sage: from sage.misc.latex import JSMathExpr
            sage: j3 = JSMathExpr('3')
            sage: jx = JSMathExpr('x')
            sage: j3 + jx
            3x
        """
        return JSMathExpr(self.__y + y)

    def __radd__(self, y):
        """
        'Add' JSMathExpr ``y`` to ``self``.  This concatenates them
        (assuming that they're strings).

        EXAMPLES::

            sage: from sage.misc.latex import JSMathExpr
            sage: j3 = JSMathExpr('3')
            sage: jx = JSMathExpr('x')
            sage: j3.__radd__(jx)
            x3
        """
        return JSMathExpr(y + self.__y)

class JSMath:
    r"""
    Render LaTeX input using JSMath.  This returns a :class:`JSMathExpr`.

    EXAMPLES::

        sage: from sage.misc.latex import JSMath
        sage: JSMath()(3)
        <html><div class="math">\newcommand{\Bold}[1]{\mathbf{#1}}3</div></html>
        sage: JSMath()(ZZ)
        <html><div class="math">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}</div></html>
    """

    def __call__(self, x):
        r"""
        Render LaTeX input using JSMath.  This returns a :class:`JSMathExpr`.

        INPUT:

        - ``x`` - a Sage object

        OUTPUT: a JSMathExpr

        EXAMPLES::

            sage: from sage.misc.latex import JSMath
            sage: JSMath()(3)
            <html><div class="math">\newcommand{\Bold}[1]{\mathbf{#1}}3</div></html>
            sage: str(JSMath().eval(ZZ[x], mode='display')) == str(JSMath()(ZZ[x]))
            True
        """
        return self.eval(x)

    def eval(self, x, globals=None, locals=None, mode='display'):
        r"""
        Render LaTeX input using JSMath.  This returns a :class:`JSMathExpr`.

        INPUT:

        - ``x`` - a Sage object

        -  ``globals`` -- a globals dictionary

        -  ``locals`` - extra local variables used when
           evaluating Sage code in ``x``.

        -  ``mode`` - string (optional, default 'display): 'display'
           for displaymath or 'inline' for inline math

        OUTPUT: a JSMathExpr

        EXAMPLES::

            sage: from sage.misc.latex import JSMath
            sage: JSMath().eval(3, mode='display')
            <html><div class="math">\newcommand{\Bold}[1]{\mathbf{#1}}3</div></html>
            sage: JSMath().eval(3, mode='inline')
            <html><span class="math">\newcommand{\Bold}[1]{\mathbf{#1}}3</span></html>
        """
        # try to get a latex representation of the object
        if hasattr(x, '_latex_'):
            x = LatexExpr(x._latex_())
        else:
            try:
                f = latex_table[type(x)]
                x = LatexExpr(f(x))
            except KeyError:
                # otherwise just get the string representation
                x = str(x)

        # in JSMath:
        # inline math: <span class="math">...</span>
        # displaymath: <div class="math">...</div>
        from sage.misc.latex_macros import sage_configurable_latex_macros
        if 'display' == mode:
            return JSMathExpr('<html><div class="math">'
                              + ''.join(sage_configurable_latex_macros)
                              + _Latex_prefs._option['macros']
                              + '%s</div></html>'%x)
        elif 'inline' == mode:
            return JSMathExpr('<html><span class="math">'
                              + ''.join(sage_configurable_latex_macros)
                              + _Latex_prefs._option['macros']
                              + '%s</span></html>'%x)
        else:
            # what happened here?
            raise ValueError, "mode must be either 'display' or 'inline'"

def jsmath(x, mode='display'):
    r"""
    Attempt to nicely render an arbitrary Sage object with jsMath typesetting.
    Tries to call ._latex_() on x. If that fails, it will render a string
    representation of x.

    .. warning::

        2009-04: This function is deprecated; use :func:`html`
        instead: replace ``jsmath('MATH', mode='display')`` with
        ``html('$$MATH$$')``, and replace ``jsmath('MATH',
        mode='inline')`` with ``html('$MATH$')``.

    INPUT:
        x -- the object to render
        mode -- 'display' for displaymath or 'inline' for inline math

    OUTPUT:
        A string of html that contains the LaTeX representation of x. In the
        notebook this gets embedded into the cell.

    EXAMPLES::

        sage: from sage.misc.latex import jsmath
        sage: f = maxima('1/(x^2+1)')
        sage: g = f.integrate()
        sage: jsmath(f)
        doctest:1: DeprecationWarning: The jsmath function is deprecated.  Use html('$math$') for inline mode or html('$$math$$') for display mode.
        <html><font color='black'><div class="math">{{1}\over{x^2+1}}</div></font></html>
        sage: jsmath(g, 'inline')
        <html><font color='black'><span class="math">\tan^{-1} x</span></font></html>
        sage: jsmath('\int' + latex(f) + '\ dx=' + latex(g))
        <html><font color='black'><div class="math">\int{{1}\over{x^2+1}}\ dx=\tan^{-1} x</div></font></html>

    AUTHORS:

    - William Stein (2006-10): general layout (2006-10)

    - Bobby Moretti (2006-10): improvements, comments, documentation
    """
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

def view(objects, title='SAGE', debug=False, sep='', tiny=False, pdflatex=None, viewer = None, tightpage = None, **kwds):
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

    -  ``pdflatex`` - bool (default: False): use pdflatex.

    -  ``viewer`` -- string or None (default: None): specify a viewer to
       use; currently the only options are None and 'pdf'.

    -  ``tightpage`` - bool (default: False): use the LaTeX package
       'preview' with the 'tightpage' option.

    OUTPUT: Display typeset objects.

    This function behaves differently depending on whether in notebook
    mode or not.

    If not in notebook mode, the output is displayed in a separate
    viewer displaying a dvi (or pdf) file, with the following: the
    title string is printed, centered, at the top. Beneath that, each
    object in ``objects`` is typeset on its own line, with the string
    ``sep`` inserted between these lines.

    The value of ``sep`` is inserted between each element of the list
    ``objects``; you can, for example, add vertical space between
    objects with ``sep='\\vspace{15mm}'``, while ``sep='\\hrule'``
    adds a horizontal line between objects, and ``sep='\\newpage'``
    inserts a page break between objects.

    If ``pdflatex`` is ``True``, then this produces a pdf file.
    Otherwise, it produces a dvi file, and if the program dvipng is
    installed, it checks the dvi file by trying to convert it to a png
    file.  If this conversion fails, the dvi file probably contains
    some postscript special commands or it has other issues which
    might make displaying it a problem; in this case, the file is
    converted to a pdf file, which is then displayed.

    Setting ``viewer`` to ``'pdf'`` forces the use of a separate
    viewer, even in notebook mode. This also sets ``pdflatex`` to
    ``True``.

    Setting the option ``tightpage`` to ``True`` tells LaTeX to use
    the package 'preview' with the 'tightpage' option. Then, each
    object is typeset in its own page, and that page is cropped to
    exactly the size of the object. This is typically useful for very
    large pictures (like graphs) generated with tikz. This only works
    when using a separate viewer. Note that the object are currently
    typeset in plain math mode rather than displaymath, because the
    latter imposes a limit on the width of the picture. Technically,
    ``tightpage`` adds ::

      \\usepackage[tightpage,active]{preview}
      \\PreviewEnvironment{page}

    to the LaTeX preamble, and replaces the ``\\[`` and ``\\]`` around
    each object by ``\\begin{page}$`` and ``$\\end{page}``.

    If in notebook mode with ``viewer`` equal to ``None``, this
    usually uses jsMath -- see the next paragraph for the exception --
    to display the output in the notebook. Only the first argument,
    ``objects``, is relevant; the others are ignored. If ``objects``
    is a list, each object is printed on its own line.

    In the notebook, this *does* *not* use jsMath if the LaTeX code
    for ``objects`` contains a string in
    :meth:`latex.jsmath_avoid_list() <Latex.jsmath_avoid_list>`.  In
    this case, it creates and displays a png file.

    EXAMPLES::

        sage: sage.misc.latex.EMBEDDED_MODE = True
        sage: view(3)
        <html><span class="math">\newcommand{\Bold}[1]{\mathbf{#1}}3</span></html>
        sage: sage.misc.latex.EMBEDDED_MODE = False
    """
    if isinstance(objects, LatexExpr):
        s = str(objects)
    else:
        if tightpage == True:
            latex_options = {'extra_preamble':'\\usepackage[tightpage,active]{preview}\\PreviewEnvironment{page}',
                             'math_left':'\\begin{page}$', 'math_right':'$\\end{page}'}
        else:
            latex_options = {}
        s = _latex_file_(objects, title=title, sep=sep, tiny=tiny, debug=debug, **latex_options)
    # notebook
    if EMBEDDED_MODE and viewer is None:
        jsMath_okay = True
        for t in latex.jsmath_avoid_list():
            if s.find(t) != -1:
                jsMath_okay = False
            if not jsMath_okay:
                break
        if jsMath_okay:
            print JSMath().eval(objects, mode='inline')  # put comma at end of line?
        else:
            base_dir = os.path.abspath("")
            png_file = graphics_filename(ext='png')
            png_link = "cell://" + png_file
            png(objects, os.path.join(base_dir, png_file),
                debug=debug, do_in_background=False, pdflatex=pdflatex)
            print '<html><img src="%s"></html>'%png_link  # put comma at end of line?
        return
    # command line
    if viewer == "pdf":
        pdflatex = True
    if pdflatex is None:
        pdflatex = _Latex_prefs._option["pdflatex"]
    tmp = tmp_dir('sage_viewer')
    tex_file = os.path.join(tmp, "sage.tex")
    open(tex_file,'w').write(s)
    suffix = _run_latex_(tex_file, debug=debug, pdflatex=pdflatex, png=False)
    if suffix == "pdf":
        from sage.misc.viewer import pdf_viewer
        viewer = pdf_viewer()
    elif suffix == "dvi":
        from sage.misc.viewer import dvi_viewer
        viewer = dvi_viewer()
    else:
        print "Latex error"
        return
    output_file = os.path.join(tmp, "sage." + suffix)
    os.system('sage-native-execute %s %s'%(viewer, output_file))
    return


def png(x, filename, density=150, debug=False,
        do_in_background=False, tiny=False, pdflatex=True):
    """
    Create a png image representation of ``x`` and save to the given
    filename.

    INPUT:

    -  ``x`` - object to be displayed

    -  ``filename`` - file in which to save the image

    -  ``density`` - integer (default: 150)

    -  ``debug`` - bool (default: False): print verbose
       output

    -  ``do_in_background`` - bool (default: False): create the
       file in the background.

    -  ``tiny`` - bool (default: False): use 'tiny' font

    -  ``pdflatex`` - bool (default: False): use pdflatex.

    EXAMPLES::

        sage: from sage.misc.latex import png
        sage: png(ZZ[x], SAGE_TMP + "zz.png", do_in_background=False) # random - error if no latex
    """
    import sage.plot.all
    if sage.plot.all.is_Graphics(x):
        x.save(filename)
        return
    # if not graphics: create a string of latex code to write in a file
    s = _latex_file_([x], math_left='$\\displaystyle', math_right='$', title='',
                     debug=debug, tiny=tiny,
                     extra_preamble='\\textheight=2\\textheight')
    # path name for permanent png output
    abs_path_to_png = os.path.abspath(filename)
    # temporary directory to store stuff
    tmp = tmp_dir('sage_viewer')
    tex_file = os.path.join(tmp, "sage.tex")
    png_file = os.path.join(tmp, "sage.png")
    # write latex string to file
    open(tex_file,'w').write(s)
    # run latex on the file, producing png output to png_file
    e = _run_latex_(tex_file, density=density, debug=debug,
                    png=True, do_in_background=do_in_background,
                    pdflatex=pdflatex)
    if e.find("Error") == -1:
        # if no errors, copy png_file to the appropriate place
        shutil.copy(png_file, abs_path_to_png)
    else:
        print "Latex error"
    if debug:
        return s
    return

def coeff_repr(c):
    r"""
    LaTeX string representing coefficients in a linear combination.

    INPUT:

    - ``c`` - a coefficient (i.e., an element of a ring)

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.latex import coeff_repr
        sage: coeff_repr(QQ(1/2))
        '\\frac{1}{2}'
        sage: coeff_repr(-x^2)
        '\\left(-x^{2}\\right)'
    """
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
    r"""
    Compute a latex representation of a linear combination of some
    formal symbols.

    INPUT:

    -  ``symbols`` - list of symbols

    -  ``coeffs`` - list of coefficients of the symbols

    OUTPUT: a string

    EXAMPLES::

        sage: t = PolynomialRing(QQ, 't').0
        sage: from sage.misc.latex import repr_lincomb
        sage: repr_lincomb(['a', 's', ''], [-t, t - 2, t^12 + 2])
        '-t\\text{a} + \\left(t - 2\\right)\\text{s} + \\left(t^{12} + 2\\right)\\text{}'
        sage: repr_lincomb(['a', 'b'], [1,1])
        '\\text{a} + \\text{b}'

    Verify that a certain corner case works (see trac 5707 and 5766)::

        sage: repr_lincomb([1,5,-3],[2,8/9,7])
        '2\\cdot 1 + \\frac{8}{9}\\cdot 5 + 7\\cdot -3'
    """
    s = ""
    first = True
    i = 0

    from sage.rings.all import CC

    for c in coeffs:
        bv = symbols[i]
        b = latex(bv)
        if c != 0:
            if c == 1:
                if first:
                    s += b
                else:
                    s += " + %s"%b
            else:
                coeff = coeff_repr(c)
                if first:
                    coeff = str(coeff)
                else:
                    coeff = " + %s"%coeff
                # this is a hack: i want to say that if the symbol
                # happens to be a number, then we should put a
                # multiplication sign in
                try:
                    if bv in CC:
                        s += "%s\cdot %s"%(coeff, b)
                except:
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
    r"""
    Try to pretty print an object in an intelligent way.  For graphics
    objects, this returns their default representation.  For other
    objects, in the notebook, this calls the :func:`view` command,
    while from the command line, this produces an html string suitable
    for processing by jsMath.

    INPUT:

    - ``object`` - a Sage object

    This function is used in the notebook when the ``Typeset`` button is
    checked.

    EXAMPLES::

        sage: from sage.misc.latex import pretty_print
        sage: pretty_print(ZZ)  # indirect doctest
        <html><span class="math">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}</span></html>
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
        if EMBEDDED_MODE:
            view(object)
        else:
            print JSMath().eval(object, mode='inline')
        return


def pretty_print_default(enable=True):
    r"""
    Enable or disable default pretty printing. Pretty printing means
    rendering things so that jsMath or some other latex-aware front end
    can render real math.

    INPUT:

    -  ``enable`` - bool (optional, default True).  If True, turn on
       pretty printing; if False, turn it off.

    EXAMPLES::

        sage: pretty_print_default(True)
        sage: sys.displayhook
        <html><span class="math">...<function pretty_print at ...></span></html>
        sage: pretty_print_default(False)
        sage: sys.displayhook == sys.__displayhook__
        True
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

def latex_varify(a, is_fname=False):
    r"""
    Convert a string ``a`` to a LaTeX string: if it's an element of
    ``common_varnames``, then prepend a backslash.  If ``a`` consists
    of a single letter, then return it.  Otherwise, return
    either "{\\rm a}" or "\\mbox{a}" if "is_fname" flag is True or False.

    INPUT:

    - ``a`` - string

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.latex import latex_varify
        sage: latex_varify('w')
        'w'
        sage: latex_varify('aleph')
        '\\mbox{aleph}'
        sage: latex_varify('aleph', is_fname=True)
        '{\\rm aleph}'
        sage: latex_varify('alpha')
        '\\alpha'
    """
    if a in common_varnames:
        return "\\" + a
    elif len(a) == 1:
        return a
    elif is_fname is True:
        return '{\\rm %s}'%a
    else:
        return '\\mbox{%s}'%a

def latex_variable_name(x, is_fname=False):
    r"""
    Return latex version of a variable name.

    Here are some guiding principles for usage of this function:

    1. If the variable is a single letter, that is the latex version.

    2. If the variable name is suffixed by a number, we put the number
       in the subscript.

    3. If the variable name contains an '_' we start the subscript at
       the underscore. Note that #3 trumps rule #2.

    4. If a component of the variable is a Greek letter, escape it
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
        sage: latex_variable_name('nothing1', is_fname=True)
        '{\\rm nothing}_{1}'
        sage: latex_variable_name('nothing_abc')
        '\\mbox{nothing}_{\\mbox{abc}}'
        sage: latex_variable_name('nothing_abc', is_fname=True)
        '{\\rm nothing}_{{\\rm abc}}'
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
            suffix = latex_variable_name(suffix, is_fname) # recurse to deal with recursive subscripts
        return '%s_{%s}'%(latex_varify(prefix, is_fname), suffix)
    else:
        return latex_varify(prefix, is_fname)

class LatexExamples():
    r"""
    A catalogue of Sage objects with complicated ``_latex_`` methods.
    Use these for testing :func:`latex`, :func:`view`, the Typeset
    button in the notebook, etc.

    The classes here only have ``__init__``, ``_repr_``, and ``_latex_`` methods.

    EXAMPLES::

        sage: from sage.misc.latex import latex_examples
        sage: K = latex_examples.knot()
        sage: K
        LaTeX example for testing display of a knot produced by xypic...
        sage: latex(K)
        \vtop{\vbox{\xygraph{!{0;/r1.5pc/:}
        [u] !{\vloop<(-.005)\khole||\vcrossneg \vunder- }
        [] !{\ar @{-}@'{p-(1,0)@+}+(-1,1)}
        [ul] !{\vcap[3]>\khole}
        [rrr] !{\ar @{-}@'{p-(0,1)@+}-(1,1)}
        }}}
    """
    class graph(SageObject):
        """
        LaTeX example for testing display of graphs.  See its string
        representation for details.

        EXAMPLES::

            sage: from sage.misc.latex import latex_examples
            sage: G = latex_examples.graph()
            sage: G
            LaTeX example for testing display of graphs...
        """

        def __init__(self):
            """
            See the string representation for complete documentation.

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: type(latex_examples.graph())
                <class 'sage.misc.latex.graph'>
            """
            pass

        def _repr_(self):
            """
            String representation

            EXAMPLES:

                sage: from sage.misc.latex import latex_examples
                sage: G = latex_examples.graph()
                sage: len(G._repr_()) > 300
                True
            """
            return r"""LaTeX example for testing display of graphs.

To use, first try calling 'view' on this object -- you should get
gibberish.  Now, make sure that you have the most recent version of
the TeX package pgf installed, along with the LaTeX package tkz-graph.
Run 'latex.add_to_preamble("\\usepackage{tkz-graph}")', and try viewing
it again.  From the command line, this should pop open a nice window
with a picture of a graph.  In the notebook, you should get a jsMath
error.  Finally, run 'latex.add_to_jsmath_avoid_list("tikzpicture")'
and try again from the notebook -- you should get a nice picture.

(LaTeX code taken from http://altermundus.com/pages/graph.html)
"""

        def _latex_(self):
            """
            LaTeX representation

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: len(latex_examples.graph()._latex_()) > 500
                True
                sage: len(latex_examples.graph()._latex_()) > 600
                False
            """
            return r"""\begin{tikzpicture}[node distance   = 4 cm]
   \GraphInit[vstyle=Shade]
   \tikzset{LabelStyle/.style =   {draw,
                                   fill  = yellow,
                                   text  = red}}
   \Vertex{A}
   \EA(A){B}
   \EA(B){C}
   \tikzset{node distance   = 8 cm}
   \NO(B){D}
   \Edge[label=1](B)(D)
   \tikzset{EdgeStyle/.append style = {bend left}}
   \Edge[label=4](A)(B)
   \Edge[label=5](B)(A)
   \Edge[label=6](B)(C)
   \Edge[label=7](C)(B)
   \Edge[label=2](A)(D)
   \Edge[label=3](D)(C)
\end{tikzpicture}"""

    class pstricks(SageObject):
        """
        LaTeX example for testing display of pstricks output.  See its
        string representation for details.

        EXAMPLES::

            sage: from sage.misc.latex import latex_examples
            sage: PS = latex_examples.pstricks()
            sage: PS
            LaTeX example for testing display of pstricks...
        """
        def __init__(self):
            """
            See the string representation for complete documentation.

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: type(latex_examples.pstricks())
                <class 'sage.misc.latex.pstricks'>
            """
            pass

        def _repr_(self):
            """
            String representation

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: len(latex_examples.pstricks()._repr_()) > 300
                True
            """
            return """LaTeX example for testing display of pstricks output.

To use, view this object.  It should work okay from the command
line but not from the notebook.  In the notebook, run
'latex.add_to_jsmath_avoid_list("pspicture")' and try again -- you
should get a picture (of forces acting on a mass on a pendulum)."""

        def _latex_(self):
            """
            LaTeX representation

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: len(latex_examples.pstricks()._latex_()) > 250
                True
            """
            return r"""\begin{pspicture}(0,-4)(14,0)
  \psline{-}(0,0)(0,-4)
  \psline[linewidth=2pt]{-}(0,0)(1,-3)
  \qdisk(1,-3){3pt}
  \psarc{-}(0,0){0.6}{270}{292}
  \psline{->}(1,-3.3)(1,-4)
  \psline{->}(1.1,-2.7)(0.85,-1.95)
  \psline{-}(5,0)(5,-4)
  \psline[linewidth=2pt]{-}(5,0)(6,-3)
  \qdisk(6,-3){3pt}
  \psarc{-}(5,0){0.6}{270}{292}
  \psarc{-}(5,0){3.2}{270}{290}
\end{pspicture}"""

    class knot(SageObject):
        """
        LaTeX example for testing display of knots.  See its string
        representation for details.

        EXAMPLES::

            sage: from sage.misc.latex import latex_examples
            sage: K = latex_examples.knot()
            sage: K
            LaTeX example for testing display of a knot...
        """
        def __init__(self):
            """
            See the string representation for complete documentation.

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: type(latex_examples.knot())
                <class 'sage.misc.latex.knot'>
            """
            pass

        def _repr_(self):
            """
            String representation

            EXAMPLES:

                sage: from sage.misc.latex import latex_examples
                sage: len(latex_examples.knot()._repr_()) > 300
                True
            """
            return r"""LaTeX example for testing display of a knot produced by xypic.

To use, try to view this object -- it won't work.  Now try
'latex.add_to_preamble("\\usepackage[graph,knot,poly,curve]{xypic}")',
and try viewing again -- it should work in the command line but not
from the notebook.  In the notebook, run
'latex.add_to_jsmath_avoid_list("xygraph")' and try again -- you
should get a nice picture.

(LaTeX code taken from the xypic manual)
"""

        def _latex_(self):
            """
            LaTeX representation

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: len(latex_examples.knot()._latex_()) > 180
                True
            """
            return r"""\vtop{\vbox{\xygraph{!{0;/r1.5pc/:}
 [u] !{\vloop<(-.005)\khole||\vcrossneg \vunder- }
 [] !{\ar @{-}@'{p-(1,0)@+}+(-1,1)}
 [ul] !{\vcap[3]>\khole}
 [rrr] !{\ar @{-}@'{p-(0,1)@+}-(1,1)}
}}}"""

    class diagram(SageObject):
        """
        LaTeX example for testing display of commutative diagrams.
        See its string representation for details.

        EXAMPLES::

            sage: from sage.misc.latex import latex_examples
            sage: CD = latex_examples.diagram()
            sage: CD
            LaTeX example for testing display of a commutative diagram...
        """
        def __init__(self):
            """
            See the string representation for complete documentation.

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: type(latex_examples.diagram())
                <class 'sage.misc.latex.diagram'>
            """
            pass

        def _repr_(self):
            """
            String representation

            EXAMPLES:

                sage: from sage.misc.latex import latex_examples
                sage: len(latex_examples.diagram()._repr_()) > 300
                True
            """
            return r"""LaTeX example for testing display of a commutative diagram produced
by xypic.

To use, try to view this object -- it won't work.  Now try
'latex.add_to_preamble("\\usepackage[matrix,arrow,curve,cmtip]{xy}")',
and try viewing again -- it should work in the command line but not
from the notebook.  In the notebook, run
'latex.add_to_jsmath_avoid_list("xymatrix")' and try again -- you
should get a picture (a part of the diagram arising from a filtered
chain complex)."""

        def _latex_(self):
            """
            LaTeX representation

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: len(latex_examples.diagram()._latex_()) > 1000
                True
            """
            return r"""\xymatrix{
& {} \ar[d] & & \ar[d] & & \ar[d] \\
\ldots \ar[r] & H_{p+q}(K^{p-2}) \ar[r] \ar[d] &
H_{p+q}(K^{p-2}/K^{p-3}) \ar[r] & H_{p+q-1}(K^{p-3}) \ar[r] \ar[d] &
H_{p+q-1}(K^{p-3}/K^{p-4}) \ar[r] & H_{p+q-2}(K^{p-4}) \ar[r] \ar[d] &
\ldots \\
\ldots \ar[r]^{k \quad \quad } & H_{p+q}(K^{p-1}) \ar[r] \ar[d]^{i} &
H_{p+q}(K^{p-1}/K^{p-2}) \ar[r] & H_{p+q-1}(K^{p-2}) \ar[r] \ar[d] &
H_{p+q-1}(K^{p-2}/K^{p-3}) \ar[r] & H_{p+q-2}(K^{p-3}) \ar[r] \ar[d] &
\ldots \\
\ldots \ar[r] & H_{p+q}(K^{p}) \ar[r]^{j} \ar[d] &
H_{p+q}(K^{p}/K^{p-1}) \ar[r]^{k} & H_{p+q-1}(K^{p-1}) \ar[r] \ar[d]^{i} &
H_{p+q-1}(K^{p-1}/K^{p-2}) \ar[r] & H_{p+q-2}(K^{p-2}) \ar[r] \ar[d] &
\ldots \\
\ldots \ar[r] & H_{p+q}(K^{p+1}) \ar[r] \ar[d] &
H_{p+q}(K^{p+1}/K^{p}) \ar[r] & H_{p+q-1}(K^{p}) \ar[r]^{j} \ar[d] &
H_{p+q-1}(K^{p}/K^{p-1}) \ar[r]^{k} & H_{p+q-2}(K^{p-1}) \ar[r] \ar[d]^{i} &
\ldots \\
& {} & {} & {} & {} & {} \\
{} \\
\save "3,1"+DL \PATH ~={**@{-}}
    '+<0pc,-1pc> '+<4pc,0pc> '+<0pc,-4pc> '+<16pc,0pc>
     '+<0pc,-3pc> '+<19pc,0pc>
     '+<0pc,-1pc>
\restore
\save "3,1"+DL \PATH ~={**@{-}}
    '+<0pc,2pc> '+<9pc,0pc> '+<0pc,-3pc> '+<18pc,0pc>
     '+<0pc,-4pc> '+<18pc,0pc>
     '+<0pc,-4pc>
\restore
}"""

latex_examples = LatexExamples()
