"""
LaTeX printing support

In order to support latex formatting, an object should define a
special method ``_latex_(self)`` that returns a string, which will be typeset
in a mathematical mode (the exact mode depends on circumstances).
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


EMBEDDED_MODE = False

COMMON_HEADER = \
r'''\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsfonts}
\usepackage{graphicx}
\usepackage{mathrsfs}
\pagestyle{empty}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
'''

LATEX_HEADER = (
r'''\documentclass{article}
''' + COMMON_HEADER +
r'''\oddsidemargin 0.0in
\evensidemargin 0.0in
\textwidth 6.45in
\topmargin 0.0in
\headheight 0.0in
\headsep 0.0in
\textheight 9.0in
''')

SLIDE_HEADER = (
r'''\documentclass[a0,8pt]{beamer}
''' + COMMON_HEADER +
r'''\textwidth=1.1\textwidth
\textheight=2\textheight
''')

import sys
import shutil, re
import os.path
import random
import subprocess
import types

from sage.misc.temporary_file import tmp_dir, graphics_filename
import sage_eval
from sage.misc.sage_ostools import have_program
from sage.misc.cachefunc import cached_function, cached_method

@cached_function
def have_latex():
    """
    Return ``True`` if this computer has the program ``latex``.

    If this computer doesn't have LaTeX installed, you may obtain it
    from http://ctan.org/.

    EXAMPLES::

        sage: from sage.misc.latex import have_latex
        sage: have_latex() # random
        True
    """
    return have_program('latex')


@cached_function
def have_pdflatex():
    """
    Return ``True`` if this computer has the program ``pdflatex``.

    If this computer doesn't have pdflatex installed, you may obtain it
    from http://ctan.org/.

    EXAMPLES::

        sage: from sage.misc.latex import have_pdflatex
        sage: have_pdflatex() # random
        True
    """
    return have_program('pdflatex')


@cached_function
def have_xelatex():
    """
    Return ``True`` if this computer has the program ``xelatex``.

    If this computer doesn't have xelatex installed, you may obtain it
    from http://ctan.org/.

    EXAMPLES::

        sage: from sage.misc.latex import have_xelatex
        sage: have_xelatex() # random
        True
    """
    return have_program('xelatex')


@cached_function
def have_dvipng():
    """
    Return ``True`` if this computer has the program ``dvipng``.

    If this computer doesn't have dvipng installed, you may obtain it
    from http://sourceforge.net/projects/dvipng/

    EXAMPLES::

        sage: from sage.misc.latex import have_dvipng
        sage: have_dvipng() # random
        True
    """
    return have_program('dvipng')


@cached_function
def have_convert():
    """
    Return ``True`` if this computer has the program ``convert``.

    If this computer doesn't have convert installed, you may obtain it
    (along with the rest of the ImageMagick suite) from
    http://www.imagemagick.org

    EXAMPLES::

        sage: from sage.misc.latex import have_convert
        sage: have_convert() # random
        True
    """
    return have_program('convert')


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


def tuple_function(x, combine_all=False):
    r"""
    Returns the LaTeX code for a tuple ``x``.

    INPUT:

    - ``x`` -- a tuple

    - ``combine_all`` -- boolean (Default: ``False``) If ``combine_all`` is
      ``True``, then it does not return a tuple and instead returns a string
      with all the elements separated by a single space. It does not collapse
      tuples which are inside tuples.

    EXAMPLES::

        sage: from sage.misc.latex import tuple_function
        sage: tuple_function((1,2,3))
        '\\left(1, 2, 3\\right)'

    Check that :trac:`11775` is fixed::

        sage: tuple_function((1,2,3), combine_all=True)
        '1 2 3'
        sage: tuple_function(((1,2),3), combine_all=True)
        '\\left(1, 2\\right) 3'
    """
    if combine_all:
        return " ".join([latex(v) for v in x])
    return "\\left(" + ", ".join([latex(v) for v in x]) + "\\right)"


def bool_function(x):
    r"""
    Returns the LaTeX code for a boolean ``x``.

    INPUT:

    - ``x`` -- boolean

    EXAMPLES::

        sage: from sage.misc.latex import bool_function
        sage: print bool_function(2==3)
        \mathrm{False}
        sage: print bool_function(3==(2+1))
        \mathrm{True}
    """
    return r"\mathrm{%s}" % bool(x)

def builtin_constant_function(x):
    r"""
    Returns the LaTeX code for a builtin constant ``x``.

    INPUT:

    - ``x`` -- builtin constant

    .. SEEALSO:: Python built-in Constants http://docs.python.org/library/constants.html

    EXAMPLES::

        sage: from sage.misc.latex import builtin_constant_function
        sage: builtin_constant_function(True)
        '\\mbox{\\rm True}'
        sage: builtin_constant_function(None)
        '\\mbox{\\rm None}'
        sage: builtin_constant_function(NotImplemented)
        '\\mbox{\\rm NotImplemented}'
        sage: builtin_constant_function(Ellipsis)
        '\\mbox{\\rm Ellipsis}'

    TESTS::

        sage: sage.misc.latex.EMBEDDED_MODE = True
        sage: builtin_constant_function(True)
        '{\\rm True}'
        sage: sage.misc.latex.EMBEDDED_MODE = False
    """
    if EMBEDDED_MODE:
        return "{\\rm %s}"%x
    return "\mbox{\\rm %s}"%x

def None_function(x):
    r"""
    Returns the LaTeX code for ``None``.

    INPUT:

    - ``x`` -- ``None``

    EXAMPLES::

        sage: from sage.misc.latex import None_function
        sage: print None_function(None)
        \mathrm{None}
    """
    assert x is None
    return r"\mathrm{None}"


def str_function(x):
    r"""
    Return a LaTeX representation of the string ``x``.

    The main purpose of this function is to generate LaTeX representation for
    classes that do not provide a customized method.

    If ``x`` contains only digits with, possibly, a single decimal point and/or
    a sign in front, it is considered to be its own representation. Otherwise
    each line of ``x`` is wrapped in a ``\texttt`` command and these lines are
    assembled in a left-justified array. This gives to complicated strings the
    closest look to their "terminal representation".

    .. WARNING::

        Such wrappers **cannot** be used as arguments of LaTeX
        commands or in command definitions. If this causes you any problems,
        they probably can be solved by implementing a suitable ``_latex_``
        method for an appropriate class.

    INPUT:

    - ``x`` -- a string.

    OUTPUT:

    A string

    EXAMPLES::

        sage: from sage.misc.latex import str_function
        sage: str_function('34')
        '34'
        sage: str_function('34.5')
        '34.5'
        sage: str_function('-34.5')
        '-34.5'
        sage: str_function('+34.5')
        '+34.5'
        sage: str_function('hello_world')
        '\\text{\\texttt{hello{\\char`\\_}world}}'
        sage: str_function('-1.00000?') # trac 12178
        '-1.00000?'
    """
    # Check if x is just a number with a possible sign, and/or decimal
    # point, and/or ends with "?"
    if re.match(r'(\+|-)?[0-9]*\.?[0-9]*\??$', x):
        return x
    # Deal with special characters
    char_wrapper = r"{\char`\%s}"
    x = "".join(char_wrapper % c if c in "#$%&\^_{}~" else c for c in x)
    # Avoid grouping spaces into one
    x = x.replace(" ", "{ }")
    # And dashes too, since it causes issues for the command line...
    x = x.replace("-", "{-}")
    # Make it work in math mode, but look like typewriter
    line_wrapper = r"\text{\texttt{%s}}"
    x = "\\\\\n".join(line_wrapper % line for line in x.split("\n"))
    # Preserve line breaks
    if "\n" in x:
        x = "\\begin{array}{l}\n%s\n\\end{array}" % x
    return x


def dict_function(x):
    r"""
    Returns the LaTeX code for a dictionary ``x``.

    INPUT:

    - ``x`` -- a dictionary

    EXAMPLES::

        sage: from sage.misc.latex import dict_function
        sage: x,y,z = var('x,y,z')
        sage: print dict_function({x/2: y^2})
        \left\{\frac{1}{2} \, x : y^{2}\right\}
        sage: d = {(1,2,x^2): [sin(z^2), y/2]}
        sage: latex(d)
        \left\{\left(1, 2, x^{2}\right) :
               \left[\sin\left(z^{2}\right), \frac{1}{2} \, y\right]\right\}
    """
    return "".join([r"\left\{",
                    ", ".join(r"%s : %s" % (latex(key), latex(value))
                              for key, value in x.iteritems()),
                    r"\right\}"])

# One can add to the latex_table in order to install latexing
# functionality for other types.  (Suggested by Robert Kerns of Enthought.)

def float_function(x):
    r"""
    Returns the LaTeX code for a python float ``x``.

    INPUT:

    - ``x`` -- a python float

    EXAMPLES::

        sage: from sage.misc.latex import float_function
        sage: float_function(float(3.14))
        3.14
        sage: float_function(float(1e-10))
        1 \times 10^{-10}
        sage: float_function(float(2e10))
        20000000000.0

    TESTS:

    Check that :trac:`7356` is fixed::

        sage: latex(float(2e-13))
        2 \times 10^{-13}
    """
    from sage.all import RDF
    return latex(RDF(x))


latex_table = {types.NoneType: None_function,
               bool: bool_function,
               dict: dict_function,
               float: float_function,
               int: str,
               list: list_function,
               long: str,
               str: str_function,
               tuple: tuple_function,
               type(None):builtin_constant_function,
               type(NotImplemented):builtin_constant_function,
               type(Ellipsis):builtin_constant_function}


class LatexExpr(str):
    r"""
    A class for LaTeX expressions.

    Normally, objects of this class are created by a :func:`latex` call. It is
    also possible to generate :class:`LatexExpr` directly from a string, which
    must contain valid LaTeX code for typesetting in math mode (without dollar
    signs). In the Sage notebook, use :func:`pretty_print` or the "Typeset"
    checkbox to actually see the typeset LaTeX code; alternatively, from
    either the command-line or the notebook, use the :func:`view` function.

    INPUT:

    - ``str`` -- a string with valid math mode LaTeX code (or something
      which can be converted to such a string).

    OUTPUT:

    - :class:`LatexExpr` wrapping the string representation of the input.

    EXAMPLES::

        sage: latex(x^20 + 1)
        x^{20} + 1
        sage: LatexExpr(r"\frac{x^2 + 1}{x - 2}")
        \frac{x^2 + 1}{x - 2}

    ``LatexExpr`` simply converts to string without doing anything
    extra, it does *not* call :func:`latex`::

        sage: latex(ZZ)
        \Bold{Z}
        sage: LatexExpr(ZZ)
        Integer Ring

    The result of :func:`latex` is of type ``LatexExpr``::

        sage: L = latex(x^20 + 1)
        sage: L
        x^{20} + 1
        sage: type(L)
        <class 'sage.misc.latex.LatexExpr'>

    A ``LatexExpr`` can be converted to a plain string::

        sage: str(latex(x^20 + 1))
        'x^{20} + 1'
    """
    def __add__(self, other):
        r"""
        Add a LatexExpr and another LatexExpr (or a string).

        EXAMPLES::

            sage: o = LatexExpr(r"\Delta\neq") + LatexExpr(r"\frac{x}{2}"); o
            \Delta\neq \frac{x}{2}
            sage: type(o)
            <class 'sage.misc.latex.LatexExpr'>
            sage: o = LatexExpr(r"\Delta\neq") + r"\frac{x}{2}"; o
            \Delta\neq \frac{x}{2}
            sage: type(o)
            <class 'sage.misc.latex.LatexExpr'>

        We add extra space only if it wasn't there yet::

            sage: LatexExpr("foo ") + LatexExpr("bar")
            foo bar
            sage: LatexExpr("foo") + LatexExpr(" bar")
            foo bar
            sage: str(LatexExpr("") + LatexExpr("bar"))
            'bar'
            sage: str(LatexExpr("foo") + LatexExpr(""))
            'foo'
        """
        left = str(self)
        right = str(other)
        # Add a space if self ends with a non-space and other starts
        # with a non-space
        try:
            if left[-1] != ' ' and right[0] != ' ':
                left += ' '
        except IndexError:
            pass
        return LatexExpr(left + right)

    def __radd__(self, other):
        r"""
        Add a string and a LatexExpr.

        EXAMPLES::

            sage: o = "a =" + LatexExpr("b")
            sage: o
            a = b
            sage: type(o)
            <class 'sage.misc.latex.LatexExpr'>
        """
        return LatexExpr(other) + self

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LatexExpr("abc").__repr__()
            'abc'
        """
        return str(self)

    def _latex_(self):
        """
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: latex(LatexExpr("abc")) # indirect doctest
            abc
        """
        return str(self)

def has_latex_attr(x):
    """
    Return ``True`` if ``x`` has a ``_latex_`` attribute, except if ``x``
    is a ``type``, in which case return ``False``.

    EXAMPLES::

        sage: from sage.misc.latex import has_latex_attr
        sage: has_latex_attr(identity_matrix(3))
        True
        sage: has_latex_attr("abc")  # strings have no _latex_ method
        False

    Types inherit the ``_latex_`` method of the class to which they refer,
    but calling it is broken::

        sage: T = type(identity_matrix(3)); T
        <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
        sage: hasattr(T, '_latex_')
        True
        sage: T._latex_()
        Traceback (most recent call last):
        ...
        TypeError: descriptor '_latex_' of 'sage.matrix.matrix0.Matrix' object needs an argument
        sage: has_latex_attr(T)
        False
    """
    return hasattr(x, '_latex_') and not isinstance(x, type)

from sage.structure.sage_object import SageObject

class _Latex_prefs_object(SageObject):
    """
    An object that holds LaTeX global preferences.
    """
    def __init__(self, bb=False, delimiters=["(", ")"],
                 matrix_column_alignment="r"):
        """
        Define an object that holds LaTeX global preferences.

        EXAMPLES::

            sage: from sage.misc.latex import _Latex_prefs_object
            sage: latex_prefs = _Latex_prefs_object()
            sage: TestSuite(latex_prefs).run(skip ="_test_pickling")
        """
        self._option = {}
        self._option["blackboard_bold"] = bb
        self._option["matrix_delimiters"] = list(delimiters)
        self._option["vector_delimiters"] = list(delimiters)
        self._option["matrix_column_alignment"] = matrix_column_alignment
        self._option["macros"] = ""
        self._option["preamble"] = ""
        self._option["engine"] = "pdflatex"
        self._option["engine_name"] = "LaTeX"
        self._option["mathjax_avoid"] = []

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
        sage: print latex_extra_preamble()
        ...
        <BLANKLINE>
        \newcommand{\ZZ}{\Bold{Z}}
        \newcommand{\NN}{\Bold{N}}
        \newcommand{\RR}{\Bold{R}}
        \newcommand{\CC}{\Bold{C}}
        \newcommand{\QQ}{\Bold{Q}}
        \newcommand{\QQbar}{\overline{\QQ}}
        \newcommand{\GF}[1]{\Bold{F}_{#1}}
        \newcommand{\Zp}[1]{\ZZ_{#1}}
        \newcommand{\Qp}[1]{\QQ_{#1}}
        \newcommand{\Zmod}[1]{\ZZ/#1\ZZ}
        \newcommand{\CDF}{\Bold{C}}
        \newcommand{\CIF}{\Bold{C}}
        \newcommand{\CLF}{\Bold{C}}
        \newcommand{\RDF}{\Bold{R}}
        \newcommand{\RIF}{\Bold{I} \Bold{R}}
        \newcommand{\RLF}{\Bold{R}}
        \newcommand{\CFF}{\Bold{CFF}}
        \newcommand{\Bold}[1]{\mathbf{#1}}
        <BLANKLINE>
    """
    from sage.misc.latex_macros import sage_latex_macros
    return "\n".join([_Latex_prefs._option['preamble'],
                     "\n".join(sage_latex_macros()),
                     _Latex_prefs._option['macros']])

def _run_latex_(filename, debug=False, density=150, engine=None, png=False, do_in_background=False):
    """
    This runs LaTeX on the TeX file "filename.tex".  It produces files
    "filename.dvi" (or "filename.pdf"` if engine is either ``pdflatex``
    or ``xelatex``) and if ``png`` is ``True``, "filename.png".  If ``png``
    is ``True`` and dvipng can't convert the dvi file to png (because of
    postscript specials or other issues), then dvips is called, and the
    PS file is converted to a png file.

    INPUT:

    -  ``filename`` -- string: file to process, including full path

    -  ``debug`` -- bool (optional, default ``False``): whether to print
       verbose debugging output

    -  ``density`` -- integer (optional, default 150): how big output
       image is.

    -  ``engine`` -- string: latex engine to use.

    -  ``png`` -- bool (optional, default ``False``): whether to produce a
       png file.

    -  ``do_in_background`` -- bool (optional, default ``False``).  Unused,
       kept for backwards compatibility.

    OUTPUT:

    A string which could be a string starting with 'Error' (if
    there was a problem), or it could be 'pdf' or 'dvi'.  If
    engine is latex or ``None``, then a dvi file is created, but if there
    appear to be problems with it (because of PS special commands, for
    example), then a pdf file is created instead.  The function
    returns 'dvi' or 'pdf' to indicate which type of file is created.
    (Detecting problems requires that dvipng be installed; if it is
    not, then the dvi file is not checked for problems and 'dvi' is
    returned.)  If engine is pdflatex or xelatex and there are no errors, then
    'pdf' is returned.

    .. WARNING::

       If ``png`` is ``True``, then when using latex (the default), you
       must have 'dvipng' (or 'dvips' and 'convert') installed on your
       operating system, or this command won't work.  When using
       pdflatex or xelatex, you must have 'convert' installed.

    EXAMPLES::

        sage: from sage.misc.latex import _run_latex_, _latex_file_
        sage: file = os.path.join(SAGE_TMP, "temp.tex")
        sage: O = open(file, 'w')
        sage: O.write(_latex_file_([ZZ[x], RR])); O.close()
        sage: _run_latex_(file) # random - depends on whether latex is installed
        'dvi'
    """
    if engine is None:
        engine = _Latex_prefs._option["engine"]

    if not engine or engine == "latex":
        if not have_latex():
            print("Error: LaTeX does not seem to be installed.  Download it from")
            print("ctan.org and try again.")
            return "Error"
        command = "latex"
        # 'suffix' is used in the 'convert' command list
        suffix = "ps"
        return_suffix = "dvi"
    elif engine == "pdflatex":
        if not have_pdflatex():
            print("Error: PDFLaTeX does not seem to be installed.  Download it from")
            print("ctan.org and try again.")
            return "Error"
        command = "pdflatex"
        suffix = "pdf"
        return_suffix = "pdf"
    elif engine == "xelatex":
        if not have_xelatex():
            print("Error: XeLaTeX does not seem to be installed.  Download it from")
            print("ctan.org and try again.")
            return "Error"
        command = "xelatex"
        suffix = "pdf"
        return_suffix = "pdf"
    else:
        raise ValueError("Unsupported LaTeX engine.")

    # if png output + latex, check to see if dvipng or convert is installed.
    if png:
        if (not engine or engine == "latex") and not (have_dvipng() or have_convert()):
            print()
            print("Error: neither dvipng nor convert (from the ImageMagick suite)")
            print("appear to be installed. Displaying LaTeX, PDFLaTeX output")
            print("requires at least one of these programs, so please install")
            print("and try again.")
            print()
            print("Go to http://sourceforge.net/projects/dvipng/ and")
            print("http://www.imagemagick.org to download these programs.")
            return "Error"
    # if png output + pdflatex, check to see if convert is installed.
        elif engine == "pdflatex" and not have_convert():
            print()
            print("Error: convert (from the ImageMagick suite) does not")
            print("appear to be installed. Displaying PDFLaTeX output")
            print("requires this program, so please install and try again.")
            print()
            print("Go to http://www.imagemagick.org to download it.")
            return "Error"
        elif engine == "xelatex" and not have_convert():
            print()
            print("Error: convert (from the ImageMagick suite) does not")
            print("appear to be installed. Displaying XeLaTeX output")
            print("requires this program, so please install and try again.")
            print()
            print("Go to http://www.imagemagick.org to download it.")
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
        raise ValueError("filename must contain no spaces")
    if not debug:
        redirect = subprocess.PIPE
    else:
        redirect = None
    # if do_in_background:
    #     background = ' &'
    # else:
    #     background = ''

    # Define the commands to be used:
    lt = ['sage-native-execute', command, r'\nonstopmode', r'\input{' + filename + '.tex}']
    # dvipng is run with the 'picky' option: this means that if
    # there are warnings, no png file is created.
    dvipng = ['sage-native-execute', 'dvipng', '--picky', '-q', '-T', 'tight',
              '-D', str(density), filename + '.dvi', '-o', filename + '.png']

    dvips = ['sage-native-execute', 'dvips', filename + '.dvi']

    ps2pdf = ['sage-native-execute', 'ps2pdf', filename + '.ps']

    # We seem to need a larger size when using convert compared to
    # when using dvipng:
    density = int(1.4 * density / 1.3)
    convert = ['sage-native-execute', 'convert', '-density',
               '{0}x{0}'.format(density), '-trim', filename + '.' + suffix,
               filename + '.png']

    e = False # it is possible to get through the following commands
              # without running a program, so in that case we force error

    # our standard way of calling programs here; change this if we want
    # finer-grained analysis of the return code. Think of the output as
    # a boolean: "the command exited normally"
    subpcall = lambda x: not subprocess.call(x, stdout=redirect,
                                             stderr=redirect, cwd=base)
    if engine == "pdflatex" or engine == "xelatex":
        if debug:
            print(lt)
            if png:
                print(convert)
        e = subpcall(lt)
        if png:
            e = e and subpcall(convert)
    else:  # latex
        if (png or check_validity):
            if have_dvipng():
                if debug:
                    print(lt)
                    print(dvipng)
                e = subpcall(lt) and subpcall(dvipng)
                dvipng_error = not os.path.exists(os.path.join(base, filename + '.png'))
                # If there is no png file, then either the latex
                # process failed or dvipng failed.  Assume that dvipng
                # failed, and try running dvips and convert.  (If the
                # latex process failed, then dvips and convert will
                # fail also, so we'll still catch the error.)
                if dvipng_error:
                    if png:
                        if have_convert():
                            if debug:
                                print("'dvipng' failed; trying 'convert' instead...")
                                print(dvips)
                                print(convert)
                            e = subpcall(dvips) and subpcall(convert)
                        else:
                            print("Error: 'dvipng' failed and 'convert' is not installed.")
                            return "Error: dvipng failed."
                    else:  # not png, i.e., check_validity
                        return_suffix = "pdf"
                        if debug:
                            print("bad dvi file; running dvips and ps2pdf instead...")
                            print(dvips)
                            print(ps2pdf)
                        e = subpcall(dvips) and subpcall(ps2pdf)
                        if not e:  # error running dvips and/or ps2pdf
                            pdflt = lt[:]
                            pdflt[1] = 'pdflatex'
                            if debug:
                                print("error running dvips and ps2pdf; trying pdflatex instead...")
                                print(pdflt)
                            e = subpcall(pdflt)
            else:  # don't have dvipng, so must have convert.  run latex, dvips, convert.
                if debug:
                    print(lt)
                    print(dvips)
                    print(convert)
                e = subpcall(lt) and subpcall(dvips) and subpcall(convert)
    if not e:
        print("An error occurred.")
        try:
            print(open(base + '/' + filename + '.log').read())
        except IOError:
            pass
        return "Error latexing slide."
    return return_suffix

class LatexCall:
    """
    Typeset Sage objects via a ``__call__`` method to this class,
    typically by calling those objects' ``_latex_`` methods.  The
    class :class:`Latex` inherits from this. This class is used in
    :mod:`~sage.misc.latex_macros`, while functions from
    :mod:`~sage.misc.latex_macros` are used in :class:`Latex`, so this
    is here primarily to avoid circular imports.

    EXAMPLES::

        sage: from sage.misc.latex import LatexCall
        sage: LatexCall()(ZZ)
        \Bold{Z}
        sage: LatexCall().__call__(ZZ)
        \Bold{Z}

    This returns an instance of the class :class:`LatexExpr`::

        sage: type(LatexCall()(ZZ))
        <class 'sage.misc.latex.LatexExpr'>
    """
    def __call__(self, x, combine_all=False):
        r"""
        Return a :class:`LatexExpr` built out of the argument ``x``.

        INPUT:

        - ``x`` -- a Sage object

        - ``combine_all`` -- boolean (Default: ``False``) If ``combine_all``
          is ``True`` and the input is a tuple, then it does not return a
          tuple and instead returns a string with all the elements separated by
          a single space.

        OUTPUT:

        A :class:`LatexExpr` built from ``x``

        EXAMPLES::

            sage: latex(Integer(3))  # indirect doctest
            3
            sage: latex(1==0)
            \mathrm{False}
            sage: print latex([x,2])
            \left[x, 2\right]

        Check that :trac:`11775` is fixed::

            sage: latex((x,2), combine_all=True)
            x 2
        """
        if has_latex_attr(x):
            return LatexExpr(x._latex_())
        try:
            f = latex_table[type(x)]
            if isinstance(x, tuple):
                return LatexExpr(f(x, combine_all=combine_all))
            return LatexExpr(f(x))
        except KeyError:
            return LatexExpr(str_function(str(x)))


class Latex(LatexCall):
    r"""nodetex
    Enter, e.g.,

    ::

        %latex
        The equation $y^2 = x^3 + x$ defines an elliptic curve.
        We have $2006 = \sage{factor(2006)}$.

    in an input cell in the notebook to get a typeset version. Use
    ``%latex_debug`` to get debugging output.

    Use ``latex(...)`` to typeset a Sage object.  Use :class:`LatexExpr`
    to typeset LaTeX code that you create by hand.

    Use ``%slide`` instead to typeset slides.

    .. WARNING::

       You must have dvipng (or dvips and convert) installed
       on your operating system, or this command won't work.

    EXAMPLES::

        sage: latex(x^20 + 1)
        x^{20} + 1
        sage: latex(FiniteField(25,'a'))
        \Bold{F}_{5^{2}}
        sage: latex("hello")
        \text{\texttt{hello}}
        sage: LatexExpr(r"\frac{x^2 - 1}{x + 1} = x - 1")
        \frac{x^2 - 1}{x + 1} = x - 1

    LaTeX expressions can be added; note that a space is automatically
    inserted::

        sage: LatexExpr(r"y \neq") + latex(x^20 + 1)
        y \neq x^{20} + 1
    """
    def __init__(self, debug=False, slide=False, density=150, pdflatex=None, engine=None):
        """
        Initialize the latex builder.

        EXAMPLES::

            sage: from sage.misc.latex import Latex
            sage: l = Latex()
            sage: TestSuite(l).run(skip ="_test_pickling")
        """
        self.__debug = debug
        self.__slide = slide
        self.__pdflatex = pdflatex
        self.__engine = engine
        self.__density = density

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
        Replace instances of ``'\sage{x}'`` in ``s`` with the LaTeX version of
        ``x`` in the running session.

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
                k = str(latex(sage_eval.sage_eval(var, locals)))
            except Exception as msg:
                print(msg)
                k = '\\mbox{\\rm [%s undefined]}'%var
            s = s[:i] + k + t[j+1:]

    def eval(self, x, globals, strip=False, filename=None, debug=None,
             density=None, pdflatex=None, engine=None, locals={}):
        r"""
        Compiles the formatted tex given by ``x`` as a png and writes the
        output file to the directory given by ``filename``.

        INPUT:

        -  ``globals`` -- a globals dictionary

        -  ``x`` -- string to evaluate.

        -  ``strip`` -- ignored

        -  ``filename`` -- output filename

        -  ``debug`` -- whether to print verbose debugging
           output

        -  ``density`` -- how big output image is.

        -  ``pdflatex`` -- whether to use pdflatex. This is deprecated. Use
           ``engine`` option instead.

        -  ``engine`` -- latex engine to use. Currently latex, pdflatex, and
           xelatex are supported.

        -  ``locals`` - extra local variables used when
           evaluating Sage code in ``x``.

        .. WARNING::

           When using latex (the default), you must have 'dvipng' (or
           'dvips' and 'convert') installed on your operating system,
           or this command won't work.  When using pdflatex or xelatex, you
           must have 'convert' installed.

        OUTPUT:

        If it compiled successfully, this returns an empty string ``''``,
        otherwise it returns ``None``.

        EXAMPLES::

            # This would generate a file named "test.png"
            sage: latex.eval("\\ZZ[x]", locals(), filename="test") # not tested
            ''
            # This would generate a file named "/path/to/test.png"
            sage: latex.eval("\\ZZ[x]", locals(), filename="/path/to/test") # not tested
            ''
            sage: latex.eval("\ThisIsAnInvalidCommand", {}) # optional -- ImageMagick
            An error
            ...
            No pages of output.
            <BLANKLINE>
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
            raise ValueError("filename must contain no spaces")
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

        from sagenb.misc.misc import encoded_str
        O.write(encoded_str(x))
        if self.__slide:
            O.write('\n\n\\end{document}')
        else:
            O.write('\n\n\\end{document}\n')

        O.close()
        if engine is None:
            if self.__engine is None:
                engine = _Latex_prefs._option["engine"]
            else:
                engine = self.__engine
        e = _run_latex_(os.path.join(base, filename + ".tex"), debug=debug,
                               density=density, engine=engine, png=True)
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

        - ``t`` -- boolean or ``None``

        OUTPUT:

        If ``t`` is ``None``, return the current setting (``True`` or
        ``False``).

        If ``t`` is ``True``, use blackboard bold (``\mathbb``); otherwise use
        boldface (``\mathbf``).

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
        from latex_macros import sage_configurable_latex_macros
        global sage_configurable_latex_macros
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
            sage_configurable_latex_macros.remove(old_macro)
            sage_configurable_latex_macros.append(macro)

    def matrix_delimiters(self, left=None, right=None):
        r"""nodetex
        Change the left and right delimiters for the LaTeX representation
        of matrices

        INPUT:

        - ``left``, ``right`` - strings or ``None``

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

        .. NOTE::

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

        - ``left``, ``right`` -- strings or ``None``

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

        .. NOTE::

           Putting aside aesthetics, you may combine these in any way
           imaginable; for example, you could set ``left`` to be a
           right-hand bracket ']' and ``right`` to be a right-hand
           brace '\\}', and it will be typeset correctly.

        EXAMPLES::

            sage: a = vector(QQ, [1,2,3])
            sage: latex(a)
            \left(1,\,2,\,3\right)
            sage: latex.vector_delimiters("[", "]")
            sage: latex(a)
            \left[1,\,2,\,3\right]
            sage: latex.vector_delimiters(right="\\}")
            sage: latex(a)
            \left[1,\,2,\,3\right\}
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

    def matrix_column_alignment(self, align=None):
        r"""nodetex
        Changes the column-alignment of the LaTeX representation of
        matrices.

        INPUT:

        - ``align`` - a string (``'r'`` for right, ``'c'`` for center,
          ``'l'`` for left) or ``None``.

        OUTPUT:

        If ``align`` is ``None``, then returns the current
        alignment-string. Otherwise, set this alignment.

        The input ``align`` can be any string which the LaTeX
        ``array``-environment understands as a parameter for
        aligning a column.

        EXAMPLES::

            sage: a = matrix(1, 1, [42])
            sage: latex(a)
            \left(\begin{array}{r}
            42
            \end{array}\right)
            sage: latex.matrix_column_alignment('c')
            sage: latex(a)
            \left(\begin{array}{c}
            42
            \end{array}\right)
            sage: latex.matrix_column_alignment('l')
            sage: latex(a)
            \left(\begin{array}{l}
            42
            \end{array}\right)

        Restore defaults::

            sage: latex.matrix_column_alignment('r')
        """
        if align is None:
            return _Latex_prefs._option['matrix_column_alignment']
        else:
            _Latex_prefs._option['matrix_column_alignment'] = align

    @cached_method
    def has_file(self, file_name):
        """
        INPUT:

        - ``file_name`` -- a string

        Tests whether the local LaTeX installation includes ``file_name``.

        EXAMPLES::

            sage: latex.has_file("article.cls")      # optional - latex
            True
            sage: latex.has_file("some_inexistent_file.sty")
            False
        """
        assert isinstance(file_name, str)
        from subprocess import call, PIPE
        try:
            retcode = call("kpsewhich %s"%file_name, shell=True, stdout=PIPE, stderr=PIPE)
            return (retcode == 0)
        except OSError:
            return False

    @cached_method
    def check_file(self, file_name, more_info = ""):
        """
        INPUT:

        - ``file_name`` -- a string

        - ``more_info`` -- a string (default: "")

        Emit a warning if the local LaTeX installation does not
        include ``file_name``. The string ``more_info`` is appended
        to the warning message. The warning is only emitted the first
        time this method is called.

        EXAMPLES::

            sage: latex.check_file("article.cls")       # optional - latex
            sage: latex.check_file("some_inexistent_file.sty")
            Warning: `some_inexistent_file.sty` is not part of this computer's TeX installation.
            sage: latex.check_file("some_inexistent_file.sty")
            sage: latex.check_file("some_inexistent_file.sty", "This file is required for blah. It can be downloaded from: http://blah.org/")
            Warning: `some_inexistent_file.sty` is not part of this computer's TeX installation.
            This file is required for blah. It can be downloaded from: http://blah.org/

        This test checks that the bug in :trac:`9091` is fixed::

            sage: latex.check_file("article.cls", "The article class is really critical.")    # optional - latex
        """
        assert isinstance(file_name, str)
        if not self.has_file(file_name):
            print("""
Warning: `{}` is not part of this computer's TeX installation.""".format(file_name))
            if more_info:
                print(more_info)


    def extra_macros(self, macros=None):
        r"""nodetex
        String containing extra LaTeX macros to use with %latex,
        %html, and %mathjax.

        INPUT:

        - ``macros`` -- string (default: ``None``)

        If ``macros`` is ``None``, return the current string.  Otherwise,
        set it to ``macros``.  If you want to *append* to the string
        of macros instead of replacing it, use
        :meth:`latex.add_macro <Latex.add_macro>`.

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
        %latex, %html, and %mathjax.

        INPUT:

        - ``macro`` -- string

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
        Anything in this string won't be processed by %mathjax.

        INPUT:

        - ``s`` -- string or ``None``

        If ``s`` is ``None``, return the current preamble.  Otherwise, set
        it to ``s``.  If you want to *append* to the current extra
        preamble instead of replacing it, use
        :meth:`latex.add_to_preamble <Latex.add_to_preamble>`.

        You will almost certainly need to use this when using the
        XeLaTeX engine; see below or the documentation for
        :func:`engine` for a suggested preamble.

        EXAMPLES::

            sage: latex.extra_preamble("\\DeclareMathOperator{\\Ext}{Ext}")
            sage: latex.extra_preamble()
            '\\DeclareMathOperator{\\Ext}{Ext}'
            sage: latex.extra_preamble("\\"+r"usepackage{fontspec,xunicode,xltxtra}\setmainfont[Mapping=tex-text]{UnBatang}\setmonofont[Mapping=tex-text]{UnDotum}")
            sage: latex.extra_preamble()
            '\\usepackage{fontspec,xunicode,xltxtra}\\setmainfont[Mapping=tex-text]{UnBatang}\\setmonofont[Mapping=tex-text]{UnDotum}'
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
        Append to the string ``s`` of extra LaTeX macros, for use with
        %latex.  Anything in this string won't be processed by
        %mathjax.

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

    def add_package_to_preamble_if_available(self, package_name):
        r"""
        Adds a ``\usepackage{package_name}`` instruction to the latex
        preamble if not yet present there, and if ``package_name.sty``
        is available in the LaTeX installation.

        INPUT:

        - ``package_name`` -- a string

        .. SEEALSO::

            - :meth:`add_to_preamble`
            - :meth:`has_file`.

        TESTS::

            sage: latex.add_package_to_preamble_if_available("xypic")
            sage: latex.add_package_to_preamble_if_available("nonexistent_package")
            sage: latex.extra_preamble()       # optional - latex
            '\\usepackage{xypic}\n'
            sage: latex.extra_preamble('')
        """
        assert isinstance(package_name, str)
        if self.has_file(package_name+".sty"):
            self.add_to_preamble("\\usepackage{%s}\n"%package_name)

    def mathjax_avoid_list(self, L=None):
        r"""nodetex
        List of strings which signal that MathJax should not
        be used when 'view'ing.

        INPUT:

        - ``L`` -- A list or ``None``

        If ``L`` is ``None``, then return the current list.
        Otherwise, set it to ``L``.  If you want to *append* to the
        current list instead of replacing it, use
        :meth:`latex.add_to_mathjax_avoid_list <Latex.add_to_mathjax_avoid_list>`.

        EXAMPLES::

            sage: latex.mathjax_avoid_list(["\\mathsf", "pspicture"])
            sage: latex.mathjax_avoid_list()  # display current setting
            ['\\mathsf', 'pspicture']
            sage: latex.mathjax_avoid_list([])  # reset to default
            sage: latex.mathjax_avoid_list()
            []
        """
        if L is None:
            return _Latex_prefs._option['mathjax_avoid']
        else:
            _Latex_prefs._option['mathjax_avoid'] = L

    # Couldn't use deprecated_function_alias for this because of circular imports.
    def jsmath_avoid_list(self, L=None):
        """
        Deprecated. Use :meth:`mathjax_avoid_list` instead.

        EXAMPLES::

            sage: latex.jsmath_avoid_list()
            doctest:...: DeprecationWarning: Use mathjax_avoid_list instead.
            See http://trac.sagemath.org/13508 for details.
            []
        """
        from superseded import deprecation
        deprecation(13508, 'Use mathjax_avoid_list instead.')
        if L is None:
            return _Latex_prefs._option['mathjax_avoid']
        else:
            _Latex_prefs._option['mathjax_avoid'] = L

    def add_to_mathjax_avoid_list(self, s):
        r"""nodetex
        Add to the list of strings which signal that MathJax should not
        be used when 'view'ing.

        INPUT:

        - ``s`` -- string; add ``s`` to the list of 'MathJax avoid' strings

        If you want to replace the current list instead of adding to
        it, use :meth:`latex.mathjax_avoid_list <Latex.mathjax_avoid_list>`.

        EXAMPLES::

            sage: latex.add_to_mathjax_avoid_list("\\mathsf")
            sage: latex.mathjax_avoid_list()  # display current setting
            ['\\mathsf']
            sage: latex.add_to_mathjax_avoid_list("tkz-graph")
            sage: latex.mathjax_avoid_list()  # display current setting
            ['\\mathsf', 'tkz-graph']
            sage: latex.mathjax_avoid_list([])  # reset to default
            sage: latex.mathjax_avoid_list()
            []
        """
        current = latex.mathjax_avoid_list()
        if s not in current:
            _Latex_prefs._option['mathjax_avoid'].append(s)

    # Couldn't use deprecated_function_alias for this because of circular imports.
    def add_to_jsmath_avoid_list(self, s):
        """
        Deprecated. Use :meth:`add_to_mathjax_avoid_list` instead.

        EXAMPLES::

            sage: latex.add_to_jsmath_avoid_list('\\text')
            doctest:...: DeprecationWarning: Use add_to_mathjax_avoid_list instead.
            See http://trac.sagemath.org/13508 for details.
            sage: latex.mathjax_avoid_list([])  # reset list to default
        """
        from superseded import deprecation
        deprecation(13508, 'Use add_to_mathjax_avoid_list instead.')
        self.add_to_mathjax_avoid_list(s)

    def engine(self, e = None):
        r"""
        Set Sage to use ``e`` as latex engine when typesetting with
        :func:`view`, in ``%latex`` cells, etc.

        INPUT:

        - ``e`` -- 'latex', 'pdflatex', 'xelatex' or ``None``

        If  ``e`` is ``None``, return the current engine.

        If using the XeLaTeX engine, it will almost always be necessary
        to set the proper preamble with :func:`extra_preamble` or
        :func:`add_to_preamble`. For example::

            latex.extra_preamble(r'''\usepackage{fontspec,xunicode,xltxtra}
            \setmainfont[Mapping=tex-text]{some font here}
            \setmonofont[Mapping=tex-text]{another font here}''')

        EXAMPLES::

            sage: latex.engine()
            'pdflatex'
            sage: latex.engine("latex")
            sage: latex.engine()
            'latex'
            sage: latex.engine("xelatex")
            sage: latex.engine()
            'xelatex'
        """
        if e is None:
            return _Latex_prefs._option["engine"]

        if e == "latex":
            _Latex_prefs._option["engine"] = "latex"
            _Latex_prefs._option["engine_name"] = "LaTeX"
        elif e == "pdflatex":
            _Latex_prefs._option["engine"] = "pdflatex"
            _Latex_prefs._option["engine_name"] = "PDFLaTeX"
        elif e == "xelatex":
            _Latex_prefs._option["engine"] = e
            _Latex_prefs._option["engine_name"] = "XeLaTeX"
        else:
            raise ValueError("%s is not a supported LaTeX engine. Use latex, pdflatex, or xelatex" % e)

# Note: latex used to be a separate function, which by default was
# only loaded in command-line mode: in the notebook, all_notebook.py
# defined (and still defines) latex by 'latex = Latex(density=130)'.
# Meanwhile, the __call__ method for Latex used to call the latex
# function.  This has been changed around so that the contents of the
# old latex function are now in Latex.__call__; thus the following
# assignment.

latex = Latex()
# Ensure that latex appear in the sphinx doc as a function
# so that the link :func:`latex` is correctly set up.
latex.__doc__  = Latex.__call__.__doc__
#########################################

def _latex_file_(objects, title='SAGE', debug=False, \
                 sep='', tiny=False, math_left='\\[',
                 math_right='\\]',
                 extra_preamble=''):
    r"""nodetex
    Produce a string to be used as a LaTeX file, containing a
    representation of each object in objects.

    INPUT:

    -  ``objects`` -- list (or object)

    -  ``title`` -- string (default: 'Sage'): title for the document

    -  ``math_left`` -- string (default: '\\['), left delimiter for math mode

    -  ``math_right`` -- string (default: '\\]'), right delimiter for math mode

    -  ``debug`` -- bool (default: False): print verbose output

    -  ``sep`` -- string (default: ''): separator between math objects

    -  ``tiny`` -- bool (default: False): use 'tiny' font.

    -  ``extra_preamble`` -- string (default: ''): extra LaTeX commands,
       inserted before "\\begin{document}"

    This creates a string intended to be a LaTeX file containing the
    LaTeX representations of objects. It contains the following:

    - a header (with documentclass and usepackage commands)

    - ``extra_preamble``

    - the title (centered)

    - a size specification if ``tiny`` is ``True``

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

    TESTS:

    This makes sure that latex is called only once on an object::

        sage: class blah():
        ...       def _latex_(x):
        ...           print "coucou"
        ...           return "x"
        sage: latex(blah())
        coucou
        x
        sage: s = sage.misc.latex._latex_file_(blah())
        coucou
    """
    process = True
    if has_latex_attr(objects):
        objects = [objects]

    if not isinstance(objects, list):
        objects = [objects]

    if tiny:
        size='\\tiny\n'
    else:
        size=''

    s = '%s\n\\begin{document}\n\\begin{center}{\\Large\\bf %s}\\end{center}\n%s'%(
        extra_preamble, title, size)

    #s += "(If something is missing it may be on the next page or there may be errors in the latex.  Use view with {\\tt debug=True}.)\\vfill"
    s += '\\vspace{40mm}'
    if process:
        for i in range(len(objects)):
            x = objects[i]
            L = latex(x)
            if not '\\begin{verbatim}' in L:
                s += '%s%s%s'%(math_left, L, math_right)
            else:
                s += '%s'%L
            if i < len(objects)-1:
                s += '\n\n%s\n\n'%sep
    else:
        s += "\n\n".join([str(x) for x in objects])

    # latex_extra_preamble() is called here and not before because some objects
    # may require additional packages to be displayed in LaTeX. Hence, the call
    # to latex(x) in the previous loop may change the result of
    # latex_extra_preamble()
    MACROS = latex_extra_preamble()
    s = LATEX_HEADER + '\n' + MACROS + s + '\n\\end{document}'

    if debug:
        print(s)

    return s

class MathJaxExpr:
    """
    An arbitrary MathJax expression that can be nicely concatenated.

    EXAMPLES::

        sage: from sage.misc.latex import MathJaxExpr
        sage: MathJaxExpr("a^{2}") + MathJaxExpr("x^{-1}")
        a^{2}x^{-1}
    """
    def __init__(self, y):
        """
        Initialize a MathJax expression.

        INPUT:

        - ``y`` - a string

        Note that no error checking is done on the type of ``y``.

        EXAMPLES::

            sage: from sage.misc.latex import MathJaxExpr
            sage: jax = MathJaxExpr(3); jax  # indirect doctest
            3
            sage: TestSuite(jax).run(skip ="_test_pickling")
        """
        self.__y = y

    def __repr__(self):
        """
        Print representation.

        EXAMPLES::

            sage: from sage.misc.latex import MathJaxExpr
            sage: jax = MathJaxExpr('3')
            sage: jax.__repr__()
            '3'
        """
        return str(self.__y)

    def __add__(self, y):
        """
        'Add' MathJaxExpr ``self`` to ``y``.  This concatenates them
        (assuming that they're strings).

        EXAMPLES::

            sage: from sage.misc.latex import MathJaxExpr
            sage: j3 = MathJaxExpr('3')
            sage: jx = MathJaxExpr('x')
            sage: j3 + jx
            3x
        """
        return MathJaxExpr(self.__y + y)

    def __radd__(self, y):
        """
        'Add' MathJaxExpr ``y`` to ``self``.  This concatenates them
        (assuming that they're strings).

        EXAMPLES::

            sage: from sage.misc.latex import MathJaxExpr
            sage: j3 = MathJaxExpr('3')
            sage: jx = MathJaxExpr('x')
            sage: j3.__radd__(jx)
            x3
        """
        return MathJaxExpr(y + self.__y)

class MathJax:
    r"""
    Render LaTeX input using MathJax.  This returns a :class:`MathJaxExpr`.

    EXAMPLES::

        sage: from sage.misc.latex import MathJax
        sage: MathJax()(3)
        <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}3</script></html>
        sage: MathJax()(ZZ)
        <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}</script></html>
    """

    def __call__(self, x, combine_all=False):
        r"""
        Render LaTeX input using MathJax.  This returns a :class:`MathJaxExpr`.

        INPUT:

        - ``x`` - a Sage object

        - ``combine_all`` - boolean (Default: ``False``): If ``combine_all`` is
          ``True`` and the input is a tuple, then it does not return a tuple
          and instead returns a string with all the elements separated by
          a single space.

        OUTPUT:

        A :calss:`MathJaxExpr`

        EXAMPLES::

            sage: from sage.misc.latex import MathJax
            sage: MathJax()(3)
            <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}3</script></html>
            sage: str(MathJax().eval(ZZ[x], mode='display')) == str(MathJax()(ZZ[x]))
            True
        """
        return self.eval(x, combine_all=combine_all)

    def eval(self, x, globals=None, locals=None, mode='display',
            combine_all=False):
        r"""
        Render LaTeX input using MathJax.  This returns a :class:`MathJaxExpr`.

        INPUT:

        - ``x`` - a Sage object

        -  ``globals`` - a globals dictionary

        -  ``locals`` - extra local variables used when
           evaluating Sage code in ``x``.

        -  ``mode`` - string (optional, default ``'display'``): ``'display'``
           for displaymath or ``'inline'`` for inline math

        - ``combine_all`` - boolean (Default: ``False``): If ``combine_all`` is
          ``True`` and the input is a tuple, then it does not return a tuple
          and instead returns a string with all the elements separated by
          a single space.

        OUTPUT:

        A :class:`MathJaxExpr`

        EXAMPLES::

            sage: from sage.misc.latex import MathJax
            sage: MathJax().eval(3, mode='display')
            <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}3</script></html>
            sage: MathJax().eval(3, mode='inline')
            <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}3</script></html>
            sage: MathJax().eval(type(3), mode='inline')
            <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\verb|<type|\phantom{\verb!x!}\verb|'sage.rings.integer.Integer'>|</script></html>
        """
        # Get a regular LaTeX representation of x
        x = latex(x, combine_all=combine_all)

        # The following block, hopefully, can be removed in some future MathJax.
        prefix = r"\text{\texttt{"
        parts = x.split(prefix)
        for i, part in enumerate(parts):
            if i == 0:
                continue    # Nothing to do with the head part
            n = 1
            for closing, c in enumerate(part):
                if c == "{" and part[closing - 1] != "\\":
                    n += 1
                if c == "}" and part[closing - 1] != "\\":
                    n -= 1
                if n == -1:
                    break
            # part should end in "}}", so omit the last two characters
            # from y
            y = part[:closing-1]
            for delimiter in """|"'`#%&,.:;?!@_~^+-/\=<>()[]{}0123456789E""":
                if delimiter not in y:
                    break
            if delimiter == "E":
                # y is too complicated
                delimiter = "|"
                y = "(complicated string)"
            wrapper = r"\verb" + delimiter + "%s" + delimiter
            spacer = r"\phantom{\verb!%s!}"
            y = y.replace("{ }", " ").replace("{-}", "-")
            for c in r"#$%&\^_{}~":
                char_wrapper = r"{\char`\%s}" % c
                y = y.replace(char_wrapper, c)
            subparts = []
            nspaces = 0
            for subpart in y.split(" "):
                if subpart == "":
                    nspaces += 1
                    continue
                if nspaces > 0:
                    subparts.append(spacer % ("x" * nspaces))
                nspaces = 1
                subparts.append(wrapper % subpart)
            # There is a bug with omitting empty lines in arrays
            if not y:
                subparts.append(spacer % "x")
            subparts.append(part[closing + 1:])
            parts[i] = "".join(subparts)
        x = "".join(parts)

        # In MathJax:
        #   inline math: <script type="math/tex">...</script>
        #   displaymath: <script type="math/tex; mode=display">...</script>
        from sage.misc.latex_macros import sage_configurable_latex_macros
        if mode == 'display':
            modecode = '; mode=display'
        elif mode == 'inline':
            modecode = ''
        else:
            # what happened here?
            raise ValueError("mode must be either 'display' or 'inline'")

        return MathJaxExpr('<html><script type="math/tex{}">'.format(modecode)
                         + ''.join(sage_configurable_latex_macros)
                         + _Latex_prefs._option['macros']
                         + '{}</script></html>'.format(x))

def view(objects, title='Sage', debug=False, sep='', tiny=False,
        pdflatex=None, engine=None, viewer = None, tightpage = None,
        mode='inline', combine_all=False, **kwds):
    r"""nodetex
    Compute a latex representation of each object in objects, compile,
    and display typeset. If used from the command line, this requires
    that latex be installed.

    INPUT:

    -  ``objects`` -- list (or object)

    -  ``title`` -- string (default: ``'Sage'``): title for the
       document

    -  ``debug`` -- bool (default: ``False``): print verbose
       output

    -  ``sep`` -- string (default: ''): separator between
       math objects

    -  ``tiny`` -- bool (default: ``False``): use tiny font.

    -  ``pdflatex`` -- bool (default: ``False``): use pdflatex. This is
       deprecated. Use ``'engine'`` option instead.

    -  ``engine`` -- string or ``None`` (default: ``None``). Can take the
       following values:

       - ``None`` -- the value defined in the LaTeX global preferences
         ``latex.engine()`` is used.

       - ``'pdflatex'`` -- compilation does tex -> pdf

       - ``'xelatex'`` -- compilation does tex -> pdf

       - ``'latex'`` -- compilation first tries tex -> dvi -> png and if an
         error occurs then tries dvi -> ps -> pdf. This is slower than
         ``'pdflatex'`` and known to be broken when overfull hbox are detected.

    -  ``viewer`` -- string or ``None`` (default: ``None``): specify a viewer
       to use; currently the only options are ``None`` and ``'pdf'``.

    -  ``tightpage`` -- bool (default: ``False``): use the LaTeX package
       'preview' with the 'tightpage' option.

    - ``mode`` -- string (default: ``'inline'``): ``'display'`` for
      displaymath or ``'inline'`` for inline math

    - ``combine_all`` -- bool (default: ``False``): If ``combine_all`` is
      ``True`` and the input is a tuple, then it does not return a tuple and
      instead returns a string with all the elements separated by a single
      space.

    OUTPUT:

    Display typeset objects.

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

    If ``pdflatex`` is ``True``, then the latex engine is set to
    pdflatex.

    If the ``engine`` is either ``pdflatex`` or ``xelatex``,  it produces
    a pdf file. Otherwise, it produces a dvi file, and if the program dvipng is
    installed, it checks the dvi file by trying to convert it to a png
    file.  If this conversion fails, the dvi file probably contains
    some postscript special commands or it has other issues which
    might make displaying it a problem; in this case, the file is
    converted to a pdf file, which is then displayed.

    Setting ``viewer`` to ``'pdf'`` forces the use of a separate
    viewer, even in notebook mode. This also sets the latex engine to be
    ``pdflatex`` if the current engine is latex.

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
    usually uses MathJax -- see the next paragraph for the exception --
    to display the output in the notebook. Only the first argument,
    ``objects``, is relevant; the others are ignored. If ``objects``
    is a list, each object is printed on its own line.

    In the notebook, this *does* *not* use MathJax if the LaTeX code
    for ``objects`` contains a string in
    :meth:`latex.mathjax_avoid_list() <Latex.mathjax_avoid_list>`.  In
    this case, it creates and displays a png file.

    EXAMPLES::

        sage: sage.misc.latex.EMBEDDED_MODE = True
        sage: view(3)
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}3</script></html>
        sage: view(3, mode='display')
        <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}3</script></html>
        sage: view((x,2), combine_all=True) # trac 11775
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}x 2</script></html>
        sage: sage.misc.latex.EMBEDDED_MODE = False

    TESTS::

        sage: from sage.misc.latex import _run_latex_, _latex_file_
        sage: g = sage.misc.latex.latex_examples.graph()
        sage: latex.add_to_preamble(r"\usepackage{tkz-graph}")
        sage: file = os.path.join(SAGE_TMP, "temp.tex")
        sage: O = open(file, 'w'); O.write(_latex_file_(g)); O.close()
        sage: _run_latex_(file, engine="pdflatex") # optional - latex
        'pdf'

        sage: latex.extra_preamble('') # reset the preamble

        sage: view(4, engine="garbage")
        Traceback (most recent call last):
        ...
        ValueError: Unsupported LaTeX engine.
        sage: sage.misc.latex.EMBEDDED_MODE = True
        sage: view(4, engine="garbage", viewer="pdf")
        Traceback (most recent call last):
        ...
        ValueError: Unsupported LaTeX engine.

    """
    if tightpage == True:
        latex_options = {'extra_preamble':'\\usepackage[tightpage,active]{preview}\\PreviewEnvironment{page}',
                         'math_left':'\\begin{page}$', 'math_right':'$\\end{page}'}
    else:
        latex_options = {}
    s = _latex_file_(objects, title=title, sep=sep, tiny=tiny, debug=debug, **latex_options)
    if engine is None:
        engine = _Latex_prefs._option["engine"]
    if pdflatex or (viewer == "pdf" and engine == "latex"):
        engine = "pdflatex"
    # notebook
    if EMBEDDED_MODE and viewer is None:
        MathJax_okay = True
        for t in latex.mathjax_avoid_list():
            if s.find(t) != -1:
                MathJax_okay = False
            if not MathJax_okay:
                break
        if MathJax_okay:  # put comma at end of line below?
            print(MathJax().eval(objects, mode=mode, combine_all=combine_all))
        else:
            base_dir = os.path.abspath("")
            png_file = graphics_filename(ext='png')
            png_link = "cell://" + png_file
            png(objects, os.path.join(base_dir, png_file),
                debug=debug, engine=engine)
            print('<html><img src="{}"></html>'.format(png_link))  # put comma at end of line?
        return
    # command line or notebook with viewer
    tmp = tmp_dir('sage_viewer')
    tex_file = os.path.join(tmp, "sage.tex")
    open(tex_file,'w').write(s)
    suffix = _run_latex_(tex_file, debug=debug, engine=engine, png=False)
    if suffix == "pdf":
        from sage.misc.viewer import pdf_viewer
        viewer = pdf_viewer()
    elif suffix == "dvi":
        from sage.misc.viewer import dvi_viewer
        viewer = dvi_viewer()
    else:
        print("Latex error")
        return
    output_file = os.path.join(tmp, "sage." + suffix)
    # this should get changed if we switch the stuff in misc.viewer to
    # producing lists
    if not viewer.startswith('sage-native-execute '):
        viewer = 'sage-native-execute ' + viewer
    if debug:
        print('viewer: "{}"'.format(viewer))
    subprocess.call('%s %s' % (viewer, output_file), shell=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return

def png(x, filename, density=150, debug=False,
        do_in_background=False, tiny=False, pdflatex=True, engine='pdflatex'):
    """
    Create a png image representation of ``x`` and save to the given
    filename.

    INPUT:

    -  ``x`` -- object to be displayed

    -  ``filename`` -- file in which to save the image

    -  ``density`` -- integer (default: 150)

    -  ``debug`` -- bool (default: ``False``): print verbose
       output

    -  ``do_in_background`` -- bool (default: ``False``): Unused,
       kept for backwards compatibility

    -  ``tiny`` -- bool (default: ``False``): use 'tiny' font

    -  ``pdflatex`` -- bool (default: ``True``): use pdflatex. This option is
       deprecated. Use ``engine`` option instead. See below.

    -  ``engine`` -- (default: ``'pdflatex'``) ``'latex'``, ``'pdflatex'``,
       or ``'xelatex'``

    EXAMPLES::

        sage: from sage.misc.latex import png
        sage: png(ZZ[x], os.path.join(SAGE_TMP, "zz.png")) # random - error if no latex
    """
    if not pdflatex:
        engine = "latex"
    import sage.plot.all
    if sage.plot.graphics.is_Graphics(x):
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
                    png=True, engine=engine)
    if e.find("Error") == -1:
        # if no errors, copy png_file to the appropriate place
        shutil.copy(png_file, abs_path_to_png)
    else:
        print("Latex error")
    if debug:
        return s
    return

def coeff_repr(c):
    r"""
    LaTeX string representing coefficients in a linear combination.

    INPUT:

    - ``c`` -- a coefficient (i.e., an element of a ring)

    OUTPUT:

    A string

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

    -  ``symbols`` -- list of symbols

    -  ``coeffs`` -- list of coefficients of the symbols

    OUTPUT:

    A string

    EXAMPLES::

        sage: t = PolynomialRing(QQ, 't').0
        sage: from sage.misc.latex import repr_lincomb
        sage: repr_lincomb(['a', 's', ''], [-t, t - 2, t^12 + 2])
        '-t\\text{\\texttt{a}} + \\left(t - 2\\right)\\text{\\texttt{s}} + \\left(t^{12} + 2\\right)'
        sage: repr_lincomb(['a', 'b'], [1,1])
        '\\text{\\texttt{a}} + \\text{\\texttt{b}}'

    Verify that a certain corner case works (see :trac:`5707` and
    :trac:`5766`)::

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
                except Exception:
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

    INPUT:

    - ``object`` -- Anything

    EXAMPLES::

        sage: sage.misc.latex.print_or_typeset(3)
        3
        sage: sage.misc.latex.EMBEDDED_MODE=True
        sage: sage.misc.latex.print_or_typeset(3)
        3
        sage: TEMP = sys.displayhook
        sage: sys.displayhook = sage.misc.latex.pretty_print
        sage: sage.misc.latex.print_or_typeset(3)
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}3</script></html>
        sage: sage.misc.latex.EMBEDDED_MODE=False
        sage: sys.displayhook = TEMP
    """
    import sys
    if EMBEDDED_MODE and sys.displayhook == pretty_print:
        view(object)
    else:
        print(object)

def pretty_print (*args):
    r"""
    Try to pretty print the arguments in an intelligent way. For graphics
    objects, this returns their default representation. For other
    objects, in the notebook, this calls the :func:`view` command,
    while from the command line, this produces an html string suitable
    for processing by MathJax.

    INPUT:

    - ``objects`` -- The input can be any Sage object, a list or tuple of
      Sage objects, or Sage objects passed in as separate arguments.

    This function is used in the notebook when the "Typeset" button is
    checked.

    EXAMPLES::

        sage: pretty_print(ZZ)  # indirect doctest
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}</script></html>
        sage: pretty_print("Integers = ", ZZ) # trac 11775
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\verb|Integers|\phantom{\verb!x!}\verb|=| \Bold{Z}</script></html>

    To typeset LaTeX code as-is, use :class:`LatexExpr`::

        sage: pretty_print(LatexExpr(r"\frac{x^2 + 1}{x - 2}"))
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\frac{x^2 + 1}{x - 2}</script></html>
    """
    # view s if it is not empty. Used twice.
    def _show_s(s):
        if s != []:
            if EMBEDDED_MODE:
                view(tuple(s), combine_all=True)
            else:
                print(MathJax().eval(tuple(s), mode='inline',
                        combine_all=True))

    import __builtin__
    in_ipython = hasattr(__builtin__, 'get_ipython')
    s = []
    for object in args:
        if object is None:
            continue
        if in_ipython:
            get_ipython().displayhook.update_user_ns(object)
        else:
            __builtin__._ = object

        from sage.plot.plot import Graphics
        from sage.plot.plot3d.base import Graphics3d
        if isinstance(object, (Graphics, Graphics3d)):
            _show_s(s)
            s = []
            print(repr(object))

        else:
            s.append(object)

    _show_s(s)
    return

def pretty_print_default(enable=True):
    r"""
    Enable or disable default pretty printing. Pretty printing means
    rendering things so that MathJax or some other latex-aware front end
    can render real math.

    This function is pretty useless without the notebook, it shoudn't
    be in the global namespace.

    INPUT:

    -  ``enable`` -- bool (optional, default ``True``).  If ``True``, turn on
       pretty printing; if ``False``, turn it off.

    EXAMPLES::

        sage: pretty_print_default(True)
        sage: 'foo'
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\verb|foo|</script></html>
        sage: pretty_print_default(False)
        sage: 'foo'
        'foo'
    """
    import sys
    sys.displayhook.set_display('typeset' if enable else 'simple')


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
                   'Omega',
                   'ast',
                   'bullet',
                   'circ',
                   'times',
                   'star']

def latex_varify(a, is_fname=False):
    r"""
    Convert a string ``a`` to a LaTeX string: if it's an element of
    ``common_varnames``, then prepend a backslash.  If ``a`` consists
    of a single letter, then return it.  Otherwise, return
    either "{\\rm a}" or "\\mbox{a}" if "is_fname" flag is ``True``
    or ``False``.

    INPUT:

    - ``a`` -- string

    OUTPUT:

    A string

    EXAMPLES::

        sage: from sage.misc.latex import latex_varify
        sage: latex_varify('w')
        'w'
        sage: latex_varify('aleph')
        '\\mathit{aleph}'
        sage: latex_varify('aleph', is_fname=True)
        '{\\rm aleph}'
        sage: latex_varify('alpha')
        '\\alpha'
        sage: latex_varify('ast')
        '\\ast'

    TESTS:

        sage: abc = var('abc')
        sage: latex((abc/(abc+1)+42)/(abc-1))  # trac #15870
        \frac{\frac{\mathit{abc}}{\mathit{abc} + 1} + 42}{\mathit{abc} - 1}
    """
    if a in common_varnames:
        return "\\" + a
    elif len(a) == 0:
        return ''
    elif len(a) == 1:
        return a
    elif is_fname is True:
        return '{\\rm %s}'%a
    else:
        return '\\mathit{%s}'%a

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
        '\\mathit{abc}'
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
        '\\mathit{nothing}_{1}'
        sage: latex_variable_name('nothing1', is_fname=True)
        '{\\rm nothing}_{1}'
        sage: latex_variable_name('nothing_abc')
        '\\mathit{nothing}_{\\mathit{abc}}'
        sage: latex_variable_name('nothing_abc', is_fname=True)
        '{\\rm nothing}_{{\\rm abc}}'
        sage: latex_variable_name('alpha_beta_gamma12')
        '\\alpha_{\\beta_{\\gamma_{12}}}'
        sage: latex_variable_name('x_ast')
        'x_{\\ast}'

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

    The classes here only have ``__init__``, ``_repr_``, and ``_latex_``
    methods.

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

            EXAMPLES::

                sage: from sage.misc.latex import latex_examples
                sage: G = latex_examples.graph()
                sage: len(G._repr_()) > 300
                True
            """
            return r"""LaTeX example for testing display of graphs.

To use, first try calling 'view' on this object -- it won't work.
Now, make sure that you have the most recent version of the TeX
package pgf installed, along with the LaTeX package tkz-graph.  Run
'latex.add_to_preamble("\\usepackage{tkz-graph}")', and try viewing it
again.  From the command line, this should pop open a nice window with
a picture of a graph.  In the notebook, it still won't work.  Finally,
run 'latex.add_to_mathjax_avoid_list("tikzpicture")' and try again
from the notebook -- you should get a nice picture.

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

To use, first try calling 'view' on this object -- it won't work. Now,
make sure that you have the most recent version of the TeX package
pstricks installed.  Run 'latex.add_to_preamble("\\usepackage{pstricks}")'
and try viewing it again. Call 'view' with the option `engine='latex'`
-- the default behavior is to use pdflatex, which doesn't work with
pstricks.  From the command line, this should pop open a nice window
with a picture of forces acting on a mass on a pendulum.  In the
notebook, it still won't work, so run
'latex.add_to_mathjax_avoid_list("pspicture")' and try again -- you
should get a nice picture."""

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
'latex.add_to_mathjax_avoid_list("xygraph")' and try again -- you
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
'latex.add_to_mathjax_avoid_list("xymatrix")' and try again -- you
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
