r"""
LaTeX macros

AUTHORS:

- John H. Palmieri (2009-03)

The code here sets up LaTeX macro definitions for use in the
documentation. To add a macro, modify the list ``macros``, near the
end of this file, and then run 'sage -b'. The entries in this list are
used to produce ``sage_latex_macros``, a list of strings of the form
'\\newcommand...', and ``sage_js_macros``, a list of strings of the
form 'jsMath.Macro...'.  The LaTeX macros are produced using the
``_latex_`` method for each Sage object listed in ``macros``, and the
jsMath macros are produced from the LaTeX macros.  The list of LaTeX
macros is used in the file ``SAGE_ROOT/devel/sage/doc/common/conf.py``
to add to the preambles of both the LaTeX file used to build the PDF
version of the documentation and the LaTeX file used to build the html
version.  The list of jsMath macros is used in the file
``sage/server/notebook/notebook.py`` to define jsMath macros for use
in the live documentation (and also in the notebook).

Any macro defined here may be used in docstrings or in the tutorial
(or other pieces of documentation).  In a docstring, for example,
"\ZZ" in backquotes (demarking math mode) will appear as "ZZ" in
interactive help, but will be typeset as "\\mathbf{Z}" in the
reference manual.

More details on the list ``macros``: the entries are lists or tuples
of the form ``[name]`` or ``[name, arguments]``, where ``name`` is a
string and ``arguments`` consists of valid arguments for the Sage
object named ``name``.  For example, ``["ZZ"]`` and ``["GF", 2]``
produce the LaTeX macros '\\newcommand{\\ZZ}{\\mathbf{Z}}' and
'\\newcommand{\\GF}[1]{\\mathbf{F}_{#1}}', respectively.  (For the
second of these, ``latex(GF(2))`` is called and the string '2' gets
replaced by '#1', so ``["GF", 17]`` would have worked just as well.
``["GF", p]`` would have raised an error, though, because ``p`` is not
defined, and ``["GF", 4]`` would have raised an error, because to
define the field with four elements in Sage, you also need to specify
the name of a generator.)

To see evidence of the results of the code here, run ``sage -docbuild
tutorial latex`` (for example), and look at the resulting LaTeX file in
``SAGE_ROOT/sage/doc/output/latex/en/tutorial/``.  The preamble should
contain '\newcommand' lines for each of the entries in ``macros``.
"""

def produce_latex_macro(name, *sample_args):
    r"""
    Produce a string defining a LaTeX macro.

    INPUT:

    -  ``name`` - name of macro to be defined, also name of corresponding Sage object
    -  ``sample_args`` - (optional) sample arguments for this Sage object

    EXAMPLES::

        sage: from sage.misc.latex_macros import produce_latex_macro
        sage: produce_latex_macro('ZZ')
        '\\newcommand{\\ZZ}{\\mathbf{Z}}'

    If the Sage object takes arguments, then the LaTeX macro will
    accept arguments as well. You must pass valid arguments, which
    will then be converted to #1, #2, etc. in the macro
    definition. The following allows the use of "\GF{p^n}", for
    example::

         sage: produce_latex_macro('GF', 37)
         '\\newcommand{\\GF}[1]{\\mathbf{F}_{#1}}'

    If the Sage object is not in the global name space, describe it
    like so::

         sage: produce_latex_macro('sage.rings.finite_field.FiniteField', 3)
         '\\newcommand{\\FiniteField}[1]{\\mathbf{F}_{#1}}'
    """
    from sage.misc.latex import latex
    names_split = name.rsplit('.', 1)
    if len(names_split) == 1:
        module = 'sage.all'
        real_name = names_split[0]
    else:
        module, real_name = names_split
    newcommand = '\\newcommand{\\' + real_name + '}'
    count = 0
    args = "("
    for x in sample_args:
        count += 1
        args += str(x) + ','
    args += ')'
    exec('from ' + module + ' import ' + real_name)
    if count > 0:
        defn = '[' + str(count) + ']{'
        defn += eval('str(latex(' + real_name + args + '))') + '}'
    else:
        defn = '{' + eval('str(latex(' + real_name + '))') + '}'
    count = 0
    for x in sample_args:
        count += 1
        defn = defn.replace(str(x), "#" + str(count))
    return newcommand + defn

def convert_latex_macro_to_jsmath(macro):
    r"""
    This converts a LaTeX macro definition (\newcommand...) to a
    jsMath macro definition (jsMath.Macro...).

    INPUT:

    -  ``macro`` - LaTeX macro definition

    See the web page
    http://www.math.union.edu/~dpvc/jsMath/authors/macros.html for a
    description of the format for jsMath macros.

    EXAMPLES::

        sage: from sage.misc.latex_macros import convert_latex_macro_to_jsmath
        sage: convert_latex_macro_to_jsmath('\\newcommand{\\ZZ}{\\mathbf{Z}}')
        "jsMath.Macro('ZZ','\\\\mathbf{Z}')"
        sage: convert_latex_macro_to_jsmath('\\newcommand{\\GF}[1]{\\mathbf{F}_{#1}}')
        "jsMath.Macro('GF','\\\\mathbf{F}_{#1}',1)"
    """
    left_bracket = macro.find('[')
    right_bracket = macro.find('[')
    if left_bracket >= 0:
        right_bracket = macro.find(']')
        num_args = macro[left_bracket+1:right_bracket]
    else:
        num_args = 0
    start_name = macro.find('{') + 1  # add one to go past the backslash
    end_name = macro.find('}')
    name = macro[start_name+1:end_name]
    start_defn = macro.find('{', end_name)
    end_defn = macro.rfind('}')
    defn = macro[start_defn+1: end_defn].replace('\\', '\\\\')
    if num_args > 0:
        args_str = "," + str(num_args)
    else:
        args_str = ""
    return "jsMath.Macro('" + name + "','" + defn + "'" + args_str + ")"

# To add a new macro for use in the Sage documentation, add a list or
# tuple to the following list.  Each list (or tuple) should have the
# form [name, arguments], which will be passed to the function
# produce_latex_macro: see that for more documentation.
#
# To see the results of this, run 'sage -docbuild tutorial latex' (for
# example -- you could replace 'tutorial' with your favorite piece of
# documentation), and look at the resulting tex file in
# SAGE_DOC/output/latex/en/tutorial.  The preamble should contain
# \newcommand's for each of the entries here.
macros = [["ZZ"],
          ["RR"],
          ["CC"],
          ["QQ"],
          ["QQbar"],
          ["GF", 2],
          ["Zp", 2],
          ["Qp", 2],
          ["Zmod", 2],
          ["CDF"],
          ["CIF"],
          ["CLF"],
          ["RDF"],
          ["RIF"],
          ["RLF"],
          ["RQDF"],
          ["CFF"],
          ]

sage_latex_macros = [produce_latex_macro(*x) for x in macros]
sage_jsmath_macros = [convert_latex_macro_to_jsmath(m) for m in sage_latex_macros]
