r"""
Sage's IPython Modifications

This module contains all of Sage's customizations to the IPython
interpreter.  These changes consist of the following major components:

  - :class:`SageTerminalApp`
  - :class:`SageInteractiveShell`
  - :class:`SageTerminalInteractiveShell`
  - :func:`interface_shell_embed`

SageTerminalApp
---------------

This is the main application object.  It is used by the
``$SAGE_LOCAL/bin/sage-ipython`` script to start the Sage
command-line.  It's primary purpose is to

  - Initialize the :class:`SageTerminalInteractiveShell`.

  - Provide default configuration options for the shell, and its
    subcomponents.  These work with (and can be overridden by)
    IPython's configuration system.

  - Load the Sage ipython extension (which does things like preparsing,
    add magics, etc.).

  - Provide a custom :class:`SageCrashHandler` to give the user
    instructions on how to report the crash to the Sage support
    mailing list.

SageInteractiveShell
--------------------

The :class:`SageInteractiveShell` object is the object responsible for
accepting input from the user and evaluating it.  From the command-line,
this object can be retrieved by running::

    sage: shell = get_ipython()   # not tested

Any input is preprocessed and evaluated inside the ``shell.run_cell``
method. If the command line processing does not do what you want it to
do, you can step through it in the debugger::

    sage: %debug shell.run_cell('?')        # not tested

The :class:`SageInteractiveShell` provides the following
customizations:

  - Modify the libraries before calling system commands. See
    :meth:`~SageInteractiveShell.system_raw`.

SageTerminalInteractiveShell
----------------------------

The :class:`SageTerminalInteractiveShell` is a close relative of
:class:`SageInteractiveShell` that is specialized for running in a
terminal. In particular, running commands like ``!ls`` will directly
write to stdout. Technically, the ``system`` attribute will point to
``system_raw`` instead of ``system_piped``.

Interface Shell
---------------

The function :func:`interface_shell_embed` takes a
:class:`~sage.interfaces.interface.Interface` object and returns an
embeddable IPython shell which can be used to directly interact with
that shell.  The bulk of this functionality is provided through
:class:`InterfaceShellTransformer`.

TESTS:

Check that Cython source code appears in tracebacks::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: print("dummy line"); shell.run_cell('1/0') # see #25320 for the reason of the `...` and the dummy line in this test
    dummy line
    ...
    ZeroDivisionError...Traceback (most recent call last)
    <ipython-input-...> in <module>...
    ----> 1 Integer(1)/Integer(0)
    .../sage/rings/integer.pyx in sage.rings.integer.Integer...div...
    ...
    -> ...                  raise ZeroDivisionError("rational division by zero")
       ....:            x = <Rational> Rational.__new__(Rational)
       ....:            mpq_div_zz(x.value, ....value, (<Integer>right).value)
    <BLANKLINE>
    ZeroDivisionError: rational division by zero
    sage: shell.quit()

Test prompt transformer::

    sage: from sage.repl.interpreter import SagePromptTransformer
    sage: spt = SagePromptTransformer
    sage: spt(["sage: 2 + 2"])
    ['2 + 2']
    sage: spt([''])
    ['']
    sage: spt(["....: 2+2"])
    ['2+2']

This should strip multiple prompts: see :trac:`16297`::

    sage: spt(["sage:   sage: 2+2"])
    ['2+2']
    sage: spt(["   sage: ....: 2+2"])
    ['2+2']

The prompt contains a trailing space. Extra spaces between the
last prompt and the remainder should not be stripped::

    sage: spt(["   sage: ....:    2+2"])
    ['   2+2']

We test that the input transformer is enabled on the Sage command
line::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell('sage: a = 123')              # single line
    sage: shell.run_cell('sage: a = [\n... 123]')      # old-style multi-line
    sage: shell.run_cell('sage: a = [\n....: 123]')    # new-style multi-line

We test that :trac:`16196` is resolved::

    sage: shell.run_cell('    sage: 1+1')
    2
    sage: shell.quit()
"""

# ****************************************************************************
#       Copyright (C) 2004-2012 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import re
from traitlets import Bool, Type

from sage.repl.preparse import preparse, containing_block
from sage.repl.prompts import InterfacePrompts
from sage.repl.configuration import sage_ipython_config, SAGE_EXTENSION

from IPython.core.interactiveshell import InteractiveShell
from IPython.core.prefilter import PrefilterTransformer
from IPython.core.crashhandler import CrashHandler
from IPython.core.inputtransformer2 import (PromptStripper,
                                            TokenTransformBase,
                                            tokenize)

from IPython.terminal.interactiveshell import TerminalInteractiveShell
from IPython.terminal.embed import InteractiveShellEmbed
from IPython.terminal.ipapp import TerminalIPythonApp, IPAppCrashHandler


# TODO: This global variable _do_preparse should be associated with an
# IPython InteractiveShell as opposed to a global variable in this
# module.
_do_preparse = True


def preparser(on=True):
    """
    Turn on or off the Sage preparser.

    :keyword on: if True turn on preparsing; if False, turn it off.
    :type on: bool

    EXAMPLES::

        sage: 2/3
        2/3
        sage: preparser(False)
        sage: 2/3  # not tested since doctests are always preparsed
        0
        sage: preparser(True)
        sage: 2^3
        8
    """
    global _do_preparse
    _do_preparse = on is True


##############################
# Sage[Terminal]InteractiveShell
##############################
class SageShellOverride(object):
    """
    Mixin to override methods in IPython's [Terminal]InteractiveShell
    classes.
    """

    def show_usage(self):
        """
        Print the basic Sage usage.

        This method ends up being called when you enter ``?`` and
        nothing else on the command line.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('?')
            Welcome to Sage ...
            sage: shell.quit()
        """
        from sage.misc.sagedoc import help
        help()

    def system_raw(self, cmd):
        """
        Run a system command.

        This is equivalent to the sage-native-execute shell script.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.system_raw('false')
            sage: shell.user_ns['_exit_code'] > 0
            True
            sage: shell.system_raw('true')
            sage: shell.user_ns['_exit_code']
            0
            sage: shell.system_raw('R --version')   # optional - r
            R version ...
            sage: shell.user_ns['_exit_code']       # optional - r
            0
            sage: shell.quit()
        """
        return super(SageShellOverride, self).system_raw(cmd)


class SageNotebookInteractiveShell(SageShellOverride, InteractiveShell):
    """
    IPython Shell for the Sage IPython Notebook

    The doctests are not tested since they would change the current
    rich output backend away from the doctest rich output backend.

    EXAMPLES::

        sage: from sage.repl.interpreter import SageNotebookInteractiveShell
        sage: SageNotebookInteractiveShell()   # not tested
        <sage.repl.interpreter.SageNotebookInteractiveShell object at 0x...>
    """

    def init_display_formatter(self):
        """
        Switch to the Sage IPython notebook rich output backend

        EXAMPLES::

            sage: from sage.repl.interpreter import SageNotebookInteractiveShell
            sage: SageNotebookInteractiveShell().init_display_formatter()   # not tested
        """
        from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
        backend = BackendIPythonNotebook()
        backend.get_display_manager().switch_backend(backend, shell=self)


class SageTerminalInteractiveShell(SageShellOverride, TerminalInteractiveShell):
    """
    IPython Shell for the Sage IPython Commandline Interface

    The doctests are not tested since they would change the current
    rich output backend away from the doctest rich output backend.

    EXAMPLES::

        sage: from sage.repl.interpreter import SageTerminalInteractiveShell
        sage: SageTerminalInteractiveShell()   # not tested
        <sage.repl.interpreter.SageNotebookInteractiveShell object at 0x...>
    """

    def init_display_formatter(self):
        """
        Switch to the Sage IPython commandline rich output backend

        EXAMPLES::

            sage: from sage.repl.interpreter import SageTerminalInteractiveShell
            sage: SageTerminalInteractiveShell().init_display_formatter()   # not tested
        """
        from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
        backend = BackendIPythonCommandline()
        backend.get_display_manager().switch_backend(backend, shell=self)


class SageTestShell(SageShellOverride, TerminalInteractiveShell):
    """
    Test Shell

    Care must be taken in these doctests to quit the test shell in
    order to switch back the rich output display backend to the
    doctest backend.

    EXAMPLES::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell();  shell
        <sage.repl.interpreter.SageTestShell object at 0x...>
        sage: shell.quit()
    """

    def init_display_formatter(self):
        """
        Switch to the Sage IPython commandline rich output backend

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell();  shell
            <sage.repl.interpreter.SageTestShell object at 0x...>
            sage: shell.quit()
            sage: shell.init_display_formatter()
            sage: shell.quit()
        """
        from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
        self._ipython_backend = backend = BackendIPythonCommandline()
        self._display_manager = backend.get_display_manager()
        self._doctest_backend = self._display_manager.switch_backend(backend, shell=self)

    def quit(self):
        """
        Quit the test shell.

        To make the test shell as realistic as possible, we switch to
        the
        :class:`~sage.repl.rich_output.backend_ipython.BackendIPythonCommandline`
        display backend. This method restores the previous display
        backend, which is the
        :class:`~sage.repl.rich_output.backend_doctest.BackendDoctest`
        during doctests.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.repl.rich_output import get_display_manager
            sage: get_display_manager()
            The Sage display manager using the doctest backend

            sage: shell = get_test_shell()
            sage: get_display_manager()
            The Sage display manager using the IPython command line backend

            sage: shell.quit()
            sage: get_display_manager()
            The Sage display manager using the doctest backend
        """
        self._display_manager.switch_backend(self._doctest_backend)

    def _restart(self):
        """
        Restart the test shell (after :meth:`quit`).

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.quit()
            sage: shell._restart()
            sage: shell.quit()
            sage: from sage.repl.rich_output import get_display_manager
            sage: get_display_manager()
            The Sage display manager using the doctest backend
        """
        self._display_manager.switch_backend(self._ipython_backend, shell=self)

    def run_cell(self, *args, **kwds):
        """
        Run IPython cell

        Starting with IPython-3.0, this returns an success/failure
        information. Since it is more convenient for doctests, we
        ignore it.


        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: rc = shell.run_cell('2^50')
            1125899906842624
            sage: rc is None
            True
            sage: shell.quit()
        """
        super(SageTestShell, self).run_cell(*args, **kwds)


###################################################################
# Transformers used in the SageInputSplitter
###################################################################

def SagePreparseTransformer(lines):
    r"""
    EXAMPLES::

        sage: from sage.repl.interpreter import SagePreparseTransformer
        sage: SagePreparseTransformer(['1+1r+2.3^2.3r\n'])
        ["Integer(1)+1+RealNumber('2.3')**2.3\n"]
        sage: preparser(False)
        sage: SagePreparseTransformer(['2.3^2\n'])
        ['2.3^2\n']
        sage: preparser(True)

    Numbers are parsed as Sage's numbers unless an ``r`` is appended::

        sage: from sage.repl.interpreter  import SagePreparseTransformer
        sage: SagePreparseTransformer(["2+2r + 4"])
        ['Integer(2)+2 + Integer(4)']
        sage: SagePreparseTransformer(["def foo():\n", "    return 2"])
        ['def foo():\n', '    return Integer(2)']

    Multiline-strings are left as they are::

        sage: SagePreparseTransformer(["def foo():\n", "    return a + '''\n", "    2 - 3r\n", "    '''"])
        ['def foo():\n', "    return a + '''\n", '    2 - 3r\n', "    '''"]

    Preparses ``[0,2,..,n]`` notation::

        sage: SagePreparseTransformer(["for i in [2 .. 5 ..a]"])
        ['for i in (ellipsis_range(Integer(2) ,Ellipsis, Integer(5) ,Ellipsis,a))']
        sage: SagePreparseTransformer(["for i in (2 .. 5r)"])
        ['for i in (ellipsis_iter(Integer(2) ,Ellipsis, 5))']

    Preparses generator access::

        sage: SagePreparseTransformer(["K = QuadraticField(2)\n", "print(K.0)"])
        ['K = QuadraticField(Integer(2))\n', 'print(K.gen(0))']

    Preparses implicit multiplication::

        sage: SagePreparseTransformer(["2a"])
        ['2a']
        sage: implicit_multiplication(True)
        sage: SagePreparseTransformer(["2a"])
        ['Integer(2)*a']
        sage: implicit_multiplication(False)

    Replaces ``^`` by exponentiation and ``^^`` as bitwise xor::

        sage: SagePreparseTransformer(["x^2"])
        ['x**Integer(2)']
        sage: SagePreparseTransformer(["x^^2"])
        ['x^Integer(2)']


    TESTS:

    Check that syntax errors in the preparser do not crash IPython,
    see :trac:`14961`. ::

        sage: preparser(True)
        sage: bad_syntax = "R.<t> = QQ{]"
        sage: preparse(bad_syntax)
        'R = QQ{]; (t,) = R._first_ngens(1)'
        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell(bad_syntax)
          File "...", line 1
            R = QQ{]; (t,) = R._first_ngens(1)
                  ^
        SyntaxError: invalid syntax
        <BLANKLINE>
        sage: shell.quit()

    Make sure the quote state is carried over across subsequent lines in order
    to avoid interfering with multi-line strings, see :trac:`30417`. ::

        sage: SagePreparseTransformer(["'''\n", 'abc-1-2\n', "'''\n"])
        ["'''\n", 'abc-1-2\n', "'''\n"]
        sage: # instead of ["'''\n", 'abc-Integer(1)-Integer(2)\n', "'''\n"]

    .. NOTE::

        IPython may call this function more than once for the same input lines.
        So when debugging the preparser, print outs may be duplicated. If using
        IPython >= 7.17, try:
        ``sage.repl.interpreter.SagePreparseTransformer.has_side_effects = True``

    """
    if _do_preparse:
        # IPython ensures the input lines end with a newline, and it expects
        # the same of the output lines.
        lines = preparse(''.join(lines)).splitlines(keepends=True)
    return lines


SagePromptTransformer = PromptStripper(prompt_re=re.compile(r'^(\s*(:?sage: |\.\.\.\.: ))+'))


class SageGenConstructionTransformer(TokenTransformBase):
    r"""
    Transform Sage's construction with generators.

    TESTS::

        sage: from IPython import get_ipython
        sage: ip = get_ipython()

    Vanilla::

        sage: ip.transform_cell('''
        ....: R.<x> = ZZ['x']''')
        "R = ZZ['x']; (x,) = R._first_ngens(1)\n"
        sage: ip.transform_cell('''
        ....: R.<x,y> = ZZ['x,y']''')
        "R = ZZ['x,y']; (x, y,) = R._first_ngens(2)\n"

    No square brackets::

        sage: ip.transform_cell('''
        ....: R.<x> = PolynomialRing(ZZ, 'x')''')
        "R = PolynomialRing(ZZ, 'x', names=('x',)); (x,) = R._first_ngens(1)\n"
        sage: ip.transform_cell('''
        ....: R.<x,y> = PolynomialRing(ZZ, 'x,y')''')
        "R = PolynomialRing(ZZ, 'x,y', names=('x', 'y',)); (x, y,) = R._first_ngens(2)\n"

    Names filled in::

        sage: ip.transform_cell('''
        ....: R.<x> = ZZ[]''')
        "R = ZZ['x']; (x,) = R._first_ngens(1)\n"
        sage: ip.transform_cell('''
        ....: R.<x,y> = ZZ[]''')
        "R = ZZ['x', 'y']; (x, y,) = R._first_ngens(2)\n"

    Names given not the same as generator names::

        sage: ip.transform_cell('''
        ....: R.<x> = ZZ['y']''')
        "R = ZZ['y']; (x,) = R._first_ngens(1)\n"
        sage: ip.transform_cell('''
        ....: R.<x,y> = ZZ['u,v']''')
        "R = ZZ['u,v']; (x, y,) = R._first_ngens(2)\n"

    Number fields::

        sage: ip.transform_cell('''
        ....: K.<a> = QuadraticField(2)''')
        "K = QuadraticField(Integer(2), names=('a',)); (a,) = K._first_ngens(1)\n"
        sage: ip.transform_cell('''
        ....: K.<a> = QQ[2^(1/3)]''')
        'K = QQ[Integer(2)**(Integer(1)/Integer(3))]; (a,) = K._first_ngens(1)\n'
        sage: ip.transform_cell('''
        ....: K.<a, b> = QQ[2^(1/3), 2^(1/2)]''')
        'K = QQ[Integer(2)**(Integer(1)/Integer(3)), Integer(2)**(Integer(1)/Integer(2))]; (a, b,) = K._first_ngens(2)\n'

    Just the ``.<>`` notation::

        sage: ip.transform_cell('''
        ....: R.<x> = ZZx''')
        'R = ZZx; (x,) = R._first_ngens(1)\n'
        sage: ip.transform_cell('''
        ....: R.<x, y> = a+b''')
        'R = a+b; (x, y,) = R._first_ngens(2)\n'
        sage: ip.transform_cell('''
        ....: A.<x,y,z>=FreeAlgebra(ZZ,3)''')
        "A=FreeAlgebra(ZZ,Integer(3), names=('x', 'y', 'z',)); (x, y, z,) = A._first_ngens(3)\n"

    Ensure we do not eat too much::

        sage: ip.transform_cell('''
        ....: R.<x, y> = ZZ;2''')
        'R = ZZ; (x, y,) = R._first_ngens(2);Integer(2)\n'
        sage: ip.transform_cell('''
        ....: R.<x, y> = ZZ['x,y'];2''')
        "R = ZZ['x,y']; (x, y,) = R._first_ngens(2);Integer(2)\n"
        sage: ip.transform_cell('''
        ....: F.<b>, f, g = S.field_extension()''')
        "F, f, g = S.field_extension(names=('b',)); (b,) = F._first_ngens(1)\n"

    Multiple lines::

        sage: ip.transform_cell('''
        ....: K.<a, b,c> = some_large_field(
        ....:     input1,
        ....:     input2,
        ....:     intput3)''')
        "K = some_large_field(\n    input1,\n    input2,\n    intput3, names=('a', 'b', 'c',)); (a, b, c,) = K._first_ngens(3)\n"
        sage: ip.transform_cell('''
        ....: def foo(L):
        ....:     K.<a,b,c>= L
        ....:     return K, a, b, c''')
        'def foo(L):\n    K= L; (a, b, c,) = K._first_ngens(3)\n    return K, a, b, c\n'
        sage: ip.transform_cell('''
        ....: def foo(L, M):
        ....:     K1.<a,b,c>, K2.<d,e> = L, M
        ....:     return a,b,c,d,e''')
        'def foo(L, M):\n    K1, K2 = L, M; (d, e,) = K2._first_ngens(2); (a, b, c,) = K1._first_ngens(3)\n    return a,b,c,d,e\n'
        sage: ip.transform_cell('''
        ....: def p in primes(2, 20):
        ....:     K.<a> = QuadraticField(p); print(a)''')
        "def p in primes(Integer(2), Integer(20)):\n    K = QuadraticField(p, names=('a',)); (a,) = K._first_ngens(1); print(a)\n"

    See :trac:`16731`::

        sage: ip.transform_cell('''
        ....: R.<x> = ''')
        'R = ; (x,) = R._first_ngens(1)\n'

    Check support for unicode characters (:trac:`29278`)::

        sage: ip.transform_cell('''
        ....: Ω.<λ,μ> = QQ[]''')
        "Ω = QQ['λ', 'μ']; (λ, μ,) = Ω._first_ngens(2)\n"

    Check that :trac:`30953` is fixed::

        sage: K.<a> = QuadraticField(2 +  # some comment
        ....:                        1); K, a
        (Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?,
        a)
    """
    priority = 25

    def __init__(self, del_start, del_end, keyword_pos, argument_pos, insert_pos, name, gens):
        """
        INPUT:

        - ``del_start`` -- start of the ``<``

        - ``del_end`` -- end of the ``>``

        - ``keyword_pos`` -- position to insert generator names as keyword or ``None``

        - ``argument_pos`` -- position to insert generator names as arguments or ``None``

        - ``names_pos`` -- position to insert the additional command to declare the generators

        - ``name`` -- name of the field or similar

        - ``gens`` -- names of the generators

        TESTS::

            sage: from IPython.core.inputtransformer2 import make_tokens_by_line
            sage: from sage.repl.interpreter import SageGenConstructionTransformer
            sage: lines = ['K.<a,b> = QQ[]']
            sage: tokens = make_tokens_by_line(lines)
            sage: S = SageGenConstructionTransformer.find(tokens)  # indirect doctest
            sage: S.transform(lines)
            ["K = QQ['a', 'b']; (a, b,) = K._first_ngens(2)"]
        """
        super().__init__(del_start)
        self.name = name
        self.gens = gens

        # Shift from 1-index to 0-indexed.
        self.del_start = (del_start[0] - 1, del_start[1])
        self.del_end = (del_end[0] - 1, del_end[1])
        self.insert_pos = (insert_pos[0] - 1, insert_pos[1])
        self.keyword_pos = (keyword_pos[0] - 1, keyword_pos[1]) if keyword_pos is not None else None
        self.argument_pos = (argument_pos[0] - 1, argument_pos[1]) if argument_pos is not None else None

    @classmethod
    def find(cls, tokens_by_line):
        r"""
        Find the first construction with generators.

        EXAMPLES::

            sage: from IPython.core.inputtransformer2 import make_tokens_by_line
            sage: from sage.repl.interpreter import SageGenConstructionTransformer
            sage: lines = ['for p in primes(2, 20):\n', '    K.<a> = QuadraticField(p)\n', '    print(a)']
            sage: tokens = make_tokens_by_line(lines)
            sage: S = SageGenConstructionTransformer.find(tokens)
            sage: S.transform(lines)
            ['for p in primes(2, 20):\n',
             "    K = QuadraticField(p, names=('a',)); (a,) = K._first_ngens(1)\n",
             '    print(a)']
        """
        for line in tokens_by_line:
            name = None
            gens = []
            for i, token in enumerate(line[2:], 2):
                if token.string == '<' and line[i-1].string == '.' and token.start == line[i-1].end and line[i-2].type == tokenize.NAME:
                    if name is None:
                        name = line[i-2].string
                        del_start = line[i-1].start
                    else:
                        # Do not transform syntax errors.
                        break
                elif token.string == '<':
                    # Do not transform syntax errors.
                    break
                elif name is None:
                    continue
                elif token.type == tokenize.NAME:
                    gens.append(token.string)
                elif token.string[:1] == '>':
                    del_end = (token.start[0], token.start[1] + 1)

                    if del_end[0] != del_start[0]:
                        # The generators must be written in one line.
                        break
                    if gens == []:
                        # At least one generator needed.
                        break

                    # Find the position to insert the declaration of the generators.
                    ix = i + 1
                    while not line[ix].string == ';' and not line[ix].type == tokenize.NEWLINE and ix + 1 < len(line):
                        ix += 1
                    while line[ix - 1].type == tokenize.COMMENT:
                        ix -= 1
                    insert_pos = line[ix].start

                    # See if the names of the generators need to be given as keyword or argument.
                    keyword_pos = None
                    argument_pos = None
                    ix -= 1
                    while line[ix].type in {tokenize.COMMENT, tokenize.NL}:
                        ix -= 1
                    if line[ix].string == ')':
                        keyword_pos = True
                    elif line[ix].string == ']':
                        while line[ix-1].string == ']':
                            ix -= 1
                        argument_pos = True

                    if argument_pos or keyword_pos:
                        ix -= 1
                        while line[ix].type in {tokenize.COMMENT, tokenize.NL}:
                            ix -= 1
                        if keyword_pos:
                            keyword_pos = line[ix].end
                        else:
                            if line[ix].string == '[':
                                argument_pos = line[ix].end
                            else:
                                argument_pos = None

                    return cls(del_start, del_end, keyword_pos, argument_pos, insert_pos, name, gens)
                elif token.string != ',':
                    # Do not transform syntax errors.
                    break

    def transform(self, lines):
        r"""
        Transform the first construction with generators.

        EXAMPLES::

            sage: from IPython.core.inputtransformer2 import make_tokens_by_line
            sage: from sage.repl.interpreter import SageGenConstructionTransformer
            sage: lines = ['F.<b>, f, g = S.field_extension()']
            sage: tokens = make_tokens_by_line(lines)
            sage: S = SageGenConstructionTransformer.find(tokens)
            sage: S.transform(lines)
            ["F, f, g = S.field_extension(names=('b',)); (b,) = F._first_ngens(1)"]
        """
        lines = [l for l in lines]

        new_command = '; (' + ', '.join(self.gens) + ',) = ' + self.name \
                      + "._first_ngens({})".format(len(self.gens))
        lines[self.insert_pos[0]] = lines[self.insert_pos[0]][:self.insert_pos[1]] \
                                    + new_command \
                                    + lines[self.insert_pos[0]][self.insert_pos[1]:]

        if self.argument_pos or self.keyword_pos:
            names = "'" + "', '".join(self.gens) + "'"
            if self.keyword_pos:
                names_pos = self.keyword_pos
                names = "names=(" + names
                if len(self.gens):
                    names += ","
                names += ")"
            else:
                names_pos = self.argument_pos
            if lines[names_pos[0]][names_pos[1] - 1] not in ('(', '['):
                names = ", " + names
            lines[names_pos[0]] = lines[names_pos[0]][:names_pos[1]] + names \
                                  + lines[names_pos[0]][names_pos[1]:]

        lines[self.del_start[0]] = lines[self.del_start[0]][:self.del_start[1]] \
                                   + lines[self.del_start[0]][self.del_end[1]:]

        return lines


class SageCalculusTransformer(TokenTransformBase):
    r"""
    Supports calculus-like function assignment, e.g., transforms::

       f(x,y,z) = sin(x^3 - 4*y) + y^x

    into::

       __tmp__=var("x,y,z")
       f = symbolic_expression(sin(x**3 - 4*y) + y**x).function(x,y,z)

    TESTS:

        sage: from IPython import get_ipython
        sage: ip = get_ipython()
        sage: ip.transform_cell('''
        ....: f(x) = x^3 - x''')
        '__tmp__ = var("x"); __tmpf__ = x**Integer(3) - x; f = symbolic_expression(__tmpf__).function(x)\n'
        sage: ip.transform_cell('''
        ....: f(u,v) = u - v''')
        '__tmp__ = var("u,v"); __tmpf__ = u - v; f = symbolic_expression(__tmpf__).function(u,v)\n'
        sage: ip.transform_cell('''
        ....: f(x) =-5''')
        '__tmp__ = var("x"); __tmpf__ =-Integer(5); f = symbolic_expression(__tmpf__).function(x)\n'
        sage: ip.transform_cell('''
        ....: f(x) -= 5''')
        'f(x) -= Integer(5)\n'
        sage: ip.transform_cell('''
        ....: f(x_1, x_2) = x_1^2 - x_2^2''')
        '__tmp__ = var("x_1,x_2"); __tmpf__ = x_1**Integer(2) - x_2**Integer(2); f = symbolic_expression(__tmpf__).function(x_1,x_2)\n'

        sage: ip.transform_cell('''
        ....: f(t,s)=t^2''')
        '__tmp__ = var("t,s"); __tmpf__=t**Integer(2); f = symbolic_expression(__tmpf__).function(t,s)\n'


        sage: ip.transform_cell('''
        ....: f(x, y) = x^3 - y''')
        '__tmp__ = var("x,y"); __tmpf__ = x**Integer(3) - y; f = symbolic_expression(__tmpf__).function(x,y)\n'
        sage: ip.transform_cell('''
        ....: μ(x,y) = (x^3 -
        ....:           y)''')
        '__tmp__ = var("x,y"); __tmpμ__ = (x**Integer(3) -\n          y); μ = symbolic_expression(__tmpμ__).function(x,y)\n'
        sage: ip.transform_cell('''
        ....: f(x,
        ....:  y) = (x^3 -
        ....:           y)''')
        '__tmp__ = var("x,y"); __tmpf__ = (x**Integer(3) -\n          y); f = symbolic_expression(__tmpf__).function(x,y)\n'
        sage: ip.transform_cell('''
        ....: def foo:
        ....:     f(x,
        ....:       y) = (x^3 -
        ....:             y); return f''')
        'def foo:\n    __tmp__ = var("x,y"); __tmpf__ = (x**Integer(3) -\n            y); f = symbolic_expression(__tmpf__).function(x,y); return f\n'

    Check that :trac:`30953` is fixed::

        sage: f(x) = (x +
        ....:         1); f
        x |--> x + 1
        sage: f(x,  # Some comment
        ....:   y,
        ....:   z) = (x +  # Some comment
        ....:         y +
        ....:         z); f
        (x, y, z) |--> x + y + z
    """
    priority = 30

    def __init__(self, del_start, del_end, insert_pos, name, variables):
        """
        INPUT:

        - ``del_start`` -- start of the name of the function

        - ``del_end`` -- end of the ``)``

        - ``insert_pos`` -- position to insert the new command

        - ``name`` -- name of the field or similar

        - ``variables`` -- names of the variables

        TESTS::

            sage: from IPython.core.inputtransformer2 import make_tokens_by_line
            sage: from sage.repl.interpreter import SageCalculusTransformer
            sage: lines = ['f(z,zz) = z + zz\n']
            sage: tokens = make_tokens_by_line(lines)
            sage: S = SageCalculusTransformer.find(tokens)  # indirect doctest
            sage: S.transform(lines)
            ['__tmp__ = var("z,zz"); __tmpf__ = z + zz; f = symbolic_expression(__tmpf__).function(z,zz)\n']
        """
        super().__init__(del_start)
        self.name = name
        self.variables= variables

        # Shift from 1-index to 0-indexed.
        self.del_end = (del_end[0] - 1, del_end[1])
        self.insert_pos = (insert_pos[0] - 1, insert_pos[1])

    @classmethod
    def find(cls, tokens_by_line):
        r"""
        Find the first calculus-like function assignment.

        EXAMPLES::

            sage: from IPython.core.inputtransformer2 import make_tokens_by_line
            sage: from sage.repl.interpreter import SageCalculusTransformer
            sage: lines = ['for i in range(2,20)):\n', '    f(x) = i*x\n']
            sage: tokens = make_tokens_by_line(lines)
            sage: S = SageCalculusTransformer.find(tokens)
            sage: S.transform(lines)
            ['for i in range(2,20)):\n',
             '    __tmp__ = var("x"); __tmpf__ = i*x; f = symbolic_expression(__tmpf__).function(x)\n']
        """
        for line in tokens_by_line:
            name = None
            for i, token in enumerate(line[:-1]):
                if token.string == ')' and line[i+1].string == '=':
                    # Find matching '('
                    variables = []
                    ix = i - 1
                    while ix >= 1:
                        if line[ix].string == '(':
                            break
                        elif line[ix].string == ',':
                            pass
                        elif line[ix].type == tokenize.NAME:
                            variables = [line[ix].string] + variables
                        elif line[ix].type not in {tokenize.NL, tokenize.COMMENT}:
                            ix = 1
                        ix -= 1
                    else:
                        # Incorrect syntax or first token is '('.
                        break

                    if line[ix - 1].type != tokenize.NAME:
                        # A tuple assignment.
                        break

                    name = line[ix - 1].string

                    del_start = line[ix - 1].start
                    del_end = token.end

                    # Find the position to insert the declaration of the generators.
                    ix = i + 1
                    while (not line[ix].string == ';'
                           and not line[ix].type == tokenize.NEWLINE
                           and ix + 1 < len(line)):
                        ix += 1
                    while line[ix - 1].type == tokenize.COMMENT:
                        ix -= 1
                    insert_pos = line[ix].start

                    return cls(del_start, del_end, insert_pos, name, variables)


    def transform(self, lines):
        r"""
        Find the first calculus-like function assignment.

        EXAMPLES::

            sage: from IPython.core.inputtransformer2 import make_tokens_by_line
            sage: from sage.repl.interpreter import SageCalculusTransformer
            sage: lines = ['f(x,y,z) = sqrt(x) - y^z\n']
            sage: tokens = make_tokens_by_line(lines)
            sage: S = SageCalculusTransformer.find(tokens)
            sage: S.transform(lines)
            ['__tmp__ = var("x,y,z"); __tmpf__ = sqrt(x) - y^z; f = symbolic_expression(__tmpf__).function(x,y,z)\n']
        """
        start_line, start_col = self.start_line, self.start_col

        lines = [l for l in lines]

        tmpf = '__tmp{}__'.format(self.name)
        variables = ','.join(self.variables)
        line2 = '; {} = symbolic_expression({}).function({})'.format(self.name, tmpf, variables)

        lines[self.insert_pos[0]] = lines[self.insert_pos[0]][:self.insert_pos[1]] \
                                    + line2 \
                                    + lines[self.insert_pos[0]][self.insert_pos[1]:]

        line1 = lines[start_line][:start_col] + '__tmp__ = var("{}")'.format(variables) + "; " + tmpf
        line1 += lines[self.del_end[0]][self.del_end[1]:]

        lines_before = lines[:start_line]
        lines_after = lines[self.del_end[0] + 1:]

        return lines_before + [line1] + lines_after


SageTokenTransformers = [SageGenConstructionTransformer,
                         SageCalculusTransformer]


###################
# Interface shell #
###################

class InterfaceShellTransformer(PrefilterTransformer):
    priority = 50

    def __init__(self, *args, **kwds):
        """
        Initialize this class.  All of the arguments get passed to
        :meth:`PrefilterTransformer.__init__`.

        .. attribute:: temporary_objects

           a list of hold onto interface objects and keep them from being
           garbage collected

        .. SEEALSO:: :func:`interface_shell_embed`

        TESTS::

            sage: from sage.repl.interpreter import interface_shell_embed
            sage: shell = interface_shell_embed(maxima)
            sage: ift = shell.prefilter_manager.transformers[0]
            sage: ift.temporary_objects
            set()
            sage: ift._sage_import_re.findall('sage(a) + maxima(b)')
            ['sage(', 'maxima(']
        """
        super(InterfaceShellTransformer, self).__init__(*args, **kwds)
        self.temporary_objects = set()
        self._sage_import_re = re.compile(r'(?:sage|%s)\('
                                          % self.shell.interface.name())

    def preparse_imports_from_sage(self, line):
        """
        Finds occurrences of strings such as ``sage(object)`` in
        *line*, converts ``object`` to :attr:`shell.interface`,
        and replaces those strings with their identifier in the new
        system.  This also works with strings such as
        ``maxima(object)`` if :attr:`shell.interface` is
        ``maxima``.

        :param line: the line to transform
        :type line: string

        EXAMPLES::

            sage: from sage.repl.interpreter import interface_shell_embed, InterfaceShellTransformer
            sage: shell = interface_shell_embed(maxima)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.shell.ex('a = 3')
            sage: ift.preparse_imports_from_sage('2 + sage(a)')
            '2 + sage0 '
            sage: maxima.eval('sage0')
            '3'
            sage: ift.preparse_imports_from_sage('2 + maxima(a)') # maxima calls set_seed on startup which is why 'sage0' will becomes 'sage4' and not just 'sage1'
            '2 +  sage4 '
            sage: ift.preparse_imports_from_sage('2 + gap(a)')
            '2 + gap(a)'

        Since :trac:`28439`, this also works with more complicated expressions
        containing nested parentheses::

            sage: shell = interface_shell_embed(gap)
            sage: shell.user_ns = locals()
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: line = '2 + sage((1+2)*gap(-(5-3)^2).sage()) - gap(1+(2-1))'
            sage: line = ift.preparse_imports_from_sage(line)
            sage: gap.eval(line)
            '-12'
        """
        new_line = []
        pos = 0
        while True:
            m = self._sage_import_re.search(line, pos)
            if not m:
                new_line.append(line[pos:])
                break
            expr_start, expr_end = containing_block(line, m.end() - 1,
                                                    delimiters=['()'])
            expr = preparse(line[expr_start + 1:expr_end - 1])
            result = self.shell.interface(eval(expr, self.shell.user_ns))
            self.temporary_objects.add(result)
            new_line += [line[pos:m.start()], result.name()]
            pos = expr_end
        return ' '.join(new_line)

    def transform(self, line, continue_prompt):
        r'''
        Evaluates *line* in :attr:`shell.interface` and returns a
        string representing the result of that evaluation.

        :param line: the line to be transformed *and evaluated*
        :type line: string
        :param continue_prompt: is this line a continuation in a sequence of multiline input?
        :type continue_prompt: bool

        EXAMPLES::

            sage: from sage.repl.interpreter import interface_shell_embed, InterfaceShellTransformer
            sage: shell = interface_shell_embed(maxima)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('2+2', False)   # note: output contains triple quotation marks
            'sage.misc.all.logstr(r"""4""")'
            sage: ift.shell.ex('a = 4')
            sage: ift.transform(r'sage(a)+4', False)
            'sage.misc.all.logstr(r"""8""")'
            sage: ift.temporary_objects
            set()
            sage: shell = interface_shell_embed(gap)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('2+2', False)
            'sage.misc.all.logstr(r"""4""")'

        TESTS:

        Check that whitespace is not stripped and that special characters are
        escaped (:trac:`28439`)::

            sage: shell = interface_shell_embed(gap)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.transform(r'Print("  -\n\\\\-  ");', False)
            'sage.misc.all.logstr(r"""  -\n\\\\-""")'

            sage: shell = interface_shell_embed(macaulay2)  # optional - macaulay2
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager) # optional - macaulay2
            sage: ift.transform('net(ZZ^2)', False)  # optional - macaulay2
            'sage.misc.all.logstr(r"""  2\nZZ""")'
        '''
        line = self.preparse_imports_from_sage(line)

        try:
            t = self.shell.interface.eval(line)
        finally:
            # Once we've evaluated the lines, we can clear the
            # temporary objects
            self.temporary_objects = set()
        # We do not strip whitespace from t here as the individual interface is
        # responsible for that
        return 'sage.misc.all.logstr(r"""%s""")' % t


def interface_shell_embed(interface):
    """
    Returns an IPython shell which uses a Sage interface on the
    backend to perform the evaluations.  It uses
    :class:`InterfaceShellTransformer` to transform the input into the
    appropriate ``interface.eval(...)`` input.

    INPUT:

    - ``interface`` -- A Sage ``PExpect`` interface instance.

    EXAMPLES::

        sage: from sage.repl.interpreter import interface_shell_embed
        sage: shell = interface_shell_embed(gap)
        sage: shell.run_cell('List( [1..10], IsPrime )')
        [ false, true, true, false, true, false, true, false, false, false ]
        <ExecutionResult object at ..., execution_count=None error_before_exec=None error_in_exec=None ...result=[ false, true, true, false, true, false, true, false, false, false ]>
    """
    cfg = sage_ipython_config.copy()
    ipshell = InteractiveShellEmbed(config=cfg,
                                    banner1='\n  --> Switching to %s <--\n\n' % interface,
                                    exit_msg='\n  --> Exiting back to Sage <--\n')
    ipshell.interface = interface
    ipshell.prompts = InterfacePrompts(interface.name())

    while ipshell.prefilter_manager.transformers:
        ipshell.prefilter_manager.transformers.pop()
    while ipshell.prefilter_manager.checkers:
        ipshell.prefilter_manager.checkers.pop()
    ipshell.ex('import sage.misc.all')

    InterfaceShellTransformer(shell=ipshell,
                              prefilter_manager=ipshell.prefilter_manager,
                              config=cfg)
    return ipshell


def get_test_shell():
    """
    Return a IPython shell that can be used in testing the functions
    in this module.

    OUTPUT:

    An IPython shell

    EXAMPLES::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell(); shell
        <sage.repl.interpreter.SageTestShell object at 0x...>
        sage: shell.parent.shell_class
        <class 'sage.repl.interpreter.SageTestShell'>
        sage: shell.parent.test_shell
        True
        sage: shell.quit()

    TESTS:

    Check that :trac:`14070` has been resolved::

        sage: from sage.tests.cmdline import test_executable
        sage: cmd = 'from sage.repl.interpreter import get_test_shell; shell = get_test_shell()'
        sage: (out, err, ret) = test_executable(["sage", "-c", cmd])
        sage: out + err
        ''
    """
    config = sage_ipython_config.default()
    config.TerminalIPythonApp.test_shell = True
    config.TerminalIPythonApp.shell_class = SageTestShell
    app = SageTerminalApp.instance(config=config)
    if app.shell is None:
        app.initialize(argv=[])
    else:
        try:
            app.shell._restart()
        except AttributeError:
            pass
    # No quit noise
    app.shell.verbose_quit = False
    return app.shell


#######################
# IPython TerminalApp #
#######################

class SageCrashHandler(IPAppCrashHandler):
    def __init__(self, app):
        """
        A custom :class:`CrashHandler` which gives the user
        instructions on how to post the problem to sage-support.

        EXAMPLES::

            sage: from sage.repl.interpreter import SageTerminalApp, SageCrashHandler
            sage: app = SageTerminalApp.instance()
            sage: sch = SageCrashHandler(app); sch
            <sage.repl.interpreter.SageCrashHandler object at 0x...>
            sage: sorted(sch.info.items())
            [('app_name', u'Sage'),
             ('bug_tracker', 'http://trac.sagemath.org'),
             ('contact_email', 'sage-support@googlegroups.com'),
             ('contact_name', 'sage-support'),
             ('crash_report_fname', u'Crash_report_Sage.txt')]
        """
        contact_name = 'sage-support'
        contact_email = 'sage-support@googlegroups.com'
        bug_tracker = 'http://trac.sagemath.org'
        CrashHandler.__init__(self,
            app, contact_name, contact_email, bug_tracker, show_crash_traceback=True)
        self.crash_report_fname = 'Sage_crash_report.txt'


class SageTerminalApp(TerminalIPythonApp):
    name = u'Sage'
    crash_handler_class = SageCrashHandler

    test_shell = Bool(False, help='Whether the shell is a test shell')
    test_shell.tag(config=True)
    shell_class = Type(InteractiveShell, help='Type of the shell')
    shell_class.tag(config=True)

    def load_config_file(self, *args, **kwds):
        r"""
        Merges a config file with the default sage config.

        .. note::

            This code is based on :meth:`Application.update_config`.

        TESTS:

        Test that :trac:`15972` has been fixed::

            sage: from sage.misc.temporary_file import tmp_dir
            sage: from sage.repl.interpreter import SageTerminalApp
            sage: d = tmp_dir()
            sage: from IPython.paths import get_ipython_dir
            sage: IPYTHONDIR = get_ipython_dir()
            sage: os.environ['IPYTHONDIR'] = d
            sage: SageTerminalApp().load_config_file()
            sage: os.environ['IPYTHONDIR'] = IPYTHONDIR
        """
        super(SageTerminalApp, self).load_config_file(*args, **kwds)
        newconfig = sage_ipython_config.default()
        # merge in the config loaded from file
        newconfig.merge(self.config)
        self.config = newconfig

    def init_shell(self):
        r"""
        Initialize the :class:`SageInteractiveShell` instance.

        .. note::

            This code is based on
            :meth:`TerminalIPythonApp.init_shell`.

        EXAMPLES::

            sage: from sage.repl.interpreter import SageTerminalApp
            sage: app = SageTerminalApp.instance()
            sage: app.shell
            <sage.repl.interpreter.SageTestShell object at 0x...>
        """
        # Shell initialization
        self.shell = self.shell_class.instance(
            parent=self,
            config=self.config,
            profile_dir=self.profile_dir,
            ipython_dir=self.ipython_dir)
        self.shell.configurables.append(self)
        self.shell.has_sage_extensions = SAGE_EXTENSION in self.extensions

        # Load the %lprun extension if available
        try:
            import line_profiler
            assert line_profiler  # silence pyflakes
        except ImportError:
            pass
        else:
            self.extensions.append('line_profiler')

        if self.shell.has_sage_extensions:
            self.extensions.remove(SAGE_EXTENSION)

            # load sage extension here to get a crash if
            # something is wrong with the sage library
            self.shell.extension_manager.load_extension(SAGE_EXTENSION)
