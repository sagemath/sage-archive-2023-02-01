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

"""

#*****************************************************************************
#       Copyright (C) 2004-2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import copy
import os
import re
import sys
from sage.misc.preparser import preparse

from IPython import Config

from sage.env import SAGE_LOCAL

SAGE_EXTENSION = 'sage'

DEFAULT_SAGE_CONFIG = Config(
    PromptManager = Config(
        in_template = 'sage: ',
        in2_template = '....: ',
        justify = False,
        out_template = ''),
    TerminalIPythonApp = Config(
        display_banner = False,
        verbose_crash = True),
    InteractiveShell = Config(
        ast_node_interactivity = 'all',
        colors = 'LightBG' if sys.stdout.isatty() else 'NoColor',
        confirm_exit = False,
        separate_in = ''),
    InteractiveShellApp = Config(extensions=[SAGE_EXTENSION]),
    )

def embedded():
    """
    Returns True if Sage is being run from the notebook.

    EXAMPLES::

        sage: from sage.repl.interpreter import embedded
        sage: embedded()
        False
    """
    import sage.server.support
    return sage.server.support.EMBEDDED_MODE

#TODO: This global variable do_preparse should be associtated with an
#IPython InteractiveShell as opposed to a global variable in this
#module.
_do_preparse=True
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
# (Terminal)InteractiveShell #
##############################
class SageNativeExecute(object):

    def system_raw(self, cmd):
        """
        Run a system command.

        If the command is not a sage-specific binary, adjust the library
        paths before calling system commands.  See :trac:`975` for a
        discussion of running system commands.

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
            sage: shell.system_raw('env | grep "^LD_LIBRARY_PATH=" | grep $SAGE_LOCAL')
            sage: shell.user_ns['_exit_code']
            1
            sage: shell.system_raw('R --version')
            R version ...
            sage: shell.user_ns['_exit_code']
            0
        """
        path = os.path.join(SAGE_LOCAL, 'bin',
                            re.split(r'[\s|;&]', cmd)[0])
        if not os.access(path, os.X_OK):
            libraries = 'LD_LIBRARY_PATH="$SAGE_ORIG_LD_LIBRARY_PATH";export LD_LIBRARY_PATH;'
            if os.uname()[0]=='Darwin':
                libraries += 'DYLD_LIBRARY_PATH="$SAGE_ORIG_DYLD_LIBRARY_PATH";export DYLD_LIBRARY_PATH;'
            cmd = libraries+cmd
        return super(SageNativeExecute, self).system_raw(cmd)


from IPython.core.interactiveshell import InteractiveShell
from IPython.terminal.interactiveshell import TerminalInteractiveShell

class SageInteractiveShell(SageNativeExecute, InteractiveShell):
    pass

class SageTerminalInteractiveShell(SageNativeExecute, TerminalInteractiveShell):
    pass


###################################################################
# Transformers used in the SageInputSplitter
###################################################################
from IPython.core.inputtransformer import (CoroutineInputTransformer,
                                           StatelessInputTransformer,
                                           _strip_prompts)

@StatelessInputTransformer.wrap
def SagePreparseTransformer(line):
    r"""
    EXAMPLES::

        sage: from sage.repl.interpreter import SagePreparseTransformer
        sage: spt = SagePreparseTransformer()
        sage: spt.push('1+1r+2.3^2.3r')
        "Integer(1)+1+RealNumber('2.3')**2.3"
        sage: preparser(False)
        sage: spt.push('2.3^2')
        '2.3^2'

    TESTS:

    Check that syntax errors in the preparser do not crash IPython,
    see :trac:`14961`. ::

        sage: preparser(True)
        sage: bad_syntax = "R.<t> = QQ{]"
        sage: preparse(bad_syntax)
        Traceback (most recent call last):
        ...
        SyntaxError: Mismatched ']'
        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell(bad_syntax)
          File "<string>", line unknown
        SyntaxError: Mismatched ']'
        <BLANKLINE>
    """
    if _do_preparse and not line.startswith('%'):
        return preparse(line)
    else:
        return line

@CoroutineInputTransformer.wrap
def SagePromptTransformer():
    r"""
    Strip the sage:/....: prompts of Sage.

    EXAMPLES::

        sage: from sage.repl.interpreter import SagePromptTransformer
        sage: spt = SagePromptTransformer()
        sage: spt.push("sage: 2 + 2")
        '2 + 2'
        sage: spt.push('')
        ''
        sage: spt.push("....: 2+2")
        '2+2'

    This should strip multiple prompts: see :trac:`16297`::

        sage: spt.push("sage:   sage: 2+2")
        '2+2'
        sage: spt.push("   sage: ....: 2+2")
        '2+2'

    The prompt contains a trailing space. Extra spaces between the
    last prompt and the remainder should not be stripped::

        sage: spt.push("   sage: ....:    2+2")
        '   2+2'

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
    """
    _sage_prompt_re = re.compile(r'^(\s*(:?sage: |\.\.\.\.: ))+')
    return _strip_prompts(_sage_prompt_re)

###################
# Interface shell #
###################
from IPython.core.prefilter import PrefilterTransformer
from IPython.terminal.embed import InteractiveShellEmbed

class InterfaceShellTransformer(PrefilterTransformer):
    priority = 50
    def __init__(self, *args, **kwds):
        """
        Initialize this class.  All of the arguments get passed to
        :meth:`PrefilterTransformer.__init__`.

        .. attribute:: temporary_objects

           a list of hold onto interface objects and keep them from being
           garbage collected

        .. seealso:: :func:`interface_shell_embed`

        EXAMPLES::

            sage: from sage.repl.interpreter import interface_shell_embed
            sage: shell = interface_shell_embed(maxima)
            sage: ift = shell.prefilter_manager.transformers[0]
            sage: ift.temporary_objects
            set()
            sage: ift._sage_import_re.findall('sage(a) + maxima(b)')
            ['a', 'b']
        """
        super(InterfaceShellTransformer, self).__init__(*args, **kwds)
        self.temporary_objects = set()
        self._sage_import_re = re.compile(r'(?:sage|%s)\((.*?)\)'%self.shell.interface.name())

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

        .. warning::

            This does not parse nested parentheses correctly.  Thus,
            lines like ``sage(a.foo())`` will not work correctly.
            This can't be done in generality with regular expressions.

        EXAMPLES::

            sage: from sage.repl.interpreter import interface_shell_embed, InterfaceShellTransformer
            sage: shell = interface_shell_embed(maxima)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.shell.ex('a = 3')
            sage: ift.preparse_imports_from_sage('2 + sage(a)')
            '2 + sage0 '
            sage: maxima.eval('sage0')
            '3'
            sage: ift.preparse_imports_from_sage('2 + maxima(a)')
            '2 +  sage1 '
            sage: ift.preparse_imports_from_sage('2 + gap(a)')
            '2 + gap(a)'
        """
        for sage_code in self._sage_import_re.findall(line):
            expr = preparse(sage_code)
            result = self.shell.interface(eval(expr, self.shell.user_ns))
            self.temporary_objects.add(result)
            line = self._sage_import_re.sub(' ' + result.name() + ' ', line, 1)
        return line

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
            'sage.misc.all.logstr("""4""")'
            sage: ift.shell.ex('a = 4')
            sage: ift.transform(r'sage(a)+4', False)
            'sage.misc.all.logstr("""8""")'
            sage: ift.temporary_objects
            set()
            sage: shell = interface_shell_embed(gap)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('2+2', False)
            'sage.misc.all.logstr("""4""")'
        '''
        line = self.preparse_imports_from_sage(line)

        try:
            t = self.shell.interface.eval(line)
        finally:
            # Once we've evaluated the lines, we can clear the
            # temporary objects
            self.temporary_objects = set()
        return 'sage.misc.all.logstr("""%s""")'%t.strip()

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
    """
    try:
        cfg = copy.deepcopy(get_ipython().config)
    except NameError:
        cfg = copy.deepcopy(DEFAULT_SAGE_CONFIG)
    cfg.PromptManager['in_template'] = interface.name() + ': '
    cfg.PromptManager['in2_template'] = len(interface.name())*'.' + ': '

    ipshell = InteractiveShellEmbed(config=cfg,
                                    banner1='\n  --> Switching to %s <--\n\n'%interface,
                                    exit_msg = '\n  --> Exiting back to Sage <--\n')
    ipshell.interface = interface

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
    Returns a IPython shell that can be used in testing the functions
    in this module.

    :returns: an IPython shell

    EXAMPLES::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell(); shell
        <sage.repl.interpreter.SageTerminalInteractiveShell object at 0x...>

    TESTS:

    Check that :trac:`14070` has been resolved::

        sage: from sage.tests.cmdline import test_executable
        sage: cmd = 'from sage.repl.interpreter import get_test_shell; shell = get_test_shell()'
        sage: (out, err, ret) = test_executable(["sage", "-c", cmd])
        sage: out + err
        ''
    """
    app = SageTerminalApp.instance(config=copy.deepcopy(DEFAULT_SAGE_CONFIG))
    if app.shell is None:
        app.initialize(argv=[])
    # overwrite the default (console + graphics) formatter with the plain text one
    import sage.repl.display.formatter as formatter
    app.shell.display_formatter.formatters['text/plain'] = (
        formatter.SagePlainTextFormatter(config=app.shell.config))
    # No quit noise
    app.shell.verbose_quit = False
    return app.shell

#######################
# IPython TerminalApp #
#######################
from IPython.terminal.ipapp import TerminalIPythonApp, IPAppCrashHandler
from IPython.core.crashhandler import CrashHandler

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
            app, contact_name, contact_email, bug_tracker, show_crash_traceback=False)
        self.crash_report_fname = 'Sage_crash_report.txt'


class SageTerminalApp(TerminalIPythonApp):
    name = u'Sage'
    crash_handler_class = SageCrashHandler
    test_shell = False

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
            sage: IPYTHONDIR = os.environ['IPYTHONDIR']
            sage: os.environ['IPYTHONDIR'] = d
            sage: SageTerminalApp().load_config_file()
            sage: os.environ['IPYTHONDIR'] = IPYTHONDIR
        """
        super(SageTerminalApp, self).load_config_file(*args, **kwds)

        newconfig = copy.deepcopy(DEFAULT_SAGE_CONFIG)

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

            sage: from sage.repl.interpreter import SageTerminalApp, DEFAULT_SAGE_CONFIG
            sage: app = SageTerminalApp(config=DEFAULT_SAGE_CONFIG)
            sage: app.initialize(argv=[])  # indirect doctest
            sage: app.shell
            <sage.repl.interpreter.SageTerminalInteractiveShell object at 0x...>
        """
        # Shell initialization
        self.shell = SageTerminalInteractiveShell.instance(config=self.config,
                        display_banner=False, profile_dir=self.profile_dir,
                        ipython_dir=self.ipython_dir)
        self.shell.configurables.append(self)
        self.shell.has_sage_extensions = SAGE_EXTENSION in self.extensions

        if self.shell.has_sage_extensions:
            self.extensions.remove(SAGE_EXTENSION)

            # load sage extension here to get a crash if
            # something is wrong with the sage library
            self.shell.extension_manager.load_extension(SAGE_EXTENSION)


