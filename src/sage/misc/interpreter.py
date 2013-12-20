r"""
Sage's IPython Modifications

This module contains all of Sage's customizations to the IPython
interpreter.  These changes consist of the following major components:

  - :class:`SageTerminalApp`
  - :class:`SageInteractiveShell`
  - :func:`interface_shell_embed`

SageTerminalApp
---------------

This is the main application object.  It is used by the
``$SAGE_LOCAL/bin/sage-ipython`` script to start the Sage
command-line.  It's primary purpose is to

  - Initialize the :class:`SageInteractiveShell`.

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

  - Cleanly deinitialize the Sage library before exiting.  See
    :meth:`~SageInteractiveShell.ask_exit`.

  - Modify the libraries before calling system commands. See
    :meth:`~SageInteractiveShell.system_raw`.

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
import os, log, re, new, sys
from IPython.utils.py3compat import cast_unicode
from IPython.utils.traitlets import (Integer, CBool, CaselessStrEnum, Enum,
                                     List, Unicode, Instance, Type)
from preparser import (preparse, preparse_file, load_wrap)

def embedded():
    """
    Returns True if Sage is being run from the notebook.

    EXAMPLES::

        sage: from sage.misc.interpreter import embedded
        sage: embedded()
        False
    """
    import sage.server.support
    return sage.server.support.EMBEDDED_MODE


#TODO: This global variable do_preparse should be associtated with an
#IPython InteractiveShell as opposed to a global variable in this
#module.
do_preparse=True
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
    global do_preparse
    if on:
        do_preparse = True
    else:
        do_preparse = False


###############################################################
# Old code for handling the sage prompt in previous verisons of
# IPython
###############################################################
def set_sage_prompt(s):
    """
    Sets the Sage prompt to the string ``s``.

    :param s: the new prompt
    :type s: string
    :returns: None

    EXAMPLES::

        sage: from sage.misc.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell('from sage.misc.interpreter import set_sage_prompt')
        sage: shell.run_cell('set_sage_prompt(u"new")')
        sage: shell.prompt_manager.in_template
        u'new: '
        sage: shell.run_cell('set_sage_prompt(u"sage")')
    """
    ipython = get_ipython()
    ipython.prompt_manager.in_template = s+': '

def sage_prompt():
    """
    Returns the current Sage prompt.

    EXAMPLES::

        sage: from sage.misc.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell('sage_prompt()')
        u'sage'
    """
    ipython = get_ipython()
    return ipython.prompt_manager.in_template[:-2]

####################
# InteractiveShell #
####################
from IPython.frontend.terminal.interactiveshell import TerminalInteractiveShell
class SageInteractiveShell(TerminalInteractiveShell):

    def system_raw(self, cmd):
        """
        Run a system command.

        If the command is not a sage-specific binary, adjust the library paths before calling
        system commands.  See Trac #975 for a discussion of running system commands.

        This is equivalent to the sage-native-execute shell script.

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell
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
        path = os.path.join(os.environ['SAGE_LOCAL'],'bin',
                            re.split('[ |\n\t;&]', cmd)[0])
        if not os.access(path, os.X_OK):
            libraries = 'LD_LIBRARY_PATH="$SAGE_ORIG_LD_LIBRARY_PATH";export LD_LIBRARY_PATH;'
            if os.uname()[0]=='Darwin':
                libraries += 'DYLD_LIBRARY_PATH="$SAGE_ORIG_DYLD_LIBRARY_PATH";export DYLD_LIBRARY_PATH;'
            cmd = libraries+cmd
        return super(SageInteractiveShell, self).system_raw(cmd)

    def ask_exit(self):
        """
        We need to run the :func:`sage.all.quit_sage` function when we
        exit the shell.

        EXAMPLES:

        We install a fake :func:`sage.all.quit_sage`::

            sage: import sage.all
            sage: old_quit = sage.all.quit_sage
            sage: def new_quit(): print "Quitter!!!"
            sage: sage.all.quit_sage = new_quit

        Now, we can check to see that this method works::

            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.ask_exit()
            Quitter!!!
            sage: shell.exit_now
            True

        Clean up after ourselves::

            sage: shell.exit_now = False
            sage: sage.all.quit_sage = old_quit
        """
        from sage.all import quit_sage
        quit_sage()
        super(SageInteractiveShell, self).ask_exit()

###################################################################
# Transformers used in the SageInputSplitter
###################################################################
# These used to be part of the older PrefilterTransformer framework,
# but that is not the modern way of doing things. For more details, see
# http://mail.scipy.org/pipermail/ipython-dev/2011-March/007334.html

class SagePreparseTransformer():
    """
    Preparse the line of code before it get evaluated by IPython.
    """
    def __call__(self, line, line_number):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import SagePreparseTransformer
            sage: spt = SagePreparseTransformer()
            sage: spt('2', 0)
            'Integer(2)'
            sage: preparser(False)
            sage: spt('2', 0)
            '2'
            sage: preparser(True)

        TESTS:

        Check that syntax errors in the preparser do not crash IPython,
        see :trac:`14961`. ::

            sage: bad_syntax = "R.<t> = QQ{]"
            sage: preparse(bad_syntax)
            Traceback (most recent call last):
            ...
            SyntaxError: Mismatched ']'
            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell(bad_syntax)
            SyntaxError: Mismatched ']'
        """
        if do_preparse and not line.startswith('%'):
            # we use preparse_file instead of just preparse because preparse_file
            # automatically prepends attached files
            try:
                return preparse(line, reset=(line_number==0))
            except SyntaxError as err:
                print "SyntaxError: {0}".format(err)
            return ''
        else:
            return line

class MagicTransformer():
    r"""
    Handle input lines that start out like ``load ...`` or ``attach
    ...``.

    Since there are objects in the Sage namespace named ``load`` and
    ``attach``, IPython's automagic will not transform these lines
    into ``%load ...`` and ``%attach ...``, respectively.  Thus, we
    have to do it manually.
    """
    deprecations = {'load': '%runfile',
                    'attach': '%attach',
                    'time': '%time'}
    def __call__(self, line, line_number):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell, MagicTransformer
            sage: mt = MagicTransformer()
            sage: mt('load /path/to/file', 0)
            doctest:1: DeprecationWarning: Use %runfile instead of load.
            See http://trac.sagemath.org/12719 for details.
            '%runfile /path/to/file'
            sage: mt('attach /path/to/file', 0)
            doctest:1: DeprecationWarning: Use %attach instead of attach.
            See http://trac.sagemath.org/12719 for details.
            '%attach /path/to/file'
            sage: mt('time 1+2', 0)
            doctest:1: DeprecationWarning: Use %time instead of time.
            See http://trac.sagemath.org/12719 for details.
            '%time 1+2'
        """
        for old,new in self.deprecations.items():
            if line.startswith(old+' '):
                from sage.misc.superseded import deprecation
                deprecation(12719, 'Use %s instead of %s.'%(new,old))
                return new+line[len(old):]
        return line

class SagePromptTransformer():
    """
    Remove Sage prompts from the input line.
    """
    _sage_prompt_re = re.compile(r'(^[ \t]*sage: |^[ \t]*\.+:? )+')

    def __call__(self, line, line_number):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import SagePromptTransformer
            sage: spt = SagePromptTransformer()
            sage: spt("sage: sage: 2 + 2", 0)
            '2 + 2'
            sage: spt('', 0)
            ''
            sage: spt("      sage: 2+2", 0)
            '2+2'
            sage: spt("      ... 2+2", 0)
            '2+2'
        """
        if not line or line.isspace() or line.strip() == '...':
            # This allows us to recognize multiple input prompts separated by
            # blank lines and pasted in a single chunk, very common when
            # pasting doctests or long tutorial passages.
            return ''
        while True:
            m = self._sage_prompt_re.match(line)
            if m:
                line = line[len(m.group(0)):]
            else:
                break
        return line

class SagePromptDedenter():
    """
    Remove leading spaces from the input line.
    """
    def __call__(self, line, line_number):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        - ``line_number`` -- integer. The line number. For a single-line input, this is always zero.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import SagePromptDedenter
            sage: spd = SagePromptDedenter()
            sage: spd('  1 + \\', 0)
            '1 + \\'
            sage: spd('  2', 1)
            '2'
            sage: spd('3', 2)   # try our best with incorrect indentation
            '3'
        """
        if line_number == 0:
            dedent_line = line.lstrip()
            self._dedent = len(line) - len(dedent_line)
            return dedent_line
        else:
            dedent = min(len(line)-len(line.lstrip()), self._dedent)
            return line[dedent:]

###################
# Interface shell #
###################
from IPython.core.prefilter import PrefilterTransformer
from IPython.frontend.terminal.embed import InteractiveShellEmbed
from IPython import Config

class InterfaceShellTransformer(PrefilterTransformer):
    priority = 50
    def __init__(self, *args, **kwds):
        """
        Initialize this class.  All of the arguments get passed to
        :meth:`PrefilterTransformer.__init__`.

        .. attribute::  lines_queue

            a list of lines to be evaluated

        .. attribute:: temporary_objects

           a list of hold onto interface objects and keep them from being
           garbage collected

        .. seealso:: :func:`interface_shell_embed`

        EXAMPLES::

            sage: from sage.misc.interpreter import interface_shell_embed
            sage: shell = interface_shell_embed(maxima)
            sage: ift = shell.prefilter_manager.transformers[0]
            sage: ift.lines_queue
            []
            sage: ift.temporary_objects
            []
            sage: ift._sage_import_re.findall('sage(a) + maxima(b)')
            ['a', 'b']
        """
        super(InterfaceShellTransformer, self).__init__(*args, **kwds)
        self.lines_queue = []
        self.temporary_objects = []
        self._sage_import_re = re.compile(r'(?:sage|%s)\((.*?)\)'%self.shell.interface.name())

    def preparse_imports_from_sage(self, line):
        """
        Finds occurrences of strings such as ``sage(object)`` in
        *line*, converts ``object`` to :attr:`shell.interface`,
        and replaces those strings with their identifier in the new
        system.  This also works with strings such as
        ``maxima(object`` if :attr:`shell.interface` is
        ``maxima``.

        :param line: the line to transform
        :type line: string

        .. warning::

            This does not parse nested parentheses correctly.  Thus,
            lines like ``sage(a.foo())`` will not work correctly.
            This can't be done in generality with regular expressions.

        EXAMPLES::

            sage: from sage.misc.interpreter import interface_shell_embed, InterfaceShellTransformer
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
        from sage_eval import sage_eval
        for sage_code in self._sage_import_re.findall(line):
            expr = preparse(sage_code)
            result = self.shell.interface(sage_eval(expr, self.shell.user_ns))
            self.temporary_objects.append(result)
            line = self._sage_import_re.sub(' ' + result.name() + ' ', line, 1)
        return line

    def transform(self, line, continue_prompt):
        """
        Evaluates *line* in :attr:`shell.interface` and returns a
        string representing the result of that evaluation.  If a line
        ends in backspace, then this method will store *line* in
        :attr:`lines_queue` until it receives a line not ending in
        backspace.

        :param line: the line to be transformed *and evaluated*
        :type line: string
        :param continue_prompt: is this line a continuation in a sequence of multiline input?
        :type continue_prompt: bool

        EXAMPLES::

            sage: from sage.misc.interpreter import interface_shell_embed, InterfaceShellTransformer
            sage: shell = interface_shell_embed(maxima)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('2+2', False)   # note: output contains triple quotation marks
            'sage.misc.all.logstr(...4...)'
            sage: ift.shell.ex('a = 4')
            sage: ift.transform(r'sage(a)+\\', False)
            sage: ift.shell.prompt_manager.in_template
            u'....: '
            sage: ift.lines_queue
            [' sage2 +']
            sage: ift.temporary_objects
            [4]
            sage: ift.transform('4', False)
            'sage.misc.all.logstr(...8...)'
            sage: ift.lines_queue
            []
            sage: ift.temporary_objects
            []

            sage: shell = interface_shell_embed(gap)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config, prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('2+2', False)
            'sage.misc.all.logstr(...4...)'
        """
        line = self.preparse_imports_from_sage(line)
        line = line.rstrip()
        if line.endswith('\\'):
            line = line.strip('\\')
            self.lines_queue.append(line)
            #TODO: This should beter handled by getting IPython to
            #switch to a continuation prompt
            self.shell.prompt_manager.in_template = self.shell.prompt_manager.in2_template
        else:
            self.lines_queue.append(line)
            line = "".join(self.lines_queue)

            if self.shell.interface.name() in ['gap', 'magma', 'kash', 'singular']:
                if not line.endswith(';'):
                    line += ';'
            elif self.shell.interface.name() == 'mathematica':
                line = 'InputForm[%s]'%line

            t = self.shell.interface.eval(line)

            #Once we've evaluated the lines, we can clear the queue
            #and temporary objects and switch the prompt back
            self.lines_queue = []
            self.temporary_objects = []
            self.shell.prompt_manager.in_template = self.shell.interface.name() + ': '
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

        sage: from sage.misc.interpreter import interface_shell_embed
        sage: shell = interface_shell_embed(gap)
        sage: shell.run_cell('List( [1..10], IsPrime )')
        [ false, true, true, false, true, false, true, false, false, false ]
    """
    try:
        cfg = Config(get_ipython().config)
    except NameError:
        cfg = Config(DEFAULT_SAGE_CONFIG)
    cfg.PromptManager['in_template'] = interface.name() + ': '

    ipshell = InteractiveShellEmbed(config=cfg,
                                    banner1='\n  --> Switching to %s <--\n\n'%interface,
                                    exit_msg = '\n  --> Exiting back to Sage <--\n')
    ipshell.interface = interface

    while ipshell.prefilter_manager.transformers:
        ipshell.prefilter_manager.transformers.pop()
    while ipshell.prefilter_manager.checkers:
        ipshell.prefilter_manager.checkers.pop()
    ipshell.ex('from sage.all import *')

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

        sage: from sage.misc.interpreter import get_test_shell
        sage: shell = get_test_shell(); shell
        <sage.misc.interpreter.SageInteractiveShell object at 0x...>

    TESTS:

    Check that :trac:`14070` has been resolved::

        sage: from sage.tests.cmdline import test_executable
        sage: cmd = 'from sage.misc.interpreter import get_test_shell; shell = get_test_shell()'
        sage: (out, err, ret) = test_executable(["sage", "-c", cmd])
        sage: out + err
        ''
    """
    app = SageTerminalApp.instance(config=DEFAULT_SAGE_CONFIG)
    if app.shell is None:
        app.initialize(argv=[])
    return app.shell


#######################
# IPython TerminalApp #
#######################
from IPython.frontend.terminal.ipapp import TerminalIPythonApp, IPAppCrashHandler
from IPython.core.crashhandler import CrashHandler
from IPython import Config

class SageCrashHandler(IPAppCrashHandler):
    def __init__(self, app):
        """
        A custom :class:`CrashHandler` which gives the user
        instructions on how to post the problem to sage-support.

        EXAMPLES::

            sage: from sage.misc.interpreter import SageTerminalApp, SageCrashHandler
            sage: app = SageTerminalApp.instance()
            sage: sch = SageCrashHandler(app); sch
            <sage.misc.interpreter.SageCrashHandler object at 0x...>
            sage: sorted(sch.info.items())
            [('app_name', u'Sage'),
             ('bug_tracker', 'http://trac.sagemath.org/sage_trac'),
             ('contact_email', 'sage-support@googlegroups.com'),
             ('contact_name', 'sage-support'),
             ('crash_report_fname', u'Crash_report_Sage.txt')]
        """
        contact_name = 'sage-support'
        contact_email = 'sage-support@googlegroups.com'
        bug_tracker = 'http://trac.sagemath.org/sage_trac'
        CrashHandler.__init__(self,
            app, contact_name, contact_email, bug_tracker, show_crash_traceback=False)
        self.crash_report_fname = 'Sage_crash_report.txt'

DEFAULT_SAGE_CONFIG = Config(
    PromptManager = Config(
        in_template = 'sage: ',
        in2_template = '....: ',
        justify = False,
        out_template = ''),
    TerminalIPythonApp = Config(
        display_banner = False,
        verbose_crash = True),
    TerminalInteractiveShell = Config(
        ast_node_interactivity = 'all',
        colors = 'LightBG' if sys.stdout.isatty() else 'NoColor',
        confirm_exit = False,
        separate_in = ''),
    # The extension is *always* loaded for SageTerminalApp
    # See the code for SageTerminalApp.init_shell
    #InteractiveShellApp = Config(extensions=['sage.misc.sage_extension']),
    )

class SageTerminalApp(TerminalIPythonApp):
    name = u'Sage'
    crash_handler_class = SageCrashHandler
    test_shell = False

    def __init__(self, **kwargs):
        self.command_line_config = kwargs.get('config', Config())
        super(SageTerminalApp, self).__init__(**kwargs)


    def load_config_file(self, *args, **kwds):
        from IPython.frontend.terminal.ipapp import default_config_file_name
        from IPython.config.loader import PyFileConfigLoader, ConfigFileNotFound
        from IPython.core.profiledir import ProfileDir
        from IPython.utils.path import get_ipython_dir

        conf = Config()
        conf._merge(DEFAULT_SAGE_CONFIG)
        conf._merge(self.command_line_config)

        # Get user config.
        sage_profile_dir = ProfileDir.find_profile_dir_by_name(
            get_ipython_dir(), 'sage').location
        try:
            cl = PyFileConfigLoader(default_config_file_name, sage_profile_dir)
            conf._merge(cl.load_config())
        except ConfigFileNotFound:
            pass
        self.update_config(conf)


    def init_shell(self):
        r"""
        Initialize the :class:`SageInteractiveShell` instance.
        Additionally, this also does the following:

          - Merges the default shell configuration with the user's.

        .. note::

            This code is based on
            :meth:`TermintalIPythonApp.init_shell`.

        EXAMPLES::

            sage: from sage.misc.interpreter import SageTerminalApp, DEFAULT_SAGE_CONFIG
            sage: app = SageTerminalApp(config=DEFAULT_SAGE_CONFIG)
            sage: app.initialize(argv=[])  # indirect doctest
            sage: app.shell
            <sage.misc.interpreter.SageInteractiveShell object at 0x...>
        """
        # We need verbose crashes for the Sage crash handler.  We set it here
        # so that we don't overwrite the traitlet attribute
        self.verbose_crash = True

        # Shell initialization
        self.shell = SageInteractiveShell.instance(config=self.config,
                        display_banner=False, profile_dir=self.profile_dir,
                        ipython_dir=self.ipython_dir)
        self.shell.configurables.append(self)
        self.shell.extension_manager.load_extension('sage.misc.sage_extension')
