r"""
Sage's IPython Modifications

This module contains all of Sage's customizations to the IPython
interpreter.  These changes consist of the following magjor components:

  - :class:`SageTerminalApp`
  - :class:`SageInteractiveShell`
  - :func:`interface_shell_embed`

SageTerminalApp
---------------

This is the main application object.  It is used by the
``$SAGE_ROOT/local/bin/sage-ipython`` script to start the Sage
command-line.  It's primary purpose is to

  - Initialize the :class:`SageInteractiveShell`.

  - Provide default configuration options for the shell, and its
    subcomponents.  These work with (and can be overrided by)
    IPython's configuration system.

  - Monkey-patch IPython in order to support Sage's customization when
    introspecting objects.

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

  - Modifying the input before it is evaluated. This is done in
    :class:`SageInputSplitter`, which uses
    :class:`SagePromptTransformer`, :class:`SagePreparseTransformer`,
    :class:`LoadAttachTransformer`, and
    :class:`InterfaceMagicTransformer` to do the actual work.

  - Provide a number of IPython magic functions that work with Sage
    and its preparser.  See :meth:`~SageInteractiveShell.magic_timeit`,
    :meth:`~SageInteractiveShell.magic_prun`,
    :meth:`~SageInteractiveShell.magic_load`,
    :meth:`~SageInteractiveShell.magic_attach`, and
    :meth:`~SageInteractiveShell.magic_iload`.

  - Adding support for attached files.  See
    :meth:`~SageInteractiveShell.run_cell`.

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
from preparser import (preparse, preparse_file, load_wrap,
                       modified_attached_files, attached_files)

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
        sage: shell.run_cell('set_sage_prompt(u"new: ")')
        sage: shell.prompt_manager.in_template
        u'new: '
        sage: shell.run_cell('set_sage_prompt(u"sage: ")')
    """
    ipython = get_ipython()
    ipython.prompt_manager.in_template = s

def sage_prompt():
    """
    Returns the current Sage prompt.

    EXAMPLES::

        sage: from sage.misc.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell('sage_prompt()')
        u'sage: '
    """
    ipython = get_ipython()
    return ipython.prompt_manager.in_template

###############
# Displayhook #
###############
from IPython.core.displayhook import DisplayHook
class SageDisplayHook(DisplayHook):
    """
    A replacement for ``sys.displayhook`` which correctly print lists
    of matrices.

    EXAMPLES::

        sage: from sage.misc.interpreter import SageDisplayHook, get_test_shell
        sage: shell = get_test_shell()
        sage: shell.displayhook
        <sage.misc.interpreter.SageDisplayHook object at 0x...>
        sage: shell.run_cell('a = identity_matrix(ZZ, 2); [a,a]')
        [
        [1 0]  [1 0]
        [0 1], [0 1]
        ]
    """
    def compute_format_data(self, result):
        r"""
        Computes the format data of ``result``.  If the
        :func:`sage.misc.displayhook.print_obj` writes a string, then
        we override IPython's :class:`DisplayHook` formatting.

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.displayhook
            <sage.misc.interpreter.SageDisplayHook object at 0x...>
            sage: shell.displayhook.compute_format_data(2)
            {u'text/plain': '2'}
            sage: a = identity_matrix(ZZ, 2)
            sage: shell.displayhook.compute_format_data([a,a])
            {u'text/plain': '[\n[1 0]  [1 0]\n[0 1], [0 1]\n]'}
        """
        format_data = super(SageDisplayHook, self).compute_format_data(result)

        from cStringIO import StringIO
        s = StringIO()
        from sage.misc.displayhook import print_obj
        print_obj(s, result)
        if s:
            format_data['text/plain'] = s.getvalue().strip()
        return format_data


####################
# InteractiveShell #
####################
from IPython.frontend.terminal.interactiveshell import TerminalInteractiveShell
class SageInteractiveShell(TerminalInteractiveShell):

    # used by init_displayhook
    displayhook_class = Type(SageDisplayHook)

    # Input splitter, to split entire cells of input into either individual
    # interactive statements or whole blocks.
    input_splitter = Instance('sage.misc.interpreter.SageInputSplitter', (), {})

    def run_cell(self, *args, **kwds):
        r"""
        This method loads all modified attached files before executing
        any code.  Any arguments and keyword arguments are passed to
        :meth:`TerminalInteractiveShell.run_cell`.

        EXAMPLES::

            sage: import os
            sage: from sage.misc.interpreter import get_test_shell
            sage: from sage.misc.misc import tmp_dir
            sage: shell = get_test_shell()
            sage: tmp = os.path.join(tmp_dir(), 'run_cell.py')
            sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
            sage: shell.run_cell('%attach ' + tmp)
            sage: shell.run_cell('a')
            2
            sage: import time; time.sleep(1)
            sage: f = open(tmp, 'w'); f.write('a = 3\n'); f.close()
            sage: shell.run_cell('a')
            3
            sage: os.remove(tmp)
        """
        for f in modified_attached_files():
            super(SageInteractiveShell, self).run_cell('%%load "%s"'%f)
        return super(SageInteractiveShell, self).run_cell(*args, **kwds)

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

    def system_raw(self, cmd):
        """
        Adjust the libraries before calling system commands.  See Trac
        #975 for a discussion of this function.

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.system_raw('false')
            sage: shell.user_ns['_exit_code']
            256
        """
        sage_commands = os.listdir(os.environ['SAGE_ROOT']+"/local/bin/")
        DARWIN_SYSTEM = os.uname()[0]=='Darwin'
        if cmd in sage_commands:
            return super(SageInteractiveShell, self).system_raw(cmd)
        else:
            libraries = 'LD_LIBRARY_PATH=$$SAGE_ORIG_LD_LIBRARY_PATH;'
            if DARWIN_SYSTEM:
                libraries += 'DYLD_LIBRARY_PATH=$$SAGE_ORIG_DYLD_LIBRARY_PATH;'
            return super(SageInteractiveShell, self).system_raw(cmd)

    def run_ast_nodes(self, nodelist, cell_name, interactivity='last_expr'):
        """
        Override ``interactivity``

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('for i in range(3): i')  # indirect doctest
            0
            1
            2
        """
        super(SageInteractiveShell, self).run_ast_nodes(nodelist, cell_name, interactivity='all')

    #######################################
    # Magic functions
    #######################################
    def magic_timeit(self, s):
        """
        Runs :func:`sage_timeit` on the code s.  This is designed to
        be used from the command line as ``%timeit 2+2``.

        :param s: code to be timed
        :type s: string

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.magic_timeit('2+2') #random output
            625 loops, best of 3: 525 ns per loop
            sage: shell.magic_timeit('2r+2r').__class__
            <class sage.misc.sage_timeit.SageTimeitResult at 0x...>
        """
        from sage.misc.sage_timeit import sage_timeit
        return sage_timeit(s, self.user_ns)

    def magic_prun(self, parameter_s='', **kwds):
        """
        Profiles the code contained in ``parameter_s``. This is
        designed to be used from the command line as ``%prun 2+2``.

        :param parameter_s: code to be profiled
        :type parameter_s: string

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.magic_prun('2+2')
                     2 function calls in 0.000 seconds
            <BLANKLINE>
               Ordered by: internal time
            <BLANKLINE>
               ncalls  tottime  percall  cumtime  percall filename:lineno(function)
                    1    0.000    0.000    0.000    0.000 <string>:1(<module>)
                    1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
        """
        return super(SageInteractiveShell, self).magic_prun(parameter_s=preparse(parameter_s),
                                                            **kwds)

    def magic_load(self, s):
        r"""
        Loads the code contained in the file ``s``. This is designed
        to be used from the command line as ``%load /path/to/file``.

        :param s: file to be loaded
        :type s: string

        EXAMPLES::

            sage: import os
            sage: from sage.misc.interpreter import get_test_shell
            sage: from sage.misc.misc import tmp_dir
            sage: shell = get_test_shell()
            sage: tmp = os.path.join(tmp_dir(), 'run_cell.py')
            sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
            sage: shell.magic_load(tmp)
            sage: shell.run_cell('a')
            2
        """
        from sage.misc.preparser import load_wrap
        return self.ex(load_wrap(s, attach=False))

    def magic_attach(self, s):
        r"""
        Attaches the code contained in the file ``s``. This is
        designed to be used from the command line as
        ``%attach /path/to/file``.

        :param s: file to be attached
        :type s: string

        EXAMPLES::

            sage: import os
            sage: from sage.misc.interpreter import get_test_shell
            sage: from sage.misc.misc import tmp_dir
            sage: shell = get_test_shell()
            sage: tmp = os.path.join(tmp_dir(), 'run_cell.py')
            sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
            sage: shell.magic_attach(tmp)
            sage: shell.run_cell('a')
            2
            sage: import time; time.sleep(1)
            sage: f = open(tmp, 'a'); f.write('a = 3\n'); f.close()
            sage: shell.run_cell('a')
            3
        """
        from sage.misc.preparser import load_wrap
        return self.ex(load_wrap(s, attach=True))

    def magic_iload(self, s):
        """
        A magic command to interactively load a file as in MAGMA.

        :param s: the file to be interactively loaded
        :type s: string

        .. note::

            Currently, this cannot be completely doctested as it
            relies on :func:`raw_input`.

        EXAMPLES::

            sage: ip = get_ipython()           # not tested: works only in interactive shell
            sage: ip.magic_iload('/dev/null')  # not tested: works only in interactive shell
            Interactively loading "/dev/null"  # not tested: works only in interactive shell
        """
        try:
            name = str(eval(s))
        except Exception:
            name = s.strip()

        try:
            F = open(name)
        except IOError:
            raise ImportError, 'could not open file "%s"'%name

        #We need to update the execution count so that the history for the
        #iload command and the history for the first line of the loaded
        #file are not written to the history database with the same line
        #number (execution count).  This happens since the execution count
        #is updated only after the magic command is run.
        self.execution_count += 1

        print 'Interactively loading "%s"'%name

        # The following code is base on IPython's
        # InteractiveShell.interact,
        more = False
        for line in F.readlines():
            prompt = self.prompt_manager.render('in' if not more else 'in2', color=True)
            raw_input(prompt.encode('utf-8') + ' ' + line.rstrip())

            self.input_splitter.push(line)
            more = self.input_splitter.push_accepts_more()
            if not more:
                source, source_raw = self.input_splitter.source_raw_reset()
                self.run_cell(source_raw, store_history=True)


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
    def __call__(self, line, line_number, block_start):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        - ``line_number`` -- integer. The line number. For a single-line input, this is always zero.

        - ``block_start`` -- boolean. Whether the line is at the
          beginning of a new block.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import SagePreparseTransformer
            sage: spt = SagePreparseTransformer()
            sage: spt('2', 0, False)
            'Integer(2)'
            sage: preparser(False)
            sage: spt('2', 0, False)
            '2'
            sage: preparser(True)
        """
        if do_preparse and not line.startswith('%'):
            return preparse(line)
        else:
            return line

class LoadAttachTransformer():
    r"""
    Handle input lines that start out like ``load ...`` or ``attach
    ...``.

    Since there are objects in the Sage namespace named ``load`` and
    ``attach``, IPython's automagic will not transform these lines
    into ``%load ...`` and ``%attach ...``, respectively.  Thus, we
    have to do it manually.
    """
    def __call__(self, line, line_number, block_start):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        - ``line_number`` -- integer. The line number. For a single-line input, this is always zero.

        - ``block_start`` -- boolean. Whether the line is at the
          beginning of a new block.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell, LoadAttachTransformer
            sage: lat = LoadAttachTransformer()
            sage: lat('load /path/to/file', 0, False)
            '%load /path/to/file'
            sage: lat('attach /path/to/file', 0, False)
            '%attach /path/to/file'
        """
        if line_number>0:
            return line
        for cmd in ['load', 'attach']:
            if line.startswith(cmd + ' '):
                return '%' + line
        return line

class SagePromptTransformer():
    """
    Remove Sage prompts from the imput line.
    """
    _sage_prompt_re = re.compile(r'(^[ \t]*sage: |^[ \t]*\.+:? )+')

    def __call__(self, line, line_number, block_start):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        - ``line_number`` -- integer. The line number. For a single-line input, this is always zero.

        - ``block_start`` -- boolean. Whether the line is at the
          beginning of a new block.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import SagePromptTransformer
            sage: spt = SagePromptTransformer()
            sage: spt("sage: sage: 2 + 2", 0, False)
            '2 + 2'
            sage: spt('', 0, False)
            ''
            sage: spt("      sage: 2+2", 0, False)
            '2+2'
            sage: spt("      ... 2+2", 0, False)
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

class InterfaceMagicTransformer():
    """
    This transformer is for handling commands like ``%maxima`` to
    switch to a Maxima shell.
    """
    _magic_command_re = re.compile(r"get_ipython\(\).magic\(u'([^\d\W]\w+)'\)", re.UNICODE)

    def interfaces(self):
        """
        Return the list of interfaces

        OUTPUT:

        A list of stings.

        EXAMPLES::

            sage: from sage.misc.interpreter import InterfaceMagicTransformer
            sage: imt = InterfaceMagicTransformer()
            sage: imt.interfaces()
            ['LiE', 'Lisp', 'MuPAD', 'axiom', 'fricas', 'gap', 'gap3', 'giac',
             'kash', 'macaulay2', 'magma', 'maple', 'mathematica', 'matlab',
             'maxima', 'mwrank', 'octave', 'pari', 'r', 'sage', 'scilab', 'singular']
        """
        if '_interfaces' in self.__dict__:
            return self._interfaces
        import sage.interfaces
        self._interfaces = sorted([ obj.name()
                                    for obj in sage.interfaces.all.__dict__.values()
                                    if isinstance(obj, sage.interfaces.interface.Interface) ])
        return self._interfaces

    def __call__(self, line, line_number, block_start):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        - ``line_number`` -- integer. The line number. For a single-line input, this is always zero.

        - ``block_start`` -- boolean. Whether the line is at the
          beginning of a new block.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import InterfaceMagicTransformer
            sage: imt = InterfaceMagicTransformer()
            sage: imt('%maxima', 0, False)
            'maxima.interact()'
            sage: imt('%prun', 0, False)
            '%prun'
        """
        if line_number>0:
            return line
        if line.startswith('%'):
            interface = line[1:].strip()
            if interface in self.interfaces():
                return interface + '.interact()'
        return line

class SagePromptDedenter():
    """
    Remove leading spaces from the imput line.
    """
    def __call__(self, line, line_number, block_start):
        """
        Transform ``line``.

        INPUT:

        - ``line`` -- string. The line to be transformed.

        - ``line_number`` -- integer. The line number. For a single-line input, this is always zero.

        - ``block_start`` -- boolean. Whether the line is at the
          beginning of a new block.

        OUTPUT:

        A string, the transformed line.

        EXAMPLES::

            sage: from sage.misc.interpreter import SagePromptDedenter
            sage: spd = SagePromptDedenter()
            sage: spd('  1 + \\', 0, False)
            '1 + \\'
            sage: spd('  2',     1, False)
            '2'
            sage: spd('3',     2, False)   # try our best with incorrect indentation
            '3'
        """
        if line_number == 0:
            dedent_line = line.lstrip()
            self._dedent = len(line) - len(dedent_line)
            return dedent_line
        else:
            dedent = min(len(line)-len(line.lstrip()), self._dedent)
            return line[dedent:]

###################################################################
# Input Splitter
###################################################################
# the main entry point for input handling in Sage

from IPython.core.inputsplitter import InputSplitter
class SageInputSplitter(InputSplitter):
    """
    An input splitter for Sage syntax
    """
    # Whether a multi-line input is complete, i.e. line = empty line at the end
    _is_complete = False

    # String with raw, untransformed input.
    source_raw = ''

    ### Private attributes
    # List with lines of raw input accumulated so far.
    _buffer_raw = None
    _magic_interfaces = []

    def __init__(self, input_mode=None):
        """
        The Python constructor

        TESTS::

            sage: from sage.misc.interpreter import SageInputSplitter
            sage: SageInputSplitter()
            <sage.misc.interpreter.SageInputSplitter object at 0x...>
        """
        InputSplitter.__init__(self, input_mode)
        self._buffer_raw = []
        self._transforms = [
            SagePromptDedenter(),
            SagePromptTransformer(),
            InterfaceMagicTransformer(),
            LoadAttachTransformer(),
            SagePreparseTransformer() ]

    def reset(self):
        """
        Reset the input buffer and associated state.

        EXAMPLES::

            sage: from sage.misc.interpreter import SageInputSplitter
            sage: sis = SageInputSplitter()
            sage: sis.push('1')
            True
            sage: sis._buffer
            [u'Integer(1)\n']
            sage: sis.reset()
            sage: sis._buffer
            []
        """
        # print 'SageInputSplitter.reset buf = '+str(self._buffer)
        InputSplitter.reset(self)
        self._buffer_raw[:] = []
        self.source_raw = ''

    def source_raw_reset(self):
        """
        Return input and raw source and perform a full reset.

        EXAMPLES::

            sage: from sage.misc.interpreter import SageInputSplitter
            sage: sis = SageInputSplitter()
            sage: sis.push('1')
            True
            sage: sis._buffer
            [u'Integer(1)\n']
            sage: sis.source_raw_reset()
            (u'Integer(1)\n', u'1\n')
            sage: sis._buffer
            []
        """
        # print 'SageInputSplitter.source_raw_reset'
        out = self.source
        out_r = self.source_raw
        self.reset()
        return out, out_r

    def push(self, lines):
        """
        Push one or more lines of Sage input.

        This is mostly copy & paste from IPython, but with different transformers

        EXAMPLES::

            sage: from sage.misc.interpreter import SageInputSplitter
            sage: sis = SageInputSplitter()
            sage: sis._buffer
            []
            sage: sis.push('1')
            True
            sage: sis._buffer
            [u'Integer(1)\n']
        """
        # print 'SageInputSplitter.push '+lines+' ('+str(len(self._buffer))+')'  # dbg

        if not lines:
            return super(SageInputSplitter, self).push(lines)

        # We must ensure all input is pure unicode
        lines = cast_unicode(lines, self.encoding)
        lines_list = lines.splitlines()

        # Transform logic
        #
        # We only apply the line transformers to the input if we have either no
        # input yet, or complete input, or if the last line of the buffer ends
        # with ':' (opening an indented block).  This prevents the accidental
        # transformation of escapes inside multiline expressions like
        # triple-quoted strings or parenthesized expressions.
        #
        # The last heuristic, while ugly, ensures that the first line of an
        # indented block is correctly transformed.
        #
        # FIXME: try to find a cleaner approach for this last bit.
        # If we were in 'block' mode, since we're going to pump the parent
        # class by hand line by line, we need to temporarily switch out to
        # 'line' mode, do a single manual reset and then feed the lines one
        # by one.  Note that this only matters if the input has more than one
        # line.
        saved_input_mode = None
        if self.input_mode != 'line':
            saved_input_mode = self.input_mode
            self.reset()
            self.input_mode = 'line'

        # Store raw source before applying any transformations to it.  Note
        # that this must be done *after* the reset() call that would otherwise
        # flush the buffer.
        self._store(lines, self._buffer_raw, 'source_raw')

        push = super(SageInputSplitter, self).push
        buf = self._buffer
        try:
            for line in lines_list:
                line_number = len(buf)
                block_start = self._is_complete or not buf or \
                    (buf and (buf[-1].rstrip().endswith(':') or
                              buf[-1].rstrip().endswith(',')) )
                for f in self._transforms:
                    line = f(line, line_number, block_start)
                out = push(line)
        finally:
            if saved_input_mode is not None:
                self.input_mode = saved_input_mode
        return out


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

    EXAMPLES::

        sage: from sage.misc.interpreter import interface_shell_embed
        sage: shell = interface_shell_embed(gap)
        sage: shell.run_cell('List( [1..10], IsPrime )')
        [ false, true, true, false, true, false, true, false, false, false ]
    """
    try:
        cfg = Config(get_ipython().config)
    except NameError:
        cfg = Config(SageTerminalApp.DEFAULT_SAGE_CONFIG)
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

    :keyword sage_ext: load the Sage extension
    :type sage_ext: bool
    :returns: an IPython shell

    EXAMPLES::

        sage: from sage.misc.interpreter import get_test_shell
        sage: shell = get_test_shell(); shell
        <sage.misc.interpreter.SageInteractiveShell object at 0x...>
    """
    app = SageTerminalApp.instance()
    app.test_shell = True
    if app.shell is None:
        app.initialize()
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

class SageTerminalApp(TerminalIPythonApp):
    name = u'Sage'
    crash_handler_class = SageCrashHandler
    verbose_crash = True   # Needed to use Sage Crash Handler
    display_banner = False
    test_shell = False

    DEFAULT_SAGE_CONFIG = Config(
        {'PromptManager':    {'in_template':  u'sage: ',
                              'in2_template': u'....: ',
                              'out_template': u'',
                              'justify':      False},
         'InteractiveShell': {'separate_in':  u'',
                              'autoindent':   True,
                              'confirm_exit': False,
                              'colors':       'NoColor'},
         })

    def init_shell(self):
        r"""
        Initialize the :class:`SageInteractiveShell` instance.
        Additionally, this also does the following:

          - Merges the default shell configuration with the user's.

          - Monkey-patch :mod:`IPython.core.oinspect` to add Sage
            introspection functions.

          - Run additional Sage startup code.

        .. note::

            This code is based on
            :meth:`TermintalIPythonApp.init_shell`.

        EXAMPLES::

            sage: from sage.misc.interpreter import SageTerminalApp
            sage: app = SageTerminalApp.instance()
            sage: app.initialize() #indirect doctest
            sage: app.shell
            <sage.misc.interpreter.SageInteractiveShell object at 0x...>
        """
        import sys
        sys.path.insert(0, '')

        # Overwrite the default Sage configuration with the user's.
        new_config = Config()
        new_config._merge(self.DEFAULT_SAGE_CONFIG)
        new_config._merge(self.config)

        # Shell initialization
        self.shell = SageInteractiveShell.instance(config=new_config,
                        display_banner=False, profile_dir=self.profile_dir,
                        ipython_dir=self.ipython_dir)
        self.shell.configurables.append(self)


        # Ideally, these would just be methods of the Inspector class
        # that we could override; however, IPython looks them up in
        # the global :class:`IPython.core.oinspect` module namespace.
        # Thus, we have to monkey-patch.
        import sagedoc, sageinspect, IPython.core.oinspect
        IPython.core.oinspect.getdoc = sageinspect.sage_getdoc #sagedoc.my_getdoc
        IPython.core.oinspect.getsource = sagedoc.my_getsource
        IPython.core.oinspect.getargspec = sageinspect.sage_getargspec

        from sage.misc.edit_module import edit_devel
        self.shell.set_hook('editor', edit_devel)

        # Additional initalization code
        preparser(True)
        os.chdir(os.environ["CUR"])

        from sage.misc.misc import branch_current_hg_notice, branch_current_hg
        branch = branch_current_hg_notice(branch_current_hg())
        if branch and self.test_shell is False:
            print branch

        try:
            self.shell.ex('from sage.all import Integer, RealNumber')
        except Exception:
            import traceback
            print "Error importing the Sage library"
            traceback.print_exc()
            print
            print "To debug this, you can run:"
            print 'sage -ipython -i -c "import sage.all"'
            print 'and then type "%debug" to enter the interactive debugger'
            sys.exit(1)

        self.shell.push(dict(sage_prompt=sage_prompt))

        if os.environ.get('SAGE_IMPORTALL', 'yes') != 'yes':
            return

        self.shell.ex('from sage.all_cmdline import *')
        startup_file = os.environ.get('SAGE_STARTUP_FILE', '')

        if os.path.exists(startup_file):
            self.shell.run_cell('%%load "%s"'%startup_file)

