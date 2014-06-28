r"""
Sage's IPython Extension

A Sage extension which adds sage-specific features:

* magics

  - ``%crun``

  - ``%runfile``

  - ``%attach``

  - ``%display``

  - ``%mode`` (like ``%maxima``, etc.)

* preparsing of input

* loading Sage library

* running init.sage

* changing prompt to Sage prompt

* Display hook

TESTS:

We test that preparsing is off for ``%runfile``, on for ``%time``::

    sage: import os, re
    sage: from sage.repl.interpreter import get_test_shell
    sage: from sage.misc.misc import tmp_dir
    sage: shell = get_test_shell()
    sage: TMP = tmp_dir()

The temporary directory should have a name of the form
``.../12345/...``, to demonstrate that file names are not
preparsed when calling ``%runfile`` ::

    sage: bool(re.search('/[0-9]+/', TMP))
    True
    sage: tmp = os.path.join(TMP, 'run_cell.py')
    sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
    sage: shell.run_cell('%runfile '+tmp)
    sage: shell.run_cell('a')
    2

In contrast, input to the ``%time`` magic command is preparsed::

    sage: shell.run_cell('%time 594.factor()')
    CPU times: user ...
    Wall time: ...
    2 * 3^3 * 11
"""

from IPython.core.magic import Magics, magics_class, line_magic

import sage.misc.preparser as preparser

from sage.env import SAGE_IMPORTALL, SAGE_STARTUP_FILE

@magics_class
class SageMagics(Magics):

    @line_magic
    def crun(self, s):
        r"""
        Profile C function calls

        INPUT:

        - ``s`` -- string. Sage command to profile.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%crun sum(1/(1+n^2) for n in range(100))')   # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
            Using local file ...
            Using local file ...
        """
        import sage.misc.gperftools
        sage.misc.gperftools.crun(s, evaluator=self.shell.ex)

    @line_magic
    def runfile(self, s):
        r"""
        Execute the code contained in the file ``s``.

        This is designed to be used from the command line as
        ``%runfile /path/to/file``.

        - ``s`` -- string. The file to be loaded.

        EXAMPLES::

            sage: import os
            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.misc.misc import tmp_dir
            sage: shell = get_test_shell()
            sage: tmp = os.path.join(tmp_dir(), 'run_cell.py')
            sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
            sage: shell.run_cell('%runfile '+tmp)
            sage: shell.run_cell('a')
            2
        """
        return self.shell.ex(preparser.load_wrap(s, attach=False))

    @line_magic
    def attach(self, s):
        r"""
        Attach the code contained in the file ``s``.

        This is designed to be used from the command line as ``%attach
        /path/to/file``.

        - ``s`` -- string. The file to be attached

        EXAMPLES::

            sage: import os
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: tmp = os.path.normpath(os.path.join(SAGE_TMP, 'run_cell.py'))
            sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
            sage: shell.run_cell('%attach ' + tmp)
            sage: shell.run_cell('a')
            2
            sage: sleep(1)  # filesystem timestamp granularity
            sage: f = open(tmp, 'w'); f.write('a = 3\n'); f.close()

        Note that the doctests are never really at the command prompt, so
        we call the input hook manually::

            sage: shell.run_cell('from sage.repl.inputhook import sage_inputhook')
            sage: shell.run_cell('sage_inputhook()')
            ### reloading attached file run_cell.py modified at ... ###
            0

            sage: shell.run_cell('a')
            3
            sage: shell.run_cell('detach(%r)'%tmp)
            sage: shell.run_cell('attached_files()')
            []
            sage: os.remove(tmp)
        """
        return self.shell.ex(preparser.load_wrap(s, attach=True))

    @line_magic
    def iload(self, args):
        """
        A magic command to interactively load a file as in MAGMA.

        - ``args`` -- string. The file to be interactively loaded

        .. note::

            Currently, this cannot be completely doctested as it
            relies on :func:`raw_input`.

        EXAMPLES::

            sage: ip = get_ipython()           # not tested: works only in interactive shell
            sage: ip.magic_iload('/dev/null')  # not tested: works only in interactive shell
            Interactively loading "/dev/null"  # not tested: works only in interactive shell
        """
        content = self.shell.find_user_code(args).splitlines()

        # we create a stack so e.g. having an iload inside of an iload
        # will process the inner iload and then resume the outer iload
        orig_readline = self.shell.pre_readline

        def pre_readline():
            if self.shell.rl_next_input is None:
                self.shell.rl_next_input = content.pop(0)
                self.shell.rl_do_indent = False
            orig_readline()
            if not content:
                # restore original hook
                self.shell.readline_startup_hook(orig_readline)
                self.shell.pre_readline = orig_readline

        self.shell.readline_startup_hook(pre_readline)
        self.shell.pre_readline = pre_readline

        print('Interactively loading "%s"'%args)

    _magic_display_status = 'simple'
    @line_magic
    def display(self, args):
        r"""
        A magic command to switch between simple display and ASCII art display.

        - ``args`` -- string.  See
          :meth:`sage.misc.display_hook.DisplayHookBase.set_display`
          for allowed values. If the mode is ``ascii_art``, it can
          optionally be followed by a width.

        How to use: if you want activate the ASCII art mod::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%display ascii_art')

        That means you don't have to use :func:`ascii_art` to get an ASCII art
        output::

            sage: shell.run_cell("i = var('i')")
            sage: shell.run_cell('sum(i^2*x^i, i, 0, 10)')
                 10       9       8       7       6       5       4      3      2
            100*x   + 81*x  + 64*x  + 49*x  + 36*x  + 25*x  + 16*x  + 9*x  + 4*x  + x

        Then when you want return in 'textual mode'::

            sage: shell.run_cell('%display simple')
            sage: shell.run_cell('sum(i^2*x^i, i, 0, 10)')
            100*x^10 + 81*x^9 + 64*x^8 + 49*x^7 + 36*x^6 + 25*x^5 + 16*x^4 + 9*x^3 + 4*x^2 + x

        Sometime you could have to use a special output width and you
        could specify it::

            sage: shell.run_cell('%display ascii_art')
            sage: shell.run_cell('StandardTableaux(4).list()')
            [
            [                                                                  1  4    1  3
            [                 1  3  4    1  2  4    1  2  3    1  3    1  2    2       2
            [   1  2  3  4,   2      ,   3      ,   4      ,   2  4,   3  4,   3   ,   4
            <BLANKLINE>
                        1 ]
                1  2    2 ]
                3       3 ]
            ,   4   ,   4 ]
            sage: shell.run_cell('%display ascii_art 50')
            sage: shell.run_cell('StandardTableaux(4).list()')
            [
            [
            [                 1  3  4    1  2  4    1  2  3
            [   1  2  3  4,   2      ,   3      ,   4      ,
            <BLANKLINE>
                                                      1 ]
                              1  4    1  3    1  2    2 ]
              1  3    1  2    2       2       3       3 ]
              2  4,   3  4,   3   ,   4   ,   4   ,   4 ]

        As yet another option, typeset mode. This is used in the emacs
        interface::

            sage: shell.run_cell('%display typeset')
            sage: shell.run_cell('1/2')
            <html><script type="math/tex">...\frac{1}{2}</script></html>

        Switch back::

            sage: shell.run_cell('%display simple')

        TESTS::

            sage: shell.run_cell('%display invalid_mode')
            ---------------------------------------------------------------------------
            ValueError                                Traceback (most recent call last)
            ...
            ValueError: invalid mode set
        """
        args = args.strip().split()
        if not args:
            # toggle
            self._magic_display_status = ('ascii_art' if
                    self._magic_display_status == 'simple' else 'simple')
        else:
            mode = args[0]
            self._magic_display_status = mode

        if self._magic_display_status == 'ascii_art' and len(args) > 1:
            try:
                max_width = int(args[1])
            except ValueError:
                max_width = 0
            if max_width <= 0:
                raise ValueError(
                        "ascii_art max width must be a positive integer")
            import sage.misc.ascii_art as ascii_art
            ascii_art.MAX_WIDTH = max_width

        self.shell.display_formatter.formatters['text/plain'].set_display(
                self._magic_display_status)

class SageCustomizations(object):

    def __init__(self, shell=None):
        """
        Initialize the Sage plugin.
        """
        self.shell = shell

        self.auto_magics = SageMagics(shell)
        self.shell.register_magics(self.auto_magics)

        import sage.misc.displayhook as displayhook
        self.shell.display_formatter.formatters['text/plain'] = (
                displayhook.SagePlainTextFormatter(config=shell.config))

        import sage.misc.edit_module as edit_module
        self.shell.set_hook('editor', edit_module.edit_devel)

        self.init_inspector()
        self.init_line_transforms()

        import inputhook
        inputhook.install()

        import sage.all # until sage's import hell is fixed

        self.shell.verbose_quit = True
        self.set_quit_hook()

        self.register_interface_magics()

        if SAGE_IMPORTALL == 'yes':
            self.init_environment()


    def register_interface_magics(self):
        """
        Register magics for each of the Sage interfaces
        """
        from sage.misc.superseded import deprecation
        import sage.interfaces.all
        interfaces = [(name, obj)
                      for name, obj in sage.interfaces.all.__dict__.items()
                      if isinstance(obj, sage.interfaces.interface.Interface)]

        for real_name, obj in interfaces:
            def tmp(line, name=real_name):
                self.shell.run_cell('%s.interact()' % name)
            tmp.__doc__ = "Interact with %s" % real_name
            self.shell.register_magic_function(tmp, magic_name=real_name)

            obj_name = obj.name()
            if real_name != obj_name:
                def tmp_deprecated(line, name=real_name, badname=obj_name):
                    deprecation(6288, 'Use %%%s instead of %%%s.' % (name,
                                                                     badname))
                    self.shell.run_cell('%s.interact()' % name)
                tmp_deprecated.__doc__ = "Interact with %s" % real_name
                self.shell.register_magic_function(tmp_deprecated, magic_name=obj_name)

    def set_quit_hook(self):
        """
        Set the exit hook to cleanly exit Sage.
        """
        def quit():
            import sage.all
            sage.all.quit_sage(self.shell.verbose_quit)
        import atexit
        atexit.register(quit)

    def init_environment(self):
        """
        Set up Sage command-line environment
        """
        # import outside of cell so we don't get a traceback
        import sage.all_cmdline as sage
        self.shell.user_ns.update((key, getattr(sage, key))
                for key in dir(sage) if key[0] != '_')
        self.run_init()

    def run_init(self):
        """
        Run Sage's initial startup file.
        """
        try:
            with open(SAGE_STARTUP_FILE, 'r') as f:
                self.shell.run_cell(f.read(), store_history=False)
        except IOError:
            pass

    def init_inspector(self):
        # Ideally, these would just be methods of the Inspector class
        # that we could override; however, IPython looks them up in
        # the global :class:`IPython.core.oinspect` module namespace.
        # Thus, we have to monkey-patch.
        import sage.misc.sagedoc as sagedoc
        import sage.misc.sageinspect as sageinspect
        import IPython.core.oinspect
        IPython.core.oinspect.getdoc = sageinspect.sage_getdoc
        IPython.core.oinspect.getsource = sagedoc.my_getsource
        IPython.core.oinspect.getargspec = sageinspect.sage_getargspec

    def init_line_transforms(self):
        """
        Set up transforms (like the preparser).
        """
        from interpreter import (SagePreparseTransformer,
                                 SagePromptTransformer)

        for s in (self.shell.input_splitter, self.shell.input_transformer_manager):
            s.physical_line_transforms.insert(1, SagePromptTransformer())
            s.python_line_transforms.append(SagePreparseTransformer())

# from http://stackoverflow.com/questions/4103773/efficient-way-of-having-a-function-only-execute-once-in-a-loop
from functools import wraps
def run_once(func):
    """
    Runs a function (successfully) only once.

    The running can be reset by setting the ``has_run`` attribute to False

    TEST::

        sage: from sage.repl.ipython_extension import run_once
        sage: @run_once
        ....: def foo(work):
        ....:     if work:
        ....:         return 'foo worked'
        ....:     raise RuntimeError("foo didn't work")
        sage: foo(False)
        Traceback (most recent call last):
        ...
        RuntimeError: foo didn't work
        sage: foo(True)
        'foo worked'
        sage: foo(False)
        sage: foo(True)
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if not wrapper.has_run:
            result = func(*args, **kwargs)
            wrapper.has_run = True
            return result
    wrapper.has_run = False
    return wrapper


@run_once
def load_ipython_extension(ip):
    """
    Load the extension in IPython.
    """
    # this modifies ip
    SageCustomizations(shell=ip)
