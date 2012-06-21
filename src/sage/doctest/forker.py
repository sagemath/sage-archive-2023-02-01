"""
This module controls the processes started by Sage that actually run
the doctests.

EXAMPLES:

The following examples are used in doctesting this file::

    sage: doctest_var = 42; doctest_var^2
    1764
    sage: R.<a> = ZZ[]
    sage: a + doctest_var
    a + 42

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.
"""

#*****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import hashlib, multiprocessing, os, sys, time, warnings, signal, subprocess, linecache
from collections import defaultdict
import doctest, pdb, traceback
from Queue import Empty
import sage.misc.randstate as randstate
from util import Timer, RecordingDict, count_noun
from sources import DictAsObject
from StringIO import StringIO
from parsing import OriginalSource, reduce_hex
from sage.structure.sage_object import SageObject
from parsing import SageOutputChecker, pre_hash, get_source
from sage.misc.misc import walltime

################
from sage.structure.sage_object import dumps, loads
################

debug_semaphore = None
output_semaphore = None

def init_sage():
    """
    Import the Sage library.

    This function is called once at the beginning of a doctest run
    (rather than once for each file).  It imports the Sage library,
    sets DOCTEST_MODE to True, and invalidates any interfaces.

    EXAMPLES::

        sage: from sage.doctest.forker import init_sage
        sage: sage.plot.plot.DOCTEST_MODE = False
        sage: init_sage()
        sage: sage.plot.plot.DOCTEST_MODE
        True
    """
    # Do this once before forking.
    import sage.all_cmdline
    sage.plot.plot; sage.plot.plot.DOCTEST_MODE=True
    sage.interfaces.quit.invalidate_all()
    import sage.misc.displayhook
    sys.displayhook = sage.misc.displayhook.DisplayHook(sys.displayhook)

def warning_function(file):
    """
    Creates a function that prints warnings to the given file.

    INPUT:

    - ``file`` -- an open file handle.

    OUPUT:

    - a function that prings warnings to the given file.

    EXAMPLES::

        sage: from sage.doctest.forker import warning_function
        sage: import os
        sage: F = os.tmpfile()
        sage: wrn = warning_function(F)
        sage: wrn("bad stuff", UserWarning, "myfile.py", 0)
        sage: F.seek(0)
        sage: F.read()
        'doctest:0: UserWarning: bad stuff\n'
    """
    def doctest_showwarning(message, category, filename, lineno, file=file, line=None):
        try:
            file.write(warnings.formatwarning(message, category, 'doctest', lineno, line))
        except IOError:
            pass # the file (probably stdout) is invalid
    return doctest_showwarning

def _test_activate_stdin():
    """
    This function is designed to test
    :meth:`SageSpoofOut.activate_stdin` and
    :meth:`SageSpoofOut.deactivate_stdin`.

    TESTS::

        sage: import os
        sage: os.environ['SAGE_PEXPECT_LOG'] = "1"
        sage: sage0.quit()
        sage: _ = sage0.eval("from sage.doctest.forker import _test_activate_stdin")
        sage: sage0._prompt = "my_prompt: "
        sage: print sage0.eval("_test_activate_stdin()")
        bad_input: input not available before activation
        sage: sage0._prompt = "sage: "
        sage: print sage0.eval("input active")
        input received: input active
        bad_input: input not available after deactivation
    """
    O = os.tmpfile()
    save_stdin = sys.stdin
    sys.stdin = open(os.devnull)
    v = None
    try:
        S = SageSpoofOut(O)
        try:
            w = raw_input("bad_input: ")
        except EOFError:
            print "input not available before activation"
        S.activate_stdin()
        v = raw_input("my_prompt: ")
    finally:
        print "input received:", v
        S.deactivate_stdin()
        try:
            w = raw_input("bad_input: ")
        except EOFError:
            print "input not available after deactivation"
        sys.stdin = save_stdin

class SageSpoofOut(SageObject):
    """
    We replace the standard :class:`doctest._SpoofOut` for two reasons:

    - we need to divert the output of C programs that don't print through sys.stdout,
    - we want the ability to recover partial output from doctest processes that segfault.

    INPUT:

    - ``outfile`` -- an open file handle (usually from os.tmpfile) to which stdout should be redirected.

    EXAMPLES::

        sage: import os
        sage: from sage.doctest.forker import SageSpoofOut
        sage: O = os.tmpfile()
        sage: S = SageSpoofOut(O)
        sage: try:
        ....:     S.start_spoofing()
        ....:     print "hello world"
        ....: finally:
        ....:     S.stop_spoofing()
        ....:
        sage: S.getvalue()
        'hello world\n'
    """
    def __init__(self, outfile):
        """
        Initialization.

        TESTS::

            sage: import os
            sage: from sage.doctest.forker import SageSpoofOut
            sage: SageSpoofOut(os.tmpfile())
            <class 'sage.doctest.forker.SageSpoofOut'>
        """
        self.outfile = outfile
        self.spoofing = False
        self.stdin_active = False #(sys.stdin.fileno() == 0)
        #self.stdin_active_original = self.stdin_active
        self.dup_stdout = os.dup(sys.stdout.fileno())
        self.dup_stderr = os.dup(sys.stderr.fileno())
        self.dup_stdin = os.dup(sys.stdin.fileno())
        self.position = 0

    def fileno(self):
        """
        This object serves as stdout in various places, so we need to
        emulate stdout's fileno in order to be able to recurse.

        EXAMPLES::

            sage: import os
            sage: from sage.doctest.forker import SageSpoofOut
            sage: O = os.tmpfile()
            sage: S = SageSpoofOut(O)
            sage: S.fileno() == 1
            True
        """
        return sys.stdout.fileno()

    def start_spoofing(self):
        """
        Set stdout to print to outfile.

        EXAMPLES::

            sage: import os
            sage: from sage.doctest.forker import SageSpoofOut
            sage: O = os.tmpfile()
            sage: S = SageSpoofOut(O)
            sage: try:
            ....:     S.start_spoofing()
            ....:     print "this is not printed"
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'this is not printed\n'

        We also catch C output::

            sage: try:
            ....:     S.start_spoofing()
            ....:     retval = os.system('''echo "Hello there"\nif [ $? -eq 0 ]; then\necho "good"\nfi''')
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'Hello there\ngood\n'
        """
        if not self.spoofing:
            sys.stdout.flush()
            os.dup2(self.outfile.fileno(), sys.stdout.fileno())
            sys.stderr.flush()
            os.dup2(self.outfile.fileno(), sys.stderr.fileno())
            self.spoofing = True

    def stop_spoofing(self):
        """
        Reset stdout to its original value.

        EXAMPLES::

            sage: import os
            sage: from sage.doctest.forker import SageSpoofOut
            sage: O = os.tmpfile()
            sage: S = SageSpoofOut(O)
            sage: try:
            ....:     S.start_spoofing()
            ....:     print "this is not printed"
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: print "this is now printed"
            this is now printed
        """
        if self.spoofing:
            sys.stdout.flush()
            os.dup2(self.dup_stdout, sys.stdout.fileno())
            sys.stderr.flush()
            os.dup2(self.dup_stderr, sys.stderr.fileno())
            self.spoofing = False

    def getvalue(self):
        """
        Gets the value that has been printed to the outfile since the last time this function was called.

        EXAMPLES::

            sage: import os
            sage: from sage.doctest.forker import SageSpoofOut
            sage: O = os.tmpfile()
            sage: S = SageSpoofOut(O)
            sage: try:
            ....:     S.start_spoofing()
            ....:     print "step 1"
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'step 1\n'
            sage: try:
            ....:     S.start_spoofing()
            ....:     print "step 2"
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'step 2\n'
        """
        sys.stdout.flush()
        self.outfile.seek(self.position)
        result = self.outfile.read()
        self.position = self.outfile.tell()
        if not result.endswith("\n"):
            result += "\n"
        return result

    def write(self, a):
        """
        When debugging we need to pretend this object is stdout.

        EXAMPLES::

            sage: import os
            sage: from sage.doctest.forker import SageSpoofOut
            sage: O = os.tmpfile()
            sage: S = SageSpoofOut(O)
            sage: try:
            ....:     S.start_spoofing()
            ....:     S.write("this is not printed\n")
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'this is not printed\n'
        """
        sys.stdout.write(a)
        sys.stdout.flush()

    def unspoofed_write(self, a):
        """
        Prints to the unspoofed standard out.

        EXAMPLES::

            sage: import os
            sage: from sage.doctest.forker import SageSpoofOut
            sage: O = os.tmpfile()
            sage: S = SageSpoofOut(O)
            sage: try:
            ....:     S.start_spoofing()
            ....:     print "spoofed"
            ....:     S.unspoofed_write("unspoofed\n")
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            unspoofed
            sage: S.getvalue()
            'spoofed\n'
        """
        spoofed = self.spoofing
        try:
            if output_semaphore is not None:
                output_semaphore.acquire()
            if spoofed:
                self.stop_spoofing()
            sys.stdout.write(a)
        finally:
            sys.stdout.flush()
            if spoofed:
                self.start_spoofing()
            if output_semaphore is not None:
                output_semaphore.release()

    def activate_stdin(self):
        """
        Turns on stdin, which is useful when user input is needed in a subprocess.

        TESTS::

            sage: import os
            sage: os.environ['SAGE_PEXPECT_LOG'] = "1"
            sage: sage0.quit()
            sage: _ = sage0.eval("from sage.doctest.forker import _test_activate_stdin")
            sage: sage0._prompt = "my_prompt: "
            sage: print sage0.eval("_test_activate_stdin()")
            bad_input: input not available before activation
            sage: sage0._prompt = "sage: "
            sage: print sage0.eval("input active") # indirect doctest
            input received: input active
            bad_input: input not available after deactivation
        """
        if not self.stdin_active:
            os.dup2(0, sys.stdin.fileno())
            self.stdin_active = True

    def deactivate_stdin(self):
        """
        Restores stdin to its original state (which ignores user input
        if this function is being called in a subprocess).

        TESTS::

            sage: import os
            sage: os.environ['SAGE_PEXPECT_LOG'] = "1"
            sage: sage0.quit()
            sage: _ = sage0.eval("from sage.doctest.forker import _test_activate_stdin")
            sage: sage0._prompt = "my_prompt: "
            sage: print sage0.eval("_test_activate_stdin()")
            bad_input: input not available before activation
            sage: sage0._prompt = "sage: "
            sage: print sage0.eval("input active") # indirect doctest
            input received: input active
            bad_input: input not available after deactivation
        """
        if self.stdin_active:
            os.dup2(self.dup_stdin, sys.stdin.fileno())
            self.stdin_active = False

class SagePdb(doctest._OutputRedirectingPdb):
    """
    We subclass :class:`doctest._OutputRedirectingPdb` since switching
    in and out of spoofing can't be done just be assigning to
    sys.stdout.

    INPUT:

    - ``fakeout`` -- a :class:`SageSpoofOut` instance

    EXAMPLES::

        sage: from sage.doctest.forker import SageSpoofOut, SagePdb
        sage: import os; O = os.tmpfile()
        sage: S = SageSpoofOut(O)
        sage: SagePdb(S)
        <sage.doctest.forker.SagePdb instance at ...>
    """
    def __init__(self, fakeout):
        """
        Initialization.

        TESTS::

            sage: from sage.doctest.forker import SageSpoofOut, SagePdb
            sage: import os; O = os.tmpfile()
            sage: S = SageSpoofOut(O)
            sage: SagePdb(S).use_rawinput
            1
        """
        self._fakeout = fakeout
        doctest._OutputRedirectingPdb.__init__(self, fakeout)

    def trace_dispatch(self, *args):
        """
        We stop spoofing before calling :meth:`pdb.Pdb.trace_dispatch` and start again afterward.

        TESTS::

            sage: import os
            sage: os.environ['SAGE_PEXPECT_LOG'] = "1"
            sage: sage0.quit()
            sage: _ = sage0.eval("import doctest, sys, os, multiprocessing, subprocess")
            sage: _ = sage0.eval("from sage.doctest.parsing import SageOutputChecker")
            sage: _ = sage0.eval("import sage.doctest.forker as sdf")
            sage: _ = sage0.eval("from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()")
            sage: _ = sage0.eval("ex = doctest.Example('E = EllipticCurve([0,0]); E', 'A singular Elliptic Curve')")
            sage: _ = sage0.eval("DT = doctest.DocTest([ex], globals(), 'singular_curve', None, 0, None)")
            sage: _ = sage0.eval("DTR = sdf.SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)")
            sage: _ = sage0.eval("sdf.debug_semaphore = multiprocessing.Semaphore()")
            sage: _ = sage0.eval("sdf.output_semaphore = multiprocessing.Semaphore()")
            sage: sage0._prompt = r"\(Pdb\) "
            sage: sage0.eval("DTR.run(DT, clear_globs=False)") # indirect doctest
            '... "Invariants %s define a singular curve."%ainvs'
            sage: sage0._prompt = "sage: "
            sage: tb = sage0.eval("quit")
        """
        if self._fakeout.spoofing:
            self._fakeout.stop_spoofing()
            try:
                return pdb.Pdb.trace_dispatch(self, *args)
            finally:
                self._fakeout.start_spoofing()
        else:
            return pdb.Pdb.trace_dispatch(self, *args)

class SageDocTestRunner(doctest.DocTestRunner):
    def __init__(self, *args, **kwds):
        """
        A customized version of DocTestRunner that tracs dependencies of doctests.

        INPUT:

        - ``stdout`` -- an open file to restore for debugging

        - ``checker`` -- None, or an instance of
          :class:`doctest.OutputChecker`

        - ``verbose`` -- boolean, determines whether verbose printing
          is enabled.

        - ``optionflags`` -- Controls the comparison with the expected
          output.  See :mod:`testmod` for more information.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: DTR
            <sage.doctest.forker.SageDocTestRunner instance at ...>
        """
        O = kwds.pop('output')
        self.options = kwds.pop('sage_options')
        doctest.DocTestRunner.__init__(self, *args, **kwds)
        self._fakeout = SageSpoofOut(O)
        self.history = []
        self.references = []
        self.setters = {}
        self.running_global_digest = hashlib.md5()
        self.delayed_output = []

    def _run(self, test, compileflags, out):
        """
        This function replaces :meth:`doctest.DocTestRunner.__run`.

        It changes the following behavior:

        - We call :meth:`SageDocTestRunner.execute` rather than just
          exec

        - We don't truncate _fakeout after each example since we want
          the output file to be readable by the calling
          :class:`SageWorker`.

        Since it needs to be able to read stdout, it should be called
        while spoofing using :class:`SageSpoofOut`.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: DTR.run(doctests[0], clear_globs=False) # indirect doctest
            TestResults(failed=0, attempted=4)
        """
        # Keep track of the number of failures and tries.
        failures = tries = 0

        # Save the option flags (since option directives can be used
        # to modify them).
        original_optionflags = self.optionflags

        SUCCESS, FAILURE, BOOM = range(3) # `outcome` state

        check = self._checker.check_output

        # Process each example.
        for examplenum, example in enumerate(test.examples):

            # If REPORT_ONLY_FIRST_FAILURE is set, then suppress
            # reporting after the first failure.
            quiet = (self.optionflags & doctest.REPORT_ONLY_FIRST_FAILURE and
                     failures > 0)

            # Merge in the example's options.
            self.optionflags = original_optionflags
            if example.options:
                for (optionflag, val) in example.options.items():
                    if val:
                        self.optionflags |= optionflag
                    else:
                        self.optionflags &= ~optionflag

            # If 'SKIP' is set, then skip this example.
            if self.optionflags & doctest.SKIP:
                continue

            # Record that we started this example.
            tries += 1
            # We print the example we're running for easier debugging
            # if this file times out or sefaults.
            with OriginalSource(example):
                print "sage: " + example.source[:-1] + " ## line %s ##"%(test.lineno + example.lineno + 1)
            # Update the position so that result comparison works
            throwaway = self._fakeout.getvalue()
            if not quiet:
                self.report_start(out, test, example)

            # Use a special filename for compile(), so we can retrieve
            # the source code during interactive debugging (see
            # __patched_linecache_getlines).
            filename = '<doctest %s[%d]>' % (test.name, examplenum)

            # Run the example in the given context (globs), and record
            # any exception that gets raised.  (But don't intercept
            # keyboard interrupts.)
            try:
                # Don't blink!  This is where the user's code gets run.
                compiled = compile(example.source, filename, "single",
                             compileflags, 1)
                self.execute(example, compiled, test.globs)
                self.debugger.set_continue() # ==== Example Finished ====
                exception = None
            except KeyboardInterrupt:
                # We check to see if it's expected
                exception = sys.exc_info()
                exc_msg = traceback.format_exception_only(*exception[:2])[-1]
                if example.exc_msg is None or not check(example.exc_msg, exc_msg, self.optionflags):
                    raise
                else:
                    # KeyboardInterrupt was expected
                    self.debugger.set_continue()
            except Exception:
                exception = sys.exc_info()
                self.debugger.set_continue() # ==== Example Finished ====

            got = self._fakeout.getvalue()  # the actual output
            outcome = FAILURE   # guilty until proved innocent or insane

            # If the example executed without raising any exceptions,
            # verify its output.
            if exception is None:
                if check(example.want, got, self.optionflags):
                    outcome = SUCCESS

            # The example raised an exception:  check if it was expected.
            else:
                exc_info = sys.exc_info()
                exc_msg = traceback.format_exception_only(*exc_info[:2])[-1]
                if not quiet:
                    got += doctest._exception_traceback(exc_info)

                # If `example.exc_msg` is None, then we weren't expecting
                # an exception.
                if example.exc_msg is None:
                    outcome = BOOM

                # We expected an exception:  see whether it matches.
                elif check(example.exc_msg, exc_msg, self.optionflags):
                    outcome = SUCCESS

                # Another chance if they didn't care about the detail.
                elif self.optionflags & doctest.IGNORE_EXCEPTION_DETAIL:
                    m1 = re.match(r'(?:[^:]*\.)?([^:]*:)', example.exc_msg)
                    m2 = re.match(r'(?:[^:]*\.)?([^:]*:)', exc_msg)
                    if m1 and m2 and check(m1.group(1), m2.group(1),
                                           self.optionflags):
                        outcome = SUCCESS

            # Report the outcome.
            if outcome is SUCCESS:
                if self.options.warn_long and walltime(self.tmp_time) > self.options.warn_long:
                    self.report_overtime(out, test, example, got)
                    failures += 1
                elif not quiet:
                    self.report_success(out, test, example, got)
            elif outcome is FAILURE:
                if not quiet:
                    self.report_failure(out, test, example, got, test.globs)
                failures += 1
            elif outcome is BOOM:
                if not quiet:
                    self.report_unexpected_exception(out, test, example,
                                                     exc_info)
                failures += 1
            else:
                assert False, ("unknown outcome", outcome)

        # Restore the option flags (in case they were modified)
        self.optionflags = original_optionflags

        # Record and return the number of failures and tries.
        self._DocTestRunner__record_outcome(test, failures, tries)
        return doctest.TestResults(failures, tries)

    def run(self, test, compileflags=None, out=None, clear_globs=True):
        """
        Runs the examples in a given doctest.

        This function replaces :class:`doctest.DocTestRunner.run`
        since it needs to handle spoofing and the output redirecting
        Pdb differently.  It also leaves the display hook in place.

        INPUT:

        - ``test`` -- an instance of :class:`doctest.DocTest`

        - ``compileflags`` -- the set of compiler flags used to
          execute examples (passed in to the :func:`compile`).  If
          None, they are filled in from the result of
          :func:`doctest._extract_future_flags` applied to
          ``test.globs``.

        - ``out`` -- a function for writing the output (defaults to
          :func:`sys.stdout.write`).

        - ``clear_globs`` -- boolean (default True): whether to clear
          the namespace after running this doctest.

        OUTPUT:

        - ``f`` -- integer, the number of examples that failed

        - ``t`` -- the number of examples tried

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: DTR.run(doctests[0], clear_globs=False)
            TestResults(failed=0, attempted=4)
        """
        self.setters = {}
        randstate.set_random_seed(long(0))
        warnings.showwarning = warning_function(sys.stdout)
        self.running_doctest_digest = hashlib.md5()
        self.test = test
        if compileflags is None:
            compileflags = doctest._extract_future_flags(test.globs)
        save_set_trace = pdb.set_trace
        self.debugger = SagePdb(self._fakeout)
        self.debugger.reset()
        pdb.set_trace = self.debugger.set_trace
        self.save_linecache_getlines = linecache.getlines
        linecache.getlines = self._DocTestRunner__patched_linecache_getlines
        if out is None and (self.options.serial or self.options.verbose or self.options.debug):
            out = self._fakeout.unspoofed_write
        elif out is None:
            def out(s):
                self.delayed_output.append(s)

        self._fakeout.start_spoofing()
        # If self.options.initial is set, we show only the first failure in each doctest block.
        self.no_failure_yet = True
        try:
            return self._run(test, compileflags, out)
        finally:
            self._fakeout.stop_spoofing()
            pdb.set_trace = save_set_trace
            linecache.getlines = self.save_linecache_getlines
            if clear_globs:
                test.globs.clear()

    def summarize(self, verbose=None, delay_list=None):
        """
        Returns a string summarizing the results of testing.

        Summarize works a bit differently for Sage testing since we
        want to continue delaying output so that it is printed at the
        same time as
        :meth:`sage.doctest.reporting.DocTestReporter.report`, and we
        want to log it.

        INPUT:

        - ``verbose`` -- whether to print lots of stuff

        - ``delay_list`` -- either None or a list of strings on which to append delayed output

        OUTPUT:

        - if ``delay_list`` is None, prints a summary of tests run.
          Otherwise, appends lines of the summary to ``delay_list``.

        - returns ``(f, t)``, a :class:`doctest.TestResults` instance
          giving the number of failures and the total number of tests
          run.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: DTR.delayed_output = ["A delayed failure\n","And another one\n"]
            sage: DTR._name2ft['sage.doctest.forker'] = (1,120)
            sage: DTR.summarize()
            A delayed failure
            And another one
            **********************************************************************
            1 item had failures:
                1 of 120 in sage.doctest.forker
            TestResults(failed=1, attempted=120)
        """
        if verbose is None:
            verbose = self._verbose
        print_now = delay_list is None
        if delay_list is None:
            delay_list = []
        delay_list.extend(self.delayed_output)
        notests = []
        passed = []
        failed = []
        totalt = totalf = 0
        for x in self._name2ft.items():
            name, (f, t) = x
            assert f <= t
            totalt += t
            totalf += f
            if not t:
                notests.append(name)
            elif not f:
                passed.append( (name, t) )
            else:
                failed.append(x)
        if verbose:
            if notests:
                delay_list.append(count_noun(len(notests), "item") + " had no tests:\n")
                notests.sort()
                for thing in notests:
                    delay_list.append("    %s\n"%thing)
            if passed:
                delay_list.append(count_noun(len(passed), "item") + " passed all tests:\n")
                passed.sort()
                for thing, count in passed:
                    delay_list.append(" %s in %s\n"%(count_noun(count, "test", pad_number=3, pad_noun=True), thing))
        if failed:
            delay_list.append(self.DIVIDER + "\n")
            delay_list.append(count_noun(len(failed), "item") + " had failures:\n")
            failed.sort()
            for thing, (f, t) in failed:
                delay_list.append(" %3d of %3d in %s\n"%(f, t, thing))
        if verbose:
            delay_list.append(count_noun(totalt, "test") + " in " + count_noun(len(self._name2ft), "item") + ".\n")
            delay_list.append("%s passed and %s failed.\n"%(totalt - totalf, totalf))
            if totalf:
                delay_list.append("***Test Failed***\n")
            else:
                delay_list.append("Test passed.\n")
        if print_now:
            delayed_output = "".join(delay_list)
            if self.options.logfile is not None:
                try:
                    with open(self.options.logfile, 'a') as logger:
                        logger.write(delayed_output)
                except IOError:
                    pass
            sys.stdout.write(delayed_output)
            sys.stdout.flush()
        return doctest.TestResults(totalf, totalt)

    def update_digests(self, example):
        """
        Update global and doctest digests.

        Sage's doctest runner tracks the state of doctests so that
        their dependencies are known.  For example, in the following
        two lines ::

            sage: R.<x> = ZZ[]
            sage: f = x^2 + 1

        it records that the second line depends on the first since the
        first INSERTS ``x`` into the global namespace and the second
        line RETRIEVES ``x`` from the global namespace.

        This function updates the hashes that record these
        dependencies.

        INPUT:

        - ``example`` -- a :class:`doctest.Example` instance

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os, hashlib
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: DTR.running_global_digest.hexdigest()
            'd41d8cd98f00b204e9800998ecf8427e'
            sage: DTR.running_doctest_digest = hashlib.md5()
            sage: ex = doctests[0].examples[0]; ex.predecessors = None
            sage: DTR.update_digests(ex)
            sage: DTR.running_global_digest.hexdigest()
            '3cb44104292c3a3ab4da3112ce5dc35c'
        """
        s = pre_hash(get_source(example))
        self.running_global_digest.update(s)
        self.running_doctest_digest.update(s)
        if example.predecessors is not None:
            digest = hashlib.md5(s)
            digest.update(reduce_hex(e.running_state for e in example.predecessors))
            example.running_state = digest.hexdigest()

    def execute(self, example, compiled, globs):
        """
        Runs the given example, recording dependencies.

        Rather than using a basic dictionary, Sage's doctest runner
        uses a :class:`sage.doctest.util.RecordingDict`, which records
        every time a value is set or retrieved.  Executing the given
        code with this recording dictionary as the namespace allows
        Sage to track dependencies between doctest lines.  For
        example, in the following two lines ::

            sage: R.<x> = ZZ[]
            sage: f = x^2 + 1

        the recording dictionary records that the second line depends
        on the first since the first INSERTS ``x`` into the global
        namespace and the second line RETRIEVES ``x`` from the global
        namespace.

        INPUT:

        - ``example`` -- a :class:`doctest.Example` instance.
        - ``compiled`` -- a code object produced by compiling ``example.source``
        - ``globs`` -- a dictionary in which to execute ``compiled``.

        OUTPUT:

        - the output of the compiled code snippet.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.util import RecordingDict
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os, hashlib
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: DTR.running_doctest_digest = hashlib.md5()
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: globs = RecordingDict(globals())
            sage: globs.has_key('doctest_var')
            False
            sage: doctests, extras = FDS.create_doctests(globs)
            sage: ex0 = doctests[0].examples[0]
            sage: compiled = compile(ex0.source, '<doctest sage.doctest.forker[0]>', 'single', 32768, 1)
            sage: DTR.execute(ex0, compiled, globs)
            1764
            sage: globs['doctest_var']
            42
            sage: globs.set
            set(['doctest_var'])
            sage: globs.got
            set(['Integer'])

        Now we can execute some more doctests to see the dependencies.

            sage: ex1 = doctests[0].examples[1]
            sage: compiled = compile(ex1.source, '<doctest sage.doctest.forker[1]>', 'single', 32768, 1)
            sage: DTR.execute(ex1, compiled, globs)
            sage: sorted(list(globs.set))
            ['R', 'a']
            sage: globs.got
            set(['ZZ'])
            sage: ex1.predecessors
            []

        ::

            sage: ex2 = doctests[0].examples[2]
            sage: compiled = compile(ex2.source, '<doctest sage.doctest.forker[2]>', 'single', 32768, 1)
            sage: DTR.execute(ex2, compiled, globs)
            a + 42
            sage: list(globs.set)
            []
            sage: sorted(list(globs.got))
            ['a', 'doctest_var']
            sage: set(ex2.predecessors) == set([ex0,ex1])
            True
        """
        if isinstance(globs, RecordingDict):
            globs.start()
        example.sequence_number = len(self.history)
        self.history.append(example)
        timer = Timer().start()
        try:
            exec compiled in globs
        finally:
            timer.stop().annotate(example)
            if isinstance(globs, RecordingDict):
                example.predecessors = []
                for name in globs.got:
                    ref = self.setters.get(name)
                    if ref is not None:
                        example.predecessors.append(ref)
                for name in globs.set:
                    self.setters[name] = example
            else:
                example.predecessors = None
            self.update_digests(example)
            example.total_state = self.running_global_digest.hexdigest()
            example.doctest_state = self.running_doctest_digest.hexdigest()

    def _failure_header(self, test, example):
        """
        We strip out ``sage:`` prompts, so we override
        :meth:`doctest.DocTestRunner._failure_header` for better
        reporting.

        INPUT:

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``.

        OUTPUT:

        - a string used for reporting that the given example failed.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: print DTR._failure_header(doctests[0], ex)
            **********************************************************************
            File ".../sage/doctest/forker.py", line 9, in sage.doctest.forker
            Failed example:
                doctest_var = 42; doctest_var^2
            <BLANKLINE>

       Without the source swapping::

            sage: import doctest
            sage: print doctest.DocTestRunner._failure_header(DTR, doctests[0], ex)
            **********************************************************************
            File ".../sage/doctest/forker.py", line 9, in sage.doctest.forker
            Failed example:
                doctest_var = Integer(42); doctest_var**Integer(2)
            <BLANKLINE>
        """
        with OriginalSource(example):
            return doctest.DocTestRunner._failure_header(self, test, example)

    def report_start(self, out, test, example):
        """
        Called when an example starts.

        INPUT:

        - ``out`` -- a function for printing (delayed or immediate)

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        OUTPUT:

        - prints a report to ``out``

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: DTR.report_start(sys.stdout.write, doctests[0], ex)
            Trying (line 9):    doctest_var = 42; doctest_var^2
            Expecting:
                1764
        """
        # We completely replace doctest.DocTestRunner.report_start so that we can include line numbers
        with OriginalSource(example):
            if self._verbose:
                start_txt = ('Trying (line %s):'%(test.lineno + example.lineno + 1)
                             + doctest._indent(example.source))
                if example.want:
                    start_txt += 'Expecting:\n' + doctest._indent(example.want)
                else:
                    start_txt += 'Expecting nothing\n'
                out(start_txt)
            if self._verbose or self.options.warn_long:
                self.tmp_time = walltime()

    def report_success(self, out, test, example, got):
        """
        Called when an example succeeds.

        INPUT:

        - ``out`` -- a function for printing (delayed or immediate)

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``got`` -- a string, the result of running ``example``

        OUTPUT:

        - prints a report to ``out``

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: from sage.misc.misc import walltime
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: DTR.tmp_time = walltime()
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: DTR.report_success(sys.stdout.write, doctests[0], ex, '1764')
            ok [...s]
        """
        # We completely replace doctest.DocTestRunner.report_success so that we can include time taken for the test
        if self._verbose:
            out("ok [%.2f s]\n"%(walltime(self.tmp_time)))

    # We name the local strangely since they get imported into
    def report_failure(self, out, test, example, got, globs):
        """
        Called when a doctest fails.

        INPUT:

        - ``out`` -- a function for printing (delayed or immediate)

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``got`` -- a string, the result of running ``example``

        - ``globs`` -- a dictionary of globals, used if in debugging mode

        OUTPUT:

        - prints a report to ``out``

        - if in debugging mode, starts an iPython prompt at the point
          of the failure

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: DTR.no_failure_yet = True
            sage: DTR.report_failure(sys.stdout.write, doctests[0], ex, 'BAD ANSWER\n', {})
            **********************************************************************
            File ".../sage/doctest/forker.py", line 9, in sage.doctest.forker
            Failed example:
                doctest_var = 42; doctest_var^2
            Expected:
                1764
            Got:
                BAD ANSWER
        """
        if not self.options.initial or self.no_failure_yet:
            self.no_failure_yet = False
            doctest.DocTestRunner.report_failure(self, out, test, example, got)

    def report_overtime(self, out, test, example, got):
        """
        Called when the ``warn_long`` option flag is set and a doctest
        runs longer than the specified time.

        INPUT:

        - ``out`` -- a function for printing (delayed or immediate)

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``got`` -- a string, the result of running ``example``

        OUTPUT:

        - prints a report to ``out``

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: from sage.misc.misc import walltime
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: DTR.tmp_time = walltime() - 1
            sage: DTR.report_overtime(sys.stdout.write, doctests[0], ex, 'BAD ANSWER\n')
            **********************************************************************
            File ".../sage/doctest/forker.py", line 9, in sage.doctest.forker
            Failed example:
                doctest_var = 42; doctest_var^2
            Test ran for ... s
        """
        out(self._failure_header(test, example) +
            "Test ran for %.2f s\n"%(walltime(self.tmp_time)))

    def report_unexpected_exception(self, out, test, example, exc_info):
        """
        Called when a doctest raises an exception that's not matched by the expected output.

        If debugging has been turned on, starts an interactive debugger.

        INPUT:

        - ``out`` -- a function for printing (delayed or immediate)

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``exc_info`` -- the result of ``sys.exc_info()``

        OUTPUT:

        - prints a report to ``out``

        - if in debugging mode, starts PDB with the given traceback

        EXAMPLES::

            sage: import os
            sage: os.environ['SAGE_PEXPECT_LOG'] = "1"
            sage: sage0.quit()
            sage: _ = sage0.eval("import doctest, sys, os, multiprocessing, subprocess")
            sage: _ = sage0.eval("from sage.doctest.parsing import SageOutputChecker")
            sage: _ = sage0.eval("import sage.doctest.forker as sdf")
            sage: _ = sage0.eval("from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()")
            sage: _ = sage0.eval("ex = doctest.Example('E = EllipticCurve([0,0]); E', 'A singular Elliptic Curve')")
            sage: _ = sage0.eval("DT = doctest.DocTest([ex], globals(), 'singular_curve', None, 0, None)")
            sage: _ = sage0.eval("DTR = sdf.SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)")
            sage: _ = sage0.eval("sdf.debug_semaphore = multiprocessing.Semaphore()")
            sage: _ = sage0.eval("sdf.output_semaphore = multiprocessing.Semaphore()")
            sage: sage0._prompt = r"\(Pdb\) "
            sage: sage0.eval("DTR.run(DT, clear_globs=False)") # indirect doctest
            '... "Invariants %s define a singular curve."%ainvs'
            sage: sage0.eval("l")
            '...if self.discriminant() == 0:...raise ArithmeticError...'
            sage: sage0.eval("u")
            '...EllipticCurve_field.__init__(self, [field(x) for x in ainvs])'
            sage: sage0.eval("p ainvs")
            '[0, 0]'
            sage: sage0._prompt = "sage: "
            sage: sage0.eval("quit")
            'TestResults(failed=1, attempted=1)'
        """
        if debug_semaphore is not None or not self.options.initial or self.no_failure_yet:
            self.no_failure_yet = False
            returnval = doctest.DocTestRunner.report_unexpected_exception(self, out, test, example, exc_info)
            if debug_semaphore is not None:
                debug_semaphore.acquire()
                output_semaphore.acquire()
                self._fakeout.stop_spoofing()
                self._fakeout.activate_stdin()
                try:
                    exc_type, exc_val, exc_tb = exc_info
                    self.debugger.reset()
                    self.debugger.interaction(None, exc_tb)
                finally:
                    self._fakeout.start_spoofing()
                    self._fakeout.deactivate_stdin()
                    sys.stdout.flush()
                    output_semaphore.release()
                    debug_semaphore.release()
            return returnval

    def update_results(self, D):
        """
        When returning results we pick out the results of interest
        since many attributes are not pickleable.

        INPUT:

        - ``D`` -- a dictionary to update with cputime and walltime

        OUTPUT:

        - the number of failures (or False if there is no failure attribute)

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource, DictAsObject
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','forker.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: from sage.doctest.util import Timer
            sage: T = Timer().start()
            sage: DTR.run(doctests[0])
            TestResults(failed=0, attempted=4)
            sage: T.stop().annotate(DTR)
            sage: D = DictAsObject({'cputime':[],'walltime':[],'err':None})
            sage: DTR.update_results(D)
            0
            sage: sorted(list(D.iteritems()))
            [('cputime', [...]), ('err', None), ('failures', 0), ('walltime', [...])]
        """
        for key in ["cputime","walltime"]:
            if not D.has_key(key):
                D[key] = []
            if hasattr(self, key):
                D[key].append(self.__dict__[key])
        if hasattr(self, 'failures'):
            D['failures'] = self.failures
            return self.failures
        else:
            return False

class DocTestDispatcher(SageObject):
    """
    Creates and dispatches doctesting tasks to workers in series or parallel.
    """
    def __init__(self, controller):
        """
        INPUT:

        - ``controller`` -- a :class:`sage.doctest.control.DocTestController` instance

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: DocTestDispatcher(DocTestController(DocTestDefaults(), []))
            <class 'sage.doctest.forker.DocTestDispatcher'>
        """
        self.controller = controller

    def _serial_dispatch(self):
        """
        Run the doctests from the controller's specified sources in series.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.util import Timer
            sage: import os
            sage: powser = os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage', 'sage', 'rings', 'homset.py')
            sage: inf = os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage', 'sage', 'rings', 'ideal.py')
            sage: DC = DocTestController(DocTestDefaults(), [powser, inf])
            sage: DC.expand_files_into_sources()
            sage: DD = DocTestDispatcher(DC)
            sage: DR = DocTestReporter(DC)
            sage: DC.reporter = DR
            sage: DC.dispatcher = DD
            sage: DC.timer = Timer().start()
            sage: DD._serial_dispatch()
            sage -t .../rings/homset.py
                [... tests, ... s]
            sage -t .../rings/ideal.py
                [... tests, ... s]
        """
        opt = self.controller.options
        sources = self.controller.sources
        try:
            for source in sources:
                output = os.tmpfile()
                result = DocTestTask(source)(None, output, opt)
                output.seek(0)
                output = output.read()
                self.controller.reporter.report(source, False, 0, result, output)
        except KeyboardInterrupt:
            print "Tesing %s interrupted!"%(source.basename)
            sys.stdout.flush()
            # We rely on the reporter to give more information back to the user.

    def _parallel_dispatch(self):
        """
        Run the doctests from the controller's specified sources in series.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.util import Timer
            sage: import os
            sage: crem = os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage', 'sage', 'databases', 'cremona.py')
            sage: bigo = os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage', 'sage', 'rings', 'big_oh.py')
            sage: DC = DocTestController(DocTestDefaults(), [crem, bigo])
            sage: DC.expand_files_into_sources()
            sage: DD = DocTestDispatcher(DC)
            sage: DR = DocTestReporter(DC)
            sage: DC.reporter = DR
            sage: DC.dispatcher = DD
            sage: DC.timer = Timer().start()
            sage: DD._parallel_dispatch()
            sage -t .../databases/cremona.py
                [... tests, ... s]
            sage -t .../rings/big_oh.py
                [... tests, ... s]
        """
        opt = self.controller.options
        sources = self.controller.sources
        nthreads = opt.nthreads
        tasks = multiprocessing.JoinableQueue(len(sources)+nthreads)
        for source in sources:
            tasks.put(DocTestTask(source))
        results = multiprocessing.Queue(len(sources))
        workers = [DocTestWorker(tasks, results, options=opt) for i in xrange(nthreads)]
        for w in workers:
            w.start()
            tasks.put(None)
        done_count = 0
        while True:
            try:
                result = results.get()
                if result is None:
                    done_count += 1
                    if done_count >= len(workers):
                        break
                else:
                    try:
                        if output_semaphore is not None:
                            output_semaphore.acquire()
                        self.controller.reporter.report(*result)
                    finally:
                        sys.stdout.flush()
                        if output_semaphore is not None:
                            output_semaphore.release()
            except Empty:
                time.sleep(.05)

    def dispatch(self): # todo, nthreads=options.nthreads, streaming=False, verbose=options.verbose, debug=options.debug, run_id=run_id
        """
        Runs the doctests for the controller's specified sources.

        It will run in series or parallel depending on the ``serial``
        option on the controller's option object.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.util import Timer
            sage: import os
            sage: freehom = os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage', 'sage', 'modules', 'free_module_homspace.py')
            sage: bigo = os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage', 'sage', 'rings', 'big_oh.py')
            sage: DC = DocTestController(DocTestDefaults(), [freehom, bigo])
            sage: DC.expand_files_into_sources()
            sage: DD = DocTestDispatcher(DC)
            sage: DR = DocTestReporter(DC)
            sage: DC.reporter = DR
            sage: DC.dispatcher = DD
            sage: DC.timer = Timer().start()
            sage: DD.dispatch()
            sage -t .../sage/modules/free_module_homspace.py
                [... tests, ... s]
            sage -t .../sage/rings/big_oh.py
                [... tests, ... s]
        """
        init_sage()
        global debug_semaphore, output_semaphore
        try:
            output_semaphore = multiprocessing.Semaphore()
            if self.controller.options.debug:
                debug_semaphore = multiprocessing.Semaphore()
            if self.controller.options.serial:
                self._serial_dispatch()
            else:
                self._parallel_dispatch()
        finally:
            output_semaphore = None
            debug_semaphore = None

class DocTestWorker(multiprocessing.Process):
    """
    The DocTestWorker process pulls tasks from its task queue and runs
    them, putting the results onto its result queue.

    It can recover gracefully from timeouts and segfaults in the
    tasks.

    INPUT:

    - ``task_queue`` -- a :class:multiprocessing.Queue instance

    - ``result_queue`` -- a :class:multiprocessing.Queue instance

    - ``timeout`` -- a positive number

    - ``verbose`` -- a boolean, controls verbosity when killing processes.

    EXAMPLES::

        sage: import multiprocessing, os, time
        sage: from sage.doctest.forker import DocTestWorker, DocTestTask
        sage: from sage.doctest.sources import FileDocTestSource
        sage: from sage.doctest.reporting import DocTestReporter
        sage: from sage.doctest.control import DocTestController, DocTestDefaults
        sage: tasks = multiprocessing.JoinableQueue(2)
        sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','sources.py')
        sage: DD = DocTestDefaults()
        sage: DC = DocTestController(DD, filename)
        sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
        sage: DTT = DocTestTask(FDS)
        sage: tasks.put(DTT)
        sage: tasks.put(None)
        sage: results = multiprocessing.Queue(1)
        sage: W = DocTestWorker(tasks, results, DD)
        sage: from Queue import Empty
        sage: W.start()
        sage: reporter = DocTestReporter(DC)
        sage: while True:
        ....:     try:
        ....:         result = results.get()
        ....:         if result is None: break
        ....:         reporter.report(*result)
        ....:     except Empty:
        ....:         time.sleep(0.05)
        ....:
        sage -t .../sage/doctest/sources.py
            [... tests, ... s]
    """
    def __init__(self, task_queue, result_queue, options):
        """
        Initialization.

        TESTS::

            sage: run_doctests(sage.rings.big_oh) # indirect doctest
            Doctesting .../sage/rings/big_oh.py
            Running doctests with ID ...
            Doctesting 1 file.
            sage -t .../sage/rings/big_oh.py
                [... tests, ... s]
            ------------------------------------------------------------------------
            All tests passed!
            ------------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
        """
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.options = options

    def run(self):
        """
        Iterates through the tasks in the task queue and runs the
        doctests in the corresponding sources.

        TESTS::

            sage: run_doctests(sage.symbolic.units) # indirect doctest
            Doctesting .../sage/symbolic/units.py
            Running doctests with ID ...
            Doctesting 1 file.
            sage -t .../sage/symbolic/units.py
                [... tests, ... s]
            ------------------------------------------------------------------------
            All tests passed!
            ------------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
        """
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                self.task_queue.task_done()
                self.result_queue.put(None)
                break
            result_pipe = multiprocessing.Queue() # pipes block
            output_file = os.tmpfile()
            p = multiprocessing.Process(target = next_task, args=(result_pipe,output_file,self.options))
            try:
                p.start()
                ## We would like to do the following:
                # results = result_pipe.get(timeout=self.options.timeout)
                # p.join(timeout=self.options.timeout)
                ## But this hangs when the underlying process segfaults.
                ## So instead we query the underlying process manually.
                maxwait = 0.1
                curwait = 0.004
                deathwait = 1
                results = None
                t = walltime()
                alive = True
                while True:
                    try:
                        results = result_pipe.get(timeout=curwait)
                    except Empty:
                        curwait *= 2
                        if curwait > maxwait:
                            curwait = maxwait
                    if results is not None or walltime(t) >= self.options.timeout:
                        break
                    if not p.is_alive():
                        if not alive:
                            # We've already waited, but results is still None
                            results = 0, DictAsObject(dict(err='noresult')), ""
                            break
                        alive = False
                        if results is None:
                            curwait = deathwait
                answer = next_task.source, walltime(t) >= self.options.timeout, p.exitcode, results
            ## This except block was used in the approch mentioned above
            #except Empty:
            #    answer = next_task.source, True, -1, None
            except KeyboardInterrupt:
                if output_semaphore is not None:
                    output_semaphore.acquire()
                try:
                    print "Testing %s interrupted!"%(next_task.source.basename)
                finally:
                    sys.stdout.flush()
                    if output_semaphore is not None:
                        output_semaphore.release()
                # We quietly die: the child process is killed in the finally block.
                break
            finally:
                output_file.seek(0)
                output = output_file.read()
                self.annihilate(p, False)
            self.task_queue.task_done()
            answer = answer + (output,)
            self.result_queue.put(answer)

    def annihilate(self, p, verbose):
        """
        Kills the process p and all of its descendents.

        INPUT:

        - ``p`` -- a process to kill

        - ``verbose`` -- boolean

        EXAMPLES::

            sage: from sage.doctest.forker import DocTestWorker
            sage: DTW = DocTestWorker(None, None, None)
            sage: from multiprocessing import Process
            sage: import time
            sage: def make_babies():
            ....:     for i in range(3):
            ....:         p = Process(target = time.sleep, args = (4,))
            ....:         p.start()
            ....:     time.sleep(6)
            ....:
            sage: p = Process(target = make_babies)
            sage: p.start()
            sage: p.is_alive()
            True
            sage: DTW.annihilate(p, True)
            killing ...
            subprocesses [...]
            sage: p.is_alive()
            False
        """
        ## The proper way to do this is to create a process group.
        ## But then control-Z doesn't get forwarded automatically....
        #        os.system("ps -fTj")
        #        os.system("./pstree -p %s" % os.getpid())
        #        os.system("kill -9 -%s" % p.pid)
        #        os.system("kill -9 %s" % p.pid)
        #        os.killpg(p.pid, signal.SIGSTOP)
        def kill_all(pids, sig):
            for pid in pids:
                try:
                    os.kill(pid, sig)
                except OSError:
                    pass
        def all_subprocesses(root_pid):
            listing = subprocess.Popen(['ps', '-o', 'pid,ppid'], stdout=subprocess.PIPE).communicate()[0].split('\n')[1:]
            children = defaultdict(list)
            for line in listing:
                line = line.strip()
                if line:
                    pid, ppid = [int(s.strip()) for s in line.split(' ') if s]
                    children[ppid].append(pid)
            all = []
            def add_children(ppid):
                for pid in children[ppid]:
                    all.append(pid)
                    add_children(pid)
            add_children(root_pid)
            return all
        if p.is_alive() and verbose:
            print "killing", p.pid
        if p.is_alive():
            subprocesses = all_subprocesses(p.pid)
            if verbose:
                print "subprocesses", subprocesses
            p.terminate()
            p.join(1)
            if subprocesses:
                kill_all(reversed(subprocesses), signal.SIGTERM)
                time.sleep(.1)
                kill_all(reversed(subprocesses), signal.SIGKILL)
                p.join(1)
        if p.is_alive():
            os.kill(p.pid, signal.SIGKILL)
            p.join(1)

class DocTestTask(object):
    """
    This class encapsulates the tests from a single source.  It can be
    called in series (via :meth:`DocTestDispatcher._serial_dispatch`)
    or in parallel through :class:`DocTestWorker` instances.

    This class does not insulate from problems in the source
    (e.g. entering an infinite loop or causing a segfault).

    INPUT:

    - ``source`` -- a :class:`sage.doctest.sources.DocTestSource` instance.

    - ``verbose`` -- boolean, controls reporting of progress by :class:`doctest.DocTestRunner`.

    EXAMPLES::

        sage: from sage.doctest.forker import DocTestTask
        sage: from sage.doctest.sources import FileDocTestSource
        sage: from sage.doctest.control import DocTestDefaults, DocTestController
        sage: import os
        sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','sources.py')
        sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
        sage: DTT = DocTestTask(FDS)
        sage: DD = DocTestDefaults()
        sage: DC = DocTestController(DD,[filename])
        sage: ntests, results, delayed_output = DTT(output = os.tmpfile(),options=DD); ntests
        276
        sage: sorted(results.keys())
        ['cputime', 'err', 'failures', 'walltime']
    """
    def __init__(self, source):
        """
        Initialization.

        TESTS::

            sage: from sage.doctest.forker import DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: import os
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','sources.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: DocTestTask(FDS)
            <sage.doctest.forker.DocTestTask object at ...>
        """
        self.source = source

    def __call__(self, result_pipe=None, output=None, options=None):
        """
        Calling the task does the actual work of running the doctests.

        INPUT:

        - ``result_pipe`` -- either None (if run in series) or a :class:`multiprocessing.Queue`

        - ``output`` -- a temporary file that's used by the doctest runner to redirect stdout.

        OUPUT:

        - If ``result_pipe`` is None, returns a pair ``(doctests,
          runner)``, where ``doctests`` is a list of
          :class:`doctest.DocTest` instances and ``runner`` is an
          annotated ``SageDocTestRunner`` instance.

        - If ``result_pipe`` is not None, puts ``(doctests, runner)``
          onto the result pipe and returns nothing.

        EXAMPLES::

            sage: from sage.doctest.forker import DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: import os
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','parsing.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: DTT = DocTestTask(FDS)
            sage: DD = DocTestDefaults()
            sage: DC = DocTestController(DD, [filename])
            sage: ntests, runner, delayed_output = DTT(None, os.tmpfile(), DD)
            sage: runner.failures
            0
            sage: ntests
            192
        """
        result = None
        try:
            file = self.source.path
            basename = self.source.basename
            import sage.all_cmdline
            runner = SageDocTestRunner(SageOutputChecker(), verbose=options.verbose, output=output, sage_options=options, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            runner.basename = basename
            N = options.file_iterations
            results = DictAsObject(dict(walltime=[],cputime=[],err=None))
            if options.verbose:
                # This will cause summaries to print immediately
                delayed_summaries = None
            else:
                delayed_summaries = []
            for it in range(N):
                sage_namespace = RecordingDict(dict(sage.all_cmdline.__dict__))
                sage_namespace['__name__'] = '__main__'
                sage_namespace['__package__'] = None
                doctests, extras = self.source.create_doctests(sage_namespace)
                timer = Timer().start()

                for test in doctests:
                    runner.run(test)
                runner.filename = file
                if output_semaphore is not None:
                    output_semaphore.acquire()
                failed, tried = runner.summarize(verbose=options.verbose, delay_list=delayed_summaries)
                if output_semaphore is not None:
                    sys.stdout.flush()
                    output_semaphore.release()
                timer.stop().annotate(runner)
                if runner.update_results(results):
                    break
            if extras['tab']:
                results.err = 'tab'
                results.tab_linenos = extras['tab']
            # We subtract 1 to remove the sig_on_count() tests
            if delayed_summaries is None:
                delayed_summaries = ""
            else:
                delayed_summaries = "".join(delayed_summaries)
            result = sum([max(0,len(test.examples) - 1) for test in doctests]), results, delayed_summaries
        except KeyboardInterrupt:
            if result_pipe is None:
                # We're in serial mode
                raise
            # Otherwise, we've been started by a worker: we'll return the result in the finally block
            result = 0, DictAsObject(dict(err='ctlC')), ""
        except IOError:
            # File doesn't exist
            result = 0, DictAsObject(dict(err='file')), ""
        except Exception:
            exc_info = sys.exc_info()
            tb = "".join(traceback.format_exception(*exc_info))
            result = 0, DictAsObject(dict(err=exc_info[0], tb=tb)), ""
        finally:
            if result_pipe is None:
                return result
            else:
                result_pipe.put(result, False)
