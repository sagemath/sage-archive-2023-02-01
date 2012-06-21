"""
This module determines how doctest results are reported to the user.

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

import sys
from sage.structure.sage_object import SageObject
from sage.doctest.util import count_noun

class DocTestReporter(SageObject):
    """
    This class reports to the users on the results of doctests.
    """
    def __init__(self, controller):
        """
        Initialize the reporter.

        INPUT:

        - ``controller`` -- a
          :class:`sage.doctest.control.DocTestController` instance.
          Note that some methods assume that appropriate tests have
          been run by the controller.

        EXAMPLES::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: import os
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','reporting.py')
            sage: DC = DocTestController(DocTestDefaults(),[filename])
            sage: DTR = DocTestReporter(DC)
        """
        self.controller = controller
        self.postscript = dict(lines=[], cputime=0, walltime=0)
        self.sources_completed = 0
        self.stats = {}
        self.error_status = 0 # partially for backward compatibility

    def report(self, source, timeout, return_code, results, output):
        """
        Report on the result of running doctests on a given source.

        INPUT:

        - ``source`` -- a source from :mod:`sage.doctest.sources`

        - ``timeout`` -- a boolean, whether doctests timed out

        - ``return_code`` -- an int, the return code of the process
          running doctests on that file.

        - ``results`` -- (irrelevant if ``timeout`` or
          ``return_code``), a tuple

          - ``ntests`` -- the number of doctests

          - ``timings`` -- a
            :class:`sage.doctest.sources.DictAsObject` instance
            storing timing data.

          - ``delayed_output`` -- delayed output to be printed now

        - ``output`` -- a string, printed if there was some kind of
          failure

        EXAMPLES::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource, DictAsObject
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.util import Timer
            sage: import os, sys, doctest
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','reporting.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: DD = DocTestDefaults()
            sage: DC = DocTestController(DD,[filename])
            sage: DTR = DocTestReporter(DC)

        You can report a timeout::

            sage: DTR.report(FDS, True, 0, None, "Output so far...")
            sage -t .../sage/doctest/reporting.py
                Timed out!
            ********************************************************************************
            Tests run before process froze:
            Output so far...
            ********************************************************************************
            sage: DTR.stats
            {'sage.doctest.reporting': {'failed': True, 'walltime': 1000000.0}}

        Or a process that returned a bad exit code::

            sage: DTR.report(FDS, False, 3, None, "Output before bad exit")
            sage -t .../sage/doctest/reporting.py
                Bad exit: 3
            ********************************************************************************
            Tests run before process failed:
            Output before bad exit
            ********************************************************************************
            sage: DTR.stats
            {'sage.doctest.reporting': {'failed': True, 'walltime': 1000000.0}}

        Or tell the user that everything succeeded::

            sage: doctests, extras = FDS.create_doctests(globals())
            sage: runner = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: Timer().start().stop().annotate(runner)
            sage: D = DictAsObject({'err':None})
            sage: runner.update_results(D)
            0
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D, ""), "Good tests")
            sage -t .../doctest/reporting.py
                [... tests, 0.0 s]
            sage: DTR.stats
            {'sage.doctest.reporting': {'walltime': ...}}

        Or inform the user that some doctests failed::

            sage: runner.failures = 1
            sage: runner.update_results(D)
            1
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D, ""), "Doctest output including the failure...")
            sage -t .../doctest/reporting.py
                [... tests, 1 failure, 0.0 s]
        """
        log = self.controller.log
        postscript = self.postscript
        stats = self.stats
        try:
            basename = source.basename
            if self.controller.options.long:
                islong = "--long "
            else:
                islong = ""
            warnlong = self.controller.options.warn_long
            if warnlong is None:
                warnlong = ""
            elif warnlong == 1.0:
                warnlong = "--warn-long "
            else:
                warnlong = "--warn-long %.1f "%(warnlong)
            cmd = "sage -t %s%s%s"%(islong, warnlong, source.printpath)
            log(cmd)
            if timeout:
                log("    Timed out!\n%s\nTests run before process froze:"%("*"*80))
                log(output)
                log("*"*80)
                postscript['lines'].append(cmd + " # Time out")
                stats[basename] = dict(failed=True, walltime=1e6)
                self.error_status |= 64
            elif return_code:
                log("    Bad exit: %s\n%s\nTests run before process failed:"%(return_code,"*"*80))
                log(output)
                log("*"*80)
                postscript['lines'].append(cmd + " # Bad exit: %s" % return_code)
                stats[basename] = dict(failed=True, walltime=1e6)
                self.error_status |= 4
            elif results is None:
                log("    Error in doctesting framework!\n%s\nTests run before error:"%("*"*80))
                log(output)
                log("*"*80)
                stats[basename] = dict(failed=True, walltime=1e6)
                postscript['lines'].append(cmd + " # Unhandled doctest exception")
                self.error_status |= 8
            else:
                ntests, result_dict, delayed_output = results
                if hasattr(result_dict, 'walltime') and hasattr(result_dict.walltime, '__len__') and len(result_dict.walltime) > 0:
                    wall = sum(result_dict.walltime) / len(result_dict.walltime)
                else:
                    wall = 1e6
                if hasattr(result_dict, 'cputime') and hasattr(result_dict.cputime, '__len__') and len(result_dict.cputime) > 0:
                    cpu = sum(result_dict.cputime) / len(result_dict.cputime)
                else:
                    cpu = 1e6
                if result_dict.err == 'file':
                    log("    File not found!")
                    postscript['lines'].append(cmd + " # File not found")
                    self.error_status |= 1
                elif result_dict.err == 'noresult':
                    log("    Error in doctesting framework (no result returned)\n%s\nTests run before error:"%("*"*80))
                    log(output)
                    log("*"*80)
                    postscript['lines'].append(cmd + " # Testing error: no result returned")
                    self.error_status |= 8
                elif result_dict.err == 'ctlC':
                    log("    Unexpected KeyboardInterrupt raised in file\n%s\nTests run before interrupt:"%("*"*80))
                    log(output)
                    log("*"*80)
                    postscript['lines'].append(cmd + " # Unhandled KeyboardInterrupt")
                    stats[basename] = dict(failed=True, walltime = wall)
                    self.error_status |= 8
                elif result_dict.err == 'tab':
                    if len(result_dict.tab_linenos) > 5:
                        result_dict.tab_linenos[3:-1] = "..."
                    tabs = " " + ",".join(result_dict.tab_linenos)
                    if len(result_dict.tab_linenos) > 1:
                        tabs = "s" + tabs
                    log("    Error: TAB character found at line%s"%(tabs))
                    postscript['lines'].append(cmd + " # Tab character found")
                    self.error_status |= 32
                elif result_dict.err is True:
                    # This case should not occur
                    log("    Error in doctesting framework")
                    postscript['lines'].append(cmd + " # Unhandled doctest exception")
                    self.error_status |= 8
                elif result_dict.err is not None:
                    if hasattr(result_dict, 'tb'):
                        log(result_dict.tb)
                    if hasattr(result_dict.err, '__name__'):
                        err = result_dict.err.__name__
                    else:
                        err = repr(result_dict.err)
                    postscript['lines'].append(cmd + " # %s in loading"%(err))
                    if hasattr(result_dict, 'walltime'):
                        stats[basename] = dict(failed=True, walltime=wall)
                    else:
                        stats[basename] = dict(failed=True, walltime=1e6)
                    self.error_status |= 16
                if result_dict.err is None or result_dict.err == 'tab':
                    f = result_dict.failures
                    if f:
                        postscript['lines'].append(cmd + " # %s failed" % (count_noun(f, "doctest")))
                        self.error_status |= 128
                    if f or result_dict.err == 'tab':
                        stats[basename] = dict(failed=True, walltime=wall)
                    else:
                        stats[basename] = dict(walltime=wall)
                    postscript['cputime'] += cpu
                    postscript['walltime'] += wall

                    if delayed_output:
                        if delayed_output[-1] == '\n':
                            delayed_output = delayed_output[:-1]
                        log(delayed_output)
                    log("    [%s, %s%.1f s]" % (count_noun(ntests, "test"), "%s, "%(count_noun(f, "failure")) if f else "", wall))
            self.sources_completed += 1

        except StandardError:
            import traceback
            traceback.print_exc()

    def finalize(self):
        """
        Print out the postcript that summarizes the doctests that were run.

        EXAMPLES:

        First we have to set up a bunch of stuff::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource, DictAsObject
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.util import Timer
            sage: import os, sys, doctest
            sage: filename = os.path.join(os.environ['SAGE_ROOT'],'devel','sage','sage','doctest','reporting.py')
            sage: FDS = FileDocTestSource(filename,True,False,set(['sage']),None)
            sage: DD = DocTestDefaults()
            sage: DC = DocTestController(DD,[filename])
            sage: DTR = DocTestReporter(DC)

        Now we pretend to run some doctests::

            sage: DTR.report(FDS, True, 0, None, "Output so far...")
            sage -t .../sage/doctest/reporting.py
                Timed out!
            ********************************************************************************
            Tests run before process froze:
            Output so far...
            ********************************************************************************
            sage: DTR.report(FDS, False, 3, None, "Output before bad exit")
            sage -t .../sage/doctest/reporting.py
                Bad exit: 3
            ********************************************************************************
            Tests run before process failed:
            Output before bad exit
            ********************************************************************************
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: runner = SageDocTestRunner(SageOutputChecker(), output = os.tmpfile(), verbose=False, sage_options=DD,optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: t = Timer().start().stop()
            sage: t.annotate(runner)
            sage: DC.timer = t
            sage: D = DictAsObject({'err':None})
            sage: runner.update_results(D)
            0
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D, ""), "Good tests")
            sage -t .../doctest/reporting.py
                [... tests, 0.0 s]
            sage: runner.failures = 1
            sage: runner.update_results(D)
            1
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D, ""), "Doctest output including the failure...")
            sage -t .../doctest/reporting.py
                [... tests, 1 failure, 0.0 s]

        Now we can show the output of finalize::

            sage: DC.sources = [None] * 4 # to fool the finalize method
            sage: DTR.finalize()
            ------------------------------------------------------------------------
            sage -t .../sage/doctest/reporting.py # Time out
            sage -t .../sage/doctest/reporting.py # Bad exit: 3
            sage -t .../sage/doctest/reporting.py # 1 doctest failed
            ------------------------------------------------------------------------
            Total time for all tests: 0.0 seconds
                cpu time: 0.0 seconds
                cumulative wall time: 0.0 seconds

        If we interrupted doctests, then the number of files tested
        will not match the number of sources on the controller::

            sage: DC.sources = [None] * 6
            sage: DTR.finalize()
            <BLANKLINE>
            ------------------------------------------------------------------------
            sage -t .../sage/doctest/reporting.py # Time out
            sage -t .../sage/doctest/reporting.py # Bad exit: 3
            sage -t .../sage/doctest/reporting.py # 1 doctest failed
            Doctests interrupted: 4/6 files tested
            ------------------------------------------------------------------------
            Total time for all tests: 0.0 seconds
                cpu time: 0.0 seconds
                cumulative wall time: 0.0 seconds
        """
        log = self.controller.log
        postscript = self.postscript
        if self.sources_completed < len(self.controller.sources) * self.controller.options.global_iterations:
            postscript['lines'].append("Doctests interrupted: %s/%s files tested"%(self.sources_completed, len(self.controller.sources)))
            self.error_status |= 2
        elif not postscript['lines']:
            postscript['lines'].append("All tests passed!")
        log('-' * 72)
        log("\n".join(postscript['lines']))
        log('-' * 72)
        log("Total time for all tests: %.1f seconds" % self.controller.timer.walltime)
        log("    cpu time: %.1f seconds" % postscript['cputime'])
        log("    cumulative wall time: %.1f seconds" % postscript['walltime'])
        sys.stdout.flush()
