"""
Various tests for the doctesting framework.

Many tests (with expected failures or crashes) are run in a
subprocess, those tests can be found in the ``tests/`` subdirectory.

EXAMPLES::

    sage: import signal
    sage: import subprocess
    sage: import time
    sage: from sage.env import SAGE_SRC
    sage: tests_dir = os.path.join(SAGE_SRC, 'sage', 'doctest', 'tests')

Check that :trac:`2235` has been fixed::

    sage: subprocess.call(["sage", "-t", "longtime.rst"], cwd=tests_dir)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t longtime.rst
    [0 tests, ...s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0
    sage: subprocess.call(["sage", "-t", "-l", "longtime.rst"], cwd=tests_dir)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t --long longtime.rst
    [1 test, ...s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0

Test the ``--initial`` option::

    sage: subprocess.call(["sage", "-t", "-i", "initial.rst"], cwd=tests_dir)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t initial.rst
    **********************************************************************
    File "initial.rst", line 4, in sage.doctest.tests.initial
    Failed example:
        a = binomiak(10,5)  # random to test that we still get the exception
    Exception raised:
        Traceback (most recent call last):
        ...
        NameError: name 'binomiak' is not defined
    **********************************************************************
    File "initial.rst", line 14, in sage.doctest.tests.initial
    Failed example:
        binomial(10,5)
    Expected:
        255
    Got:
        252
    **********************************************************************
    ...
    ----------------------------------------------------------------------
    sage -t initial.rst  # 5 doctests failed
    ----------------------------------------------------------------------
    ...
    1

Test a timeout using the ``SAGE_TIMEOUT`` environment variable::

    sage: env = dict(os.environ)
    sage: env['SAGE_TIMEOUT'] = "3"
    sage: subprocess.call(["sage", "-t", "99seconds.rst"], cwd=tests_dir, env=env)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t 99seconds.rst
        Time out
    **********************************************************************
    Tests run before process timed out:
    ...
    ----------------------------------------------------------------------
    sage -t 99seconds.rst  # Time out
    ----------------------------------------------------------------------
    ...
    4

Test handling of ``KeyboardInterrupt``s in doctests::

    sage: subprocess.call(["sage", "-t", "keyboardinterrupt.rst"], cwd=tests_dir)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t keyboardinterrupt.rst
    **********************************************************************
    File "keyboardinterrupt.rst", line 11, in sage.doctest.tests.keyboardinterrupt
    Failed example:
        raise KeyboardInterrupt
    Exception raised:
        Traceback (most recent call last):
        ...
        KeyboardInterrupt
    **********************************************************************
    ...
    ----------------------------------------------------------------------
    sage -t keyboardinterrupt.rst  # 1 doctest failed
    ----------------------------------------------------------------------
    ...
    1

Interrupt the doctester::

    sage: subprocess.call(["sage", "-t", "interrupt.rst"], cwd=tests_dir)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t interrupt.rst
    Killing test interrupt.rst
    ----------------------------------------------------------------------
    Doctests interrupted: 0/1 files tested
    ----------------------------------------------------------------------
    ...
    128

Interrupt the doctester (while parallel testing) when a doctest cannot
be interrupted. We also test that passing a ridiculous number of threads
doesn't hurt::

    sage: F = tmp_filename()
    sage: env = dict(os.environ)
    sage: env['DOCTEST_TEST_PID_FILE'] = F  # Doctester will write its PID in this file
    sage: subprocess.call(["sage", "-tp", "1000000", "--timeout=120",  # long time
    ....:     "99seconds.rst", "interrupt_diehard.rst"], cwd=tests_dir, env=env)
    Running doctests...
    Doctesting 2 files using 1000000 threads.
    Killing test 99seconds.rst
    Killing test interrupt_diehard.rst
    ----------------------------------------------------------------------
    Doctests interrupted: 0/2 files tested
    ----------------------------------------------------------------------
    ...
    128

Even though the doctester master process has exited, the child process
is still alive, but it should be killed automatically
in 120 * 0.05 = 6 seconds::

    sage: pid = int(open(F).read())    # long time
    sage: time.sleep(2)                # long time
    sage: os.kill(pid, signal.SIGHUP)  # long time; 2 seconds passed => still alive
    sage: time.sleep(6)                # long time
    sage: os.kill(pid, signal.SIGHUP)  # long time; 8 seconds passed => dead
    Traceback (most recent call last):
    ...
    OSError: ...

Test a doctest failing with ``abort()``::

    sage: subprocess.call(["sage", "-t", "abort.rst"], cwd=tests_dir)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t abort.rst
        Killed due to abort
    **********************************************************************
    Tests run before process failed:
    ...
    ------------------------------------------------------------------------
    Unhandled SIGABRT: An abort() occurred in Sage.
    This probably occurred because a *compiled* component of Sage has a bug
    in it and is not properly wrapped with sig_on(), sig_off(). You might
    want to run Sage under gdb with 'sage -gdb' to debug this.
    Sage will now terminate.
    ------------------------------------------------------------------------
    ...
    ----------------------------------------------------------------------
    sage -t abort.rst  # Killed due to abort
    ----------------------------------------------------------------------
    ...
    16

A different kind of crash::

    sage: subprocess.call(["sage", "-t", "fail_and_die.rst"], cwd=tests_dir)  # long time
    Running doctests...
    Doctesting 1 file.
    sage -t fail_and_die.rst
    **********************************************************************
    File "fail_and_die.rst", line 5, in sage.doctest.tests.fail_and_die
    Failed example:
        this_gives_a_NameError
    Exception raised:
        Traceback (most recent call last):
        ...
        NameError: name 'this_gives_a_NameError' is not defined
        Killed due to kill signal
    **********************************************************************
    Tests run before process failed:
    ...
    ----------------------------------------------------------------------
    sage -t fail_and_die.rst  # Killed due to kill signal
    ----------------------------------------------------------------------
    ...
    16

Test running under gdb, without and with a timeout::

    sage: subprocess.call(["sage", "-t", "--gdb", "1second.rst"], cwd=tests_dir, stdin=open(os.devnull))  # long time, optional: gdb
    exec gdb ...
    Running doctests...
    Doctesting 1 file.
    sage -t 1second.rst
        [2 tests, 1.0 s]
    ----------------------------------------------------------------------
    All tests passed!
    ----------------------------------------------------------------------
    ...
    0
    sage: subprocess.call(["sage", "-t", "--gdb", "-T" "5", "99seconds.rst"], cwd=tests_dir, stdin=open(os.devnull))  # long time, optional: gdb
    exec gdb ...
    Running doctests...
    Doctesting 1 file.
        Time out
    4
"""
