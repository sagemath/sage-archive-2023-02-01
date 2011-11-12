.. _chapter-doctesting:

===========================
Doctesting the Sage Library
===========================

Doctesting a function ensures that the function performs as claimed by
its documentation. Testing can be performed using one thread or
multiple threads. After compiling a source version of Sage, doctesting
can be run on the whole Sage library, on all modules under a given
directory, or on a specified module only. For the purposes of this
chapter, suppose we have compiled Sage 4.8 from source and the top
level Sage directory is

::

    [jdemeyer@sage sage-4.8]$ pwd
    /scratch/jdemeyer/build/sage-4.8

See the section :ref:`chapter-testing` for information on Sage's
automated testing process. The general syntax for doctesting is as
follows. To doctest a module in the library of a version of Sage, use
this syntax::

    /path/to/sage-x.y.z/sage -t [--long] /path/to/sage-x.y.z/path/to/module.py[x]

where ``--long`` is an optional argument. The version of ``sage`` used must
match the version of Sage containing the module we want to doctest. A
Sage module can be either a Python script (with the file extension
".py") or it can be a Cython script, in which case it has the file
extension ".pyx".


Testing a module
================

Say we want to run all tests in the sudoku module
``sage/games/sudoku.py``. In a terminal window, first we ``cd`` to the
top level Sage directory of our local Sage installation. Now  we can
start doctesting as demonstrated in the following terminal session::

    [jdemeyer@sage sage-4.8]$ ./sage -t devel/sage-main/sage/games/sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [7.3 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 7.3 seconds

The numbers output by the test show that testing the sudoku module
takes about six seconds, while testing all specified modules took the
same amount of time. In this case, we only tested one module so it is
not surprising that the total testing time is approximately the same
as the time required to test only that one module. Notice that the
syntax is ::

    [jdemeyer@sage sage-4.8]$ ./sage -t devel/sage-main/sage/games/sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [7.3 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 7.3 seconds
    [jdemeyer@sage sage-4.8]$ ./sage -t "devel/sage-main/sage/games/sudoku.py"
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [7.5 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 7.6 seconds

but not

::

    [jdemeyer@sage sage-4.8]$ ./sage -t sage/games/sudoku.py
    ERROR: File ./sage/games/sudoku.py is missing
    exit code: 1

    ----------------------------------------------------------------------
    The following tests failed:

    ./sage/games/sudoku.py
    Total time for all tests: 0.0 seconds
    [jdemeyer@sage sage-4.8]$ ./sage -t "sage/games/sudoku.py"
    ERROR: File ./sage/games/sudoku.py is missing

    ----------------------------------------------------------------------
    The following tests failed:


            ./sage/games/sudoku.py # File not found
    Total time for all tests: 0.0 seconds

We can also first ``cd`` to the directory containing the module
``sudoku.py`` and doctest that module as follows::

    [jdemeyer@sage sage-4.8]$ cd devel/sage-main/sage/games/
    [jdemeyer@sage games]$ ls
    __init__.py  hexad.py       sudoku.py           sudoku_backtrack.pyx
    all.py       quantumino.py  sudoku_backtrack.c
    [jdemeyer@sage games]$ ../../../../sage -t sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [7.1 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 7.1 seconds

In all of the above terminal sessions, we used a local installation of
Sage to test its own modules. Even if we have a system-wide Sage
installation, using that version to doctest the modules of a local
installation is a recipe for confusion.


Troubleshooting
===============

To doctest modules of a Sage installation, from a terminal window we
first ``cd`` to the top level directory of that Sage installation,
otherwise known as the ``SAGE_ROOT`` of that installation. When we
run tests, we use that particular Sage installation via the syntax
``./sage``; notice the "dot-forward-slash" at the front of
``sage``. This is a precaution against confusion that can arise when
our system has multiple Sage installations. For example, the following
syntax is acceptable because we explicitly specify the Sage
installation in the current ``SAGE_ROOT``::

    [jdemeyer@sage sage-4.8]$ ./sage -t devel/sage-main/sage/games/sudoku.py
    ./sage -t "devel/sage-main/sage/games/sudoku.py"
             [6.9 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 6.9 seconds
    [jdemeyer@sage sage-4.8]$ ./sage -t "devel/sage-main/sage/games/sudoku.py"
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [7.7 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 7.7 seconds

The following syntax is not recommended as we are using a system-wide
Sage installation (if it exists):

.. skip

::

    [jdemeyer@sage sage-4.8]$ sage -t devel/sage-main/sage/games/sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
    **********************************************************************
    File "/home/jdemeyer/sage/sage-4.8/devel/sage-main/sage/games/sudoku.py", line 515:
        sage: h.solve(algorithm='backtrack').next()
    Exception raised:
        Traceback (most recent call last):
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1231, in run_one_test
            self.run_one_example(test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/sagedoctest.py", line 38, in run_one_example
            OrigDocTestRunner.run_one_example(self, test, example, filename, compileflags)
          File "/usr/local/sage/local/bin/ncadoctest.py", line 1172, in run_one_example
            compileflags, 1) in test.globs
          File "<doctest __main__.example_13[4]>", line 1, in <module>
            h.solve(algorithm='backtrack').next()###line 515:
        sage: h.solve(algorithm='backtrack').next()
          File "/home/jdemeyer/.sage/tmp/sudoku.py", line 607, in solve
            for soln in gen:
          File "/home/jdemeyer/.sage/tmp/sudoku.py", line 719, in backtrack
            from sudoku_backtrack import backtrack_all
        ImportError: No module named sudoku_backtrack
    **********************************************************************
    [...more errors...]
    2 items had failures:
       4 of  15 in __main__.example_13
       2 of   8 in __main__.example_14
    ***Test Failed*** 6 failures.
    For whitespace errors, see the file /home/jdemeyer/.sage//tmp/.doctest_sudoku.py
             [21.1 s]

    ----------------------------------------------------------------------
    The following tests failed:


            sage -t  "devel/sage-main/sage/games/sudoku.py"
    Total time for all tests: 21.3 seconds

In this case, we received an error because the system-wide Sage
installation is a different (older) version than the one we are
using for Sage development.  Make sure you always test the files
with the correct version of Sage.

Parallel testing many modules
=============================

So far we have used a single thread to doctest a module in the Sage
library. There are hundreds, even thousands of modules in the Sage
library. Testing them all using one thread would take a few
hours. Depending on our hardware, this could take up to six hours or
more. On a multi-core system, parallel doctesting can significantly
reduce the testing time. Unless we also want to use our computer
while doctesting in parallel, we can choose to devote all the cores
of our system for parallel testing.

Let us doctest all modules in a directory, first using a single thread
and then using four threads. For this example, suppose we want to test
all the modules under ``sage/crypto/``. We can use a syntax similar to
that shown above to achieve this::

    [jdemeyer@sage sage-4.8]$ ./sage -t devel/sage-main/sage/crypto/
    sage -t  "devel/sage-main/sage/crypto/block_cipher/__init__.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/block_cipher/miniaes.py"
             [5.5 s]
    sage -t  "devel/sage-main/sage/crypto/block_cipher/all.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/block_cipher/sdes.py"
             [4.2 s]
    sage -t  "devel/sage-main/sage/crypto/__init__.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/stream.py"
             [3.7 s]
    sage -t  "devel/sage-main/sage/crypto/classical_cipher.py"
             [5.1 s]
    sage -t  "devel/sage-main/sage/crypto/boolean_function.pyx"
             [7.3 s]
    sage -t  "devel/sage-main/sage/crypto/lattice.py"
             [3.7 s]
    sage -t  "devel/sage-main/sage/crypto/util.py"
             [3.4 s]
    sage -t  "devel/sage-main/sage/crypto/cryptosystem.py"
             [3.6 s]
    sage -t  "devel/sage-main/sage/crypto/all.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/mq/__init__.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/mq/sbox.py"
             [4.4 s]
    sage -t  "devel/sage-main/sage/crypto/mq/mpolynomialsystem.py"
             [12.8 s]
    sage -t  "devel/sage-main/sage/crypto/mq/sr.py"
             [10.6 s]
    sage -t  "devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py"
             [3.4 s]
    sage -t  "devel/sage-main/sage/crypto/cipher.py"
             [3.4 s]
    sage -t  "devel/sage-main/sage/crypto/classical.py"
             [13.8 s]
    sage -t  "devel/sage-main/sage/crypto/public_key/blum_goldwasser.py"
             [3.5 s]
    sage -t  "devel/sage-main/sage/crypto/public_key/__init__.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/public_key/all.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/stream_cipher.py"
             [3.4 s]
    sage -t  "devel/sage-main/sage/crypto/lfsr.py"
             [3.5 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 96.1 seconds

Now we do the same thing, but this time we also use the optional
argument ``--long``::

    [jdemeyer@sage sage-4.8]$ ./sage -t --long devel/sage-main/sage/crypto/
    sage -t --long "devel/sage-main/sage/crypto/block_cipher/__init__.py"
             [0.1 s]
    sage -t --long "devel/sage-main/sage/crypto/block_cipher/miniaes.py"
             [4.1 s]
    sage -t --long "devel/sage-main/sage/crypto/block_cipher/all.py"
             [0.1 s]
    sage -t --long "devel/sage-main/sage/crypto/block_cipher/sdes.py"
             [3.9 s]
    sage -t --long "devel/sage-main/sage/crypto/__init__.py"
             [0.0 s]
    sage -t --long "devel/sage-main/sage/crypto/stream.py"
             [3.3 s]
    sage -t --long "devel/sage-main/sage/crypto/classical_cipher.py"
             [3.9 s]
    sage -t --long "devel/sage-main/sage/crypto/boolean_function.pyx"
             [7.2 s]
    sage -t --long "devel/sage-main/sage/crypto/lattice.py"
             [3.4 s]
    sage -t --long "devel/sage-main/sage/crypto/util.py"
             [3.3 s]
    sage -t --long "devel/sage-main/sage/crypto/cryptosystem.py"
             [3.4 s]
    sage -t --long "devel/sage-main/sage/crypto/all.py"
             [0.1 s]
    sage -t --long "devel/sage-main/sage/crypto/mq/__init__.py"
             [0.1 s]
    sage -t --long "devel/sage-main/sage/crypto/mq/sbox.py"
             [3.5 s]
    sage -t --long "devel/sage-main/sage/crypto/mq/mpolynomialsystem.py"
             [11.8 s]
    sage -t --long "devel/sage-main/sage/crypto/mq/sr.py"
             [96.8 s]
    sage -t --long "devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py"
             [2.9 s]
    sage -t --long "devel/sage-main/sage/crypto/cipher.py"
             [3.2 s]
    sage -t --long "devel/sage-main/sage/crypto/classical.py"
             [13.6 s]
    sage -t --long "devel/sage-main/sage/crypto/public_key/blum_goldwasser.py"
             [3.2 s]
    sage -t --long "devel/sage-main/sage/crypto/public_key/__init__.py"
             [0.1 s]
    sage -t --long "devel/sage-main/sage/crypto/public_key/all.py"
             [0.1 s]
    sage -t --long "devel/sage-main/sage/crypto/stream_cipher.py"
             [3.4 s]
    sage -t --long "devel/sage-main/sage/crypto/lfsr.py"
             [3.0 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 174.3 seconds

Notice the time difference between the first set of tests and the
second set, which uses the optional argument ``--long``. Many tests in the
Sage library are flagged with ``# long time`` because these are known to
take a long time to run through. Without using the optional ``--long``
argument, the module ``sage/crypto/mq/sr.py`` took about ten
seconds. With this optional argument, it required 97 seconds to run
through all tests in that module. Here is a snippet of a function in
the module ``sage/crypto/mq/sr.py`` with a doctest that has been flagged
as taking a long time::

    def test_consistency(max_n=2, **kwargs):
        r"""
        Test all combinations of ``r``, ``c``, ``e`` and ``n`` in ``(1,
        2)`` for consistency of random encryptions and their polynomial
        systems. `\GF{2}` and `\GF{2^e}` systems are tested. This test takes
        a while.

        INPUT:

        - ``max_n`` -- maximal number of rounds to consider (default: 2)
        - ``kwargs`` -- are passed to the SR constructor

        TESTS:

        The following test called with ``max_n`` = 2 requires a LOT of RAM
        (much more than 2GB).  Since this might cause the doctest to fail
        on machines with "only" 2GB of RAM, we test ``max_n`` = 1, which
        has a more reasonable memory usage. ::

            sage: from sage.crypto.mq.sr import test_consistency
            sage: test_consistency(1)  # long time (80s on sage.math, 2011)
            True
        """

Now we doctest the same directory in parallel using 4 threads::

    [jdemeyer@sage sage-4.8]$ ./sage -tp 4 devel/sage-main/sage/crypto/
    Global iterations: 1
    File iterations: 1
    Using cached timings to run longest doctests first.
    Doctesting 24 files doing 4 jobs in parallel
    sage -t  devel/sage-main/sage/crypto/__init__.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/lattice.py
             [3.3 s]
    sage -t  devel/sage-main/sage/crypto/stream.py
             [3.5 s]
    sage -t  devel/sage-main/sage/crypto/classical_cipher.py
             [4.0 s]
    sage -t  devel/sage-main/sage/crypto/all.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/util.py
             [3.4 s]
    sage -t  devel/sage-main/sage/crypto/cryptosystem.py
             [3.4 s]
    sage -t  devel/sage-main/sage/crypto/boolean_function.pyx
             [6.9 s]
    sage -t  devel/sage-main/sage/crypto/cipher.py
             [3.3 s]
    sage -t  devel/sage-main/sage/crypto/block_cipher/__init__.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/lfsr.py
             [3.3 s]
    sage -t  devel/sage-main/sage/crypto/stream_cipher.py
             [3.4 s]
    sage -t  devel/sage-main/sage/crypto/block_cipher/all.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/mq/__init__.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/block_cipher/miniaes.py
             [4.0 s]
    sage -t  devel/sage-main/sage/crypto/block_cipher/sdes.py
             [3.6 s]
    sage -t  devel/sage-main/sage/crypto/mq/sbox.py
             [4.0 s]
    sage -t  devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py
             [3.2 s]
    sage -t  devel/sage-main/sage/crypto/public_key/blum_goldwasser.py
             [3.4 s]
    sage -t  devel/sage-main/sage/crypto/public_key/__init__.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/classical.py
             [14.3 s]
    sage -t  devel/sage-main/sage/crypto/public_key/all.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/mq/sr.py
             [9.3 s]
    sage -t  devel/sage-main/sage/crypto/mq/mpolynomialsystem.py
             [12.0 s]

    ----------------------------------------------------------------------
    All tests passed!
    Timings have been updated.
    Total time for all tests: 23.7 seconds
    [jdemeyer@sage sage-4.8]$ ./sage -tp 4 --long devel/sage-main/sage/crypto/
    Global iterations: 1
    File iterations: 1
    Using long cached timings to run longest doctests first.
    Doctesting 24 files doing 4 jobs in parallel
    sage -t --long devel/sage-main/sage/crypto/__init__.py
             [0.1 s]
    sage -t --long devel/sage-main/sage/crypto/stream.py
             [3.2 s]
    sage -t --long devel/sage-main/sage/crypto/lattice.py
             [3.3 s]
    sage -t --long devel/sage-main/sage/crypto/classical_cipher.py
             [4.1 s]
    sage -t --long devel/sage-main/sage/crypto/all.py
             [0.1 s]
    sage -t --long devel/sage-main/sage/crypto/util.py
             [3.1 s]
    sage -t --long devel/sage-main/sage/crypto/cryptosystem.py
             [3.3 s]
    sage -t --long devel/sage-main/sage/crypto/boolean_function.pyx
             [7.0 s]
    sage -t --long devel/sage-main/sage/crypto/cipher.py
             [3.2 s]
    sage -t --long devel/sage-main/sage/crypto/block_cipher/__init__.py
             [0.1 s]
    sage -t --long devel/sage-main/sage/crypto/stream_cipher.py
             [3.2 s]
    sage -t --long devel/sage-main/sage/crypto/block_cipher/all.py
             [0.1 s]
    sage -t --long devel/sage-main/sage/crypto/lfsr.py
             [3.4 s]
    sage -t --long devel/sage-main/sage/crypto/mq/__init__.py
             [0.1 s]
    sage -t --long devel/sage-main/sage/crypto/block_cipher/miniaes.py
             [4.2 s]
    sage -t --long devel/sage-main/sage/crypto/block_cipher/sdes.py
             [4.0 s]
    sage -t --long devel/sage-main/sage/crypto/mq/sbox.py
             [3.8 s]
    sage -t --long devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py
             [3.1 s]
    sage -t --long devel/sage-main/sage/crypto/classical.py
             [13.8 s]
    sage -t --long devel/sage-main/sage/crypto/public_key/__init__.py
             [0.0 s]
    sage -t --long devel/sage-main/sage/crypto/public_key/all.py
             [0.0 s]
    sage -t --long devel/sage-main/sage/crypto/public_key/blum_goldwasser.py
             [3.1 s]
    sage -t --long devel/sage-main/sage/crypto/mq/mpolynomialsystem.py
             [11.3 s]
    sage -t --long devel/sage-main/sage/crypto/mq/sr.py
             [95.4 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 109.4 seconds

As the number of threads increases, the total testing time
decreases. To minimize confusion, it is also a good idea to explicitly
specify the path name of the directory we want to doctest and not a
symbolic link to that directory. In the above examples, the symbolic
link ``devel/sage`` points to the directory ``devel/sage-main``, but the
actual path to the directory has been specified instead of its
symbolic link.

.. _section-parallel-test-whole-library:

Parallel testing the whole Sage library
=======================================

The main Sage library resides in the directory
``SAGE_ROOT/devel/sage-main/``. We can use the syntax described above
to doctest the main library using multiple threads. When doing release
management or patching the main Sage library, a release manager would
parallel test the library using 10 threads with the following command::

    [jdemeyer@sage sage-4.8]$ ./sage -tp 10 -long devel/sage-main/

Another way is run ``make ptestlong``, which builds Sage (if necessary),
builds the Sage documentation (if necessary), and then runs parallel
doctests.  This determines the number of threads by reading the
environment variable :envvar:`MAKE`: if it is set to ``make -j12``, then
use 12 threads.  If :envvar:`MAKE` is not set, then by default it uses
the number of CPU cores (as determined by the Python function
``multiprocessing.cpu_count()``) with a minimum of 2 and a maximum of 8.

In any case, this will test the Sage library with multiple threads::

    [jdemeyer@sage sage-4.8]$ make ptestlong

Any of the following commands would also doctest the Sage library or
one of its clones::

    make test
    make check
    make testlong
    make ptest
    make ptestlong

In each case, testing is performed on the directory that is pointed to
by the symbolic link ``devel/sage``.

* ``make test`` and ``make check`` --- These two commands run the same
  set of tests. First the Sage standard documentation is tested,
  i.e. the documentation that resides in

  * ``SAGE_ROOT/devel/sage/doc/common``
  * ``SAGE_ROOT/devel/sage/doc/en``
  * ``SAGE_ROOT/devel/sage/doc/fr``

  Finally, the commands doctest the Sage library. For more details on
  these command, see the files ``SAGE_ROOT/Makefile`` and
  ``SAGE_ROOT/local/bin/sage-maketest``.

* ``make testlong`` --- This command doctests the standard
  documentation:

  * ``SAGE_ROOT/devel/sage/doc/common``
  * ``SAGE_ROOT/devel/sage/doc/en``
  * ``SAGE_ROOT/devel/sage/doc/fr``

  and then the Sage library. Doctesting is run with the optional
  argument ``-long``. See the file ``SAGE_ROOT/Makefile`` for further
  details.

* ``make ptest`` --- Similar to the commands ``make test`` and ``make
  check``. However, doctesting is run with the number of threads as
  described above for ``make ptestlong``.

* ``make ptestlong`` --- Similar to the command ``make ptest``, but
  using the optional argument ``-long`` for doctesting.


Beyond the Sage library
=======================

The doctesting scripts of a Sage installation currently have limited
support for doctesting of modules outside of the Sage library for
that version of Sage. We cannot use the doctesting scripts of Sage
4.1.1 to doctest modules in, say, Sage 4.1. Doing so would result
in errors::

    [mvngu@sage sage-4.1.1]$ ./sage -t ../sage-4.1/devel/sage-main/sage/games/sudoku.py
    sage -t  "../sage-4.1/devel/sage-main/sage/games/sudoku.py"
      File "./sudoku.py", line 18
        from ../sage-4.1/devel/sage-main/sage/games/sudoku import *
               ^
    SyntaxError: invalid syntax

             [0.2 s]
    exit code: 1024

    ----------------------------------------------------------------------
    The following tests failed:


           sage -t  "../sage-4.1/devel/sage-main/sage/games/sudoku.py"
    Total time for all tests: 0.2 seconds

However, suppose we have a Python script called ``my_python_script.py``
that uses the Sage library. Our Python script has the following
content::

    [mvngu@sage build]$ cat my_python_script.py
    from sage.all_cmdline import *   # import sage library

    def square(n):
        """
    	Return the square of n.

	EXAMPLES::

            sage: square(2)
            4
        """
        return n**2

We can use any version of Sage to doctest our Python script, so long
as that version of Sage has features that are used in our script. For
example, we can use both Sage 4.1.1 and 4.1 to doctest the above
Python script::

    [mvngu@sage build]$ sage-4.1/sage -t my_python_script.py
    sage -t  "my_python_script.py"
             [1.3 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 1.3 seconds
    [mvngu@sage build]$ sage-4.1.1/sage -t my_python_script.py
    sage -t  "my_python_script.py"
             [1.4 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 1.4 seconds

Doctesting can also be performed on Sage scripts. Say we have a Sage
script called ``my_sage_script.sage`` with the following content::

    [mvngu@sage build]$ cat my_sage_script.sage
    def cube(n):
        r"""
        Return the cube of n.

        EXAMPLES::

            sage: cube(2)
            8
        """
        return n**3

This must be converted to an equivalent Python script prior to
doctesting. First, we use Sage to convert ``my_sage_script.sage`` to
an equivalent Python script called ``my_sage_script.py``::

    [mvngu@sage build]$ sage-4.1.1/sage my_sage_script.sage
    [mvngu@sage build]$ cat my_sage_script.py
    # This file was *autogenerated* from the file my_sage_script.sage.
    from sage.all_cmdline import *   # import sage library
    _sage_const_3 = Integer(3)
    def cube(n):
        r"""
        Return the cube of n.

        EXAMPLES::

            sage: cube(2)
            8
        """
        return n**_sage_const_3

Doctesting is then performed on that equivalent Python script::

    [mvngu@sage build]$ sage-4.1.1/sage -t my_sage_script.py
    sage -t  "my_sage_script.py"
             [1.5 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 1.5 seconds
