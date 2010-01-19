.. _chapter-doctesting:

=================================
Parallel Testing the Sage Library
=================================

Doctesting a function ensures that the function performs as claimed by
its documentation. Testing can be performed using one thread or
multiple threads. After compiling a source version of Sage, doctesting
can be run on the whole Sage library, on all modules under a given
directory, or on a specified module only. For the purposes of this
chapter, suppose we have compiled Sage 4.1.1 from source and the top
level Sage directory is

::

    [mvngu@sage sage-4.1.1]$ pwd
    /scratch/mvngu/build/sage-4.1.1

See the section :ref:`chapter-testing` for information on Sage's
automated testing process. The general syntax for doctesting is as
follows. To doctest a module in the library of a version of Sage, use
this syntax::

    /path/to/sage-x.y.z/sage -t [-long] /path/to/sage-x.y.z/path/to/module.py[x]

where ``-long`` is an optional argument. The version of ``sage`` used must
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

    [mvngu@sage sage-4.1.1]$ ./sage -t devel/sage-main/sage/games/sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [6.0 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 6.0 seconds

The numbers output by the test show that testing the sudoku module
takes about six seconds, while testing all specified modules took the
same amount of time. In this case, we only tested one module so it is
not surprising that the total testing time is approximately the same
as the time required to test only that one module. Notice that the
syntax is

::

    [mvngu@sage sage-4.1.1]$ ./sage -t devel/sage-main/sage/games/sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [5.7 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 5.7 seconds
    [mvngu@sage sage-4.1.1]$ ./sage -t "devel/sage-main/sage/games/sudoku.py"
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [5.4 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 5.4 seconds

but not

::

    [mvngu@sage sage-4.1.1]$ ./sage -t sage/games/sudoku.py
    ERROR: File ./sage/games/sudoku.py is missing
    exit code: 1

    ----------------------------------------------------------------------
    The following tests failed:

    ./sage/games/sudoku.py
    Total time for all tests: 0.0 seconds
    [mvngu@sage sage-4.1.1]$ ./sage -t "sage/games/sudoku.py"
    ERROR: File ./sage/games/sudoku.py is missing
    exit code: 1

    ----------------------------------------------------------------------
    The following tests failed:

    ./sage/games/sudoku.py
    Total time for all tests: 0.0 seconds

We can also first ``cd`` to the directory containing the module
``sudoku.py`` and doctest that module as follows::

    [mvngu@sage sage-4.1.1]$ cd devel/sage-main/sage/games/
    [mvngu@sage games]$ ls
    all.py    __init__.py         sudoku_backtrack.pyx
    hexad.py  sudoku_backtrack.c  sudoku.py
    [mvngu@sage games]$ ../../../../sage -t sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [5.1 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 5.1 seconds

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

    [mvngu@sage sage-4.1.1]$ ./sage -t devel/sage-main/sage/games/sudoku.py
    sage -t  "devel/sage-main/sage/games/sudoku.py"
             [5.1 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 5.1 seconds
    [mvngu@sage sage-4.1.1]$ ./sage -t "devel/sage-main/sage/games/sudoku.py"
    sage -t  "devel/sage-main/sage/games/sudoku.py"
            [5.0 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 5.0 seconds

With a regular user account, the following syntax is not recommended
as we are using a system-wide Sage installation (if it exists)::

    [mvngu@sage sage-4.1.1]$ sage -t devel/sage-main/sage/games/sudoku.py
    Traceback (most recent call last):
      File "/usr/local/sage/local/bin/sage-test", line 49, in
        os.makedirs(TMP)
      File "/usr/local/sage/local/lib/python/os.py", line 157, in makedirs
        mkdir(name, mode)
    OSError: [Errno 13] Permission denied: '/usr/local/sage/tmp/tmp'
    [mvngu@sage sage-4.1.1]$ sage -t "devel/sage-main/sage/games/sudoku.py"
    Traceback (most recent call last):
      File "/usr/local/sage/local/bin/sage-test", line 49, in
        os.makedirs(TMP)
      File "/usr/local/sage/local/lib/python/os.py", line 157, in makedirs
        mkdir(name, mode)
    OSError: [Errno 13] Permission denied: '/usr/local/sage/tmp/tmp'

In this case, we received a permission error because the system-wide
Sage installation attempts to write some data to a system-wide
directory using our login privileges. The system-wide directory is a
temporary directory under the system-wide ``SAGE_ROOT``. Most likely a
system-wide Sage installation was performed by a system administrator
using an account with more privileges than a regular user. As a
regular user, we cannot write to directories where we do not have
write permission. The following syntax is also discouraged when we
login as a regular user::

    [mvngu@sage sage-4.1.1]$ cd
    [mvngu@sage ~]$ sage -t devel/sage-main/sage/games/sudoku.py
    Traceback (most recent call last):
      File "/usr/local/sage/local/bin/sage-test", line 49, in
        os.makedirs(TMP)
      File "/usr/local/sage/local/lib/python/os.py", line 157, in makedirs
        mkdir(name, mode)
    OSError: [Errno 13] Permission denied: '/usr/local/sage/tmp/tmp'
    [mvngu@sage ~]$ sage -t "devel/sage-main/sage/games/sudoku.py"
    Traceback (most recent call last):
      File "/usr/local/sage/local/bin/sage-test", line 49, in
        os.makedirs(TMP)
      File "/usr/local/sage/local/lib/python/os.py", line 157, in makedirs
        mkdir(name, mode)
    OSError: [Errno 13] Permission denied: '/usr/local/sage/tmp/tmp'

This is exactly the same as the previous syntax because in both cases
we attempted to doctest modules in a system-wide Sage installation
using privileges of a regular user. If we do not have permission to
read or write somewhere on a system, we cannot read or write
there. As a regular user, we do not usually have privileges to read
the directory ``/root`` nor do we have privileges to write to the root
directory::

    [mvngu@sage sage-4.1.1]$ ls /root/
    ls: cannot open directory /root/: Permission denied
    [mvngu@sage sage-4.1.1]$ cd /
    [mvngu@sage /]$ touch demo.txt
    touch: cannot touch `demo.txt': Permission denied


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
and then using two threads. For this example, suppose we want to test
all the modules under ``sage/crypto/``. We can use a syntax similar to
that shown above to achieve this::

    [mvngu@sage sage-4.1.1]$ ./sage -t devel/sage-main/sage/crypto/
    sage -t  "devel/sage-main/sage/crypto/lfsr.py"
             [2.5 s]
    sage -t  "devel/sage-main/sage/crypto/cryptosystem.py"
             [1.9 s]
    sage -t  "devel/sage-main/sage/crypto/block_cipher/miniaes.py"
             [2.5 s]
    sage -t  "devel/sage-main/sage/crypto/block_cipher/all.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/block_cipher/__init__.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/classical.py"
             [2.7 s]
    sage -t  "devel/sage-main/sage/crypto/mq/mpolynomialsystem.py"
             [8.7 s]
    sage -t "devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py"
             [1.9 s]
    sage -t  "devel/sage-main/sage/crypto/mq/__init__.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/mq/sbox.py"
             [2.8 s]
    sage -t  "devel/sage-main/sage/crypto/mq/sr.py"
             [4.9 s]
    sage -t  "devel/sage-main/sage/crypto/stream_cipher.py"
             [1.9 s]
    sage -t  "devel/sage-main/sage/crypto/all.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/stream.py"
             [1.9 s]
    sage -t  "devel/sage-main/sage/crypto/__init__.py"
             [0.1 s]
    sage -t  "devel/sage-main/sage/crypto/classical_cipher.py"
             [1.9 s]
    sage -t  "devel/sage-main/sage/crypto/cipher.py"
             [1.9 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 35.7 seconds

Now we do the same thing, but this time we also use the optional
argument ``-long``::

    [mvngu@sage sage-4.1.1]$ ./sage -t -long devel/sage-main/sage/crypto/
    sage -t -long "devel/sage-main/sage/crypto/lfsr.py"
                  [1.9 s]
    sage -t -long "devel/sage-main/sage/crypto/cryptosystem.py"
                  [2.0 s]
    sage -t -long "devel/sage-main/sage/crypto/block_cipher/miniaes.py"
                  [2.6 s]
    sage -t -long "devel/sage-main/sage/crypto/block_cipher/all.py"
                  [0.1 s]
    sage -t -long "devel/sage-main/sage/crypto/block_cipher/__init__.py"
                  [0.1 s]
    sage -t -long "devel/sage-main/sage/crypto/classical.py"
                  [2.7 s]
    sage -t -long "devel/sage-main/sage/crypto/mq/mpolynomialsystem.py"
                  [8.7 s]
    sage -t -long "devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py"
                  [2.2 s]
    sage -t -long "devel/sage-main/sage/crypto/mq/__init__.py"
                  [0.1 s]
    sage -t -long "devel/sage-main/sage/crypto/mq/sbox.py"
                  [2.9 s]
    sage -t -long "devel/sage-main/sage/crypto/mq/sr.py"
                  [56.6 s]
    sage -t -long "devel/sage-main/sage/crypto/stream_cipher.py"
                  [2.5 s]
    sage -t -long "devel/sage-main/sage/crypto/all.py"
                  [0.1 s]
    sage -t -long "devel/sage-main/sage/crypto/stream.py"
                  [1.9 s]
    sage -t -long "devel/sage-main/sage/crypto/__init__.py"
                  [0.1 s]
    sage -t -long "devel/sage-main/sage/crypto/classical_cipher.py"
                  [1.9 s]
    sage -t -long "devel/sage-main/sage/crypto/cipher.py"
                  [1.9 s]

    ----------------------------------------------------------------------
    All tests passed!
    Total time for all tests: 88.0 seconds

Notice the time difference between the first set of tests and the
second set, which uses the optional argument ``-long``. Many tests in the
Sage library are flagged with ``# long time`` because these are known to
take a long time to run through. Without using the optional ``-long``
argument, the module ``sage/crypto/mq/sr.py`` took about five
seconds. With this optional argument, it required 57 seconds to run
through all tests in that module. Here is a snippet of a function in
the module ``sage/crypto/mq/sr.py`` with a doctest that has been flagged
as taking a long time::

    def test_consistency(max_n=2, **kwargs):
        r"""
        Test all combinations of ``r``, ``c``, ``e`` and ``n`` in ``(1,
	2)`` for consistency of random encryptions and their polynomial
        systems. `\GF{2}` and `\GF{2^e}` systems are tested. This test
        takes
        a while.

        INPUT:

        - ``max_n`` - maximal number of rounds to consider (default: 2)
        - ``kwargs`` - are passed to the SR constructor

        TESTS::

            sage: from sage.crypto.mq.sr import test_consistency
            sage: test_consistency(1) # long time -- calling w/ max_n = 2 requires a LOT of RAM (>> 2GB, evidently).  Calling w/ max_n = 1 is far more manageable.
            True

        The above doctest used to fail on a machine with "only" 2GB RAM.
        Using ``max_n = 1`` appears to be a more reasonable memory usage.
        """

Now we doctest the same directory in parallel using two threads::

    [mvngu@sage sage-4.1.1]$ ./sage -tp 2 devel/sage-main/sage/crypto/
    Global iterations: 1
    File iterations: 1
    Using cached timings to run longest doctests first.
    Doctesting 17 files doing 2 jobs in parallel
    sage -t  devel/sage-main/sage/crypto/lfsr.py
             [2.7 s]
    sage -t  devel/sage-main/sage/crypto/cryptosystem.py
             [2.0 s]
    sage -t  devel/sage-main/sage/crypto/mq/mpolynomialsystem.py
             [9.4 s]
    sage -t  devel/sage-main/sage/crypto/mq/sr.py
             [5.2 s]
    sage -t  devel/sage-main/sage/crypto/classical.py
             [2.8 s]
    sage -t  devel/sage-main/sage/crypto/mq/sbox.py
             [3.2 s]
    sage -t  devel/sage-main/sage/crypto/block_cipher/miniaes.py
             [2.6 s]
    sage -t  devel/sage-main/sage/crypto/stream_cipher.py
             [2.0 s]
    sage -t  devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py
             [2.0 s]
    sage -t  devel/sage-main/sage/crypto/classical_cipher.py
             [2.1 s]
    sage -t  devel/sage-main/sage/crypto/cipher.py
             [2.1 s]
    sage -t  devel/sage-main/sage/crypto/__init__.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/block_cipher/__init__.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/mq/__init__.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/block_cipher/all.py
             [0.1 s]
    sage -t  devel/sage-main/sage/crypto/stream.py
             [2.0 s]
    sage -t  devel/sage-main/sage/crypto/all.py
             [0.1 s]

    ----------------------------------------------------------------------
    All tests passed!
    Timings have been updated.
    Total time for all tests: 19.3 seconds

    [mvngu@sage sage-4.1.1]$ ./sage -tp 2 -long devel/sage-main/sage/crypto/
    Global iterations: 1
    File iterations: 1
    No long cached timings exist; will create upon successful finish.
    Doctesting 17 files doing 2 jobs in parallel
    sage -t -long devel/sage-main/sage/crypto/cryptosystem.py
             [2.7 s]
    sage -t -long devel/sage-main/sage/crypto/lfsr.py
             [2.7 s]
    sage -t -long devel/sage-main/sage/crypto/stream_cipher.py
             [2.2 s]
    sage -t -long devel/sage-main/sage/crypto/all.py
             [0.1 s]
    sage -t -long devel/sage-main/sage/crypto/classical.py
             [3.0 s]
    sage -t -long devel/sage-main/sage/crypto/__init__.py
             [0.1 s]
    sage -t -long devel/sage-main/sage/crypto/stream.py
             [2.1 s]
    sage -t -long devel/sage-main/sage/crypto/classical_cipher.py
             [2.1 s]
    sage -t -long devel/sage-main/sage/crypto/cipher.py
             [2.1 s]
    sage -t -long devel/sage-main/sage/crypto/block_cipher/all.py
             [0.1 s]
    sage -t -long devel/sage-main/sage/crypto/block_cipher/__init__.py
             [0.1 s]
    sage -t -long devel/sage-main/sage/crypto/block_cipher/miniaes.py
             [2.8 s]
    sage -t -long devel/sage-main/sage/crypto/mq/mpolynomialsystemgenerator.py
             [2.0 s]
    sage -t -long devel/sage-main/sage/crypto/mq/__init__.py
             [0.1 s]
    sage -t -long devel/sage-main/sage/crypto/mq/sbox.py
             [3.1 s]
    sage -t -long devel/sage-main/sage/crypto/mq/mpolynomialsystem.py
             [9.1 s]
    sage -t -long devel/sage-main/sage/crypto/mq/sr.py
             [56.0 s]

    ----------------------------------------------------------------------
    All tests passed!
    Timings have been updated.
    Total time for all tests: 71.8 seconds

As the number of threads increases, the total testing time
decreases. To minimize confusion, it is also a good idea to explicitly
specify the path name of the directory we want to doctest and not a
symbolic link to that directory. In the above examples, the symbolic
link ``devel/sage`` points to the directory ``devel/sage-main``, but the
actual path to the directory has been specified instead of its
symbolic link.


Parallel testing the whole Sage library
=======================================

The main Sage library resides in the directory
``SAGE_ROOT/devel/sage-main/``. We can use the syntax described above
to doctest the main library using multiple threads. When doing release
management or patching the main Sage library, a release manager would
parallel test the library using ten or more threads::

    [mvngu@sage sage-4.1.1]$ ./sage -tp 10 -long devel/sage-main/

Another way is to edit the file ``makefile`` in the top level Sage
directory so that the variable ``NUM_THREADS`` is set to ``10``::

    # How many threads should be used when doing parallel testing (and
    # sometime in the future, parallel building)?
    NUM_THREADS=10

After saving all changes to ``makefile``, we can parallel test with the
``-long`` option using ten threads::

    [mvngu@sage sage-4.1.1]$ make ptestlong

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
  these command, see the files ``SAGE_ROOT/makefile`` and
  ``SAGE_ROOT/local/bin/sage-maketest``.

* ``make testlong`` --- This command doctests the standard
  documentation:

  * ``SAGE_ROOT/devel/sage/doc/common``
  * ``SAGE_ROOT/devel/sage/doc/en``
  * ``SAGE_ROOT/devel/sage/doc/fr``

  and then the Sage library. Doctesting is run with the optional
  argument ``-long``. See the file ``SAGE_ROOT/makefile`` for further
  details.

* ``make ptest`` --- Similar to the commands ``make test`` and ``make
  check``. However, doctesting is run with the number of threads as
  specified by the variable ``NUM_THREADS``. See the file
  ``SAGE_ROOT/makefile`` for further details.

* ``make ptestlong`` --- Similar to the command ``make ptest``, but
  using the optional argument ``-long`` for doctesting.


Beyond the Sage library
=======================

The doctesting scripts of a Sage installation currently have limited
support for doctesting of modules outside of that Sage library. We
cannot use the doctesting scripts of Sage 4.1.1 to doctest modules in,
say, Sage 4.1. Doing so would result in errors::

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
