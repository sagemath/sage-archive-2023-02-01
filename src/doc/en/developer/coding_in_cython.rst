.. _chapter-cython:

================
Coding in Cython
================

This chapter discusses Cython, which is a compiled language based on
Python.  The major advantage it has over Python is that code can be
much faster (sometimes orders of magnitude) and can directly call
C and C++ code.  As Cython is essentially a superset of the Python
language, one often doesn’t make a distinction between Cython and 
Python code in Sage (e.g. one talks of the “Sage Python Library”
and “Python Coding Conventions”).

Python is an interpreted language and has no declared data types for
variables. These features make it easy to write and debug, but Python
code can sometimes be slow. Cython code can look a lot like Python,
but it gets translated into C code (often very efficient C code) and
then compiled. Thus it offers a language which is familiar to Python
developers, but with the potential for much greater speed. Cython also
allows Sage developers to interface with C and C++ much easier than
using the Python C API directly.

Cython is a compiled version of Python. It was originally based on
Pyrex but has changed based on what Sage's developers needed; Cython
has been developed in concert with Sage. However, it is an independent
project now, which is used beyond the scope of Sage. As such, it is a
young, but developing language, with young, but developing
documentation. See its web page, http://www.cython.org/, for the most
up-to-date information.



Writing Cython Code in Sage
===========================

There are several ways to create and build Cython code in Sage.

#. In the Sage Notebook, begin any cell with ``%cython``. When you
   evaluate that cell,

   #. It is saved to a file.

   #. Cython is run on it with all the standard Sage libraries
      automatically linked if necessary.

   #. The resulting shared library file (``.so`` / ``.dll`` /
      ``.dylib``) is then loaded into your running instance of Sage.

   #. The functionality defined in that cell is now available for you
      to use in the notebook. Also, the output cell has a link to the C
      program that was compiled to create the ``.so`` file.

   #. A ``cpdef`` or ``def`` function, say ``testfunction``, defined in
      a ``%cython`` cell in a worksheet can be imported and made available
      in a different ``%cython`` cell within the same worksheet by
      importing it as shown below::

          %cython
          from __main__ import testfunction

#. Create an ``.spyx`` file and attach or load it from the command
   line. This is similar to creating a ``%cython`` cell in the
   notebook but works completely from the command line (and not from
   the notebook).

#. Create a ``.pyx`` file and add it to the Sage library.

   #. First, add a listing for the Cython extension to the variable
      ``ext_modules`` in the file
      ``SAGE_ROOT/src/module_list.py``. See the
      ``distutils.extension.Extension`` class for more information on
      creating a new Cython extension.

   #. Run ``sage -b`` to rebuild Sage.

   For example, the file
   ``SAGE_ROOT/src/sage/graphs/chrompoly.pyx`` has the lines::

       Extension('sage.graphs.chrompoly',
                 sources = ['sage/graphs/chrompoly.pyx']),

   in ``module_list.py``. In addition, ``sage.graphs`` is included in
   the ``packages`` list under the Distutils section of ``setup.py``
   since ``chrompoly.pyx`` is contained in the directory
   ``sage/graphs``.


Special Pragmas
===============

If Cython code is either attached or loaded as a ``.spyx`` file or
loaded from the notebook as a ``%cython`` block, the following
pragmas are available:

* clang --- may be either c or c++ indicating whether a C or C++
  compiler should be used.

* clib --- additional libraries to be linked in, the space separated
  list is split and passed to distutils.

* cinclude --- additional directories to search for header files. The
  space separated list is split and passed to distutils.

* cfile -- additional C or C++ files to be compiled

* cargs -- additional parameters passed to the compiler

For example::

    #clang C++
    #clib givaro
    #cinclude /usr/local/include/
    #cargs -ggdb
    #cfile foo.c


Attaching or Loading .spyx Files
================================

The easiest way to try out Cython without having to learn anything
about distutils, etc., is to create a file with the extension
``spyx``, which stands for "Sage Pyrex":

#. Create a file ``power2.spyx``.

#. Put the following in it::

       def is2pow(n):
           while n != 0 and n%2 == 0:
               n = n >> 1
           return n == 1

#. Start the Sage command line interpreter and load the ``spyx`` file
   (this will fail if you do not have a C compiler installed).

   .. skip

   ::

       sage: load "power2.spyx"
       Compiling power2.spyx...
       sage: is2pow(12)
       False

Note that you can change ``power2.spyx``, then load it again and it
will be recompiled on the fly. You can also attach ``power2.spyx`` so
it is reloaded whenever you make changes:

.. skip

::

    sage: attach "power2.spyx"

Cython is used for its speed. Here is a timed test on a 2.6 GHz
Opteron:

.. skip

::

    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 0.60 s, sys: 0.00 s, total: 0.60 s
    Wall time: 0.60 s

Now, the code in the file ``power2.spyx`` is valid Python, and if we
copy this to a file ``powerslow.py`` and load that, we get the
following:

.. skip

::

    sage: load "powerslow.py"
    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 1.01 s, sys: 0.04 s, total: 1.05 s
    Wall time: 1.05 s

By the way, we could gain even a little more speed with the Cython
version with a type declaration, by changing ``def is2pow(n):`` to
``def is2pow(unsigned int n):``.


.. _section_sig_on:

Interrupt and Signal Handling
=============================

When writing Cython code for Sage, special care must be taken to
ensure the code can be interrupted with ``CTRL-C``.  Since Cython
optimizes for speed, Cython normally does not check for interrupts.
For example, code like the following cannot be interrupted:

.. skip

::

    sage: cython('while True: pass')  # DON'T DO THIS

While this is running, pressing ``CTRL-C`` has no effect.  The only
way out is to kill the Sage process.  On certain systems, you can
still quit Sage by typing ``CTRL-\`` (sending a Quit signal) instead
of ``CTRL-C``.

Using ``sig_on()`` and ``sig_off()``
------------------------------------

.. highlight:: cython

To enable interrupt handling, use the ``sig_on()`` and ``sig_off()`` functions.
You should put ``sig_on()`` *before* and ``sig_off()`` *after* any Cython code
which could potentially take a long time.
These two *must always* be called in **pairs**, i.e. every
``sig_on()`` must be matched by a closing ``sig_off()``.

In practice your function will probably look like::

    def sig_example():
        # (some harmless initialization)
        sig_on()
        # (a long computation here, potentially calling a C library)
        sig_off()
        # (some harmless post-processing)
        return something

You can put ``sig_on()`` and ``sig_off()`` in all kinds of Cython
functions: ``def``, ``cdef`` or ``cpdef``.  You cannot put them in
pure Python code (i.e. files with extension ``.py``).

It is possible to put ``sig_on()`` and ``sig_off()`` in different
functions, provided that ``sig_off()`` is called before the function
which calls ``sig_on()`` returns.  The following code is *invalid*::

    # INVALID code because we return from function foo()
    # without calling sig_off() first.
    cdef foo():
        sig_on()

    def f1():
        foo()
        sig_off()

But the following is valid since you cannot call ``foo``
interactively::

    cdef int foo():
        sig_off()
        return 2+2

    def f1():
        sig_on()
        return foo()

For clarity however, it is best to avoid this.  One good example where
the above makes sense is the ``new_gen()`` function in
:ref:`section-pari-library`.

A common mqistake is to put ``sig_off()`` towards the end of a
function (before the ``return``) when the function has multiple
``return`` statements.  So make sure there is a ``sig_off()`` before
*every* ``return`` (and also before every ``raise``).

.. WARNING::

    The code inside ``sig_on()`` should be pure C or Cython code. If
    you call Python code, an interrupt is likely to mess up Python.

    Also, when an interrupt occurs inside ``sig_on()``, code execution
    immediately stops without cleaning up.  For example, any memory
    allocated inside ``sig_on()`` is lost.  See :ref:`advanced-sig`
    for ways to deal with this.

When the user presses ``CTRL-C`` inside ``sig_on()``, execution will
jump back to ``sig_on()`` (the first one if there is a stack) and
``sig_on()`` will raise ``KeyboardInterrupt``.  These can be handled
just like other Python exceptions::

    def catch_interrupts():
        try:
            sig_on()  # This MUST be inside the try
            # (some long computation)
            sig_off()
        except KeyboardInterrupt:
            # (handle interrupt)

Certain C libraries in Sage are written in a way that they will raise
Python exceptions: NTL and PARI are examples of this.  NTL can raise
``RuntimeError`` and PARI can raise ``PariError``.  Since these use
the ``sig_on()`` mechanism, these exceptions can be caught just like
the ``KeyboardInterrupt`` in the example above.

It is possible to stack ``sig_on()`` and ``sig_off()``.  If you do
this, the effect is exactly the same as if only the outer
``sig_on()``/``sig_off()`` was there.  The inner ones will just change
a reference counter and otherwise do nothing.  Make sure that the
number of ``sig_on()`` calls equal the number of ``sig_off()`` calls::

    def stack_sig_on():
        sig_on()
        sig_on()
        sig_on()
        # (some code)
        sig_off()
        sig_off()
        sig_off()

Extra care must be taken with exceptions raised inside ``sig_on()``.
The problem is that, if you do not do anything special, the
``sig_off()`` will never be called if there is an exception.  If you
need to *raise* an exception yourself, call a ``sig_off()`` before
it::

    def raising_an_exception():
        sig_on()
        # (some long computation)
        if (something_failed):
            sig_off()
            raise RuntimeError("something failed")
        # (some more computation)
        sig_off()
        return something

Alternatively, you can use ``try``/``finally`` which will also catch
exceptions raised by subroutines inside the ``try``::

    def try_finally_example():
        sig_on()
        try:
            # (some long computation, potentially raising exceptions)
        finally:
            sig_off()
        return something


Other Signals
-------------

Apart from handling interrupts, ``sig_on()`` provides more general
signal handling.  Indeed, if the code inside ``sig_on()`` would
generate a segmentation fault or call the C function ``abort()`` (or
more generally, raise any of SIGSEGV, SIGILL, SIGABRT, SIGFPE,
SIGBUS), this is caught by the interrupt framework and a
``RuntimeError`` is raised::

    cdef extern from 'stdlib.h':
        void abort()

    def abort_example():
        sig_on()
        abort()
        sig_off()

.. code-block:: python

    sage: abort_example()
    Traceback (most recent call last):
    ...
    RuntimeError: Aborted

This exception can then be caught as explained above.  This means that
``abort()`` can be used as an alternative to exceptions within
``sig_on()``/``sig_off()``.  A segmentation fault unguarded by
``sig_on()`` would simply terminate Sage.

Instead of ``sig_on()``, there is also a function ``sig_str(s)``,
which takes a C string ``s`` as argument.  It behaves the same as
``sig_on()``, except that the string ``s`` will be used as a string
for the exception.  ``sig_str(s)`` should still be closed by
``sig_off()``.  Example Cython code::

    cdef extern from 'stdlib.h':
        void abort()

    def abort_example_with_sig_str():
        sig_str("custom error message")
        abort()
        sig_off()

Executing this gives:

.. code-block:: python

    sage: abort_example_with_sig_str()
    Traceback (most recent call last):
    ...
    RuntimeError: custom error message

With regard to ordinary interrupts (i.e. SIGINT), ``sig_str(s)``
behaves the same as ``sig_on()``: a simple ``KeyboardInterrupt`` is
raised.


.. _advanced-sig:

Advanced Functions
------------------

There are several more specialized functions for dealing with
interrupts.  The function ``sig_check()`` behaves exactly as
``sig_on(); sig_off()`` (except that ``sig_check()`` is faster since
it does not involve a ``setjmp()`` call).

``sig_check()`` can be used to check for pending interrupts.  If an
interrupt happens outside of a ``sig_on()``/``sig_off()`` block, it
will be caught by the next ``sig_check()`` or ``sig_on()``.

The typical use case for ``sig_check()`` is within tight loops doing
complicated stuff (mixed Python and Cython code, potentially raising
exceptions).  It gives more control, because a ``KeyboardInterrupt``
can *only* be raised during ``sig_check()``::

    def sig_check_example():
        for x in foo:
            # (one loop iteration which does not take a long time)
            sig_check()

As mentioned above, ``sig_on()`` makes no attempt to clean anything up
(restore state or freeing memory) when an interrupt occurs.  In fact,
it would be impossible for ``sig_on()`` to do that.  If you want to
add some cleanup code, use ``sig_on_no_except()`` for this. This
function behaves *exactly* like ``sig_on()``, except that any
exception raised (either ``KeyboardInterrupt`` or ``RuntimeError``) is
not yet passed to Python. Essentially, the exception is there, but we
prevent Cython from looking for the exception.  Then
``cython_check_exception()`` can be used to make Cython look for the
exception.

Normally, ``sig_on_no_except()`` returns 1.  If a signal was caught
and an exception raised, ``sig_on_no_except()`` instead returns 0.
The following example shows how to use ``sig_on_no_except()``::

    def no_except_example():
        if not sig_on_no_except():
            # (clean up messed up internal state)

            # Make Cython realize that there is an exception.
            # It will look like the exception was actually raised
            # by cython_check_exception().
            cython_check_exception()
        # (some long computation, messing up internal state of objects)
        sig_off()

There is also a function ``sig_str_no_except(s)`` which is analogous
to ``sig_str(s)``.

.. NOTE::

    See the file :file:`SAGE_ROOT/src/sage/tests/interrupt.pyx`
    for more examples of how to use the various ``sig_*()`` functions.

Testing Interrupts
------------------

.. highlight:: python

When writing :ref:`section-docstrings`,
one sometimes wants to check that certain code can be interrupted in a clean way.
In the module ``sage.tests.interrupt``, there is a function
``interrupt_after_delay(ms_delay = 500)`` which can be used to test interrupts.
That function simulates a ``CTRL-C`` (by sending SIGINT)
after ``ms_delay`` milliseconds.

The following is an example of a doctest demonstrating that
the function ``factor()`` can be interrupted::

    sage: import sage.tests.interrupt
    sage: try:
    ...     sage.tests.interrupt.interrupt_after_delay()
    ...     factor(10^1000 + 3)
    ... except KeyboardInterrupt:
    ...     print "ok!"
    ok!


Unpickling Cython Code
======================

Pickling for python classes and extension classes, such as cython, is
different.  This is discussed in the `python pickling
documentation`_. For the unpickling of extension classes you need to
write a :meth:`__reduce__` method which typically returns a tuple
``(f, args,...)`` such that ``f(*args)`` returns (a copy of) the
original object. As an example, the following code snippet is the
:meth:`~sage.rings.integer.Integer.__reduce__` method from
:class:`sage.rings.integer.Integer`::

    def __reduce__(self):
        '''
        This is used when pickling integers.

        EXAMPLES::

            sage: n = 5
            sage: t = n.__reduce__(); t
            (<built-in function make_integer>, ('5',))
            sage: t[0](*t[1])
            5
            sage: loads(dumps(n)) == n
            True
        '''
        # This single line below took me HOURS to figure out.
        # It is the *trick* needed to pickle Cython extension types.
        # The trick is that you must put a pure Python function
        # as the first argument, and that function must return
        # the result of unpickling with the argument in the second
        # tuple as input. All kinds of problems happen
        # if we don't do this.
        return sage.rings.integer.make_integer, (self.str(32),)


.. _python pickling documentation: http://docs.python.org/library/pickle.html#pickle-protocol

