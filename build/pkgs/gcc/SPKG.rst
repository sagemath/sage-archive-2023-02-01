gcc: The GNU Compiler Collection or other suitable C and C++ compilers
======================================================================

Description
-----------

This package represents the required C and C++ compilers.

- GCC (GNU Compiler Collection) versions 8.x to 12.x are supported.

- Clang (LLVM) is also supported.

The required Fortran compiler is represented by the package ``gfortran``.

You can pass the names of compilers to use to ``./configure`` using
the environment variables :envvar:`CC`, :envvar:`CXX`, and
:envvar:`FC`, for C, C++, and Fortran compilers, respectively.

For example, if your C compiler is ``clang``, your C++ compiler is
``clang++``, and your Fortran compiler is ``flang``, then you would
need to run::

    $ ./configure CC=clang CXX=clang++ FC=flang

Vendor and versions of the C and C++ compilers should match.

This package uses the non-standard default
``configure --with-system-gcc=force``, giving an error at ``configure``
time when no suitable system compilers are configured.

You can override this using ``./configure --without-system-gcc``.  In
this case, Sage builds and installs the GNU Compiler Collection,
including the C, C++ and Fortran compiler. This is not recommended.
You will need suitable C and C++ compilers from which GCC can
bootstrap itself. There are some known problems with old assemblers,
in particular when building the ``ecm`` and ``fflas_ffpack``
packages. You should ensure that your assembler understands all
instructions for your processor. On Linux, this means you need a
recent version of ``binutils`` (not provided by an SPKG); on macOS
you need a recent version of Xcode.

(Installing the
``gfortran`` SPKG becomes a no-op in this case.)

License
-------

GPL version 2 or version 3


Upstream Contact
----------------

https://gcc.gnu.org/
