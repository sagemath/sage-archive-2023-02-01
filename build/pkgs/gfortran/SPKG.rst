gfortran: Fortran compiler from the GNU Compiler Collection
===========================================================

Description
-----------

This package represents the required Fortran compiler.

Officially we support ``gfortran`` from `GNU Compiler Collection (GCC)
<https://gcc.gnu.org/>`_.  It has also been reported that using ``flang``
(from LLVM) might work.

You can pass the names of compilers to use to ``./configure`` using
the environment variables :envvar:`CC`, :envvar:`CXX`, and
:envvar:`FC`, for C, C++, and Fortran compilers, respectively.

For example, if your C compiler is ``clang``, your C++ compiler is
``clang++``, and your Fortran compiler is ``flang``, then you would
need to run::

    $ ./configure CC=clang CXX=clang++ FC=flang

Building Sage from source on Apple Silicon (M1/M2) requires the use of
​the `Homebrew package manager <https://brew.sh>`_ (recommended) or
conda-forge, which package versions of GCC 12.x (including
``gfortran``) with the necessary changes for this platform.  These
changes are not in a released upstream version of GCC, and hence
the ``gfortran`` SPKG is not suitable for the M1/M2.

License
-------

GPL version 2 or version 3


Upstream Contact
----------------

http://gcc.gnu.org/

Special Update/Build Instructions
---------------------------------

None.
