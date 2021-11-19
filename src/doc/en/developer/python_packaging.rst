
.. _chapter-modularization:

===============================================
 Packaging Portions of the Sage Python Library
===============================================


Modules, packages, distribution packages
========================================

The Sage Python library consists of a large number of Python modules,
organized into a hierarchical set of packages that fill the namespace
`sage`.  All source files are located in a subdirectory of the
directory ``SAGE_ROOT/src/sage/``.

For example,

- the file ``SAGE_ROOT/src/sage/coding/code_bounds.py`` provides the
  module ``sage.coding.code_bounds``;

- the directory containing this file, ``SAGE_ROOT/src/sage/coding/``,
  thus provides the package ``sage.coding``.  This directory contains
  an ``__init__.py`` file, which makes this package an "ordinary"
  package.  (Later, we will also discuss package directories without
  ``__init__.py`` files, which provide "implicit" (or "native")
  "namespace" packages.)

There is another notion of "package" in Python, the "distribution
package" (also known as a "distribution" or a "pip-installable
package").  Currently, the entire Sage Python library is provided by a
single distribution, https://pypi.org/project/sagemath-standard/,
which is generated from the directory
``SAGE_ROOT/pkgs/sagemath-standard``.

Note that the distribution name is not required to be a Python
identifier. In fact, using dashes (``-``) is preferred to underscores in
distribution names; ``setuptools`` and other parts of Python's packaging
infrastructure normalize underscores to dashes. (Using dots in
distribution names, to indicate ownership by organizations, still
mentioned in https://www.python.org/dev/peps/pep-0423/, appears to
have largely fallen out of favor.)


Source directories of distribution packages
===========================================

The source directory of a distribution package, such as
``SAGE_ROOT/pkgs/sagemath-standard``, contains the following files:

- ``sage`` -- a relative symbolic link to the monolithic Sage library
  source tree ``SAGE_ROOT/src/sage/``

- ``MANIFEST.in`` -- controls which files and directories of the
  monolithic Sage library source tree are included in the distribution

- ``pyproject.toml.m4``, ``setup.cfg.m4``, ``requirements.txt.m4`` --
  Python packaging metadata, declaring the distribution name, dependencies,
  etc.  (These files are run through the ``m4`` macro processor by the
  ``SAGE_ROOT/bootstrap`` script to generate standard files
  ``pyproject.toml`` etc.)

- ``setup.py`` -- a ``setuptools``-based installation script

- ``tox.ini`` -- testing infrastructure


Dependencies and distribution packages
======================================

When preparing a portion of the Sage library as a distribution
package, dependencies matter.

Build-time dependencies
-----------------------

If the portion of the library contains any Cython modules, these
modules are compiled during the wheel-building phase of the
distribution package. If the Cython module uses ``cimport`` to pull in
anything from ``.pxd`` files, these files must be either part of the
portion shipped as the distribution being built, or the distribution
that provides these files must be installed in the build
environment. Also, any C/C++ libraries that the Cython module uses
must be accessible from the build environment.

Declaring build-time dependencies: Modern Python packaging provides a
mechanism to declare build-time dependencies on other distribution
packages via the file ``pyproject.toml`` ("build-system requires"); this
has superseded the older ``setup_requires`` declaration. (There is no
mechanism to declare anything regarding the C/C++ libraries.)

Module-level runtime dependencies
---------------------------------

Any ``import`` statements at the top level of a Python or Cython
module are executed when the module is imported. Hence, the imported
modules must be part of the distribution, or provided by another
distribution -- which then must be declared as a run-time dependency.

Declaring run-time dependencies: These dependencies are declared in
``setup.cfg`` (generated from ``setup.cfg.m4``) as ``install_requires``.

Other runtime dependencies
--------------------------

If ``import`` statements are used within a method, the imported module
is loaded the first time that the method is called. Hence the module
defining the method can still be imported even if the module needed by
the method is not present.

It is then a question whether a run-time dependency should be
declared. If the method needing that import provides core
functionality, then probably yes. But if it only provides what can be
considered "optional functionality", then probably not, and in this
case it will be up to the user to install the distribution enabling
this optional functionality.

Declaring optional run-time dependencies: It is possible to declare
such optional dependencies as ``extra_requires`` in ``setup.cfg``
(generated from ``setup.cfg.m4``).  This is a very limited mechanism
-- in particular it does not affect the build phase of the
distribution in any way. It basically only provides a way to give a
nickname to a distribution that can be installed as an add-on. It
allows users to write, for example, ``pip install
sagemath-polyhedra[normaliz]`` instead of ``pip install
sagemath-polyhedra pynormaliz``.

Lazy module-level imports
-------------------------

Sage provides the ``lazy_import`` mechanism. Lazy imports can be
declared at the module level, but the actual importing is only done on
demand. It is a runtime error at that time if the imported module is
not present. This can be convenient compared to local imports in
methods when the same imports are needed in several methods.

Doctest-only dependencies
-------------------------

Doctests often use examples constructed using functionality provided
by other portions of the Sage library.  This kind of integration
testing as one of the strengths of Sage; but it also creates extra
dependencies.

Fortunately, these dependencies are very mild, and we can deal with
them using the same mechanism that we use for making doctests
conditional on the presence of optional libraries: using ``# optional -
FEATURE`` directives in the doctests.  Adding these directives will
allow developers to test the distribution separately, without
requiring all of Sage to be present.

Declaring doctest-only dependencies: The ``extra_requires`` mechanism
mentioned above can also be used for this.
