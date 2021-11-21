
.. _chapter-modularization:

===============================================
 Packaging Portions of the Sage Python Library
===============================================


Modules, packages, distribution packages
========================================

The Sage Python library consists of a large number of Python modules,
organized into a hierarchical set of packages that fill the namespace
``sage``.  All source files are located in a subdirectory of the
directory ``SAGE_ROOT/src/sage/``.

For example,

- the file ``SAGE_ROOT/src/sage/coding/code_bounds.py`` provides the
  module ``sage.coding.code_bounds``;

- the directory containing this file, ``SAGE_ROOT/src/sage/coding/``,
  thus provides the package ``sage.coding``.  This directory contains
  an ``__init__.py`` file, which makes this package an "ordinary"
  package.

There is another notion of "package" in Python, the "distribution
package" (also known as a "distribution" or a "pip-installable
package").  Currently, the entire Sage Python library is provided by a
single distribution, https://pypi.org/project/sagemath-standard/,
which is generated from the directory
``SAGE_ROOT/pkgs/sagemath-standard``.

Note that the distribution name is not required to be a Python
identifier. In fact, using dashes (``-``) is preferred to underscores in
distribution names; **setuptools** and other parts of Python's packaging
infrastructure normalize underscores to dashes. (Using dots in
distribution names, to indicate ownership by organizations, still
mentioned in https://www.python.org/dev/peps/pep-0423/, appears to
have largely fallen out of favor.)

Ordinary vs. implicit namespace packages
----------------------------------------

Each module of the Sage library must be packaged in exactly one
distribution package.

An important constraint is that an "ordinary" package (directory with
``__init__.py`` file) cannot be split into more than one distribution
package.

By removing the ``__init__.py`` file, however, we can make
the package an "implicit" (or "native") "namespace" package, following
https://www.python.org/dev/peps/pep-0420/; in this case,
``__init__.py`` cannot be used any more for initializing the package.

Whenever there are two distribution packages that provide modules with
a common prefix of Python packages, that prefix needs to be a native
namespace package, i.e., there cannot be an ``__init__.py`` file.
For example,

- **sagemath-tdlib** will provide ``sage.graphs.graph_decompositions.tdlib``,

- **sagemath-rw** will provide ``sage.graphs.graph_decompositions.rankwidth``,

- **sagemath-graphs** will provide all of the rest of
  ``sage.graphs.graph_decompositions`` (and most of ``sage.graphs``).

Then, none of

- ``sage``,

- ``sage.graphs``,

- ``sage.graphs.graph_decomposition``

can be an ordinary package (with an ``__init__.py`` file), but rather
each of them has to be an implicit namespace package (no
``__init__.py`` file).


In the Sage 9.6 development cycle, we still use ordinary packages by
default, but several packages are converted to implicit namespace
packages to support modularization.



Source directories of distribution packages
===========================================

The development of the SageMath library uses a monorepo strategy for
all distribution packages that fill the `sage.*` namespace.  This
means that the source trees of these distributions are included in a
single ``git`` repository, in a subdirectory of ``SAGE_ROOT/pkgs``.

All these distribution packages have matching version numbers.  From
the viewpoint of a single distribution, this means that sometimes
there will be a new release of some distribution where the only thing
changing is the version number.

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

- ``setup.py`` -- a **setuptools**-based installation script

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

*Declaring build-time dependencies:* Modern Python packaging provides a
mechanism to declare build-time dependencies on other distribution
packages via the file ``pyproject.toml`` ("build-system requires"); this
has superseded the older ``setup_requires`` declaration. (There is no
mechanism to declare anything regarding the C/C++ libraries.)

While the namespace ``sage.*`` is organized roughly according to
mathematical fields or categories, how we partition the implementation
modules into distribution packages has to respect the hard constraints
that are imposed by the build-time dependencies.

We can define some meaningful small distributions that just consist of
a single or a few Cython modules. For example, **sagemath-tdlib**
(https://trac.sagemath.org/ticket/29864) would just package the single
Cython module that must be linked with ``tdlib``,
:mod:`sage.graphs.graph_decompositions.tdlib`. Starting with the Sage
9.6 development cycle, as soon as namespace packages are activated, we
can start to create these distributions. This is quite a mechanical
task.

Module-level runtime dependencies
---------------------------------

Any ``import`` statements at the top level of a Python or Cython
module are executed when the module is imported. Hence, the imported
modules must be part of the distribution, or provided by another
distribution -- which then must be declared as a run-time dependency.

*Declaring run-time dependencies:* These dependencies are declared in
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

*Declaring optional run-time dependencies:* It is possible to declare
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

*Declaring doctest-only dependencies:* The ``extra_requires`` mechanism
mentioned above can also be used for this.


Version constraints of dependencies
-----------------------------------

The version information for dependencies comes from the files
``build/pkgs/*/install-requires.txt`` and
``build/pkgs/*/package-version.txt``.  We use the ``m4`` macro
processor to insert the version information in the generated files
``pyproject.toml``, ``setup.cfg``, ``requirements.txt``.


Testing distribution packages
=============================

Of course, we need tools for testing modularized distributions of
portions of the Sage library.

1. Modularized distributions must be testable separately!

2. But we want to keep integration testing with other portions of Sage too!

Preparing doctests
------------------

Enter ``# optional``, the doctest annotation that we also use whenever
an optional package is needed for a particular test.

This mechanism can also be used for making a doctest conditional on
the presence of a portion of the Sage library.  The available tags
take the form of package or module names such as ``sage.combinat``,
``sage.graphs``, ``sage.plot``, ``sage.rings.number_field``,
``sage.rings.real_double``, and ``sage.symbolic``.  They are defined
via "features" in a single file,
``SAGE_ROOT/src/sage/features/sagemath.py``, which also provides the
mapping from features to the distributions providing them (actually,
to SPKG names).  Using this mapping, Sage can issue installation hints
to the user.

For example, the package ``sage.tensor`` is purely algebraic and has
no dependency on symbolics. However, there are a small number of
doctests that depend on the Symbolic Ring for integration
testing. Hence, these doctests are marked ``# optional -
sage.symbolic``.

Testing the distribution in virtual environments with tox
---------------------------------------------------------

So how to test that this works?

Sure, we could go into the installation directory
``SAGE_VENV/lib/python3.9/site-packages/`` and do ``rm -rf
sage/symbolic`` and test that things still work. But that's not a good
way of testing.

Instead, we use a virtual environment in which we only install the
distribution to be tested (and its Python dependencies).

Let's try it out first with the entire Sage library, represented by
the distribution **sagemath-standard**.  Note that after Sage has been
built normally, a set of wheels for all installed Python packages is
available in ``SAGE_VENV/var/lib/sage/wheels/``::

  $ ls venv/var/lib/sage/wheels
  Babel-2.9.1-py2.py3-none-any.whl
  Cython-0.29.24-cp39-cp39-macosx_11_0_x86_64.whl
  Jinja2-2.11.2-py2.py3-none-any.whl
  ...
  sage_conf-9.5b6-py3-none-any.whl
  ...
  scipy-1.7.2-cp39-cp39-macosx_11_0_x86_64.whl
  setuptools-58.2.0-py3-none-any.whl
  ...
  wheel-0.37.0-py2.py3-none-any.whl
  widgetsnbextension-3.5.1-py2.py3-none-any.whl
  zipp-3.5.0-py3-none-any.whl

Note in particular the wheel for **sage-conf**, which provides
configuration variable settings and the connection to the non-Python
packages installed in ``SAGE_LOCAL``.

We can now set up a separate virtual environment, in which we install
these wheels and our distribution to be tested.  This is where ``tox``
comes into play: It is the standard Python tool for creating
disposable virtual environments for testing.  Every distribution in
``SAGE_ROOT/pkgs/`` provides a configuration file ``tox.ini``.

Following the comments in the file
``SAGE_ROOT/pkgs/sagemath-standard/tox.ini``, we can try the following
command::

  $ ./sage -sh -c '(cd pkgs/sagemath-standard && SAGE_NUM_THREADS=16 tox -v -v -v -e py39-sagewheels-nopypi)'

This command does not make any changes to the normal installation of
Sage. The virtual environment is created in a subdirectory of
``pkgs/sagemath-standard-no-symbolics/.tox/``. After the command
finishes, we can start the separate installation of the Sage library
in its virtual environment::

  $ pkgs/sagemath-standard/.tox/py39-sagewheels-nopypi/bin/sage

We can also run parts of the testsuite::

  $ pkgs/sagemath-standard/.tox/py39-sagewheels-nopypi/bin/sage -tp 4 src/sage/graphs/

The whole ``.tox`` directory can be safely deleted at any time.
