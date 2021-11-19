
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
