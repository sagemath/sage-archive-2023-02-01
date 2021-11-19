
.. _chapter-modularization:

===============================================
 Packaging Portions of the Sage Python Library
===============================================

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
  ``__init__.py`` files, implicit (or "native") namespace packages.)

There is another notion of "package" in Python, the "distribution
package" (also known as a "distribution" or a "pip-installable
package").  Currently, the entire Sage Python library is provided by a
single distribution, https://pypi.org/project/sagemath-standard/,
which is generated from the directory
``SAGE_ROOT/pkgs/sagemath-standard``.
