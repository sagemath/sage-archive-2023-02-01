.. nodoctest

.. highlight:: shell-session

.. _chapter-tools:

========================================
Additional development and testing tools
========================================

Tox
===

`Tox <https://tox.readthedocs.io/en/latest/>`_ is a popular package that is
used by a large number of Python projects as the standard entry point
for testing and linting.

Sage includes tox as a standard package and uses it for three purposes:

- For portability testing of the Sage distribution, as we explain in
  :ref:`chapter-portability_testing`.  This is configured in the file
  ``SAGE_ROOT/tox.ini``.

- For testing modularized distributions of the Sage library, as we explain
  in :ref:`chapter-modularization`.  This is configured in ``tox.ini`` files in subdirectories
  of ``SAGE_ROOT/pkgs/``, such as ``SAGE_ROOT/pkgs/sagemath-standard/tox.ini``.

- As an entry point for testing and linting of the Sage library, as we describe below.
  This is configured in the file ``SAGE_ROOT/src/tox.ini``.

The tox configuration ``SAGE_ROOT/src/tox.ini`` can be invoked by using the command
``./sage --tox``.  (If ``tox`` is available in your system installation,
you can just type ``tox`` instead.)

This configuration provides an entry point for various testing/linting methods,
known as "tox environments".  We can type ``./sage --advanced`` so see what is
available::

  $ ./sage --advanced
  SageMath version 9.2
  ...
  Testing files:
  ...
  --tox [options] <files|dirs> -- general entry point for testing
                                  and linting of the Sage library
     -e <envlist>     -- run specific test environments
                         (default: run all except full pycodestyle)
        doctest                -- run the Sage doctester
                                  (same as "sage -t")
        coverage               -- give information about doctest coverage of files
                                  (same as "sage --coverage[all]")
        startuptime            -- display how long each component of Sage takes to start up
                                  (same as "sage --startuptime")
        pycodestyle-minimal    -- check against Sage's minimal style conventions
        relint                 -- check whether some forbidden patterns appear
                                  (includes all patchbot pattern-exclusion plugins)
        codespell              -- check for misspelled words in source code
        pycodestyle            -- check against the Python style conventions of PEP8
     -p auto          -- run test environments in parallel
     --help           -- show tox help

Pycodestyle
===========
`Pycodestyle <https://pycodestyle.pycqa.org/en/latest/>`_ (formerly known as pep8)
checks Python code against the style conventions of `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_.

Sage defines two configurations for pycodestyle.  The command ``./sage -tox -e pycodestyle-minimal`` uses
pycodestyle in a minimal configuration.
As of Sage 9.5, the entire Sage library conforms to this configuration::

  $ ./sage -tox -e pycodestyle-minimal -- src/sage/
  pycodestyle-minimal installed: pycodestyle==2.8.0
  pycodestyle-minimal run-test-pre: PYTHONHASHSEED='28778046'
  pycodestyle-minimal run-test: commands[0] | pycodestyle --select E401,E70,W605,E711,E712,E721 sage
  ___________ summary ____________
    pycodestyle-minimal: commands succeeded
    congratulations :)

The second configuration is used with the command ``./sage -tox -e pycodestyle`` and runs a
more thorough check::

  $ ./sage -tox -e pycodestyle -- src/sage/quadratic_forms/quadratic_form.py
  pycodestyle installed: pycodestyle==2.8.0
  pycodestyle run-test-pre: PYTHONHASHSEED='2520226550'
  pycodestyle run-test: commands[0] | pycodestyle sage/quadratic_forms/quadratic_form.py
  sage/quadratic_forms/quadratic_form.py:135:9: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:163:64: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:165:52: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:173:42: E228 missing whitespace around modulo operator
  ...
  sage/quadratic_forms/quadratic_form.py:1620:9: E266 too many leading '#' for block comment
  sage/quadratic_forms/quadratic_form.py:1621:9: E266 too many leading '#' for block comment
  25      E111 indentation is not a multiple of 4
  2       E117 over-indented
  129     E127 continuation line over-indented for visual indent
  1       E128 continuation line under-indented for visual indent
  4       E201 whitespace after '['
  4       E202 whitespace before ']'
  2       E222 multiple spaces after operator
  7       E225 missing whitespace around operator
  1       E228 missing whitespace around modulo operator
  25      E231 missing whitespace after ','
  1       E262 inline comment should start with '# '
  3       E265 block comment should start with '# '
  62      E266 too many leading '#' for block comment
  2       E272 multiple spaces before keyword
  2       E301 expected 1 blank line, found 0
  17      E303 too many blank lines (2)
  ERROR: InvocationError for command .../pycodestyle sage/quadratic_forms/quadratic_form.py (exited with code 1)
  ___________ summary ____________
  ERROR:   pycodestyle: commands failed

By passing the options ``--count -qq`` we can reduce the output to
only show the number of style violation warnings::

  $ ./sage -tox -e pycodestyle -- --count -qq src/sage
  pycodestyle installed: pycodestyle==2.8.0
  pycodestyle run-test-pre: PYTHONHASHSEED='3166223974'
  pycodestyle run-test: commands[0] | pycodestyle --count -qq sage
  557     E111 indentation is not a multiple of 4
  1       E112 expected an indented block
  194     E114 indentation is not a multiple of 4 (comment)
  ...
  7       E743 ambiguous function definition 'l'
  335     W291 trailing whitespace
  4       W292 no newline at end of file
  229     W293 blank line contains whitespace
  459     W391 blank line at end of file
  97797
  ERROR: InvocationError for command .../pycodestyle --count -qq sage (exited with code 1)
  ___________ summary ____________
  ERROR:   pycodestyle: commands failed

*Installation:* (for manual use:) ``pip install -U pycodestyle --user``

*Usage:*

- With tox: See above.

- Manual: Run ``pycodestyle path/to/the/file.py``.

- VS Code: Activate by adding the setting ``"python.linting.pycodestyleEnabled": true``, see `official VS Code documentation <https://code.visualstudio.com/docs/python/linting>`__ for details.

*Configuration:* ``[pycodestyle]`` block in ``SAGE_ROOT/src/tox.ini``

*Documentation:* https://pycodestyle.pycqa.org/en/latest/index.html



Pytest
======
`Pytest <https://docs.pytest.org/en/stable/>`_ is a testing framework.
It is included in the Sage distribution as an optional package.

Currently, Sage only makes very limited use of pytest, for testing the
package :mod:`sage.numerical.backends`.

*Installation:*

- (for use with the Sage doctester:) ``./sage -i pytest``.

- (for manual use:) ``pip install -U pytest``, see `documentation <https://docs.pytest.org/en/stable/getting-started.html#installation-and-getting-started>`__ for details.

*Usage:*

- Tox, Sage doctester: At the end of ``./sage -t`` (or ``./sage --tox -e doctest``), Pytest is automatically invoked.

- Manual: Run ``pytest path/to/the/test_file.py`` or ``pytest`` to run all tests (from a virtual environment with the Sage library installed)

- VS Code: Install the `Python extension <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`_ and follow the `offical VS Code documentation <https://code.visualstudio.com/docs/python/testing>`__.

*Configuration:* ``SAGE_ROOT/src/conftest.py``

*Documentation:* https://docs.pytest.org/en/stable/index.html

Pyright 
===============================
`Pyright <https://github.com/microsoft/pyright>`_ is static type checker.

*Installation:* ``npm install -g pyright``, see `documentation <https://github.com/microsoft/pyright#installation>`__ for details.

*Usage:*

- Manual: Run ``pyright path/to/the/file.py``

- VS Code: Install the `Pylance <https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance>`__ extension.

*Configuration:* ``SAGE_ROOT/pyrightconfig.json``

*Note*: Currently, only the package :mod:`sage.manifolds` is checked. Further packages can be added in the ``pyrightconfig.json`` file.

*Documentation:* https://github.com/microsoft/pyright#documentation

Pyflakes
===============================
`Pyflakes <https://github.com/PyCQA/pyflakes>`_ checks for common coding errors.
