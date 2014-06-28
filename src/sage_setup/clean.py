"""
Clean the Install Dir
"""
#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os

from sage_setup.find import installed_files_by_module



def _remove(file_set, module_base, to_remove):
    """
    Helper to remove files from a set of filenames.

    INPUT:

    - ``file_set`` -- a set of filenames.

    - ``module_base`` -- string. Name of a Python package/module.

    - ``to_remove`` -- list/tuple/iterable of strings. Either
      filenames or extensions (starting with ``'.'``)

    OUTPUT:

    This function does not return anything. The ``file_set`` parameter
    is modified in place.

    EXAMPLES::

        sage: files = set(['a/b/c.py', 'a/b/d.py', 'a/b/c.pyx'])
        sage: _remove(files, 'a.b', ['c.py', 'd.py'])
        sage: files
        set(['a/b/c.pyx'])

        sage: files = set(['a/b/c.py', 'a/b/d.py', 'a/b/c.pyx'])
        sage: _remove(files, 'a.b.c', ['.py', '.pyx'])
        sage: files
        set(['a/b/d.py'])
    """
    path = os.path.join(*module_base.split('.'))
    for filename in to_remove:
        if filename.startswith('.'):
            filename = path + filename
        else:
            filename = os.path.join(path, filename)
        file_set.difference_update([filename])


def _find_stale_files(site_packages, python_packages, python_modules, ext_modules):
    """
    Find stale files

    This method lists all files installed and then subtracts the ones
    which are intentionally being installed.

    EXAMPLES:

    It is crucial that only truly stale files are being found, of
    course. We check that when the doctest is being run, that is,
    after installation, there are no stale files::

        sage: from sage.env import SAGE_SRC, SAGE_LIB
        sage: from sage_setup.find import find_python_sources
        sage: python_packages, python_modules = find_python_sources(
        ....:     SAGE_SRC, ['sage', 'sage_setup'])
        sage: from sage_setup.clean import _find_stale_files

    TODO: move ``module_list.py`` into ``sage_setup`` and also check
    extension modules::

        sage: stale_iter = _find_stale_files(SAGE_LIB, python_packages, python_modules, [])
        sage: for f in stale_iter:
        ....:     if f.endswith('.so'): continue
        ....:     print('Found stale file: ' + f)
    """
    PYMOD_EXTS = (os.path.extsep + 'py', os.path.extsep + 'pyc')
    CEXTMOD_EXTS = (os.path.extsep + 'so',)
    INIT_FILES= map(lambda x: '__init__' + x, PYMOD_EXTS)

    module_files = installed_files_by_module(site_packages, ['sage', 'sage_setup'])

    for mod in python_packages:
        try:
            files = module_files[mod]
        except KeyError:
            # the source module "mod" has not been previously installed, fine.
            continue
        _remove(files, mod, INIT_FILES)
    for mod in python_modules:
        try:
            files = module_files[mod]
        except KeyError:
            continue
        _remove(files, mod, PYMOD_EXTS)
    for ext in ext_modules:
        mod = ext.name
        try:
            files = module_files[mod]
        except KeyError:
            continue
        _remove(files, mod, CEXTMOD_EXTS)
    for files in module_files.values():
        for f in files:
            yield f


def clean_install_dir(site_packages, python_packages, python_modules, ext_modules):
    """
    Delete all modules that are **not** being installed

    If you switch branches it is common to (re)move the source for an
    already installed module. Subsequent rebuilds will leave the stale
    module in the install directory, which can break programs that try
    to import all modules. In particular, the Sphinx autodoc builder
    does this and is susceptible to hard-to-reproduce failures that
    way. Hence we must make sure to delete all stale modules.

    INPUT:

    - ``site_packages`` -- the root Python path where the Sage library
      is being installed.

    - ``python_packages`` -- list of pure Python packages (directories
      with ``__init__.py``).

    - ``python_modules`` -- list of pure Python modules.

    - ``ext_modules`` -- list of distutils ``Extension`` classes. The
      output of ``cythonize``.
    """
    stale_file_iter = _find_stale_files(
        site_packages, python_packages, python_modules, ext_modules)
    for f in stale_file_iter:
        f = os.path.join(site_packages, f)
        print('Cleaning up stale file: {0}'.format(f))
        os.unlink(f)
