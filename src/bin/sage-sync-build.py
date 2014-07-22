#!/usr/bin/env python
"""
Author: Craig Citro
        R. Andrew Ohana
Date: 2013 April 10

This script will walk through the output trees of the sage
library and delete any files for which it cannot find a corresponding
source file in $SAGE_SRC/sage. If any directory is empty
after this process, we delete it.

Assumptions made by this code:

  * we want to be aggressive, i.e. delete anything that we don't
    recognize
"""

from __future__ import print_function

import os, sys, getopt, glob
try:
    from sage.env import SAGE_SRC, SAGE_LIB
except ImportError:
    SAGE_SRC = os.environ.get('SAGE_SRC',
            os.path.join(os.environ['SAGE_ROOT'], 'src'))
    SAGE_LIB = os.path.join(
            os.environ.get('SAGE_LOCAL',
                os.path.join(os.environ['SAGE_ROOT'], 'local')),
            'lib', 'python', 'site-packages')

dot = os.path.extsep
out_dirs = [SAGE_LIB]+glob.glob(os.path.join(SAGE_SRC, 'build', '*'))

def check_file(relative_path, ext_module_bases):
    """
    Takes a path relative to $SAGE_SRC, and decides whether or not the
    corresponding source file still exists.

    Returns True if the corresponding source file still exists,
    and False otherwise.
    """
    # look for the appropriate source file.
    relative_base, build_ext = os.path.splitext(relative_path)
    if build_ext in (dot+'o', dot+'so'):
        return relative_base in ext_module_bases
    elif build_ext == dot+'pyc':
        source_ext = dot+'py'
    else:
        source_ext = build_ext
    return os.path.exists(os.path.join(SAGE_SRC, relative_base+source_ext))

def clean_tree(prune_directories=True, dry_run=False):
    os.chdir(SAGE_SRC)

    # we need to monkey with the path a bit to load the
    # appropriate list of extension modules
    sys.path.append(os.getcwd())
    from module_list import ext_modules
    from sage.ext.gen_interpreters import modules as gen_modules
    sys.path.pop()

    # finally build the list of extension modules
    ext_module_bases = set()
    for x in ext_modules + gen_modules:
        for y in x.sources:
            y_base, y_ext = os.path.splitext(y)
            ext_module_bases.add(y_base)

    # start walking the directory tree. we explicitly set
    # topdown=False, so that we can prune empty directories as we
    # encounter them.
    def walk():
        for dir in out_dirs:
            for path, dirs, files in os.walk(
                    os.path.join(dir, 'sage'), topdown=False):
                yield dir, path, dirs, files
    for build_dir, path, dirs, files in walk():
        for filename in files:
            file_path = os.path.join(path, filename)
            found = check_file(os.path.relpath(file_path, build_dir),
                    ext_module_bases=ext_module_bases)
            if not found:
                if dry_run:
                    print("Remove file: %s"%file_path)
                else:
                    os.remove(file_path)
        for directory in dirs:
            dir_path = os.path.join(path, directory)
            file_ls = os.listdir(dir_path)
            if prune_directories and not file_ls:
                if dry_run:
                    print("Remove empty directory: %s"%dir_path)
                else:
                    os.rmdir(dir_path)
            elif not file_ls:
                print("Remove empty directory: %s"%dir_path)

dry_run = False
prune_directories = True

# process command-line arguments
args, extra_args = getopt.getopt(sys.argv[1:], 'dp', ['dry-run', 'prune-directories'])

for arg, opt in args:
    if arg in ['-d', '--dry-run']:
        dry_run = True
    elif arg in ['-p', '--prune-directories']:
        prune_directories = False

clean_tree(prune_directories=prune_directories, dry_run=dry_run)
