"""
Command-line script for checking an installed SPKG in $SAGE_LOCAL.
"""

# ****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.m.bray@gmail.com>
#                     2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function

import glob
import json
import os
import shutil
import subprocess
import sys
import argparse
from warnings import warn

from .env import SAGE_ROOT

pth = os.path
PKGS = pth.join(SAGE_ROOT, 'build', 'pkgs')
"""Directory where all spkg sources are found."""


def installcheck(spkg_name, sage_local, verbose=False):
    """
    Given a package name and path to SAGE_LOCAL, check the installation of the package
    in SAGE_LOCAL.
    """

    # The default path to this directory; however its value should be read
    # from the environment if possible
    spkg_inst = pth.join(sage_local, 'var', 'lib', 'sage', 'installed')

    # Find all stamp files for the package; there should be only one, but if
    # there is somehow more than one we'll work with the most recent and delete
    # the rest
    pattern = pth.join(spkg_inst, '{0}-*'.format(spkg_name))
    stamp_files = sorted(glob.glob(pattern), key=pth.getmtime)

    if stamp_files:
        stamp_file = stamp_files[-1]
    else:
        stamp_file = None

    spkg_meta = {}
    if stamp_file:
        try:
            with open(stamp_file) as f:
                spkg_meta = json.load(f)
        except (OSError, ValueError):
            pass

    if 'files' not in spkg_meta:
        if stamp_file:
            print("Old-style or corrupt stamp file '{0}'"
                  .format(stamp_file), file=sys.stderr)
        else:
            print("Package '{0}' is currently not installed in '{1}'"
                  .format(spkg_name, sage_local), file=sys.stderr)
    else:
        files = spkg_meta['files']

        for f in files:
            f = os.path.join(sage_local, f)
            if f.endswith(('.so', '.dylib')):
                if verbose:
                    print("Checking shared library file '{0}'"
                          .format(f), file=sys.stderr)
                try:
                    from delocate.libsana import _tree_libs_from_libraries, _filter_system_libs
                except ImportError:
                    warnings.warn('delocate is not available, so nothing is actually checked')
                else:
                    _tree_libs_from_libraries([f],
                                              lib_filt_func=_filter_system_libs,
                                              copy_filt_func=lambda path: True)
            elif f.endswith('-any.whl'):
                # pure Python wheel, nothing to check
                pass
            elif f.endswith('.whl'):
                if verbose:
                    print("Checking wheel file '{0}'"
                          .format(f), file=sys.stderr)
                try:
                    from delocate import wheel_libs
                except ImportError:
                    warnings.warn('delocate is not available, so nothing is actually checked')
                else:
                    wheel_libs(f)


def dir_type(path):
    """
    A custom argument 'type' for directory paths.
    """

    if path and not pth.isdir(path):
        raise argparse.ArgumentTypeError(
            "'{0}' is not a directory".format(path))

    return path


def spkg_type(pkg):
    """
    A custom argument 'type' for spkgs--checks whether the given package name
    is a known spkg.
    """
    pkgbase = pth.join(PKGS, pkg)

    if not pth.isdir(pkgbase):
        raise argparse.ArgumentTypeError(
            "'{0}' is not a known spkg".format(pkg))

    return pkg


def make_parser():
    """Returns the command-line argument parser for sage-spkg-installcheck."""

    doc_lines = __doc__.strip().splitlines()

    parser = argparse.ArgumentParser(
        description=doc_lines[0],
        epilog='\n'.join(doc_lines[1:]).strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('spkg', type=spkg_type, help='the spkg to check')
    parser.add_argument('sage_local', type=dir_type, nargs='?',
                        default=os.environ.get('SAGE_LOCAL'),
                        help='the SAGE_LOCAL path (default: the $SAGE_LOCAL '
                             'environment variable if set)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose output showing all files removed')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)

    return parser


def run(argv=None):
    parser = make_parser()

    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    if args.sage_local is None:
        print('Error: SAGE_LOCAL must be specified either at the command '
              'line or in the $SAGE_LOCAL environment variable',
              file=sys.stderr)
        sys.exit(1)

    try:
        installcheck(args.spkg, args.sage_local,
                  verbose=args.verbose)
    except Exception as exc:
        print("Error during installcheck of '{0}': {1}".format(
            args.spkg, exc), file=sys.stderr)

        if args.debug:
            raise

        sys.exit(1)


if __name__ == '__main__':
    run()
