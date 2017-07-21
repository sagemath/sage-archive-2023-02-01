"""
Command-line script for uninstalling an existing Sage spkg from $SAGE_LOCAL.

This performs two types of uninstallation:

    1) Old-style uninstallation: This is close to what existed before this
       script, where *some* packages had uninstall steps (mostly consisting of
       some broad `rm -rf` commands) that were run before installing new
       versions of the packages.  This convention was applied inconsistently,
       but for those packages that did have old-style uninstall steps, those
       steps should be in a script called `spkg-legacy-uninstall` under the
       spkg directory (build/pkgs/<pkg_name>).  If this script is found, it is
       run for backwards-compatibility support.

    2) New-style uninstallation: More recently installed packages that were
       installed with staged installation have a record of all files installed
       by that package.  That file is stored in the $SAGE_SPKG_INST directory
       (typically $SAGE_LOCAL/var/lib/sage/installed) and is created when the
       spkg is installed.  This is a JSON file containing some meta-data about
       the package, including the list of all files associated with the
       package.  This script removes all these files, including the record
       file.  Any directories that are empty after files are removed from them
       are also removed.
"""

#*****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.m.bray@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

import glob
import json
import os
import shutil
import subprocess
import sys

from os import path as pth

try:
    import argparse
except ImportError:
    from sage_bootstrap.compat import argparse

from .env import SAGE_ROOT


PKGS = pth.join(SAGE_ROOT, 'build', 'pkgs')
"""Directory where all spkg sources are found."""


def uninstall(spkg_name, sage_local):
    """
    Given a package name and path to SAGE_LOCAL, uninstall that package from
    SAGE_LOCAL if it is currently installed.
    """

    # TODO: I don't like that this is hard-coded here; there should be better
    # centralization for this without having to load sage-env or sage.env;
    # https://trac.sagemath.org/ticket/22652 would help here
    spkg_inst = pth.join(sage_local, 'var', 'lib', 'sage', 'installed')

    # Find all stamp files for the package; there should be only one, but if
    # there is somehow more than one we'll work with the most recent and delete
    # the rest
    stamp_files = glob.glob(pth.join(spkg_inst, '{0}-*'.format(spkg_name)))

    if len(stamp_files) > 1:
        stamp_files.sort(key=pth.getmtime)
    elif not stamp_files:
        print("No record that '{0}' was ever installed; skipping "
              "uninstall".format(spkg_name), file=sys.stderr)
        return

    stamp_file = stamp_files[-1]

    try:
        with open(stamp_file) as f:
            spkg_meta = json.load(f)
    except (OSError, ValueError):
        spkg_meta = {}

    if 'files' not in spkg_meta:
        print("Old-style or corrupt stamp file '{0}'".format(stamp_file),
              file=sys.stderr)

        legacy_uninstall(spkg_name)
    else:
        files = spkg_meta['files']
        if not files:
            print("Warning: No files to uninstall for "
                  "'{0}'".format(spkg_name), file=sys.stderr)

        modern_uninstall(spkg_name, sage_local, files)

    # Finally, if all went well, delete all the stamp files.
    for stamp_file in stamp_files:
        os.remove(stamp_file)


def legacy_uninstall(spkg_name):
    """
    Run the spkg's legacy uninstall script, if one exists; otherwise do
    nothing.
    """

    spkg_dir = pth.join(PKGS, spkg_name)
    legacy_uninstall = pth.join(spkg_dir, 'spkg-legacy-uninstall')

    if not (pth.isfile(legacy_uninstall) and
            os.access(legacy_uninstall, os.X_OK)):
        print("No legacy uninstaller found for '{0}'; nothing to "
              "do".format(spkg_name), file=sys.stderr)
        return

    print("Uninstalling '{0}' with legacy uninstaller".format(spkg_name),
          file=sys.stderr)

    # Any errors from this, including a non-zero return code will
    # bubble up and exit the uninstaller
    subprocess.check_call([legacy_uninstall])


def modern_uninstall(spkg_name, sage_local, files):
    """
    Remove all listed files from the given SAGE_LOCAL (all file paths should
    be assumed relative to SAGE_LOCAL).

    This is otherwise (currently) agnostic about what package is actually
    being uninstalled--all it cares about is removing a list of files.

    If the directory containing any of the listed files is empty after all
    files are removed then the directory is removed as well.
    """

    spkg_scripts = pth.join(sage_local, 'var', 'lib', 'sage', 'scripts',
                            spkg_name)

    # Sort the given files first by the highest directory depth, then by name,
    # so that we can easily remove a directory once it's been emptied
    files.sort(key=lambda f: (-f.count(os.sep), f))
    cur_dir = None

    print("Uninstalling existing '{0}'".format(spkg_name), file=sys.stderr)

    for filename in files:
        filename = pth.join(sage_local, filename)
        dirname = pth.dirname(filename)

        if cur_dir != dirname:
            cur_dir = dirname

            # New directory; see if the previous directory has been emptied out
            if cur_dir and not os.listdir(cur_dir):
                os.rmdir(cur_dir)

        if os.path.exists(filename):
            os.remove(filename)
        else:
            print("Warning: File '{0}' not found".format(filename),
                  file=sys.stderr)

    # Remove the last directory, if it was empty
    if cur_dir and not os.listdir(cur_dir):
        os.rmdir(cur_dir)

    # Run the package's postrm script, if it exists
    postrm = pth.join(spkg_scripts, 'spkg-postrm')
    if pth.exists(postrm):
        print("Running post-uninstall script for '{0}'".format(spkg_name),
              file=sys.stderr)
        # Any errors from this, including a non-zero return code will
        # bubble up and exit the uninstaller
        subprocess.check_call([postrm])
        shutil.rmtree(spkg_scripts)


def dir_type(path):
    """
    A custom argument 'type' for directory paths.
    """

    if path and not pth.isdir(path):
        raise argparse.ArgumentError(
            "'{0}' is not a directory".format(path))

    return path


def spkg_type(pkg):
    """
    A custom argument 'type' for spkgs--checks whether the given package name
    is a known spkg.
    """

    pkgbase = pth.join(PKGS, pkg)

    if not pth.isdir(pkgbase):
        raise argparse.ArgumentError(
                "'{0}' is not a known spkg".format(pkg))

    return pkg


def make_parser():
    """Returns the command-line argument parser for sage-spkg-uninstall."""

    doc_lines = __doc__.strip().splitlines()

    parser = argparse.ArgumentParser(
            description=doc_lines[0],
            epilog='\n'.join(doc_lines[1:]).strip(),
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('spkg', type=spkg_type, help='the spkg to uninstall')
    parser.add_argument('sage_local', type=dir_type, nargs='?',
                        default=os.environ.get('SAGE_LOCAL'),
                        help='the SAGE_LOCAL path (default: the $SAGE_LOCAL '
                             'environment variable if set)')
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
        uninstall(args.spkg, args.sage_local)
    except Exception as exc:
        print("Error during uninstallation of '{0}': {1}".format(
            args.spkg, exc), file=sys.stderr)

        if args.debug:
            raise

        sys.exit(1)


if __name__ == '__main__':
    run()
