#!/usr/bin/env python

##
## sync-build.py
##
## Author: Craig Citro
## Date: 2009 May 26
##
## This script will walk through the tree in $SAGE_DEVEL/sage/build/
## and delete any files for which it cannot find a corresponding
## source file in $SAGE_DEVEL/sage/sage/. If any directory is empty
## after this process, we delete it.
##
## Assumptions made by this code:
##
##  * no filename in the sage source tree contains a copy of
##    os.path.extsep (which is the generic name for '.') other than to
##    separate the filename from the extension
##
##  * we want to be aggressive, i.e. delete anything that we don't
##    recognize
##
## CAVEAT: we need to look at the subdirectories in
## $SAGE_DEVEL/sage/build/, and find corresponding source files in
## $SAGE_DEVEL/sage/sage/. We know that distutils makes several copies
## of this in build/, and we need to be able to go back and forth
## between the trees in build/ and sage/. This means we need to make
## some assumption about the filenames; in particular, we assume that
## the first occurrence of the string '/sage' is the root of a copy of
## the sage source tree. This is unlikely to cause any issue, as far
## as I know.

import os, sys, getopt

dot = os.path.extsep
slash = os.path.sep

try:
    sage_root      = os.environ['SAGE_ROOT']
    sage_main_repo = sage_root + slash + 'devel' + slash + 'sage'
    sage_source    = sage_main_repo + slash + 'sage'
    sage_build     = sage_main_repo + slash + 'build'
except KeyError:
    raise RuntimeError, "no SAGE_ROOT defined!"

def check_file(relative_path, filename, special_source_files=None, verbose=False):
    """
    Takes a relative path and a filename, where the path
    is relative to $SAGE_ROOT/devel/sage/build, and
    decides whether or not the corresponding source file
    still exists in $SAGE_ROOT/devel/sage/sage.

    We assume that all .o and .so files are originally created from
    Cython source files, and look for the corresponding .pyx
    file. However, if the file does not have a corresponding .pyx, but
    a source file with the same name occurs in special_source_files,
    we return True.

    Returns True if the corresponding source file still exists,
    and False otherwise. If verbose is True, it explains what
    file it couldn't find when returning False.
    """
    # set the stub of the path to the source file
    source_path = relative_path[relative_path.find(slash + 'sage')+1:]
    file_to_check = sage_main_repo + slash + source_path

    # process the special source files
    files_to_ignore = []
    for f in special_source_files:
        p, n = os.path.split(f)
        files_to_ignore.append(sage_main_repo + slash + p + slash + n[:n.find(dot)])

    # look for the appropriate source file.
    if filename.endswith(dot + 'pyc'):
        file_to_check += filename[:-1]
        if not os.path.exists(file_to_check):
            if verbose:
                print "Cannot find source file %s for file %s."%(file_to_check, relative_path + filename)
            return False
        else:
            return True
    elif filename.endswith(dot + 'py') or filename.endswith(dot + 'pyx'):
        file_to_check += filename
        if not os.path.exists(file_to_check):
            if verbose:
                print "Cannot find source file %s for file %s."%(file_to_check, relative_path + filename)
            return False
        else:
            return True
    elif filename.endswith(dot + 'o') or filename.endswith(dot + 'so'):
        file_to_check += filename.split(dot)[0]

        # first, look for a cython source file
        if os.path.exists(file_to_check + dot + 'pyx'):
            return True
        # next, try for one of the special source files
        elif (file_to_check in files_to_ignore):
            return True
        else:
            if verbose:
                print
                print file_to_check
                print "Cannot find source file for file %s."%(relative_path + filename)
            return False

    else:
        # If we wanted to be more conservative, we could raise an error here
        # instead of simply deciding the file needs to be deleted.
        #raise RuntimeError, "unknown file type for file %s"%(relative_path + filename)

        if verbose:
            print "Could not determine file type for file %s"%(relative_path + filename)
        return False


def clean_tree(prune_directories=True, dry_run=False, verbose=False):
    # remember where we started
    start_dir = os.getcwd()

    # move to the root of the sage_main repository
    os.chdir(sage_main_repo)

    # there are several files that need to be treated
    # differently, since they correspond to .o or .so files
    # which are *not* created from .pyx files.
    non_cython_source_ls = []

    # we need to monkey with the path a bit to load the
    # appropriate lists of modules
    sys.path.append(os.getcwd())
    from module_list import ext_modules
    return_dir = os.getcwd()
    os.chdir('sage/ext')
    sys.path.append(os.getcwd())
    from gen_interpreters import modules as gen_modules
    os.chdir(return_dir)
    sys.path.pop()
    sys.path.pop()

    # finally build the list of special source files to
    # keep track of
    for x in ext_modules + gen_modules:
        for y in x.sources:
            if not y.endswith(dot + 'pyx'):
                non_cython_source_ls.append(y)


    # start walking the directory tree. we explicitly set
    # topdown=False, so that we can prune empty directories as we
    # encounter them.
    w = os.walk('build', topdown=False)
    for path, dirs, files in w:
        for filename in files:
            found = check_file(path + os.path.sep, filename,
                               special_source_files=non_cython_source_ls,
                               verbose=verbose)
            if not found:
                if dry_run:
                    print "Remove file: %s"%(path + slash + filename)
                else:
                    os.remove(path + slash + filename)
        for directory in dirs:
            file_ls = os.listdir(path + slash + directory)
            if prune_directories and (len(file_ls) == 0):
                if dry_run:
                    print "Remove empty directory: %s"%(path + slash + directory)
                else:
                    os.rmdir(path + slash + directory)
            elif len(file_ls) == 0:
                print "Remove empty directory: %s"%(path + slash + directory)

    # go back to where we started
    os.chdir(start_dir)

dry_run = False
prune_directories = True
verbose = False

# process command-line arguments
args, extra_args = getopt.getopt(sys.argv[1:], 'dpv', ['dry-run', 'prune-directories', 'verbose'])
if extra_args:
    # We could choose to throw an error if extra arguments are passed,
    # but this seems unnecessary.
    # raise ValueError, "unknown additional options: %s"%extra_args
    pass

for arg, opt in args:
    if opt:
        # Again, we could raise an error here, because this incorrect usage,
        # but it seems unnecessary.
        # raise ValueError, "unknown additional option: %s"%opt
        pass
    if arg in ['-d', '--dry-run']:
        dry_run = True
    elif arg in ['-p', '--prune-directories']:
        prune_directories = False
    elif arg in ['-v', '--verbose']:
        verbose = True

# call the function that does the work
clean_tree(prune_directories=prune_directories, dry_run=dry_run, verbose=verbose)
