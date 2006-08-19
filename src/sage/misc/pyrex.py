"""
Preparse and compile Pyrex file

AUTHORS:
    -- William Stein, 2006-01-18

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os

from misc import SPYX_TMP, SAGE_ROOT

include_dirs = ['%s/local/include'%SAGE_ROOT,  \
                '%s/local/include/python2.4'%SAGE_ROOT, \
                '%s/devel/sage/ext'%SAGE_ROOT]

standard_libs = ['mpfr', 'gmp', 'gmpxx', 'stdc++', 'pari', 'm', 'mwrank']

offset = 0

def pyx_preparse(s):
    t = """import sage.all; sage = sage.all
include 'interrupt.pxi'
include 'cdefs.pxi'
"""
    # TODO: the ZZ, QQ, etc. could be replaced by C-level constructors for efficiency?
    global offset
    offset = len(t.split('\n')) - 1
    s = t + s
    return s, standard_libs, include_dirs

################################################################
# If the user attaches a .spyx file and changes it, we have
# to reload an .so.
#
# PROBLEM: Python does not allow one to reload an .so extension module.
# Solution, we create a different .so file and load that one,
# overwriting the definitions of everythin in the original .so file.
#
# HOW: By using a sequence_number for each .spyx file; we keep
# these sequence numbers in a dict.
#
################################################################

sequence_number = {}

def pyrex(filename, verbose=False, compile_message=False):
    if filename[-5:] != '.spyx':
        print "File (=%s) must have extension .spyx"%filename

    base = os.path.split(os.path.splitext(filename)[0])[1]

    build_dir = '%s/%s'%(SPYX_TMP, base)
    if os.path.exists(build_dir):
        # There is already a module here.  Maybe we do not have to rebuild?
        # Find the name.
        prev_so = [F for F in os.listdir(build_dir) if F[-3:] == '.so']
        if len(prev_so) > 0:
            prev_so = prev_so[0]     # should have length 1 because of deletes below
            if os.path.getmtime(filename) <= os.path.getmtime('%s/%s'%(build_dir, prev_so)):
                # We do not have to rebuild.
                return prev_so[:-3], build_dir
    else:
        os.makedirs(build_dir)
    for F in os.listdir(build_dir):
        G = '%s/%s'%(build_dir,F)
        if not os.path.isdir(G):
            os.unlink(G)

    os.system('cd "%s"; ln -s "%s/devel/sage/sage/ext/"*.pxi .'%(build_dir,
                   SAGE_ROOT))
    os.system('cd "%s"; ln -s "%s/devel/sage/sage/ext/interrupt.c" .'%(build_dir,
                   SAGE_ROOT))

    if compile_message:
        print "Compiling %s..."%filename

    F = open(filename).read()

    F, libs, includes = pyx_preparse(F)

    global sequence_number
    if not sequence_number.has_key(base):
        sequence_number[base] = 0
    name = '%s_%s'%(base, sequence_number[base])

    # increment the sequence number so will use a different one next time.
    sequence_number[base] += 1

    pyx = '%s/%s.pyx'%(build_dir, name)
    open(pyx,'w').write(F)
    setup="""
import distutils.sysconfig, os, sys
from distutils.core import setup, Extension

ext_modules = [Extension('%s', sources=['%s.c', 'interrupt.c'], libraries=%s,
                     extra_compile_args = ['-w'])]

setup(ext_modules = ext_modules,
      include_dirs = %s)
    """%(name, name, libs, includes)
    open('%s/setup.py'%build_dir,'w').write(setup)
    cmd = 'cd %s; pyrexc -I "%s/devel/sage" -I "%s/devel/sage/ext" %s.pyx 1>log 2>err'%(build_dir, SAGE_ROOT, SAGE_ROOT, name)
    if verbose:
        print cmd
    if os.system(cmd):
        log = open('%s/log'%build_dir).read()
        err = subtract_from_line_numbers(open('%s/err'%build_dir).read(), offset)
        raise RuntimeError, "Error converting %s to C:\n%s\n%s"%(filename, log, err)

    cmd = 'cd %s; python setup.py build 1>log 2>err'%build_dir
    if verbose: print cmd
    if os.system(cmd):
        log = open('%s/log'%build_dir).read()
        err = open('%s/err'%build_dir).read()
        raise RuntimeError, "Error compiling %s:\n%s\n%s"%(filename, log, err)

    # Move from lib directory.
    cmd = 'mv %s/build/lib.*/* %s'%(build_dir, build_dir)
    if verbose: print cmd
    if os.system(cmd):
        raise RuntimeError, "Error copying extension module for %s"%filename


    return name, build_dir



def subtract_from_line_numbers(s, n):
    ans = []
    for X in s.split('\n'):
        i = X.find(':')
        j = i+1 + X[i+1:].find(':')
        try:
            ans.append('%s:%s:%s\n'%(X[:i], int(X[i+1:j]) - n, X[j+1:]))
        except ValueError:
            ans.append(X)
    return '\n'.join(ans)
