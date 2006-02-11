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
from preparser import preparse_file

include_dirs = ['%s/local/include'%SAGE_ROOT, '%s/local/include/python2.4'%SAGE_ROOT]
standard_libs = ['mpfr', 'gmp', 'gmpxx', 'stdc++', 'pari', 'm', 'mwrank']

offset = 0

def pyx_preparse(s):
    t = """
import sage.all
sage = sage.all
ZZ = sage.ZZ
QQ = sage.QQ
RR = sage.RR

include 'interrupt.pxi'
"""
    # TODO: the ZZ, QQ, etc. could be replaced by C-level constructors for efficiency?
    global offset
    offset = len(t.split('\n')) - 1
    s = t + s
    return s, standard_libs, include_dirs

sequence_number = 0

def pyrex(filename, verbose=False, compile_message=False):
    if filename[-5:] != '.spyx':
        print "File (=%s) must have extension .spyx"%filename

    global sequence_number
    name = '%s_%s'%(filename[:-5],sequence_number)
    sequence_number += 1

    build_dir = '%s/%s'%(SPYX_TMP, filename[:-5])
    if os.path.exists(build_dir):
        # There is already a module here.  Maybe we do not have to rebuild?
        # Find the name
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

    os.system('cp %s/devel/sage/sage/ext/interrupt* %s/'%(SAGE_ROOT,
                                                          build_dir))

    if compile_message:
        print "***************************************************"
        print "                Recompiling %s"%filename
        print "***************************************************"

    F = open(filename).read()

    if F.find('__no_preparse__') == -1:
        F = preparse_file(F)
    F, libs, includes = pyx_preparse(F)
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
    cmd = 'cd %s; pyrexc %s.pyx 1>log 2>err'%(build_dir, name)
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
