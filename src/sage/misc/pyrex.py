"""
Preparse and compile Pyrex file

AUTHORS:
    -- William Stein, 2006-01-18

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

from misc import SPYX_TMP, SAGE_ROOT

include_dirs = ['%s/local/include'%SAGE_ROOT,  \
                '%s/local/include/python'%SAGE_ROOT, \
                '%s/devel/sage/sage/ext'%SAGE_ROOT, \
                '%s/devel/sage/'%SAGE_ROOT, \
                '%s/devel/sage/sage/gsl'%SAGE_ROOT]

standard_libs = ['mpfr', 'gmp', 'gmpxx', 'stdc++', 'pari', 'm', 'mwrank', 'gsl', 'gslcblas', 'ntl', 'csage']

offset = 0

def parse_keywords(kwd, s):
    j = 0
    v = []
    while True:
        i = s[j:].find(kwd)
        if i == -1: break
        j = i + j
        s = s[:j] + '#' + s[j:]
        j += len(kwd) + 1
        k = s[j:].find('\n')
        if k == -1:
            k = len(s)
        for X in s[j:j+k].split():
            if X[0] == '#':   # skip rest of line
                break
            v.append(X)
    return v, s

def environ_parse(s):
    i = s.find('$')
    if i == -1:
        return s
    j = s[i:].find('/')
    if j == -1:
        j = len(s)
    else:
        j = i + j
    name = s[i+1:j]
    if os.environ.has_key(name):
        s = s[:i] + os.environ[name] + s[j:]
    else:
        return s
    return environ_parse(s)

def pyx_preparse(s):
    lang = parse_keywords('clang', s)
    if lang[0]:
        lang = lang[0][0]
    else:
        lang = "c"

    v, s = parse_keywords('clib', s)
    libs = v + standard_libs

    additional_source_files, s = parse_keywords('cfile', s)

    v, s = parse_keywords('cinclude', s)
    inc = [environ_parse(x.replace('"','').replace("'","")) for x in v] + include_dirs
    s = """
include "cdefs.pxi"
include "interrupt.pxi"  # ctrl-c interrupt block support
include "stdsage.pxi"  # ctrl-c interrupt block support
""" + s
    return s, libs, inc, lang, additional_source_files

################################################################
# If the user attaches a .spyx file and changes it, we have
# to reload an .so.
#
# PROBLEM: Python does not allow one to reload an .so extension module.
# Solution, we create a different .so file and load that one,
# overwriting the definitions of everything in the original .so file.
#
# HOW: By using a sequence_number for each .spyx file; we keep
# these sequence numbers in a dict.
#
################################################################

sequence_number = {}

def pyrex(filename, verbose=False, compile_message=False,
          use_cache=False):
    if filename[-5:] != '.spyx':
        print "File (=%s) must have extension .spyx"%filename

    base = os.path.split(os.path.splitext(filename)[0])[1]

    build_dir = '%s/%s'%(SPYX_TMP, base)
    if os.path.exists(build_dir):
        # There is already a module here.  Maybe we do not have to rebuild?
        # Find the name.
        if use_cache:
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

    os.system('cd "%s"; ln -s "%s"/devel/sage/sage/ext/*.c .'%(build_dir, SAGE_ROOT))

    if compile_message:
        print "Compiling %s..."%filename

    F = open(filename).read()

    F, libs, includes, language, additional_source_files = pyx_preparse(F)

    # add the working directory to the includes so custom headers etc. work
    includes.append(os.path.split(os.path.splitext(filename)[0])[0])

    if language == 'c++':
        extension = "cpp"
    else:
        extension = "c"

    global sequence_number
    if not sequence_number.has_key(base):
        sequence_number[base] = 0
    name = '%s_%s'%(base, sequence_number[base])

    # increment the sequence number so will use a different one next time.
    sequence_number[base] += 1

    additional_source_files = ",".join(["'"+os.path.abspath(os.curdir)+"/"+filename+"'" \
                                        for filename in additional_source_files])

    pyx = '%s/%s.pyx'%(build_dir, name)
    open(pyx,'w').write(F)
    setup="""
# Build using 'python setup.py'
import distutils.sysconfig, os, sys
from distutils.core import setup, Extension

if not os.environ.has_key('SAGE_ROOT'):
    print "    ERROR: The environment variable SAGE_ROOT must be defined."
    sys.exit(1)
else:
    SAGE_ROOT  = os.environ['SAGE_ROOT']
    SAGE_LOCAL = SAGE_ROOT + '/local/'

extra_link_args =  ['-L' + SAGE_LOCAL + '/lib']
extra_compile_args = ['-w']

ext_modules = [Extension('%s', sources=['%s.%s', %s],
                     libraries=%s,
                     library_dirs=[SAGE_LOCAL + '/lib/'],
                     extra_compile_args = extra_compile_args,
                     extra_link_args = extra_link_args,
                     language = '%s' )]

setup(ext_modules = ext_modules,
      include_dirs = %s)
    """%(name, name, extension, additional_source_files, libs, language, includes)
    open('%s/setup.py'%build_dir,'w').write(setup)

    pyrex_include = ' '.join(['-I %s'%x for x in includes])

    target_c = '%s/_%s.c'%(os.path.abspath(os.curdir), base)

    if language == 'c++':
        target_c = target_c + "pp"



    cmd = 'cd %s && pyrexc -p %s %s.pyx 1>log 2>err && cp %s.c %s'%(build_dir, pyrex_include, name,
                                                                  name, target_c)

    if verbose:
        print cmd
    if os.system(cmd):
        log = open('%s/log'%build_dir).read()
        err = subtract_from_line_numbers(open('%s/err'%build_dir).read(), offset)
        raise RuntimeError, "Error converting %s to C:\n%s\n%s"%(filename, log, err)

    if language=='c++':
        os.system("cd %s && mv %s.c %s.cpp"%(build_dir,name,name))

##     if make_c_file_nice and os.path.exists(target_c):
##         R = open(target_c).read()
##         R = "/* THIS IS A PARSED TO MAKE READABLE VERSION OF THE C FILE. */" + R

##         # 1. Get rid of the annoying __pyx_'s before variable names.
##         # R = R.replace('__pyx_v_','').replace('__pyx','')
##         # 2. Replace the line number references by the actual code from the file,
##         #    since it is very painful to go back and forth, and the philosophy
##         #    of SAGE is that everything that can be very easy *is*.

##         pyx_file = os.path.abspath('%s/%s.pyx'%(build_dir,name))
##         S = '/* "%s":'%pyx_file
##         n = len(S)
##         last_i = -1
##         X = F.split('\n')
##         stars = '*'*80
##         while True:
##             i = R.find(S)
##             if i == -1 or i == last_i: break
##             last_i = i
##             j = R[i:].find('\n')
##             if j == -1: break
##             line_number = int(R[i+n: i+j])
##             try:
##                 line = X[line_number-1]
##             except IndexError:
##                 line = '(missing code)'
##             R = R[:i+2] + '%s\n\n Line %s: %s\n\n%s'%(stars, line_number, line, stars) + R[i+j:]

##         open(target_c,'w').write(R)


    cmd = 'cd %s && python setup.py build 1>log 2>err'%build_dir
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


