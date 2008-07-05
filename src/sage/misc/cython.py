"""
Cython -- C-Extensions for Python

AUTHORS:
    -- William Stein (2006-01-18): initial version
    -- William Stein (2007-07-28): update from sagex to cython
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os, sys

from misc import SPYX_TMP, SAGE_ROOT

def cblas():
    if os.environ.has_key('SAGE_CBLAS'):
        return os.environ['SAGE_CBLAS']
    elif os.path.exists('/usr/lib/libcblas.dylib') or \
         os.path.exists('/usr/lib/libcblas.so'):
        return 'cblas'
    elif os.path.exists('/usr/lib/libblas.dll.a'):   # untested.
        return 'blas'
    else:
        # This is very slow  (?), but *guaranteed* to be available.
        return 'gslcblas'

# In case of ATLAS we need to link against cblas as well as atlas
# In the other cases we just return the same library name as cblas()
# which is fine for the linker
def atlas():
    if os.environ.has_key('SAGE_CBLAS'):
        return os.environ['SAGE_CBLAS']
    elif os.path.exists('/usr/lib/libatlas.dylib') or \
        os.path.exists('/usr/lib/libatlas.so'):
        return 'atlas'
    elif os.path.exists('/usr/lib/libblas.dll.a'):   # untested
        return 'blas'
    else:
        # This is very slow  (?), but *guaranteed* to be available.
        return 'gslcblas'

include_dirs = ['%s/local/include/csage/'%SAGE_ROOT,
                '%s/local/include/'%SAGE_ROOT,  \
                '%s/local/include/python%s/'%(SAGE_ROOT, sys.version[:3]), \
                '%s/devel/sage/sage/ext/'%SAGE_ROOT, \
                '%s/devel/sage/'%SAGE_ROOT, \
                '%s/devel/sage/sage/gsl/'%SAGE_ROOT]


standard_libs = ['mpfr', 'gmp', 'gmpxx', 'stdc++', 'pari', 'm', 'curvesntl', \
                 'g0nntl', 'jcntl', 'rankntl', 'gsl', cblas(), atlas(), 'ntl', 'csage']

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
""" + s
    if lang != "c++": # has issues with init_csage()
        s = """
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

def cython(filename, verbose=False, compile_message=False,
          use_cache=False, create_local_c_file=False, annotate=True, sage_namespace=True):
    if not filename.endswith('pyx'):
        print "File (=%s) should have extension .pyx"%filename

    # base is the name of the .so module that we create.
    # The main constraint is that it is unique and determined by the file
    # that we're running Cython on, so that in some cases we can cache
    # the result (e.g., recompiling the same pyx file during the same session).
    base = sanitize(os.path.abspath(filename))

    # This is the *temporary* directory where we build the pyx file.
    # This is deleted when sage exits, which means pyx files must be
    # rebuilt every time Sage is restarted at present.
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
        try:
            if not os.path.isdir(G):
                os.unlink(G)
        except OSError:
            pass

    # Get the absolute path to the directory that contains the pyx file.
    # We will use this only to make some convenient symbolic links.
    abs_base = os.path.split(os.path.abspath(filename))[0]

    cmd = 'cd "%s"; ln -sf "%s"/* .'%(build_dir, abs_base)
    os.system(cmd)

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

    additional_source_files = ",".join(["'"+os.path.abspath(os.curdir)+"/"+fname+"'" \
                                        for fname in additional_source_files])

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
extra_compile_args = ['-w','-O2']

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

    cython_include = ' '.join(["-I '%s'"%x for x in includes if len(x.strip()) > 0 ])

    options = ['-p', '--incref-local-binop']
    if annotate:
        options.append('-a')
    if sage_namespace:
        options.append('--pre-import sage.all')

    cmd = "cd '%s' && cython %s %s '%s.pyx' 1>log 2>err " % (build_dir, ' '.join(options), cython_include, name)

    if create_local_c_file:
        target_c = '%s/_%s.c'%(os.path.abspath(os.curdir), base)
        if language == 'c++':
            target_c = target_c + "pp"
        cmd += " && cp '%s.c' '%s'"%(name, target_c)
        if annotate:
            target_html = '%s/_%s.html'%(os.path.abspath(os.curdir), base)
            cmd += " && cp '%s.html' '%s'"%(name, target_html)

    if verbose:
        print cmd
    if os.system(cmd):
        log = open('%s/log'%build_dir).read()
        err = subtract_from_line_numbers(open('%s/err'%build_dir).read(), offset)
        raise RuntimeError, "Error converting %s to C:\n%s\n%s"%(filename, log, err)

    if language=='c++':
        os.system("cd '%s' && mv '%s.c' '%s.cpp'"%(build_dir,name,name))

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


################################################################
# COMPILE
################################################################
def cython_lambda(vars, expr,
                 verbose=False,
                 compile_message=False,
                 use_cache=False):
    """
    Create a compiled function which evaluates expr assuming machine
    values for vars.

    WARNING: This implementation is not well tested.

    INPUT:
        vars -- list of pairs (variable name, c-data type), where
                the variable names and data types are strings.
            OR -- a string such as
                         'double x, int y, int z'
        expr -- an expression involving the vars and constants;
                You can access objects defined in the current
                module scope globals() using sagobject_name.
                See the examples below.

    EXAMPLES:
    We create a Lambda function in pure Python (using the r to make sure the 3.2
    is viewed as a Python float):
        sage: f = lambda x,y: x*x + y*y + x + y + 17r*x + 3.2r

    We make the same Lambda function, but in a compiled form.
        sage: g = cython_lambda('double x, double y', 'x*x + y*y + x + y + 17*x + 3.2')
        sage: g(2,3)
        55.200000000000003
        sage: g(0,0)
        3.2000000000000002

    We access a global function and variable.
        sage: a = 25
        sage: f = cython_lambda('double x', 'sage.math.sin(x) + sage.a')
        sage: f(10)
        24.455978889110629
        sage: a = 50
        sage: f(10)
        49.455978889110632
    """
    if isinstance(vars, str):
        v = vars
    else:
        v = ', '.join(['%s %s'%(typ,var) for typ, var in vars])

    s = """
class _s:
   def __getattr__(self, x):
       return globals()[x]

sage = _s()

def f(%s):
 return %s
    """%(v, expr)
    if verbose:
        print s
    import sage.misc.misc
    tmpfile = sage.misc.misc.tmp_filename() + ".spyx"
    open(tmpfile,'w').write(s)

    import sage.server.support
    d = {}
    sage.server.support.cython_import_all(tmpfile, d,
                                         verbose=verbose, compile_message=compile_message,
                                         use_cache=use_cache,
                                         create_local_c_file=False)
    return d['f']



def sanitize(f):
    """
    Given a filename f, replace it by a filename that is a valid Python
    module name.

    This means that the characters are all alphanumeric or _'s and
    doesn't begin with a numeral.
    """
    s = ''
    if f[0].isdigit():
        s += '_'
    for a in f:
        if a.isalnum():
            s += a
        else:
            s += '_'
    return s



