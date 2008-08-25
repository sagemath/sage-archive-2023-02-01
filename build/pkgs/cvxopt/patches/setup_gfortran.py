from distutils.core import setup, Extension
from os import listdir
import os

SAGE_LIB=os.environ['SAGE_LOCAL']+'/lib'

ATLAS_LIB_DIR=SAGE_LIB

if os.uname()[0]=="Darwin":
    libraries = ['m','lapack','blas','gfortran']
else:
    libraries = ['m','lapack','cblas','f77blas','atlas','gfortran']


# Set to 1 if you are installing the fftw module.
BUILD_FFTW = 0

# directory containing libfftw3 (used only when BUILD_FFTW = 1)
FFTW_LIB_DIR = '/usr/lib'

# directory containing fftw.h (used only when BUILD_FFTW = 1)
FFTW_INC_DIR = '/usr/include'

# Set to 1 if you are installing the glpk module.
BUILD_GLPK = 0

# directory containing libglpk (used only when BUILD_GLPK = 1)
GLPK_LIB_DIR = '/usr/lib'

# directory containing glpk.h (used only when BUILD_GLPK = 1)
GLPK_INC_DIR = '/usr/include'

# Set to 1 if you are installing the MOSEK 4.0 module.
BUILD_MOSEK = 0

# directory containing libmosek (used only when BUILD_MOSEK = 1)
MOSEK_LIB_DIR = '/usr/local/mosek/4/tools/platform/linux32x86/bin'

# directory containing mosek.h (used only when BUILD_MOSEK = 1)
MOSEK_INC_DIR = '/usr/local/mosek/4/tools/platform/linux32x86/h'

# Set to 1 if you are installing the DSDP module.
BUILD_DSDP = 0

# directory containing libdsdp.a (used only when BUILD_DSDP = 1)
DSDP_LIB_DIR = '/usr/lib'

# directory containing dsdp5.h (used only when BUILD_DSDP = 1)
DSDP_INC_DIR = '/usr/include'

extmods = []

# optional modules

if BUILD_FFTW:
    fftw = Extension('fftw', libraries = ['fftw3', 'blas'],
        include_dirs = [ FFTW_INC_DIR ],
        library_dirs = [ FFTW_LIB_DIR, ATLAS_LIB_DIR ],
        sources = ['C/fftw.c'] )
    extmods += [fftw];

if BUILD_GLPK:
    glpk = Extension('glpk', libraries = ['glpk'],
        include_dirs = [ GLPK_INC_DIR ],
        library_dirs = [ GLPK_LIB_DIR ],
        sources = ['C/glpk.c'] )
    extmods += [glpk];

# If the MOSEK module fails to build, modify the script according to
# http://www.mosek.com/products/4_0/tools/doc/html/tools/node20.html,
# section 18.13.
if BUILD_MOSEK:
    mosek = Extension('mosek',
        include_dirs = [ MOSEK_INC_DIR ],
        library_dirs = [ MOSEK_LIB_DIR ],
        libraries = ['mosek', 'pthread', 'c', 'dl', 'm'],
        sources = ['C/mosek.c'] )
    extmods += [mosek];

if BUILD_DSDP:
    dsdp = Extension('dsdp', libraries = ['dsdp', 'blas', 'lapack'],
        include_dirs = [ DSDP_INC_DIR ],
        library_dirs = [ DSDP_LIB_DIR, ATLAS_LIB_DIR ],
        sources = ['C/dsdp.c'] )
    extmods += [dsdp];


# required modules

# Modify this for compilation on Windows.
# Set to True if your BLAS/LAPACK do not use trailing underscores
# (eg, on Windows).
BLAS_NOUNDERSCORES = False
if BLAS_NOUNDERSCORES:
    MACROS = [('BLAS_NO_UNDERSCORE','')]
else:
    MACROS = []

base = Extension('base', libraries = libraries,
    library_dirs = [ ATLAS_LIB_DIR ],
    define_macros = MACROS,
    sources = ['C/base.c','C/dense.c','C/sparse.c'])

random = Extension('random',
    sources = ['C/random.c', 'C/rngs/rngs.c', 'C/rngs/rvgs.c'])

blas = Extension('blas', libraries = libraries,
    library_dirs = [ ATLAS_LIB_DIR ],
    define_macros = MACROS,
    sources = ['C/blas.c'] )

lapack = Extension('lapack', libraries = libraries,
    library_dirs = [ ATLAS_LIB_DIR ],
    define_macros = MACROS,
    sources = ['C/lapack.c'] )

umfpack = Extension('umfpack',
    include_dirs = [ 'C/SuiteSparse/UMFPACK/Include',
        'C/SuiteSparse/AMD/Include', 'C/SuiteSparse/AMD/Source',
        'C/SuiteSparse/UFconfig' ],
    library_dirs = [ ATLAS_LIB_DIR ],
    define_macros = MACROS,
    libraries = libraries,
    sources = [ 'C/umfpack.c',
        'C/SuiteSparse/UMFPACK/Source/umfpack_global.c',
        'C/SuiteSparse/UMFPACK/Source/umfpack_tictoc.c' ] +
        ['C/SuiteSparse_cvxopt_extra/umfpack/' + s for s in
            listdir('C/SuiteSparse_cvxopt_extra/umfpack')])

# Build for int or long?
import sys
if sys.maxint > 2**31: MACROS += [('DLONG','')]

cholmod = Extension('cholmod',
    library_dirs = [ ATLAS_LIB_DIR ],
    libraries = libraries,
    include_dirs = [ 'C/SuiteSparse/CHOLMOD/Include',
        'C/SuiteSparse/COLAMD', 'C/SuiteSparse/AMD/Include',
        'C/SuiteSparse/UFconfig', 'C/SuiteSparse/COLAMD/Include' ],
    define_macros = MACROS + [('NPARTITION', '1')],
    sources = [ 'C/cholmod.c' ] +
        ['C/SuiteSparse/AMD/Source/' + s for s in ['amd_global.c',
            'amd_postorder.c', 'amd_post_tree.c', 'amd_2.c']] +
        ['C/SuiteSparse/COLAMD/Source/' + s for s in ['colamd.c',
            'colamd_global.c']] +
        ['C/SuiteSparse/CHOLMOD/Core/' + s for s in
            listdir('C/SuiteSparse/CHOLMOD/Core') if s[-2:] == '.c' and
            s[0] == 'c'] +
        ['C/SuiteSparse/CHOLMOD/Cholesky/' + s for s in
            listdir('C/SuiteSparse/CHOLMOD/Cholesky') if s[-2:] == '.c'
            and s[0] == 'c'] +
        ['C/SuiteSparse/CHOLMOD/Check/cholmod_check.c'] +
        ['C/SuiteSparse/CHOLMOD/Supernodal/' + s for s in
            listdir('C/SuiteSparse/CHOLMOD/Supernodal') if
            s[-2:] == '.c' and s[0] == 'c'] )

amd = Extension('amd',
    include_dirs = [ 'C/SuiteSparse/AMD/Include',
        'C/SuiteSparse/UFconfig' ],
    define_macros = MACROS,
    sources = [ 'C/amd.c' ] + [ 'C/SuiteSparse/AMD/Source/' + s for s in
        listdir('C/SuiteSparse/AMD/Source') if s[-2:] == '.c' ])

extmods += [base, blas, lapack, random, umfpack, cholmod, amd]

setup (name = 'cvxopt',
    description = 'Convex optimization package',
    version = '0.9',
    long_description = '''
CVXOPT is a free software package for convex optimization based on the
Python programming language. It can be used with the interactive Python
interpreter, on the command line by executing Python scripts, or
integrated in other software via Python extension modules. Its main
purpose is to make the development of software for convex optimization
applications straightforward by building on Python's extensive standard
library and on the strengths of Python as a high-level programming
language.''',
    author='J. Dahl and L. Vandenberghe',
    author_email='joachim@es.aau.dk, vandenbe@ee.ucla.edu',
    url='http://abel.ee.ucla.edu/cvxopt',
    license='GNU GPL version 3',
    ext_package = "cvxopt",
    ext_modules = extmods,
    package_dir = {"cvxopt": "python"},
    packages = ["cvxopt"])
