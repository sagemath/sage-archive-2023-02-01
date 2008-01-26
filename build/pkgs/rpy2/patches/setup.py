"""
This file builds rpy.

It has been adapted to allow multiple versions of the rpy shared
library to be included in the build.  To enable this, set the
environment variable "RHOMES" to a colon-delimited
(semicolon-delimited on MS-Windows) list of R installation directories.

For linux csh derivatives (including tcsh) use something like

  setenv RHOMES  /usr/local/R-2.1.0/lib/R:/usr/local/R-2.1.1/lib/R
  python setup.py bdist

For linux sh derivatives (including bash) use something like:

  RHOMES=/usr/local/R-2.1.0/lib/R:/usr/local/R-2.1.1/lib/R
  export RHOMES
  python setup.py bdist

For windows, edit setup.32 to include the correct paths, then use
something like:

  copy setup.Win32 setup.cfg
  python setup.py build --compiler=mingw32 bdist_wininst

setup.py will automatically look in C:\Program Files\R to determine
which versions of R are installed, and will build an installer that
can be used for each of these R versions.

See the files INSTALL.UNIX and INSTALL.WINDOWS for more details.
"""

DEBUG=True

import os, os.path, sys, shutil
from distutils.core import setup, Extension
from distutils.sysconfig import *
from distutils.errors import *
import rpy_tools
if sys.platform=="win32":
    import rpy_wintools

# Get list of R Home directories that we will be processing
try:
    if sys.platform=="win32":
        RHOMES = os.environ['RHOMES'].split(';')
    else:
        RHOMES = os.environ['RHOMES'].split(':')
except:
    RHOMES = []

print "RHOMES=", RHOMES
print "DEBUG=", DEBUG

if not RHOMES:
    if sys.platform=="win32":
        RHOMES=rpy_wintools.get_RHOMES()
    else:
        RHOMES = [rpy_tools.get_R_HOME(force_exec=False)]
    print "Setting RHOMES to ", RHOMES

# ensure RHOMES is in ascii, since most command line tools are unhappy with unicode...
RHOMES = map( lambda x: x.encode('ascii'), RHOMES )

# On windows, check for/create the python link library
#if sys.platform=="win32":
#    rpy_wintools.CreatePythonWinLib()


#  For sage
import os
SAGE_LOCAL=os.environ['SAGE_LOCAL']

# Find the directory where the g95 libraries are; it looks like this:
#        GCC_LIB_DIR= SAGE_LOCAL+'/lib/gcc-lib/i686-pc-linux-gnu/4.0.3'

GCC_LIB_DIR=SAGE_LOCAL+"/lib/"
if os.path.exists(GCC_LIB_DIR + "gcc-lib"):
    GCC_LIB_DIR += "gcc-lib/"
    GCC_LIB_DIR += os.listdir(GCC_LIB_DIR)[0] + "/"
    GCC_LIB_DIR += os.listdir(GCC_LIB_DIR)[0] + "/"

###

modules = []
for RHOME in RHOMES:

    RHOME = RHOME.strip()

    if DEBUG:
        # to avoid strict prototypes errors from R includes
        get_config_vars()['OPT'] = '-g -Wall'
    else:
        # to avoid strict prototypes errors from R includes
        get_config_vars()['OPT'] = '-DNDEBUG -g -O3 -Wall'

    # get the Python version
    if sys.version[:3] >= '2.2':
        DEFINE = []
        UNDEF = ['PRE_2_2']
    else:
        DEFINE = [('PRE_2_2', None)]
        UNDEF = []

    # configure the R paths
    RVERSION = rpy_tools.get_R_VERSION(RHOME, force_exec=True)
    RVER     = rpy_tools.get_R_VERSION_CODE(RVERSION)

    print "### Using R verion %s installed at %s ###" % (RVERSION, RHOME)

    r_libs = [ # Different verisons of R put .so/.dll in different places
              os.path.join(RHOME, 'bin'),  # R 2.0.0+
              os.path.join(RHOME, 'lib'),  # Pre 2.0.0
             ]

    print "RHOME=",RHOME

    base_source_files = ["src/rpymodule.c", "src/R_eval.c",
                         "src/io.c"]

    exists = os.path.exists
    mtime = os.path.getmtime

    # Make one copy of the source files for this R version
    source_files = []
    for f in base_source_files:
        # file.c => file2010.c
        nfile = f[0:-2] + RVER + '.c'
        print "copying %s -> %s" % (f, nfile)
        if (not exists(nfile)) or  ( mtime(f) > mtime(nfile) ):
            shutil.copy(f, nfile)
        source_files.append(nfile)

    if sys.platform=='win32':
        include_dirs = [ os.path.join(RHOME.strip(), 'include'),
                         'src' ]
        libraries= ['R']
        library_dirs = r_libs
        runtime_libs = []
        extra_compile_args= []
        source_files = source_files + ["src\\setenv.c"]
    elif sys.platform=='darwin':
        include_dirs = [ os.path.join(RHOME.strip(), 'include'),
                         'src' ]
        libraries=['R']
        library_dirs= r_libs
        runtime_libs = r_libs
        extra_compile_args=[]
    elif sys.platform=='osf1V5':
        include_dirs = [ os.path.join(RHOME.strip(), 'include'),
                         'src' ]
        libraries=['R','Rlapack']
        library_dirs = r_libs
        runtime_libs = r_libs
        extra_compile_args=["-shared"]
        source_files = source_files + ["src/setenv.c"]
    else: # unix-like systems, this is known to work for Linux and Solaris
        include_dirs = [ os.path.join(RHOME.strip(), 'include'),
                         'src', '/usr/share/R/include' ]
        libraries=['R'] # changed for Sage: since we build R with ATLAS we don't need RLapack
        library_dirs = r_libs
        runtime_libs = r_libs
        extra_compile_args=["-shared"]
        source_files = source_files + ["src/setenv.c"]

    # Discover which array packages are present
    try:
        import numpy
        DEFINE.append(('WITH_NUMERIC', '3'))
        DEFINE.append(('PY_ARRAY_TYPES_PREFIX', 'PyArray_'))
        include_dirs.append(numpy.get_include())
    except ImportError:
        # fall back to Numeric
        try:
            import Numeric
            DEFINE.append(('WITH_NUMERIC', '1'))
        except ImportError:
            UNDEF.append('WITH_NUMERIC')

    # Added for Sage
    libraries.append('lapack')
    libraries.append('cblas')
    libraries.append('atlas')
    if os.popen2('which_fortran')[1].read().startswith('g95'):
        libraries.append('f95')
    else:
        libraries.append('gfortran')
    library_dirs.append(GCC_LIB_DIR)



    # get the RPy version
    from rpy_version import rpy_version

    LONG_DESC = """RPy provides a robust Python interface to the R
    programming language.  It can manage all kinds of R objects and can
    execute arbitrary R functions. All the errors from the R language are
    converted to Python exceptions."""

    # use R version specific shared library name
    shlib_name = "_rpy%s" % RVER
    DEFINE.append( ('RPY_SHNAME', shlib_name))
    DEFINE.append( ('INIT_RPY', 'init_rpy%s' % RVER ) )


    # add debugging

    modules.append( Extension(
        shlib_name,
        source_files,
        include_dirs=include_dirs,
        libraries=libraries,
        library_dirs=library_dirs,
        define_macros=DEFINE,
        undef_macros=UNDEF,
        extra_compile_args=extra_compile_args,
        runtime_library_dirs = runtime_libs,
        ) )


setup(name="rpy",
      version=rpy_version,
      description="Python interface to the R language",
      maintainer="Gregory R. Warnes",
      maintainer_email="warnes@bst.rochester.edu",
      url="http://rpy.sourceforge.net",
      license="GPL",
      long_description=LONG_DESC,
      py_modules=['rpy', 'rpy_io', 'rpy_version', 'rpy_tools', 'rpy_options',
                  'rpy_wintools'],
      ext_modules=modules
      )
