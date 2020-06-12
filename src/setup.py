#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import time
from distutils import log
from distutils.core import setup

#########################################################
### Set source directory
#########################################################

import sage.env
sage.env.SAGE_SRC = os.getcwd()
from sage.env import *

from sage_setup import excepthook
sys.excepthook = excepthook

# This import allows instancemethods to be pickable
import sage_setup.fpickle_setup

#########################################################
### List of Extensions
###
### The list of extensions resides in module_list.py in
### the same directory as this file
#########################################################

from module_list import ext_modules

#########################################################
### Configuration
#########################################################

if len(sys.argv) > 1 and sys.argv[1] == "sdist":
    sdist = True
else:
    sdist = False

#########################################################
### Testing related stuff
#########################################################

# Remove (potentially invalid) star import caches
import sage.misc.lazy_import_cache
if os.path.exists(sage.misc.lazy_import_cache.get_cache_file()):
    os.unlink(sage.misc.lazy_import_cache.get_cache_file())


from sage_setup.command.sage_build import sage_build
from sage_setup.command.sage_build_cython import sage_build_cython
from sage_setup.command.sage_build_ext import sage_build_ext


#########################################################
### Discovering Sources
#########################################################

# TODO: This should be quiet by default
print("Discovering Python/Cython source code....")
t = time.time()
from sage_setup.find import find_python_sources
python_packages, python_modules = find_python_sources(
    SAGE_SRC, ['sage', 'sage_setup'])

log.debug('python_packages = {0}'.format(python_packages))

print("Discovered Python/Cython sources, time: %.2f seconds." % (time.time() - t))


from sage_setup.command.sage_install import sage_install

#########################################################
### Distutils
#########################################################

code = setup(name = 'sage',
      version     =  SAGE_VERSION,
      description = 'Sage: Open Source Mathematics Software',
      license     = 'GNU Public License (GPL)',
      author      = 'William Stein et al.',
      author_email= 'https://groups.google.com/group/sage-support',
      url         = 'https://www.sagemath.org',
      packages    = python_packages,
      package_data = {
          'sage.libs.gap': ['sage.gaprc'],
          'sage.doctest':  ['tests/*'],
          'sage': ['ext_data/*',
                   'ext_data/kenzo/*',
                   'ext_data/singular/*',
                   'ext_data/singular/function_field/*',
                   'ext_data/images/*',
                   'ext_data/doctest/*',
                   'ext_data/doctest/invalid/*',
                   'ext_data/doctest/rich_output/*',
                   'ext_data/doctest/rich_output/example_wavefront/*',
                   'ext_data/gap/*',
                   'ext_data/gap/joyner/*',
                   'ext_data/mwrank/*',
                   'ext_data/notebook-ipython/*',
                   'ext_data/nbconvert/*',
                   'ext_data/graphs/*',
                   'ext_data/pari/*',
                   'ext_data/pari/dokchitser/*',
                   'ext_data/pari/buzzard/*',
                   'ext_data/pari/simon/*',
                   'ext_data/magma/*',
                   'ext_data/magma/latex/*',
                   'ext_data/magma/sage/*',
                   'ext_data/valgrind/*',
                   'ext_data/threejs/*']
      },
      scripts = [## The sage script
                 'bin/sage',
                 ## Other scripts that should be in the path also for OS packaging of sage:
                 'bin/sage-eval',
                 ## Helper scripts invoked by sage script
                 ## (they would actually belong to something like libexec)
                 'bin/sage-cachegrind',
                 'bin/sage-callgrind',
                 'bin/sage-massif',
                 'bin/sage-omega',
                 'bin/sage-valgrind',
                 'bin/sage-version.sh'
                 ## Only makes sense in sage-the-distribution
                 'bin/sage-list-experimental',
                 'bin/sage-list-optional',
                 'bin/sage-list-packages',
                 'bin/sage-list-standard',
                 ## Uncategorized scripts in alphabetical order
                 'bin/math-readline',
                 'bin/sage-cleaner',
                 'bin/sage-clone-source',
                 'bin/sage-coverage',
                 'bin/sage-coverageall',
                 'bin/sage-cython',
                 'bin/sage-download-upstream',
                 'bin/sage-env',
                 'bin/sage-env-config',
                 # sage-env-config.in -- not to be installed',
                 'bin/sage-fix-pkg-checksums',
                 'bin/sage-fixdoctests',
                 'bin/sage-gdb-commands',
                 'bin/sage-grep',
                 'bin/sage-grepdoc',
                 'bin/sage-inline-fortran',
                 'bin/sage-ipynb2rst',
                 'bin/sage-ipython',
                 'bin/sage-location',
                 'bin/sage-maxima.lisp',
                 'bin/sage-native-execute',
                 'bin/sage-notebook',
                 'bin/sage-num-threads.py',
                 'bin/sage-open',
                 'bin/sage-preparse',
                 'bin/sage-pypkg-location',
                 'bin/sage-python',
                 'bin/sage-rebase.bat',
                 'bin/sage-rebase.sh',
                 'bin/sage-rebaseall.bat',
                 'bin/sage-rebaseall.sh',
                 'bin/sage-rst2sws',
                 'bin/sage-rst2txt',
                 'bin/sage-rsyncdist',
                 'bin/sage-run',
                 'bin/sage-run-cython',
                 'bin/sage-runtests',
                 'bin/sage-sdist',
                 'bin/sage-startuptime.py',
                 'bin/sage-sws2rst',
                 'bin/sage-unzip',
                 'bin/sage-update-src',
                 'bin/sage-update-version',
                 'bin/sage-upgrade',
                 ],
      cmdclass = dict(build=sage_build,
                      build_cython=sage_build_cython,
                      build_ext=sage_build_ext,
                      install=sage_install),
      ext_modules = ext_modules)
