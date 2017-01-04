"""Cross-platform compatibility routines and wrappers."""

#*****************************************************************************
#       Copyright (C) 2017 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import subprocess
import sys

from sage.env import SAGE_LOCAL
from sage.misc.decorators import sage_wraps


#################################################################
# Replacements (as needed) for Python stdlib functions to provide
# better platform compatibility
#################################################################
from ctypes.util import find_library
if sys.platform == 'cygwin':
    # find_library that works in cygwin adapted from
    # http://cygwin-ports.svn.sourceforge.net/viewvc/cygwin-ports/ports/trunk/lang/python/2.5.2-ctypes-util-find_library.patch?revision=8245&view=markup
    @sage_wraps(find_library)
    def find_library(name):
        for libdir in [os.path.join(SAGE_LOCAL, 'lib'),
                       '/usr/local/lib', '/usr/lib']:
            for libext in ['dll.a', 'a']:
                implib = os.path.join(libdir,
                                      'lib{0}.{1}'.format(name, libext))
                if not os.path.exists(implib):
                    continue

                cmd = ['dlltool', '-I', implib]

                p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE,
                                          universal_newlines=True)

                stdout, stderr = p.communicate()

                if p.returncode == 0:
                    return stdout.strip()
