"""
Sage Runtime Environment

AUTHORS:

- \R. Andrew Ohana (2012): Initial version.

"""

########################################################################
#       Copyright (C) 2013 R. Andrew Ohana <andrew.ohana@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, socket
import version

opj = os.path.join

# set default values for sage environment variables
# every variable can be overwritten by os.environ

# $VAR will be expanded as the appropriate variable,
# so do not make recursive definitions
SAGE_ENV = {
        # system info
        'UNAME'            : os.uname()[0],
        'HOSTNAME'         : socket.gethostname().replace('-','_').replace('/','_').replace('\\','_'),
        'LOCAL_IDENTIFIER' : '$HOSTNAME.%s'%os.getpid(),

        # bunch of sage directories and files
        'SAGE_ROOT'        : None,
        'SAGE_LOCAL'       : opj('$SAGE_ROOT', 'local'),
        'SAGE_SHARE'       : opj('$SAGE_LOCAL', 'share'),
        'SAGE_EXTCODE'     : opj('$SAGE_SHARE', 'sage', 'ext'),
        'SAGE_LOGS'        : opj('$SAGE_ROOT', 'logs', 'pkgs'),
        'SAGE_SPKG_INST'   : opj('$SAGE_LOCAL', 'var', 'lib', 'sage', 'installed'),
        'SAGE_DOC'         : opj('$SAGE_SRC', 'doc'),
        'DOT_SAGE'         : opj(os.environ.get('HOME','$SAGE_ROOT'), '.sage'),
        # SAGE_LIB is the site-packages directory if the sage library
        # has been installed, otherwise it is the same of SAGE_SRC
        'SAGE_SRC'         : opj('$SAGE_ROOT', 'src'),
        'SAGE_LIB'         : os.path.dirname(os.path.dirname(__file__)),

        # misc
        'SAGE_URL'         : 'http://sage.math.washington.edu/sage/',
        'SAGE_VERSION'     : version.version,
        'SAGE_DATE'        : version.date,
        }

# set any variables that are already in os.environ
for var in SAGE_ENV:
    try:
        exec(var + ' = os.environ[var]')
    except KeyError:
        pass

# create the dictionary of variables that need to
# be set to their default value
_tmp_env = {var:val for var,val in SAGE_ENV.items() if var not in globals()}

# end once everything has been set
while _tmp_env:
    for var,val in _tmp_env.items():
        try:
            # expand $VAR expressions
            if isinstance(val, str):
                for var2 in SAGE_ENV:
                    if '$'+var2 in val:
                        val = val.replace('$'+var2, eval(var2))
            exec(var + ' = val')
        except NameError:
            # the value depends upon something that hasn't been set yet
            # so try again later
            continue
        except TypeError:
            # depended on a variable that was set to None
            # so set this variable to None as well
            exec(var + ' = None')

    # remove things that have been set
    _tmp_env = {var:val for var,val in _tmp_env.items() if var not in globals()}

# repopulate SAGE_ENV with the corrected values
for var in SAGE_ENV:
    SAGE_ENV[var] = eval(var)

# post process
if ' ' in DOT_SAGE:
    if UNAME[:6] == 'CYGWIN':
        # on windows/cygwin it is typical for the home directory
        # to have a space in it.  Fortunately, users also have
        # write privileges to c:\cygwin\home, so we just put
        # .sage there.
        DOT_SAGE = "/home/.sage"
    else:
        print("Your home directory has a space in it.  This")
        print("will probably break some functionality of Sage.  E.g.,")
        print("the GAP interface will not work. A workaround")
        print("is to set the environment variable HOME to a")
        print("directory with no spaces that you have write")
        print("permissions to before you start sage.")

# delete temporary variables used for setting up sage.env
del _tmp_env, var, val, var2, opj, os, socket, version
