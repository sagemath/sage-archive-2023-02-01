######################################################################
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
######################################################################

import platform, os, sys, time, shutil, glob, subprocess


# this dictionary will hold all configuration information
conf = dict()

# Here is a list of keys and (some of) their possible values. First,
# strings:
#
# system:  Linux, SunOS, Darwin, FreeBSD, CYGWIN
# machine: i386, x86_64,   # Linux, Darwin, *BSD
#          sun4u, i86pc    # SunOS
# processor: i386, x86_64, powerpc, sparc
# bits:    32bit, 64bit
# fortran: g95, gfortran
# ld:      GNU, Solaris

# The following pre-defined boolean values are determined from the
# strings. The keys are distinguished by a question mark at the
# end. If possible use these keys as they guard against typos.
#
# Linux?, Solaris?, Darwin?, FreeBSD?, CYGWIN?  # system
# Intel?, PPC?, SPARC?                          # processor
# 64bit?, 32bit?                                # bit width
# fortran_g95?, fortran_GNU?                    # sage_fortran
# linker_GNU?, linker_Solaris?                  # ld


######################################################################
### Sanity check
######################################################################

if not os.environ.has_key('SAGE_LOCAL'):
    print "SAGE_LOCAL undefined ... exiting"
    print "Maybe run 'sage -sh'?"
    sys.exit(1)


######################################################################
### Functions
######################################################################

def try_run(command, ignore=False):
    """
    Try to execute ``command`` and return its output as string.

    Return ``None`` if command not found. Return ``None`` if execution
    fails and ``ignore=False`` (default). Otherwise, return output of
    ``command`` as string. Here, output always means the concatenation
    of stdout and stderr.
    """
    f = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = f.communicate()
    rc = f.wait()
    if (not ignore) and (rc!=0):
        return None
    # concatenate stdout and stderr
    return result[0].strip() + result[1].strip()


def cp(source_pattern, destination):
    """
    Portable implementation of "cp -p"
    """
    for filename in glob.iglob(source_pattern):
        print 'Copying', filename, 'to', destination
        shutil.copy2(filename, destination)


class edit_in_place(object):
    """
    Edit a file in-place by search/replacing regular expressions.
    """
    def __init__(self, filename):
        self.filename = filename
        f = open(filename, 'r')
        self.data = f.read()
        f.close()

    def replace(self, find_regex, subs, count=0):
        self.data = glob.re.sub(find_regex, subs, self.data, count)
        return self

    def close(self):
        f = open(self.filename, 'w')
        f.write(self.data)
        f.close()


######################################################################
### uname
######################################################################

try:
    conf['system'] = os.environ['UNAME']
except KeyError:
    conf['system'] = platform.system()

conf['Linux?'  ] = (conf['system'] == 'Linux')
conf['Solaris?'] = (conf['system'] == 'SunOS')
conf['Darwin?' ] = (conf['system'] == 'Darwin')
conf['FreeBSD?'] = (conf['system'] == 'FreeBSD')
conf['CYGWIN?' ] = (conf['system'] == 'CYGWIN')

conf['machine'] = platform.machine()

conf['processor'] = platform.processor()

conf['Intel?'] = (platform.machine() in ('i386', 'i486', 'i586', 'i686', 'x86_64',
                                         'AMD64', 'i86pc'))
conf['PPC?']   = (platform.processor() == 'powerpc')
conf['SPARC?'] = (platform.processor() == 'sparc')


######################################################################
### bit width
######################################################################

conf['bits'] = platform.architecture()[0]

if os.environ.get('SAGE64', 'no') == 'yes':
    assert conf['bits'] == '64bit', 'SAGE64=yes on a 32-bit system!'
    conf['64bit?'] = True
else:
    conf['64bit?'] = (conf['bits'] == '64bit')

conf['32bit?'] = not conf['64bit?']


######################################################################
### fortran compiler
######################################################################

fortran_version = try_run('sage_fortran --version')
if fortran_version is None:
    print 'Cannot execute fortran compiler (sage_fortran)!'
    sys.exit(3)
if 'G95' in fortran_version:
    conf['fortran'] = 'g95'
elif 'GNU Fortran' in fortran_version:
    conf['fortran'] = 'gfortran'
else:
    print 'Unknown fortran compiler version: '+fortran_version
    conf['fortran'] = None


conf['fortran_g95?'] = (conf['fortran'] == 'g95')
conf['fortran_GNU?'] = (conf['fortran'] == 'gfortran')


if conf['fortran_g95?']:
    g95_dir = glob.glob(os.environ['SAGE_LOCAL']+'/lib/gcc-lib/*/*')
    assert len(g95_dir)==1, 'Could not find G95 dir: '+str(g95_dir)
    conf['fortran_g95_dir'] = g95_dir[0]

######################################################################
### linker
######################################################################

# Note: Solaris linker does not accept --version
# ignore error code as 'ld -V' likes to error on some Solaris versions
ld_version = try_run('ld  -V', ignore=True)
if ld_version is None:
    print 'Cannot execute ld!'
    sys.exit(3)
if 'GNU' in ld_version:
    conf['ld'] = 'GNU'
elif 'Solaris' in ld_version:
    conf['ld'] = 'Solaris'
else:
    print 'Unknown linker: '+ld_version
    conf['ld'] = None

conf['linker_GNU?'] = (conf['ld'] == 'GNU')
conf['linker_Solaris?'] = (conf['ld'] == 'Solaris')


if conf['Solaris?'] and conf['linker_GNU?']:
    print "WARNING: You are using the GNU linker from 'binutils'"
    print "Generally it is considered better to use the Sun linker"
    print "but Sage has been built on Solaris using the GNU linker"
    print "although that was a very old version of Sage, which"
    print "never passes all the Sage test-suite."



######################################################################
### paths, files, and environment variables
######################################################################

conf['SPKG_DIR'] = os.getcwd()
conf['SAGE_LOCAL'] = os.environ['SAGE_LOCAL']


######################################################################
### The end: print configuration
######################################################################

print "Configuration:"
for key, value in conf.items():
    print '    '+str(key)+': '+str(value)




