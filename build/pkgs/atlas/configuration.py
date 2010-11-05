######################################################################
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
######################################################################

import platform, os, sys, time, shutil, glob


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

def try_run(command):
    """
    Try to execute `command` and return its output as string. Return
    `None` if execution fails.
    """
    try:
        f = os.popen(command, 'r')
    except NameError:
        return None
    result = f.read()
    if f.close() is None:
        return result.strip()
    else:
        return None


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
        self.data = re.sub(find_regex, subs, self.data, count)
        return self

    def close(self):
        f = open(filename, 'w')
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

conf['Intel?'] = (platform.processor() in ('i386', 'x86_64'))
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

ld_version = try_run('ld  --version')
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
### The end: print configuration
######################################################################

print "Configuration:"
for key, value in conf.items():
    print '    '+str(key)+': '+str(value)




