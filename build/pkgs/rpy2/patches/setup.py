
import os, os.path, sys, shutil, re, itertools
from distutils.core import setup, Extension


pack_name = 'rpy2'
pack_version = __import__('rpy').__version__

RHOMES = os.getenv('RHOMES')

if RHOMES is None:

    RHOMES = os.popen("R RHOME").readlines()
    if len(RHOMES) == 0:
        raise RuntimeError(
            "R_HOME not defined, and no R command in the PATH."
            )

    #Twist if 'R RHOME' spits out a warning
    if RHOMES[0].startswith("WARNING"):
        RHOMES = RHOMES[1]
    else:
        RHOMES = RHOMES[0]
    RHOMES = [RHOMES, ]
else:
    RHOMES = RHOMES.split(os.pathsep)


def get_rversion(RHOME):
    r_exec = os.path.join(RHOME, 'bin', 'R')
    # Twist if Win32
    if sys.platform == "win32":
        rp = os.popen3('"'+r_exec+'" --version')[2]
    else:
        rp = os.popen('"'+r_exec+'" --version')
    rversion = rp.readline()
    #Twist if 'R RHOME' spits out a warning
    if rversion.startswith("WARNING"):
        rversion = rp.readline()
    m = re.match('^R version ([^ ]+) .+$', rversion)
    rversion = m.groups()[0]
    rversion = rversion.split('.')
    rversion[0] = int(rversion[0])
    rversion[1] = int(rversion[1])
    return rversion

def cmp_version(x, y):
    if (x[0] < y[0]):
        return -1
    if (x[0] > y[0]):
        return 1
    if (x[0] == y[0]):
        if len(x) == 1 or len(y) == 1:
            return 0
        return cmp_version(x[1:], y[1:])

def get_rconfig(RHOME, about, allow_empty = False):
    r_exec = os.path.join(RHOME, 'bin', 'R')
    cmd = '"'+r_exec+'" CMD config '+about
    rp = os.popen(cmd)
    rconfig = rp.readline()
    #Twist if 'R RHOME' spits out a warning
    if rconfig.startswith("WARNING"):
        rconfig = rp.readline()
    rconfig = rconfig.strip()
    #sanity check of what is returned into rconfig
    rconfig_m = re.match('^(-L.+) (-l.+)$', rconfig)
    #cheap fix for the case -lblas is returned
    #FIXME: clean/unify that at one point
    if rconfig_m is None:
        rconfig_m = re.match('^(-l.+)$', rconfig)
    if rconfig_m is None:
        # MacOSX
        rconfig_m = re.match('^(-F.+) (-framework.+)$', rconfig)
    if rconfig_m is None:
        # MacOSX
        rconfig_m = re.match('^(-framework.+)$', rconfig)
    if rconfig_m is None:
        rconfig_m = re.match('^(-I.+)$', rconfig)

    if (rconfig_m is None):
        if allow_empty:
            return ()
        else:
            raise Exception(cmd + '\nreturned\n' + rconfig)
    else:
        return rconfig_m.groups()

rnewest = [0, 0, 0]
rversions = []
for RHOME in RHOMES:
    RHOME = RHOME.strip()
    rversion = get_rversion(RHOME)
    if cmp_version(rversion[:2], [2, 7]) == -1:
        raise Exception("R >= 2.7 required.")
    rversions.append(rversion)

def getRinterface_ext(RHOME, r_packversion):
    r_libs = [os.path.join(RHOME, 'lib'), os.path.join(RHOME, 'modules')]

    #FIXME: crude way (will break in many cases)
    #check how to get how to have a configure step
    define_macros = []

    if sys.platform == 'win32':
        define_macros.append(('Win32', 1))
    else:
        define_macros.append(('R_INTERFACE_PTRS', 1))

    define_macros.append(('CSTACK_DEFNS', 1))
    define_macros.append(('RIF_HAS_RSIGHAND', 1))

    # defines for debugging
    #define_macros.append(('RPY_DEBUG_PRESERVE', 1))
    #define_macros.append(('RPY_DEBUG_PROMISE', 1))
    #define_macros.append(('RPY_DEBUG_OBJECTINIT', 1))
    #define_macros.append(('RPY_DEBUG_CONSOLE', 1))

    include_dirs = get_rconfig(RHOME, '--cppflags')[0].split()
    for i, d in enumerate(include_dirs):
        if d.startswith('-I'):
           include_dirs[i] = d[2:]

    rinterface_ext = Extension(
            pack_name + '.rinterface.rinterface',
            [os.path.join('rpy', 'rinterface', 'array.c'),
             os.path.join('rpy', 'rinterface', 'r_utils.c'),
             os.path.join('rpy', 'rinterface', 'rinterface.c')],
            include_dirs = include_dirs +
                            [os.path.join('rpy', 'rinterface'),],
            libraries = ['R', ],
            library_dirs = r_libs,
            define_macros = define_macros,
            runtime_library_dirs = r_libs,
            #extra_compile_args=['-O0', '-g'],
            extra_link_args = get_rconfig(RHOME, '--ldflags') +\
                              get_rconfig(RHOME, 'LAPACK_LIBS',
                                          allow_empty = True)
                              #+\
                              #get_rconfig(RHOME, 'BLAS_LIBS'),
            )

    return rinterface_ext


rinterface_exts = []
rinterface_rversions = []

for rversion, RHOME in itertools.izip(rversions, RHOMES):
    RHOME = RHOME.strip()
    #if (cmp_version(rversion, rnewest) == 0):
    #r_packversion = None
    r_packversion = '%i%02i%s' %(rversion[0], rversion[1], rversion[2])
    ri_ext = getRinterface_ext(RHOME, r_packversion)
    rinterface_exts.append(ri_ext)
    rinterface_rversions.append(r_packversion)

pack_dir = {pack_name: 'rpy'}

import distutils.command.install
for scheme in distutils.command.install.INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

setup(name = pack_name,
      version = pack_version,
      description = "Python interface to the R language",
      url = "http://rpy.sourceforge.net",
      license = "(L)GPL",
      author = "Laurent Gautier <lgautier@gmail.com>",
      ext_modules = rinterface_exts,
      package_dir = pack_dir,
      packages = [pack_name,
                  pack_name+'.rlike',
                  pack_name+'.rlike.tests',
                  pack_name+'.robjects',
                  pack_name+'.robjects.tests'] + \
                 [pack_name + '.rinterface', pack_name + '.rinterface.tests'],
                 #[pack_name + '.rinterface_' + x for x in rinterface_rversions] + \
                 #[pack_name + '.rinterface_' + x + '.tests' for x in rinterface_rversions]
      classifiers = ['Programming Language :: Python',
                     'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
                     'Intended Audience :: Developers',
                     'Intended Audience :: Science/Research',
                     'Development Status :: 5 - Production/Stable'
                    ],
      data_files = [(os.path.join('rpy2','images'), [os.path.join('doc', 'source', 'rpy2_logo.png')])]
      )

