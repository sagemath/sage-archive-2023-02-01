#!/usr/bin/env python

import os, sys, time
from distutils.core import setup
from distutils.extension import Extension

#########################################################
### Configuration
#########################################################

if len(sys.argv) > 1 and sys.argv[1] == "sdist":
    sdist = True
else:
    sdist = False

if not os.environ.has_key('SAGE_ROOT'):
    print "    ERROR: The environment variable SAGE_ROOT must be defined."
    sys.exit(1)
else:
    SAGE_ROOT  = os.environ['SAGE_ROOT']
    SAGE_LOCAL = SAGE_ROOT + '/local/'
    SAGE_DEVEL = SAGE_ROOT + '/devel/'

if not os.environ.has_key('SAGE_VERSION'):
    SAGE_VERSION=0
else:
    SAGE_VERSION = os.environ['SAGE_VERSION']

SITE_PACKAGES = '%s/lib/python/site-packages/'%SAGE_LOCAL
if not os.path.exists(SITE_PACKAGES):
    SITE_PACKAGES = '%s/lib/python2.5/site-packages/'%SAGE_LOCAL
    if not os.path.exists(SITE_PACKAGES):
        SITE_PACKAGES = '%s/lib/python2.4/site-packages/'%SAGE_LOCAL
        if not os.path.exists(SITE_PACKAGES) and os.environ.has_key('SAGE_DEBIAN'):
            SITE_PACKAGES = '%s/lib/python2.5/site-packages/'%SAGE_LOCAL
            os.system('mkdir -p "%s"'%SITE_PACKAGES)
        if not os.path.exists(SITE_PACKAGES):
            raise RuntimeError, "Unable to find site-packages directory (see setup.py file in sage python code)."

if not os.path.exists('build/sage'):
    os.makedirs('build/sage')

sage_link = SITE_PACKAGES + '/sage'
if not os.path.islink(sage_link) or not os.path.exists(sage_link):
    os.system('rm -rf "%s"'%sage_link)
    os.system('cd %s; ln -sf ../../../../devel/sage/build/sage .'%SITE_PACKAGES)

include_dirs = ['%s/include'%SAGE_LOCAL, \
                '%s/include/csage'%SAGE_LOCAL, \
                ## this is included, but doesn't actually exist
                ## '%s/include/python'%SAGE_LOCAL, \
                '%s/sage/sage/ext'%SAGE_DEVEL]

extra_compile_args = [ ]

# comment these four lines out to turn on warnings from gcc
import distutils.sysconfig
NO_WARN = True
if NO_WARN and distutils.sysconfig.get_config_var('CC').startswith("gcc"):
    extra_compile_args.append('-w')

DEVEL = False
if DEVEL:
    extra_compile_args.append('-ggdb')



######################################################################
# CODE for generating C/C++ code from Cython and doing dependency
# checking, etc.  In theory distutils would run Cython, but I don't
# trust it at all, and it won't have the more sophisticated dependency
# checking that we need.
######################################################################


from module_list import ext_modules

for m in ext_modules:
    m.libraries = ['csage'] + m.libraries + ['stdc++', 'ntl']
    m.extra_compile_args += extra_compile_args
    if os.environ.has_key('SAGE_DEBIAN'):
        m.library_dirs += ['/usr/lib','/usr/lib/eclib','/usr/lib/singular','/usr/lib/R/lib','%s/lib' % SAGE_LOCAL]
    else:
        m.library_dirs += ['%s/lib' % SAGE_LOCAL]




#############################################
###### Parallel Cython execution
#############################################

def execute_list_of_commands_in_serial(command_list):
    """
    INPUT:
        command_list -- a list of commands (strings) to run at the shell using os.system

    OUTPUT:
        the given list of commands are all executed in serial
    """
    for cmd in command_list:
        print cmd
        r = os.system(cmd)
        if r != 0:
            print "Error running command, exited with status %s."%r
            sys.exit(1)

def run_command(cmd):
    """
    INPUT:
        cmd -- a string; a command to run

    OUTPUT:
        prints cmd to the console and then runs os.system
    """
    print cmd
    return os.system(cmd)

def execute_list_of_commands_in_parallel(command_list, ncpus):
    """
    INPUT:
        command_list -- a list of strings (commands)
        ncpus -- integer; number of cpus to use

    OUTPUT:
        Executes the given list of commands, possibly in parallel,
        using ncpus cpus.  Terminates setup.py with an exit code of 1
        if an error occurs in any subcommand.

    WARNING: commands are run roughly in order, but of course successive
    commands may be run at the same time.
    """
    print "Execute %s commands (using %s cpus)"%(len(command_list), min(len(command_list),ncpus))
    from processing import Pool
    p = Pool(ncpus)
    for r in p.imap(run_command, command_list):
        if r:
            print "Parallel build failed with status %s."%r
            sys.exit(1)

def number_of_cpus():
    """
    Try to determine the number of cpu's on this system.
    If successful return that number.  Otherwise return 0
    to indicate failure.

    OUTPUT:
        int
    """
    if hasattr(os, "sysconf") and os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"): # Linux and Unix
        n = os.sysconf("SC_NPROCESSORS_ONLN")
        if isinstance(n, int) and n > 0:
            return n
    try:
        return int(os.popen2("sysctl -n hw.ncpu")[1].read().strip())
    except:
        return 0

def execute_list_of_commands(command_list):
    """
    INPUT:
        command_list -- a list of strings (commands)
    OUTPUT:
        runs the given list of commands using os.system; on machines
        with more than 1 cpu the commands are run in parallel.
    """
    t = time.time()
    if not os.environ.has_key('MAKE'):
        ncpus = 1
    else:
        MAKE = os.environ['MAKE']
        z = [w[2:] for w in MAKE.split() if w.startswith('-j')]
        if len(z) == 0:  # no command line option
            ncpus = 1
        else:
            # Determine number of cpus from command line argument.
            # Also, use the OS to cap the number of cpus, in case
            # user annoyingly makes a typo and asks to use 10000
            # cpus at once.
            try:
                ncpus = int(z[0])
                n = 2*number_of_cpus()
                if n:  # prevent dumb typos.
                    ncpus = min(ncpus, n)
            except ValueError:
                ncpus = 1

    if ncpus > 1:
        # parallel version
        execute_list_of_commands_in_parallel(command_list, ncpus)
    else:
        # non-parallel version
        execute_list_of_commands_in_serial(command_list)
    print "Time to execute %s commands: %s seconds"%(len(command_list), time.time() - t)


#############################################
###### Dependency checking
#############################################


# matches any dependency
import re
dep_regex = re.compile(r'^ *(?:cimport +(\S+))|(?:from +(\S+) *cimport)|(?:include *[\'"]([^\'"]+)[\'"])', re.M)

class DependencyTree:
    """
    This class stores all the information about the dependencies of a set of
    Cython files. It uses a lot of caching so information only needs to be
    looked up once per build.
    """
    def __init__(self):
        self._last_parse = {}
        self._timestamps = {}
        self._deps = {}
        self._deps_all = {}
        self.root = "%s/devel/sage/" % SAGE_ROOT

    def __getstate__(self):
        """
        Used for pickling.

        Timestamps and deep dependencies may change between builds,
        so we don't want to save those.
        """
        state = dict(self.__dict__)
        state['_timestamps'] = {}
        state['_deps_all'] = {}
        return state

    def __setstate__(self, state):
        """
        Used for unpickling.
        """
        self.__dict__.update(state)
        self._timestamps = {}
        self._deps_all = {}
        self.root = "%s/devel/sage/" % SAGE_ROOT

    def timestamp(self, filename):
        """
        Look up the last modified time of a file, with caching.
        """
        if filename not in self._timestamps:
            try:
                self._timestamps[filename] = os.path.getmtime(filename)
            except OSError:
                self._timestamps[filename] = 0
        return self._timestamps[filename]

    def parse_deps(self, filename, verify=True):
        """
        Open a Cython file and extract all of its dependencies.

        INPUT:
            filename -- the file to parse
            verify   -- only return existing files (default True)

        OUTPUT:
            list of dependency files
        """
        # only parse cython files
        if filename[-4:] not in ('.pyx', '.pxd', '.pxi'):
            return []

        dirname = os.path.split(filename)[0]
        deps = set()
        if filename.endswith('.pyx'):
            pxd_file = filename[:-4] + '.pxd'
            if os.path.exists(pxd_file):
                deps.add(pxd_file)

        f = open(filename)
        for m in dep_regex.finditer(open(filename).read()):
            groups = m.groups()
            module = groups[0] or groups[1] # cimport or from ... cimport
            if module is not None:
                if '.' in module:
                    path = module.replace('.', '/') + '.pxd'
                else:
                    path = "%s/%s.pxd" % (dirname, module)
            else:
                path = '%s/%s'%(dirname, groups[2])
                if not os.path.exists(path):
                    path = groups[2]
            deps.add(os.path.normpath(path))
        f.close()
        return list(deps)

    def immediate_deps(self, filename):
        """
        Returns a list of files directly referenced by this file.
        """
        if (filename not in self._deps
                or self.timestamp(filename) < self._last_parse[filename]):
            self._deps[filename] = self.parse_deps(filename)
            self._last_parse[filename] = self.timestamp(filename)
        return self._deps[filename]

    def all_deps(self, filename, path=None):
        """
        Returns all files directly or indirectly referenced by this file.

        A recursive algorithm is used here to maximize caching, but it is
        still robust for circular cimports (via the path parameter).
        """
        if filename not in self._deps_all:
            circular = False
            deps = set([filename])
            if path is None:
                path = set([filename])
            else:
                path.add(filename)
            for f in self.immediate_deps(filename):
                if f not in path:
                    deps.update(self.all_deps(f, path))
                else:
                    circular = True
            path.remove(filename)
            if circular:
                return deps # Don't cache, as this may be incomplete
            else:
                self._deps_all[filename] = deps
        return self._deps_all[filename]

    def newest_dep(self, filename):
        """
        Returns the most recently modified file that filename depends on,
        along with its timestamp.
        """
        nfile = filename
        ntime = self.timestamp(filename)
        for f in self.all_deps(filename):
            if self.timestamp(f) > ntime:
                nfile = f
                ntime = self.timestamp(f)
        return nfile, ntime


#############################################
###### Build code
#############################################

def process_filename(f, m):
    base, ext = os.path.splitext(f)
    if ext == '.pyx':
        if m.language == 'c++':
            return base + '.cpp'
        else:
            return base + '.c'
    else:
        return f

def compile_cmd(f, m):
    """
    Given a .pyx file f, which is a part of module m, copy the
    file to SITE_PACKAGES, and return a string which will call
    Cython on it.
    """
    if f.endswith('.pyx'):
        # process cython file
        pyx_inst_file = '%s/%s'%(SITE_PACKAGES, f)
        retval = os.system('cp %s %s 2>/dev/null'%(f, pyx_inst_file))
        # we could do this more elegantly -- load the files, use
        # os.path.exists to check that they exist, etc. ... but the
        # *vast* majority of the time, the copy just works. so this is
        # just specializing for the most common use case.
        if retval:
            dirname, filename = os.path.split(pyx_inst_file)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            retval = os.system('cp %s %s 2>/dev/null'%(f, pyx_inst_file))
            if retval:
                raise OSError, "cannot copy %s to %s"%(f,pyx_inst_file)
        print "%s --> %s"%(f, pyx_inst_file)
        outfile = f[:-4]
        if m.language == 'c++':
            outfile += ".cpp"
        else:
            outfile += ".c"
        cmd = "python2.5 `which cython` --embed-positions --incref-local-binop -I%s -o %s %s"%(os.getcwd(), outfile, f)

    elif f.endswith(('.c','.cc','.cpp')):
        # process C/C++ file
        cmd = "touch %s"%f

    return cmd #"NEED TO COMPILE file " + f + " in module " + m.name

def compile_command_list(ext_modules, deps):
    """
    Computes a list of commands needed to compile and link the
    extension modules given in 'ext_modules'
    """
    queue_compile_high = []
    queue_compile_med = []
    queue_compile_low = []

    for m in ext_modules:
        new_sources = []
        for f in m.sources:
            if f.endswith('.pyx'):
                dep_file, dep_time = deps.newest_dep(f)
                dest_file = "%s/%s"%(SITE_PACKAGES, f)
                dest_time = deps.timestamp(dest_file)
                if dest_time < dep_time:
                    if dep_file == f:
                        print "Building modified file %s."%f
                        cmd = compile_cmd(f, m)
                        queue_compile_high.append(cmd)
                    elif dep_file == (f[:-4] + '.pxd'):
                        print "Building %s because it depends on %s."%(f, dep_file)
                        cmd = compile_cmd(f, m)
                        queue_compile_med.append(cmd)
                    else:
                        print "Building %s because it depends on %s."%(f, dep_file)
                        cmd = compile_cmd(f, m)
                        queue_compile_low.append(cmd)
            new_sources.append(process_filename(f, m))
        m.sources = new_sources
    # print "# compile high = ", len(queue_compile_high)
    # print queue_compile_high
    # print "# compile med =  ", len(queue_compile_med)
    # print queue_compile_med
    # print "# compile low =  ", len(queue_compile_low)
    # print queue_compile_low
    return queue_compile_high + queue_compile_med + queue_compile_low


import cPickle as pickle
CYTHON_DEPS_FILE='.cython_deps'

if not sdist:
    print "Updating Cython code...."
    t = time.time()
    try:
        f = open(CYTHON_DEPS_FILE)
        deps = pickle.load(open(CYTHON_DEPS_FILE))
        f.close()
    except:
        deps = DependencyTree()
    queue = compile_command_list(ext_modules, deps)
    execute_list_of_commands(queue)
    f = open(CYTHON_DEPS_FILE, 'w')
    pickle.dump(deps, f)
    f.close()
    print "Finished compiling Cython code (time = %s seconds)"%(time.time() - t)


#########################################################
### Distutils
#########################################################

code = setup(name = 'sage',

      version     =  SAGE_VERSION,

      description = 'Sage: Open Source Mathematics Software',

      license     = 'GNU Public License (GPL)',

      author      = 'William Stein et al.',

      author_email= 'http://groups.google.com/group/sage-support',

      url         = 'http://www.sagemath.org',

      packages    = ['sage',

                     'sage.algebras',

                     'sage.calculus',

                     'sage.catalogue',

                     'sage.categories',

                     'sage.coding',

                     'sage.combinat',
                     'sage.combinat.crystals',
                     'sage.combinat.designs',
                     'sage.combinat.sf',
                     'sage.combinat.root_system',
                     'sage.combinat.matrices',
                     'sage.combinat.posets',
                     'sage.combinat.species',

                     'sage.crypto',
                     'sage.crypto.mq',

                     'sage.databases',

                     'sage.ext',

                     'sage.finance',

                     'sage.functions',

                     'sage.geometry',

                     'sage.games',

                     'sage.gsl',

                     'sage.graphs',
                     'sage.graphs.base',

                     'sage.groups',
                     'sage.groups.abelian_gps',
                     'sage.groups.matrix_gps',
                     'sage.groups.perm_gps',
                     'sage.groups.perm_gps.partn_ref',

                     'sage.interfaces',

                     'sage.lfunctions',

                     'sage.libs',
                     'sage.libs.fplll',
                     'sage.libs.linbox',
                     'sage.libs.mwrank',
                     'sage.libs.ntl',
                     'sage.libs.flint',
                     'sage.libs.pari',
                     'sage.libs.singular',
                     'sage.libs.symmetrica',
                     'sage.libs.cremona',

                     'sage.logic',

                     'sage.matrix',
                     'sage.media',
                     'sage.misc',

                     'sage.modules',

                     'sage.modular',
                     'sage.modular.abvar',
                     'sage.modular.hecke',
                     'sage.modular.modform',
                     'sage.modular.modsym',
                     'sage.modular.ssmod',

                     'sage.monoids',

                     'sage.numerical',

                     'sage.plot',
                     'sage.plot.plot3d',

                     'sage.probability',

                     'sage.quadratic_forms',
                     'sage.quadratic_forms.genera',

                     'sage.rings',
                     'sage.rings.number_field',
                     'sage.rings.padics',
                     'sage.rings.polynomial',
                     'sage.rings.polynomial.padics',

                     'sage.tests',

                     'sage.sets',

                     'sage.stats',

                     'sage.stats.hmm',

                     'sage.symbolic',

                     'sage.parallel',

                     'sage.schemes',
                     'sage.schemes.generic',
                     'sage.schemes.jacobians',
                     'sage.schemes.plane_curves',
                     'sage.schemes.plane_quartics',
                     'sage.schemes.elliptic_curves',
                     'sage.schemes.hyperelliptic_curves',

                     'sage.server',
                     'sage.server.simple',
                     'sage.server.notebook',
                     'sage.server.notebook.compress',
                     'sage.server.wiki',
                     'sage.server.trac',

                     'sage.structure',
                     'sage.structure.proof',

                     'sage.dsage',
                     'sage.dsage.tests',
                     'sage.dsage.database',
                     'sage.dsage.database.tests',
                     'sage.dsage.server',
                     'sage.dsage.server.tests',
                     'sage.dsage.interface',
                     'sage.dsage.interface.tests',
                     'sage.dsage.errors',
                     'sage.dsage.twisted',
                     'sage.dsage.twisted.tests',
                     'sage.dsage.dist_functions',
                     'sage.dsage.dist_functions.tests',
                     'sage.dsage.misc',
                     'sage.dsage.misc.tests',
                     'sage.dsage.web',
                     'sage.dsage.scripts',
                     ],

      scripts = ['sage/dsage/scripts/dsage_worker.py',
                 'sage/dsage/scripts/dsage_setup.py',
                 'spkg-debian-maybe',
                ],

      data_files = [('dsage/web/static',
                    ['sage/dsage/web/static/dsage_web.css',
                     'sage/dsage/web/static/dsage_web.js',
                     'sage/dsage/web/static/jquery.js',
                     'sage/dsage/web/static/jquery.tablesorter.pack.js',
                     'sage/dsage/web/static/jquery.history.js',
                     'sage/dsage/web/static/jquery.form.js',
                     'sage/dsage/web/static/asc.gif',
                     'sage/dsage/web/static/desc.gif',
                     'sage/dsage/web/static/bg.gif',
                     'sage/dsage/README.html']),
                    ('dsage/web/',
                    ['sage/dsage/web/index.html'])],

      ext_modules = ext_modules,
      include_dirs = include_dirs)

