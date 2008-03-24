import os, signal, sys, time, thread, threading, tempfile

from build.all import *
from sagebuild import buildsage, build_sage_clean
from clib import buildclib, build_clib_clean
argv = sys.argv
numthreads = 2
try:
    numthreads = int(os.environ['SAGE_BUILD_THREADS'])
except:
    pass
verbose=1
OM = optionmanager.OM
taskmanager.TM = taskmanager.taskmanager(numthreads)
TM = taskmanager.TM
env = None

def handle_argv():
    if verbose>20:
        print "Handling Argv"
    del argv[0]

class sage_env(enviroment):
    def __init__(self):
        enviroment.__init__(self)
        self.choose_blas()
        self.setup_options()
        self.setup_defaults()

    def choose_blas(self):
        ## Choose cblas library -- note -- make sure to update sage/misc/cython.py
        ## if you change this!!
        if os.environ.has_key('SAGE_BLAS'):
            self.options['BLAS']=os.environ['SAGE_BLAS']
            self.options['BLAS2']=os.environ['SAGE_BLAS']
        elif os.path.exists('%s/lib/libatlas.so'%os.environ['SAGE_LOCAL']):
            self.options['BLAS']='cblas'
            self.options['BLAS2']='atlas'
        elif os.path.exists('/usr/lib/libcblas.dylib') or \
             os.path.exists('/usr/lib/libcblas.so'):
            self.options['BLAS']='cblas'
            self.options['BLAS2']='atlas'
        elif os.path.exists('/usr/lib/libblas.dll.a'):
            self.options['BLAS']='gslcblas'
            self.options['BLAS2']='gslcblas'
        else:
            # This is very slow  (?), but *guaranteed* to be available.
            self.options['BLAS']='gslcblas'
            self.options['BLAS2']='gslcblas'

    def setup_options(self):
        if not os.environ.has_key('SAGE_ROOT'):
            print "    ERROR: The environment variable SAGE_ROOT must be defined."
            sys.exit(1)
        else:
            SAGE_ROOT  = os.environ['SAGE_ROOT']
            self.options['SAGE_ROOT'] = SAGE_ROOT
            SAGE_LOCAL = self.options['SAGE_ROOT'] + '/local/'
            self.options['SAGE_LOCAL'] = SAGE_LOCAL
            SAGE_DEVEL = self.options['SAGE_ROOT'] + '/devel/'
            self.options['SAGE_DEVEL'] = SAGE_DEVEL
        #Change to SAGE_ROOT
        os.chdir(SAGE_ROOT)
        if not os.environ.has_key('SAGE_VERSION'):
            self.options['SAGE_VERSION']=0
        else:
            self.options['SAGE_VERSION'] = os.environ['SAGE_VERSION']
        self.options['SITE_PACKAGES'] = '%slib/python/site-packages/'%self.options['SAGE_LOCAL']
        if not os.path.exists(self.options['SITE_PACKAGES']):
            self.options['SITE_PACKAGES'] = '%slib/python2.5/site-packages/'%self.options['SAGE_LOCAL']
            if not os.path.exists(self.options['SITE_PACKAGES']):
                self.options['SITE_PACKAGES'] = '%slib/python2.4/site-packages/'%self.options['SAGE_LOCAL']
                if not os.path.exists(self.options['SITE_PACKAGES']) and os.environ.has_key('SAGE_DEBIAN'):
                    self.options['SITE_PACKAGES'] = '/usr/lib/python2.5/site-packages/'
                if not os.path.exists(self.options['SITE_PACKAGES']):
                    raise RuntimeError, "Unable to find site-packages directory (see setup.py file in sage python code)."
        if not self.options.has_key('SAGE_DEBIAN'):
            self.options['SAGE_DEBIAN'] = False
        self.options['SITE_PACKAGES_REL'] = self.options['SITE_PACKAGES'][len(self.options['SAGE_LOCAL']):]
        try:
            self.options['PLATFORM'] = os.environ['PLATFORM']
        except:
            pass
    def setup_defaults(self):
        self.defaults['include_dirs'] = ['%sinclude'%self.options['SAGE_LOCAL'], \
        '%sinclude/csage'%self.options['SAGE_LOCAL'], \
        '%ssage/sage/ext'%self.options['SAGE_DEVEL'], \
        '%sinclude/python2.5'%self.options['SAGE_LOCAL']]

        if self.options['SAGE_DEBIAN']:
            self.debian_include_dirs=["/usr/include","/usr/include/numpy","/usr/include/FLINT","/usr/include/givaro", "/usr/include/gsl","/usr/include/fplll","/usr/include/eclib","/usr/include/gmp++","/usr/include/linbox","/usr/include/NTL","/usr/include/pari","/usr/include/qd","/usr/include/singular","/usr/include/singular/singular","/usr/include/symmetrica","/usr/include/polybori"]
            self.defaults['include_dirs'].expand(self.debian_include_dirs)
        else:
            self.debian_include_dirs=[]
        self.defaults['prelibraries'] = ['csage']
        self.defaults['postlibraries'] = ['stdc++', 'ntl']
        if self.options['SAGE_DEBIAN']:
            self.defaults['library_dirs'] =  ['/usr/lib','/usr/lib/eclib','/usr/lib/singular','/usr/lib/R/lib','%s/lib' % self.options['SAGE_LOCAL'] ]
        else:
            self.defaults['library_dirs'] = ['%slib' % self.options['SAGE_LOCAL'] ]

class debian_env(sage_env):
    def __init__(self):
        self.options = { }
        self.options['SAGE_DEBIAN'] = True
        sage_env.__init__(self)



def default_option(env,args):
    os.system('sage')

def help_option(env,args):
    outtext = """
  --coverage <files>  -- give info about doctest coverage of files
  --help, -h            -- print this help message
  --version, -v  -- print the SAGE version
  --notebook [options], -n -- start the SAGE notebook (options are
                   the same as to the notebook command in SAGE)
  --inotebook [options] -- start the *insecure* SAGE notebook
  --install [packages], -i -- install the given SAGE packages
  --min ...      -- do not populate global namespace (must be first option)
  --optional     -- list all optional packages that can be downloaded
  --test <files|dir>, -t -- test examples in .py, .pyx, .sage or .tex files
        options to '-t':   --verbose, --gdb, --long, --optional,
                  --valgrind, --memcheck, --massif, --cachegrind,
                  --callgrind, --omega
  --upgrade      -- download, build and install latest non-optional SAGE packages
  --advanced     -- list more options
  --command cmd, -c        -- Evaluates cmd as sage code
  file.[sage|py|spyx] [options] -- run given .sage, .py or .spyx files"""
    print outtext

def advanced_option(env, args):
    pass


def put_header(env, args):
    outtext = """
-----------------------------------------------------------
| SAGE: Software for Algebra and Geometry Experimentation |
-----------------------------------------------------------
"""
    print outtext

opt_sage_build = False
opt_sage_clib = False
opt_sage_all = False
opt_TM_go = False


def build_option(env,args):
    global opt_sage_build, opt_sage_clib, opt_sage_all, opt_TM_go
    if verbose>10:
        print "Building sage"
    opt_sage_build=True
    opt_TM_go=True

def build_clib_option(env, args):
    global opt_sage_build, opt_sage_clib, opt_sage_all, opt_TM_go
    if verbose>10:
        print "Building clib"
    opt_sage_clib = True
    opt_TM_go = True

def all_option(env, args):
    global opt_sage_build, opt_sage_clib, opt_sage_all, opt_TM_go
    if verbose>10:
        print "Building all"
    opt_sage_build = True
    opt_sage_clib = True
    opt_sage_all = True
    opt_TM_go = True
    build_clib_clean(env)
    build_sage_clean(env)
    pass
#def __init__(self, longstr, args, func, helpstr, shortstr = None, help_levels = [0], greedy = False):
def init():
    OM.set_default(default_option)
    OM.add_option(option( "help", [], help_option, "print this help message", "h", help_levels = [0,1]))
    OM.add_option(option("build", [], build_option, "switch to and build Sage branch in devel/sage-branch", "b", help_levels = [1]))
    OM.add_option(option("all", [], all_option, "rebuild all Cython code and clib", "a", help_levels = [1]))
    OM.add_option(option("buildclib", [], build_clib_option, "rebuild clib", help_levels = [1]))
#    OM.add_option('t',test_option,0,True)
#    OM.add_option('test',test_option,0,True)
    handle_argv()
    if 'SAGE_DEBIAN' in os.environ:
        env = debian_env()
    else:
        env = sage_env()
    if verbose > 10:
        env.put()
    put_header(env, argv)
    OM.execute_options(argv, env)
    global opt_sage_build, opt_sage_clib, opt_sage_all, opt_TM_go
    if opt_sage_clib:
        gccc = GCC_compiler(env, define_macros = [ ("NDEBUG",None) ], libraries = ['pthread'], options = { "-O2":None, "-g":None, "-fno-strict-aliasing":None } )
        buildclib(env, gccc)
    if opt_sage_build:
        gccc = GCC_compiler(env, define_macros = [ ("NDEBUG",None) ], libraries = ['pthread'], options = { "-O2":None, "-g":None, "-fno-strict-aliasing":None } )
        buildsage(env, gccc)
    if opt_TM_go:
        TM.go()

init()





