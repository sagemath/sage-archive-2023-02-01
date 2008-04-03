import os, signal, sys, time, thread, threading, tempfile

from build.all import *
from sageenv import create_env
from sagebuild import buildsage, build_sage_clean
from clib import buildclib, build_clib_clean
argv = sys.argv
numthreads = 2
try:
    numthreads = int(os.environ['SAGE_BUILD_THREADS'])
except:
    pass
print 'Building with %s threads.' % numthreads
verbose=1
OM = optionmanager.OM
taskmanager.TM = taskmanager.taskmanager(numthreads)
TM = taskmanager.TM
env = None

def handle_argv():
    if verbose>20:
        print "Handling Argv"
    del argv[0]

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
    env = create_env()
    if verbose > 10:
        env.put()
    #Change to SAGE_ROOT
    os.chdir(env.options['SAGE_ROOT'])
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





