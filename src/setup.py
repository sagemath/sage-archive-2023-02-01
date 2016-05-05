#!/usr/bin/env python
from __future__ import print_function

import os, sys, time, errno, platform, subprocess, glob
from distutils.core import setup


def excepthook(*exc):
    """
    When an error occurs, display an error message similar to the error
    messages from ``sage-spkg``.

    In particular, ``build/make/install`` will recognize "sage" as a failed
    package, see :trac:`16774`.
    """
    stars = '*' * 72

    print(stars, file=sys.stderr)
    import traceback
    traceback.print_exception(*exc, file=sys.stderr)
    print(stars, file=sys.stderr)
    print("Error building the Sage library", file=sys.stderr)
    print(stars, file=sys.stderr)

    try:
        logfile = os.path.join(os.environ['SAGE_LOGS'],
                "sagelib-%s.log" % os.environ['SAGE_VERSION'])
    except:
        pass
    else:
        print("Please email sage-devel (http://groups.google.com/group/sage-devel)", file=sys.stderr)
        print("explaining the problem and including the relevant part of the log file", file=sys.stderr)
        print("  " + logfile, file=sys.stderr)
        print("Describe your computer, operating system, etc.", file=sys.stderr)
        print(stars, file=sys.stderr)

sys.excepthook = excepthook


#########################################################
### List of Extensions
###
### The list of extensions resides in module_list.py in
### the same directory as this file
#########################################################

from module_list import ext_modules, library_order, aliases
from sage.env import *

#########################################################
### Configuration
#########################################################

if len(sys.argv) > 1 and sys.argv[1] == "sdist":
    sdist = True
else:
    sdist = False

try:
    compile_result_dir = os.environ['XML_RESULTS']
    keep_going = True
except KeyError:
    compile_result_dir = None
    keep_going = False

# search for dependencies and add to gcc -I<path>
include_dirs = sage_include_directories(use_sources=True)

# Manually add -fno-strict-aliasing, which is needed to compile Cython
# and disappears from the default flags if the user has set CFLAGS.
extra_compile_args = [ "-fno-strict-aliasing" ]
extra_link_args = [ ]

DEVEL = False
if DEVEL:
    extra_compile_args.append('-ggdb')

#########################################################
### Testing related stuff
#########################################################

class CompileRecorder(object):

    def __init__(self, f):
        self._f = f
        self._obj = None

    def __get__(self, obj, type=None):
        # Act like a method...
        self._obj = obj
        return self

    def __call__(self, *args):
        t = time.time()
        try:
            if self._obj:
                res = self._f(self._obj, *args)
            else:
                res = self._f(*args)
        except Exception as ex:
            print(ex)
            res = ex
        t = time.time() - t

        errors = failures = 0
        if self._f is compile_command0:
            name = "cythonize." + args[0][1].name
            failures = int(bool(res))
        else:
            name = "gcc." + args[0][1].name
            errors = int(bool(res))
        if errors or failures:
            type = "failure" if failures else "error"
            failure_item = """<%(type)s/>""" % locals()
        else:
            failure_item = ""
        output = open("%s/%s.xml" % (compile_result_dir, name), "w")
        output.write("""
            <?xml version="1.0" ?>
            <testsuite name="%(name)s" errors="%(errors)s" failures="%(failures)s" tests="1" time="%(t)s">
            <testcase classname="%(name)s" name="compile">
            %(failure_item)s
            </testcase>
            </testsuite>
        """.strip() % locals())
        output.close()
        return res

if compile_result_dir:
    record_compile = CompileRecorder
else:
    record_compile = lambda x: x

# Remove (potentially invalid) star import caches
import sage.misc.lazy_import_cache
if os.path.exists(sage.misc.lazy_import_cache.get_cache_file()):
    os.unlink(sage.misc.lazy_import_cache.get_cache_file())


######################################################################
# CODE for generating C/C++ code from Cython and doing dependency
# checking, etc.  In theory distutils would run Cython, but I don't
# trust it at all, and it won't have the more sophisticated dependency
# checking that we need.
######################################################################

# Do not put all, but only the most common libraries and their headers
# (that are likely to change on an upgrade) here:
# [At least at the moment. Make sure the headers aren't copied with "-p",
# or explicitly touch them in the respective spkg's spkg-install.]
lib_headers = { "gmp":     [ os.path.join(SAGE_INC, 'gmp.h') ],   # cf. #8664, #9896
                "gmpxx":   [ os.path.join(SAGE_INC, 'gmpxx.h') ],
                "ntl":     [ os.path.join(SAGE_INC, 'NTL', 'config.h') ]
              }

# In the loop below, don't append to any list, since many of these
# lists are actually identical Python objects. For every list, we need
# to write (at least the first time):
#
#   list = list + [foo]
#
for m in ext_modules:
    # Make everything depend on *this* setup.py file
    m.depends = m.depends + [__file__]

    # Add dependencies for the libraries
    for lib in lib_headers:
        if lib in m.libraries:
            m.depends += lib_headers[lib]

    m.extra_compile_args = m.extra_compile_args + extra_compile_args
    m.extra_link_args = m.extra_link_args + extra_link_args
    m.library_dirs = m.library_dirs + [os.path.join(SAGE_LOCAL, "lib")]
    m.include_dirs = m.include_dirs + include_dirs


#############################################
###### Parallel Cython execution
#############################################

def run_command(cmd):
    """
    INPUT:

    - ``cmd`` -- a string; a command to run

    OUTPUT: prints ``cmd`` to the console and then runs
    ``os.system(cmd)``.
    """
    print(cmd)
    sys.stdout.flush()
    return os.system(cmd)

def apply_func_progress(p):
    """
    Given a triple p consisting of a function, value and a string,
    output the string and apply the function to the value.

    The string could for example be some progress indicator.

    This exists solely because we can't pickle an anonymous function
    in execute_list_of_commands_in_parallel below.
    """
    sys.stdout.write(p[2])
    sys.stdout.flush()
    return p[0](p[1])

def execute_list_of_commands_in_parallel(command_list, nthreads):
    """
    Execute the given list of commands, possibly in parallel, using
    ``nthreads`` threads.  Terminates ``setup.py`` with an exit code
    of 1 if an error occurs in any subcommand.

    INPUT:

    - ``command_list`` -- a list of commands, each given as a pair of
       the form ``[function, argument]`` of a function to call and its
       argument

    - ``nthreads`` -- integer; number of threads to use

    WARNING: commands are run roughly in order, but of course successive
    commands may be run at the same time.
    """
    # Add progress indicator strings to the command_list
    N = len(command_list)
    progress_fmt = "[{:%i}/{}] " % len(str(N))
    for i in range(N):
        progress = progress_fmt.format(i+1, N)
        command_list[i] = command_list[i] + (progress,)

    from multiprocessing import Pool
    import fpickle_setup #doing this import will allow instancemethods to be pickable
    # map_async handles KeyboardInterrupt correctly if an argument is
    # given to get().  Plain map() and apply_async() do not work
    # correctly, see Trac #16113.
    pool = Pool(nthreads)
    result = pool.map_async(apply_func_progress, command_list, 1).get(99999)
    pool.close()
    pool.join()
    process_command_results(result)

def process_command_results(result_values):
    error = None
    for r in result_values:
        if r:
            print("Error running command, failed with status %s."%r)
            if not keep_going:
                sys.exit(1)
            error = r
    if error:
        sys.exit(1)

def execute_list_of_commands(command_list):
    """
    INPUT:

    - ``command_list`` -- a list of strings or pairs

    OUTPUT:

    For each entry in command_list, we attempt to run the command.
    If it is a string, we call ``os.system()``. If it is a pair [f, v],
    we call f(v).

    If the environment variable :envvar:`SAGE_NUM_THREADS` is set, use
    that many threads.
    """
    t = time.time()
    # Determine the number of threads from the environment variable
    # SAGE_NUM_THREADS, which is set automatically by sage-env
    try:
        nthreads = int(os.environ['SAGE_NUM_THREADS'])
    except KeyError:
        nthreads = 1

    # normalize the command_list to handle strings correctly
    command_list = [ [run_command, x] if isinstance(x, str) else x for x in command_list ]

    # No need for more threads than there are commands, but at least one
    nthreads = min(len(command_list), nthreads)
    nthreads = max(1, nthreads)

    def plural(n,noun):
        if n == 1:
            return "1 %s"%noun
        return "%i %ss"%(n,noun)

    print("Executing %s (using %s)"%(plural(len(command_list),"command"), plural(nthreads,"thread")))
    execute_list_of_commands_in_parallel(command_list, nthreads)
    print("Time to execute %s: %.2f seconds."%(plural(len(command_list),"command"), time.time() - t))


########################################################################
##
## Parallel gcc execution
##
## This code is responsible for making distutils dispatch the calls to
## build_ext in parallel. Since distutils doesn't seem to do this by
## default, we create our own extension builder and override the
## appropriate methods.  Unfortunately, in distutils, the logic of
## deciding whether an extension needs to be recompiled and actually
## making the call to gcc to recompile the extension are in the same
## function. As a result, we can't just override one function and have
## everything magically work. Instead, we split this work between two
## functions. This works fine for our application, but it means that
## we can't use this modification to make the other parts of Sage that
## build with distutils call gcc in parallel.
##
########################################################################

from distutils.command.build_ext import build_ext
from distutils.command.install import install
from distutils.dep_util import newer_group
from distutils import log

class sage_build_ext(build_ext):

    def build_extensions(self):

        from distutils.debug import DEBUG

        if DEBUG:
            print("self.compiler.compiler:")
            print(self.compiler.compiler)
            print("self.compiler.compiler_cxx:")
            print(self.compiler.compiler_cxx) # currently not used
            print("self.compiler.compiler_so:")
            print(self.compiler.compiler_so)
            print("self.compiler.linker_so:")
            print(self.compiler.linker_so)
            # There are further interesting variables...


        # At least on MacOS X, the library dir of the *original* Sage
        # installation is "hard-coded" into the linker *command*, s.t.
        # that directory is always searched *first*, which causes trouble
        # after the Sage installation has been moved (or its directory simply
        # been renamed), especially in conjunction with upgrades (cf. #9896).
        # (In principle, the Python configuration should be modified on
        # Sage relocations as well, but until that's done, we simply fix
        # the most important.)
        # Since the following is performed only once per call to "setup",
        # and doesn't hurt on other systems, we unconditionally replace *any*
        # library directory specified in the (dynamic) linker command by the
        # current Sage library directory (if it doesn't already match that),
        # and issue a warning message:

        if True or sys.platform[:6]=="darwin":

            sage_libdir = os.path.realpath(SAGE_LOCAL+"/lib")
            ldso_cmd = self.compiler.linker_so # a list of strings, like argv

            for i in range(1, len(ldso_cmd)):

                if ldso_cmd[i][:2] == "-L":
                    libdir = os.path.realpath(ldso_cmd[i][2:])
                    self.debug_print(
                      "Library dir found in dynamic linker command: " +
                      "\"%s\"" % libdir)
                    if libdir != sage_libdir:
                        self.compiler.warn(
                          "Replacing library search directory in linker " +
                          "command:\n  \"%s\" -> \"%s\"\n" % (libdir,
                                                              sage_libdir))
                        ldso_cmd[i] = "-L"+sage_libdir

        if DEBUG:
            print("self.compiler.linker_so (after fixing library dirs):")
            print(self.compiler.linker_so)


        # First, sanity-check the 'extensions' list
        self.check_extensions_list(self.extensions)

        import time
        t = time.time()

        compile_commands = []
        for ext in self.extensions:
            need_to_compile, p = self.prepare_extension(ext)
            if need_to_compile:
                compile_commands.append((record_compile(self.build_extension), p))

        execute_list_of_commands(compile_commands)

        print("Total time spent compiling C/C++ extensions: %.2f seconds." % (time.time() - t))

    def prepare_extension(self, ext):
        sources = ext.sources
        if sources is None or not isinstance(sources, (list, tuple)):
            raise DistutilsSetupError(("in 'ext_modules' option (extension '%s'), " +
                   "'sources' must be present and must be " +
                   "a list of source filenames") % ext.name)
        sources = list(sources)

        fullname = self.get_ext_fullname(ext.name)
        if self.inplace:
            # ignore build-lib -- put the compiled extension into
            # the source tree along with pure Python modules

            modpath = string.split(fullname, '.')
            package = string.join(modpath[0:-1], '.')
            base = modpath[-1]

            build_py = self.get_finalized_command('build_py')
            package_dir = build_py.get_package_dir(package)
            ext_filename = os.path.join(package_dir,
                                        self.get_ext_filename(base))
            relative_ext_filename = self.get_ext_filename(base)
        else:
            ext_filename = os.path.join(self.build_lib,
                                        self.get_ext_filename(fullname))
            relative_ext_filename = self.get_ext_filename(fullname)

        # while dispatching the calls to gcc in parallel, we sometimes
        # hit a race condition where two separate build_ext objects
        # try to create a given directory at the same time; whoever
        # loses the race then seems to throw an error, saying that
        # the directory already exists. so, instead of fighting to
        # fix the race condition, we simply make sure the entire
        # directory tree exists now, while we're processing the
        # extensions in serial.
        relative_ext_dir = os.path.split(relative_ext_filename)[0]
        prefixes = ['', self.build_lib, self.build_temp]
        for prefix in prefixes:
            path = os.path.join(prefix, relative_ext_dir)
            try:
                os.makedirs(path)
            except OSError as e:
                assert e.errno==errno.EEXIST, 'Cannot create %s.' % path
        depends = sources + ext.depends
        if not (self.force or newer_group(depends, ext_filename, 'newer')):
            log.debug("skipping '%s' extension (up-to-date)", ext.name)
            need_to_compile = False
        elif getattr(ext, "skip_build", False):
            log.debug("skipping '%s' extension (optional)", ext.name)
            need_to_compile = False
        else:
            log.info("building '%s' extension", ext.name)
            need_to_compile = True

        # If we need to compile, adjust the given extension
        if need_to_compile:
            libs = ext.libraries
            if ext.language == 'c++' and 'stdc++' not in libs:
                libs = libs + ['stdc++']

            # Sort libraries according to library_order
            ext.libraries = sorted(libs, key=lambda x: library_order.get(x, 0))

        return need_to_compile, (sources, ext, ext_filename)

    def build_extension(self, p):

        sources, ext, ext_filename = p

        # First, scan the sources for SWIG definition files (.i), run
        # SWIG on 'em to create .c files, and modify the sources list
        # accordingly.
        sources = self.swig_sources(sources, ext)

        # Next, compile the source code to object files.

        # XXX not honouring 'define_macros' or 'undef_macros' -- the
        # CCompiler API needs to change to accommodate this, and I
        # want to do one thing at a time!

        # Two possible sources for extra compiler arguments:
        #   - 'extra_compile_args' in Extension object
        #   - CFLAGS environment variable (not particularly
        #     elegant, but people seem to expect it and I
        #     guess it's useful)
        # The environment variable should take precedence, and
        # any sensible compiler will give precedence to later
        # command line args.  Hence we combine them in order:
        extra_args = ext.extra_compile_args or []

        macros = ext.define_macros[:]
        for undef in ext.undef_macros:
            macros.append((undef,))

        objects = self.compiler.compile(sources,
                                        output_dir=self.build_temp,
                                        macros=macros,
                                        include_dirs=ext.include_dirs,
                                        debug=self.debug,
                                        extra_postargs=extra_args,
                                        depends=ext.depends)

        # XXX -- this is a Vile HACK!
        #
        # The setup.py script for Python on Unix needs to be able to
        # get this list so it can perform all the clean up needed to
        # avoid keeping object files around when cleaning out a failed
        # build of an extension module.  Since Distutils does not
        # track dependencies, we have to get rid of intermediates to
        # ensure all the intermediates will be properly re-built.
        #
        self._built_objects = objects[:]

        # Now link the object files together into a "shared object" --
        # of course, first we have to figure out all the other things
        # that go into the mix.
        if ext.extra_objects:
            objects.extend(ext.extra_objects)
        extra_args = ext.extra_link_args or []

        # Detect target language, if not provided
        language = ext.language or self.compiler.detect_language(sources)

        self.compiler.link_shared_object(
            objects, ext_filename,
            libraries=self.get_libraries(ext),
            library_dirs=ext.library_dirs,
            runtime_library_dirs=ext.runtime_library_dirs,
            extra_postargs=extra_args,
            export_symbols=self.get_export_symbols(ext),
            debug=self.debug,
            build_temp=self.build_temp,
            target_lang=language)




#############################################
###### Cythonize
#############################################

def run_cythonize():
    from Cython.Build import cythonize
    import Cython.Compiler.Options
    import Cython.Compiler.Main

    # Sage uses these directives (mostly for historical reasons).
    Cython.Compiler.Options.embed_pos_in_docstring = True
    Cython.Compiler.Options.directive_defaults['autotestdict'] = False
    Cython.Compiler.Options.directive_defaults['cdivision'] = True
    Cython.Compiler.Options.directive_defaults['fast_getattr'] = True
    # The globals() builtin in Cython was fixed to return to the current scope,
    # but Sage relies on the broken behavior of returning to the nearest
    # enclosing Python scope (e.g. to perform variable injection).
    Cython.Compiler.Options.old_style_globals = True

    debug = False
    if os.environ.get('SAGE_DEBUG', None) != 'no':
        print('Enabling Cython debugging support')
        debug = True
        Cython.Compiler.Main.default_options['gdb_debug'] = True
        Cython.Compiler.Main.default_options['output_dir'] = 'build'

    profile = False
    if os.environ.get('SAGE_PROFILE', None) == 'yes':
        print('Enabling Cython profiling support')
        profile = True

    # Disable Cython caching (the cache is stored in ~/.cycache which is
    # Cython's default).
    # Cycache is currently too broken to use reliably:
    # http://trac.sagemath.org/ticket/17851
    Cython.Compiler.Main.default_options['cache'] = False

    force = True
    version_file = os.path.join(os.path.dirname(__file__), '.cython_version')
    version_stamp = '\n'.join([
        'cython version: ' + str(Cython.__version__),
        'debug: ' + str(debug),
        'profile: ' + str(profile),
    ""])
    if os.path.exists(version_file) and open(version_file).read() == version_stamp:
        force = False

    global ext_modules
    ext_modules = cythonize(
        ext_modules,
        nthreads=int(os.environ.get('SAGE_NUM_THREADS', 0)),
        build_dir=SAGE_CYTHONIZED,
        force=force,
        aliases=aliases,
        compiler_directives={
            'embedsignature': True,
            'profile': profile,
        })

    open(version_file, 'w').write(version_stamp)


print("Updating Cython code....")
t = time.time()
run_cythonize()
print("Finished Cythonizing, time: %.2f seconds." % (time.time() - t))


#########################################################
### Discovering Sources
#########################################################

print("Discovering Python/Cython source code....")
t = time.time()
from sage_setup.find import find_python_sources, find_extra_files
python_packages, python_modules = find_python_sources(
    SAGE_SRC, ['sage', 'sage_setup'])
python_data_files = find_extra_files(python_packages,
    ".", SAGE_CYTHONIZED, SAGE_LIB, ["ntlwrap.cpp"])
print("Discovered Python/Cython sources, time: %.2f seconds." % (time.time() - t))


#########################################################
### Clean
#########################################################

print('Cleaning up stale installed files....')
t = time.time()
from sage_setup.clean import clean_install_dir
output_dirs = SITE_PACKAGES + glob.glob(os.path.join(SAGE_SRC, 'build', 'lib*'))
for output_dir in output_dirs:
    print('- cleaning {0}'.format(output_dir))
    clean_install_dir(output_dir, python_packages, python_modules,
            ext_modules, python_data_files)
print('Finished cleaning, time: %.2f seconds.' % (time.time() - t))


#########################################################
### Install also Jupyter kernel spec
#########################################################

# We cannot just add the installation of the kernel spec to data_files
# since the file is generated, not copied.
class sage_install(install):
    def run(self):
        install.run(self)
        self.install_kernel_spec()

    def install_kernel_spec(self):
        from sage.repl.ipython_kernel.install import SageKernelSpec
        SageKernelSpec.update()


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
      packages    = python_packages,
      data_files  = python_data_files,
      scripts = [],
      cmdclass = dict(build_ext=sage_build_ext, install=sage_install),
      ext_modules = ext_modules)
