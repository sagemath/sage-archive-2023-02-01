#!/usr/bin/env python

import shutil
from distutils.core import setup, Extension
import distutils.sysconfig
import os
import stat
import commands
import sys

gnofract4d_version = "3.6"

if float(sys.version[:3]) < 2.4:
    print "Sorry, you need Python 2.4 or higher to run Gnofract 4D."
    print "You have version %s. Please upgrade." % sys.version
    sys.exit(1)

# hack to use a different Python for building if an env var is set
# I use this to build python-2.2 and 2.4 RPMs.
build_version = os.environ.get("BUILD_PYTHON_VERSION")
build_python = os.environ.get("BUILD_PYTHON")

if build_version and build_python and sys.version[:3] != build_version:
    print sys.version[:3], build_version
    args = ["/usr/bin/python"] + sys.argv
    print "running other Python version %s with args: %s" % (build_python,args)
    os.execv(build_python, args)

from buildtools import my_bdist_rpm, my_build, my_build_ext, my_install_lib

# Extensions need to link against appropriate libs
# We use pkg-config to find the appropriate set of includes and libs

pkgs = "gconf-2.0"

def call_package_config(package,option,optional=False):
    '''invoke pkg-config, if it exists, to find the appropriate
    arguments for a library'''
    cmd = "pkg-config %s %s" % (package, option)
    (status,output) = commands.getstatusoutput(cmd)
    if status != 0:
        if optional:
            print >>sys.stderr, "Can't find '%s'" % package
            print >>sys.stderr, "Some functionality will be disabled"
            return []
        else:
            print >>sys.stderr, "Can't set up. Error running '%s'." % cmd
            print >>sys.stderr, output
            print >>sys.stderr, "Possibly you don't have one of these installed: '%s'." % package
            sys.exit(1)

    return output.split()

gconf_flags = call_package_config(pkgs,"--cflags",True)
gconf_libs =  call_package_config(pkgs,"--libs",True)

extra_macros = []
if gconf_flags != []:
    extra_macros.append(('GCONF_ENABLED',1))

png_flags = call_package_config("libpng", "--cflags", True)
if png_flags != []:
    extra_macros.append(('PNG_ENABLED', 1))
else:
    print "NO PNG HEADERS FOUND"

png_libs = call_package_config("libpng", "--libs", True)

jpg_lib = "jpeg"
if os.path.isfile("/usr/include/jpeglib.h"):
    extra_macros.append(('JPG_ENABLED', 1))
    jpg_libs = [ jpg_lib ]
else:
    print "NO JPEG HEADERS FOUND"
    jpg_libs = []

#not ready yet.
have_gmp = False # os.path.isfile("/usr/include/gmp.h")

# use currently specified compilers, not ones from when Python was compiled
# this is necessary for cross-compilation
compiler = os.environ.get("CC","gcc")
cxxcompiler = os.environ.get("CXX","g++")

fract4d_sources = [
    'fract4d/c/fract4dmodule.cpp',
    'fract4d/c/cmap.cpp',
    'fract4d/c/pointFunc.cpp',
    'fract4d/c/fractFunc.cpp',
    'fract4d/c/STFractWorker.cpp',
    'fract4d/c/MTFractWorker.cpp',
    'fract4d/c/image.cpp',
    'fract4d/c/imageIO.cpp',
    'fract4d/c/fract_stdlib.cpp'
    ]

# this is a hack to build 2 versions of the same extension.
# we want to create a standard fract4dc which doesn't depend on gmp
# and a second fract4dcgmp which does. These are both built from the
# same source files but in the latter case with USE_GMP defined
# This is so I can ship a single binary for gnofract4d which supports
# both users with GMP and users without it, by conditionally loading
# the appropriate extension
fract4d_gmp_sources = []
if(have_gmp):
    for sourcefile in fract4d_sources:
        # this particular part of the hack is so that each file gets
        # compiled twice
        gmp_sourcefile = sourcefile.replace(".cpp","_gmp.cpp")
        os.system("cp %s %s" % (sourcefile, gmp_sourcefile))
        fract4d_gmp_sources.append(gmp_sourcefile)

module_gmp = Extension(
    'fract4d.gmpy',
    sources = [
    'fract4d/gmpy/gmpy.c'
    ],
    libraries = ['gmp']
    )

defines = [ ('_REENTRANT',1),
            #('NO_CALC', 1),  # set this to not calculate the fractal
            #('DEBUG_CREATION',1), # debug spew for allocation of objects
            #('DEBUG_ALLOCATION',1), # debug spew for array handling
            ]
module_fract4dgmp = Extension(
    'fract4d.fract4dcgmp',
    sources = fract4d_gmp_sources,
    include_dirs = [
    'fract4d/c'
    ],
    libraries = [
    'stdc++', 'gmp'
    ] + jpg_libs,
    extra_compile_args = [
    '-Wall',
    ] + png_flags,
    extra_link_args = png_libs,
    define_macros = defines + [('USE_GMP',1)] + extra_macros,
    undef_macros = [ 'NDEBUG']
    )

module_fract4dc = Extension(
    'fract4d.fract4dc',
    sources = fract4d_sources,
    include_dirs = [
    'fract4d/c'
    ],
    libraries = [
    'stdc++'
    ] + jpg_libs,
    extra_compile_args = [
    '-Wall', '-O0'
    ] + png_flags,
    extra_link_args = png_libs,
    define_macros = defines + extra_macros,
    undef_macros = [ 'NDEBUG']
    )

module_cmap = Extension(
    'fract4d.fract4d_stdlib',
    sources = [
    'fract4d/c/cmap.cpp',
    'fract4d/c/image.cpp',
    'fract4d/c/fract_stdlib.cpp'
    ],
    include_dirs = [
    'fract4d/c'
    ],
    libraries = [
    'stdc++'
    ],
    define_macros = [ ('_REENTRANT', 1)]
    )

module_gui = Extension(
    'fract4dgui.fract4dguic',
    sources = [
    'fract4dgui/c/guicmodule.cpp',
    ],
    include_dirs = [
    'fract4dgui/c',
    'fract4d/c/'
    ],
    libraries = [
    'stdc++'
    ],
    extra_compile_args = gconf_flags,
    extra_link_args = gconf_libs,
    define_macros = [ ('_REENTRANT',1),
                      #('DEBUG_CREATION',1)
                      ] + extra_macros,
    undef_macros = [ 'NDEBUG']
    )

modules = [module_fract4dc, module_cmap, module_gui]

# SAGE -- no GUI #
modules = [module_fract4dc, module_cmap]


if have_gmp:
    modules.append(module_fract4dgmp)
    modules.append(module_gmp)

def get_files(dir,ext):
    return [ os.path.join(dir,x) for x in os.listdir(dir) if x.endswith(ext)]

setup (name = 'gnofract4d',
       version = gnofract4d_version,
       description = 'A program to draw fractals',
       long_description = \
'''Gnofract 4D is a fractal browser. It can generate many different fractals,
including some which are hybrids between the Mandelbrot and Julia sets,
and includes a Fractint-compatible parser for your own fractal formulas.''',
       author = 'Tim Whidbey',
       author_email = 'catenary@users.sourceforge.net',
       maintainer = 'Tim whidbey',
       maintainer_email = 'catenary@users.sourceforge.net',
       keywords = "fractal Mandelbrot Julia fractint chaos",
       url = 'http://gnofract4d.sourceforge.net/',
       packages = ['fract4d', 'fract4dgui', 'fractutils'],
       ext_modules = modules,
       scripts = ['gnofract4d'],
       data_files = [
           # color maps
           ('share/gnofract4d/maps',
            get_files("maps",".map") +
            get_files("maps",".cs") +
            get_files("maps", ".ugr")),

           # formulas
           ('share/gnofract4d/formulas',
            get_files("formulas","frm") +
            get_files("formulas", "ucl") +
            get_files("formulas", "uxf")),

           # documentation
           ('share/gnome/help/gnofract4d/C',
            get_files("doc/gnofract4d-manual/C", "xml")),
           ('share/gnome/help/gnofract4d/C/figures',
            get_files("doc/gnofract4d-manual/C/figures",".png")),
           ('share/gnome/help/gnofract4d/C',
            get_files("doc/gnofract4d-manual/C", "html")),
           ('share/gnome/help/gnofract4d/C/figures',
            get_files("doc/gnofract4d-manual/C/figures",".css")),

           #internal pixmaps
           ('share/pixmaps/gnofract4d',
            ['pixmaps/deepen_now.png',
             'pixmaps/explorer_mode.png']),

           # icon
           ('share/pixmaps',
            ['pixmaps/gnofract4d-logo.png']),

           # .desktop file
           ('share/applications', ['gnofract4d.desktop']),

           # MIME type registration
           ('share/mime/packages', ['gnofract4d-mime.xml']),

           # doc files
           ('share/doc/gnofract4d-%s/' % gnofract4d_version,
            ['COPYING', 'README']),
           ],
       cmdclass={
           "my_bdist_rpm": my_bdist_rpm.my_bdist_rpm,
           "build" : my_build.my_build,
           "my_build_ext" : my_build_ext.my_build_ext,
           "install_lib" : my_install_lib.my_install_lib
           }
       )

# I need to find the file I just built and copy it up out of the build
# location so it's possible to run without installing. Can't find a good
# way to extract the actual target directory out of distutils, hence
# this egregious hack

so_extension = distutils.sysconfig.get_config_var("SO")

lib_targets = {
    "fract4dguic" + so_extension : "fract4dgui",
    "fract4dc" + so_extension : "fract4d",
    "fract4d_stdlib" + so_extension : "fract4d",
    "fract4dcgmp" + so_extension : "fract4d",
    "gmpy" + so_extension: "fract4d"
    }

def copy_libs(dummy,dirpath,namelist):
     for name in namelist:
         target = lib_targets.get(name)
         if target != None:
             shutil.copy(os.path.join(dirpath, name), target)

os.path.walk("build",copy_libs,None)

