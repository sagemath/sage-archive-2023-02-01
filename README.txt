    Sage: Open Source Mathematical Software

       "Creating a Viable Open Source Alternative to
          Magma, Maple, Mathematica, and MATLAB"

    Copyright (C) 2005-2011 William Stein and the Sage Development Team

        http://www.sagemath.org

    Over 200 people have contributed code to Sage. Please see the
    following web page for a list:

        http://www.sagemath.org/development-map.html

    In many cases, documentation for modules and functions list the
    authors.


GETTING STARTED
---------------
This README.txt contains build instructions for Sage. If you downloaded
a binary, you do not need to do anything; just execute:

    ./sage

from the command line. If you downloaded the sources, please read below
on how to build Sage and work around common issues.

If you have questions or encounter problems, please do not hesitate
to email the sage-support mailing list:

    http://groups.google.com/group/sage-support


SUPPORTED PLATFORMS
-------------------

Sage fully supports several Linux distributions, recent versions of
Mac OS X, as well as a number of Solaris and OpenSolaris releases.

There should be no serious bugs in an officially released version of
Sage on any of the fully supported platforms, but any major issues with
a particular release will be documented on an errata page:

    http://wiki.sagemath.org/errata

Ports are in progress to some other, less common platforms. The list of
supported platforms and their current statuses are given at the
following web page:

    http://wiki.sagemath.org/SupportedPlatforms

If you are interested in helping port Sage to a new platform, please let
us know at the sage-devel mailing list:

    http://groups.google.com/group/sage-devel


QUICK INSTRUCTIONS TO BUILD FROM SOURCE
---------------------------------------

The following steps briefly outline the process of building Sage from
source. More detailed instructions, including how to build faster on
multicore machines are contained later in this README and in the
Installation Guide:

    http://www.sagemath.org/doc/installation

1. Make sure you have the dependencies and 2.5 GB of free disk space.

   Linux: GCC, g++, gfortran, make, m4, perl, ranlib, and tar.
   (install these using your package manager)

   OS X: XCode.

       WARNING: If "gcc -v" outputs 4.0.0, you *must* upgrade XCode
       (free from Apple). That version of GCC is very broken.

   Other platforms: See detailed instructions below.

2. Extract the tarball:

       tar -xvf sage-*.tar

3. cd into the Sage directory and type make:

       cd sage-*/
       make

   That's it! Everything is automatic and non-interactive. The build
   should work fine on all fully supported platforms. If it does not, we
   want to know!


FORTRAN
-------

To build Sage on any platform except OS X, you must use a system-wide
gfortran compiler. Sometimes you need to explicitly tell the Sage build
process about the Fortran compiler and library locations. Do this by
typing:

    export SAGE_FORTRAN=/exact/path/to/gfortran
    export SAGE_FORTRAN_LIB=/path/to/fortran/libs/libgfortran.so

Note that the SAGE_FORTRAN environment variable is supposed to impact
*only* the Fortran Sage package, otherwise known as the Fortran spkg
("spkg" stands for "Sage package").  Apart from that, this variable is
*not* designed to do anything at all to other spkg's that use
Fortran. For example, the Lapack spkg uses Fortran, but the
compilation process of Lapack should ignore the SAGE_FORTRAN
environment variable. The SAGE_FORTRAN environment variable does not
mean "build any spkg that uses Fortran using this Fortran". It means
"when installing the Fortran spkg, setup the sage_fortran script to
run the Fortran specified by the SAGE_FORTRAN variable".


SELINUX
--------

On Linux, if you get this error message:

    " restore segment prot after reloc: Permission denied "

the problem is probably related to SELinux. See the following URL for
further information:

    http://www.ittvis.com/services/techtip.asp?ttid=3092


IMPLEMENTATION
--------------

Sage has significant components written in the following languages:
C/C++, Python, Cython, Lisp, and Fortran. Lisp (ECL), Python, and Cython
are built as part of Sage and a GNU Fortran (gfortran) binary is
included (OS X only), so you do not need them in order to build Sage.


MORE DETAILED INSTRUCTIONS TO BUILD FROM SOURCE
-----------------------------------------------

1. Make sure you have about 2.5 GB of free disk space.

2. Install build dependencies.

   Linux: See quick instructions above.

   OS X: Make sure you have XCode version >= 2.4, i.e. "gcc -v" should
   output build >= 5363. If you don't, go to:

       http://developer.apple.com/

   sign up, and download the free XCode package. Only OS X >= 10.4 is
   supported.

   Solaris and OpenSolaris: Building Sage on these platforms is more
   tricky than on Linux or OS X. For details on how to build Sage on
   these platforms, see:

       http://wiki.sagemath.org/solaris

   Windows: Not supported. A solution is to download and install
   VirtualBox, install Linux into it, etc.

   NOTE: On some operating systems, it might be necessary to install
   gas/as, gld/ld, gnm/nm. On most platforms, these are automatically
   installed when you install the programs listed above.

3. Extract the Sage source tarball and cd into a directory with no
   spaces in it. If you have a machine with 4 processors, say, type the
   following to configure the build script to perform a parallel
   compilation of Sage using 4 jobs:

       export MAKE="make -j4"

   By default the above will only parallelize a small part of the build
   process. To make the spkg build part of the installation process
   parallel do:

       export SAGE_PARALLEL_SPKG_BUILD="yes"

   The amount of processors used for the building of the spkg's is
   determined by what you assigned to MAKE. So the above is useless if
   you did export MAKE="make -j1".

   If you want to run the test suite for each individual spkg as it is
   installed, type:

       export SAGE_CHECK="yes"

   before starting the Sage build. This will run each test suite and
   will raise an error if any failures occur.  Warning: on many
   platforms, this will cause failures in the installation of the
   Python spkg.

   To start the build, type:

       make


4. Wait about 1 hour to 14 days, depending on your computer (it took
   about 2 weeks to build Sage on the T-Mobile G1 Android cell phone).

5. Type "./sage" to try it out.

6. OPTIONAL: Start Sage and run the command

       install_scripts("/usr/local/bin/")   # change /usr/local/bin/

   Type "install_scripts?" in Sage for more details about what this
   command does.

7. OPTIONAL: Type "make test" to test all examples in the documentation
   (over 93,000 lines of input!) -- this takes from 30 minutes to
   several hours. Don't get too disturbed if there are 2 to 3 failures,
   but always feel free to email the section of test.log that contains
   errors to the sage-support mailing list. If there are numerous
   failures, there was a serious problem with your build.

8. OPTIONAL: If you want to (try to) build the documentation, run:

       sage -docbuild --help

   for instructions. The HTML version of the documentation is built
   during the compilation process of Sage and resides in the directory:

       $SAGE_ROOT/devel/sage/doc/output/html/

   LaTeX is required to build the PDF version of the documentation.

9. OPTIONAL: It is highly recommended that you install the optional GAP
   database by typing:

       ./sage -optional

   then installing (with "./sage -i") the package whose name begins with
   database_gap. This will download the package from
   sage.math.washington.edu and install it. While you're at it, you
   might install other databases of interest to you.

10. OPTIONAL: It is recommended that you have both LaTeX and the
    ImageMagick tools (e.g. the "convert" command) installed since some
    plotting functionality benefits from it.


SUPPORTED COMPILERS
-------------------

 * Sage needs a version of GCC that is at least version 4.0.1 that is
   configured with support for at least the C, C++, and Fortran
   languages, although a Fortran compiler is not needed on OS X, as
   noted above.

 * The versions of the C compiler (gcc), C++ compiler (g++), and Fortran
   compiler (gfortran) must all be identical.

 * Sage has never been built without using the GCC compiler.


RELOCATION
----------

You *should* be able to move the sage-x.y.z/ directory anywhere you
want. If you copy the sage script or make a symbolic link to it, you
should modify the script to reflect this (as instructed at the top of
the script). It is best if the path to Sage does not have any spaces in
it.

For a system-wide installation, as root you can move the sage-x.y.z/
directory to a system-wide directory. Afterwards, you need to start up
Sage as root at least once prior to using the system-wide Sage as a
normal user. See the Installation Guide for further information on
performing a system-wide installation:

    http://www.sagemath.org/doc/installation/source.html#installation-in-a-multiuser-environment

If you find anything that doesn't work correctly after you moved the
directory, please email the sage-support mailing list.


REDISTRIBUTION
--------------

Your local Sage install is almost exactly the same as any "developer"
install. You can make changes to documentation, source, etc., and very
easily package the complete results up for redistribution just like we
do.

1. To make your own source tarball (sage-x.y.z.tar) of Sage, type:

       sage -sdist x.y.z

   where the version is whatever you want.

2. To make a binary distribution with your currently installed packages,
   type:

       sage -bdist x.y.z

3. To make a binary that will run on the widest range of target
   machines, set the SAGE_FAT_BINARY environment variable to "yes"
   before building Sage:

       export SAGE_FAT_BINARY="yes"
       make
       ./sage -bdist x.y.z-fat

In all cases, the result is placed in the directory "$SAGE_ROOT/dist/".


CHANGES TO INCLUDED SOFTWARE
----------------------------

All software included with Sage is copyrighted by the respective authors
and released under an open source license that is "GPL version 3 or
later" compatible. See the file COPYING.txt for more details.

Each spkg in $SAGE_ROOT/spkg/standard/ is a bzip'd tarball. You can
extract it with:

    tar -jxvf name-*.spkg

Inside the spkg, there is a file SPKG.txt that details all changes made
to the given package for inclusion with Sage. The inclusion of such a
file detailing changes is specifically required by some of the packages
included with Sage (e.g. for GAP).
