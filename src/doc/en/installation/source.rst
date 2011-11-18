.. comment:
   ****************************
   If you alter this document, please change the last line ("This page
   was last updated in ...")
   ****************************

Install from Source Code
========================

More familiarity with computers may be required to build Sage from
the `source code <http://en.wikipedia.org/wiki/Source_code>`_. If you do have all the
pre-requisite tools, the process should
be completely painless. It will take your computer a while to
compile Sage from the source code, although you don't have to watch. Compiling
Sage from the source code has the major advantage that you have the latest
version of Sage with which you can change absolutely any part
or the programs on which Sage depends. You can also recompile Sage.
Also, some parts of Sage will be optimised for your particular computer,
so will run faster than a binary that you have downloaded.

Sage is supported on a number of
`Linux <http://en.wikipedia.org/wiki/Linux>`_
, Mac `OS X <http://www.apple.com/macosx/>`_ , Sun/Oracle `Solaris <http://www.oracle.com/solaris>`_ and
`OpenSolaris <http://en.wikipedia.org/wiki/OpenSolaris>`_
releases, but Sage is not supported on all versions of Linux, OS X,
Solaris or OpenSolaris. Depending on the `operating system <http://en.wikipedia.org/wiki/Operating_system>`_, Sage works
with `x86 <http://en.wikipedia.org/wiki/X86>`_, `x64 <http://en.wikipedia.org/wiki/X86-64>`_, `PowerPC <http://en.wikipedia.org/wiki/PowerPC>`_ or `SPARC <http://en.wikipedia.org/wiki/SPARC>`_ processors. There is no native version of Sage which
installs on `Microsoft Windows <http://en.wikipedia.org/wiki/Microsoft_Windows>`_, although Sage can be used on Windows
with the aid of a  `virtual machine <http://en.wikipedia.org/wiki/Virtual_machine>`_ .
Go to http://www.sagemath.org/download-windows.html
to download a version of Sage for Windows. See http://wiki.sagemath.org/SupportedPlatforms
for the list of platforms on which Sage is supported and the level of support
for these systems. You will also find details about `ports <http://en.wikipedia.org/wiki/Computer_port_%28software%29>`_
to other operating systems or processors which may be taking place.

Assumptions: You have a computer with at least 3 GB of free
disk space running one of the supported version of an
operating system listed at
http://wiki.sagemath.org/SupportedPlatforms
The following standard
command-line development tools must be installed on your computer.
(Under OS X they all come with `Xcode <http://developer.apple.com/xcode/>`_, with the exception of the
`Fortran <http://en.wikipedia.org/wiki/Fortran>`_ compiler ``gfortran``.
However, Sage includes an `executable <http://en.wikipedia.org/wiki/Executable>`_
Fortran compiler for OS X, so there is no need to install
a Fortran compiler on OS X):

::

       gcc        (Version 4.0.1 or later)
       g++        (Version 4.0.1 or later)
       gfortran   (Version 4.0.1 or later)
       make       (For Solaris or OpenSolaris, GNU make, version 3.80 or later)
       perl       (Version 5.8.0 or later)
       ranlib
       tar        (For Solaris or OpenSolaris, GNU tar, version 1.17 or later)
       ssh-keygen (Needed to run the notebook in secure mode)
       latex      (Highly recommended, though not strictly required)
       ImageMagick -- recommended
       ffmpeg -- recommended
       dvipng -- recommended

The programs ``gcc``, ``g++`` and ``gfortran`` are all part of the `GNU Compiler Collection (GCC) <http://gcc.gnu.org/>`_.
You are generally advised to use a recent version of GCC, though
some obscure bugs have stopped specific versions of GCC
compiling Sage on particular platforms.

To check if you have ``perl`` installed, for example, type

::

       command -v perl


on the command line. If it gives an error (or returns nothing), then
either ``perl`` is not installed, or it is installed but not in your
`PATH <http://en.wikipedia.org/wiki/PATH_%28variable%29>`_.
It is highly recommended that you have `Latex <http://en.wikipedia.org/wiki/LaTeX>`_
installed, but it is not required. If you don't have ``ssh-keygen`` on your
local system, then you cannot run the notebook in secure mode, which the uses
encrypted `HTTPS <http://en.wikipedia.org/wiki/HTTP_Secure>`_ protocol. To run the notebook in secure mode, type the command
``notebook(secure=True)`` instead of ``notebook()``. Unless ``notebook(secure=True)``
is used, the notebook uses the less secure `HTTP <http://en.wikipedia.org/wiki/HTTP>`_ protocol

If you don't have either ImageMagick or ffmpeg, you won't be able to
view animations.  ffmpeg can produce animations in more different
formats than ImageMagick, and seems to be faster than ImageMagick when
creating animated GIFs.  Either ImageMagick or dvipng is used for
displaying some LaTeX output in the Sage notebook.

In OS X, make sure you have a recent version of `Xcode <http://developer.apple.com/xcode/>`_.
See http://wiki.sagemath.org/SupportedPlatforms to find out what
version(s) of Xcode are supported. You can get the latest Xcode
from http://developer.apple.com/xcode/, but may have to pay a small
fee in order to download this.

On Linux systems (e.g., Ubuntu, Redhat etc), ``ranlib`` is in the
`Binutils <http://www.gnu.org/software/binutils/>`_ package.
Assuming you have sufficient privileges,
you can install the ``binutils`` and other necessary components. If you
do not have the privileges to do this, ask your system
administrator to do this, or build the components from source
code. The method of installing additional software varies from
distribution to distribution
but on a `Debian <http://www.debian.org/>`_ based system (e.g. `Ubuntu <http://www.ubuntu.com/>`_ or `Mint <http://www.linuxmint.com/>`_), you would use:

::

     sudo apt-get install build-essential gfortran

(this was tested on Ubuntu 9.04).

On other Linux systems you might use `rpm <http://en.wikipedia.org/wiki/RPM_Package_Manager>`_,
`yum <http://en.wikipedia.org/wiki/Yellowdog_Updater,_Modified>`_ or other package manager. On
Solaris you would use ``pkgadd`` and on OpenSolaris use ``ipf``. Check
the documentation for your particular operating system.

The LaTeX package and a PDF previewer are optional but they can be
installed using

::

    sudo apt-get install texlive xpdf evince xdvi

On other systems it might be necessary to install TeX Live from source code,
which is quite easy, though a rather time-consuming process.

On Solaris or OpenSolaris, you must have the GNU version of ``make``
installed and it must be the first
``make`` in your PATH. On Solaris 10, a version of GNU ``make`` may be found
at ``/usr/sfw/bin/gmake`` but you will need to copy it somewhere else
and rename it to ``make``. The same is true for GNU ``tar`` - there is a version
called ``gtar`` in ``/usr/sfw/bin`` but it will need to be copied somewhere
else and renamed to ``tar``). If you attempt to build Sage on AIX or HP-UX,
you will need to install both GNU tar and GNU make. On OS X, the BSD tar
suppied will build Sage, so there is no need to install GNU tar.

For Solaris, it is recommended you create a directory ``$HOME/bins-for-sage`` and
put the GNU versions of ``tar`` and ``make`` in that directory. Then ensure that
``$HOME/bins-for-sage`` is first in your PATH. That's because Sage also needs
``/usr/ccs/bin`` in your PATH to execute programs like ``ar`` and ``ranlib`,
but ``/usr/ccs/bin`` has the Sun/Oracle versions of ``make`` and ``tar``
which are unsuitable for building Sage. For more information on
building Sage on Solaris, see http://wiki.sagemath.org/solaris

Although some of Sage is written in `Python <http://www.python.org/>`_, you do not need Python
pre-installed on your computer, since the Sage installation
includes virtually everything you need. When the Sage installation program is run,
it will check that you have each of the above-listed prerequisites,
and inform you of any that are missing, or have unsuitable verisons.

-  If you want to use `Tcl/Tk <http://www.tcl.tk/>`_ libraries in Sage,
   do the following before compiling Sage.
   Sage's Python will automatically recognize your system's
   install of Tcl/Tk if it exists. You need to install the
   Tcl/Tk development libraries though, not just the Tck/Tk base.

   On `Ubuntu <http://www.ubuntu.com/>`_, this is the command::

       sudo apt-get install tk8.5-dev    # or the latest version available

   Now you can install Sage, If you forgot
   and installed Sage first anyway, all is not lost.
   Just issue the command::

       sage -f  python-2.6.4.p9    # or the latest version available

   after installing the Tcl/Tk development libraries as above.
   If

   .. skip

   ::

       sage: import _tkinter
       sage: import Tkinter

   does not raise an ``ImportError`` then it worked.

-  Sage developers tend to use fairly recent versions of gcc, but
   Sage should compile with almost any gcc of at least version 4.0.1

   If you are interested in working on support for commerical compilers
   from `HP <http://docs.hp.com/en/5966-9844/ch01s03.html>`_,
   `IBM <http://www-01.ibm.com/software/awdtools/xlcpp/>`_,
   `Intel <http://software.intel.com/en-us/articles/intel-compilers/>`_,
   `Sun/Oracle <http://www.oracle.com/technetwork/server-storage/solarisstudio/overview/index.html>`_ etc,
   or the open-source `Clang <http://clang.llvm.org/>`_,
   please email the sage-devel mailing list, otherwise known as the
   sage-devel Google group at
   http://groups.google.com/group/sage-devel

After extracting the Sage tarball, the subdirectory ``spkg`` contains
the source distributions for everything on which Sage depends. We
emphasize that all of this software is included with Sage, so you
do not have to worry about trying to download and install any one
of these packages (such as GAP, for example) yourself.

Fortran
-------

On Linux, Solaris and OpenSolaris systems, a working Fortran compiler is required
for building Sage from source. If you are using Fortran on a platform
for which Sage does not include g95 binaries, you must use
``gfortran``, which may be installed system wide (e.g in /usr or
/usr/local) or your own private copy.  You need to explicitly
tell the Sage build process
about the Fortran compiler and library location. Do this by typing ::

    export SAGE_FORTRAN=/exact/path/to/gfortran
    export SAGE_FORTRAN_LIB=/path/to/fortran/libs/libgfortran.so

Note that the :envvar:`SAGE_FORTRAN` environment variable is supposed to
impact *only* the Fortran Sage package, otherwise known as the Fortran
spkg. Apart from that, this variable is *not* designed to do anything
at all to other spkg's that use Fortran. For example, the Lapack spkg
uses Fortran, but the compilation process of Lapack should ignore the
:envvar:`SAGE_FORTRAN` environment variable. The :envvar:`SAGE_FORTRAN`
environment variable does not mean "build any spkg that uses Fortran
using this Fortran". It means "when installing the Fortran spkg, setup
the ``sage_fortran`` script to run the Fortran compiler specified by
the :envvar:`SAGE_FORTRAN` variable".

On Mac OS X, you are not required to have a Fortran compiler on your
system. The Sage source distribution is shipped with a Fortran
compiler for Mac OS X. This Fortran compiler is used, unless you
specify another Fortran compiler via the variable :envvar:`SAGE_FORTRAN`.

On operating systems such as `AIX <http://en.wikipedia.org/wiki/IBM_AIX>`_,
`HP-UX <http://en.wikipedia.org/wiki/HP-UX>`_, Solaris and OpenSolaris, where both 32-bit and
64-bit builds are supported, the library path variable
:envvar:`SAGE_FORTRAN_LIB` must point to the 32-bit library if you are
building Sage in 32-bit. Also, :envvar:`SAGE_FORTRAN_LIB` must point to a
64-bit library if you are building Sage in 64-bit. For example, on
Solaris & OpenSolaris, the variables :envvar:`SAGE_FORTRAN`,
:envvar:`SAGE_FORTRAN_LIB` and :envvar:`SAGE64` could be set as follows::

    # SPARC, x86 and x64.
    SAGE_FORTRAN=/path/to/gcc/install/directory/bin/gfortran

    # 32-bit SPARC
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/libgfortran.so

    # 64-bit SPARC
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/sparcv9/libgfortran.so
    SAGE64=yes

    # 32-bit x86
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/libgfortran.so

    # 64-bit x64
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/amd64/libgfortran.so
    SAGE64=yes

(It should be noted that Sage is not supported on AIX or HP-UX, although some
efforts have been made to `port Sage to AIX <http://wiki.sagemath.org/AIX>`_ and
to `port Sage to HP-UX <http://wiki.sagemath.org/HP-UX>`_.)

Steps to Install from Source
----------------------------

Installation from source is (potentially) very easy, because the
distribution contains (essentially) everything on which Sage
depends.

Make sure there are **no spaces** in the path name for the directory
in which you build: several of Sage's components will not build if
there are spaces in the path.  Running Sage from a directory with
spaces in its name will also fail.

#. Go to http://www.sagemath.org/download-source.html, select a mirror,
   and download the file ``sage-x.y.z.tar``.

   This tarfile contains the source code for Sage and the source for
   all programs on which Sage depends. Download it into a subdirectory
   of your home directory into which you want to install Sage. Note
   that this file is not compressed; it's just a plain tarball (which
   happens to be full of compressed files).

#. Extract:

   ::

             tar xvf sage-x.y.z.tar

#. This creates a directory ``sage-x.y.z``.

#. Change into that directory

   ::

             cd sage-x.y.z

   This is Sage's home directory. It is also referred to as
   ``SAGE_ROOT`` or the top level Sage directory.

#. Optional (but highly recommended): Read the ``README.txt`` file
   there.

#. On OSX 10.4, OS 10.5, Solaris 10 and OpenSolaris, if you wish to
   build a 64-bit version of Sage, then assuming your computer and
   operating system are 64-bit, type

   ::

           SAGE64=yes
           export SAGE64

   It should be noted that at the time of writing (April 2011), 64-bit
   builds of Sage on both Solaris 10 and OpenSolaris are not very stable,
   so you are advised not to set ``SAGE64`` to ``yes``. This will then
   create stable 32-bit versions of Sage.
   See http://wiki.sagemath.org/SupportedPlatforms  and
   http://wiki.sagemath.org/solaris for the latest information, as
   work is ongoing to resolve the 64-bit Solaris & OpenSolaris problems.

#. Type

   ::

             make

   This compiles Sage and all dependencies. Note that you do not need
   to be logged in as root, since no files are changed outside of the
   ``sage-x.y.z`` directory (with one exception -- the ``.ipythonrc``
   directory is created in your ``HOME`` directory if it doesn't exist).
   In fact, **it is inadvisable to build Sage as root**, as the root account
   should only be used when absolutely necessary, as mis-typed commands
   can have serious consequences if you are logged in as root.  There has been a bug
   `reported <http://trac.sagemath.org/sage_trac/ticket/9551/>`_ in Sage
   which would have overwritten a system file had the user been logged in
   as root.

   Typing ``make`` does the usual steps for each of the packages, but puts
   all the results in the local build tree. Depending on the architecture of your system (e.g.,
   Celeron, Pentium Mobile, Pentium 4, SPARC, etc.), it can take over three hours
   to build Sage from source. On slower older hardware it can take over
   a day to build Sage. If the build is successful, you will not see
   the word ERROR in the last 3-4 lines of output.

   Each component of Sage has its own build log, saved in
   ``SAGE_ROOT/spkg/logs``.  In particular,
   if the build of Sage fails, then you can type the following from the directory
   where you typed ``make``.

   ::

            grep -li "^Error installing" spkg/logs/*

   Then paste the contents of the log file(s) with errors to the Sage
   support newsgroup http://groups.google.com/group/sage-support
   If the log files are very large (and many are), then don't paste
   the whole file, but make sure to include any error messages.

   The directory where you built Sage is NOT hardcoded. You should
   be able to safely move or rename that directory. (It's a bug if
   this is not the case)

#. To start Sage, change into the Sage home directory and type:

   ::

             ./sage

   You should see the Sage prompt, which will look something like this
   (starting the first time should take well under a minute, but can
   take several minutes if the file system is slow or busy. Since Sage
   opens a lot of files, it is preferable to install Sage on a fast file
   system if this is possible.):

   ::

       $ sage
       ----------------------------------------------------------------------
       | Sage Version 4.7, Release Date: 2011-05-23                         |
       | Type notebook() for the GUI, and license() for information.        |
       ----------------------------------------------------------------------
       sage:

   Just starting successfully tests that many of the components built
   correctly. If the above is not displayed (e.g., if you get a
   massive traceback), please report the problem, e.g., to
   http://groups.google.com/group/sage-support .
   It would also be helpful to
   include the type of operating system (Linux, OS X, Solaris or OpenSolaris),
   the version and date of that operating system and the version
   number of the copy of Sage you are using. (There are no
   formal requirements for bug reports - just send them; we appreciate
   everything.)

   After Sage starts, try a command:

   ::

       sage: 2 + 2
       4

   Try something more complicated, which uses the PARI C library:

   ::

       sage: factor(2005)
       5 * 401

   Try something simple that uses the Gap, Singular, Maxima and
   PARI/GP interfaces:

   ::

       sage: gap('2+2')
       4
       sage: gp('2+2')
       4
       sage: maxima('2+2')
       4
       sage: singular('2+2')
       4
       sage: pari('2+2')
       4

   (For those familiar with GAP: Sage automatically builds a GAP
   "workspace" during installation, so the response time from this GAP
   command is relatively fast. For those familiar with GP/PARI, the
   ``gp`` command creates an object in the GP interpreter, and the
   ``pari`` command creates an object directly in the PARI C-library.)

   Try running Gap, Singular or GP from Sage:

   .. skip

   ::

       sage: gap_console()
       GAP4, Version: 4.4.12 of 17-Dec-2008, i386-pc-solaris2.11-gcc
       gap> 2+2;
       4
       [ctrl-d]

   .. skip

   ::

       sage: gp_console()
       ...
       [ctrl-d]

   .. skip

   ::

       sage: singular_console()
                            SINGULAR                             /  Development
        A Computer Algebra System for Polynomial Computations   /   version 3-1-1
                                                              0<
            by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   Feb 2010
       FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
       [ctrl-d]
       > Auf Wiedersehen.
       sage:

#. Optional: Check the interfaces to any other software that
   you have available. Note that each interface calls its
   corresponding program by a particular name:
   `Mathematica <http://www.wolfram.com/mathematica/>`_ is invoked
   by calling ``math``, `Maple <http://www.maplesoft.com/>`_ by calling ``maple``, etc. The
   easiest way to change this name or perform other customizations is
   to create a redirection script in ``$SAGE_ROOT/local/bin``. Sage
   inserts this directory at the front of your PATH, so your script
   may need to use an absolute path to avoid calling itself; also,
   your script should use ``$*`` to pass along all of its arguments.
   For example, a ``maple`` script might look like:

   ::

       #!/bin/sh

       /etc/maple10.2/maple.tty $*

#. Optional: Different possibilities to make using Sage a little
   easier:

   - Make a symbolic link from ``/usr/local/bin/sage`` (or another
     directory in your :envvar:`PATH`) to ``$SAGE_ROOT/sage``::

         ln -s /path/to/sage-x.y.z/sage /usr/local/bin/sage

     Now simply typing ``sage`` should be sufficient to run Sage.

   - Copy ``$SAGE_ROOT/sage`` to a location in your ``PATH``. If you do
     this, make sure you edit the line ``#SAGE_ROOT=/path/to/sage-version``
     at the top of the copied ``sage`` script. It is best to edit only
     the copy, not the original.

   -  For KDE users, create a bash script {sage} containing the lines

      ::

          #!/bin/bash
          konsole -T "sage" -e <SAGE_ROOT>/sage

      which you make executable (``chmod a+x sage``) and put it somewhere in
      your path. (Note that you have to change ``$SAGE_ROOT`` above!) You
      can also make a KDE desktop icon with this as the command (under
      the Application tab of the Properties of the icon, which you get my
      right clicking the mouse on the icon).

   - On Linux and OS X systems, you can make an alias to ``$SAGE_ROOT/sage``.
     For example, put something similar to the following line in your
     ``.bashrc`` file::

         alias 'sage'='/home/username/sage-4.8/sage'

     Having done so, quit your terminal emulator and restart it again.
     Now typing ``sage`` within your terminal emulator should start
     Sage.

#. Optional, but highly recommended: Test the install by typing ``./sage -testall``. This
   runs most examples in the source code and makes sure that they run
   exactly as claimed. To test all examples, use
   ``./sage -testall -optional -long``; this will run examples that take
   a long time, and those that depend on optional packages and
   software, e.g., Mathematica or Magma. Some (optional) examples will
   likely fail because they assume that a database is installed.
   Alternatively, from within ``$SAGE_ROOT``, you can type
   ``make test`` to run all the standard test code.  This can take
   from 25 minutes to several hours, depending on your hardware. On
   very old hardware building and testing Sage can take several days!

#. Optional: Install optional Sage packages and databases. Type
   ``sage -optional`` to see a list or visit
   http://www.sagemath.org/packages/optional/, and
   ``sage -i <package name>`` to automatically download and install a
   given package.

#. Optional: Run the ``install_scripts`` command from within Sage to create
   gp, singular, gap, etc., scripts in your ``PATH``. Type
   ``install_scripts?`` in Sage for details.


Have fun! Discover some amazing conjectures!

Environment variables
---------------------

Sage uses several environment variables to control its build process.
Most users won't need to set any of these: the build process just
works on many platforms.  (Note though that setting :envvar:`MAKE`, as
described below, can significantly speed up the process.)  Building
Sage involves building about 100 packages, each of which has its own
compilation instructions.

Here are some of the more commonly used variables affecting the build
process:

- :envvar:`MAKE` - one useful setting for this variable when building
  Sage is ``MAKE='make -jNUM'`` to tell the "make" program to
  run NUM jobs in parallel when building.  Some people advise using
  more jobs than there are CPU cores, at least if the system is not
  heavily loaded and has plenty of RAM; for example, a good setting
  for NUM might be between 1 and 1.5 times the number of cores.  In
  addition, the "-l" option sets a load limit: ``MAKE='make -j4
  -l5.5``, for example, tells "make" to try to use four jobs, but to
  not start more than one job if the system load average is above 5.5.
  See the manual page for GNU make: `Command-line options
  <http://www.gnu.org/software/make/manual/make.html#Options-Summary>`_
  and `Parallel building
  <http://www.gnu.org/software/make/manual/make.html#Parallel>`_.

  .. warning::

     Some users on single-core OS X machines have reported problems
     when building Sage with ``MAKE='make -jNUM'`` with NUM greater
     than one.

- :envvar:`SAGE_PARALLEL_SPKG_BUILD` - if this is set to "no", then
  build spkgs serially rather than in parallel.  If this is "no", then
  each spkg may still take advantage of the setting of :envvar:`MAKE`
  to build using multiple jobs, but the spkgs will be built one at a
  time.  Alternatively, run "make build-serial" which sets this
  environment variable for you.

- :envvar:`SAGE_CHECK` - if this is set to "yes", then during the
  build process, run the test suite for each package which has one.

- :envvar:`SAGE64` - Set this to "yes" to build a 64-bit binary on platforms
  which default to 32-bit, even though they can build 64-bit binaries.
  It adds the compiler flag
  -m64 when compiling programs.  The SAGE64 variable is mainly of use
  is on OS X (pre 10.6), Solaris and OpenSolaris, though it will add
  the -m64 on any operating system. If you are running version 10.6 of
  OS X on a 64-bit machine, then Sage will automatically build a
  64-bit binary, so this variable does not need setting.

- :envvar:`CFLAG64` - default value "-m64".  If Sage detects that it
  should build a 64-bit binary, then it uses this flag when compiling
  C code.  Modify it if necessary for your system and C compiler.
  This should not be necessary on most systems -- this flag will
  typically be set automatically, based on the setting of
  :envvar:`SAGE64`, for example.

- :envvar:`SAGE_FORTRAN` - see above, the "Fortran" section.

- :envvar:`SAGE_FORTRAN_LIB` - see above, the "Fortran" section.

- :envvar:`SAGE_DEBUG` - about half a dozen Sage packages use this
  variable.  If it is unset (the default) or set to "yes", then
  debugging is turned on.  If it is set to anything else, then
  debugging is turned off.

- :envvar:`SAGE_SPKG_LIST_FILES` - Set this to "yes" to enable
  verbose extraction of tar files, i.e. Sage's spkg files. Since
  some spkgs contain a huge number of files such that the log files
  get very large and harder to search (and listing the contained
  files is usually less valuable), we decided to turn this off
  by default. This variable affects builds of Sage with ``make``
  (and ``sage -upgrade``) as well as the manual installation of
  individual spkgs with e.g. ``sage -i``.

- :envvar:`SAGE_SPKG_INSTALL_DOCS` - Set this to "yes" to install
  package-specific documentation to
  :file:`$SAGE_ROOT/local/share/doc/PACKAGE_NAME/` when an spkg is
  installed.  This option may not be supported by all spkgs.  Some
  spkgs might also assume that certain programs are available on the
  system (for example, ``latex`` or ``pdflatex``).

- :envvar:`SAGE_FAT_BINARY` - to prepare a binary distribution that
  will run on the widest range of target machines, set this variable
  to "yes" before building Sage::

      export SAGE_FAT_BINARY="yes"
      make
      ./sage -bdist x.y.z-fat

Variables to set if you're trying to build Sage with an unusual setup,
e.g., an unsupported machine or an unusual compiler:

- :envvar:`SAGE_PORT` - if you try to build Sage on a platform which
  is recognized as being unsupported (e.g. AIX, or
  HP-UX), or with a compiler which is unsupported (anything except
  gcc), you will see a message saying something like ::

        You are attempting to build Sage on IBM's AIX operating system,
        which is not a supported platform for Sage yet. Things may or
        may not work. If you would like to help port Sage to AIX,
        please join the sage-devel discussion list - see
        http://groups.google.com/group/sage-devel
        The Sage community would also appreciate any patches you submit.

        To get past this message, export the variable SAGE_PORT to
        something non-empty.

  If this is the situation, follow the directions: set
  :envvar:`SAGE_PORT` to something non-empty (and expect to run into
  problems).

 :envvar:`SAGE_USE_OLD_GCC` - the Sage build process requires
  gcc with a version number of at least 4.0.1.
  If the most recent version of gcc on your system is the older 3.4.x series and you
  want to try building anyway, then set :envvar:`SAGE_USE_OLD_GCC` to
  something non-empty. Expect the build to fail in this case.

Environment variables dealing with specific Sage packages:

- :envvar:`SAGE_ATLAS_ARCH` - if you are compiling ATLAS (in
  particular, if :envvar:`SAGE_ATLAS_LIB` is not set), you can use
  this environment variable to set a particular architecture and
  instruction set architecture. The syntax is
  ``SAGE_ATLAS_ARCH=arch[,isaext1][,isaext2]...[,isaextN]``. While
  ATLAS comes with precomputed timings for a variety of CPUs, it only
  uses them if it finds an exact match. Otherwise, ATLAS runs through
  a lengthy automated tuning process in order to optimize performance
  for your particular system. You drastically reduce the total Sage
  compile time if you manually select a suitable architecture. It is
  recommended to specify a suitable architecture on laptops or other
  systems with CPU throttling or if you want to distribute the
  binaries. Available architectures are

    ``POWER3``, ``POWER4``, ``POWER5``, ``PPCG4``, ``PPCG5``, ``P5``,
    ``P5MMX``, ``PPRO``, ``PII``, ``PIII``, ``PM``, ``CoreSolo``,
    ``CoreDuo``, ``Core2Solo``, ``Core2``, ``Corei7``, ``P4``,
    ``P4E``, ``Efficeon``, ``K7``, ``HAMMER``, ``AMD64K10h``,
    ``IA64Itan``, ``IA64Itan2``, ``USI``, ``USII``, ``USIII``,
    ``USIV``, ``UnknownUS``, ``MIPSR1xK``, ``MIPSICE9``

  and instruction set extensions are

    ``AltiVec``, ``SSE3``, ``SSE2``, ``SSE1``, ``3DNow``.

  In addition, you can also set

  - ``SAGE_ATLAS_ARCH=fast`` picks defaults for a modern (2-3 year old)
    CPU of your processor line, and

  - ``SAGE_ATLAS_ARCH=base`` picks defaults that should work for a ~10
    year old CPU.

  For example,

    ``SAGE_ATLAS_ARCH=Corei7,SSE3,SSE2,SSE1``

  would be appropriate for a Core i7 CPU.

- :envvar:`SAGE_ATLAS_LIB` - if you have an installation of ATLAS on
  your system and you want Sage to use it instead of building and
  installing its own version of ATLAS, set this variable to be the
  directory containing your ATLAS installation. It should contain the
  files :file:`libatlas`, :file:`liblapack`, :file:`libcblas`, and
  :file:`libf77blas` with extensions ``.a``, ``.so``, or
  ``.dylib``. For backward compatibility, the libraries may also be in
  the subdirectory ``SAGE_ATLAS_LIB/lib/``.

- :envvar:`SAGE_MATPLOTLIB_GUI` - set this to anything non-empty except
  "no", and Sage will attempt to build the graphical backend when it
  builds the matplotlib package.

- :envvar:`INCLUDE_MPFR_PATCH` - This is used to add a patch to MPFR
  to bypass a bug in the memset function affecting sun4v machines with
  versions of Solaris earlier than Solaris 10 update 8
  (10/09). Earlier versions of Solaris 10 can be patched by applying
  Sun patch 142542-01.  Recognized values are:

  - ``INCLUDE_MPFR_PATCH=0`` - never include the patch - useful if you
    know all sun4v machines Sage will be used are running Solaris
    10 update 8 or later, or have been patched with Sun patch
    142542-01.

  - ``INCLUDE_MPFR_PATCH=1`` - always include the patch, so the binary
    will work on a sun4v machine, even if created on an older sun4u
    machine.

  If this variable is unset, include the patch on sun4v machines only.

- :envvar:`SAGE_BINARY_BUILD` - used by the pil package.  If set to
  "yes", then force Sage to use the versions of libjpeg, libtiff and
  libpng from :file:`$SAGE_ROOT/local/lib`.  Otherwise, allow the use
  of the system's versions of these libraries.

- :envvar:`SAGE_PIL_NOTK` - used by the pil package.  If set to "yes",
  then disable building TK.  If this is not set, then this should be
  dealt with automatically: Sage tries to build the pil package with
  TK support enabled, but if it runs into problems, it tries building
  again with TK disabled.  So only use this variable to force TK to be
  disabled.  (Building the pil package is pretty fast -- less than a
  minute on many systems -- so allowing it to build twice is not a
  serious issue.)

Some standard environment variables which you should probably **not**
set:

- :envvar:`CC` - while some programs allow you to use this to specify
  your C compiler, the Sage packages do **not** all recognize this.
  In fact, setting this variable for building Sage is likely to cause
  the build process to fail.

- :envvar:`CXX` - similarly, this will set the C++ complier for some
  Sage packages, and similarly, using it is likely quite risky.

- :envvar:`CFLAGS`, :envvar:`CXXFLAGS` - the flags for the C compiler
  and the C++ compiler, respectively.  The same comments apply to
  these: setting them may cause problems, because they are not
  universally respected among the Sage packages.

Sage uses the following environment variables when it runs:

- :envvar:`DOT_SAGE` - this is the directory, to which the user has
  read and write access, where Sage stores a number of files.  The
  default location is ``~/.sage/``, but you can change that by setting
  this variable.

- :envvar:`SAGE_STARTUP_FILE` - a file including commands to be
  executed every time Sage starts.  The default value is
  ``$DOT_SAGE/init.sage``.

- :envvar:`SAGE_SERVER` - if you want to install a Sage package using
  ``sage -i PKG_NAME``, Sage downloads the file from the web, using
  the address ``http://www.sagemath.org/`` by default, or the address
  given by :envvar:`SAGE_SERVER` if it is set.  If you wish to set up
  your own server, then note that Sage will search the directories
  ``SAGE_SERVER/packages/standard/``,
  ``SAGE_SERVER/packages/optional/``,
  ``SAGE_SERVER/packages/experimental/``, and
  ``SAGE_SERVER/packages/archive/`` for packages.  See the script
  :file:`$SAGE_ROOT/local/bin/sage-download_package` for the
  implementation.

- :envvar:`SAGE_PATH` - a colon-separated list of directories which
  Sage searches when trying to locate Python libraries.

- :envvar:`SAGE_BROWSER` - on most platforms, Sage will detect the
  command to run a web browser, but if this doesn't seem to work on
  your machine, set this variable to the appropriate command.

- :envvar:`SAGE_ORIG_LD_LIBRARY_PATH_SET` - set this to something
  non-empty to force Sage to set the :envvar:`LD_LIBRARY_PATH` before
  executing system commands.

- :envvar:`SAGE_ORIG_DYLD_LIBRARY_PATH_SET` - similar, but only used
  on Mac OS X to set the :envvar:`DYLD_LIBRARY_PATH`.

- :envvar:`SAGE_CBLAS` - used in the file
  :file:`SAGE_ROOT/devel/sage/sage/misc/cython.py`.  Set this to the
  base name of the BLAS library file on your system if you want to
  override the default setting.  That is, if the relevant file is
  called :file:`libcblas_new.so` or :file:`libcblas_new.dylib`, then
  set this to "cblas_new".

Sage overrides the user's settings of the following variables:

- :envvar:`MPLCONFIGDIR` - ordinarily, this variable lets the user set
  their matplotlib config directory.  Due to incompatibilies in the
  contents of this directory among different versions of matplotlib,
  Sage overrides the user's setting, defining it instead to be
  ``$DOT_SAGE/matplotlib-VER``,   with "VER" replaced by the
  current matplotlib version number.

Variables dealing with doctesting:

- :envvar:`SAGE_TESTDIR` - a temporary directory used during Sage's
  doctesting.  The default is to use the directory ``$DOT_SAGE/tmp``,
  but you can override that by setting this variable.

- :envvar:`SAGE_TIMEOUT` - used for Sage's doctesting: the number of
  seconds to allow a doctest before timing it out.  If this isn't set,
  the default is 360 seconds (6 minutes).

- :envvar:`SAGE_TIMEOUT_LONG` - used for Sage's doctesting: the number
  of seconds to allow a doctest before timing it out, if tests are run
  using ``sage -t --long``.  If this isn't set, the default is 1800
  seconds (30 minutes).

- :envvar:`SAGE_PICKLE_JAR` - if you want to update the the standard
  pickle jar, set this to something non-empty and run the doctest
  suite.  See the documentation for the functions :func:`picklejar`
  and :func:`unpickle_all` in
  :file:`SAGE_ROOT/devel/sage/sage/structure/sage_object.pyx`, online
  `here (picklejar)
  <http://sagemath.org/doc/reference/sage/structure/sage_object.html#sage.structure.sage_object.picklejar>`_
  and `here (unpickle_all)
  <http://sagemath.org/doc/reference/sage/structure/sage_object.html#sage.structure.sage_object.unpickle_all>`_.

..
  THIS INDENTED BLOCK IS A COMMENT.  FIX IT ONCE WE UNDERSTAND
  THESE VARIABLES.

  Variables dealing with valgrind and friends:

  - :envvar:`SAGE_TIMEOUT_VALGRIND` - used for Sage's doctesting: the
    number of seconds to allow a doctest before timing it out, if tests
    are run using ``??``.  If this isn't set, the default is 1024*1024
    seconds.

  - :envvar:`SAGE_VALGRIND` - ?

  - :envvar:`SAGE_MEMCHECK_FLAGS`, :envvar:`SAGE_MASSIF_FLAGS`,
    :envvar:`SAGE_CACHEGRIND_FLAGS`, :envvar:`SAGE_OMEGA_FLAGS` - flags
    used when using valgrind and one of the tools "memcheck", "massif",
    "cachegrind", or "omega"

Installation in a Multiuser Environment
---------------------------------------

This section addresses the question of how a system administrator
can install a single copy of Sage in a multi-user computer
network.

System-wide install
~~~~~~~~~~~~~~~~~~~

#. After you build Sage, you may optionally copy or move the entire
   build tree to ``/usr/local`` or another location.  If you do this,
   then you must run ``./sage`` once so that various hard-coded
   locations will get updated.  For this reason, it might be easier to
   simply build Sage in its final location.

#. Make a symbolic link to the ``sage`` script in ``/usr/local/bin``::

       ln -s /path/to/sage-x.y.z/sage /usr/local/bin/sage

   Alternatively, copy the Sage script::

       cp /path/to/sage-x.y.z/sage /usr/local/bin/sage

   and edit the file ``/usr/local/bin/sage``: ``SAGE_ROOT`` should be
   set to the directory ``/path/to/sage-x.y.z/`` where Sage is
   installed.  It is recommended not to edit the original ``sage``
   script, only the copy in ``/usr/local/bin/sage``.

#. Make sure that all files in the Sage tree are readable by all::

       chmod a+rX -R /usr/local/sage-4.8

#. Optionally, you can test Sage by running::

       make testlong

   or ``make ptestlong`` which tests files in parallel using multiple
   processes. You can also omit ``long`` to skip tests which take a long
   time.

Some common problems
--------------------

ATLAS
~~~~~

Sometimes the ATLAS spkg can fail to build.  Some things to check for:

- Make sure that CPU throttling mode (= power-saving mode) is turned off
  when building ATLAS.

- Also, the ATLAS build can fail if the system load is too high, and in
  particular this has been known to happen when building with
  ``MAKE='make -jNUM'`` with NUM large.  If this happens, just try
  running "make" again.  If "make" fails after five or six attempts,
  report your problem to the sage-devel mailing list.

Special Notes
-------------


-  (Found by Dorian Raymer) Sage will not build if you have only
   bison++. You should uninstall bison++ and install `bison <http://www.gnu.org/software/bison/>`_.

-  (Found by Peter Jipsen) If you get an error like

   ::

       ImportError: /home/jipsen/Desktop/sage-1.3.3.1/local/lib/libpari-gmp.so.2:
            cannot restore segment prot after reloc:
       Permission denied

   then your `SELinux <http://fedoraproject.org/wiki/SELinux>`_ configuration is preventing Sage from launching. To
   rectify this issue, you can either change the default security
   context for Sage (??) or disable SELinux altogether by setting the
   line ``SELINUX=disabled`` in your ``/etc/sysconfig/selinux`` file.

- To make SageTeX available to your users, see the instructions for
  :ref:`installation in a multiuser environment
  <sagetex_installation_multiuser>`.

  **This page was last updated in November 2011**
