
Install from Source Code
========================

More familiarity with computers may be required to build Sage from
source. If you do have all the pre-requisite tools, the process should
be completely painless. It would take your computer a while to
compiled Sage from source, though you don't have to watch. Compiling
Sage from source has the major advantage that you have the latest
version of Sage with which you can change absolutely any part
or the programs on which Sage depends. You can also recompile Sage.

As of this writing, Sage is known to work on Linux (32-bit x86, 64-bit
x86-64, IA64, or 32-bit PPC) and OS X (10.4, 10.5, 10.6, PPC or
x86, 32-bit only). (See http://wiki.sagemath.org/SupportedPlatforms
for the latest information.)

    **Solaris? FreeBSD? OS X 10.5 in 64 bit mode?**: Complete compilation
    of Sage is currently not supported on Solaris or \*BSD. It is
    possible to compile most of Sage on Solaris machines and to fill in
    the extra parts using standard packages; please email sage-devel if
    you desperately need to run Sage on Solaris. We do plan to fully
    support Solaris - it's a very important platform. Work is ongoing.

    **Sage on FreeBSD**: The binaries can be run with the help of Linux
    emulation. We are working on a fully native port and the number of
    issues that need to be fixed are relatively small compared to the
    other ports.

    We hope to support OS X 10.5 in 64-bit mode in our next
    release. You can find some instructions to build Sage on OS X 10.5
    in 64-bit mode at
    http://mvngu.wordpress.com/2009/09/02/compile-sage-4-1-in-64-bit-mode-on-os-x-10-5-8/


Assumptions: You have a computer with about 2 GB of free
disk space running Linux (32-bit or 64-bit), Mac OS X 10.4, 10.5, or
10.6 with XCode. In particular, under Linux the following standard
command-line development tools must be installed on your computer
(under OS X they all come with XCode):

::

       gcc
       g++
       gfortran
       make
       m4
       perl
       ranlib
       tar
       readline and its development headers
       ssh-keygen -- needed to run the notebook in secure mode.
       latex -- highly recommended, though not strictly required

To check if you have ``m4`` installed, for example, type ``which m4``
at a command line. If it gives an error (or returns nothing), then
it is not installed. It is highly recommended that you have LaTeX
installed, but not required. If you don't have ``ssh-keygen`` on your
local system, then you cannot run the notebook in secure mode. To
run it in insecure mode, run the command ``notebook(secure=False)``
instead of ``notebook()``.

In OS X, make sure you have XCode version at least 2.4, i.e., ``gcc -v``
should output build at least 5363. If you don't, go to
http://developer.apple.com/ sign up, and download the free Xcode
package. Only OS X :math:`$\geq 10.4$` is supported. This will give
you all of the above commands.

On a Debian-based system (e.g., Ubuntu), ranlib is in the binutils
package. On a newly installed Ubuntu system (this was tested on
Ubuntu 9.04), you can install the above commands as follows:

::

     sudo apt-get install build-essential m4 gfortran

It is recommended that you install the readline package and its
corresponding development headers. These packages make it easier to
work with the Sage command line interface by providing text editing
features at the command line level. On a Debian or Ubuntu system, use
the following commands to install the readline library and its
development headers:

::

    sudo apt-get install readline-common libreadline-dev

The LaTeX package and a PDF previewer are optional but they can be
installed using

::

    sudo apt-get install texlive xpdf evince xdvi

(You must have the GNU version of ``make`` installed.
For example, Sage won't build on a FreeBSD install that doesn't
have the optional GNU version of ``make`` installed as well
(and named ``make``).)

Although some of Sage is written in Python, you do not need Python
pre-installed on your computer, since the Sage installation
includes everything you need. When the installation program is run,
it will check that you have each of the above-listed prerequisites,
and inform you of any that are missing.

-  If you want to use Tcl/Tk libraries in Sage,
   do the following preferably before compilation.
   Sage's Python will automatically recognize your system's
   install of Tcl/Tk if it exists. You need to install the
   Tcl/Tk development libraries though, not just the Tck/Tk base.

   On Ubuntu, this is the command::

       sudo apt-get install tk8.5-dev    # or the latest version available

   Now you can install Sage and Sage's Python will automatically
   recognize your system's install of Tcl/Tk. If you forgot
   and installed Sage first anyway, all is not lost.
   Just issue the command::

       sage -f python-2.5.2.p8    # or the latest version available

   after installing the Tcl/Tk development libraries as above.
   If

   .. skip

   ::

       sage: import _tkinter
       sage: import Tkinter

   does not raise an ``ImportError`` then it worked.

-  Sage is currently being developed using GCC version 4.3.x, and
   is likely to compile fine with other GCC versions in the 4.x
   series. It does not work with older GCC releases. If you are
   interested in working on support for Intel or Sun's CC compiler,
   please email the sage-devel mailing list, otherwise known as the
   sage-devel Google group at
   http://groups.google.com/group/sage-devel

-  One reason ``perl`` is required is that both the NTL and PARI
   configuration scripts are written in Perl.



After extracting the Sage tarball, the subdirectory ``spkg`` contains
the source distributions for everything on which Sage depends. We
emphasize that all of this software is included with Sage, so you
do not have to worry about trying to download and install any one
of these packages (such as GAP, for example) yourself.

On tests using various Linux computer systems, the known problems
are:


-  Does not build with gcc 4.3.0 yet, but work is ongoing to fix
   that.

-  Moving the build after compiling breaks the PARI Galois fields
   database, which appears to be hardcoded into the PARI binary.
   (Somebody help fix this!)


Fortran
-------

On Linux and Solaris systems, a working Fortran compiler is required
for building Sage from source. If you are using Fortran on a platform
for which Sage does not include g95 binaries, you must use a
system-wide gFortran. For example, Solaris 10 does not ship with any
Fortran binaries. You need to explicitly tell the Sage build process
about the Fortran compiler and library location. Do this by typing ::

    export SAGE_FORTRAN=/exact/path/to/gfortran
    export SAGE_FORTRAN_LIB=/path/to/fortran/libs/libgfortran.so

Note that the ``SAGE_FORTRAN`` environment variable is supposed to
impact *only* the Fortran Sage package, otherwise known as the Fortran
spkg. Apart from that, this variable is *not* designed to do anything
at all to other spkg's that use Fortran. For example, the Lapack spkg
uses Fortran, but the compilation process of Lapack should ignore the
``SAGE_FORTRAN`` environment variable. The ``SAGE_FORTRAN``
environment variable does not mean "build any spkg that uses Fortran
using this Fortran". It means "when installing the Fortran spkg, setup
the ``sage_fortran`` script to run the Fortran compiler specified by
the ``SAGE_FORTRAN`` variable".

On Mac OS X, you are not required to have a Fortran compiler on your
system. The Sage source distribution is shipped with a Fortran
compiler for Mac OS X. This Fortran compiler is used, unless you
specify another Fortran compiler via the variable ``SAGE_FORTRAN``.

On platforms such as AIX, HP-UX, and Solaris, where both 32- and
64-bit builds are supported, the library path variable
``SAGE_FORTRAN_LIB`` must point to the 32-bit library if you are
building Sage in 32-bit. Also, ``SAGE_FORTRAN_LIB`` must point to a
64-bit library if you are building Sage in 64-bit. For example, on
Solaris both of the variables ``SAGE_FORTRAN`` and
``SAGE_FORTRAN_LIB`` could be set as follows (you need to check
this)::

    # SPARC and x86
    SAGE_FORTRAN=/path/to/gcc/install/directory/bin/gfortran

    # 32-bit SPARC
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/libgfortran.so

    # 64-bit SPARC
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/sparcv9/libgfortran.so

    # 32-bit x86
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/libgfortran.so

    # 64-bit x64
    SAGE_FORTRAN_LIB=/path/to/gcc/install/directory/lib/amd64/libgfortran.so


Steps to Install from Source
----------------------------

Installation from source is (potentially) very easy, because the
distribution contains (essentially) everything on which Sage
depends.

Make sure there are no spaces in the directory name under which you
build. Running from a directory with spaces in its name is supported but
discouraged. Building is not possible, since several of the
components do not build if there are spaces in the path.



#. Go to http://www.sagemath.org/download-source.html , select a mirror,
   and download the file sage-\*.tar.

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

#. Type

   ::

             make

   This compiles Sage and all dependencies. Note that you do not need
   to be logged in as root, since no files are changed outside of the
   ``sage-x.y.z`` directory. [1]_ This command does the usual steps for
   each of the packages, but puts all the results in the local build
   tree. This can take close to an hour on some machines. Depending on the
   architecture of your system (e.g., Celeron, Pentium Mobile, Pentium 4,
   etc.), it can take over three hours to build Sage from source.  If the
   build is successful, you will not see the word ERROR in the last 3-4 lines
   of output.

.. [1]
   There is one exception--the ``.ipythonrc`` directory is created in
   your ``HOME`` directory if it doesn't exist.


       The directory where you built Sage is NOT hardcoded. You should
       be able to safely move or rename that directory. (It's a bug if
       this is not the case --- unfortunately there is one
       bug which hasn't yet been fixed along these lines, namely the PARI
       install hard-codes the location of the "galois data" files. Fixes
       welcome!)


   After you build Sage, you may optionally copy or move the entire
   build tree to ``/usr/local``. You might also copy the ``sage-*/sage``
   script to ``/usr/local/bin/`` and edit ``ROOT="....."`` at the top of
   that file.

#. To start Sage, change into the Sage home directory and type:

   ::

             ./sage

   You should see the Sage prompt, which will look something like this
   (starting the first time can take a few seconds):

   ::

       $ sage
       ----------------------------------------------------------------------
       | SAGE Version 3.1, Release Date: 2008-08-16                         |
       | Type notebook() for the GUI, and license() for information.        |
       ----------------------------------------------------------------------
       sage:

   Just starting successfully tests that many of the components built
   correctly. If the above is not displayed (e.g., if you get a
   massive traceback), please report the problem, e.g., to
   http://groups.google.com/group/sage-support . Please include in
   your email the file ``install.log``. It would also be helpful to
   include the type of operating system you have and the version
   number (and date) of the copy of Sage you are using. (There are no
   formal requirements for bug reports - just send them; we appreciate
   everything.)

   After starts, try a command:

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
       GAP4, Version: 4.4.6 of 02-Sep-2005, x86_64-unknown-linux-gnu-gcc
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
        A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                              0<
            by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
       FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
       // ** executing /usr/local/sage/sage-0.8.2/bin/LIB/.singularrc
       [ctrl-d]
       > Auf Wiedersehen.
       sage:

#. Optional: Check the interfaces to any non-included software that
   you have available. Note that each interface calls its
   corresponding program by a particular name: Mathematica is invoked
   by calling ``math``, Maple by calling ``maple``, et cetera. The
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


   -  Copy ``$SAGE_ROOT/sage`` to a location in your ``PATH``. If you do
      this, make sure you edit the line with the ``....``'s at the top of
      the ``sage`` script.

   -  For KDE users, create a bash script {sage} containing the lines

      ::

          #!/bin/bash
          konsole -T "sage" -e <SAGE_ROOT>/sage

      which you make executable (``chmod a+x sage``) and put it somewhere in
      your path. (Note that you have to change ``$SAGE_ROOT`` above!) You
      can also make a KDE desktop icon with this as the command (under
      the Application tab of the Properties of the icon, which you get my
      right clicking the mouse on the icon).

   -  For bash shell users, type ``echo $PATH`` and
      ``cp sage <your-path-dir>`` into one of these directories, or else
      add this ``bin`` directory to your ``PATH`` variable, e.g., if you use
      the bash shell, add the line

      ::

          PATH="<sage-home-dir>/bin":$PATH
          export PATH

      in your .bashrc file (if it exists; if not, make one). After doing
      this and logging out and in again, typing ``sage`` at a shell prompt
      should start Sage.

   - On Linux and OS X systems, you can make an alias to ``$SAGE_ROOT/sage``.
     For example, put something similar to the following line in your
     ``.bashrc`` file:

     ::

         alias 'sage'='/home/username/sage-3.1.2/sage'

     Having done so, quit your terminal emulator and restart it again.
     Now typing ``sage`` within your terminal emulator should start
     Sage.

#. Optional: Test the install by typing ``./sage -testall``. This
   runs most examples in the source code and makes sure that they run
   exactly as claimed. To test all examples, use
   ``./sage -testall -optional -long``; this will run examples that take
   a long time, and those that depend on optional packages and
   software, e.g., Mathematica or Magma. Some (optional) examples will
   likely fail because they assume that a database is installed.
   Alternatively, from within ``$SAGE_ROOT``, you can type
   ``make test`` to run all the standard test code.  This can take
   from 30 minutes to an hour or longer.

#. Optional: The directory ``spkg/build`` contains intermediate code
   that is used to build sage. Type ``make clean`` to delete it and a
   few other directories (e.g., ``spkg/archive`` and ``devel/old``). This
   is safe and will save you about 500 MB of disk space. You may wish to
   type this periodically.

#. Optional: Install optional Sage packages and databases. Type
   ``sage -optional`` to see a list or visit
   http://www.sagemath.org/packages/optional/, and
   ``sage -i <package name>`` to automatically download and install a
   given package.

#. Optional: Run the ``install_scripts`` command from within Sage to create
   gp, singular, gap, etc., scripts in your ``PATH``. Type
   ``install_scripts?`` in Sage for details.


Have fun! Discover some amazing conjectures!

Installation in a Multiuser Environment
---------------------------------------

This section addresses the question of how a system administrator
can install a single copy of Sage in a multi-user computer
network.

System-wide install
~~~~~~~~~~~~~~~~~~~

This is a compilation of posts to the Sage support list (in
particular those of Luis Finotti).


#. Unpack the current Sage tarball (we shall assume it is
   ``sage-2.5.2.tar``) at, e.g., ``/usr/local/`` and compile it as root.
   Assuming you are in a root shell and the tarball is in your current
   directory, type:

   ::

       cp sage-2.5.2.tar /usr/local
       cd /usr/local
       tar xvf sage-2.5.2.tar
       cd sage-2.5.2/
       make

    (Comment: It's better to build in place.  It's a bug if anything goes
    wrong when relocating the entire tarball -- unfortunately there
    is one bug I haven't fixed along these lines, namely the
    PARI install hard-codes the location of the "galois data" files.
    (Fixes welcome!))

#. Make sure to modify the line with the ``.....``"'s at the top of the
   ``sage`` script. In other words, edit ``SAGE_ROOT="....."`` to say
   ``SAGE_ROOT="/usr/local/sage-2.5.2"``.

#. There are some initial files that have to be created during the
   first run of Sage. Try starting up Sage once as root (or, to be
   more thorough, try ``make test`` as root to run all the standard test
   code). You can stop the tests by pressing ``ctrl-z`` followed by
   typing ``kill %1`` (assuming you had no other jobs in the
   background of that shell).

#. Make a copy of the ``sage`` script in ``/usr/local/bin``:

   ::

       cp /usr/local/sage-2.5.2/sage /usr/local/bin/

   You make a copy instead of a symlink, since upgrading with
   ``sage -upgrade`` overwrites ``/usr/local/sage-2.5.2/sage``, hence
   deleting the ``ROOT=...`` part of that file.

   Make sure that all files in ``/usr/local/sage-2.5.2`` are readable by
   all:

   ::

       chmod a+rX -R /usr/local/sage-2.5.2


Special Notes
-------------


-  (Found by Dorian Raymer) Sage will not build if you have only
   bison++. You should uninstall bison++ and install bison.

-  (Found by Peter Jipsen) If you get an error like

   ::

       ImportError: /home/jipsen/Desktop/sage-1.3.3.1/local/lib/libpari-gmp.so.2:
            cannot restore segment prot after reloc:
       Permission denied

   then your SELinux configuration is preventing Sage from launching. To
   rectify this issue, you can either change the default security
   context for Sage (??) or disable SELinux altogether by setting the
   line ``SELINUX=disabled`` in your ``/etc/sysconfig/selinux`` file.

- To make SageTeX available to your users, see the instructions for
  :ref:`installation in a multiuser environment
  <sagetex_installation_multiuser>`.
