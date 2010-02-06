.. _chapter-spkg:

===========================
Producing New Sage Packages
===========================

If you are producing code to add new functionality to Sage, you might
consider turning it into a package (an ``spkg`` file) instead of a
patch file. If your code is very large (for instance) and should be
offered as an optional download, a package is the right
choice. Similarly, if your code depends on some other optional
component of Sage, you should produce a package.

If you are not sure whether to build an spkg file or a patch file,
ask for advice on ``sage-devel``.


Creating a new spkg file
========================

Sage packages are distributed as ``.spkg`` files, but an ``.spkg``
file is just a ``.tar.bz2`` file (or a tar file), but named with the
extension ``.spkg`` to discourage confusion. In particular, you can
type

::

         tar jxvf mypackage-version.spkg

to extract an spkg and see what is inside.  If you want to create a
new ``.spkg`` file, you should start by looking at some existing
ones. In a source code distribution of Sage, the standard spkg's can
be found in ``SAGE_ROOT/spkg/standard``, or you can download ``.spkg``
files from http://www.sagemath.org/download-packages.html.

Here is how to make your own ``.spkg`` file. First, make a directory,
e.g. ``mypackage-0.1``. The name of the directory should be a
lower-case string with no dashes, followed by a dash, followed by a
version number.


Directory structure
-------------------

Put your files in the directory ``mypackage-0.1``.  If you are porting
another software package, then the directory should contain a
subdirectory "src", containing an unaltered copy of the package.
Every file not in "src" should be under version control, i.e. checked
into an hg repository.

More precisely, the directory should contain the following:

- ``src``: this directory contains vanilla upstream code, with a few
  exceptions, e.g. when the spkg shipped with Sage is in effect
  upstream, and development on that code base is happening in close
  coordination with Sage.  See John Cremona's  eclib spkg, for
  instance.

- ``.hg``, ``.hgignore``, and ``.hgtags``: The Sage project uses
  Mercurial for its revision control system (see
  :ref:`chapter-mercurial`).  The hidden directory ``.hg`` is part
  of the standard Sage spkg layout.  It contains the Mercurial
  repository for all files not in the ``src`` directory.

  The files ``.hgignore`` and ``.hgtags`` also belong to the
  Mercurial repository.  The file ``.hgtags`` is optional, and is
  frequently omitted.  You should make sure that the file
  ``.hgignore`` contains "src", since we are not tracking its
  content.  Indeed, frequently this file contains only a single line,

  ::

      src

- ``spkg-install``: this file contains the install script. See below
  for more information and a template.

- ``SPKG.txt``: this file describes the spkg in wiki format.  Each
  new revision needs an updated changelog entry or the spkg will
  get an automatic "needs work" at review time.  See below for a
  template.

- ``spkg-check``: this file runs the test suite.  This is somewhat
  optional since not all spkg's have test suites. If possible, do
  create such a script since it helps isolate bugs in upstream
  packages

- ``patches``: this directory contains patches against ``src``. Each
  file requiring changes (e.g. ``foo.c``) needs to have a diff
  against the original file (e.g. ``foo.c.patch``) for easy rebases
  against new upstream. Updated files should be copied into the
  right place in ``src`` at the start of ``spkg-install``. Please
  document all patches in ``SPKG.txt``, i.e. what they do, if they
  are platform specific, if they should be pushed upstream, etc.


The file spkg-install
---------------------

The script ``spkg-install`` is run during installation of the Sage
package. In this script, you may make the following assumptions:

- The PATH has the locations of ``sage`` and ``python`` (from the Sage
  installation) at the front. Thus the command

  ::

      python setup.py install

  will run the correct version of Python with everything set up
  correctly. Also, running ``gap`` or ``Singular``, for example, will
  run the correct version.

- The environment variable ``SAGE_ROOT`` points to the root directory
  of the Sage installation.

- The environment variable ``SAGE_LOCAL`` points to the
  ``SAGE_ROOT/local`` directory of the Sage installation.

- The environment variables ``LD_LIBRARY_PATH`` and
  ``DYLD_LIBRARY_PATH`` both have ``SAGE_ROOT/local/lib`` at the
  front.

The ``spkg-install`` script should copy your files to the appropriate
place after doing any build that is necessary.  Here is a template::

       #!/usr/bin/env bash

       if [ "$SAGE_LOCAL" = "" ]; then
          echo "SAGE_LOCAL undefined ... exiting";
          echo "Maybe run 'sage -sh'?"
          exit 1
       fi

       cd src

       ./configure --prefix="$SAGE_LOCAL"
       if [ $? -ne 0 ]; then
          echo "Error configuring PACKAGE_NAME."
          exit 1
       fi

       make
       if [ $? -ne 0 ]; then
          echo "Error building PACKAGE_NAME."
          exit 1
       fi

       make install
       if [ $? -ne 0 ]; then
          echo "Error installing PACKAGE_NAME."
          exit 1
       fi

Note that the first line is ``/usr/bin/env bash``; this is important
for portability.  Next, the script checks that ``SAGE_LOCAL`` is
defined to make sure that the Sage environment has been set.  After
this, the script may simply run ``cd src`` and then call either
``python setup.py install`` or the autotools sequence
``./configure && make && make install``, or something else along these
lines.

Often, though, it can be more complicated. For example, it is often
necessary to apply the patches from the ``patches`` directory. Also,
you should first build (e.g. with ``python setup.py build``,  exiting
if there is an error, before installing (e.g. with ``python setup.py
install``). In this way, you would not overwrite a working older
version with a non-working newer version of the spkg.


The file SPKG.txt
-----------------

The ``SPKG.txt`` file should follow this pattern::

     = name of spkg =

     == Description ==

     Describe the package here.

     == License ==

     Describe the package's license here.

     == SPKG Maintainers ==

     List the maintainers here

     == Upstream Contact ==

     Provide information for upstream contact.

     == Dependencies ==

     List the dependencies here

     == Special Update/Build Instructions ==

     List patches that need to be applied and what they do

     == Changelog ==

     Provide a changelog of the spkg here.

When the directory (say, ``mypackage-0.1``) is ready, the command

::

    sage -pkg mypackage-0.1

will create the file ``mypackage-0.1.spkg``.  As noted above, this
creates a compressed tar file. Running ``sage -pkg_nc mypackage-0.1``
creates an uncompressed tar file.

When your spkg is ready, you should post about it on ``sage-devel``.
If people there think it is a good idea, then post a link to the spkg
on the Sage trac server (see :ref:`chapter-trac`) so it can be
refereed.  Do not post the spkg itself to the trac server. You only
need to provide a link to your spkg.  If your spkg gets a positive
review, it might be included into the core Sage library, or it might
become an optional download from the Sage website, so anybody can
automatically install it by typing ``sage -i mypackage-version.spkg``.

.. note::

   There are usually a number of things to do for all spkgs:

   - Make sure that the hg repository contains every file outside the
     ``src`` directory, and that these are all up-to-date and committed
     into the repository.

   - Ensure that ``make install`` is non-parallel, i.e. do
     ``export MAKE=make``.

   - Include an ``spkg-check`` file if possible (see `trac ticket #299`_).

   - Include md5sums for spkgs (see `trac ticket #329`_).

   - Set ``LDFLAGS`` on Mac OS X (see `trac ticket #3349`_).

   .. _trac ticket #299: http://trac.sagemath.org/sage_trac/ticket/299

   .. _trac ticket #329: http://trac.sagemath.org/sage_trac/ticket/329

   .. _trac ticket #3349: http://trac.sagemath.org/sage_trac/ticket/3349

.. note::

   - If your package depends on another package, say boehmgc, then you
     should check that this other package has been installed. Your
     ``spkg-install`` script should check that it exists, with code
     like the following:

     ::

       BOEHM_GC=`cd $SAGE_ROOT/spkg/standard/; ./newest_version boehm_gc`
       if [ $? -ne 0 ]; then
           echo "Failed to find boehm_gc.  Please install the boehm_gc spkg"
           exit 1
       fi

   - *Caveat*: Do not just copy to e.g. ``SAGE_ROOT/local/lib/gap*/``
     since that will copy your package to the lib directory of the old
     version of GAP if GAP is upgraded.

   - External Magma code goes in ``SAGE_ROOT/data/extcode/magma/user``,
     so if you want to redistribute Magma code with Sage as a package
     that Magma-enabled users can use, that is where you would put
     it. You would also want to have relevant Python code to make the
     Magma code easily usable.
