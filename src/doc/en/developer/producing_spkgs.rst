.. _chapter-spkg:

===========================
Producing New Sage Packages
===========================

If you are producing code to add new functionality to Sage, you might
consider turning it into a package (an "spkg") instead of a patch
file. If your code is very large (for instance) and should be offered
as an optional download, a package is the right choice. Similarly, if
your code depends on some other optional component of Sage, you should
produce a package. When in doubt, ask for advice on the ``sage-devel``
mailing list.

This chapter covers issues relevant to producing a package. The
directory structure of a package is discussed along with scripts for
installing a package and running the test suite (if any) contained in
an upstream project's source distribution. For guidelines on patching
an existing Sage package, see the chapter
:ref:`chapter-patching-spkgs`.


Creating a new spkg
===================

The abbreviation "spkg" stands for "Sage package". The directory
``SAGE_ROOT/spkg/standard`` contains spkg's. In a source install,
these are all Sage spkg files (actually ``.tar`` or ``.tar.bz2``
files), which are the source code that defines Sage. In a binary
install, some of these may be small placeholder files to save space.

Sage packages are distributed as ``.spkg`` files, which are
``.tar.bz2`` files (or ``tar`` files) but have the extension ``.spkg``
to discourage confusion. Although Sage packages are packed using tar
and/or bzip2, note that ``.spkg`` files contain control information
(installation scripts and metadata) that are necessary for building
and installing them.  When you compile Sage from a source distribution
(or when you run ``sage -i <pkg>`` or ``sage -f <pkg>``),
the file ``SAGE_ROOT/spkg/bin/sage-spkg`` takes care of the unpacking,
compilation, and installation of Sage packages for you. You can
type

::

    tar -jxvf mypackage-version.spkg

to extract an spkg and see what is inside.  If you want to create a
new Sage package, it is recommended that you start by examining some
existing spkg's. In a source distribution of Sage, the standard spkg's
can be found under ``SAGE_ROOT/spkg/standard/``. The URL
http://www.sagemath.org/download-packages.html lists standard spkg's
available for download.


Naming your spkg
----------------

Each Sage spkg has a name of the following form:

::

   BASENAME-VERSION.spkg

``BASENAME`` is the name of the package; it may contain lower-case
letters, numbers, and underscores, but no hyphens.  ``VERSION`` is the
version number; it should start with a number and may contain numbers,
letters, dots, and hyphens; it may end in a string of the form
"pNUM", where "NUM" is a non-negative integer.  If your spkg is a
"vanilla" (unmodified) version of some piece of software, say version
5.3 of "my-python-package", then ``BASENAME`` would be
"my_python_package" -- note the change from hyphens to underscores,
because ``BASENAME`` should not contain any hyphens -- and ``VERSION``
would be "5.3".  If you need to modify the software to use it with
Sage (as described below and in the chapter
:ref:`chapter-patching-spkgs`), then ``VERSION`` would be "5.3.p0",
the "p0" indicating a patch-level of 0.  If someone adds more patches,
later, this would become "p1", then "p2", etc.

The string ``VERSION`` must be present.  If you are using a piece
software with no obvious version number, use a date: you can see
several such names among the standard Sage packages:
http://www.sagemath.org/packages/standard/.

To give your spkg a name like this, create a directory called
``BASENAME-VERSION`` and put your files in that directory -- the
next section describes the directory structure.


Directory structure
-------------------

Put your files in a directory with a name like ``mypackage-0.1``, as
described above.  If you are porting
another software package, then the directory should contain a
subdirectory ``src/``, containing an unaltered copy of the package.
Every file not in ``src/`` should be under version control, i.e. checked
into an hg repository.

More precisely, the directory should contain the following:

- ``src/``: this directory contains vanilla upstream code, with a few
  exceptions, e.g. when the spkg shipped with Sage is in effect
  upstream, and development on that code base is happening in close
  coordination with Sage.  See John Cremona's  eclib spkg, for
  instance. The directory ``src/`` must not be under revision control.

- ``.hg``, ``.hgignore``, and ``.hgtags``: The Sage project uses
  Mercurial for its revision control system (see
  :ref:`chapter-mercurial`).  The hidden directory ``.hg`` is part
  of the standard Sage spkg layout.  It contains the Mercurial
  repository for all files not in the ``src/`` directory.
  To create this Mercurial repository from scratch, you should do

  ::

      hg init

  The files ``.hgignore`` and ``.hgtags`` also belong to the
  Mercurial repository.  The file ``.hgtags`` is optional, and is
  frequently omitted.  You should make sure that the file
  ``.hgignore`` contains "src/", since we are not tracking its
  content.  Indeed, frequently this file contains only a single line,

  ::

      src/

- ``spkg-install``: this file contains the install script.
  See :ref:`section-spkg-install` for more information and a template.

- ``SPKG.txt``: this file describes the spkg in wiki format.  Each
  new revision needs an updated changelog entry or the spkg will
  get an automatic "needs work" at review time.  See
  :ref:`section-SPKG-txt` for a template.

- ``spkg-check``: this file runs the test suite.  This is somewhat
  optional since not all spkg's have test suites. If possible, do
  create such a script since it helps isolate bugs in upstream
  packages.

- ``patches/``: this directory contains patches to
  source files in ``src/``.  See :ref:`chapter-patching-spkgs`.
  Patches to files in ``src/`` should be applied in
  ``spkg-install``, and all patches must be documented in
  ``SPKG.txt``, i.e. what they do, if they are platform
  specific, if they should be pushed upstream, etc. To ensure that all
  patched versions of upstream source files under ``src/`` are under
  revision control, the whole directory ``patches/`` must be under
  revision control.

**Never** apply patches to upstream source files under ``src/`` and
then package up an spkg. Such a mixture of upstream source with Sage
specific patched versions is a recipe for confusion. There must be a
**clean separation** between the source provided by the upstream
project and the patched versions that the Sage project generates based
on top of the upstream source.

The only exception to this rule is for *removals* of unused
files or directories.  Some packages contain parts which are not needed
for Sage.  To save space, these may be removed directly from ``src/``.
But be sure to document this in the "Special Update/Build Instructions"
section in ``SPKG.txt``!


.. _section-spkg-install:

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

       if [ -z "$SAGE_LOCAL" ]; then
           echo >&2 "SAGE_LOCAL undefined ... exiting"
           echo >&2 "Maybe run 'sage --sh'?"
           exit 1
       fi

       cd src

       # Apply patches.  See SPKG.txt for information about what each patch
       # does.
       for patch in ../patches/*.patch; do
           patch -p1 <"$patch"
           if [ $? -ne 0 ]; then
               echo >&2 "Error applying '$patch'"
               exit 1
           fi
       done

       ./configure --prefix="$SAGE_LOCAL"
       if [ $? -ne 0 ]; then
           echo >&2 "Error configuring PACKAGE_NAME."
           exit 1
       fi

       $MAKE
       if [ $? -ne 0 ]; then
           echo >&2 "Error building PACKAGE_NAME."
           exit 1
       fi

       $MAKE install
       if [ $? -ne 0 ]; then
           echo >&2 "Error installing PACKAGE_NAME."
           exit 1
       fi

       if [ "$SAGE_SPKG_INSTALL_DOCS" = yes ] ; then
          # Before trying to build the documentation, check if any
          # needed programs are present. In the example below, we
          # check for 'latex', but this will depend on the package.
          # Some packages may need no extra tools installed, others
          # may require some.  We use 'command -v' for testing this,
          # and not 'which' since 'which' is not portable, whereas
          # 'command -v' is defined by POSIX.

          # if [ `command -v latex` ] ; then
          #    echo "Good, latex was found, so building the documentation"
          # else
          #    echo "Sorry, can't build the documentation for PACKAGE_NAME as latex is not installed"
          #    exit 1
          # fi


          # make the documentation in a package-specific way
          # for example, we might have
          # cd doc
          # $MAKE html

          if [ $? -ne 0 ]; then
              echo >&2 "Error building PACKAGE_NAME docs."
              exit 1
          fi
          mkdir -p $SAGE_ROOT/local/share/doc/PACKAGE_NAME
          # assuming the docs are in doc/*
          cp -r doc/* $SAGE_ROOT/local/share/doc/PACKAGE_NAME/
       fi


Note that the first line is ``#!/usr/bin/env bash``; this is important
for portability.  Next, the script checks that ``SAGE_LOCAL`` is
defined to make sure that the Sage environment has been set.  After
this, the script may simply run ``cd src`` and then call either
``python setup.py install`` or the autotools sequence
``./configure && make && make install``, or something else along these
lines.

Sometimes, though, it can be more complicated. For example, you might need
to apply the patches from the ``patches`` directory in a particular order. Also,
you should first build (e.g. with ``python setup.py build``, exiting
if there is an error), before installing (e.g. with ``python setup.py
install``). In this way, you would not overwrite a working older
version with a non-working newer version of the spkg.

When copying documentation to
``$SAGE_ROOT/local/share/doc/PACKAGE_NAME``, it may be necessary to
check that only the actual documentation files intended for the
user are copied.  For example, if the documentation is built from
``.tex`` files, you may just need to copy the resulting pdf files,
rather than copying the entire doc directory.  When generating
documentation using Sphinx, copying the ``build/html`` directory
generally will copy just the actual output intended for the user.


.. _section-SPKG-txt:

The file SPKG.txt
-----------------

The ``SPKG.txt`` file should follow this pattern::

     = name of spkg =

     == Description ==

     Describe the package here.

     == License ==

     Describe the package's license here.

     == SPKG Maintainers ==

     * Mary Smith
     * Bill Jones
     * Leonhard Euler

     == Upstream Contact ==

     Provide information for upstream contact.

     == Dependencies ==

     Put a bulleted list of dependencies here:

     * python
     * readline

     == Special Update/Build Instructions ==

     List patches that need to be applied and what they do

     == Changelog ==

     Provide a changelog of the spkg here, where the entries have this format:

     === mypackage-0.1.p0 (Mary Smith, 1 Jan 2012) ===

      * Patch src/configure so it builds on Solaris. See Sage trac #137.

     === mypackage-0.1 (Leonhard Euler, 17 September 1783) ===

      * Initial release.  See Sage trac #007.

When the directory (say, ``mypackage-0.1``) is ready, the command

::

    sage --pkg mypackage-0.1

will create the file ``mypackage-0.1.spkg``.  As noted above, this
creates a compressed tar file. Running ``sage --pkg_nc mypackage-0.1``
creates an uncompressed tar file.

When your spkg is ready, you should post about it on ``sage-devel``.
If people there think it is a good idea, then post a link to the spkg
on the Sage trac server (see :ref:`chapter-trac`) so it can be
refereed.  Do not post the spkg itself to the trac server: you only
need to provide a link to your spkg.  If your spkg gets a positive
review, it might be included into the core Sage library, or it might
become an optional download from the Sage website, so anybody can
automatically install it by typing ``sage -i mypackage-version.spkg``.

.. note::

   For any spkg:

   - Make sure that the hg repository contains every file outside the
     ``src`` directory, and that these are all up-to-date and committed
     into the repository.

   - Include an ``spkg-check`` file if possible (see `trac ticket #299`_).

   .. _trac ticket #299: http://trac.sagemath.org/sage_trac/ticket/299

.. note::

   - If your package is intended to be a standard Sage spkg, then you
     should make sure that any dependencies for your package are
     recorded in the makefile ``SAGE_ROOT/spkg/standard/deps``.  Also
     add a line for your package to the script
     ``SAGE_ROOT/spkg/install``.  For example, the relevant line for
     the readline package is ::

       READLINE=`newest_version readline`

   - If your package is not a standard package and depends on another
     non-standard package, say ``fricas-1.0.9.spkg``, then
     your package's ``spkg-install`` script should check that the
     other package has been installed, with code like the following::

        if [ ! -f "$SAGE_ROOT/spkg/installed/fricas-1.0.9" ]; then
            echo >&2 "The fricas spkg, version 1.0.9 is required; please install it."
            exit 1
        fi

     If you don't care which version of the fricas spkg is installed,
     you could instead use ::

        if ! ls -1 "$SAGE_ROOT/spkg/installed/" | grep '^fricas-.*' > /dev/null ; then
            echo >&2 "The fricas spkg is required; please install it."
            exit 1
        fi

     (The regular expression matches the package name followed by a
     hyphen and then other characters; in particular, it was chosen so
     that it wouldn't match a package like ``fricasaldor-1.0.9`` whose
     name also starts with "fricas".)

     This could be made more sophisticated, for example testing which
     version of fricas is installed vs. which version is required,
     etc. You could, instead of or in addition to checking the
     existence of the appropriate file in
     ``$SAGE_ROOT/spkg/installed/``, check for the required
     functionality somehow. For instance, the ``spkg-install`` script
     for the ``p_group_cohomology`` package checks whether
     ``database_gap`` is installed using the following::

         SMALL_GROUPS=`echo "SmallGroup(13,1); quit;" | $SAGE_ROOT/sage -gap -b -T | grep "13"`
         if [ "$SMALL_GROUPS" = "" ]; then
             echo "It seems that GAP's SmallGroups library is missing."
             echo "One way to install it is by doing"
             echo "    sage: install_package('database_gap')"
             echo "in a Sage session."
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


.. _section-spkg-avoiding-troubles:

Avoiding troubles
=================

This section contains some guidelines on what an spkg must never do to
a Sage installation. You are encouraged to produce an spkg that is as
self-contained as possible.

#. An spkg must not modify an existing source file in the Sage
   library.
#. Do not allow an spkg to modify another spkg. One spkg can depend on
   other spkg -- see above. You need to first test for the existence of the
   prerequisite spkg before installing an spkg that depends on it.
