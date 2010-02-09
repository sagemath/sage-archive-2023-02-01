.. _chapter-patching-spkgs:

=======================
Patching a Sage Package
=======================

This chapter provides guidelines on patching an existing spkg. Also
covered are steps for upgrading an upstream project's source
distribution, as contained under the subdirectory ``src/``, to the
latest upstream release. For information on creating a new spkg, see
the chapter :ref:`chapter-spkg`.


.. _section-spkg-patching-overview:

Overview of patching spkg's
===========================

Make sure you are familiar with the structure and conventions relating
to spkg's; see the chapter :ref:`chapter-spkg` for details. Patching
an spkg involves patching the installation script of the spkg and/or
patching the upstream source code contained in the spkg. Say you want
to patch the GAP package ``gap-4.4.10.p17``. Note that "p17" denotes
the patch level of the spkg, while "4.4.10" refers to the upstream
version of GAP as contained under ``gap-4.4.10.p17/src/``. The
installation script of that spkg is ::

    gap-4.4.10.p17/spkg-install

In general, a script with the name ``spkg-install``  is an
installation script for an spkg. To patch the installation script, use
a text editor to edit that script. Then in the log file ``SPKG.txt``,
provide a high-level description of your changes. Once you are
satisfied with your changes in the installation script and the log
file ``SPKG.txt``, use Mercurial to check in your changes and make
sure to provide a meaningful commit message. See the section
:ref:`section-submitting-change` for guidelines relating to commit
messages.

The directory ``src/`` contains the source code provided by the
upstream project. For example, the source code of GAP 4.4.10 is
contained under ::

    gap-4.4.10.p17/src/

To patch the upstream source code, you should first copy the relevant
file under ``src/`` over to the directory ``patches/``. Then edit that
copied file rather than the relevant file under ``src/``. Once you are
satisfied with your changes to the copied file under ``patches/``,
generate a unified diff between the original file under ``src/`` and
the corresponding copied (and edited) file under ``patches/``. Save
the unified diff to a file with the same name as the source file you
patched, but using the file extension ".patch", and place the diff
file under ``patches/``. Note that the directory ``src/`` should not
be under revision control, whereas ``patches/`` must be under revision
control. The Mercurial configuration file ``.hgignore`` should contain
the following line ::

    src/

Ensure that the installation script ``spkg-install`` contains code to
copy over the patched file under ``patches/`` to the relevant place
under ``src/``. For example, the file ::

    gap-4.4.10.p17/patches/configure.out-ia64

is a patched version of the file ::

    gap-4.4.10.p17/src/cnf/configure.out

The installation script ``gap-4.4.10.p17/spkg-install`` contains the
following conditional to test under which conditions to copy over the
patched file to the relevant place under ``src/``::

    # If we're on an Itanium Linux box, we overwrite configure.out
    # with a slightly modified version.  The modified version has
    # all -O2's replaced by -O0. See
    # http://www.gap-system.org/Faq/Hardware-OS/hardware-os8.html
    # On the San Diego Super computer `uname -p` is unknown, but
    # uname -a includes ia64 in the output.  So this is a better
    # detection method. Note that GAP was "fixed" to work fine on
    # Itanium without this -O0 hack, but with GCC-4.4.0, GAP
    # mysterious stopped working again.  So we revert to using -O0.
    if [ `uname` = "Linux" ]; then
       uname -a | grep ia64
       if [ $? = 0 ] || [ `uname -p` = "ia64" ]; then
           cp patches/configure.out-ia64 src/cnf/configure.out
           echo "The file configure.out was patched for SAGE!" > src/cnf/configure.out.README
       fi
    fi

Now provide a high-level explanation of your changes in
``SPKG.txt``. Once you are satisfied with your changes, use Mercurial
to check in your changes with a meaningful commit message. Next,
increment the patch level of the spkg by one, e.g. rename
``gap-4.4.10.p17`` to ``gap-4.4.10.p18``. If you want your patched
spkg to be reviewed and available in the next release of Sage, you
need to compress it using tar and bzip2. Under Linux and Mac OS X, you
can compress your patched spkg as follows::

    tar -jcf gap-4.4.10.p18.spkg gap-4.4.10.p18

If your system does not have Sage installed, the above command can be
used for compressing your patched spkg. In case you have a local Sage
installation, you can achieve the equivalent effect using ::

    /path/to/sage-x.y.z/sage -pkg gap-4.4.10.p18

which would produce a tarball compressed using bzip2. Sage also
accepts the command line option ``-pkg_nc``, which is useful when your
spkg mostly contains binary files. Both of the command line options
``-pkg`` and ``-pkg_nc`` automatically appends the extension ``.spkg``
to the resulting (compressed) spkg. However, using tar and bzip2, you
need to manually put in the extension ``.spkg`` instead of
``.tar.bz2``.

As part of your request for your patched spkg to be reviewed, you need
to provide a URL to your spkg on the relevant trac ticket and/or in an
email to the relevant mailing list. Usually, you should not upload
your spkg to the relevant trac ticket.


Use cp for patching
===================

A main message of this chapter is: use the command ``cp`` to copy a
patched source file over to the relevant place under ``src/``. This
copy operation is how you should patch an upstream source and it must
be implemented in ``spkg-install``. In the installation script, do not
use the command ``patch`` for patching.

Say you have the following two lines in ``spkg-install`` for patching
the upstream source::

    cp patches/allfaces.c src/src/
    patch -p0 < patches/cdd_both_reps-make.patch

The first line is the preferred way to patch an upstream source
because it copies the patched file ``patches/allfaces.c`` over to
``src/src/``. Even before running the installation script
``spkg-install``, the original upstream source is the file
``src/src/allfaces.c`` and its corresponding patched version for Sage
is ``patches/allfaces.c``. And you patch the upstream source using
``cp``.

When you are packaging an spkg, ensure that your patch(es) are not
already applied to the relevant files under ``src/``. Similarly, a
patched version of an upstream source file must not already be copied
over to the relevant place under ``src/``. The patching process should
be implemented in ``spkg-install``, which is run during the
(re)installation of an spkg. That is, commands required for patching
and installing an spkg must be implemented in
``spkg-install``. **Never** apply patches to upstream source under
``src/`` (or copy patched files over to ``src/``) and then package up
your patched spkg. There must be a **clear separation** between the
pristine source as provided by the upstream project, and any patched
version of the upstream project's source files. The patching process
is executed at run time by the script ``spkg-install``. Implement all
the necessary patching logic in that script.

The process of copying patched files under ``patches/`` over to
``src/`` should be implemented in ``spkg-install`` using a function
named ``patch()``. This function must be called prior to running the
configuration and installation scripts of the upstream project. For
example, the file ``spkg-install`` of the cddlib spkg contains this
Bash function definition::

    # Patches to apply on top of the clean upstream source under src/.
    patch() {
        # We currently do not use this. To be removed at a later date.
	cp patches/allfaces.c src/src/

	# Required by sage.geometry.polyhedra
	cp patches/cdd_both_reps.c src/src/
	cp patches/cdd_both_reps.c src/src-gmp/
    }

    <SNIP>

    # apply patches on top of pristine upstream release under src/
    patch

In other words, at run time, patched versions of upstream source files
are to be copied over to the appropriate place under ``src/``. In this
way, there is no reliance on the Unix command ``patch``. Here are some
reasons why ``cp`` is used in preference to ``patch``:

* A version of GNU ``patch`` that can apply unified diffs does not
  exist on some platforms by default, e.g. Solaris.

* The number of changes to the sources in the ``src/`` directory are
  usually small. In many cases, changes that the Sage project makes to
  an upstream source can be reported to the relevant upstream project.

* Copying patched source files is dead simple as it follows the KISS
  principle.

* Mercurial's ``patch`` requires the sources to be under revision
  control, which unnecessarily increase the size of a Sage source
  tarball.


.. _section-unconditional-patching:

Making an unconditional patch
=============================

The section :ref:`section-spkg-patching-overview` describes how to
produce a conditional patch. That is, at run time you only want to
apply the patch to the spkg if some condition is met. This section
covers how to produce an unconditional patch, i.e. a patch that is to
be applied to an spkg no matter what.

For an unconditional patch, in the installation script
``spkg-install`` you do not need to surround the copy command with a
conditional. But immediately after the copy line, you need to test
that the copy process is successful. For example, here is a block in
``python-2.6.4.p5/spkg-install`` that copies over an unconditional
patch::

    cp patches/locale.py src/Lib/locale.py
    if [ $? -ne 0 ]; then
        echo "Error copying patched locale.py"
    	exit 1
    fi

For both conditional and unconditional patches, you must put the
patched file under ``patches/`` and leave the original file where it
is under ``src/``. The idea is to have a line in ``spkg-install`` that
copies the patched file from ``patches/`` over to the relevant place
under ``src/``. In the case of the above Python spkg, the file
``src/Lib/locale.py`` is the original upstream source and
``patches/locale.py`` is a patched version. In this way, every file in
the upstream release is kept untouched under ``src/`` and all patched
files are placed under ``patches/``. In the installation script
``spkg-install``, provide relevant code (e.g. Bash, Python, and so on)
to copy the relevant patched file under ``patches/`` over to the
appropriate place under ``src/``.


Bumping up an spkg's version
============================

If you want to bump up the version of an spkg, you need to follow some
naming conventions. Use the name and version number as given by the
upstream project, e.g. ``gap-4.4.12``. If the upstream package is
taken from some revision other than a stable version, you need to
append the date at which the revision is made, e.g. the Singular
package ``singular-3-1-0-4-20090818.p3.spkg`` is made with the
revision as of 2009-08-18. If you start afresh from an upstream
release without any patches to its source code, the resulting spkg
need not have any patch-level labels. For example, ``sagenb-0.6.spkg``
is taken from the upstream stable version ``sagenb-0.6`` without any
patches applied to its source code. So you do not see any patch-level
numbering such as ``.p0`` or ``.p1``.

Say you start with ``gap-4.4.10.p17`` and you want to replace GAP
4.4.10 with version 4.4.12. This entails replacing GAP 4.4.10 under
``gap-4.4.10.p17/src/``  with GAP 4.4.12. To start with, follow the
naming conventions as described in the section
:ref:`section-spkg-patching-overview`. If necessary, produce any
relevant patched source files against upstream source files and place
the resulting patched files under ``patches/``. Ensure to patch
``spkg-install`` as well to copy the patched files over to the
relevant place under ``src/``. To get the upgraded upstream source to
build, you might also need to patch the installation script
``spkg-install``. Then package your replacement spkg using tar and
bzip2, or the Sage command line options ``-pkg`` or ``-pkg_nc``.

To install your replacement spkg, you use ::

    sage -f /URL/to/package-x.y.z.spkg

or ::

    sage -f /path/to/package-x.y.z.spkg

To compile Sage from source with the replacement (standard) spkg,
untar a Sage source tarball, remove the existing spkg under
``SAGE_ROOT/spkg/standard/``. In its place, put your replacement
spkg. Then execute ``make`` from ``SAGE_ROOT``.
