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
to patch the Matplotlib package ``matplotlib-1.0.1.p0``. Note that "p0" denotes
the patch level of the spkg, while "1.0.1" refers to the upstream
version of Matplotlib as contained under ``matplotlib-1.0.1.p0/src/``. The
installation script of that spkg is ::

    matplotlib-1.0.1.p0/spkg-install

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
upstream project. For example, the source code of Matplotlib 1.0.1 is
contained under ::

    matplotlib-1.0.1.p0/src/

To patch the upstream source code, you should edit a copy of the
relevant file -- files in the ``src/`` directory should be untouched,
"vanilla" versions of the source code.  For example, you might copy
the entire ``src/`` directory::

    $ pwd
    matplotlib-1.0.1.p0
    $ cp -pR src src-patched

Then edit files in ``src-patched/``.  Once you are satisfied with your
changes, generate a unified diff between the original file and the
edited one, and save it in ``patches/``::

    $ diff -u src/configure src-patched/configure > patches/configure.patch

Save the unified diff to a file with the same name as the source file
you patched, but using the file extension ".patch". Note that the
directory ``src/`` should not be under revision control, whereas
``patches/`` must be under revision control. The Mercurial
configuration file ``.hgignore`` should contain the following line::

    src/

Ensure that the installation script ``spkg-install`` contains code to
apply the patches to the relevant files under ``src/``. For example,
the file ::

    matplotlib-1.0.1.p0/patches/finance.py.patch

is a patch for the file ::

    matplotlib-1.0.1.p0/src/lib/matplotlib/finance.py.patch

The installation script ``matplotlib-1.0.1.p0/spkg-install`` contains the
following code to install the relevant patches::

    # Apply patches.  See SPKG.txt for information about what each patch
    # does.
    for patch in ../patches/*.patch; do
        patch -p1 <"$patch"
        if [ $? -ne 0 ]; then
            echo >&2 "Error applying '$patch'"
            exit 1
        fi
    done

Of course, this could be modified if the order in which the patches
are applied is important, or if some patches were platform-dependent.
For example::

    if [ "$UNAME" = "Darwin" ]; then
        for patch in ../patches/darwin/*.patch; do
            patch -p1 <"$patch"
            if [ $? -ne 0 ]; then
                echo >&2 "Error applying '$patch'"
                exit 1
            fi
        done
    fi

(The environment variable :envvar:`UNAME` is defined by the script
``sage-env``, and is available when ``spkg-install`` is run.)

Now provide a high-level explanation of your changes in
``SPKG.txt``. Once you are satisfied with your changes, use Mercurial
to check in your changes with a meaningful commit message. Next,
increment the patch level of the spkg by one, e.g. rename the
directory ``matplotlib-1.0.1.p0`` to ``matplotlib-1.0.1.p1``.  To
produce the actual spkg file, change to the parent directory of
``matplotlib-1.0.1.p1`` and execute ::

    $ /path/to/sage-x.y.z/sage --pkg matplotlib-1.0.1.p1
    Creating Sage package matplotlib-1.0.1.p1

    Created package matplotlib-1.0.1.p1.spkg.

        NAME: matplotlib
     VERSION: 1.0.1.p1
        SIZE: 11.8M
     HG REPO: Good
    SPKG.txt: Good

Spkg files are either bzipped tar files or just plain tar files; the
command ``sage --pkg ...`` produces the bzipped version.  If your spkg
contains mostly binary files which will not compress well, you can use
``sage --pkg_nc ...`` to produce an uncompressed version, i.e., a
plain tar file::

    $ sage --pkg_nc matplotlib-1.0.1.p0/
    Creating Sage package matplotlib-1.0.1.p0/ with no compression

    Created package matplotlib-1.0.1.p0.spkg.

        NAME: matplotlib
     VERSION: 1.0.1.p0
        SIZE: 32.8M
     HG REPO: Good
    SPKG.txt: Good

Note that this is almost three times the size of the compressed
version, so we should use the compressed version!

At this point, you might want to submit your patched spkg for review.
So provide a URL to your spkg on the relevant trac ticket and/or in an
email to the relevant mailing list. Usually, you should not upload
your spkg itself to the relevant trac ticket -- don't post large
binary files to the trac server.


Use patch for patching
======================

A main message of this section is: use the GNU program ``patch`` to
apply patches to files in ``src/``.  GNU patch is distributed with
Sage, so if you are writing an spkg which is not part of the standard
Sage distribution, you may use ``patch`` in the ``spkg-install``
script freely.  If you are working on an spkg which is (or will be) a
standard spkg in Sage, then you should make sure that ``patch`` is
listed as a dependency for your spkg in the makefile
``SAGE_ROOT/spkg/standard/deps``.

See the section :ref:`section-spkg-patching-overview` for information
about how to produce patch files in the directory ``patches/``, and
how to apply them in ``spkg-install``.


Bumping up an spkg's version
============================

If you want to bump up the version of an spkg, you need to follow some
naming conventions. Use the name and version number as given by the
upstream project, e.g. ``matplotlib-1.0.1``. If the upstream package is
taken from some revision other than a stable version, you need to
append the date at which the revision is made, e.g. the Singular
package ``singular-3-1-0-4-20090818.p3.spkg`` is made with the
revision as of 2009-08-18. If you start afresh from an upstream
release without any patches to its source code, the resulting spkg
need not have any patch-level labels (appending ".p0" is allowed, but
is optional). For example, ``sagenb-0.6.spkg``
is taken from the upstream stable version ``sagenb-0.6`` without any
patches applied to its source code. So you do not see any patch-level
numbering such as ``.p0`` or ``.p1``.

Say you start with ``matplotlib-1.0.1.p0`` and you want to replace Matplotlib
1.0.1 with version 1.0.2. This entails replacing the source code for Matplotlib 1.0.1 under
``matplotlib-1.0.1.p0/src/`` with the new source code. To start with, follow the
naming conventions as described in the section
:ref:`section-spkg-patching-overview`. If necessary, remove any
obsolete patches and create any new ones, placing them
in the ``patches/`` directory.  Modify the script
``spkg-install`` to take any changes to the patches into account; you
might also have to deal with changes to how the new version of the
source code builds. Then package your replacement spkg using
the Sage command line options ``--pkg`` or ``--pkg_nc`` (or tar and
bzip2).

To install your replacement spkg, you use ::

    sage -f /URL/to/package-x.y.z.spkg

or ::

    sage -f /path/to/package-x.y.z.spkg

To compile Sage from source with the replacement (standard) spkg,
untar a Sage source tarball, remove the existing spkg under
``SAGE_ROOT/spkg/standard/``. In its place, put your replacement
spkg. Then execute ``make`` from ``SAGE_ROOT``.
