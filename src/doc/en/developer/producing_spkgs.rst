.. _chapter-spkg:

============================
Producing New Sage Packages
============================

If you are producing code to add new functionality to Sage, you
might consider turning it into a package (an ``spkg`` file)
instead of a patch file. If your code is very large (for instance)
and should be offered as an optional download, a package is the
right choice; similarly, if your code depends on some other
optional component of Sage, you should produce a package.

If you're not sure whether to build an spkg file or a patch file,
ask for advice on ``sage-devel``.

Creating a New spkg File
========================

Sage packages are distributed as .spkg files, but an .spkg file is
just a .tar.bz2 file (or a tar file), but named with .spkg to
discourage confusion. In particular, you can type

::

         tar jxvf mypackage-version.spkg

to extract one and see what is inside.

Here's how to make your own spkg file:


-  Make a directory, e.g., ``mypackage-0.1``. The name of the
   directory should be a lower-case string with no dashes, followed by
   a dash, followed by a version number.

-  Put your files in that directory.

-  Create an executable shell script
   ``mypackage-0.1/spkg-install`` following the skeleton given
   below.

   ::

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

-  The script ``spkg-install`` is run during installation of
   the Sage package. In this script you may make the following
   assumptions:


   -  The PATH has the locations of ``sage`` and
      ``python`` (from the Sage installation) at the front. Thus
      the command

      ::

               python setup.py install

      will run the correct version of Python with everything set up
      correctly. Also, running ``gap`` or ``Singular``, for
      example, will run the correct version.

   -  The environment variable ``SAGE_ROOT`` points to the root
      directory of the Sage installation.

   -  The environment variable ``SAGE_LOCAL`` points to the
      ``SAGE_ROOT/local`` directory of the Sage installation.

   -  The environment variables ``LD_LIBRARY_PATH`` and
      ``DYLD_LIBRARY_PATH`` both have
      ``SAGE_ROOT/local/lib`` at the front.


-  The ``spkg-install`` script should copy your files to the
   appropriate place after doing any build that is necessary. That's
   it!

-  (Optional) Post a copy on the Sage trac server (:ref:`chapter-trac`)
   so it can be refereed.  If it gets a positive review, it might be
   included into the core Sage library, or it might become an optional
   download from the Sage web site, so anybody can automatically
   install it by typing ``sage -i mypackage-version.spkg``.


If your package depends on another package, say boehmgc, then you
should check that this other package has been installed. Your
``spkg-install`` script should check that it exists, with the
code like the following:

::

    BOEHM_GC=`cd $SAGE_ROOT/spkg/standard/; ./newest_version boehm_gc`
    if [ $? -ne 0 ]; then
        echo "Failed to find boehm_gc.  Please install the boehm_gc spkg"
        exit 1
    fi

*Caveat*: Do not just copy to e.g.
``SAGE_ROOT/local/lib/gap\*/`` since that will copy your
package to the lib directory of the old version of GAP if GAP is
upgraded.

External Magma code goes in ``SAGE_ROOT/data/extcode/magma/user``,
so if you want to redistribute Magma code with Sage as a package
that Magma-enabled users can use, that's where you'd put it. You'd
also want to have relevant Python code to make the Magma code
easily usable.
