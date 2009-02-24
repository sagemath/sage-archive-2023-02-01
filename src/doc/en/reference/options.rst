Sage: Command Line Arguments
============================

There are several flags you can pass to Sage via the command line:


-  ``-h, -?`` Prints a help message.

-  ``-notebook`` Starts the Sage notebook running in default mode.
   This does not open a browser; it just starts a Sage Notebook server
   running on http://localhost:8000 (or a subsequent port if that port
   is not available). Pressing {Control-c} stops the server.

-  ``-i [packages]`` Installs the given Sage packages. If a package
   is listed at
   http://modular.math.washington.edu/sage/packages/optional/ and not
   available in your current directory, then Sage will try to download
   the package and install it.

-  ``-optional`` Lists all optional packages that can be
   downloaded.

-  ``-t <files|dir>`` Tests examples in .py, .pyx, .sage or .tex
   files.

-  ``-testall`` Tests all examples in Sage.

-  ``-upgrade`` Downloads the latest non-optional Sage packages
   from http://modular.math.washington.edu/sage/packages/standard/,
   then builds and installs them. In many cases this requires that
   your computer has the necessary software development tools (listed
   in the installation documentation). Optional packages have to
   upgraded manually using ``sage -i``. (There are no
   binary packages yet.)

-  ``file1.sage file2.sage file.py ...`` Starts Sage, compiles any
   .sage files to .py files, and runs files.

-  ``-advanced`` Displays additional options that can be passed to
   Sage.




