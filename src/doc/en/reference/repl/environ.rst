Environment variables used by Sage
==================================

Sage uses several environment variables when running.  These all have
sensible default values, so many users won't need to set any of these.
(There are also variables used to compile Sage; see the Sage
Installation Guide for more about those.)

- :envvar:`DOT_SAGE` -- this is the directory, to which the user has
  read and write access, where Sage stores a number of files.  The
  default location is ``~/.sage/``, but you can change that by setting
  this variable.

- :envvar:`SAGE_RC_FILE` -- a shell script which is sourced after
  Sage has determined its environment variables.  This script is
  executed before starting Sage or any of its subcommands (like
  ``sage -i <package>``).  The default value is ``$DOT_SAGE/sagerc``.

- :envvar:`SAGE_STARTUP_FILE` -- a file including commands to be
  executed every time Sage starts.  The default value is
  ``$DOT_SAGE/init.sage``.

- :envvar:`SAGE_SERVER` -- if you want to install a Sage package using
  ``sage -i PKG_NAME``, Sage downloads the file from the web, using
  the address ``http://www.sagemath.org/`` by default, or the address
  given by :envvar:`SAGE_SERVER` if it is set.  If you wish to set up
  your own server, then note that Sage will search the directories
  ``SAGE_SERVER/packages/standard/``,
  ``SAGE_SERVER/packages/optional/``,
  ``SAGE_SERVER/packages/experimental/``, and
  ``SAGE_SERVER/packages/archive/`` for packages.  See the script
  :file:`$SAGE_ROOT/spkg/bin/sage-spkg` for the implementation.

- :envvar:`SAGE_PATH` -- a colon-separated list of directories which
  Sage searches when trying to locate Python libraries.

- :envvar:`SAGE_BROWSER` -- on most platforms, Sage will detect the
  command to run a web browser, but if this doesn't seem to work on
  your machine, set this variable to the appropriate command.

- :envvar:`SAGE_ORIG_LD_LIBRARY_PATH_SET` -- set this to something
  non-empty to force Sage to set the :envvar:`LD_LIBRARY_PATH` before
  executing system commands.

- :envvar:`SAGE_ORIG_DYLD_LIBRARY_PATH_SET` -- similar, but only used
  on Mac OS X to set the :envvar:`DYLD_LIBRARY_PATH`.

- :envvar:`SAGE_CBLAS` -- used in the file
  :file:`SAGE_ROOT/src/sage/misc/cython.py`.  Set this to the
  base name of the BLAS library file on your system if you want to
  override the default setting.  That is, if the relevant file is
  called :file:`libcblas_new.so` or :file:`libcblas_new.dylib`, then
  set this to "cblas_new".
