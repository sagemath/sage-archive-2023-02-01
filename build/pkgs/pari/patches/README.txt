Details of which source files are patched in spkg-install, and why.
See below for information on files which used to be patched but are no
longer, or are still but now for a different reason. (Do not delete that!)

======================================================================
Files patched as of pari-2.4.3.svn-12577.p5:
======================================================================

Configuration files:
* Configure: Use "#!/usr/bin/env bash" instead of "#!/bin/sh".  Note
  that this is not strictly necessary, but it hopefully makes the
  script less system-dependent. Since Sage assumes the existence of
  bash anyway, it doesn't hurt either. -- Jeroen Demeyer
* config/get_ld: cwitty: Disable -rpath.
* config/get_tests: John Cremona: Disable testing of ellglobalred in
                    "make test-all" in spkg-check, since it requires
                    the elldata database which we do not include.
* config/get_dlcflags: Add -fno-common to DLCFLAGS on Darwin.
* config/get_config_options: Catch invalid arguments tp "--graphic".
* config/get_fltk; Add libstdc++; check presence of FLTK headers, too.
* config/get_X11: Search X11 library in */lib64/*, too (not just */lib/*).

Documentation:
* doc/gphelp.in: cwitty: Disable TeX; allow bz2 compression.

Header files:
* src/headers/paripriv.h: *After* building, on {OS X, SunOS, CYGWIN},
                          rename "ECHO".

C files:
* src/kernel/gmp/mp.c: Needed so that Sage can catch PARI's error signals.
                       Also allow disabling use of "GMP internals".
* src/language/init.c: Needed so that Sage can catch PARI's error signals.
* src/basemath/base2.c and src/basemath/polarit3.c: Fix PARI bug 1079.

======================================================================
Files previously patched:
======================================================================

* config/get_cc: On SunOS only, add "-fPIC" to compiler flags (David Kirkby)
  This happens not with a patch file, but with a sed command in
  spkg-install.  Now fixed upstream (-fPIC added on all platforms).

* config/get_dlcflags: mabshoff: To get around problem in PPC 32-bit
  Linux build.  Now fixed upstream (-fPIC added on all platforms).
  NB: We still patch this file, but for a slightly different reason.

* config/get_kernel: pjeremy: Fix for FreeBSD: #7825.  Supposedly fixed
  upstream.

* config/get_dlld: Undocumented patch for Darwin.  Removed to see what
  happens...

* config/get_cc: Changed OPTFLAG from "-O3" to "-O1" on Linux because of
  problems on Fedora 11 (32-bit) with one gcc version (ticket #7092).

  # Minh Van Nguyen: copy over patched get_cc (see ticket #7092). It's
  # reported that 32-bit Fedora 11 would fail to build otherwise.
  if [ `uname` = "Linux" ]; then
      cp "$TOP"/patches/get_cc config/get_cc
  fi

* src/headers/paridecl.h: Used to need a dummy variable changed from B0
  to N; now fixed upstream.
* src/headers/paripriv.h: Used to need a dummy variable changed from B0
  to N; now fixed upstream.  NB: There's another patch on this file still
  in place!
* config/Makefile.SH: Change "test -e" to "test -f" for Solaris.  Fixed
  upstream
