Details of which source files are patched in spkg-install, and why.
See below for information on files which used to be patched but are no
longer, or are still but now for a different reason. (Do not delete that!)

======================================================================
Current patches to PARI in Sage:
======================================================================

Patches to configuration files:
* get_ld.patch: cwitty: Disable -rpath.
* get_tests.patch: John Cremona: Disable testing of ellglobalred in
  "make test-all" in spkg-check, since it requires the elldata database
  which we do not include since it is about 18MB.
* get_config_options.patch: leif: Catch invalid arguments to "--graphic"
  (and treat such as an error) since otherwise strange compilation errors
  might occur (cf. #9722, too).
* get_fltk.patch: leif: Explicitly link against libstdc++ when using FLTK
  (for plotting) to support Fedora 13, and do an extra check for the FLTK
  include dir (cf. #9722).
* get_X11.patch: leif: Also search */lib64/* for X11 libraries (since on
  some systems this is really a separate directory, i.e. neither a
  symbolic link to */lib/* nor the target of a symbolic link */lib/*; cf.
  #9722, too).
* get_dlcflags.patch: jdemeyer: Add -fno-common to DLCFLAGS on Darwin.
  Fixed upstream, but only for PowerPC. Since this doesn't break anything
  and only improves performance, add the flag unconditionally.
* install_doc_no_make.patch: jdemeyer: Do not *build* the documentation
  when doing install-doc or install-docpdf.  We must not build the
  documentation because that requires tex.  On the other hand, to have ?
  and ?? work within gp, we must install the .tex files (but not .dvi
  files).  So simply not doing install-doc doesn't work.
* perl_path.patch: jdemeyer: change first line of all perl scripts
  to "#!/usr/bin/env perl" (#10559).  Note that this patch will always
  apply with fuzz 2 because of the svn Id.

C files:
* src/kernel/gmp/mp.c:
  Do not override GMP's memory functions.
  In addition, let PARI use "GMP internals" (access members of GMP
  structures directly) *conditionally*. (We *don't* disable that by
  default, since at least currently this is compatible with both GMP
  *and* MPIR. To disable the use, add "-DPARI_DONT_USE_GMP_INTERNALS"
  to CFLAGS. This is a trivial patch to mp.c only; see also the comment
  there.):
* reorder_init_opts.patch: call pari_init_defaults() *before* calling
  gp_expand_path(), such that warnings during gp_expand_path() do not
  cause a dereference of the NULL pointer pariErr (#12158).
  Reported upstream at
  http://pari.math.u-bordeaux.fr/cgi-bin/bugreport.cgi?bug=1264

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
* Configure: First line was changed from "#! /bin/sh" to
  "#!/usr/bin/env bash".  Now we directly call `bash Configure` instead
  in spkg-install.
* src/basemath/base2.c and src/basemath/polarit3.c: Fix PARI bug 1079.
  Now fixed upstream.
* doc/gphelp.in: cwitty: allow bz2 compression.  Now upstream.
* src/headers/paripriv.h: *After* building, on {OS X, SunOS, CYGWIN},
  rename "ECHO".  Upstream changed ECHO to gpd_ECHO, so the patch is
  not needed anymore.
* Makefile_mv.patch: Fix race condition in parallel "make install",
  see http://pari.math.u-bordeaux.fr/cgi-bin/bugreport.cgi?bug=1148
  Removed because race conditions still remain, instead make -j1 install.
* pari_1084.patch: exotic branch cut convention (#9620).  jdemeyer:
  removed the parts of the patch involving the COMPAT and CHANGES files
  as they are documentation only and give patch conflicts.  Rebased
  a hunk in src/basemath/trans1.c to make the patch apply cleanly.
* pari_1132.patch: nffactor() returns reducible factor (#10279)
* pari_1141.patch: factoring non-square-free polynomial over number
  fields (#10369)
* pari_1143.patch: rnfisnorm failing for non-integral elements (#2329)
* pari_1144.patch: rnfisnorminit requires leading coefficient 1 to be
  in ZZ instead of the base field (#2329)
* src/language/init.c:
  These was needed so that Sage can catch PARI's error signals.
  Turns out this is totally not needed since we use err_catch()
  in devel/sage/sage/libs/pari/pari_err.h
* osx_13318_13330.patch: fix linking on certain older Mac OS X systems.
  Upstreamed in PARI 2.5.1.
