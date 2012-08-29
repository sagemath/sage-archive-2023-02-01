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
  to "#!/usr/bin/env perl" (#10559).

C files:
* src/kernel/gmp/mp.c:
  Do not override GMP's memory functions.
  In addition, let PARI use "GMP internals" (access members of GMP
  structures directly) *conditionally*. (We *don't* disable that by
  default, since at least currently this is compatible with both GMP
  *and* MPIR. To disable the use, add "-DPARI_DONT_USE_GMP_INTERNALS"
  to CFLAGS. This is a trivial patch to mp.c only; see also the comment
  there.):
* GCC_PR49330.patch: in pari_init_functions(), reorder code to work
  around a GCC bug concerning pointer arithmetic:
  http://gcc.gnu.org/bugzilla/show_bug.cgi?id=49330
  This bug manifests itself as a Bus Error on OS X 10.4 PPC with
  gcc-4.6.3.
* rootpol.patch: fixes a Segmentation Fault in Sage "roots" (#13314)
