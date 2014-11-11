======================================================================
Current patches to PARI in Sage:
======================================================================

Patches to configuration files:
* get_ld.patch (Carl Witty): Disable -rpath.
* get_config_options.patch (Leif Leonhardy): Catch invalid arguments to
  "--graphic" (and treat such as an error) since otherwise strange
  compilation errors might occur (cf. #9722, too).
* get_fltk.patch (Leif Leonhardy): do an extra check for the FLTK
  include dir (cf. #9722).
* get_dlcflags.patch (Jeroen Demeyer): Add -fno-common to DLCFLAGS on
  Darwin. Submitted upstream, but upstream only applied it for PowerPC.
  Since this doesn't break anything and only improves performance, add
  the flag unconditionally.
* perl_path.patch (Jeroen Demeyer): change first line of all perl
  scripts to "#!/usr/bin/env perl" (#10559).
* KERNELCFLAGS.patch (Jeroen Demeyer): when SAGE_DEBUG=yes, compile
  kernel files with -O1 instead of -funroll-loops; -O0 gives a
  segmentation fault on some OS X systems when doing
  factor(10356613*10694706299664611221)
  See #13921, also reported upstream:
  - http://pari.math.u-bordeaux.fr/archives/pari-dev-1301/msg00000.html

C files:
* GCC_PR49330.patch (Jeroen Demeyer): in pari_init_functions(), reorder
  code to work around a GCC bug concerning pointer arithmetic:
  http://gcc.gnu.org/bugzilla/show_bug.cgi?id=49330
  This bug manifests itself as a Bus Error on OS X 10.4 PPC with
  gcc-4.6.3. Discussed with upstream:
  - http://pari.math.u-bordeaux.fr/archives/pari-dev-1203/msg00001.html
  - http://pari.math.u-bordeaux.fr/archives/pari-dev-1202/msg00006.html
* det_garbage.patch (Jeroen Demeyer, #15654): When computing a
  determinant(), only collect garbage once per outer loop iteration.
  Better increase PARI stack size instead of collecting garbage too
  often.
* nffactor.patch (Jeroen Demeyer, #16894): Fix an nffactor() bug, taken
  from PARI git commit 7630f584371c3db43c5fd6f57900c70a2c832b8e.
