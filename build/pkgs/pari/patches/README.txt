======================================================================
Current patches to PARI in Sage:
======================================================================

Patches to configuration files:
* get_ld.patch (Carl Witty): Disable -rpath.
* get_config_options.patch (Leif Leonhardy): Catch invalid arguments to
  "--graphic" (and treat such as an error) since otherwise strange
  compilation errors might occur (cf. #9722, too).
* get_dlcflags.patch (Jeroen Demeyer): Add -fno-common to DLCFLAGS on
  Darwin. Submitted upstream, but upstream only applied it for PowerPC.
  Since this doesn't break anything and only improves performance, add
  the flag unconditionally.

C files:
* public_memory_functions.patch (Jeroen Demeyer, #16997): Make some of
  PARI's private memory functions public to improve interface with Sage.
