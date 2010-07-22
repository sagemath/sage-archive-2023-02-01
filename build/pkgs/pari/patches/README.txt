Details of which source files are patched in spkg-install, and why.
See below for information on files which used to be patched but are no
longer (do not delete that!)

======================================================================
Files patched as of pari-2.4.3.svn.p6:
======================================================================

* Configuration files:

config/get_cc: on SunOS only, add "-fPIC" to compiler flags (David Kirkby)
    This happens not with a patch file, but with a sed command in
    spkg-install.

config/get_dlld: on all platforms, but affects Darwin only

config/get_ld: cwitty: disable -rpath

config/get_dlcflags: mabshoff: To get around problem in PPC 32-bit Linux build

config/get_kernel: pjeremy: fix for FreeBSD: #7825

config/get_tests: John Cremona: disable testing of ellglobalred in
                  "make test-all" in spkg-check, since it requires the
                  elldata database which we do not include.

* Documentation:

doc/gphelp.in: cwitty: disable TeX; allow bz2 compression

* Header files:

src/headers/paripriv.h: *after* building, on {OS X, SunOS, SYGWIN},
                        rename "ECHO"

* C files:
src/kernel/gmp/mp.c: needed so that Sage can catch pari's error signals.
src/language/init.c: needed so that Sage can catch pari's error signals.

======================================================================
Files previously patched:
======================================================================

config/get_cc: changed OPTFLAG from "-O3" to "-)1" on Linux because of
problems on Fedora 11 (32-bit) with one gcc version (ticket #7092).

# Minh Van Nguyen: copy over patched get_cc (see ticket #7092). It's
# reported that 32-bit Fedora 11 would fail to build otherwise.
if [ `uname` = "Linux" ]; then
    cp "$TOP"/patches/get_cc config/get_cc
fi

src/headers/paridecl.h: used to need a dummy variable changed from B0
to N; now fixed upstream.
src/headers/paripriv.h: used to need a dummy variable changed from B0
to N; now fixed upstream.  NB There's another patch on this file still
in place!

