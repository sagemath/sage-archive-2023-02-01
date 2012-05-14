This directory contains script that are run right at the beginning
of the SAGE source install.   The bzip2 source code is also included
here and built right off, since it is used for everything else.

NOTE:
* BZIP2: I modified bzip2 so it would build libbz2.a correctly.
  This involved adding -fPIC to CFLAGS in the Makefile.
  I made the contents of bzip2-1.0.2/words0 empty, since the
  useless scary message appears right at the beginning of
  SAGE compilation, and would likely confuse users.

