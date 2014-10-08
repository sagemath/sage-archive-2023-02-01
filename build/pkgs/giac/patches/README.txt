*   the patch 1.1.1._index.h.patch has been merged in giac 1.1.2 (2014 Aug 01 version)

*   the nofltk is because giac is compiled without fltk, so some outputs differ
    in the make check

*   Up to giac 1.1.1 the giac library initializes the PARI library
    and Sage also does the same thing, leading to conflicts with the cython
    interface giacpy.
    The patch 1.1.1.pari_init_opt.patch avoids this extra initialization.
    A typical symptom of this conflict was control-c after loading giacpy quits sage
    this patch has been merged in giac 1.1.2 (2014 Oct 8)
