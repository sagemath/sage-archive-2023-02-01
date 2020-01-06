SAGE_SPKG_CONFIGURE([bzip2], [
    AC_CHECK_HEADER(bzlib.h, [
      AC_SEARCH_LIBS([BZ2_bzCompress], [bz2], [
        AC_PATH_PROG([bzip2prog], [bzip2], [sage_spkg_install_bzip2=yes])
      ], [sage_spkg_install_bzip2=yes])
    ], [sage_spkg_install_bzip2=yes])
])
