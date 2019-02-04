SAGE_SPKG_CONFIGURE([gmp], [
    AC_CHECK_HEADER(gmp.h, [], [sage_spkg_install_gmp=yes])
    AC_CHECK_HEADER(gmpxx.h, [], [sage_spkg_install_gmp=yes])
dnl mpn_divexact_1 appeared in GMP 6.1.0
    AC_SEARCH_LIBS([__gmpn_divexact_1], [gmp], [break], [sage_spkg_install_gmp=yes])
])
