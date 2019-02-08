SAGE_SPKG_CONFIGURE([gf2x], [
    AC_CHECK_HEADER(gf2x.h, [], [sage_spkg_install_gf2x=yes])
dnl gf2x_mul_r appeared in version 1.2 of GF2X 
    AC_SEARCH_LIBS([gf2x_mul_r], [gf2x], [break], [sage_spkg_install_gf2x=yes])
])
