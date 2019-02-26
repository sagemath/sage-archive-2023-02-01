SAGE_SPKG_CONFIGURE([arb], [
    AC_CHECK_HEADER(arb.h, [
dnl below function added in version 2.16 of arb 
       AC_SEARCH_LIBS([acb_mat_eig_simple], [arb], [break], [sage_spkg_install_arb=yes])
    ], [sage_spkg_install_arb=yes])
])
