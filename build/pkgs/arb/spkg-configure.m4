SAGE_SPKG_CONFIGURE([arb], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_FLINT])
    SAGE_ARB_LIBRARY="arb"
    AC_MSG_CHECKING([installing flint? ])
    if test x$sage_spkg_install_flint = xyes; then
        AC_MSG_RESULT([yes; install arb as well])
        sage_spkg_install_arb=yes
    else
        AC_CHECK_HEADER(arb.h, [
        dnl below function added in version 2.16 of arb
            AC_SEARCH_LIBS([acb_mat_eig_simple], [arb], [
	      ],[dnl in Debian the name of dylib is different.
              AC_SEARCH_LIBS([acb_mat_eig_simple], [flint-arb], [
                  SAGE_ARB_LIBRARY="flint-arb"
	        ], [sage_spkg_install_arb=yes])
	    ])
        ], [sage_spkg_install_arb=yes])
    fi
    AC_SUBST(ARB_LIBRARY,[$SAGE_ARB_LIBRARY])
])
