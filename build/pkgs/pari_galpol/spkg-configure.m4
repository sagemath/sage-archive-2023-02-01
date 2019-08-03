SAGE_SPKG_CONFIGURE([pari_galpol], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install pari_galpol as well])
        sage_spkg_install_pari_galpol=yes
    else
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([is pari_galpol installed? ])
        gp_ell_check=`echo "galoisgetname(12,1)" | $GP -qf`
        if test "x$gp_ell_check = xC3\ \:\ C4"; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without galpol package])
            AC_MSG_NOTICE([Install galpol package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_galpol=yes
            sage_spkg_install_pari=yes
        fi
    fi
])
