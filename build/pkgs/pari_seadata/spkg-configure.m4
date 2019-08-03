SAGE_SPKG_CONFIGURE([pari_seadata], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install pari_seadata as well])
        sage_spkg_install_pari_seadata=yes
    else
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([is pari_seadata installed? ])
        gp_seadat_check=`echo "poldegree(ellmodulareqn(211)[[1]])" | $GP -qf`
        if test x$gp_seadat_check = x212; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without seadata package])
            AC_MSG_NOTICE([Install seadata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_seadata=yes
            sage_spkg_install_pari=yes
        fi
    fi
])
