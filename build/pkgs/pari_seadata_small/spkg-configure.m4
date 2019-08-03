SAGE_SPKG_CONFIGURE([pari_seadata_small], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install pari_seadata_small as well])
        sage_spkg_install_pari_seadata_small=yes
    else
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([is pari_seadata_small installed? ])
        gp_seadat_check=`echo "poldegree(ellmodulareqn(11)[[1]])" | $GP -qf`
        if test x$gp_seadat_check = x12; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without seadata package])
            AC_MSG_NOTICE([Install seadata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_seadata_small=yes
            sage_spkg_install_pari=yes
        fi
    fi
])
