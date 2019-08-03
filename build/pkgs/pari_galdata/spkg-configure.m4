SAGE_SPKG_CONFIGURE([pari_galdata], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install pari_galdata as well])
        sage_spkg_install_pari_galdata=yes
    else
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([is pari_galdata installed? ])
        gp_gal_check=`echo "polgalois(x^8-2)[[1]]" | $GP -qf`
        if test x$gp_gal_check = x16; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without gapdata package])
            AC_MSG_NOTICE([Install galdata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_galdata=yes
            sage_spkg_install_pari=yes
        fi
    fi
])
