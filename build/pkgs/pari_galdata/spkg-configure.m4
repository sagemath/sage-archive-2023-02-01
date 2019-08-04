SAGE_SPKG_CONFIGURE([pari_galdata], [
    AC_PATH_PROG([GP], [gp])
    if test x$GP != x; then
        AC_MSG_CHECKING([is pari_galdata installed? ])
        gp_gal_check=`echo "polgalois(x^8-2)[[1]]" | $GP -qf`
        if test x$gp_gal_check = x16; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without gapdata package])
            AC_MSG_NOTICE([Install galdata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_galdata=yes
        fi
    else
        AC_MSG_NOTICE([gp is not found])
        sage_spkg_install_pari_elldata=yes
    fi
])
