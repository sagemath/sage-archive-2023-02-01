SAGE_SPKG_CONFIGURE([pari_seadata_small], [
    AC_PATH_PROG([GP], [gp])
    if test x$GP != x; then
        AC_MSG_CHECKING([is pari_seadata_small installed? ])
        gp_seadat_check=`echo "poldegree(ellmodulareqn(11)[[1]])" | $GP -qf`
        if test x$gp_seadat_check = x12; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without seadata package])
            AC_MSG_NOTICE([Install seadata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_seadata_small=yes
        fi
    else
        AC_MSG_NOTICE([gp is not found])
        sage_spkg_install_pari_seadata_small=yes
    fi
])
