SAGE_SPKG_CONFIGURE([pari_elldata], [
    AC_PATH_PROG([GP], [gp])
    if test x$GP != x; then
        AC_MSG_CHECKING([is pari_elldata installed? ])
        gp_ell_check=`echo "r=ellinit(\"11a1\"); r[[11]]" | $GP -qf`
        if test x$gp_ell_check = x20008; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without elldata package])
            AC_MSG_NOTICE([Install elldata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_elldata=yes
        fi
    else
        AC_MSG_NOTICE([gp is not found])
        sage_spkg_install_pari_elldata=yes
    fi
])
