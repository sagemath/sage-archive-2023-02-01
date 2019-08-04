SAGE_SPKG_CONFIGURE([pari_galpol], [
    AC_PATH_PROG([GP], [gp])
    if test x$GP != x; then
        AC_MSG_CHECKING([is pari_galpol installed? ])
        gp_ell_check=`echo "galoisgetname(12,1)" | $GP -qf`
        if test "x$gp_ell_check = xC3\ \:\ C4"; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without galpol package])
            AC_MSG_NOTICE([Install galpol package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_galpol=yes
        fi
    else
        AC_MSG_NOTICE([gp is not found])
        sage_spkg_install_pari_galpol=yes
    fi
])
