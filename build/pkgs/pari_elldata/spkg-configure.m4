SAGE_SPKG_CONFIGURE([pari_elldata], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install pari_elldata as well])
        sage_spkg_install_pari_elldata=yes
    else
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([is pari_elldata installed? ])
        gp_ell_check=`echo "r=ellinit(\"11a1\"); r[[11]]" | $GP -qf`
        if test x$gp_ell_check = x20008; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without elldata package])
            AC_MSG_NOTICE([Install elldata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari_elldata=yes
            sage_spkg_install_pari=yes
        fi 
    fi
])
