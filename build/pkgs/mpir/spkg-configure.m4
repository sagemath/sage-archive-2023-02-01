SAGE_SPKG_CONFIGURE([mpir], [
    AC_ARG_WITH([mp],
    [AS_HELP_STRING([--with-mp=system],
        [use the system GMP as multiprecision library, if possible (default)])]
    [AS_HELP_STRING([--with-mp=mpir],
        [use the Sage SPKG for MPIR as multiprecision library])]
    [AS_HELP_STRING([--with-mp=gmp],
        [use the Sage SPKG for GMP as multiprecision library])])

dnl Just part the options here
    case "$with_mp" in
        system) ;;
        MPIR|mpir) with_mp=mpir;;
        GMP|gmp) with_mp=gmp;;
        "") with_mp=system;;
        *)
            AC_MSG_ERROR([allowed values for --with-mp are system, mpir, or gmp]);;
    esac

dnl Implement cases for what to do on different options here
    case "$with_mp" in
        system)
            AC_CHECK_HEADER(gmp.h, [], [sage_spkg_install_mpir=yes])
            AC_CHECK_HEADER(gmpxx.h, [], [sage_spkg_install_mpir=yes])
            dnl mpq_cmp_z appeared in GMP 6.1.0 and is used by pynac
            AC_SEARCH_LIBS([__gmpq_cmp_z], [gmp], [break],
                [sage_spkg_install_mpir=yes])
            SAGE_MP_LIBRARY=mpir
            ;;
        mpir)
            sage_spkg_install_mpir=yes
            SAGE_MP_LIBRARY=mpir
            ;;
        gmp)
            sage_spkg_install_gmp=yes
            SAGE_MP_LIBRARY=gmp
            ;;
    esac

    if test x$sage_spkg_install_mpir = xyes -o x$sage_spkg_install_gmp = xyes; then
        AC_SUBST(SAGE_GMP_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using $SAGE_MP_LIBRARY SPKG (via --with-mp=$SAGE_MP_LIBRARY)])
    else
        AC_SUBST(SAGE_GMP_PREFIX, [''])
        AC_MSG_RESULT([using GMP-compatible library from the system])
    fi

    AC_SUBST([SAGE_MP_LIBRARY], [$SAGE_MP_LIBRARY])
])

