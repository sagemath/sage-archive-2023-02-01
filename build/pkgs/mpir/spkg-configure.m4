SAGE_SPKG_CONFIGURE([mpir], [
    AC_ARG_WITH([mp],
    [AS_HELP_STRING([--with-mp=system],
        [use the system GMP as multiprecision library, if possible (default)])]
    [AS_HELP_STRING([--with-mp=mpir],
        [use the Sage SPKG for MPIR as multiprecision library])]
    [AS_HELP_STRING([--with-mp=gmp],
        [use the Sage SPKG for GMP as multiprecision library])])

    if test -n "$with_mp"; then
        SAGE_MP_LIBRARY="$with_mp"
    fi

dnl Just part the options here
    case "$SAGE_MP_LIBRARY" in
        system) ;;
        MPIR|mpir) SAGE_MP_LIBRARY=mpir;;
        GMP|gmp) SAGE_MP_LIBRARY=gmp;;
        "") SAGE_MP_LIBRARY=system;;
        *)
            AC_MSG_ERROR([allowed values for --with-mp are system, mpir, or gmp]);;
    esac

dnl Implement cases for what to do on different options here
    case "$SAGE_MP_LIBRARY" in
        system)
            AC_CHECK_HEADER(gmp.h, [], [sage_spkg_install_mpir=yes])
            AC_CHECK_HEADER(gmpxx.h, [], [sage_spkg_install_mpir=yes])
            dnl mpq_cmp_z appeared in GMP 6.1.0 and is used by pynac
            AC_SEARCH_LIBS([__gmpq_cmp_z], [gmp], [break],
                [sage_spkg_install_mpir=yes])
            AC_MSG_RESULT([using GMP-compatible library from the system])
            ;;
        mpir)
            sage_spkg_install_mpir=yes
            SAGE_MP_LIBRARY=mpir
            AC_MSG_RESULT([using mpir SPKG (via --with-mp=mpir)])
            ;;
        gmp)
            SAGE_MP_LIBRARY=gmp
            sage_spkg_install_gmp=yes
            AC_MSG_RESULT([using gmp SPKG (via --with-mp=gmp)])
            ;;
    esac

    if test $SAGE_MP_LIBRARY != system; then
        AC_SUBST(SAGE_GMP_PREFIX, ['$SAGE_LOCAL'])
    else
        AC_SUBST(SAGE_GMP_PREFIX, [''])
    fi

    AC_SUBST([SAGE_MP_LIBRARY], [$SAGE_MP_LIBRARY])
])

