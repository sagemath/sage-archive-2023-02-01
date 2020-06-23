SAGE_SPKG_CONFIGURE([gfortran], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GCC])

    AC_REQUIRE([AC_PROG_FC])
    AC_SUBST(FC)

    # Check that the Fortran compiler accepts free-format source code (as
    # opposed to the older fixed-format style from Fortran 77).
    # This helps verify the compiler works too, so if some idiot sets FC to
    # /usr/bin/ls, we will at least know it's not a working Fortran
    # compiler.
    AC_FC_FREEFORM([SAGE_HAVE_FC_FREEFORM=yes], [
	AC_MSG_NOTICE([Your Fortran compiler does not accept free-format source code])
        AC_MSG_NOTICE([which means the compiler is either seriously broken, or])
        AC_MSG_NOTICE([is too old to build Sage.])
        SAGE_HAVE_FC_FREEFORM=no])

    AS_VAR_IF(SAGE_HAVE_FC_FREEFORM, [no], [
        AS_VAR_SET(sage_spkg_install_gfortran, [yes])
    ])

    # Special case: If we are already installing gcc then don't install
    # gfortran since it's included
    if test "x$sage_spkg_install_gcc" = "xyes" -o x$SAGE_INSTALL_GCC = xexists; then
        sage_spkg_install_gfortran=no
    fi

    # Copy CFLAGS to FCFLAGS if this is not set.
    if test "x$FCFLAGS" = "x"; then
        AC_SUBST(FCFLAGS, "$CFLAGS")
        AC_SUBST(FCFLAGS_O3, "$CFLAGS_O3")
        AC_SUBST(FCFLAGS_NON_NATIVE, "$CFLAGS_NON_NATIVE")
        AC_SUBST(FCFLAGS_O3_NON_NATIVE, "$CFLAGS_O3_NON_NATIVE")
    else
        AC_SUBST(FCFLAGS_03, "$FCFLAGS")
        AC_SUBST(FCFLAGS_NON_NATIVE, "$FCFLAGS")
        AC_SUBST(FCFLAGS_O3_NON_NATIVE, "$FCFLAGS")
    fi

    # Copy FCFLAGS to F77FLAGS if this is not set.
    if test "x$F77FLAGS" = "x"; then
        AC_SUBST(F77FLAGS, "$FCFLAGS")
        AC_SUBST(F77FLAGS_O3, "$FCFLAGS_O3")
        AC_SUBST(F77FLAGS_NON_NATIVE, "$FCFLAGS_NON_NATIVE")
        AC_SUBST(F77FLAGS_O3_NON_NATIVE, "$FCFLAGS_O3_NON_NATIVE")
    else
        AC_SUBST(F77FLAGS_03, "$F77FLAGS")
        AC_SUBST(F77FLAGS_NON_NATIVE, "$F77FLAGS")
        AC_SUBST(F77FLAGS_O3_NON_NATIVE, "$F77FLAGS")
    fi

])
