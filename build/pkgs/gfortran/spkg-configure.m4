SAGE_SPKG_CONFIGURE([gfortran], [
	AC_REQUIRE([AC_PROG_FC])

    # Check that the Fortran compiler accepts free-format source code
    # (as opposed to the older fixed-format style from Fortran 77).
    # This helps verify the compiler works too, so if some idiot
    # sets FC to /usr/bin/ls, we will at least know it's
    # not a working Fortran compiler.
    if test -z "$FC"; then
        sage_spkg_install_gfortran=yes
        SAGE_MUST_INSTALL_GCC([a Fortran compiler is missing])
    fi

    # see http://www.gnu.org/software/hello/manual/autoconf/Fortran-Compiler.html
    AC_FC_FREEFORM([], [
        AC_MSG_NOTICE([Your Fortran compiler does not accept free-format source code])
        AC_MSG_NOTICE([which means the compiler is either seriously broken, or])
        AC_MSG_NOTICE([is too old to build Sage.])
        sage_spkg_install_gfortran=yes
    ])

    # Check compiler versions
    if test x$GFC != xyes; then
        sage_spkg_install_gfortran=yes
    fi
])
