dnl Process this file with autoconf to produce a configure script.
AC_INIT
AC_CONFIG_MACRO_DIR([m4])

dnl Init automake.
AM_INIT_AUTOMAKE(cddlib, 0.94)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
LT_INIT

dnl Checks for libraries.
dnl Replace `main' with a function in -lg:
dnl AC_CHECK_LIB(g, main)
dnl Replace `main' with a function in -lstdc:
dnl AC_CHECK_LIB(stdc, main)
dnl look for gmp:
AC_CHECK_LIB(gmp, main)

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

dnl Checks for library functions.

AC_CONFIG_FILES([lib-src/Makefile src/Makefile lib-src-gmp/Makefile src-gmp/Makefile Makefile])
AC_OUTPUT
