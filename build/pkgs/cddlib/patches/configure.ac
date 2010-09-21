dnl Process this file with autoconf to produce a configure script.
AC_PREREQ([2.63])
AC_INIT([cddlib],[0.94f])
AM_INIT_AUTOMAKE([1.11])
AC_CONFIG_MACRO_DIR([m4])

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
LT_INIT

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

dnl Checks for functions
AC_FUNC_MALLOC
AC_CHECK_FUNCS([strchr])


dnl Find the GNU Multiple Precision Arithmetic Library
AC_ARG_WITH(gmp,
  [AS_HELP_STRING([--with-gmp=<path>],
     	          [prefix of GMP installation. e.g. /usr/local or /usr])],
  [GMP_PREFIX="$with_gmp"],
  [AC_CHECK_LIB([gmp], [__gmpz_init], ,
    AC_MSG_ERROR([GNU Multiple Precision Arithmetic Library not found.
    Maybe you want to call configure with the --with-gmp=<path> option?
    This tells configure where to find the GMP library and headers.
    e.g. --with-gmp=/usr/local or --with-gmp=/usr]))]
)
if test "x$GMP_PREFIX" != x ; then
  GMP_LIBS=" -L${GMP_PREFIX}/lib"
  GMP_CFLAGS=" -I${GMP_PREFIX}/include"
else
  GMP_LIBS=""
  GMP_CFLAGS=""
fi
GMP_LIBS="${GMP_LIBS} -lgmp"
AC_SUBST(GMP_LIBS)
AC_SUBST(GMP_CFLAGS)


AC_CONFIG_FILES([
  lib-src/Makefile
  src/Makefile
  lib-src-gmp/Makefile
  src-gmp/Makefile
  Makefile
])
AC_OUTPUT
