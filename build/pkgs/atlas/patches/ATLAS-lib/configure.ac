AC_INIT([ATLAS],[3.10.1])
AM_INIT_AUTOMAKE(foreign)

LT_INIT

AC_PROG_CC
AC_PROG_F77

AC_CONFIG_FILES([Makefile])
AC_CONFIG_MACRO_DIR([m4])

AC_CHECK_LIB(pthread, pthread_create, [PTHREAD_LIB="-lpthread"])
AC_SUBST(PTHREAD_LIB)

AC_ARG_ENABLE([static],
    [AS_HELP_STRING([--enable-static],
                    [also install static libtool library (default: no)])])
AS_IF([test "x$enable_static" = "xyes"],
    [libtool_type=""],
    [libtool_type=-shared])
AC_SUBST(LIBTOOL_TYPE, [$libtool_type])

AC_OUTPUT
