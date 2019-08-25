SAGE_SPKG_CONFIGURE([lcalc], [
    m4_pushdef([SAGE_LCALC_MINVER],["1.22"])
    SAGE_SPKG_DEPCHECK([pari mpfr], [
        AC_PATH_PROG([LCALC], [lcalc])
        AS_IF([test x$LCALC = x], [
           AC_MSG_NOTICE([lcalc not found. Installing lcalc])
           sage_spkg_install_lcalc=yes], [
           AC_MSG_CHECKING([is lcalc's version good enough? ])
           lcalc_ver=`$LCALC --version 2>>/dev/null | $SED -e 's/lcalc\ //' | $SED -e 's/\ .*//g'`
           AX_COMPARE_VERSION([$lcalc_ver], [ge], [$SAGE_LCALC_MINVER], [
               AC_MSG_RESULT([yes.])
               AC_CHECK_HEADERS([Lfunction/L.h libLfunction/L.h],
                  [AS_IF([test x$ac_cv_header_libLfunction_L_h = x],
                         [SAGE_LCALC_INCDIR_NOLIBPREFIX='yes'],
                         [SAGE_LCALC_INCDIR_NOLIBPREFIX=''])
                   sage_spkg_install_lcalc=no
                   break;],
                  [sage_spkg_install_lcalc=yes])
          AC_MSG_CHECKING([whether we can link and run a program using libLfunction])
          LCALC_SAVED_LIBS=$LIBS
          LIBS="$LIBS -lLfunction"
          AC_RUN_IFELSE([
            AC_LANG_PROGRAM([[#ifdef HAVE_LFUNCTION_L_H
                                 #include <Lfunction/L.h>
                              #else
                                 #include <libLfunction/L.h>
                              #endif]],
                      [[initialize_globals();
                        Complex x;
                        x = Pi*I;
                        L_function<int> L4;
                        return 0;]]
            )], [AC_MSG_RESULT([yes; use lcalc from the system])], [
            AC_MSG_RESULT([no; install lcalc])
            sage_spkg_install_lcalc=yes
            LIBS=$LCALC_SAVED_LIBS
          ])
               ], [
               AC_MSG_RESULT([no. Install lcalc])
               sage_spkg_install_lcalc=yes])
        ])
    ])
    m4_popdef([SAGE_LCALC_MINVER])
], [], [], [
    AS_IF([test x$sage_spkg_install_lcalc = xyes],
          [AC_SUBST(SAGE_LCALC_INCDIR_NOLIBPREFIX,[''])],
          [AC_SUBST(SAGE_LCALC_INCDIR_NOLIBPREFIX,[$SAGE_LCALC_INCDIR_NOLIBPREFIX])])
])
