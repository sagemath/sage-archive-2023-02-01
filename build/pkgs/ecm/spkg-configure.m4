SAGE_SPKG_CONFIGURE([ecm], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GMP])
    AC_MSG_CHECKING([installing gmp/mpir? ])
    if test x$sage_spkg_install_mpir = xyes -o x$sage_spkg_install_gmp = xyes; then
        AC_MSG_RESULT([yes; install ecm as well])
        sage_spkg_install_ecm=yes
    else

        AC_CHECK_HEADER(ecm.h, [
            AX_ABSOLUTE_HEADER([ecm.h])
            if test x$gl_cv_absolute_ecm_h = x; then
                AC_MSG_ERROR(m4_normalize([
                    failed to find absolute path to ecm.h despite it being reported found
                ]))
                sage_spkg_install_ecm=yes
            else
    dnl check that the version is at least 7.0.4
                ecm_version=`grep ECM_VERSION $gl_cv_absolute_ecm_h |
                  $SED -n -e 's/\#define ECM_VERSION "*\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\)"/\1/p'`
                AS_IF([test -n "$ecm_version"], [
                    AX_COMPARE_VERSION([$ecm_version], [ge], [7.0.4], [
                        ac_cv_ECM="$ecm_version"
                        AC_SEARCH_LIBS([ecm_factor], [ecm], [], [sage_spkg_install_ecm=yes])
                    ])
                ])
            fi
            AS_IF([test -z "$ac_cv_ECM"], [sage_spkg_install_ecm=yes])
        ], [sage_spkg_install_ecm=yes])
    fi
])
