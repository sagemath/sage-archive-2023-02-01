SAGE_SPKG_CONFIGURE([polymake], [
    dnl Need polymake >= 3.5 for "improve callable library compatibility with threads"
    m4_pushdef([POLYMAKE_VERSION_MIN], [3.5])
    AC_CACHE_CHECK([for polymake-config >= ]POLYMAKE_VERSION_MIN, [ac_cv_path_POLYMAKE_CONFIG], [
        AC_PATH_PROGS_FEATURE_CHECK([POLYMAKE_CONFIG], [polymake-config], [
            polymake_config_version=$($ac_path_POLYMAKE_CONFIG --version 2>&1)
            AS_IF([test -n "$polymake_config_version"], [
                AX_COMPARE_VERSION([$polymake_config_version], [ge], POLYMAKE_VERSION_MIN, [
                    ac_cv_path_POLYMAKE_CONFIG="$ac_path_POLYMAKE_CONFIG"
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_POLYMAKE_CONFIG"],
          [sage_spkg_install_polymake=yes],
          [sage_require_perl_cpan_polymake_prereq=no
           sage_require_perl_term_readline_gnu=no])
    m4_popdef([POLYMAKE_VERSION_MIN])
])
