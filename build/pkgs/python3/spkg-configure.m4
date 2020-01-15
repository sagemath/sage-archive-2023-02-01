SAGE_SPKG_CONFIGURE(
    [python3], [
        AC_CACHE_CHECK([for python3 >= 3.7.3, < 3.8], [ac_cv_path_PYTHON3], [
            AC_PATH_PROGS_FEATURE_CHECK([PYTHON3], [python3], [
                python3_version=`$ac_path_PYTHON3 --version 2>&1 \
                    | $SED -n -e 's/\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\).*/\1/p'`
                AS_IF([test -n "$python3_version"], [
                    AX_COMPARE_VERSION([$python3_version], [ge], [3.7.3], [
                        AX_COMPARE_VERSION([$python3_version], [lt], [3.8.0], [
                            ac_cv_path_PYTHON3="$ac_path_PYTHON3"
                        ])
                    ])
                ])
            ])
        ])
        AS_IF([test -z "$ac_cv_path_PYTHON3"], [sage_spkg_install_python3_build=yes])
    ]
)
