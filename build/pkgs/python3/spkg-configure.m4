SAGE_SPKG_CONFIGURE([python3], [
    SAGE_SPKG_DEPCHECK([sqlite libpng bzip2 xz libffi], [
        AC_CACHE_CHECK([for python3 >= 3.7.3, < 3.8 with sqlite3 module], [ac_cv_path_PYTHON3], [
            AC_MSG_RESULT([])
            AC_PATH_PROGS_FEATURE_CHECK([PYTHON3], [python3.7 python3], [
                AC_MSG_CHECKING([... whether $ac_path_PYTHON3 is good])
                python3_version=`"$ac_path_PYTHON3" --version 2>&1 \
                    | $SED -n -e 's/\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\).*/\1/p'`
                AS_IF([test -n "$python3_version"], [
                    AX_COMPARE_VERSION([$python3_version], [ge], [3.7.3], [
                        AX_COMPARE_VERSION([$python3_version], [lt], [3.8.0], [
                            dnl Because the system python is not used directly but rather in a venv without site-packages,
                            dnl we test whether the module will be available in a venv.
                            dnl Otherwise, some system site-package may be providing this module to the system python.
                            rm -rf config_venv
                            AS_IF(["$ac_path_PYTHON3" -m venv --without-pip --symlinks config_venv && config_venv/bin/python3 -c "import sqlite3"], [
                                ac_cv_path_PYTHON3="$ac_path_PYTHON3"
                                ac_path_PYTHON3_found=:
                                AC_MSG_RESULT([yes])
                                dnl introduction for AC_MSG_RESULT printed by AC_CACHE_CHECK
                                AC_MSG_CHECKING([for python3 >= 3.7.3, < 3.8 with sqlite3 module])
                            ], [
                                AC_MSG_RESULT([no, the version is in the supported range but cannot import sqlite3 in a venv])
                            ])
                        ], [
                            AC_MSG_RESULT([no, $python3_version is too recent])
                        ])
                    ], [
                        AC_MSG_RESULT([no, $python3_version is too old])
                    ])
                ], [
                    AC_MSG_RESULT([no, "$ac_path_PYTHON3 --version" does not work])
                ])
            ])
        ])
        AS_IF([test -z "$ac_cv_path_PYTHON3"],
              [sage_spkg_install_python3=yes])
    ])
],, [
    dnl PRE
], [
    dnl POST
    AS_IF([test "$sage_spkg_install_python3" = yes],
          [PYTHON_FOR_VENV=python3],
          [PYTHON_FOR_VENV="$ac_cv_path_PYTHON3"])
    AC_SUBST([PYTHON_FOR_VENV])
])
