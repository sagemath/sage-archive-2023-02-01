SAGE_SPKG_CONFIGURE([python3], [
   SAGE_SPKG_DEPCHECK([sqlite libpng bzip2 xz libffi], [
      AS_IF([test $SAGE_PYTHON_VERSION = 2], [
        dnl If we use Python 2 for Sage, we install Python 3 too and do NOT attempt to do
        dnl venv using system python3 over SAGE_LOCAL.
        dnl (In particular, the setuptools and pip install scripts are not prepared for
        dnl handling this situation.)
        sage_spkg_install_python3=yes
      ], [
        dnl Using Python 3 for Sage.  Check if we can do venv with a system python3
        dnl instead of building our own copy.
        check_modules="sqlite3, ctypes, math, hashlib, crypt, readline, socket, zlib, distutils.core"
        m4_pushdef([MIN_VERSION], [3.7.0])
        m4_pushdef([LT_VERSION],  [3.9.0])
        AC_CACHE_CHECK([for python3 >= ]MIN_VERSION[, < ]LT_VERSION[ with modules $check_modules], [ac_cv_path_PYTHON3], [
            AC_MSG_RESULT([])
            AC_PATH_PROGS_FEATURE_CHECK([PYTHON3], [python3.8 python3.7 python3], [

                SAGE_CHECK_PYTHON_FOR_VENV([$ac_path_PYTHON3],
                                           MIN_VERSION, LT_VERSION,
                                           $check_modules, [
                    dnl It is good
                    ac_cv_path_PYTHON3="$ac_path_PYTHON3"
                    ac_path_PYTHON3_found=:
                    AC_MSG_RESULT([yes])
                    dnl introduction for AC_MSG_RESULT printed by AC_CACHE_CHECK
                    AC_MSG_CHECKING([for python3 >= ]MIN_VERSION[, < ]LT_VERSION[ with modules $check_modules])
                ])
            ])
        ])
        m4_popdef([MIN_VERSION])
        m4_popdef([LT_VERSION])
        AS_IF([test -z "$ac_cv_path_PYTHON3"],
              [sage_spkg_install_python3=yes])
      ])
   ])
],, [
    dnl PRE
], [
    dnl POST
    AS_IF([test x$sage_spkg_install_python3 = xno], [PYTHON_FOR_VENV="$ac_cv_path_PYTHON3"])
    AC_SUBST([PYTHON_FOR_VENV])

    dnl These temporary directories are created by the check above
    dnl and need to be cleaned up to prevent the "rm -f conftest*"
    dnl (that a bunch of other checks do) from emitting warnings about
    dnl conftest.dir and conftest_venv being directories.
    rm -rf conftest.dir conftest_venv
])
