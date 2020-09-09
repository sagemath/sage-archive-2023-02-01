SAGE_SPKG_CONFIGURE([tox], [
       dnl Accept a version of tox that can auto-provision
       dnl https://tox.readthedocs.io/en/latest/changelog.html#v3-8-0-2019-03-27
       m4_pushdef([TOX_MIN_VERSION], [3.8.0])
       AC_CACHE_CHECK([for tox >= ]TOX_MIN_VERSION, [ac_cv_path_TOX], [
         AC_PATH_PROGS_FEATURE_CHECK([TOX], [tox], [
            tox_version=$($ac_path_TOX --version 2> /dev/null | tail -1)
            AS_IF([test -n "$tox_version"], [
                AX_COMPARE_VERSION([$tox_version], [ge], TOX_MIN_VERSION, [
                    ac_cv_path_TOX="$ac_path_TOX"
                ])
            ])
         ])
       ])
       AS_IF([test -z "$ac_cv_path_TOX"],
             [sage_spkg_install_tox=yes])
])
