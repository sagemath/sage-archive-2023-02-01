SAGE_SPKG_CONFIGURE([xz], [
        AC_CACHE_CHECK([for xz >= 4.999.0], [ac_cv_path_XZ], [
            AC_PATH_PROGS_FEATURE_CHECK([XZ], [xz], [
              xz_version=`$ac_path_XZ --version 2>&1 | cut -d' ' -f4 | $SED -n 1p`
              AS_IF([test -n "$xz_version"], [
                  AX_COMPARE_VERSION([$xz_version], [ge], [4.999.0], [
                      ac_cv_path_XZ="$ac_path_XZ"
                      ac_path_XZ_found=:
                  ])
              ])
            ])
            AS_IF([test -z "$ac_cv_path_XZ"], [sage_spkg_install_xz=yes])
        ])
])
