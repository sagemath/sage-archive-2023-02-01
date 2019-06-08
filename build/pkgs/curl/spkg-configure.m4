SAGE_SPKG_CONFIGURE([curl], [
    # Only install curl (and libcurl) if we cannot find version 7.22 or
    # later, since that is what R needs.
    AC_CACHE_CHECK([for curl 7.22], [ac_cv_path_CURL], [
        AC_PATH_PROGS_FEATURE_CHECK([CURL], [curl], [
          ${ac_path_CURL}-config --checkfor 7.22 >/dev/null 2>/dev/null && ac_cv_path_CURL=${ac_path_CURL}
        ])
    ])
    AS_IF([test -z "$ac_cv_path_CURL"], [sage_spkg_install_curl=yes])
    LIBCURL_CHECK_CONFIG(, 7.22, , [sage_spkg_install_curl=yes])
])
