SAGE_SPKG_CONFIGURE([suitesparse], [
  SAGE_SPKG_DEPCHECK([openblas], [
        AC_SEARCH_LIBS([cholmod_speye], [cholmod], [
         AC_SEARCH_LIBS([umfpack_di_solve], [umfpack], [
          AC_SEARCH_LIBS([SuiteSparse_version], [suitesparseconfig], [
            AC_CHECK_HEADERS([SuiteSparse_config.h amd.h],
                             [suispa_header_found=yes])
            AC_CHECK_HEADERS([suitesparse/SuiteSparse_config.h suitesparse/amd.h],
                             [suispa_header_found=yes])
            AS_IF([test x$suispa_header_found = xyes], [],
                 [sage_spkg_install_suitesparse=yes])
          ], [sage_spkg_install_suitesparse=yes])
         ], [
         sage_spkg_install_suitesparse=yes])
        ], [
        sage_spkg_install_suitesparse=yes])
      ], [
      sage_spkg_install_suitesparse=yes])
], [], [], [
    AS_IF([test x$sage_spkg_install_suitesparse = xyes], [
       AC_SUBST(SAGE_SUITESPARSE_PREFIX, ['$SAGE_LOCAL'])
    ], [
       AC_SUBST(SAGE_SUITESPARSE_PREFIX, [''])
    ])
])
