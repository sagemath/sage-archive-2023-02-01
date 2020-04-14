SAGE_SPKG_CONFIGURE([suitesparse], [
  SAGE_SPKG_DEPCHECK([openblas], [
     AC_CHECK_HEADER([suitesparse/SuiteSparse_config.h], [
        AC_SEARCH_LIBS([cholmod_speye], [cholmod],[],
         [sage_spkg_install_suitesparse=yes])
	], [
        sage_spkg_install_suitesparse=yes
     ])
  ])
], [], [], [
    AS_IF([test x$sage_spkg_install_suitesparse = xyes], [
       AC_SUBST(SAGE_SUITESPARSE_PREFIX, ['$SAGE_LOCAL'])
    ], [
       AC_SUBST(SAGE_SUITESPARSE_PREFIX, [''])
    ])
])
