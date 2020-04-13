SAGE_SPKG_CONFIGURE([suitesparse], [
  SAGE_SPKG_DEPCHECK([openblas], [
     AC_CHECK_HEADER([suitesparse/SuiteSparse_config.h], [
        AC_SEARCH_LIBS([cholmod_speye], [cholmod],[],
         [sage_spkg_install_suitesparse=yes])
	], [
        sage_spkg_install_suitesparse=yes
     ])
  ])
])
