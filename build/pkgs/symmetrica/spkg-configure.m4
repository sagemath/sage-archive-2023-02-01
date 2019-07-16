SAGE_SPKG_CONFIGURE([symmetrica], [
    AC_LANG_PUSH(C)
    AC_CHECK_HEADER([symmetrica/def.h], [
dnl check for one of its many functions 
     AC_SEARCH_LIBS([zykelind_tetraeder_edges_extended], [symmetrica], [break],
        [sage_spkg_install_symmetrica=yes])
    ], [sage_spkg_install_symmetrica=yes])
    AC_LANG_POP(C)
])
