SAGE_SPKG_CONFIGURE([giac], [
    SAGE_SPKG_DEPCHECK([pari], [
       AC_CHECK_HEADER([giac/giac.h], [
        AC_SEARCH_LIBS([ConvertUTF16toUTF8], [giac], [
        ], [sage_spkg_install_giac=yes])
       ], [sage_spkg_install_giac=yes])
    ])
])
