SAGE_SPKG_CONFIGURE([zlib], [
    AX_CHECK_ZLIB([], [sage_spkg_install_zlib=yes])
])
