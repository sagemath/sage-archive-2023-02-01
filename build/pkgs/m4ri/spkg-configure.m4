SAGE_SPKG_CONFIGURE([m4ri], [
    SAGE_SPKG_DEPCHECK([m4ri], [libpng], [
        PKG_CHECK_MODULES([M4RI], [m4ri >= 20140914], [], [
           sage_spkg_install_m4ri=yes])
    ])
])
