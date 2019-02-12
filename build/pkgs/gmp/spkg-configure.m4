SAGE_SPKG_CONFIGURE([gmp], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_MPIR])
    if test x"$SAGE_MP_LIBRARY" = xgmp; then
        sage_spkg_install_gmp=yes
    fi
])
