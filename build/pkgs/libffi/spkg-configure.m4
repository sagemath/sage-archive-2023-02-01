SAGE_SPKG_CONFIGURE([libffi], [
    AC_SEARCH_LIBS([ffi_call], [ffi], [], [sage_spkg_install_libffi=yes])
    AC_CHECK_HEADERS([ffi.h ffi/ffi.h], [break], [sage_spkg_install_libffi=yes])
])

