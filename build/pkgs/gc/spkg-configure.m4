SAGE_SPKG_CONFIGURE([gc], [
    SAGE_SPKG_DEPCHECK([libatomic_ops], [
        dnl  checking with pkg-config
        PKG_CHECK_MODULES([GC], [bdw-gc-threaded >= 7.6.4], [], [
          PKG_CHECK_MODULES([GC], [bdw-gc >= 7.6.4], [], [
            sage_spkg_install_gc=yes])])
    ])
])
