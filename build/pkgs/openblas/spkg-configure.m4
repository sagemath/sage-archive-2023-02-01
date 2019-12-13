SAGE_SPKG_CONFIGURE([openblas], [
  PKG_CHECK_MODULES([OPENBLAS], [openblas >= 0.3.5], [
    PKG_CHECK_VAR([OPENBLASPCDIR], [openblas], [pcfiledir], [
       AC_CONFIG_LINKS([
         $SAGE_LOCAL/lib/pkgconfig/blas.pc:$OPENBLASPCDIR/openblas.pc
         $SAGE_LOCAL/lib/pkgconfig/lapack.pc:$OPENBLASPCDIR/openblas.pc])
       ], [
       AC_MSG_WARN([Unable to locate the directory of openblas.pc. This should not happen!])
       sage_spkg_install_openblas=yes
       ])
    ], [sage_spkg_install_openblas=yes])
  ], [
  AS_IF([test "x$with_blas" = xopenblas], [sage_require_openblas=yes])
  ]
)
