SAGE_SPKG_CONFIGURE([openblas], [
  PKG_CHECK_MODULES([OPENBLAS], [openblas >= 0.2.20], [
    PKG_CHECK_VAR([OPENBLASPCDIR], [openblas], [pcfiledir], [
       AC_CONFIG_LINKS([
         $SAGE_LOCAL/lib/pkgconfig/blas.pc:$OPENBLASPCDIR/openblas.pc
         $SAGE_LOCAL/lib/pkgconfig/cblas.pc:$OPENBLASPCDIR/openblas.pc
         $SAGE_LOCAL/lib/pkgconfig/lapack.pc:$OPENBLASPCDIR/openblas.pc])
       ], [
       AC_MSG_WARN([Unable to locate the directory of openblas.pc. This should not happen!])
       sage_spkg_install_openblas=yes
       ])
    ], [sage_spkg_install_openblas=yes])
  ], [
  AS_IF([test "x$with_blas" = xopenblas], [
     sage_require_openblas=yes
     sage_require_atlas=no])
  ], [
  AC_MSG_CHECKING([BLAS library])
  AC_ARG_WITH([blas],
  [AS_HELP_STRING([--with-blas=openblas],
    [use OpenBLAS as BLAS library (default)])]
  [AS_HELP_STRING([--with-blas=atlas],
    [use ATLAS as BLAS library])],,
    [with_blas=openblas]  # default
  )
  AS_CASE(["$with_blas"],
    [openblas], [],
    [atlas],    [sage_spkg_install_openblas=no],
                [AC_MSG_ERROR([allowed values for --with-blas are 'atlas' and 'openblas'])])
  AC_MSG_RESULT([$with_blas])
  AC_SUBST([SAGE_BLAS], [$with_blas])
  ]
)
