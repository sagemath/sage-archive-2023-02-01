SAGE_SPKG_CONFIGURE([openblas], [
 dnl CHECK
 SAGE_SPKG_DEPCHECK([gfortran], [
  SAVE_LIBS="$LIBS"
  SAVE_CFLAGS="$CFLAGS"
  PKG_CHECK_MODULES([OPENBLAS], [openblas >= 0.2.20], [
    LIBS="$OPENBLAS_LIBS $LIBS"
    CFLAGS="$OPENBLAS_CFLAGS $CFLAGS"
    PKG_CHECK_VAR([OPENBLASPCDIR], [openblas], [pcfiledir], [
       sage_install_blas_pc=yes
       AC_CHECK_FUNC([cblas_dgemm], [dnl openblas works as cblas
             sage_install_cblas_pc=yes
             ], [
             dnl openblas does not work as cblas; try to use system cblas as is
             PKG_CHECK_MODULES([CBLAS], [cblas], [], [sage_spkg_install_openblas=yes])
          ])
       dnl Check all name manglings that AC_FC_FUNC could check based on the
       dnl characteristics of the Fortran compiler
       m4_foreach([dgeqrf_mangled], [dgeqrf, dgeqrf_, DGEQRF, DGEQRF_], [
          AC_CHECK_FUNC(dgeqrf_mangled, [
             AS_VAR_SET([HAVE_DGEQRF], [yes])
          ])
       ])
       AS_IF([test x$HAVE_DGEQRF = xyes], [dnl openblas works as lapack
             sage_install_lapack_pc=yes
             ], [
             dnl openblas does not work as lapack; try to use system lapack as is
             PKG_CHECK_MODULES([LAPACK], [lapack], [], [sage_spkg_install_openblas=yes])
          ])
       ], [
       AC_MSG_WARN([Unable to locate the directory of openblas.pc. This should not happen!])
       sage_spkg_install_openblas=yes
       ])
     AS_IF([test x$sage_spkg_install_openblas != xyes], [
        AC_SUBST([SAGE_SYSTEM_FACADE_PC_FILES])
        AC_SUBST([SAGE_OPENBLAS_PC_COMMAND], ["\$(LN) -sf \"$OPENBLASPCDIR/openblas.pc\" \"\$(@)\""])
        m4_foreach([blaslibnam], [blas, cblas, lapack], [
         AS_IF([test x$sage_install_]blaslibnam[_pc = xyes], [
            AS_VAR_APPEND([SAGE_SYSTEM_FACADE_PC_FILES], [" \$(SAGE_PKGCONFIG)/]blaslibnam[.pc"])
         ])
        ])
     ])
    ], [
      dnl No openblas.pc
      AC_SEARCH_LIBS([cblas_dgemm], [openblas cblas blas], [
        AS_VAR_SET([HAVE_CBLAS_DGEMM], [yes])
        AS_IF([test x"$ac_cv_search_cblas_dgemm" != x"none required"], [
          AS_VAR_APPEND([OPENBLAS_LIBS], ["$ac_cv_search_cblas_dgemm "])
        ])
      ], [], [-lgfortran])
      m4_foreach([dgeqrf_mangled], [dgeqrf, dgeqrf_, DGEQRF, DGEQRF_], [
         AC_SEARCH_LIBS(dgeqrf_mangled, [openblas lapack], [
           AS_VAR_SET([HAVE_DGEQRF], [yes])
           AS_IF([test x"$ac_cv_search_]dgeqrf_mangled[" != x"none required"], [
             AS_VAR_APPEND([OPENBLAS_LIBS], ["$ac_cv_search_]dgeqrf_mangled[ "])
           ])
         ], [], [-lgfortran])
      ])
      AS_IF([test x"$HAVE_CBLAS_DGEMM" = xyes -a x"$HAVE_DGEQRF" = xyes], [
        AC_SUBST([OPENBLAS_LIBS])
        AC_SUBST([SAGE_SYSTEM_FACADE_PC_FILES])
        AC_SUBST([SAGE_OPENBLAS_PC_COMMAND], ["  (echo \"Name: openblas\"; echo \"Description: OpenBLAS\"; echo \"Version: 0.3\"; echo \"Libs: $OPENBLAS_LIBS\") > \"\$(@)\""])
        m4_foreach([blaslibnam], [openblas, blas, cblas, lapack], [
          AS_VAR_APPEND([SAGE_SYSTEM_FACADE_PC_FILES], [" \$(SAGE_PKGCONFIG)/]blaslibnam[.pc"])
        ])
      ], [
        dnl No system BLAS found
        sage_spkg_install_openblas=yes
      ])
    ])
  LIBS="$SAVE_LIBS"
  CFLAGS="$SAVE_CFLAGS"
 ])
 ], [
  dnl REQUIRED-CHECK
  AS_IF([test "x$with_blas" = xopenblas], [
     sage_require_openblas=yes
     sage_require_atlas=no])
  ], [
  dnl PRE
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
