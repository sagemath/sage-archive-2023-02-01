SAGE_SPKG_CONFIGURE([fflas_ffpack], [
  dnl A check for a system package is not implemented yet.
  sage_spkg_install_fflas_ffpack=yes

  dnl https://github.com/linbox-team/fflas-ffpack/blob/master/macros/instr_set.m4
  dnl discovers these flags from the processor but fails to check whether
  dnl compiler (and assembler) actually support these instruction sets.

  AX_CHECK_COMPILE_FLAG([-mavx512f -mavx512vl -mavx512dq], [], [
    AS_VAR_APPEND([SAGE_CONFIGURE_FFLAS_FFPACK], [" --disable-avx512f --disable-avx512vl --disable-avx512dq"])
  ])
  m4_foreach([ISFLAG], [fma, fma4], [
    AX_CHECK_COMPILE_FLAG([-m]ISFLAG, [], [AS_VAR_APPEND]([SAGE_CONFIGURE_FFLAS_FFPACK], [" --disable-]ISFLAG[ "]))
  ])
  AC_SUBST([SAGE_CONFIGURE_FFLAS_FFPACK])
])
