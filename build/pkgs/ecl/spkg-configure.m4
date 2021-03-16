SAGE_SPKG_CONFIGURE([ecl], [
  SAGE_SPKG_DEPCHECK([gc gmp mpir], [
    AC_PATH_PROG([ECL_CONFIG], [ecl-config])
      AS_IF([test x$ECL_CONFIG = x], [
        AC_MSG_NOTICE([ecl-config not found. Installing ecl])
        sage_spkg_install_ecl=yes])

      # "CPPFLAGS" is not a typo, the --cflags output from
      # ecl-config typically contains -D and -I flags.
      CPPFLAGS="${CPPFLAGS} $($ECL_CONFIG --cflags)"

      AC_LANG_PUSH([C])
      AC_RUN_IFELSE([AC_LANG_PROGRAM([[
        #include <ecl/config.h>
      ]],[[
	if (ECL_VERSION_NUMBER < 210201) { return 1; }
      ]])],
      [], # if it works
      [
        AC_MSG_NOTICE([ecl found but too old])
        sage_spkg_install_ecl=yes
      ])
      AC_LANG_POP([C])
  ])
])
