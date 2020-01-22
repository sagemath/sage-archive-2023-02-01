AC_DEFUN([SAGE_SPKG_ENABLE], [
  m4_pushdef([SPKG_NAME], [$1])
  m4_pushdef([SPKG_TYPE], [$2])
  AC_ARG_ENABLE(SPKG_NAME,
      AS_HELP_STRING([--enable-]SPKG_NAME={no|if_installed|yes},
                     [enable build and use of the $2 package ]SPKG_NAME[ (default: "if_installed")])
AS_HELP_STRING([], [package information: ./sage -info ]SPKG_NAME)
AS_HELP_STRING([--disable-]SPKG_NAME,
                     [disable build and uninstall if previously installed by Sage in PREFIX; same as --enable-]SPKG_NAME[=no]),
      AS_VAR_SET([sage_enable_]SPKG_NAME, [$enableval])
      AS_VAR_SET([sage_enable_]SPKG_NAME, [if_installed])
  )
  AS_IF([test "$sage_enable_]SPKG_NAME[" = if_installed],
        AS_IF([test -r "$SAGE_SPKG_INST/$SPKG_NAME"],
              AS_VAR_SET([sage_enable_]SPKG_NAME, [yes]),
              AS_VAR_SET([sage_enable_]SPKG_NAME, [no]))
  )
])
