AC_ARG_ENABLE([editable],
              [AS_HELP_STRING([--enable-build-as-root],
                              [use an editable install of the Sage library])],
              [AC_SUBST([SAGE_EDITABLE], [yes])])
