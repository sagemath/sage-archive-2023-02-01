SAGE_SPKG_CONFIGURE([sqlite], [
  m4_define([sqlite3_min_version_major], [3])
  m4_define([sqlite3_min_version_minor], [8])
  m4_define([sqlite3_min_version_micro], [7])
  m4_define([sqlite3_min_version], [sqlite3_min_version_major.sqlite3_min_version_minor.sqlite3_min_version_micro])
  PKG_CHECK_MODULES([SQLITE],
                    [sqlite3 >= sqlite3_min_version],
                    [],
                    [AC_MSG_CHECKING([sqlite3 >= sqlite3_min_version without pkg-config])
                     dnl https://www.sqlite.org/c3ref/libversion.html
                     dnl https://www.sqlite.org/c3ref/c_source_id.html
                     SQLITE_SAVED_LIBS="$LIBS"
                     LIBS="$LIBS -lsqlite3"
                     AC_RUN_IFELSE([
                       AC_LANG_PROGRAM([[
                                         #include <sqlite3.h>
                                         #include <assert.h>
                                         #include <stdlib.h>
                                         #include <string.h>
                                       ]],
                                       [[
                                         assert( strcmp(sqlite3_libversion(),SQLITE_VERSION)==0 );
                                         if (SQLITE_VERSION_NUMBER < ]]sqlite3_min_version_major[[*1000000 + ]]sqlite3_min_version_minor[[*1000 + ]]sqlite3_min_version_micro[[)
                                           exit(1);
                                         else
                                           exit(0);
                                       ]])
                       ],
                       [AC_MSG_RESULT([yes])],
                       [AC_MSG_RESULT([no])
                        LIBS="$SQLITE_SAVED_LIBS"
                        sage_spkg_install_sqlite=yes])
                    ])
])
