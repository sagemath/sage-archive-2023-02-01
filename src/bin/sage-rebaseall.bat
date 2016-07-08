@echo off
rem Authors:
rem * Dmitrii Pasechnik <dimpase@gmail.com>
rem * Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>
rem
rem Rebase all dlls in the SAGE_LOCAL directory (and its subdirectories),
rem as well as the ones already stored in the system database,
rem and update the database.
rem This system-wide database is located in '/etc/rebase.db.i386' and
rem includes the Cygwin system dlls.
rem
rem Invoke this script from a Windows command prompt,
rem after adjusting SAGE_LOCAL to the Windows location of the Sage directory,
rem and, if Cygwin is installed in a non-standard location,
rem adjusting CYGWIN_ROOT.
rem Ensure that no other Cygwin processes are currently running.
rem Note that you need write access to the system-wide rebase database
rem (which usually means admin rights).

set CYGWIN_ROOT=C:\cygwin\
set SAGE_LOCAL=C:\cygwin\usr\local\sage\

rem Make sure dash can be called from MSDOS prompt:
path %CYGWIN_ROOT%\bin;%path%
rem Suppress warning about MSDOS-style path:
set CYGWIN=%CYGWIN% nodosfilewarning
rem Call the dash script to do the real work:
cd %SAGE_LOCAL%
dash bin\sage-rebaseall.sh
