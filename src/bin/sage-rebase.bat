@echo off
rem Author:
rem * Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>
rem
rem Rebase all dlls in the SAGE_LOCAL directory (and its subdirectories),
rem but do not touch the ones already stored in the system database,
rem and do not update it.
rem Note that subsequent calls to 'rebaseall' will not update the Sage dlls.
rem
rem Invoke this script from a Windows command prompt,
rem after adjusting SAGE_LOCAL to the Windows location of the Sage directory,
rem and, if Cygwin is installed in a non-standard location,
rem adjusting CYGWIN_ROOT.

set CYGWIN_ROOT=C:\cygwin\
set SAGE_LOCAL=C:\cygwin\usr\local\sage\

rem Make sure bash can be called from MSDOS prompt:
path %CYGWIN_ROOT%\bin;%path%
rem Suppress warning about MSDOS-style path:
set CYGWIN=%CYGWIN% nodosfilewarning
rem Call the bash script to do the real work:
cd %SAGE_LOCAL%
bash bin\sage-rebase.sh
