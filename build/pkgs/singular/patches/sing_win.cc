/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: sing_win.cc,v 1.7 2008/06/10 14:39:43 wienand Exp $ */

/*
* ABSTRACT: Windows specific routines
*/

#include "mod2.h"
#ifdef ix86_Win
#include <windows.h>
#include <winuser.h>
#include <htmlhelp.h>
#include <sys/cygwin.h>
#include <stdio.h>
#ifndef MAXPATHLEN
#define MAXPATHLEN 1024
#endif

void heOpenWinHtmlHelp(const char* keyw, char* helppath )
// API Call Sequence for Microsoft HTML Help System
{
  char path[MAXPATHLEN];
#ifdef TEST
  printf("keyw:%s\n", keyw);
#endif
  cygwin_conv_to_full_win32_path(helppath, path);
#ifdef TEST
  printf("path:%s\n", path);
#endif
HH_ALINKA link;
   link.cbStruct =     sizeof(HH_ALINKA) ;
   link.fReserved =    FALSE ;
   link.pszKeywords =  keyw;
   link.pszUrl =       NULL ;
   link.pszMsgText =   NULL ;
   link.pszMsgTitle =  NULL ;
   link.pszWindow =    NULL ;
   link.fIndexOnFail = TRUE ;
   /* Commented out, since this breaks building on *Cygwin* for SAGE, which evidently
      isn't setup to open the standard Windows help system.  Since we use Singular
      as a library anyways, this is reasonable to comment out for now. */
   /* HtmlHelpA(NULL, path, HH_KEYWORD_LOOKUP, (DWORD)&link); */
}

void heOpenWinntHlp(const char* keyw, char* helppath )
{
  char path[MAXPATHLEN];
#ifdef TEST
  printf("keyw:%s\n", keyw);
#endif
  cygwin_conv_to_full_win32_path(helppath, path);
#ifdef TEST
  printf("path:%s\n", path);
#endif
      WinHelp(NULL, path, HELP_PARTIALKEY,(DWORD)keyw);
}

void heOpenWinntUrl(const char* url, int local)
{
#ifdef TEST
  printf("url:%s:local:%d\n", url, local);
#endif
  if (local)
  {
    char path[MAXPATHLEN];
    char *p;
    cygwin_conv_to_full_win32_path(url, path);
    /* seems like I can not open url's wit # at the end */
    if ((p=strchr(path, '#')) != NULL) *p = '\0';
#ifdef TEST
    printf("path:%s:local:%d\n", path, local);
#endif
    ShellExecute(NULL, "open", path, 0, 0, SW_SHOWNORMAL);
  }
  else
  {
    // need to check whether this works
    ShellExecute(NULL, "open", url, 0, 0, SW_SHOWNORMAL);
  }
}

void *dynl_open(char *filename)
{
  char path[MAXPATHLEN];
  cygwin_conv_to_full_win32_path(filename, path);
  HINSTANCE hLibrary = LoadLibrary( TEXT(path));

  return(hLibrary);
}

void *dynl_sym(void *handle, char *symbol)
{
  FARPROC f;
  f = GetProcAddress((HINSTANCE)handle, TEXT (symbol));
  return((void*) f);
}

int dynl_close (void *handle)
{
  FreeLibrary((HINSTANCE)handle);
  return(0);
}

const char *dynl_error()
{
  static char errmsg[] = "support for dynamic loading not implemented";

  return errmsg;
}

#endif /*ix86_Win */
