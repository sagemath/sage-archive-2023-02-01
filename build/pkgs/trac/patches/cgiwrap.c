/*
 * Copyright 2001-2004 Brandon Long
 * All Rights Reserved.
 *
 * ClearSilver Templating System
 *
 * This code is made available under the terms of the ClearSilver License.
 * http://www.clearsilver.net/license.hdf
 *
 */

#include "cs_config.h"

/* #include <features.h> */
#include "features.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util/neo_misc.h"
#include "util/neo_err.h"
#include "cgi/cgiwrap.h"

typedef struct _cgiwrapper
{
  int argc;
  char **argv;
  char **envp;
  int env_count;

  READ_FUNC read_cb;
  WRITEF_FUNC writef_cb;
  WRITE_FUNC write_cb;
  GETENV_FUNC getenv_cb;
  PUTENV_FUNC putenv_cb;
  ITERENV_FUNC iterenv_cb;

  void *data;
  int emu_init;
} CGIWRAPPER;

static CGIWRAPPER GlobalWrapper = {0, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0};

void cgiwrap_init_std (int argc, char **argv, char **envp)
{
  /* Allow setting of these even after cgiwrap_init_emu is called */
  GlobalWrapper.argc = argc;
  GlobalWrapper.argv = argv;
  GlobalWrapper.envp = envp;
  GlobalWrapper.env_count = 0;
  while (envp[GlobalWrapper.env_count] != NULL) GlobalWrapper.env_count++;

  /* so you can compile the same code for embedded without mods.
   * Note that this setting is global for the lifetime of the program, so
   * you can never reset these values after calling cgiwrap_init_emu by
   * calling cgiwrap_init_std, you'll have to call cgiwrap_init_emu with NULL
   * values to reset */
  if (GlobalWrapper.emu_init) return;

  GlobalWrapper.read_cb = NULL;
  GlobalWrapper.writef_cb = NULL;
  GlobalWrapper.write_cb = NULL;
  GlobalWrapper.getenv_cb = NULL;
  GlobalWrapper.putenv_cb = NULL;
  GlobalWrapper.iterenv_cb = NULL;
  GlobalWrapper.data = NULL;
}

void cgiwrap_init_emu (void *data, READ_FUNC read_cb,
    WRITEF_FUNC writef_cb, WRITE_FUNC write_cb, GETENV_FUNC getenv_cb,
    PUTENV_FUNC putenv_cb, ITERENV_FUNC iterenv_cb)
{
  /* leave argc, argv, envp, env_count alone since we either don't use them, or
   * they are used by the default versions if any of these are passed as NULL.
   * Note that means that if you pass NULL for anything here, you'd better
   * have called cgiwrap_init_std first! */
  GlobalWrapper.data = data;
  GlobalWrapper.read_cb = read_cb;
  GlobalWrapper.writef_cb = writef_cb;
  GlobalWrapper.write_cb = write_cb;
  GlobalWrapper.getenv_cb = getenv_cb;
  GlobalWrapper.putenv_cb = putenv_cb;
  GlobalWrapper.iterenv_cb = iterenv_cb;
  GlobalWrapper.emu_init = 1;
}

NEOERR *cgiwrap_getenv (const char *k, char **v)
{
  if (GlobalWrapper.getenv_cb != NULL)
  {
    *v = GlobalWrapper.getenv_cb (GlobalWrapper.data, k);
  }
  else
  {
    char *s = getenv(k);

    if (s != NULL)
    {
      *v = strdup(s);
      if (*v == NULL)
      {
	return nerr_raise (NERR_NOMEM, "Unable to duplicate env var %s=%s",
	    k, s);
      }
    }
    else
    {
      *v = NULL;
    }
  }
  return STATUS_OK;
}

NEOERR *cgiwrap_putenv (const char *k, const char *v)
{
  if (GlobalWrapper.putenv_cb != NULL)
  {
    if (GlobalWrapper.putenv_cb(GlobalWrapper.data, k, v))
      return nerr_raise(NERR_NOMEM, "putenv_cb says nomem when %s=%s", k, v);
  }
  else
  {
    char *buf;
    int l;
    l = strlen(k) + strlen(v) + 2;
    buf = (char *) malloc(sizeof(char) * l);
    if (buf == NULL)
      return nerr_raise(NERR_NOMEM, "Unable to allocate memory for putenv %s=%s", k, v);
    snprintf (buf, l, "%s=%s", k, v);
    if (putenv (buf))
      return nerr_raise(NERR_NOMEM, "putenv says nomem when %s", buf);
  }
  return STATUS_OK;
}

NEOERR *cgiwrap_iterenv (int num, char **k, char **v)
{
  *k = NULL;
  *v = NULL;
  if (GlobalWrapper.iterenv_cb != NULL)
  {
    int r;

    r = GlobalWrapper.iterenv_cb(GlobalWrapper.data, num, k, v);
    if (r)
      return nerr_raise(NERR_SYSTEM, "iterenv_cb returned %d", r);
  }
  else if (GlobalWrapper.envp != NULL && num < GlobalWrapper.env_count)
  {
    char *c, *s = GlobalWrapper.envp[num];

    c = strchr (s, '=');
    if (c == NULL) return STATUS_OK;
    *c = '\0';
    *k = strdup(s);
    *c = '=';
    if (*k == NULL)
      return nerr_raise(NERR_NOMEM, "iterenv says nomem for %s", s);
    *v = strdup(c+1);
    if (*v == NULL)
    {
      free(*k);
      *k = NULL;
      return nerr_raise(NERR_NOMEM, "iterenv says nomem for %s", s);
    }
  }
  return STATUS_OK;
}

NEOERR *cgiwrap_writef (const char *fmt, ...)
{
  va_list ap;

  va_start (ap, fmt);
  cgiwrap_writevf (fmt, ap);
  va_end (ap);
  return STATUS_OK;
}

NEOERR *cgiwrap_writevf (const char *fmt, va_list ap)
{
  int r;

  if (GlobalWrapper.writef_cb != NULL)
  {
    r = GlobalWrapper.writef_cb (GlobalWrapper.data, fmt, ap);
    if (r)
      return nerr_raise_errno (NERR_IO, "writef_cb returned %d", r);
  }
  else
  {
    vprintf (fmt, ap);
    /* vfprintf(stderr, fmt, ap); */
  }
  return STATUS_OK;
}

NEOERR *cgiwrap_write (const char *buf, int buf_len)
{
  int r;

  if (GlobalWrapper.write_cb != NULL)
  {
    r = GlobalWrapper.write_cb (GlobalWrapper.data, buf, buf_len);
    if (r != buf_len)
      return nerr_raise_errno (NERR_IO, "write_cb returned %d<%d", r, buf_len);
  }
  else
  {
    /* r = fwrite(buf, sizeof(char), buf_len, stderr);  */
    r = fwrite(buf, sizeof(char), buf_len, stdout);
    if (r != buf_len)
      return nerr_raise_errno (NERR_IO, "fwrite returned %d<%d", r, buf_len);
  }
  return STATUS_OK;
}

void cgiwrap_read (char *buf, int buf_len, int *read_len)
{
  if (GlobalWrapper.read_cb != NULL)
  {
    *read_len = GlobalWrapper.read_cb (GlobalWrapper.data, buf, buf_len);
  }
  else
  {
#ifdef __UCLIBC__
    /* According to
     * http://cvs.uclinux.org/cgi-bin/cvsweb.cgi/uClibc/libc/stdio/stdio.c#rev1.28
     * Note: there is a difference in behavior between glibc and uClibc here
     * regarding fread() on a tty stream.  glibc's fread() seems to return
     * after reading all _available_ data even if not at end-of-file, while
     * uClibc's fread() continues reading until all requested or eof or error.
     * The latter behavior seems correct w.r.t. the standards.
     *
     * So, we use read on uClibc.  This may be required on other platforms as
     * well.  Using raw and buffered i/o interchangeably can be problematic,
     * but everyone should be going through the cgiwrap interfaces which only
     * provide this one read function.
     */
     *read_len = read (fileno(stdin), buf, buf_len);
#else
     *read_len = fread (buf, sizeof(char), buf_len, stdin);
#endif
  }
}
