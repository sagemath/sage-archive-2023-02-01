/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: feResource.cc,v 1.7 2005/07/27 09:46:19 Singular Exp $ */
/*
* ABSTRACT: management of resources
*/

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "mod2.h"
#ifdef AIX_4
#define HAVE_PUTENV 1
#endif

#if defined(HPUX_10) || defined(HPUX_9)
extern "C" int setenv(const char *name, const char *value, int overwrite);
#endif


#include "distrib.h"
#include "dError.h"
#if !defined(ESINGULAR) && !defined(TSINGULAR)
#include "febase.h"
#include "omalloc.h"
#else
char* feResource(const char id, int warn = -1);
char* feResource(const char* key, int warn = -1);
#include "dError.c"
#endif

// define RESOURCE_DEBUG for chattering about resource management
// #define RESOURCE_DEBUG

#if defined(MAKE_DISTRIBUTION)
#if defined(ix86_Win) && ! defined(__CYGWIN__)
#define SINGULAR_DEFAULT_DIR "/Singular/"S_VERSION1
#else // unix
#define SINGULAR_DEFAULT_DIR "/usr/local/Singular/"S_VERSION1
#endif
#else // ! defined(MAKE_DISTRIBUTION)
#define SINGULAR_DEFAULT_DIR S_ROOT_DIR
#endif // defined(MAKE_DISTRIBUTION)

/*****************************************************************
 *
 * Declarations: Data  structures
 *
 *****************************************************************/
typedef enum {feResUndef = 0, feResBinary, feResDir, feResFile, feResUrl, feResPath} feResourceType;

typedef struct feResourceConfig_s
{
  char*           key;   // key to identify resource
  char            id;    // char id to identify resource
  feResourceType  type;  // type of Resource
  char*           env;   // env variable to look for
  char*           fmt;   // format string -- see below for epxlaination
  char*                        value; // what it was set to
} feResourceConfig_s;
typedef feResourceConfig_s * feResourceConfig;

// feSprintf transforms format strings as follows:
// 1.) substrings of the form %c (c being a letter) are replaced by respective resource value
// 2.) substrings of the form $string are replaced by value of resp. env variable

// feCleanResource makes furthermore  the following transformations (except for URL resources)
// 1.) '/' characters are replaced by respective directory - separators
// 2.) ';' characters are replaced by respective path separators
static feResourceConfig_s feResourceConfigs[] =
{
  {"SearchPath",    's', feResPath,  NULL,
   "$SINGULARPATH;"
   "%b/LIB;"
   "%b/MOD;"
   "%r/LIB;"
   "%r/../LIB;"
   "%d/LIB;"
   "%d/../LIB"
   ""},
  {"Singular",  'S',    feResBinary,"SINGULAR_EXECUTABLE",  "%d/"S_UNAME"/Singular",""},
  {"BinDir",    'b',    feResDir,   "SINGULAR_BIN_DIR",     "%d/"S_UNAME            ""},
  {"RootDir",   'r',    feResDir,   "SINGULAR_ROOT_DIR",    "%b/..",                ""},
  {"DefaultDir",'d',    feResDir,   "SINGULAR_DEFAULT_DIR",  SINGULAR_DEFAULT_DIR,  ""},
  {"InfoFile",  'i',    feResFile,  "SINGULAR_INFO_FILE",   "%r/info/singular.hlp", ""},
  {"IdxFile",   'x',    feResFile,  "SINGULAR_IDX_FILE",    "%r/doc/singular.idx",  ""},
  {"HtmlDir",   'h',    feResDir,   "SINGULAR_HTML_DIR",    "%r/html",              ""},
#ifdef ix86_Win
  {"HtmlHelpFile",   'C',    feResFile,   "SINGULAR_CHM_FILE",    "%r/doc/Manual.chm",              ""},
#endif
  {"ManualUrl", 'u',    feResUrl,   "SINGULAR_URL",         "http://www.singular.uni-kl.de/Manual/"S_VERSION1,    ""},
  {"ExDir",      'm',   feResDir,   "SINGULAR_EXAMPLES_DIR","%r/examples",              ""},
  {"Path",      'p',    feResPath,  NULL,                   "%b;$PATH",         ""},

#ifdef ESINGULAR
  {"emacs",    'E',    feResBinary, "ESINGULAR_EMACS",      "%b/emacs",              ""},
  {"xemacs",    'A',    feResBinary, "ESINGULAR_EMACS",      "%b/xemacs",              ""},
  {"SingularEmacs",'M',feResBinary, "ESINGULAR_SINGULAR",    "%b/Singular",           ""},
  {"EmacsLoad",'l',    feResFile,   "ESINGULAR_EMACS_LOAD",  "%e/.emacs-singular",             ""},
  {"EmacsDir",  'e',    feResDir,   "ESINGULAR_EMACS_DIR",   "%r/emacs",             ""},
#elif defined(TSINGULAR)
  {"SingularXterm",'M',feResBinary, "TSINGULAR_SINGULAR",    "%b/Singular",           ""},
#ifdef ix86_Win
  {"rxvt",     'X',    feResBinary,"RXVT",                "%b/rxvt",             ""},
#else
  {"xterm",     'X',    feResBinary,"XTERM",                "%b/xterm",             ""},
#endif
#else
  {"EmacsDir",  'e',    feResDir,   "SINGULAR_EMACS_DIR",   "%r/emacs",             ""},
#endif
  {NULL, 0, feResUndef, NULL, NULL, NULL}, // must be the last record
};


/*****************************************************************
 *
 * Declarations: Local variables / functions
 *
 *****************************************************************/
char* feArgv0=NULL;
#define MAXRESOURCELEN 5*MAXPATHLEN

char fePathSep =
#if defined(ix86_Win)
';'
#else
':'
#endif
;

static feResourceConfig feGetResourceConfig(const char id);
static feResourceConfig feGetResourceConfig(const char* key);
static char* feResource(feResourceConfig config, int warn);
static char* feResourceDefault(feResourceConfig config);
static char* feInitResource(feResourceConfig config, int warn);
static char* feGetExpandedExecutable();
static BOOLEAN feVerifyResourceValue(feResourceType type, char* value);
static char* feCleanResourceValue(feResourceType type, char* value);
static char* feCleanUpFile(char* fname);
static char* feCleanUpPath(char* path);
static void mystrcpy(char* d, char* s);
static char* feSprintf(char* s, const char* fmt, int warn = -1);
#if defined(ix86_Win) && defined(__GNUC__)
// utility function of Cygwin32:
extern "C" int cygwin32_posix_path_list_p (const char *path);
#endif

/*****************************************************************
 *
 * Public functions
 *
 *****************************************************************/
char* feResource(const char* key, int warn)
{
  return feResource(feGetResourceConfig(key), warn);
}

char* feGetResource(const char id)
{
  return feResource(feGetResourceConfig(id), -1);
}

char* feResource(const char id, int warn)
{
  return feResource(feGetResourceConfig(id), warn);
}

char* feResourceDefault(const char id)
{
  return feResourceDefault(feGetResourceConfig(id));
}

char* feResourceDefault(const char* key)
{
  return feResourceDefault(feGetResourceConfig(key));
}

void feInitResources(char* argv0)
{
#if defined(ix86_Win) && defined(__GNUC__)
  if (cygwin32_posix_path_list_p (getenv("PATH")))
    fePathSep = ':';
#endif
  feArgv0 = omStrDup(argv0);
#ifdef RESOURCE_DEBUG
  printf("feInitResources: entering with argv0=%s=\n", argv0);
#endif
  // init some Resources
  feResource('b');
  feResource('r');
  // don't complain about stuff when initializing SingularPath
  feResource('s',0);

#if defined(HAVE_SETENV) || defined(HAVE_PUTENV)
  char* path = feResource('p');
#ifdef RESOURCE_DEBUG
  printf("feInitResources: setting path with argv0=%s=\n", path);
#endif
#ifdef HAVE_PUTENV
  if (path != NULL) { char *s=(char *)omAlloc0(strlen(path)+6);
                      sprintf(s,"PATH=%s",path);
                      putenv(s);
                    }
#else
  if (path != NULL) setenv("PATH", path, 1);
#endif
#endif
}

void feReInitResources()
{
  int i = 0;
  while (feResourceConfigs[i].key != NULL)
  {
    if (feResourceConfigs[i].value != "")
    {
      if (feResourceConfigs[i].value != NULL)
        omFree(feResourceConfigs[i].value);
      feResourceConfigs[i].value = "";
    }
    i++;
  }
#ifdef RESOURCE_DEBUG
  printf("feInitResources: entering with argv0=%s=\n", feArgv0);
#endif
  // init some Resources
  feResource('b');
  feResource('r');
  // don't complain about stuff when initializing SingularPath
  feResource('s',0);
}

/*****************************************************************
 *
 * Local functions
 *
 *****************************************************************/
static feResourceConfig feGetResourceConfig(const char id)
{
  int i = 0;
  while (feResourceConfigs[i].key != NULL)
  {
    if (feResourceConfigs[i].id == id) return &(feResourceConfigs[i]);
    i++;
  }
  return NULL;
}

static feResourceConfig feGetResourceConfig(const char* key)
{
  int i = 0;
  while (feResourceConfigs[i].key != NULL)
  {
    if (strcmp(feResourceConfigs[i].key, key) == 0)
      return &(feResourceConfigs[i]);
    i++;
  }
  return NULL;
}

static char* feResource(feResourceConfig config, int warn)
{
  if (config == NULL) return NULL;
  if (config->value != NULL && *(config->value) != '\0') return config->value;
  return feInitResource(config, warn);
}

static char* feResourceDefault(feResourceConfig config)
{
  if (config == NULL) return NULL;
  char* value = (char*) omAlloc(MAXRESOURCELEN);
  feSprintf(value, config->fmt, -1);
  return value;
}

static char* feInitResource(feResourceConfig config, int warn)
{
  assume(config != NULL);
#ifdef RESOURCE_DEBUG
  printf("feInitResource: entering for %s\n", config->key);
#endif

  char value[MAXRESOURCELEN];
  // now we have to work
  // First, check Environment variable
  if (config->env != NULL)
  {
    char* evalue = getenv(config->env);
    if (evalue != NULL)
    {
#ifdef RESOURCE_DEBUG
      printf("feInitResource: Found value from env:%s\n", evalue);
#endif
      strcpy(value, evalue);
      if (config->type == feResBinary  // do not verify binaries
          ||
          feVerifyResourceValue(config->type,
                                feCleanResourceValue(config->type, value)))
      {
#ifdef RESOURCE_DEBUG
        printf("feInitResource: Set value of %s to =%s=\n", config->key, value);
#endif
        config->value = omStrDup(value);
        return config->value;
      }
    }
  }

  *value = '\0';
  // Special treatment of executable
  if (config->id == 'S')
  {
    char* executable = feGetExpandedExecutable();
    if (executable != NULL)
    {
#ifdef RESOURCE_DEBUG
      printf("exec:%s\n", executable);
#endif
      strcpy(value, executable);
#ifdef RESOURCE_DEBUG
      printf("value:%s\n", value);
#endif
      omFree(executable);
    }
  }
  // and bindir
  else if (config->id == 'b')
  {
    char* executable = feResource('S');
#ifdef RESOURCE_DEBUG
      printf("feInitResource: Get %s from %s\n", config->key, executable);
#endif
    if (executable != NULL)
    {
      strcpy(value, executable);
      executable = strrchr(value, DIR_SEP);
      if (executable != NULL) *executable = '\0';
    }
  }

#ifdef RESOURCE_DEBUG
  printf("value:%s\n", value);
#endif

  if (*value == '\0' && config->fmt != NULL )
  {
    feSprintf(value, config->fmt, warn);
  }
  else if (config->fmt == NULL)
  {
    sprintf(value, "Wrong Resource Specification of %s", config->key);
    dReportBug(value);
    return NULL;
  }

  // Clean and verify
  if (feVerifyResourceValue(config->type,
                            feCleanResourceValue(config->type, value)))
  {
#ifdef RESOURCE_DEBUG
    printf("feInitResource: Set value of %s to =%s=\n", config->key, value);
#endif
    config->value = omStrDup(value);
    return config->value;
  }
  else if (config->type == feResBinary)
  {
    // for binaries, search through PATH once more
    char* executable = omFindExec(config->key, value);
    if (executable != NULL)
    {
      if (feVerifyResourceValue(config->type,
                                feCleanResourceValue(config->type, value)))
      {
        config->value = omStrDup(value);
#ifdef RESOURCE_DEBUG
        printf("feInitResource: Set value of %s to =%s=\n", config->key, config->value);
#endif
        return config->value;
      }
    }
  }

  // issue warning if explicitely requested, or if
  // this value is gotten for the first time
  if (warn > 0 || (warn < 0 && config->value != NULL))
  {
    Warn("Could not get %s. ", config->key);
    Warn("Either set environment variable %s to %s,",
         config->env, config->key);
    feSprintf(value, config->fmt, warn);
    Warn("or make sure that %s is at %s", config->key, value);
  }
#ifdef RESOURCE_DEBUG
      printf("feInitResource: Set value of %s to NULL", config->key);
#endif
  config->value = NULL;
  return NULL;
}

static char* feGetExpandedExecutable()
{
  if (feArgv0 == NULL || *feArgv0 == '\0')
  {
    if (feArgv0 == NULL) dReportBug("feArgv0 == NULL");
    else dReportBug("feArgv0 == ''");
    return NULL;
  }
#ifdef ix86_Win // stupid WINNT sometimes gives you argv[0] within ""
  if (*feArgv0 == '"')
  {
    int l = strlen(feArgv0);
    if (feArgv0[l-1] == '"')
    {
      feArgv0[l-1] = '\0';
      feArgv0++;
    }
  }
#endif
#ifdef RESOURCE_DEBUG
  printf("feGetExpandedExecutable: calling find_exec with =%s=\n", feArgv0);
#endif
  char executable[MAXRESOURCELEN];
  char* value = omFindExec(feArgv0, executable);
#ifdef RESOURCE_DEBUG
  printf("feGetExpandedExecutable: find_exec exited with =%s=%d\n", executable, access(executable, X_OK));
#endif
  if (value == NULL)
  {
    char message[MAXRESOURCELEN];
    sprintf(message, "Could not get expanded executable from %s", feArgv0);
    dReportBug(message);
    return NULL;
  }
  return omStrDup(value);
}


static BOOLEAN feVerifyResourceValue(feResourceType type, char* value)
{
#ifdef RESOURCE_DEBUG
  printf("feVerifyResourceValue: entering with =%s=\n", value);
  printf("%d:%d\n", access(value, R_OK), access(value, X_OK));
#endif
  switch(type)
  {
      case feResUrl:
      case feResPath:
        return TRUE;

      case feResFile:
        return ! access(value, R_OK);

      case feResBinary:
      case feResDir:
        return ! access(value, X_OK);

      default:
        return FALSE;
  }
}

/*****************************************************************
 *
 * Cleaning/Transformations of resource values
 *
 *****************************************************************/

static char* feCleanResourceValue(feResourceType type, char* value)
{
  if (value == NULL || *value == '\0') return value;
#ifdef RESOURCE_DEBUG
      printf("Clean value:%s\n", value);
#endif
#ifdef ix86_Win
#ifdef RESOURCE_DEBUG
      printf("Clean WINNT value:%s\n", value);
#endif
  if (type == feResBinary)
  {
    int l = strlen(value);
    if (l < 4 || (strcmp(&value[l-4], ".exe") != 0 &&
                  strcmp(&value[l-4], ".EXE") != 0))
      strcat(value, ".exe");
  }
#endif
  if (type == feResFile || type == feResBinary || type == feResDir)
    return feCleanUpFile(value);
  if (type == feResPath)
    return feCleanUpPath(value);
  return value;
}

static char* feCleanUpFile(char* fname)
{
  char* fn, *s;

#ifdef RESOURCE_DEBUG
  printf("feCleanUpFile: entering with =%s=\n", fname);
#endif
  // Remove unnecessary .. and //
  for (fn = fname; *fn != '\0'; fn++)
  {
    if (*fn == '/')
    {
      if (*(fn+1) == '\0')
      {
        if (fname != fn) *fn = '\0';
        break;
      }
      if (*(fn + 1) == '/' && (fname != fn))
      {
        mystrcpy(fn, fn+1);
        fn--;
      }
      else if (*(fn+1) == '.')
      {
        if (*(fn+2) == '.' && (*(fn + 3) == '/' || *(fn + 3) == '\0'))
        {
          *fn = '\0';
          s = strrchr(fname, '/');
          if (s != NULL)
          {
            mystrcpy(s+1, fn + (*(fn + 3) != '\0' ? 4 : 3));
            fn = s-1;
          }
          else
          {
            *fn = '/';
          }
        }
        else if (*(fn+2) == '/' || *(fn+2) == '\0')
        {
          mystrcpy(fn+1, fn+3);
          fn--;
        }
      }
    }
  }

#ifdef RESOURCE_DEBUG
  printf("feCleanUpFile: leaving with =%s=\n", fname);
#endif
  return fname;
}

// remove duplicates dir resp. those which do not exist
static char* feCleanUpPath(char* path)
{
#ifdef RESOURCE_DEBUG
  printf("feCleanUpPath: entering with: =%s=\n", path);
#endif
  if (path == NULL) return path;

  int n_comps = 1, i, j;
  char* opath = path;
  char** path_comps;

  for (; *path != '\0'; path++)
  {
    if (*path == fePathSep) n_comps++;
    if (*path == ';')
    {
      *path = fePathSep;
      n_comps++;
    }
  }

  path_comps = (char**) omAlloc(n_comps*sizeof(char*));
  path_comps[0]=opath;
  path=opath;
  i = 1;

  if (i < n_comps)
  {
    while (1)
    {
      if (*path == fePathSep)
      {
        *path = '\0';
        path_comps[i] = path+1;
        i++;
        if (i == n_comps) break;
      }
      path++;
    }
  }

  for (i=0; i<n_comps; i++)
    path_comps[i] = feCleanUpFile(path_comps[i]);
#ifdef RESOURCE_DEBUG
  PrintS("feCleanUpPath: after CleanUpName: ");
  for (i=0; i<n_comps; i++)
    Print("%s:", path_comps[i]);
  Print("\n");
#endif

  for (i=0; i<n_comps;)
  {
#ifdef RESOURCE_DEBUG
    if (access(path_comps[i], X_OK | R_OK))
      Print("feCleanUpPath: remove %d:%s -- can not access\n", i, path_comps[i]);
#endif
    if ( ! access(path_comps[i], X_OK | R_OK))
    {
      // x- permission is granted -- we assume that it is a dir
      for (j=0; j<i; j++)
      {
        if (strcmp(path_comps[j], path_comps[i]) == 0)
        {
          // found a duplicate
#ifdef RESOURCE_DEBUG
          Print("feCleanUpPath: remove %d:%s -- equal to %d:%s\n", j, path_comps[j], i, path_comps[i]);
#endif
          j = i+1;
          break;
        }
      }
      if (j == i)
      {
        i++;
        continue;
      }
    }
    // now we can either not access or found a duplicate
    path_comps[i] = NULL;
    for (j=i+1; j<n_comps; j++)
        path_comps[j-1] = path_comps[j];
    n_comps--;
  }


  // assemble everything again
  for (path=opath, i=0;i<n_comps-1;i++)
  {
    mystrcpy(path, path_comps[i]);
    path += strlen(path);
    *path = fePathSep;
    path++;
  }
  if (n_comps)
  {
    mystrcpy(path, path_comps[i]);
  }
  else
  {
    *opath = '\0';
  }
  omFree(path_comps);
#ifdef RESOURCE_DEBUG
  Print("feCleanUpPath: leaving with path=%s=\n", opath);
#endif
  return opath;
}

// strcpy where source and destination may overlap
static void mystrcpy(char* d, char* s)
{
  assume(d != NULL && s != NULL);
  while (*s != '\0')
  {
    *d = *s;
    d++;
    s++;
  }
  *d = '\0';
}

/*****************************************************************
 *
 * feSprintf
 *
 *****************************************************************/
static char* feSprintf(char* s, const char* fmt, int warn)
{
  char* s_in = s;
  if (fmt == NULL) return NULL;

  while (*fmt != '\0')
  {
    *s = *fmt;

    if (*fmt == '%' && *(fmt + 1) != '\0')
    {
      fmt++;
      char* r = feResource(*fmt, warn);
      if (r != NULL)
      {
        strcpy(s, r);
        s += strlen(r) - 1;
      }
      else
      {
        s++;
        *s = *fmt;
      }
    }
    else if (*fmt == '$' && *(fmt + 1) != '\0')
    {
      fmt++;
      char* v = s + 1;
      while (*fmt == '_' ||
             (*fmt >= 'A' && *fmt <= 'Z') ||
             (*fmt >= 'a' && *fmt <= 'z'))
      {
        *v = *fmt;
        v++;
        fmt++;
      }
      fmt--;
      *v = '\0';
      v = getenv(s + 1);
      if (v != NULL) strcpy(s, v);
      s += strlen(s) - 1;
    }
    s++;
    fmt++;
  }
  *s = '\0';
  return s_in;
}

void feStringAppendResources(int warn)
{
  int i = 0;
  char* r;
  StringAppend("%-10s:\t%s\n", "argv[0]", feArgv0);
  while (feResourceConfigs[i].key != NULL)
  {
    r = feResource(feResourceConfigs[i].key, warn);
    StringAppend("%-10s:\t%s\n", feResourceConfigs[i].key,
                 (r != NULL ? r : ""));
    i++;
  }
}
