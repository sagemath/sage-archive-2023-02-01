import os, signal, sys, time, thread, threading, tempfile

from build.all import *
verbose=999
TM = taskmanager.TM

def abs_sage_path(env, path):
    return os.path.realpath(env.options['SAGE_ROOT'] + '/devel/sage/' + path)
def abs_sage_local(env, path):
    return os.path.realpath(env.options['SAGE_LOCAL'] + path)

def build_clib_clean(env):
    try:
        os.remove(env.options['SAGE_ROOT']+'/config.h')
    except:
        pass
    try:
        os.remove(abs_sage_path(env, "c_lib/libcsage.so"))
    except:
        pass
    try:
        os.remove(env.options['SAGE_ROOT'] + "config.chc")
    except:
        pass

def buildclib(env, gccc):
    _build_config(env, gccc)

    SAGE_LOCAL = env.options['SAGE_LOCAL']
    efw = extfilewalker()
    efw.addcallback('.c',lambda x: True)
    efw.addcallback('.cc',lambda x: True)
    efw.addcallback('.cpp',lambda x: True)
    c_list = efw.walk('devel/sage/c_lib/src')
    c_list.remove(abs_sage_path(env,'c_lib/src/gthread-win32.c'))
    c_list.remove(abs_sage_path(env,'c_lib/src/gthread-posix.c'))
    c_list.remove(abs_sage_path(env,'c_lib/src/gthread-none.c'))
    c_ext_dict = { }
    for x in c_list:
        ext = os.path.splitext(x)[1]
        fdir = os.path.split(x)[0]
        gcceo = GCC_extension_object(gccc, env, [x], fdir, options = { '-fPIC':None } )
        gcceo.generate_action(env).enable()
        c_ext_dict[x] = gcceo
        if ext == ".cpp":
            language = "C++"
        else:
            language = "C"
        config_gcc_file(env, c_ext_dict, x, language = language , include_dirs = [abs_sage_local(env, 'include/NTL/'), abs_sage_path(env, 'devel/sage/c_lib/include')], libraries = ['ntl', 'gmp', 'pari'])
    gcceso = GCC_extension_shared_object(gccc,env, c_ext_dict.values(), fdir, outfile = abs_sage_path(env, "c_lib/libcsage.so"))
    for x in c_ext_dict.values():
        gcceso.generate_action(env).dep_register(x.generate_action(env))
    gcceso.generate_action(env).enable()

def _build_config(env, gccc):
    f = file("config.h", mode="w")
    f.write('#include "gmacros.h"\n')
    confdict = read_config_cache("config.chc")
    if confdict == None:
        confdict = { }
        create_config_file(gccc,f, confdict)
        sage_dict_tests(f, confdict, gccc)
        output_config_file(f, confdict)
        sage_write_config(f, confdict, gccc)
        write_config_cache("config.chc", confdict)
    else:
        output_config_file(f, confdict)
        sage_write_config(f, confdict, gccc)
    f.close()
    cmd = 'cp %s %s' % (env.options['SAGE_ROOT']+'/config.h', env.options['SAGE_ROOT']+'/devel/sage/c_lib/include')
    os.system(cmd)



def sage_dict_tests(f, dict, gccc):
    if dict["SIZEOF_VOID_P"] == dict["SIZEOF_INT"]:

        compiledata = """
typedef union _GSystemThread GSystemThread;
union _GSystemThread
{
char   data[4];
double dummy_double;
void  *dummy_pointer;
long   dummy_long;
};
int main() { return sizeof(GSystemThread); }"""
        ret = gccc._try_code(compiledata)
        dict["GLIB_SIZEOF_SYSTEM_THREAD"] = int(ret)



        compiledata = """
typedef struct _GStaticMutex GStaticMutex;
struct _GStaticMutex
{
struct _GMutex *runtime_mutex;
union {
    char   pad[24];
    double dummy_double;
    void  *dummy_pointer;
    long   dummy_long;
} static_mutex;
};
int main() { return sizeof(GStaticMutex); }"""
        ret = gccc._try_code(compiledata)
        dict["GLIB_SIZEOF_GMUTEX"] = int(ret)

    else:
        compiledata = """
typedef union _GSystemThread GSystemThread;
union _GSystemThread
{
char   data[8];
double dummy_double;
void  *dummy_pointer;
long   dummy_long;
};
int main() { return sizeof(GSystemThread); }"""
        ret = gccc._try_code(compiledata)
        dict["GLIB_SIZEOF_SYSTEM_THREAD"] = int(ret)



        compiledata = """
typedef struct _GStaticMutex GStaticMutex;
struct _GStaticMutex
{
struct _GMutex *runtime_mutex;
union {
    char   pad[40];
    double dummy_double;
    void  *dummy_pointer;
    long   dummy_long;
} static_mutex;
};
int main() { return sizeof(GStaticMutex); }"""
        ret = gccc._try_code(compiledata)
        dict["GLIB_SIZEOF_GMUTEX"] = int(ret)

def sage_write_config(f, dict, gccc):
    if dict["SIZEOF_VOID_P"] == dict["SIZEOF_INT"]:
        outtxt = """
#define GPOINTER_TO_INT(p) ((gint) (p))
#define GPOINTER_TO_UINT(p) ((guint) (p))
#define GINT_TO_POINTER(i) ((gpointer) (i))
#define GUINT_TO_POINTER(u) ((gpointer) (u))
typedef signed long long gint64;
typedef unsigned long long guint64;
#define G_GINT64_CONSTANT(val)	(val##LL)
#define G_GUINT64_CONSTANT(val)	(val##ULL)
#define G_GINT64_MODIFIER "ll"
#define G_GINT64_FORMAT "lli"
#define G_GUINT64_FORMAT "lul"
typedef unsigned long gsize;
#define G_GSIZE_MODIFIER "l"
#define G_GSSIZE_FORMAT "li"
#define G_GSIZE_FORMAT "lu"
typedef struct _GStaticMutex GStaticMutex;
struct _GStaticMutex
{
  struct _GMutex *runtime_mutex;
  union {
    char   pad[24];
    double dummy_double;
    void  *dummy_pointer;
    long   dummy_long;
  } static_mutex;
};
#define	G_STATIC_MUTEX_INIT	{ NULL, { { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} } }
#define	g_static_mutex_get_mutex(mutex) \
  (g_thread_use_default_impl ? ((GMutex*) ((mutex)->static_mutex.pad)) : \
   g_static_mutex_get_mutex_impl_shortcut (&((mutex)->runtime_mutex)))
/* This represents a system thread as used by the implementation. An
 * alien implementaion, as loaded by g_thread_init can only count on
 * "sizeof (gpointer)" bytes to store their info. We however need more
 * for some of our native implementations. */
typedef union _GSystemThread GSystemThread;
union _GSystemThread
{
  char   data[4];
  double dummy_double;
  void  *dummy_pointer;
  long   dummy_long;
};
"""
    else:
        outtxt = """
#define GPOINTER_TO_INT(p) ((gint)  (glong) (p))
#define GPOINTER_TO_UINT(p) ((guint) (gulong) (p))
#define GINT_TO_POINTER(i) ((gpointer) (glong) (i))
#define GUINT_TO_POINTER(u) ((gpointer) (gulong) (u))
typedef signed long gint64;
typedef unsigned long guint64;
#define G_GINT64_CONSTANT(val)	(val##L)
#define G_GUINT64_CONSTANT(val)	(val##UL)
#define G_GINT64_MODIFIER "l"
#define G_GINT64_FORMAT "li"
#define G_GUINT64_FORMAT "lu"
typedef signed long gssize;
typedef unsigned long gsize;
#define G_GSIZE_MODIFIER ""
#define G_GSSIZE_FORMAT "i"
#define G_GSIZE_FORMAT "u"
typedef struct _GStaticMutex GStaticMutex;
struct _GStaticMutex
{
struct _GMutex *runtime_mutex;
union {
    char   pad[40];
    double dummy_double;
    void  *dummy_pointer;
    long   dummy_long;
} static_mutex;
};
#define	G_STATIC_MUTEX_INIT	{ NULL, { { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} } }
#define	g_static_mutex_get_mutex(mutex) \
(g_thread_use_default_impl ? ((GMutex*) ((mutex)->static_mutex.pad)) : \
g_static_mutex_get_mutex_impl_shortcut (&((mutex)->runtime_mutex)))
/* This represents a system thread as used by the implementation. An
* alien implementaion, as loaded by g_thread_init can only count on
* "sizeof (gpointer)" bytes to store their info. We however need more
* for some of our native implementations. */
typedef union _GSystemThread GSystemThread;
union _GSystemThread
{
char   data[8];
double dummy_double;
void  *dummy_pointer;
long   dummy_long;
};
"""
    f.write(outtxt)
    outtxt = """#include <limits.h>
#define G_THREADS_ENABLED
#define G_THREADS_IMPL_POSIX
#define GLIB_SYSDEF_POLLIN =1
#define GLIB_SYSDEF_POLLOUT =4
#define GLIB_SYSDEF_POLLPRI =2
#define GLIB_SYSDEF_POLLHUP =16
#define GLIB_SYSDEF_POLLERR =8
#define GLIB_SYSDEF_POLLNVAL =32
typedef int GPid;
#if defined(__SUNPRO_C) && (__SUNPRO_C >= 0x550)
#define G_GNUC_INTERNAL __hidden
#elif defined (__GNUC__) && defined (G_HAVE_GNUC_VISIBILITY)
#define G_GNUC_INTERNAL __attribute__((visibility("hidden")))
#else
#define G_GNUC_INTERNAL
#endif
typedef signed char gint8;
typedef unsigned char guint8;
typedef signed short gint16;
typedef unsigned short guint16;
#define G_GINT16_MODIFIER "h"
#define G_GINT16_FORMAT "hi"
#define G_GUINT16_FORMAT "hu"
typedef signed int gint32;
typedef unsigned int guint32;
#define G_GINT32_MODIFIER ""
#define G_GINT32_FORMAT "i"
#define G_GUINT32_FORMAT "u"
#define G_HAVE_GINT64 1          /* deprecated, always true */

#define G_MINFLOAT	FLT_MIN
#define G_MAXFLOAT	FLT_MAX
#define G_MINDOUBLE	DBL_MIN
#define G_MAXDOUBLE	DBL_MAX
#define G_MINSHORT	SHRT_MIN
#define G_MAXSHORT	SHRT_MAX
#define G_MAXUSHORT	USHRT_MAX
#define G_MININT	INT_MIN
#define G_MAXINT	INT_MAX
#define G_MAXUINT	UINT_MAX
#define G_MINLONG	LONG_MIN
#define G_MAXLONG	LONG_MAX
#define G_MAXULONG	ULONG_MAX

#define g_memmove(dest,src,len) G_STMT_START { memmove ((dest), (src), (len)); } G_STMT_END

/* Define to 1 if you have the `memalign' function. */
#define HAVE_MEMALIGN 1

/* Define to 1 if you have the `memmove' function. */
#define HAVE_MEMMOVE 1

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Have a monotonic clock */
#define HAVE_MONOTONIC_CLOCK 1

/* Define to 1 if you have the `nanosleep' function. */
#define HAVE_NANOSLEEP 1
/* Maximum POSIX RT priority */
#define POSIX_MAX_PRIORITY sched_get_priority_max(SCHED_OTHER)

/* define if posix_memalign() can allocate any size */
#define POSIX_MEMALIGN_WITH_COMPLIANT_ALLOCS 1

/* The POSIX RT yield function */
#define POSIX_YIELD_FUNC sched_yield()

/* Minimum POSIX RT priority */
#define POSIX_MIN_PRIORITY sched_get_priority_min(SCHED_OTHER)

#define G_THREAD_SOURCE "gthread-posix.c"

#define g_return_if_fail(expr)			G_STMT_START{ (void)0; }G_STMT_END
#define g_return_val_if_fail(expr,val)		G_STMT_START{ (void)0; }G_STMT_END
#define g_return_if_reached()			G_STMT_START{ return; }G_STMT_END
#define g_return_val_if_reached(val)		G_STMT_START{ return (val); }G_STMT_END
#define g_assert(expr)		G_STMT_START{ (void)0; }G_STMT_END
#define g_assert_not_reached()	G_STMT_START{ (void)0; }G_STMT_END

#define G_OS_UNIX
"""
    f.write(outtxt)




    outtxt = """
#endif
"""
    f.write(outtxt)