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
        os.remove(abs_sage_path(env, "c_lib/libcsage.so"))
    except:
        pass

def buildclib(env, gccc):

    SAGE_LOCAL = env.options['SAGE_LOCAL']
    efw = extfilewalker()
    efw.addcallback('.c',lambda x: True)
    efw.addcallback('.cc',lambda x: True)
    efw.addcallback('.cpp',lambda x: True)
    c_list = efw.walk('devel/sage/c_lib/src')
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
    if env.options["UNAME"]=="Darwin":
        outfile = abs_sage_path(env, "c_lib/libcsage.dylib")
    else:
	outfile = abs_sage_path(env, "c_lib/libcsage.so")
    gcceso = GCC_extension_shared_object(gccc,env, c_ext_dict.values(), fdir, outfile = outfile)
    gcceso.libraries = ['ntl', 'gmp', 'pari']
    for x in c_ext_dict.values():
        gcceso.generate_action(env).dep_register(x.generate_action(env))
    gcceso.generate_action(env).enable()
