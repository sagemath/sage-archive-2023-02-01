import sys, pickle
from compilers.gcc import *

def get_type_size(gccc, ctype):
    compiledata = 'int main() { return sizeof(%s); }' % ctype
    ret = gccc._try_code(compiledata)
    if ret == -1:
        return 0
    else:
        return ret

def has_header(gccc, header):
    compiledata = '#include <%s> \n int main() { return 0; }' % header
    return gccc._try_code(compiledata) == 0

def confdict_has_header(gccc, confdict, header):
    has_hdr = has_header(gccc,header+".h")
    if has_hdr:
        macro = "HAVE_%s_H" % header.upper().replace("/","_")
        confdict[macro] = None

def create_config_file(gccc, config, confdict = { }):
    #We use stdio as the stereotypical example of a stdc header
    confdict_has_header(gccc, confdict, "stdio")
    if confdict.has_key("stdio"):
        confdict["HAVE_STDIO_H"] = None
    #Header defines
    confdict_has_header(gccc, confdict, "alloca")
    confdict_has_header(gccc, confdict, "complex")
    confdict_has_header(gccc, confdict, "crt_externs")
    confdict_has_header(gccc, confdict, "dirent")
    confdict_has_header(gccc, confdict, "dlfcn")
    confdict_has_header(gccc, confdict, "fenv")
    confdict_has_header(gccc, confdict, "float")
    confdict_has_header(gccc, confdict, "inttypes")
    confdict_has_header(gccc, confdict, "langinfo")
    confdict_has_header(gccc, confdict, "limits")
    confdict_has_header(gccc, confdict, "locale")
    confdict_has_header(gccc, confdict, "malloc")
    confdict_has_header(gccc, confdict, "memory")
    confdict_has_header(gccc, confdict, "pwd")
    confdict_has_header(gccc, confdict, "sched")
    confdict_has_header(gccc, confdict, "stddef")
    confdict_has_header(gccc, confdict, "stdint")
    confdict_has_header(gccc, confdict, "stdlib")
    confdict_has_header(gccc, confdict, "strings")
    confdict_has_header(gccc, confdict, "string")
    confdict_has_header(gccc, confdict, "sys/param")
    confdict_has_header(gccc, confdict, "sys/poll")
    confdict_has_header(gccc, confdict, "sys/resource")
    confdict_has_header(gccc, confdict, "sys/select")
    confdict_has_header(gccc, confdict, "sys/stat")
    confdict_has_header(gccc, confdict, "sys/times")
    confdict_has_header(gccc, confdict, "sys/time")
    confdict_has_header(gccc, confdict, "sys/types")
    confdict_has_header(gccc, confdict, "sys/wait")
    confdict_has_header(gccc, confdict, "unistd")
    confdict_has_header(gccc, confdict, "values")

    #size of types
    sizeof_char = get_type_size(gccc, "char")
    sizeof_int = get_type_size(gccc, "int")
    sizeof_long = get_type_size(gccc, "long")
    sizeof_long_long = get_type_size(gccc, "long long")
    sizeof_short = get_type_size(gccc, "short")
    sizeof_void_p = get_type_size(gccc, "void*")
    sizeof___int64 = get_type_size(gccc, "__int64")

    #For size_t we first try to calculate without sys/types, then with if necessary
    sizeof_size_t = get_type_size(gccc, "size_t")
    if sizeof_size_t == 0 and confdict.has_key("HAVE_SYS_TYPES_H"):
        compiledata = '#include <sys/types.h> \n int main() { return sizeof(size_t); }'
        sizeof_size_t = gccc._try_code(compiledata)
        if sizeof_size_t==-1:
            sizeof_size_t = 0

    confdict["SIZEOF_CHAR"] = sizeof_char
    confdict["SIZEOF_INT"] = sizeof_int
    confdict["SIZEOF_LONG"] = sizeof_long
    confdict["SIZEOF_LONG_LONG"] = sizeof_long_long
    confdict["SIZEOF_SHORT"] = sizeof_short
    confdict["SIZEOF_VOID_P"] = sizeof_void_p
    confdict["SIZEOF_SIZE_T"] = sizeof_size_t

    #Do we have inline capeabilities?
    compiledata = 'inline int foo() { return 0; } \n int main() { return sizeof(size_t); }'
    ret = gccc._try_code(compiledata) == 0
    if ret!= -1:
        confdict["HAVE_INLINE"] = 1

    compiledata = '__inline int foo() { return 0; } \n int main() { return sizeof(size_t); }'
    ret = gccc._try_code(compiledata) == 0
    if ret!= -1:
        confdict["HAVE___INLINE"] = 1

    compiledata = '__inline__ int foo() { return 0; } \n int main() { return sizeof(size_t); }'
    ret = gccc._try_code(compiledata) == 0
    if ret!= -1:
        confdict["HAVE___INLINE__"] = 1

    compiledata = 'int main() { int x = 1; if (*(char*)&x==1) return 0; else return 1; }'
    ret = gccc._try_code(compiledata) == 0
    if ret == 0:
        confdict["BIG_ENDIAN"] = None
    else:
        confdict["LITTLE_ENDIAN"] = None

def output_config_file(config, confdict):
    config.write("#ifndef _CONFIG_H_\n")
    config.write("#define _CONFIG_H_\n")
    printlist = []
    for x in confdict.iteritems():
        outstr = "#define %s" % x[0]
        if x[1]!=None:
            outstr = outstr + " %s\n" % x[1]
        else:
            outstr = outstr + "\n"
        printlist.append(outstr)
    printlist.sort()
    for x in printlist:
        config.write(x)

def write_config_cache(filenm, confdict):
    output = open(filenm, 'wb')
    pickle.dump(confdict, output, 2)
    output.close()

def read_config_cache(filenm):
    confdict = None
    try:
        pkl_file = open(filenm, 'rb')
        confdict = pickle.load(pkl_file)
        pkl_file.close()
    except IOError:
        return None
    if _verify_dict(confdict)==False:
        raise EnvironmentError, "Fatal configuration error, rebuild all necessary"
    return confdict

def _verify_dict(confdict):
    if confdict["SIZEOF_VOID_P"] == confdict["SIZEOF_INT"]:
        #32-bit
        if sys.maxint == 9223372036854775807:
            return False
    else:
        #64-bit
        if sys.maxint != 9223372036854775807:
            return False
    if confdict.has_key("BIG_ENDIAN"):
        if sys.byteorder!='big':
            return False
    elif confdict.has_key("LITTLE_ENDIAN"):
        if sys.byteorder!='little':
            return False
    return True