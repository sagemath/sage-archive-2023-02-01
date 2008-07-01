##########################################################################
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#
#    * Neither the name of Gary Furnish nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##########################################################################
import os, signal, sys, time, thread, threading, tempfile, subprocess
from build.compiler import Compiler, Extension

def config_gcc_file(env, dict, path, language = 'C', define_macros = list(), include_dirs = list(), library_dirs = list(), libraries = list(), prop_options = { }):
    """
        This is a utility function for quickly modifying an existing gcc extension in a dictionary by pathname.
        It extends existing values (for lists) and updates (for dictionaries)
        INPUT:
            env -- <Enviroment> The enviroment used for this Extension
            dict -- <dict> The dictionary used for pathname lookup
            path -- <string> The pathname of the file we wish to modify
            language -- <string> (default: 'C') A value of either 'C' or 'C++' dictates what the output of this extension should be.
            define_macros -- <list> (default: list() ) A list of define macros that should be passed to the C/C++ compiler
            libraries -- <list> (default: list() ) A list of libraries that should be passed to the C/C++ linker
            include_dirs -- <list> (default: list() ) A list of include directories that should be passed to the C/C++ compiler
            library_dirs -- <list> (default: list() ) A list of library directories that should be passed to the C/C++ linker
            prop_options -- <dict> A dictionary that is used as a base dictionary for any Extensions that use this Extension as a source file
    """
    target = dict[path]
    target.language = language
    target.define_macros.extend(define_macros)
    target.include_dirs.extend(include_dirs)
    target.library_dirs.extend(library_dirs)
    target.libraries.extend(libraries)
    target.prop_options.update(prop_options)

class GCC_compiler(Compiler):
    def __init__(self, env, define_macros = list(), libraries=list(), include_dirs = list(), library_dirs = list(), options = { }):
        Compiler.__init__(self)
        self.options = dict(options)
        try:
            self.include_dirs = env.get_default('include_dirs') + include_dirs
        except:
            self.include_dirs = list() + include_dirs
        try:
            self.library_dirs = env.get_default('library_dirs') + library_dirs
        except:
            self.library_dirs = list() + library_dirs
        try:
            self.libraries = env.get_default('libraries') + libraries
        except:
            self.libraries = list() + libraries

        self.define_macros = list(define_macros)

    def get_file_extensions(self):
        return list('.c','.cc','cpp')

    def _init_extension(self, extension):
        self.mutex.acquire()
        extension.extmutex.acquire()
        newdict = dict(self.options)
        newdict.update(extension.options)
        extension.options = newdict
        extension.libraries = self.libraries+extension.libraries
        extension.define_macros=self.define_macros+extension.define_macros
        extension.library_dirs = self.library_dirs + extension.library_dirs
        extension.include_dirs = self.include_dirs + extension.include_dirs
        extension.extmutex.release()
        self.mutex.release()

    def _try_code(self, code):
        """
            This function attempts to compile and run the given code.
            OUTPUT:
                <integer> -- 0 on success
        """
        code = code + '\n'
        outfile = tempfile.NamedTemporaryFile(mode='r')
        infile = tempfile.NamedTemporaryFile(suffix = '.c')
        dumpfile = tempfile.TemporaryFile()
        infilenm = infile.name
        outfilenm = outfile.name
        outfile.close()
        infile.write(code)
        infile.flush()
        gccarg = "%s -o %s" % (infilenm,'%s' % outfilenm)
        ret = subprocess.call(["gcc", infilenm, "-o", outfilenm],stdout=dumpfile,stderr=dumpfile)
        if ret!=0:
            return -1
        ret = int(os.system("%s" % outfilenm)/256)
        return ret

    def get_ptr_len(self):
        """
            This function returns the length of a pointer on the given platform.
            It uses sizeof and is guarenteed correct
            OUTPUT:
                <integer> -- sizeof a pointer in bytes
        """
        compiledata = 'int main() { return sizeof(void*); }'
        return self._try_code(compiledata)

    def is_64bit(self):
        """
            This function returns if the platform is 64 bit
            It uses sizeof and is guarenteed correct
            OUTPUT:
                <bool> -- is the platform 64 bit
        """
        return self.get_ptr_len()==8

    def is_32bit(self):
        """
            This functoin returns if the platform is 32 bit
            It uses sizeof and is guarenteed correct
            OUTPUT:
                <bool> -- is the platform 32 bit
        """
        return self.get_ptr_len()==4


class GCC_extension(Extension):
    def __init__(self, compiler, env, sources=list(), outfile = None, language='C', define_macros = list(), libraries=list(), include_dirs = list(), library_dirs = list(), options = { }, prop_options = { },):
        Extension.__init__(self, compiler, env, sources, outfile, options, prop_options)
        self.language = str(language)
        self.define_macros=list(define_macros)
        self.libraries = list(libraries)
        self.include_dirs = list(include_dirs)
        self.library_dirs = list(library_dirs)
        compiler._init_extension(self)
    def _get_gcc_flags(self,env):
        ret = ""
        return ret

    def generate_command(self):
        self.extmutex.acquire()
        if self.language == 'C':
            cmd = 'gcc'
        elif self.language == 'C++':
            cmd = 'g++'
        else:
            raise TypeError
        for x in self.options.iteritems():
            if x[1]==None:
                cmd+=' '+x[0]
            else:
                cmd+=' '+x[0]+'='+x[1]
        for src in self.sources:
            cmd = cmd + ' '  + src
        cmd +=self._get_gcc_flags(self.env)
        if self.outfile != None:
            cmd+= ' -o %s' %self.outfile
        self.extmutex.release()
        return cmd

class GCC_extension_object(GCC_extension):
    def __init__(self, compiler, env, sources, outdir, language='C', define_macros = list(), libraries=list(), include_dirs = list(), library_dirs = list(), options = { }, prop_options = { }, outfile = None):
        newoptions = dict(options)
        newoptions["-c"] = None
        options = dict(options)
        language = str(language)
        include_dirs = list(include_dirs)
        define_macros = list(define_macros)
        libraries = list(libraries)
        library_dirs = list(library_dirs)
        if isinstance(sources[0], Extension):
            relfile = os.path.split(sources[0].outfile)[1]
            try:
                options.update(sources[0].prop_options[GCC_Extension_Object])
            except:
                pass
            try:
                language = sources[0].language
            except:
                pass
            try:
                include_dirs.extend(sources[0].include_dirs)
            except:
                pass
            try:
                define_macros.extend(sources[0].define_macros)
            except:
                pass
            try:
                libraries.extend(sources[0].libraries)
            except:
                pass
            try:
                library_dirs.extend(source[0].library_dirs)
            except:
                pass
        else:
            relfile = os.path.split(sources[0])[1]
        if outfile==None:
            if outdir[len(outdir)-1]!='/':
                outdir = outdir+'/'
            self.outdir = str(outdir)
            outfile = outdir+relfile
            outfile = outfile.replace(os.path.splitext(relfile)[1],".o")
        GCC_extension.__init__(self,compiler,env,sources,outfile,language,define_macros,libraries,include_dirs,library_dirs, newoptions, prop_options)
    def _get_gcc_flags(self,env):
        ret = ""
        for macro in self.define_macros:
            if macro[1]!=None:
                ret+= (" -D"+macro[0]+"="+macro[1]+" ")
            else:
                ret+= (" -D"+macro[0]+" ")
        for dir in self.include_dirs:
            ret+= (" -I"+dir+" ")
        if env.options['UNAME']=="Darwin" and sys.maxint == 9223372036854775807:
            ret+= (" -m64 ")
        return ret
    def get_out_dir(self):
        self.extmutex.acquire()
        ret = self.outdir
        self.extmutex.release()
        return ret

class GCC_extension_shared_object(GCC_extension):
    def __init__(self, compiler, env, sources, outdir, language='C', define_macros = list(), libraries=list(), include_dirs = list(), library_dirs = list(), options = { }, prop_options = { }, outfile = None):
        newoptions = dict(options)
        if env.options['UNAME']=='Darwin':
            newoptions["-dynamiclib"] = None
        else:
            newoptions["-shared"] = None
        options = dict(options)
        libraries = list(libraries)
        library_dirs = list(library_dirs)
        if isinstance(sources[0], Extension):
            relfile = os.path.split(sources[0].outfile)[1]
            try:
                options.update(sources[0].prop_options[GCC_Extension_Linker])
            except:
                pass
            try:
                libraries.extend(sources[0].libraries)
            except:
                pass
            try:
                library_dirs.extend(source[0].library_dirs)
            except:
                pass
        else:
            relfile = os.path.split(sources[0])[1]
        try:
            libraries = env.get_default('prelibraries') + libraries + env.get_default('postlibraries')
        except:
            pass
        if env.options['UNAME']=="Darwin":
            newoptions.update({"-single_module":None,  "-flat_namespace":None,  "-undefined dynamic_lookup":None })
            if sys.maxint == 9223372036854775807:
                newoptions.update({"-m64":None})
        if outfile==None:
            if outdir[len(outdir)-1]!='/':
                outdir = outdir+'/'
            self.outdir = str(outdir)
            outfile = outdir+relfile
            if env.options['UNAME']=="Darwin":
                outfile = outfile.replace(os.path.splitext(relfile)[1],".dylib")
            else:
                outfile = outfile.replace(os.path.splitext(relfile)[1],".so")
        GCC_extension.__init__(self,compiler,env,sources,outfile,language,define_macros,libraries,include_dirs,library_dirs, newoptions, prop_options)
    def _get_gcc_flags(self,env):
        ret = ""
        for dir in self.library_dirs:
            ret+= (" -L"+dir+" ")
        for lib in self.libraries:
            ret+= (" -l"+lib+" ")
        return ret
    def get_out_dir(self):
        self.extmutex.acquire()
        ret = self.outdir
        self.extmutex.release()
        return ret

class GCC_extension_executable(GCC_extension):
    def __init__(self, compiler, env, sources, outfile, language='C', define_macros = list(), libraries=list(), include_dirs = list(), library_dirs = list(), options = { }, prop_options = { }):
        newoptions = dict(options)
        options = dict(options)
        libraries = list(libraries)
        library_dirs = list(library_dirs)
        if isinstance(sources[0], Extension):
            relfile = os.path.split(sources[0].outfile)[1]
            try:
                options.update(sources[0].prop_options[GCC_Extension_Linker])
            except:
                pass
            try:
                libraries.extend(sources[0].libraries)
            except:
                pass
            try:
                library_dirs.extend(source[0].library_dirs)
            except:
                pass
        else:
            relfile = os.path.split(sources[0])[1]
        try:
            libraries = env.get_default('prelibraries') + libraries + env.get_default('postlibraries')
        except:
            pass
        GCC_extension.__init__(self,compiler,env,sources,outfile,language,define_macros,libraries,include_dirs,library_dirs, newoptions, prop_options)
    def _get_gcc_flags(self,env):
        ret = ""
        for dir in self.library_dirs:
            ret+= (" -L"+dir+" ")
        for lib in self.libraries:
            ret+= (" -l"+lib+" ")
        return ret
    def get_out_dir(self):
        self.extmutex.acquire()
        ret = self.outdir
        self.extmutex.release()
        return ret
