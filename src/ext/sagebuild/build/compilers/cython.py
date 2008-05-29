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

import os, signal, sys, time, thread, threading, tempfile
from build.compiler import Compiler, Extension

def _get_cython_file(curfile, modulename,includelist):
    """
        Get the pxd file for a given cimported module name
        OUTPUT:
            <string> -- pathname of the pxd
    """
    splitname = modulename.split(".")
    dirname = os.path.dirname(curfile)
    if len(splitname)==1:
        newpath = dirname+'/'+splitname[0]+'.pxd'
        return newpath

    newpath = includelist[0]
    for x in range(0,len(splitname)):
        newpath = newpath+'/'+splitname[x]
    newpath = newpath + '.pxd'
    return newpath

def _parse_file(filenm, includelist):
    """
        For a given cython file, get a list of all compile time dependencies
        INPUT:
            filenm -- <string> pathname for the file to analyze
            includelist -- <list> a list of directories to look for dependencies in
        OUTPUT:
            <list> -- A list of the pathnames of all dependencies
    """
    f = open(filenm, "r")
    deps = [ ]
    for line in f:
        splitline = line.replace(",","").split()
        if len(splitline)==0:
            continue
        elif splitline[0] == "cimport":
            for i in range(1,len(splitline)):
                deps.append(_get_cython_file(filenm,splitline[i],includelist))
        elif len(splitline)==2 and splitline[0]=="include":
            splitline[1] = splitline[1].replace('"','')
            if splitline[1].endswith(".pxi")  or splitline[1].endswith(".pxd"):
                curdir = os.path.dirname(filenm)
                curdirfile = curdir+'/'+splitline[1]
                if os.path.exists(curdirfile):
                    deps.append(curdirfile)
                else:
                    found = False
                    for incdir in includelist:
                        curdirfile = incdir+'/'+splitline[1]
                        if os.path.exists(curdirfile):
                            deps.append(curdirfile)
                            found = True
                            break
                    if found == False:
                        print splitline[1] + " not found"
                        raise Exception
        elif len(splitline)>3 and splitline[0]=="from" and splitline[2]=="cimport":
            deps.append(_get_cython_file(filenm,splitline[1],includelist))
    f.close()
    return deps

def get_cython_file_deps(filenm, includelist):
    """
        Given a cython file, also get any dependencies for an associated pxd
        INPUT:
            filenm -- <string> pathname for the file to analyze
            includelist -- <list> a list of directories to look for dependencies in
        OUTPUT:
            <list> -- A list of the pathnames of all dependencies
    """
    deps = [ ]
    deps.extend(_parse_file(filenm, includelist))
    try:
        pxdfilenm = filenm.replace("pyx","pxd")
        if pxdfilenm != filenm and os.path.exists(pxdfilenm):
            deps.extend(_parse_file(pxdfilenm, includelist))
            deps.append(pxdfilenm)
    except:
        pass
    return list(set(deps))



class CythonExtension(Extension):
    """
        This class, modeled after the distutils extension, represents a
        single Cython compilation
    """
    def __init__(self, compiler, env, sources=list(), language='C', define_macros = list(), libraries=list(), include_dirs = list(), library_dirs = list(), options = { }, prop_options = { }, cwd = os.getcwd()):
        """
            Create a Cython extension class
            INPUT:
                compiler -- <Compiler> The compiler for this Extension
                env -- <Enviroment> The enviroment used for this Extension
                sources -- <list> A list of source files (either pathnames or other extensions) for this extension
                language -- <string> (default: 'C') A value of either 'C' or 'C++' dictates what the output of this extension should be.
                define_macros -- <list> (default: list() ) A list of define macros that should be passed to the C/C++ compiler
                libraries -- <list> (default: list() ) A list of libraries that should be passed to the C/C++ linker
                include_dirs -- <list> (default: list() ) A list of include directories that should be passed to the C/C++ compiler
                library_dirs -- <list> (default: list() ) A list of library directories that should be passed to the C/C++ linker
                options -- <dict> A dictionary giving a mapping between options and their values for the Extension
                prop_options -- <dict> A dictionary that is used as a base dictionary for any Extensions that use this Extension as a source file
                cwd -- <string> (default: os.getcwd()) A path that is used as the working directory for Cython
        """
#        outfile = os.path.split(sources[0])[1]
        outfile = sources[0].replace(".pyx",".c")
        Extension.__init__(self, compiler, env, sources, outfile, options, prop_options)
        self.define_macros=list(define_macros)
        self.libraries = list(libraries)
        self.include_dirs = list(include_dirs)
        self.library_dirs = list(library_dirs)
        self.cwd = str(cwd)

    def _get_cython_flags(self,env):
        self.extmutex.acquire()
        ret = ""
        try:
            include_dirs = env.get_default('include_dirs')+self.include_dirs
        except:
            include_dirs = self.include_dirs
        for dir in include_dirs:
            ret+= ("-I"+dir+" ")
        self.extmutex.release()
        return ret

    def generate_command(self):
        """
            Create the command to be executed to call Cython
        """
        self.extmutex.acquire()
        file = (self.sources[0])[len(self.cwd)+1:]
        fileout = (self.outfile)[len(self.cwd)+1:]
        cmd = 'cython'
        for x in self.options.iteritems():
            if x[1]==None:
                cmd+=' '+x[0]
            else:
                cmd+=' '+x[0]+'='+x[1]
        if True:
            cmd+= ' -o %s' %fileout
        cmd = cmd + ' '  + file
        cmd = cmd + ' -w ' + self.cwd
        self.extmutex.release()
        return cmd
    def get_cwd(self):
        """
            Get the current working directory of Cython
        """
        self.extmutex.acquire()
        ret = self.cwd
        self.extmutex.release()
        return ret


class Cython_compiler(Compiler):
    """
        This class represents the Cython compiler
    """
    def __init__(self):
        Compiler.__init__(self)
        #self.options = { "--embed-positions":None, "--incref-local-binop":None }
        self.options = { "--embed-positions":None }
    def get_file_extensions(self):
        """
            Get the file extensions that this compiler can handle
        """
        return '.pyx'
    def get_build_extensions(self,env, filelist, cwd = os.getcwd()):
        """
            Given a list of files, return a list of cython extensions for it
            INPUT:
                env -- <Enviroment> enviroment for the build extensions
                filelist -- <list>
            OUTPUT:
                <list> - A list of CythonExtensions
        """
        self.mutex.acquire()
        ret = { }
        for file in filelist:
            ret[file] = CythonExtension(self, env, [file], options = self.options, cwd=cwd)
        self.mutex.release()
        return ret

