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
from build.action import Action

class Compiler:
    def __init__(self):
        self.action = None
        self.options = { }
        self.mutex = thread.allocate_lock()
    def get_file_extensions():
        return None
    def get_option_val(self, val):
        self.mutex.acquire()
        try:
            ret = self.options[key]
            self.mutex.release()
            return ret
        except:
            self.mutex.release()
            return None

    def set_option_val(self, key, val=True):
        self.mutex.acquire()
        self.options[key]=val
        self.mutex.release()

class Compiler_action(Action):
    def __init__(self, extension, env):
        Action.__init__(self)
        self.extension = extension
        self.env = env
    def execute(self):
        try:
            filedir = os.path.realpath(self.extension.get_cwd())
            os.chdir(filedir)
        except:
            pass
        oldcmd = cmd = self.extension.generate_command()
        outfile = tempfile.NamedTemporaryFile()
        errfile = tempfile.NamedTemporaryFile()
        cmd = 'bash -c "%s > %s 2> %s" ' %(cmd, outfile.name, errfile.name)
        ret = os.system(cmd)
        ol = outfile.read()
        el = errfile.read()
        import build.taskmanager
        TM = build.taskmanager.TM
        TM.startput()
        if ret!=0:
            TM.putsafe(el)
            TM.putsafe(ol)
        TM.putsafe(oldcmd)
        TM.endput()
        return ret

def ext_list_to_primary_src(extnlist):
    srclist = [ ]
    for x in extnlist:
        srclist.append(x.sources[0])
    return srclist

class Extension:
    def __init__(self, compiler, env, sources, outfile, options, prop_options):
        self.extmutex = thread.allocate_lock()
        self.compiler = compiler
        self.env = env
        self.sources = list(sources)
        for i in range(0,len(self.sources)):
            if isinstance(self.sources[i], Extension):
                self.sources[i] = self.sources[i].outfile
        self.outfile = str(outfile)
        self.options = dict(options)
        try:
            self.options.update(sources[0].prop_options[str(self.__class__)])
        except:
            pass
        self.prop_options = dict(prop_options)
    def generate_action(self,env):
        try:
            return self._action
        except:
            self._action = Compiler_action(self, env)
            return self._action
    def get_option_val(self, val):
        self.extmutex.acquire()
        try:
            ret = self.options[key]
            self.extmutex.release()
            return ret
        except:
            self.extmutex.release()
            return None

    def set_option_val(self, key, val=True):
        self.extmutex.acquire()
        self.options[key]=val
        self.extmutex.release()