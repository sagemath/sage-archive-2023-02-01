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

class filewalker:
    def __init__(self):
        pass
    def walk(self, path):
        if os.path.isfile(path):
            if _test(str(path)):
                return [path]
            else:
                return list()
        curlist = list()
        if os.path. isabs(path):
            walkdir = path
        else:
            walkdir = os.path.join(os.getcwd(),path)
        for root, dirs, lfiles in os.walk(walkdir):
            for F in lfiles:
                appendstr = os.path.realpath(os.path.join(root,F))
                if(self._test(appendstr)):
                    curlist.append(appendstr)
            for D in dirs:
                appendstr = os.path.realpath(os.path.join(root,appendstr))
                if(self._test(appendstr)):
                    curlist.append(appendstr)
        return curlist

    def _test(self,x):
        return False

class extfilewalker(filewalker):
    def __init__(self):
        filewalker.__init__(self)
        self.callbacks = { }
        self.callbackmutex = thread.allocate_lock()
    def addcallback(self, ext, callback):
        self.callbackmutex.acquire()
        self.callbacks[ext]=callback
        self.callbackmutex.release()
    def _test(self,x):
        if(os.path.isdir(x)):
            return False
        self.callbackmutex.acquire()
        try:
            ext = os.path.splitext(x)[1]
            ret = self.callbacks[ ext ]
            self.callbackmutex.release()
            return ret(x)
        except:
            pass
        self.callbackmutex.release()
        return False
