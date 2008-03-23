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

TM = None

class taskthread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        try:
            wasexecuting = False
            while True:
                act = TM._get_action(wasexecuting)
                if act == None:
                    wasexecuting=False
                else:
                    act.execute()
                    act.post_process()
                    wasexecuting=True
        except KeyboardInterrupt:
           TM.exit()

class taskmanager:
    def __init__(self, numthreads):
        TM = self
        self.numaction = 0
        self.numthreads = numthreads
        self.numexecuting = 0
        self.queue = list()
        self.queuemutex = thread.allocate_lock()
        self.printmutex = thread.allocate_lock()
        self.available = threading.Event()
        self._threads = list()
        self._finishtasks = list()

    def _action_become_ready(self, act):
        self.queuemutex.acquire()
        self.queue.insert(0,act)
#        if(len(self.queue)==1):
#            self.available.set()
        self.queuemutex.release()
    def _get_action(self, wasexecuting):
#        self.available.wait(.01)
        self.queuemutex.acquire()
        if wasexecuting:
            self.numexecuting-=1
        if len(self.queue)>0:
            ret = self.queue.pop()
            self.numexecuting+=1
        else:
            ret = None
        if self.numexecuting<=0:
            if ret!=None:
                raise TypeError
            self.queuemutex.release()
            self.available.set()
            thread.exit()
#        if (len(self.queue)==0) and (ret==None) and (self.numaction!=0):
 #           self.available.clear()
        self.queuemutex.release()
        return ret
    def _init_action(self,act):
        self.numaction+=1
    def _destroy_action(self,act):
        self.numaction-=1
    def _init_threads(self):
        for i in range(0,self.numthreads):
            newthread = taskthread()
            self._threads.append(newthread)
            newthread.start()

    def put(self, x):
        self.printmutex.acquire()
        print x
        self.printmutex.release()
    def putsafe(self, x):
        print x
    def startput(self):
        self.printmutex.acquire()
    def endput(self):
        self.printmutex.release()
    def addfinishtask(self,task):
        self.queuemutex.acquire()
        self._finishtasks.append(task)
        self.queuemutex.release()
    def go(self):
        self._init_threads()
        try:
            for t in self._threads:
                t.join()
        except:
            sys.exit(0)
        for x in self._finishtasks:
            x()
        sys.exit(0)

    def exit(self,x=0):
        thread.exit()

