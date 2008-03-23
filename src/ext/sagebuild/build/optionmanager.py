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
import os, signal, sys, time, thread, threading, tempfile, getopt

class argument:
    def __init__(self, helpstr):
        self.helpstr = helpstr

class option:
    def __init__(self, longstr, args, func, helpstr, shortstr = None, help_levels = [0], greedy = False):
        self.longstr = longstr
        self.args = args
        self.func = func
        self.helpstr = helpstr
        self.shortstr = shortstr
        self.help_levels = help_levels
        self.greedy = greedy

class optionmanager:
    def __init__(self):
        self.targets = { }
        self.abriv = { }
        self.default = None
        self.optionmutex = thread.allocate_lock()
    def set_default(self, target):
        self.optionmutex.acquire()
        self.default = target
        self.optionmutex.release()
    def add_option(self, opt):
        if not isinstance(opt, option):
            raise TypeError, "Option must be of type option"
        self.optionmutex.acquire()
        self.targets[opt.longstr] = opt
        if opt.shortstr != None:
            self.abriv[opt.shortstr] = opt
        self.optionmutex.release()

    def execute_options(self, lst, env):
        self.optionmutex.acquire()
        shortargs = ''
        for target in self.abriv.values():
            shortargs = shortargs+target.shortstr
        optlist, args = getopt.getopt(lst, shortargs, self.targets.keys())
        argnum = 0
        anyfound = False
        for x in optlist:
            try:
                if self.abriv.has_key((x[0])[1:]):
                    target = self.abriv[ (x[0])[1:] ]
                elif self.targets.has_key((x[0])[2:]):
                    target = self.targets[ (x[0])[2:] ]
                else:
                    self.optionmutex.release()
                    return False
                localargs = []
                if x[1]!="":
                    localargs.append(x[1])
                if len(args)>len(localargs):
                    grabargs = target.len(args)-len(localargs)
                    localargs.extend(args[argnum:argnum+grabargs ] )
                    argnum=argnum+grabargs
                if target.greedy:
                    localargs.extend( args[argnum:] )
                    argnum = len(args)
                f = target.func
                self.optionmutex.release()
                ret =f(env,localargs)
                if ret == False:
                    print "Bad arguments to " + str(x)
                    return ret
                anyfound = True
                self.optionmutex.acquire()
            except None:
                if self.optionmutex.locked():
                    self.optionmutex.release()
                return False
        if anyfound == False:
            f = self.default
            self.optionmutex.release()
            f(env,args[argnum:])
        else:
            self.optionmutex.release()

        return True

    def helpstr(self, level = 0):
        outtxt = ""
        for target in self.targets:
            if level in target.help_levels:
                outtxt = outtxt + target.longstr + " "
                if target.shortstr != None:
                    outtxt = outtxt + shortstr + " "
                for arg in target.args:
                    outtxt = outtxt + arg.helpstr + " "
                outtxt = outtxt + target.helpstr + " "
                outtxt = outtxt + "\n"
        return outtxt


OM = optionmanager()



