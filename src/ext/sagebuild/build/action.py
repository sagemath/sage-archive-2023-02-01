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

import taskmanager

class Action:
    def __init__(self):
        self._dependent = list()
        self._children = list()
    def __del__(self):
        #taskmanager.TM._destroy_action(self)
        pass
    def dep_register(self, dep):
        self._dependent.append(dep)
        dep._child_register(self)

    def _child_register(self, child):
        self._children.append(child)

    def _child_unregister(self, child):
        self._children.remove(child)

    def dep_succeed(self, dep):
        self._dependent.remove(dep)
        if len(self._dependent) == 0:
            taskmanager.TM._action_become_ready(self)

    def dep_fail(self, dep):
        for x in self._dependent:
            if x!=dep:
                x._child_unregister(self)
        taskmanager.TM._action_unregister(self)
    def enable(self):
        taskmanager.TM._init_action(self)
        if len(self._dependent)==0:
            taskmanager.TM._action_become_ready(self)
    def execute(self):
        pass
    def post_process(self):
        for act in self._children:
            act.dep_succeed(self)

